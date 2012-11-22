//
// Copyright (c) 2012,
// Florent Lamiraux
//
// CNRS
//
// This file is part of sot-dynamic.
// sot-dynamic is free software: you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.
// sot-dynamic is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.  You should
// have received a copy of the GNU Lesser General Public License along
// with sot-dynamic.  If not, see <http://www.gnu.org/licenses/>.
//

#include <dynamic-graph/linear-algebra.h>
#include <dynamic-graph/factory.h>
#include <dynamic-graph/signal-time-dependent.h>
#include <dynamic-graph/signal.h>
#include <dynamic-graph/signal-ptr.h>
#include <dynamic-graph/command-setter.h>
#include <dynamic-graph/command-getter.h>
#include <dynamic-graph/command-direct-setter.h>
#include <dynamic-graph/command-direct-getter.h>
#include <dynamic-graph/command-bind.h>

#include <sot/core/task-abstract.hh>
#include <sot/core/matrix-rotation.hh>
#include <sot/core/matrix-homogeneous.hh>
#include <sot/core/multi-bound.hh>

namespace sot {
  namespace dynamic {
    using dynamicgraph::sot::TaskAbstract;
    using dynamicgraph::Signal;
    using dynamicgraph::SignalPtr;
    using dynamicgraph::SignalTimeDependent;
    using dynamicgraph::Vector;
    using dynamicgraph::Matrix;
    using dynamicgraph::Entity;
    using dynamicgraph::sot::VectorMultiBound;
    using dynamicgraph::command::makeCommandVoid0;
    using dynamicgraph::command::docCommandVoid0;
    using dynamicgraph::command::docDirectSetter;
    using dynamicgraph::command::makeDirectSetter;
    using dynamicgraph::command::docDirectGetter;
    using dynamicgraph::command::makeDirectGetter;
    using dynamicgraph::sot::MatrixHomogeneous;
    using dynamicgraph::sot::MatrixRotation;

    /// Dynamic balance stabilizer
    ///
    /// This task takes as input four signals
    /// \li comSIN, the position of the center of mass (COM)
    /// \li comDesSIN, the the desired position of the center of mass,
    /// \li zmpSIN, the position of the center of pressure (ZMP),
    /// \li zmpDesSIN, the desired position of the center of pressure,
    /// \li jacobianSIN, the jacobian of the center of mass,
    /// \li comdotSIN, reference velocity of the center of mass,
    /// and provides as output two signals
    /// \li taskSOUT, the desired time derivative of the center of mass,
    /// \li jacobianSOUT, the jacobian of the center of mass
    class Stabilizer : public TaskAbstract
    {
      DYNAMIC_GRAPH_ENTITY_DECL ();
    public:
      // Constant values
      static double m_;
      static double g_;
      static double zeta_;

      /// Constructor by name
      Stabilizer(const std::string& inName) :
	TaskAbstract(inName),
	deltaComSIN_ (NULL, "Stabilizer("+inName+")::input(vector)::deltaCom"),
	jacobianSIN_ (NULL, "Stabilizer("+inName+")::input(matrix)::Jcom"),
	comdotSIN_ (NULL, "Stabilizer("+inName+")::input(vector)::comdot"),
	leftFootPositionSIN_
	(NULL, "Stabilizer("+inName+")::input(vector)::leftFootPosition"),
	rightFootPositionSIN_
	(NULL, "Stabilizer("+inName+")::input(vector)::rightFootPosition"),
	forceRightFootSIN_
	(NULL, "Stabilizer("+inName+")::input(vector)::forceRLEG"),
	forceLeftFootSIN_
	(NULL, "Stabilizer("+inName+")::input(vector)::forceLLEG"),
	stateFlexRfxSIN_
	(NULL, "Stabilizer("+inName+")::input(vector)::stateFlex_rfx"),
	stateFlexRfySIN_
	(NULL, "Stabilizer("+inName+")::input(vector)::stateFlex_rfy"),
	stateFlexLfxSIN_
	(NULL, "Stabilizer("+inName+")::input(vector)::stateFlex_lfx"),
	stateFlexLfySIN_
	(NULL, "Stabilizer("+inName+")::input(vector)::stateFlex_lfy"),
	controlGainSIN_
	(NULL, "Stabilizer("+inName+")::input(double)::controlGain"),
	sideGainSIN_
	(NULL, "Stabilizer("+inName+")::input(double)::sideGain"),
	d2comSOUT_ ("Stabilizer("+inName+")::output(vector)::d2com"),
	debugSOUT_ ("Stabilize("+inName+")::debug"),
	gain1_ (4), gain2_ (4),
	prevCom_(3), flexAngle_ (2), flexDeriv_ (2),
	dcom_ (3), timePeriod_ (.005), on_ (false),
	forceThreshold_ (.036), angularStiffness_ (425.), d2com_ (3)
      {
	// Register signals into the entity.
	signalRegistration (deltaComSIN_);
	signalRegistration (jacobianSIN_);
	signalRegistration (comdotSIN_);
	signalRegistration (leftFootPositionSIN_);
	signalRegistration (rightFootPositionSIN_);
	signalRegistration (forceRightFootSIN_);
	signalRegistration (forceLeftFootSIN_);
	signalRegistration (stateFlexRfxSIN_);
	signalRegistration (stateFlexRfySIN_);
	signalRegistration (stateFlexLfxSIN_);
	signalRegistration (stateFlexLfySIN_);
	signalRegistration (controlGainSIN_);
	signalRegistration (sideGainSIN_);
	signalRegistration (d2comSOUT_);
	signalRegistration (debugSOUT_);

	taskSOUT.addDependency (deltaComSIN_);
	taskSOUT.addDependency (comdotSIN_);
	taskSOUT.addDependency (leftFootPositionSIN_);
	taskSOUT.addDependency (rightFootPositionSIN_);
	taskSOUT.addDependency (forceRightFootSIN_);
	taskSOUT.addDependency (forceLeftFootSIN_);
	taskSOUT.addDependency (stateFlexRfxSIN_);
	taskSOUT.addDependency (stateFlexRfySIN_);
	taskSOUT.addDependency (stateFlexLfxSIN_);
	taskSOUT.addDependency (stateFlexLfySIN_);
	taskSOUT.addDependency (controlGainSIN_);
	taskSOUT.addDependency (sideGainSIN_);

	jacobianSOUT.addDependency (jacobianSIN_);

	taskSOUT.setFunction (boost::bind(&Stabilizer::computeControlFeedback,
					  this,_1,_2));
	jacobianSOUT.setFunction (boost::bind(&Stabilizer::computeJacobianCom,
					      this,_1,_2));

	d2com_.fill (0.);
	dcom_.fill (0.);
	d2comSOUT_.setConstant (d2com_);
	sideGainSIN_.setConstant (1.);

	std::string docstring;
	docstring =
	  "\n"
	  "    Set sampling time period task\n"
	  "\n"
	  "      input:\n"
	  "        a floating point number\n"
	  "\n";
	addCommand("setTimePeriod",
		   new dynamicgraph::command::Setter<Stabilizer, double>
		   (*this, &Stabilizer::setTimePeriod, docstring));
	docstring =
	  "\n"
	  "    Get sampling time period task\n"
	  "\n"
	  "      return:\n"
	  "        a floating point number\n"
	  "\n";
	addCommand("getTimePeriod",
		   new dynamicgraph::command::Getter<Stabilizer, double>
		   (*this, &Stabilizer::getTimePeriod, docstring));

	addCommand ("start",
		    makeCommandVoid0 (*this, &Stabilizer::start,
				      docCommandVoid0 ("Start stabilizer")));

	addCommand ("setGain1",
		    makeDirectSetter (*this, &gain1_,
				      docDirectSetter
				      ("Set gains single support",
				       "vector")));

	addCommand ("getGain1",
		    makeDirectGetter (*this, &gain1_,
				      docDirectGetter
				      ("Get gains single support",
				       "vector")));

	addCommand ("setGain2",
		    makeDirectSetter (*this, &gain2_,
				      docDirectSetter
				      ("Set gains double support",
				       "vector")));

	addCommand ("getGain2",
		    makeDirectGetter (*this, &gain2_,
				      docDirectGetter
				      ("Get gains double support",
				       "vector")));

	prevCom_.fill (0.);
	flexAngle_.fill (0.);
	flexDeriv_.fill (0.);
	debug_.resize (11);

	gain1_ (0) = 177.57303317647063;
	gain1_ (1) = -29.735033684033631;
	gain1_ (2) = 54.413552941176476;
	gain1_ (3) = -27.530842352941178;

	gain2_ (0) = 82.655266588235293;
	gain2_ (1) = 36.712572443697468;
	gain2_ (2) = 27.206776470588238;
	gain2_ (3) = -5.76542117647059;
      }

      ~Stabilizer() {}

      /// Documentation of the entity
      virtual std::string getDocString ()	const
      {
	std::string doc =
	  "Dynamic balance humanoid robot stabilizer\n"
	  "\n"
	  "This task aims at controlling balance for a walking legged humanoid robot.\n"
	  "The entity takes 6 signals as input:\n"
	  "  - deltaCom: the difference between the position of the center of mass and the\n"
	  " reference,\n"
	  "  - Jcom: the Jacobian of the center of mass wrt the robot configuration,\n"
	  "  - comdot: the reference velocity of the center of mass \n"
	  "  \n"
	  "As any task, the entity provide two output signals:\n"
	  "  - task: the velocity of the center of mass so as to cope with\n"
	  "          perturbations,\n"
	  "  - jacobian: the Jacobian of the center of mass with respect to robot\n"
	  "              configuration.\n";
	return doc;
      }

      /// Start stabilizer
      void start () {
	on_ = true;
      }

      /// @}
      /// \name Sampling time period
      /// @{

      /// \brief Set sampling time period
      void setTimePeriod(const double& inTimePeriod)
      {
	timePeriod_ = inTimePeriod;
      }
      /// \brief Get sampling time period
      double getTimePeriod() const
      {
	return timePeriod_;
      }
      /// @}

    private:
      /// Compute flexibility state from both feet
      void computeFlexibility (const int& time)
      {
	const Vector& flexRfx = stateFlexRfxSIN_.access (time);
	const Vector& flexRfy = stateFlexRfySIN_.access (time);
	const Vector& flexLfx = stateFlexLfxSIN_.access (time);
	const Vector& flexLfy = stateFlexLfySIN_.access (time);
	const MatrixHomogeneous& Mr = rightFootPositionSIN_.access (time);
	const MatrixHomogeneous& Ml = leftFootPositionSIN_.access (time);
	const Vector& fr = forceRightFootSIN_.access (time);
	const Vector& fl = forceLeftFootSIN_.access (time);

	// Express vertical component of force in global basis
	double flz = Ml (2,0) * fl (0) + Ml(2,1) * fl (1) + Ml (2,2) * fl (2);
	double frz = Mr (2,0) * fr (0) + Mr(2,1) * fr (1) + Mr (2,2) * fr (2);

	nbSupport_ = 0;
	if (on_) {
	  if (frz >= forceThreshold_) nbSupport_++;
	  if (flz >= forceThreshold_) nbSupport_++;
	}
	if (frz < 0) frz = 0;
	if (flz < 0) flz = 0;
	double fz = flz + frz;
	if (fz == 0) {
	  flexAngle_ (0) = 0;
	  flexAngle_ (1) = 0;
	  flexDeriv_ (0) = 0;
	  flexDeriv_ (1) = 0;
	  return;
	}
	// Extract yaw from right foot position
	double nx = Mr (0,0);
	double ny = Mr (1,0);
	double norm = sqrt (nx*nx + ny*ny);
	double cth = nx/norm;
	double sth = ny/norm;
	double flexAngleRfx = cth * flexRfx (1) + sth * flexRfy (1);
	double flexAngleRfy = -sth * flexRfx (1) + cth * flexRfy (1);
	double flexDerivRfx = cth * flexRfx (3) + sth * flexRfy (3);
	double flexDerivRfy = -sth * flexRfx (3) + cth * flexRfy (3);
	// Extract yaw from left foot position
	nx = Ml (0,0);
	ny = Ml (1,0);
	norm = sqrt (nx*nx + ny*ny);
	cth = nx/norm;
	sth = ny/norm;
	double flexAngleLfx = cth * flexLfx (1) + sth * flexLfy (1);
	double flexAngleLfy = -sth * flexLfx (1) + cth * flexLfy (1);
	double flexDerivLfx = cth * flexLfx (3) + sth * flexLfy (3);
	double flexDerivLfy = -sth * flexLfx (3) + cth * flexLfy (3);

	flexAngle_ (0) = (frz * flexAngleRfx + flz * flexAngleLfx)/fz;
	flexAngle_ (1) = (frz * flexAngleRfy + flz * flexAngleLfy)/fz;
	flexDeriv_ (0) = (frz * flexDerivRfx + flz * flexDerivLfx)/fz;
	flexDeriv_ (1) = (frz * flexDerivRfy + flz * flexDerivLfy)/fz;
      }

      /// Compute the control law
      VectorMultiBound& computeControlFeedback(VectorMultiBound& comdot,
					       const int& time)
      {
	const Vector& deltaCom = deltaComSIN_ (time);
	const Vector& comdotRef = comdotSIN_ (time);
	const MatrixHomogeneous& leftFootPosition =
	  leftFootPositionSIN_.access (time);
	const MatrixHomogeneous& rightFootPosition =
	  rightFootPositionSIN_.access (time);
	const double& gain = controlGainSIN_.access (time);
	const double& sideGain = sideGainSIN_.access (time);

	computeFlexibility (time);

	double x = deltaCom (0);
	double y = deltaCom (1);
	double z = deltaCom (2);

	double theta0, dtheta0, norm, u2x, u2y, u1x, u1y,
	  lat, dlat, ddlat;
	double theta1, dtheta1, delta_x, delta_y, theta, dtheta, xi, dxi, ddxi;

	for (unsigned int i=0; i<11; ++i) {debug_ (i) = 0;}
	debug_ (10) = nbSupport_;

	switch (nbSupport_) {
	case 0:
	  dcom_ (0) = -gain * x;
	  dcom_ (1) = -gain * y;
	  break;
	case 1: //single support
	  //along x
	  theta0 = flexAngle_ (0);
	  dtheta0 = flexDeriv_ (0);
	  d2com_ (0)= -(gain1_ (0)*x + gain1_ (1)*theta0 + gain1_ (2)*dcom_ (0) +
			  gain1_ (3)*dtheta0);
	  debug_ (0) = x;
	  debug_ (1) = theta0;
	  debug_ (2) = dcom_ (0);
	  debug_ (3) = dtheta0;
	  debug_ (4) = d2com_ (0);
	  dcom_ (0) += timePeriod_ * d2com_ (0);
	  // along y
	  theta1 = flexAngle_ (1);
	  dtheta1 = flexDeriv_ (1);
	  d2com_ (1) = - (gain1_ (0)*y + gain1_ (1)*theta1 +
			  gain1_ (2)*dcom_ (1) + gain1_ (3)*dtheta1);
	  debug_ (5) = y;
	  debug_ (6) = theta1;
	  debug_ (7) = dcom_ (1);
	  debug_ (8) = dtheta1;
	  debug_ (9) = d2com_ (1);
	  dcom_ (1) += timePeriod_ * d2com_ (1);

	  break;
	case 2: //double support
	  // compute component of angle orthogonal to the line joining the feet
	  delta_x = leftFootPosition (0, 3) - rightFootPosition (0, 3);
	  delta_y = leftFootPosition (1, 3) - rightFootPosition (1, 3);
	  norm = sqrt (delta_x*delta_x+delta_y*delta_y);
	  u2x = delta_x/norm;
	  u2y = delta_y/norm;
	  u1x = u2y;
	  u1y = -u2x;
	  theta = - (u2x * flexAngle_ (0) + u2y * flexAngle_ (1));
	  dtheta = - (u2x * flexDeriv_ (0) + u2y * flexDeriv_ (1));
	  xi = u1x*x + u1y*y;
	  dxi = u1x*dcom_ (0) + u1y*dcom_ (1);
	  ddxi = - (gain2_ (0)*xi + gain2_ (1)*theta + gain2_ (2)*dxi +
		    gain2_ (3)*dtheta);
	  lat = u2x*x + u2y*y;
	  dlat = u2x*dcom_ (0) + u2y*dcom_ (1);
	  ddlat = -2*sideGain*dlat - sideGain*sideGain*lat;

	  d2com_ (0) = ddxi * u1x;
	  d2com_ (1) = ddxi * u1y;
	  dcom_ (0) += timePeriod_ * (d2com_ (0) + ddlat*u2x);
	  dcom_ (1) += timePeriod_ * (d2com_ (1) + ddlat*u2y);
	  break;
	default:
	  break;
	};
	dcom_ (2) = -gain * z;

	comdot.resize (3);
	comdot [0].setSingleBound (comdotRef (0) + dcom_ (0));
	comdot [1].setSingleBound (comdotRef (1) + dcom_ (1));
	comdot [2].setSingleBound (comdotRef (2) + dcom_ (2));

	d2comSOUT_.setConstant (d2com_);
	debugSOUT_.setConstant (debug_);
	debugSOUT_.setTime (time);
	return comdot;
      }

      Matrix& computeJacobianCom(Matrix& jacobian, const int& time)
      {
	typedef unsigned int size_t;
	jacobian = jacobianSIN_ (time);
	return jacobian;
      }

      /// Position of center of mass
      SignalPtr < dynamicgraph::Vector, int> deltaComSIN_;
      /// Position of center of mass
      SignalPtr < dynamicgraph::Matrix, int> jacobianSIN_;
      /// Reference velocity of the center of mass
      SignalPtr < dynamicgraph::Vector, int> comdotSIN_;
      // Position of left foot force sensor in global frame
      SignalPtr <MatrixHomogeneous, int> leftFootPositionSIN_;
      // Position of right foot force sensor in global frame
      SignalPtr <MatrixHomogeneous, int> rightFootPositionSIN_;
      // Force in right foot sensor
      SignalPtr <dynamicgraph::Vector, int> forceRightFootSIN_;
      // Force in left foot sensor
      SignalPtr <dynamicgraph::Vector, int> forceLeftFootSIN_;
      // Right foot flexibility state along local x axis
      SignalPtr <dynamicgraph::Vector, int> stateFlexRfxSIN_;
      // Right foot flexibility state along local y axis
      SignalPtr <dynamicgraph::Vector, int> stateFlexRfySIN_;
      // Left foot flexibility state along local x axis
      SignalPtr <dynamicgraph::Vector, int> stateFlexLfxSIN_;
      // Left foot flexibility state along local y axis
      SignalPtr <dynamicgraph::Vector, int> stateFlexLfySIN_;
      // Gain of center task of mass when stabilizer is off
      SignalPtr <double, int> controlGainSIN_;
      // Sideways gain in double support
      SignalPtr <double, int> sideGainSIN_;
      // Acceleration of center of mass
      SignalTimeDependent <dynamicgraph::Vector, int> d2comSOUT_;
      // Debug signal
      Signal <dynamicgraph::Vector, int> debugSOUT_;

      /// Gains single support
      Vector gain1_;      /// Gains double support
      Vector gain2_;
      /// Store center of mass for finite-difference evaluation of velocity
      Vector prevCom_;
      /// Angle of the flexibility
      Vector flexAngle_;
      Vector flexDeriv_;
      /// coordinates of center of mass velocity in moving frame
      Vector dcom_;
      // Time sampling period
      double timePeriod_;
      // Whether stabilizer is on
      bool on_;
      // Threshold on normal force above which the foot is considered in contact
      double forceThreshold_;
      // Angular stiffness of flexibility of each foot
      double angularStiffness_;
      // Number of feet in support
      unsigned int nbSupport_;
      // Acceleration of the center of mass computed by stabilizer
      Vector d2com_;
      // Debug
      Vector debug_;
    }; // class Stabilizer

    double Stabilizer::m_ = 56.;
    double Stabilizer::g_ = 9.81;
    double Stabilizer::zeta_ = .80;

    DYNAMICGRAPH_FACTORY_ENTITY_PLUGIN (Stabilizer, "Stabilizer");

    // This namespace contains a few classes related to the dynamics and
    // observation of the ankle flexibility of a humanoid robot
    //
    // The state of the flexibility is represented by vector
    //                   .    .
    //  x = (xi, theta, xi, theta, k     )
    //                              theta
    // where xi is the position of the center of mass in a moving frame rotating
    // around the contact foot,
    // theta is the rotation angle of the moving frame with respect to the
    // world frame
    // k      is the angular stiffness of the flexibility
    //  theta

    namespace flexibility {
      class Function : public Entity
      {
      protected:
	SignalPtr <double, int> stiffnessSIN_;
	SignalPtr < Vector, int > stateSIN_;
	double dt_;

      public:
	Function (const std::string& name) : Entity (name),
	  stiffnessSIN_ (0, "flexibility_f(" + name +
		       ")::input(double)::stiffness"),
	  stateSIN_ (0, "flexibility::Function(" + name +
		     ")::input(vector)::state"),
	  dt_ (.005)
	{
	  signalRegistration (stiffnessSIN_ << stateSIN_);
	  addCommand ("setTimePeriod",
		      makeDirectSetter (*this, &dt_,
					docDirectSetter ("time period",
							 "float")));
	  addCommand ("getTimePeriod",
		      makeDirectGetter (*this, &dt_,
					docDirectGetter ("time period",
							 "float")));
	}
      }; // class Function
      //
      // State transition of time discretized system
      //
      //   x    =  f (x , u )
      //    k+1        k   k
      //
      //   with
      //        ..
      //   u  = xi,   x  = x (k dt)
      //    k          k
      //
      class f : public Function
      {
	SignalPtr <Vector, int> controlSIN_;
	Signal <Vector, int> newStateSOUT_;
	Signal <Matrix, int> jacobianSOUT_;

	DYNAMIC_GRAPH_ENTITY_DECL();
	f (const std::string& name) :
	  Function (name),
	  controlSIN_ (0, "flexibility_f(" + name +
		       ")::input(vector)::control"),
	  newStateSOUT_ ("flexibility_f(" + name +
			 ")::output(vector)::newState"),
	  jacobianSOUT_ ("flexibility_f(" + name +
			 ")::output(vector)::jacobian")

	{
	  signalRegistration (controlSIN_ << newStateSOUT_
			      << jacobianSOUT_);
	  newStateSOUT_.setFunction (boost::bind (&f::computeNewState, this, _1, _2));
	  jacobianSOUT_.setFunction (boost::bind (&f::computeJacobian, this, _1, _2));
	}

	Vector& computeNewState (Vector& x, const int& time)
	{
	  double m = Stabilizer::m_;
	  double g = Stabilizer::g_;
	  double zeta = Stabilizer::zeta_;

	  const Vector& state = stateSIN_.accessCopy ();
	  const Vector& control = controlSIN_.access (time);
	  const double& kth = stiffnessSIN_.accessCopy ();

	  double xi = state (0);
	  double th = state (1);
	  double dxi = state (2);
	  double dth = state (3);

	  double u = control (0);
	  double d2 = (xi*xi+zeta*zeta);
	  x.resize (4);

	  x (0) = xi + dt_ * dxi;
	  x (1) = th + dt_ * dth;
	  x (2) = dxi + dt_ * u;
	  x (3) = dth + dt_* (-kth*th - m*g*(cos (th)*xi - sin (th)*zeta) +
			     m*(zeta*u -2*dth*xi*dxi))/(m*d2);
	  return x;
	}

	Matrix& computeJacobian (Matrix& J, const int&)
	{
	  double m = Stabilizer::m_;
	  double g = Stabilizer::g_;
	  double zeta = Stabilizer::zeta_;

	  const Vector& state = stateSIN_.accessCopy ();
	  const Vector& control = controlSIN_.accessCopy ();
	  const double& kth = stiffnessSIN_.accessCopy ();

	  double xi = state (0);
	  double th = state (1);
	  double dxi = state (2);
	  double dth = state (3);

	  double u = control (0);

	  double d2 = (xi*xi+zeta*zeta);
	  double d4 = d2*d2;

	  J.resize (4, 4);
	  J.fill (0.);
	  J (0, 0) = 1.;
	  J (0, 2) = dt_;
	  J (1, 1) = 1.;
	  J (1, 3) = dt_;
	  J (2, 2) = 1.;
	  J (3, 0) = dt_*((-g*cos (th) -2*dth*dxi)/d2
			  +(2*xi*(kth*th + m*g*(cos (th)*xi-sin (th)*zeta)
				  -m*(zeta*u-2*dth*xi*dxi)))
			  /(m*d4));
	  J (3, 1) = dt_*(-kth + m*g*(sin (th)*xi + cos (th)*zeta))/
	    (m*d2);
	  J (3, 2) = -2*dt_*dth*xi/d2;
	  J (3, 3) = 1. - 2.*dt_*xi*dxi/d2;

	  return J;
	}
      };  // class f

      //
      // Observation function
      //
      class h : public Function
      {
	SignalTimeDependent <Vector, int> observationSOUT_;
	SignalTimeDependent <Matrix, int> jacobianSOUT_;

	DYNAMIC_GRAPH_ENTITY_DECL();
	h (const std::string& name) :
	  Function (name),
	  observationSOUT_ ("flexibility_h(" + name +
			    ")::output(vector)::observation"),
	  jacobianSOUT_ ("flexibility_H(" + name +
			 ")::output(vector)::jacobian")
	{
	  signalRegistration (observationSOUT_ << jacobianSOUT_);
	  observationSOUT_.setFunction
	    (boost::bind (&h::computeObservation, this, _1, _2));
	  observationSOUT_.addDependency (stateSIN_);
	  jacobianSOUT_.setFunction (boost::bind (&h::computeJacobian, this, _1, _2));
	  jacobianSOUT_.addDependency (stateSIN_);
	}

	Vector& computeObservation (Vector& obs, const int& time)
	{
	  const Vector& state = stateSIN_.access (time);
	  const double& kth = stiffnessSIN_.accessCopy ();

	  double xi = state (0);
	  double th = state (1);

	  obs.resize (2);
	  obs (0) = xi;
	  obs (1) = kth * th;

	  return obs;
	}
	Matrix& computeJacobian (Matrix& J, const int& time)
	{
	  const Vector& state = stateSIN_.access (time);
	  const double& kth = stiffnessSIN_.accessCopy ();

	  double th = state (1);

	  J.resize (2, 4);
	  J.fill (0.);
	  J (0, 0) = 1.;
	  J (1, 1) = kth;

	  return J;
	}
      }; // class h

      DYNAMICGRAPH_FACTORY_ENTITY_PLUGIN (f, "flexibility_f");
      DYNAMICGRAPH_FACTORY_ENTITY_PLUGIN (h, "flexibility_h");

    } // namespace flexibility

    class MatrixHomoToYawOrientation
      :public Entity
    {
    public: /* --- CONSTRUCTION --- */

      static const std::string CLASS_NAME;

      virtual const std::string& getClassName  () const
      {
	return CLASS_NAME;
      }

      std::string getDocString () const
      {
	return
	  "Extract yaw orientation of matrix homogeneous as rotation matrix\n"
	  "\n"
	  "  Export also inverse of rotation matrix in signal inverse\n";
      }

      MatrixHomoToYawOrientation( const std::string& name ) :
	Entity(name)
	,positionSIN_ (NULL, CLASS_NAME+"(" + name +
		       ")::input(MatrixHomogeneous)::position")
	,yawOrientationSOUT_ (boost::bind (&MatrixHomoToYawOrientation::
					   computeOrientation,this,_1,_2),
			      positionSIN_, CLASS_NAME+"(" + name +
			      ")::output(MatrixRotation)::yawOrientation")
	,inverseSOUT_ (boost::bind (&MatrixHomoToYawOrientation::
				    computeInverse,this, _1, _2),
		       yawOrientationSOUT_,
		       CLASS_NAME+"(" + name +
		       ")::output(MatrixRotation)::inverse")
      {
	signalRegistration (positionSIN_<<yawOrientationSOUT_
			    <<inverseSOUT_);
      }

      virtual ~MatrixHomoToYawOrientation( void ) {};

    public: /* --- SIGNAL --- */

      SignalPtr<MatrixHomogeneous,int> positionSIN_;
      SignalTimeDependent<MatrixRotation,int> yawOrientationSOUT_;
      SignalTimeDependent<MatrixRotation,int> inverseSOUT_;

    protected:
      MatrixRotation& computeOrientation (MatrixRotation& res,int time)
      {
	const MatrixHomogeneous& m = positionSIN_(time);
	const double& nx = m (0,0);
	const double& ny = m (1,0);
	double norm = sqrt (nx*nx + ny*ny);
	double cth = nx/norm;
	double sth = ny/norm;
	res.setIdentity ();
	res (0,0) = cth;
	res (1,0) = sth;
	res (0,1) = -sth;
	res (1,1) = cth;

	return res;
      }

      MatrixRotation& computeInverse (MatrixRotation& res,int time)
      {
	const MatrixRotation& r = yawOrientationSOUT_ (time);
	res.setIdentity ();
	res (0,0) = r (0,0);
	res (1,0) = r (0,1);
	res (0,1) = r (1,0);
	res (1,1) = r (1,1);

	return res;
      }
    };
    DYNAMICGRAPH_FACTORY_ENTITY_PLUGIN (MatrixHomoToYawOrientation,
					"MatrixHomoToYawOrientation");
  } // namespace dynamic
} // namespace sot
