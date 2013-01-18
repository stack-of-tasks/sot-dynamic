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

#ifndef SOT_DYNAMIC_FLEXIBILITY_KALMAN_HH
# define SOT_DYNAMIC_FLEXIBILITY_KALMAN_HH

#include <dynamic-graph/command-setter.h>
#include <dynamic-graph/command-getter.h>
#include <dynamic-graph/command-direct-setter.h>
#include <dynamic-graph/command-direct-getter.h>
#include <dynamic-graph/command-bind.h>

# include "stabilizer.hh"

namespace sot {
  namespace dynamic {
    using dynamicgraph::command::makeCommandVoid0;
    using dynamicgraph::command::docCommandVoid0;
    using dynamicgraph::command::docDirectSetter;
    using dynamicgraph::command::makeDirectSetter;
    using dynamicgraph::command::docDirectGetter;
    using dynamicgraph::command::makeDirectGetter;

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
	SignalPtr < Vector, int > stateSIN_;
	double dt_;

      public:
	Function (const std::string& name) :
	  Entity (name),
	  stateSIN_ (0, "flexibility::Function(" + name +
		     ")::input(vector)::state"),
	  dt_ (.005)
	{
	  signalRegistration (stateSIN_);
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
      // State transition of time discretized system (angular flexibility)
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
	// Cosine of the angle between the line linking both support and
	// the axis of the flexibility measured by Extended Kalman filter
	SignalPtr <double, int> cosineFootSIN_;
	SignalPtr <unsigned int, int> nbSupportSIN_;
	Signal <Vector, int> newStateSOUT_;
	Signal <Matrix, int> jacobianSOUT_;

	DYNAMIC_GRAPH_ENTITY_DECL();
	f (const std::string& name) :
	  Function (name),
	  controlSIN_ (0, "flexibility_f(" + name +
		       ")::input(vector)::control"),
	  cosineFootSIN_ (0, "flexibility_f(" + name +
			  ")::input(double)::cosineFoot"),
	  nbSupportSIN_ (0, "flexibility_f(" + name +
			 ")::input(double)::nbSupport"),
	  newStateSOUT_ ("flexibility_f(" + name +
			 ")::output(vector)::newState"),
	  jacobianSOUT_ ("flexibility_f(" + name +
			 ")::output(vector)::jacobian")

	{
	  signalRegistration (controlSIN_ << cosineFootSIN_ << newStateSOUT_
			      << jacobianSOUT_ << nbSupportSIN_);
	  newStateSOUT_.setFunction (boost::bind (&f::computeNewState, this,
						  _1, _2));
	  newStateSOUT_.addDependency (controlSIN_);
	  newStateSOUT_.addDependency (cosineFootSIN_);
	  newStateSOUT_.addDependency (nbSupportSIN_);

	  jacobianSOUT_.setFunction (boost::bind (&f::computeJacobian, this,
						  _1, _2));
	  jacobianSOUT_.addDependency (controlSIN_);
	  jacobianSOUT_.addDependency (cosineFootSIN_);
	  jacobianSOUT_.addDependency (nbSupportSIN_);
	}

	Vector& computeNewState (Vector& x, const int& time)
	{
	  double m = Stabilizer::m_;
	  double g = Stabilizer::g_;
	  double zeta = Stabilizer::zeta_;

	  const Vector& state = stateSIN_.accessCopy ();
	  const Vector& control = controlSIN_.access (time);
	  const double& cosineFoot = cosineFootSIN_.access (time);
	  const unsigned int& nbSupport = nbSupportSIN_.access (time);

	  double xi = state (0);
	  double th = state (1);
	  double dxi = state (2);
	  double dth = state (3);
	  double kth = state (4);

	  double u = control (0);
	  double d2 = (xi*xi+zeta*zeta);
	  x.resize (5);

	  x (0) = xi + dt_ * dxi;
	  x (1) = th + dt_ * dth;
	  x (2) = dxi + dt_ * u;
	  switch (nbSupport) {
	  case 0:
	    x (3) = 0;
	    break;
	  case 1:
	    x (3) = dth + dt_* (-kth*th - m*g*(cos (th)*xi - sin (th)*zeta) +
				m*(zeta*u -2*dth*xi*dxi))/(m*d2);
	    break;
	  case 2:
	    x (3) = dth + dt_* cosineFoot *
	      (-2*kth*th - m*g*(cos (th)*xi - sin (th)*zeta) +
	       m*(zeta*u -2*dth*xi*dxi))/(m*d2);
	  default:
	    break;
	  }

	  x (4) = kth;

	  return x;
	}

	Matrix& computeJacobian (Matrix& J, const int& t)
	{
	  double m = Stabilizer::m_;
	  double g = Stabilizer::g_;
	  double zeta = Stabilizer::zeta_;

	  const Vector& state = stateSIN_.accessCopy ();
	  const Vector& control = controlSIN_.access (t);
	  const double& cosineFoot = cosineFootSIN_.access (t);
	  const unsigned int& nbSupport = nbSupportSIN_.access (t);

	  double xi = state (0);
	  double th = state (1);
	  double dxi = state (2);
	  double dth = state (3);
	  double kth = state (4);

	  double u = control (0);

	  double d2 = (xi*xi+zeta*zeta);
	  double d4 = d2*d2;

	  J.resize (5, 5);
	  J.fill (0.);
	  J (0, 0) = 1.;
	  J (0, 2) = dt_;
	  J (1, 1) = 1.;
	  J (1, 3) = dt_;
	  J (2, 2) = 1.;
	  switch (nbSupport) {
	  case 0:
	    J (3, 3) = 1.;
	    break;
	  case 1:
	    J (3, 0) = dt_*((-g*cos (th) -2*dth*dxi)/d2
			    +(2*xi*(kth*th + m*g*(cos (th)*xi-sin (th)*zeta)
				    -m*(zeta*u-2*dth*xi*dxi)))
			    /(m*d4));
	    J (3, 1) = dt_*(-kth + m*g*(sin (th)*xi + cos (th)*zeta))/
	      (m*d2);
	    J (3, 2) = -2*dt_*dth*xi/d2;
	    J (3, 3) = 1. - 2.*dt_*xi*dxi/d2;
	    J (3, 4) = -dt_*th/(m*d2);
	    break;
	  case 2:
	    J (3, 0) = dt_*cosineFoot*
	      ((-g*cos (th) -2*dth*dxi)/d2
	       +(2*xi*(kth*th + m*g*(cos (th)*xi-sin (th)*zeta)
		       -m*(zeta*u-2*dth*xi*dxi)))/(m*d4));
	    J (3, 1) = dt_*cosineFoot*
	      (-2*kth + m*g*(sin (th)*xi + cos (th)*zeta))/(m*d2);
	    J (3, 2) = -2*dt_*cosineFoot*dth*xi/d2;
	    J (3, 3) = 1. - 2.*dt_*cosineFoot*xi*dxi/d2;
	    J (3, 4) = -2*dt_*cosineFoot*th/(m*d2);
	  default:
	    break;
	  }
	  J (4, 4) = 1.;

	  return J;
	}
	std::string getDocString () const
	{
	  return
	    "State transition for foot flexibility along one local coordinate "
	    "axis\n"
	    "\n"
	    "  Compute expected state at time k from state and control at time"
	    " k-1.\n"
	    "  State is defined by center of mass deviation (wrt reference), "
	    "flexibility\n"
	    "  angular deviation and derivatives of these two values, along a "
	    "local\n"
	    "  axis of one foot.\n"
	    "  \n"
	    "  Signal \"cosineFoot\" contains the cosine of angle between the "
	    "line linking\n"
	    "  both ankles (in horizontal plane) in double support and the "
	    "local coordinate\n"
	    "  axis of the foot.\n";
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
	  jacobianSOUT_ ("flexibility_h(" + name +
			 ")::output(vector)::jacobian")
	{
	  signalRegistration (observationSOUT_ << jacobianSOUT_);
	  observationSOUT_.setFunction
	    (boost::bind (&h::computeObservation, this, _1, _2));
	  observationSOUT_.addDependency (stateSIN_);
	  jacobianSOUT_.setFunction (boost::bind (&h::computeJacobian,
						  this, _1, _2));
	  jacobianSOUT_.addDependency (stateSIN_);
	}

	Vector& computeObservation (Vector& obs, const int& time)
	{
	  const Vector& state = stateSIN_.access (time);

	  double xi = state (0);
	  double th = state (1);
	  double kth = state (4);

	  obs.resize (2);
	  obs (0) = xi;
	  obs (1) = kth * th;

	  return obs;
	}
	Matrix& computeJacobian (Matrix& J, const int& time)
	{
	  const Vector& state = stateSIN_.access (time);

	  double th = state (1);
	  double kth = state (4);

	  J.resize (2, 5);
	  J.fill (0.);
	  J (0, 0) = 1.;
	  J (1, 1) = kth;
	  J (1, 4) = th;

	  return J;
	}
	std::string getDocString () const
	{
	  return
	    "Observation function for foot flexibility along one local "
	    "coordinate axis\n"
	    "\n"
	    "  Compute expected observation at time k from state at time k.\n"
	    "  State is defined by center of mass deviation (wrt reference), "
	    "flexibility\n"
	    "  angular deviation and derivatives of these two values, along a "
	    "local\n"
	    "  axis of one foot.\n"
	    "  Observation is defined by center of mass deviation and moment "
	    "at the ankle\n";
	}
      }; // class h

      DYNAMICGRAPH_FACTORY_ENTITY_PLUGIN (f, "flexibility_f");
      DYNAMICGRAPH_FACTORY_ENTITY_PLUGIN (h, "flexibility_h");

      //
      // State transition of time discretized system (linear vertical
      //                                              flexibility)
      //
      //   x    =  f (x , u )
      //    k+1        k   k
      //
      //   with
      //        ..
      //   u  = xi,   x  = x (k dt)
      //    k          k
      //
      class fz : public Function
      {
	SignalPtr <Vector, int> controlSIN_;
	Signal <Vector, int> newStateSOUT_;
	Signal <Matrix, int> jacobianSOUT_;

	DYNAMIC_GRAPH_ENTITY_DECL();
	fz (const std::string& name) :
	  Function (name),
	  controlSIN_ (0, "flexibility_fz(" + name +
		       ")::input(vector)::control"),
	  newStateSOUT_ ("flexibility_fz(" + name +
			 ")::output(vector)::newState"),
	  jacobianSOUT_ ("flexibility_fz(" + name +
			 ")::output(vector)::jacobian")

	{
	  signalRegistration (controlSIN_ << newStateSOUT_ << jacobianSOUT_);
	  newStateSOUT_.setFunction (boost::bind (&fz::computeNewState, this,
						  _1, _2));
	  newStateSOUT_.addDependency (controlSIN_);
	  jacobianSOUT_.setFunction (boost::bind (&fz::computeJacobian, this,
						  _1, _2));
	  jacobianSOUT_.addDependency (controlSIN_);
	}

	Vector& computeNewState (Vector& x, const int& time)
	{
	  double m = Stabilizer::m_;

	  const Vector& state = stateSIN_.accessCopy ();
	  const Vector& control = controlSIN_.access (time);

	  double zeta = state (0);
	  double th = state (1);
	  double dzeta = state (2);
	  double dth = state (3);
	  double kz = state (4);

	  double u = control (0);
	  x.resize (5);

	  x (0) = zeta + dt_ * dzeta;
	  x (1) = th + dt_ * dth;
	  x (2) = dzeta + dt_ * u;
	  x (3) = dth + dt_ * (-kz/m*th - u);
	  x (4) = kz;

	  return x;
	}

	Matrix& computeJacobian (Matrix& J, const int&)
	{
	  double m = Stabilizer::m_;

	  const Vector& state = stateSIN_.accessCopy ();

	  double kz = state (4);

	  J.resize (5, 5);
	  J.setIdentity ();
	  J (0, 2) = dt_;
	  J (1, 3) = dt_;
	  J (3, 1) = dt_ * (-kz/m);

	  return J;
	}
	std::string getDocString () const
	{
	  return
	    "State transition for foot flexibility along z axis\n"
	    "\n"
	    "  Compute expected state at time k from state and control at time"
	    " k-1.\n"
	    "  State is defined by center of mass deviation (wrt reference), "
	    "flexibility\n"
	    "  linear deviation and derivatives of these two values.\n";
	}
      };  // class f

      //
      // Observation function
      //
      class hz : public Function
      {
	SignalTimeDependent <Vector, int> observationSOUT_;
	SignalTimeDependent <Matrix, int> jacobianSOUT_;

	DYNAMIC_GRAPH_ENTITY_DECL();
	hz (const std::string& name) :
	  Function (name),
	  observationSOUT_ ("flexibility_hz(" + name +
			    ")::output(vector)::observation"),
	  jacobianSOUT_ ("flexibility_hz(" + name +
			 ")::output(vector)::jacobian")
	{
	  signalRegistration (observationSOUT_ << jacobianSOUT_);
	  observationSOUT_.setFunction
	    (boost::bind (&hz::computeObservation, this, _1, _2));
	  observationSOUT_.addDependency (stateSIN_);
	  jacobianSOUT_.setFunction (boost::bind (&hz::computeJacobian,
						  this, _1, _2));
	  jacobianSOUT_.addDependency (stateSIN_);
	}

	Vector& computeObservation (Vector& obs, const int& time)
	{
	  const Vector& state = stateSIN_.access (time);

	  double zeta = state (0);
	  double th = state (1);
	  double kz = state (4);

	  obs.resize (2);
	  obs (0) = zeta;
	  obs (1) = - kz * th;

	  return obs;
	}
	Matrix& computeJacobian (Matrix& J, const int& time)
	{
	  const Vector& state = stateSIN_.access (time);

	  double th = state (1);
	  double kz = state (4);

	  J.resize (2, 5);
	  J.setZero ();
	  J (0, 0) = 1.;
	  J (1, 1) = -kz;
	  J (1, 4) = -th;

	  return J;
	}
	std::string getDocString () const
	{
	  return
	    "Observation function for foot flexibility along z axis\n"
	    "\n"
	    "  Compute expected observation at time k from state at time k.\n"
	    "  State is defined by center of mass deviation (wrt reference), "
	    "flexibility\n"
	    "  angular deviation and derivatives of these two values\n"
	    "  Observation is defined by center of mass deviation and vertical\n"
	    "  component of the forces in the feet\n";
	}
      }; // class hz

      DYNAMICGRAPH_FACTORY_ENTITY_PLUGIN (fz, "flexibility_fz");
      DYNAMICGRAPH_FACTORY_ENTITY_PLUGIN (hz, "flexibility_hz");

    } // namespace flexibility

    class VarianceDoubleSupport
      :public Entity
    {
      SignalPtr <Matrix, int> varianceSIN_;
      SignalPtr <unsigned int, int> nbSupportSIN_;
      SignalPtr <double, int> cosineFootSIN_;
      SignalPtr <double, int> sigmaSIN_;
      SignalTimeDependent <Matrix, int> varianceSOUT_;

      DYNAMIC_GRAPH_ENTITY_DECL();
      virtual ~VarianceDoubleSupport () {};

      std::string getDocString () const
      {
	return
	  "Add variance on flexibility evolution in double support\n";
      }

      VarianceDoubleSupport (const std::string& name) :
	Entity (name),
	varianceSIN_ (0, "VarianceDoubleSupport(" + name +
	 	      ")::input(Matrix)::varianceIn"),
	nbSupportSIN_ (0, "VarianceDoubleSupport(" + name +
	 	       ")::input(unsigned int)::nbSupport"),
	cosineFootSIN_ (0, "VarianceDoubleSupport(" + name +
			")::input(double)::cosineFoot"),
	sigmaSIN_ (0, "VarianceDoubleSupport(" + name +
	 	   ")::input(double)::sigma"),
	varianceSOUT_ ("VarianceDoubleSupport(" + name +
	 	       ")::input(double)::varianceOut")
      {
	signalRegistration (varianceSIN_  << nbSupportSIN_ << cosineFootSIN_
			    << sigmaSIN_ << varianceSOUT_ );
	varianceSOUT_.addDependency (varianceSIN_);
	varianceSOUT_.addDependency (nbSupportSIN_);
	varianceSOUT_.addDependency (cosineFootSIN_);
	varianceSOUT_.addDependency (sigmaSIN_);

	varianceSOUT_.setFunction (boost::bind
	 			   (&VarianceDoubleSupport::computeVariance,
	 			    this, _1, _2));
      }

      Matrix& computeVariance (Matrix& sout, const int& t)
      {
       	const unsigned int& nbSupport = nbSupportSIN_.access (t);
       	const Matrix& varianceIn = varianceSIN_.access (t);
       	const double& sigma = sigmaSIN_.access (t);
	const double& cosine = cosineFootSIN_.access (t);

 	sout = varianceIn;
      	if (nbSupport == 2) {
       	  sout (3, 3) += sigma * (1 - cosine*cosine);
       	}
       	return sout;
      }
    }; // class VarianceDoubleSupport

    DYNAMICGRAPH_FACTORY_ENTITY_PLUGIN (VarianceDoubleSupport,
    					"VarianceDoubleSupport");
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

#endif // SOT_DYNAMIC_FLEXIBILITY_KALMAN_HH
