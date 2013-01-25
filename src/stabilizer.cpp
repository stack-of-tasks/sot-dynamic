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

#include <dynamic-graph/command-setter.h>
#include <dynamic-graph/command-getter.h>
#include <dynamic-graph/command-direct-setter.h>
#include <dynamic-graph/command-direct-getter.h>
#include <dynamic-graph/command-bind.h>

#include "stabilizer.hh"

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
    using dynamicgraph::sot::VectorUTheta;

    double Stabilizer::m_ = 56.;
    double Stabilizer::g_ = 9.81;
    double Stabilizer::zeta_ = .80;

    DYNAMICGRAPH_FACTORY_ENTITY_PLUGIN (Stabilizer, "Stabilizer");

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
    Stabilizer::Stabilizer(const std::string& inName) :
      TaskAbstract(inName),
      deltaComSIN_ (NULL, "Stabilizer("+inName+")::input(vector)::deltaCom"),
      jacobianSIN_ (NULL, "Stabilizer("+inName+")::input(matrix)::Jcom"),
      comdotSIN_ (NULL, "Stabilizer("+inName+")::input(vector)::comdot"),
      leftFootPositionSIN_
      (NULL, "Stabilizer("+inName+")::input(vector)::leftFootPosition"),
      rightFootPositionSIN_
      (NULL, "Stabilizer("+inName+")::input(vector)::rightFootPosition"),
      forceLeftFootSIN_
      (NULL, "Stabilizer("+inName+")::input(vector)::force_lf"),
      forceRightFootSIN_
      (NULL, "Stabilizer("+inName+")::input(vector)::force_rf"),
      forceLeftFootRefSIN_
      (NULL, "Stabilizer("+inName+")::input(vector)::forceRef_lf"),
      forceRightFootRefSIN_
      (NULL, "Stabilizer("+inName+")::input(vector)::forceRef_rf"),
      stateFlexRfxSIN_
      (NULL, "Stabilizer("+inName+")::input(vector)::stateFlex_rfx"),
      stateFlexRfySIN_
      (NULL, "Stabilizer("+inName+")::input(vector)::stateFlex_rfy"),
      stateFlexLfxSIN_
      (NULL, "Stabilizer("+inName+")::input(vector)::stateFlex_lfx"),
      stateFlexLfySIN_
      (NULL, "Stabilizer("+inName+")::input(vector)::stateFlex_lfy"),
      stateFlexZSIN_
      (NULL, "Stabilizer("+inName+")::input(vector)::stateFlex_z"),
      stateFlexLatSIN_
      (0, "Stabilizer("+inName+")::input(vector)::stateFlex_lat"),
      controlGainSIN_
      (NULL, "Stabilizer("+inName+")::input(double)::controlGain"),
      d2comSOUT_ ("Stabilizer("+inName+")::output(vector)::d2com"),
      cosineRightFootXSOUT_
      ("Stabilizer("+inName+")::output(double)::cosine_rfx"),
      cosineRightFootYSOUT_
      ("Stabilizer("+inName+")::output(double)::cosine_rfy"),
      cosineLeftFootXSOUT_
      ("Stabilizer("+inName+")::output(double)::cosine_lfx"),
      cosineLeftFootYSOUT_
      ("Stabilizer("+inName+")::output(double)::cosine_lfy"),
      nbSupportSOUT_
      ("Stabilizer("+inName+")::output(unsigned int)::nbSupport"),
      flexPositionSOUT_
      ("Stabilizer("+inName+")::output(MatrixHomo)::flexPosition"),
      flexPositionLfSOUT_
      ("Stabilizer("+inName+")::output(MatrixHomo)::flexPosition_lf"),
      flexPositionRfSOUT_
      ("Stabilizer("+inName+")::output(MatrixHomo)::flexPosition_rf"),
      flexPositionLatSOUT_
      ("Stabilizer("+inName+")::output(MatrixHomo)::flexPosition_lat"),
      flexVelocitySOUT_
      ("Stabilizer("+inName+")::output(MatrixHomo)::flexVelocity"),
      flexVelocityLfSOUT_
      ("Stabilizer("+inName+")::output(MatrixHomo)::flexVelocity_lf"),
      flexVelocityRfSOUT_
      ("Stabilizer("+inName+")::output(MatrixHomo)::flexVelocity_rf"),
      flexVelocityLatSOUT_
      ("Stabilizer("+inName+")::output(MatrixHomo)::flexVelocity_lat"),
      flexZobsSOUT_ ("Stabilizer("+inName+")::output(Vector)::flexZobs"),
      stepLengthSOUT_ ("Stabilizer("+inName+")::output(double)::stepLength"),
      flexLatControlSOUT_
      ("Stabilizer("+inName+")::output(Vector)::flexLatControl"),
      flexLatObsSOUT_
      ("Stabilizer("+inName+")::output(Vector)::flexLatObs"),
      debugSOUT_ ("Stabilizer("+inName+")::output(vector)::debug"),
      gain1_ (4), gain2_ (4), gainz_ (4), gainLat_ (4),
      prevCom_(3), flexValue_ (3), flexDeriv_ (3),
      dcom_ (3), dt_ (.005), on_ (false),
      forceThreshold_ (.036*m_*g_), angularStiffness_ (425.), d2com_ (3),
      deltaCom_ (3),
      cosineLeftFootX_ (0.), cosineLeftFootY_ (0.),
      cosineRightFootX_ (0.), cosineRightFootY_ (0.),
      flexPosition_ (), flexPositionLf_ (), flexPositionRf_ (),
      flexPositionLat_ (),
      flexVelocity_ (6), flexVelocityLf_ (6), flexVelocityRf_ (6),
      flexVelocityLat_ (6),
      flexZobs_ (2), flexLatControl_ (1), flexLatObs_ (2),
      timeBeforeFlyingFootCorrection_ (.1),
      iterationsSinceLastSupportLf_ (0), iterationsSinceLastSupportRf_ (0),
      supportCandidateLf_ (0), supportCandidateRf_ (0),
      uth_ (),
      R_ (), translation_ (3), zmp_ (3), debug_ (4)
    {
      // Register signals into the entity.
      signalRegistration (deltaComSIN_);
      signalRegistration (jacobianSIN_);
      signalRegistration (comdotSIN_);
      signalRegistration (leftFootPositionSIN_ << rightFootPositionSIN_
			  << forceRightFootSIN_ << forceLeftFootSIN_
			  << forceLeftFootRefSIN_ << forceRightFootRefSIN_);
      signalRegistration (stateFlexRfxSIN_ << stateFlexRfySIN_
			  << stateFlexLfxSIN_ << stateFlexLfySIN_
			  << stateFlexZSIN_ << stateFlexLatSIN_);
      signalRegistration (controlGainSIN_);
      signalRegistration (d2comSOUT_);
      signalRegistration (cosineRightFootXSOUT_ << cosineRightFootYSOUT_
			  << cosineLeftFootXSOUT_ << cosineLeftFootYSOUT_);
      signalRegistration (nbSupportSOUT_);
      signalRegistration (flexPositionSOUT_ << flexPositionLfSOUT_
			  << flexPositionRfSOUT_ << flexPositionLatSOUT_
			  << flexZobsSOUT_);
      signalRegistration (flexVelocitySOUT_ << flexVelocityLfSOUT_
			  << flexVelocityRfSOUT_ << flexVelocitySOUT_);
      signalRegistration (stepLengthSOUT_ << flexLatControlSOUT_
			  << flexLatObsSOUT_);
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

      jacobianSOUT.addDependency (jacobianSIN_);

      taskSOUT.setFunction (boost::bind(&Stabilizer::computeControlFeedback,
					this,_1,_2));
      jacobianSOUT.setFunction (boost::bind(&Stabilizer::computeJacobianCom,
					    this,_1,_2));
      cosineRightFootXSOUT_.addDependency (taskSOUT);
      cosineRightFootYSOUT_.addDependency (taskSOUT);
      cosineLeftFootXSOUT_.addDependency (taskSOUT);
      cosineLeftFootYSOUT_.addDependency (taskSOUT);
      nbSupportSOUT_.addDependency (taskSOUT);

      d2com_.setZero ();
      dcom_.setZero ();
      deltaCom_.setZero ();
      d2comSOUT_.setConstant (d2com_);
      debug_.setZero ();
      debugSOUT_.setConstant (debug_);
      flexVelocity_.setZero ();
      flexVelocityLat_.setZero ();

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

      addCommand ("setGainz",
		  makeDirectSetter (*this, &gainz_,
				    docDirectSetter
				    ("Set gains of vertical flexibility",
				     "vector")));

      addCommand ("getGainz",
		  makeDirectGetter (*this, &gainz_,
				    docDirectGetter
				    ("Get gains of vertical flexibility",
				     "vector")));

      addCommand ("setGainLateral",
		  makeDirectSetter (*this, &gainLat_,
				    docDirectSetter
				    ("Set gains of lateral flexibility",
				     "vector")));

      addCommand ("getGainLateral",
		  makeDirectGetter (*this, &gainLat_,
				    docDirectGetter
				    ("Get gains of lateral flexibility",
				     "vector")));

      prevCom_.fill (0.);
      flexValue_.fill (0.);
      flexDeriv_.fill (0.);
      flexZobs_.setZero ();
      flexZobsSOUT_.setConstant (flexZobs_);
      stepLengthSOUT_.setConstant (0.);

      flexLatControl_.setZero ();
      flexLatControlSOUT_.setConstant (flexLatControl_);
      flexLatObs_.setZero ();
      flexLatObsSOUT_.setConstant (flexLatObs_);

      // Single support gains for
      //  - kth = 510,
      //  - zeta = .8
      //  - m = 56
      //  - eigen values = 3, 3, 6, 6.
      gain1_ (0) = 121.27893435294121;
      gain1_ (1) = -19.899754625210093;
      gain1_ (2) = 38.133760000000009;
      gain1_ (3) = -17.707008000000005;

      // Double support gains for
      //  - kth = 2*510,
      //  - zeta = .8
      //  - m = 56
      //  - eigen values = 4, 4, 8, 8.
      gain2_ (0) = 118.62268196078432;
      gain2_ (1) = 58.543997288515406;
      gain2_ (2) = 37.326305882352941;
      gain2_ (3) = -10.661044705882359;

      // Vectical gains for
      //  - kz = 150000
      //  - m = 56
      //  - eigen values = 21, 21, 21, 21.
      gainz_ (0) = 25.23;
      gainz_ (1) = -124.77;
      gainz_ (2) = 10.19;
      gainz_ (3) = -9.81;

      // Lateral gains for
      //  - kth = 3000 (kz = 160000, h=.19)
      //  - m = 56
      //  - eigen values = 8, 8, 8, 8.
      gainLat_ (0) = 94.721917866666672;
      gainLat_ (1) = 174.26817999238096;
      gainLat_ (2) = 29.154645333333335;
      gainLat_ (3) = 2.2762837333333308;

      zmp_.setZero ();
    }

    /// Compute flexibility state from both feet
    void Stabilizer::computeFlexibility (const int& time)
    {
      const Vector& flexRfx = stateFlexRfxSIN_.access (time);
      const Vector& flexRfy = stateFlexRfySIN_.access (time);
      const Vector& flexLfx = stateFlexLfxSIN_.access (time);
      const Vector& flexLfy = stateFlexLfySIN_.access (time);
      const Vector& flexZ = stateFlexZSIN_.access (time);
      const Vector& flexLat = stateFlexLatSIN_.access (time);
      const MatrixHomogeneous& Mr = rightFootPositionSIN_.access (time);
      const MatrixHomogeneous& Ml = leftFootPositionSIN_.access (time);
      const Vector& fr = forceRightFootSIN_.access (time);
      const Vector& fl = forceLeftFootSIN_.access (time);
      const Vector& forceRefLf = forceLeftFootRefSIN_.access (time);
      const Vector& forceRefRf = forceRightFootRefSIN_.access (time);
      double deltaComRfx, deltaComRfy, deltaComLfx, deltaComLfy;
      double dcomRfx, dcomRfy, dcomLfx, dcomLfy;

      // compute component of angle orthogonal to the line joining the feet
      double delta_x = Ml (0, 3) - Mr (0, 3);
      double delta_y = Ml (1, 3) - Mr (1, 3);
      double stepLength = sqrt (delta_x*delta_x+delta_y*delta_y);
      stepLengthSOUT_.setConstant (stepLength);

      u2x_ = delta_x/stepLength;
      u2y_ = delta_y/stepLength;
      u1x_ = u2y_;
      u1y_ = -u2x_;

      // Express vertical component of force in global basis
      double flz = Ml (2,0) * fl (0) + Ml(2,1) * fl (1) + Ml (2,2) * fl (2);
      double frz = Mr (2,0) * fr (0) + Mr(2,1) * fr (1) + Mr (2,2) * fr (2);

      kth_ = flexRfx (4) + flexLfx (4);
      nbSupport_ = 0;
      if (on_) {
	if (frz >= forceThreshold_) {
	  nbSupport_++;
	  supportCandidateRf_++;
	  if (supportCandidateRf_ >= 3) {
	    iterationsSinceLastSupportRf_ = 0;
	  }
	} else {
	  supportCandidateRf_ = 0;
	  iterationsSinceLastSupportRf_ ++;
	}
	if (flz >= forceThreshold_) {
	  nbSupport_++;
	  supportCandidateLf_++;
	  if (supportCandidateLf_ >= 3) {
	    iterationsSinceLastSupportLf_ = 0;
	  }
	} else {
	  supportCandidateLf_ = 0;
	  iterationsSinceLastSupportLf_++;
	}
      }
      if (nbSupport_ == 2) {
	// Compute reference moment from reference forces
	double kthLat = .5*stepLength*stepLength*flexLat (4);
	double Mu1Ref;
	Mu1Ref = .5*(forceRefRf (2) - forceRefLf (2))*stepLength;
	theta1Ref_ = Mu1Ref/kthLat;
      }

      if (frz < 0) frz = 0;
      if (flz < 0) flz = 0;
      double Fz = flz + frz;
      flexZobs_ (1) = Fz - m_ * g_;
      flexLatObs_ (1) = .5*stepLength*(frz - flz);
      if (Fz == 0) {
	flexValue_ (0) = 0;
	flexValue_ (1) = 0;
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
      // Flexibility right foot in global frame
      double flexAngleRfx = cth * flexRfx (1) + sth * flexRfy (1);
      double flexAngleRfy = -sth * flexRfx (1) + cth * flexRfy (1);
      double flexDerivRfx = cth * flexRfx (3) + sth * flexRfy (3);
      double flexDerivRfy = -sth * flexRfx (3) + cth * flexRfy (3);
      // Compute deviation of center of mass
      deltaComRfx = cth * flexRfx (0) - sth * flexRfy (0);
      deltaComRfy = sth * flexRfx (0) + cth * flexRfy (0);
      dcomRfx = cth * flexRfx (2) - sth * flexRfy (2);
      dcomRfy = sth * flexRfx (2) + cth * flexRfy (2);
      // Extract yaw from left foot position
      nx = Ml (0,0);
      ny = Ml (1,0);
      norm = sqrt (nx*nx + ny*ny);
      cth = nx/norm;
      sth = ny/norm;
      // Flexibility left foot in global frame
      double flexAngleLfx = cth * flexLfx (1) + sth * flexLfy (1);
      double flexAngleLfy = -sth * flexLfx (1) + cth * flexLfy (1);
      double flexDerivLfx = cth * flexLfx (3) + sth * flexLfy (3);
      double flexDerivLfy = -sth * flexLfx (3) + cth * flexLfy (3);

      flexValue_ (0) = (frz * flexAngleRfx + flz * flexAngleLfx)/Fz;
      flexValue_ (1) = (frz * flexAngleRfy + flz * flexAngleLfy)/Fz;
      flexValue_ (2) = flexZ (1);
      flexDeriv_ (0) = (frz * flexDerivRfx + flz * flexDerivLfx)/Fz;
      flexDeriv_ (1) = (frz * flexDerivRfy + flz * flexDerivLfy)/Fz;
      flexDeriv_ (2) = flexZ (3);
      // Compute deviation of center of mass
      deltaComLfx = cth * flexLfx (0) - sth * flexLfy (0);
      deltaComLfy = sth * flexLfx (0) + cth * flexLfy (0);
      dcomLfx = cth * flexLfx (2) - sth * flexLfy (2);
      dcomLfy = sth * flexLfx (2) + cth * flexLfy (2);

      deltaCom_ (0) = (frz * deltaComRfx + flz * deltaComLfx)/Fz;
      deltaCom_ (1) = (frz * deltaComRfy + flz * deltaComLfy)/Fz;
      deltaCom_ (2) = flexZ (0);
      dcom_ (0) = (frz * dcomRfx + flz * dcomLfx)/Fz;
      dcom_ (1) = (frz * dcomRfy + flz * dcomLfy)/Fz;
      dcom_ (2) = flexZ (2);
      // Compute flexibility transformations and velocities
      if (Fz != 0) {
	zmp_ (0) = (frz * Mr (0, 3) + flz * Ml (0, 3))/Fz;
	zmp_ (1) = (frz * Mr (1, 3) + flz * Ml (1, 3))/Fz;
	zmp_ (2) = (frz * Mr (2, 3) + flz * Ml (2, 3))/Fz;
	uth_ (0) = flexValue_ (1);
	uth_ (1) = -flexValue_ (0);
	translation_ = zmp_ - R_ * zmp_;
	uth_.toMatrix (R_);
	for (std::size_t row = 0; row < 3; ++row) {
	  for (std::size_t col = 0; col < 3; ++col) {
	    flexPosition_ (row, col) = R_ (row, col);
	  }
	  flexPosition_ (row, 3) = translation_ (row);
	}
	// Lateral flexibility
	double theta = flexLat (1);
	uth_ (0) = u1x_ * theta;
	uth_ (1) = u1y_ * theta;
	uth_.toMatrix (R_);
	translation_ = zmp_ - R_ * zmp_;
	for (std::size_t row = 0; row < 3; ++row) {
	  for (std::size_t col = 0; col < 3; ++col) {
	    flexPositionLat_ (row, col) = R_ (row, col);
	  }
	  flexPositionLat_ (row, 3) = translation_ (row);
	}

	flexVelocity_ (0) = zmp_ (2) * flexDeriv_ (0);
	flexVelocity_ (1) = zmp_ (2) * flexDeriv_ (1);
	flexVelocity_ (2) = -zmp_ (0)*flexDeriv_ (0)-zmp_ (1)*flexDeriv_ (1);
	flexVelocity_ (3) = flexDeriv_ (1);
	flexVelocity_ (4) = -flexDeriv_ (0);

	if (iterationsSinceLastSupportLf_ * dt_ >
	    timeBeforeFlyingFootCorrection_) {
	  flexPositionLf_ = flexPosition_;
	  flexVelocityLf_ = flexVelocity_;
	} else {
	  flexPositionLf_.setIdentity ();
	  flexVelocityLf_.setZero ();
	}
	if (iterationsSinceLastSupportRf_ * dt_ >
	    timeBeforeFlyingFootCorrection_) {
	  flexPositionRf_ = flexPosition_;
	  flexVelocityRf_ = flexVelocity_;
	} else {
	  flexPositionRf_.setIdentity ();
	  flexVelocityRf_.setZero ();
	}
	flexPositionSOUT_.setConstant (flexPosition_);
	flexPositionLfSOUT_.setConstant (flexPositionLf_);
	flexPositionRfSOUT_.setConstant (flexPositionRf_);
	flexPositionLatSOUT_.setConstant (flexPositionLat_);
	flexVelocitySOUT_.setConstant (flexVelocity_);
	flexVelocityLfSOUT_.setConstant (flexVelocityLf_);
	flexVelocityRfSOUT_.setConstant (flexVelocityRf_);
	flexVelocityLatSOUT_.setConstant (flexVelocityLat_);
      }
    }

    /// Compute the control law
    VectorMultiBound&
    Stabilizer::computeControlFeedback(VectorMultiBound& comdot,
				       const int& time)
    {
      const Vector& deltaCom = deltaComSIN_ (time);
      const Vector& comdotRef = comdotSIN_ (time);
      const MatrixHomogeneous& leftFootPosition =
	leftFootPositionSIN_.access (time);
      const MatrixHomogeneous& rightFootPosition =
	rightFootPositionSIN_.access (time);
      const double& gain = controlGainSIN_.access (time);
      const Vector& forceLf = forceLeftFootSIN_.access (time);
      const Vector& forceRf = forceRightFootSIN_.access (time);
      const double& stepLength = stepLengthSOUT_.accessCopy ();

      computeFlexibility (time);

      double x = deltaCom_ (0);
      double y = deltaCom_ (1);
      double z = deltaCom_ (2);

      // z-component of center of mass deviation in global frame
      flexZobs_ (0) = deltaCom (2);
      flexZobsSOUT_.setConstant (flexZobs_);
      flexLatObs_ (0) = 0;
      flexLatObsSOUT_.setConstant (flexLatObs_);

      double theta0, dtheta0;
      double theta1, dtheta1, ddxi;
      double xi, dxi, lat, dlat, ddlat;
      double thetaz = flexValue_ (2);
      double dthetaz = flexDeriv_ (2);
      double fzRef, Zrefx, Zrefy, fz, Zx, Zy;

      cosineLeftFootX_ = leftFootPosition (0, 0)*u1x_ +
	leftFootPosition (1, 0)*u1y_;
      cosineLeftFootY_ = leftFootPosition (0, 1)*u1x_ +
	leftFootPosition (1, 1)*u1y_;;
      cosineRightFootX_ = rightFootPosition (0, 0)*u1x_ +
	rightFootPosition (1, 0)*u1y_;
      cosineRightFootY_ = rightFootPosition (0, 1)*u1x_ +
	rightFootPosition (1, 1)*u1y_;;

      switch (nbSupport_) {
      case 0:
	dcom_ (0) = -gain * x;
	dcom_ (1) = -gain * y;
	dcom_ (2) = -gain * z;
	break;
      case 1: //single support
	//along x
	theta0 = flexValue_ (0);
	dtheta0 = flexDeriv_ (0);
	d2com_ (0)= -(gain1_ (0)*x + gain1_ (1)*theta0 +
		      gain1_ (2)*dcom_ (0) + gain1_ (3)*dtheta0);
	dcom_ (0) += dt_ * d2com_ (0);
	// along y
	theta1 = flexValue_ (1);
	dtheta1 = flexDeriv_ (1);
	d2com_ (1) = - (gain1_ (0)*y + gain1_ (1)*theta1 +
			gain1_ (2)*dcom_ (1) + gain1_ (3)*dtheta1);
	dcom_ (1) += dt_ * d2com_ (1);
	// along z
	d2com_ (2) = - (gainz_ (0)*z + gainz_ (1)*thetaz +
			gainz_ (2)*dcom_ (2) + gainz_ (3)*dthetaz);
	dcom_ (2) += dt_ * d2com_ (2);
	break;
      case 2: //double support
	theta0 = u1x_ * flexValue_ (0) + u1y_ * flexValue_ (1);
	dtheta0 = u1x_ * flexDeriv_ (0) + u1y_ * flexDeriv_ (1);
	xi = u1x_*x + u1y_*y;
	dxi = u1x_*dcom_ (0) + u1y_*dcom_ (1);
	ddxi = - (gain2_ (0)*xi + gain2_ (1)*theta0 + gain2_ (2)*dxi +
		  gain2_ (3)*dtheta0);

	theta1 = u2x_ * flexValue_ (0) + u2y_ * flexValue_ (1);
	dtheta1 = u2x_ * flexDeriv_ (0) + u2y_ * flexDeriv_ (1);
	lat = u2x_*x + u2y_*y;
	dlat = u2x_*dcom_ (0) + u2y_*dcom_ (1);
	ddlat = - (gainLat_ (0)*lat + gainLat_ (1)*(theta1-theta1Ref_)
		   + gainLat_ (2)*dlat + gainLat_ (3)*dtheta1);

	d2com_ (0) = ddxi * u1x_ + ddlat*u2x_;
	d2com_ (1) = ddxi * u1y_ + ddlat*u2y_;
	dcom_ (0) += dt_ * d2com_ (0);
	dcom_ (1) += dt_ * d2com_ (1);

	// along z
	d2com_ (2) = - (gainz_ (0)*z + gainz_ (1)*thetaz +
			gainz_ (2)*dcom_ (2) + gainz_ (3)*dthetaz);
	dcom_ (2) += dt_ * d2com_ (2);
	break;
      default:
	break;
      };

      comdot.resize (3);
      comdot [0].setSingleBound (comdotRef (0) + dcom_ (0));
      comdot [1].setSingleBound (comdotRef (1) + dcom_ (1));
      comdot [2].setSingleBound (comdotRef (2) + dcom_ (2));

      debug_ (0) = deltaCom (0);
      debug_ (1) = deltaCom (1);
      debug_ (2) = dcom_ (0);
      debug_ (3) = dcom_ (1);

      d2comSOUT_.setConstant (d2com_);
      d2comSOUT_.setTime (time);

      cosineLeftFootXSOUT_.setConstant (cosineLeftFootX_);
      cosineLeftFootXSOUT_.setTime (time);
      cosineLeftFootYSOUT_.setConstant (cosineLeftFootY_);
      cosineLeftFootYSOUT_.setTime (time);
      cosineRightFootXSOUT_.setConstant (cosineRightFootX_);
      cosineRightFootXSOUT_.setTime (time);
      cosineRightFootYSOUT_.setConstant (cosineRightFootY_);
      cosineRightFootYSOUT_.setTime (time);
      nbSupportSOUT_.setConstant (nbSupport_);
      nbSupportSOUT_.setTime (time);
      debugSOUT_.setConstant (debug_);
      debugSOUT_.setTime (time);

      return comdot;
    }

    Matrix& Stabilizer::computeJacobianCom(Matrix& jacobian, const int& time)
    {
      typedef unsigned int size_t;
      jacobian = jacobianSIN_ (time);
      return jacobian;
    }

  } // namespace dynamic
} // namespace sot
