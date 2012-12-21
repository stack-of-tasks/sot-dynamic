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
      stateFlexInvSOUT_
      ("Stabilizer("+inName+")::output(MatrixHomo)::stateFlex_inv"),
      correctionLfSOUT_
      ("Stabilizer("+inName+")::output(MatrixHomo)::correction_lf"),
      correctionRfSOUT_
      ("Stabilizer("+inName+")::output(MatrixHomo)::correction_rf"),
      flexZobsSOUT_ ("Stabilizer("+inName+")::output(Vector)::flexZobs"),
      debugSOUT_ ("Stabilizer("+inName+")::output(vector)::debug"),
      gain1_ (4), gain2_ (4),
      prevCom_(3), flexAngle_ (2), flexDeriv_ (2),
      dcom_ (3), dt_ (.005), on_ (false),
      forceThreshold_ (.036*m_*g_), angularStiffness_ (425.), d2com_ (3),
      deltaCom_ (3),
      cosineLeftFootX_ (0.), cosineLeftFootY_ (0.),
      cosineRightFootX_ (0.), cosineRightFootY_ (0.),
      stateFlex_ (), stateFlexInv_ (), correctionLf_ (), correctionRf_ (),
      flexZobs_ (2),
      timeBeforeFlyingFootCorrection_ (.1),
      iterationsSinceLastSupportLf_ (0), iterationsSinceLastSupportRf_ (0),
      uth_ (),
      R_ (), translation_ (3), zmp_ (3), debug_ (4)
    {
      // Register signals into the entity.
      signalRegistration (deltaComSIN_);
      signalRegistration (jacobianSIN_);
      signalRegistration (comdotSIN_);
      signalRegistration (leftFootPositionSIN_ << rightFootPositionSIN_
			  << forceRightFootSIN_ << forceLeftFootSIN_);
      signalRegistration (stateFlexRfxSIN_ << stateFlexRfySIN_
			  << stateFlexLfxSIN_ << stateFlexLfySIN_);
      signalRegistration (controlGainSIN_);
      signalRegistration (sideGainSIN_);
      signalRegistration (d2comSOUT_);
      signalRegistration (cosineRightFootXSOUT_ << cosineRightFootYSOUT_
			  << cosineLeftFootXSOUT_ << cosineLeftFootYSOUT_);
      signalRegistration (nbSupportSOUT_);
      signalRegistration (stateFlexInvSOUT_ << correctionLfSOUT_
			  << correctionRfSOUT_ << flexZobsSOUT_);
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
      cosineRightFootXSOUT_.addDependency (taskSOUT);
      cosineRightFootYSOUT_.addDependency (taskSOUT);
      cosineLeftFootXSOUT_.addDependency (taskSOUT);
      cosineLeftFootYSOUT_.addDependency (taskSOUT);
      nbSupportSOUT_.addDependency (taskSOUT);

      d2com_.setZero ();
      dcom_.setZero ();
      deltaCom_.setZero ();
      d2comSOUT_.setConstant (d2com_);
      sideGainSIN_.setConstant (1.);
      debug_.setZero ();
      debugSOUT_.setConstant (debug_);

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

      // Single support gains for
      //  - kth = 510,
      //  - zeta = .8
      //  - m = 56
      //  - eigen values = 3, 3, 6, 6.
      gain1_ (0) = 121.89726;
      gain1_ (1) = -5.4917365714285769;
      gain1_ (2) = 38.280282352941178;
      gain1_ (3) = -16.224225882352943;

      // Double support gains for
      //  - kth = 2*510,
      //  - zeta = .8
      //  - m = 56
      //  - eigen values = 4, 4, 8, 8.
      gain2_ (0) = 118.62268196078432;
      gain2_ (1) = 58.543997288515406;
      gain2_ (2) = 37.326305882352941;
      gain2_ (3) = -10.661044705882359;

      zmp_.setZero ();
    }

    /// Compute flexibility state from both feet
    void Stabilizer::computeFlexibility (const int& time)
    {
      const Vector& flexRfx = stateFlexRfxSIN_.access (time);
      const Vector& flexRfy = stateFlexRfySIN_.access (time);
      const Vector& flexLfx = stateFlexLfxSIN_.access (time);
      const Vector& flexLfy = stateFlexLfySIN_.access (time);
      const MatrixHomogeneous& Mr = rightFootPositionSIN_.access (time);
      const MatrixHomogeneous& Ml = leftFootPositionSIN_.access (time);
      const Vector& fr = forceRightFootSIN_.access (time);
      const Vector& fl = forceLeftFootSIN_.access (time);
      double deltaComRfx, deltaComRfy, deltaComLfx, deltaComLfy;
      double dcomRfx, dcomRfy, dcomLfx, dcomLfy;

      // Express vertical component of force in global basis
      double flz = Ml (2,0) * fl (0) + Ml(2,1) * fl (1) + Ml (2,2) * fl (2);
      double frz = Mr (2,0) * fr (0) + Mr(2,1) * fr (1) + Mr (2,2) * fr (2);

      nbSupport_ = 0;
      if (on_) {
	if (frz >= forceThreshold_) {
	  nbSupport_++;
	  iterationsSinceLastSupportRf_ = 0;
	} else {
	  iterationsSinceLastSupportRf_ ++;
	}
	if (flz >= forceThreshold_) {
	  nbSupport_++;
	  iterationsSinceLastSupportLf_ = 0;
	} else {
	  iterationsSinceLastSupportLf_++;
	}
      }
      if (frz < 0) frz = 0;
      if (flz < 0) flz = 0;
      double fz = flz + frz;
      flexZobs_ (1) = fz;
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

      flexAngle_ (0) = (frz * flexAngleRfx + flz * flexAngleLfx)/fz;
      flexAngle_ (1) = (frz * flexAngleRfy + flz * flexAngleLfy)/fz;
      flexDeriv_ (0) = (frz * flexDerivRfx + flz * flexDerivLfx)/fz;
      flexDeriv_ (1) = (frz * flexDerivRfy + flz * flexDerivLfy)/fz;
      // Compute deviation of center of mass
      deltaComLfx = cth * flexLfx (0) - sth * flexLfy (0);
      deltaComLfy = sth * flexLfx (0) + cth * flexLfy (0);
      dcomLfx = cth * flexLfx (2) - sth * flexLfy (2);
      dcomLfy = sth * flexLfx (2) + cth * flexLfy (2);

      deltaCom_ (0) = (frz * deltaComRfx + flz * deltaComLfx)/fz;
      deltaCom_ (1) = (frz * deltaComRfy + flz * deltaComLfy)/fz;
      dcom_ (0) = (frz * dcomRfx + flz * dcomLfx)/fz;
      dcom_ (1) = (frz * dcomRfy + flz * dcomLfy)/fz;
      // Compute inverses of flexibility transformations
      if (fz != 0) {
	zmp_ (0) = (frz * Mr (0, 3) + flz * Ml (0, 3))/fz;
	zmp_ (1) = (frz * Mr (1, 3) + flz * Ml (1, 3))/fz;
	uth_ (0) = flexAngle_ (1);
	uth_ (1) = -flexAngle_ (0);
	translation_ = zmp_ - R_ * zmp_;
	uth_.toMatrix (R_);
	for (std::size_t row = 0; row < 3; ++row) {
	  for (std::size_t col = 0; col < 3; ++col) {
	    stateFlex_ (row, col) = R_ (row, col);
	  }
	  stateFlex_ (row, 3) = translation_ (row);
	  stateFlex_.inverse (stateFlexInv_);
	}
	if (iterationsSinceLastSupportLf_ * dt_ >
	    timeBeforeFlyingFootCorrection_) {
	  correctionLf_ = stateFlexInv_;
	} else {
	  correctionLf_.setIdentity ();
	}
	if (iterationsSinceLastSupportRf_ * dt_ >
	    timeBeforeFlyingFootCorrection_) {
	  correctionRf_ = stateFlexInv_;
	} else {
	  correctionRf_.setIdentity ();
	}
	stateFlexInvSOUT_.setConstant (stateFlexInv_);
	correctionLfSOUT_.setConstant (correctionLf_);
	correctionRfSOUT_.setConstant (correctionRf_);
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

      computeFlexibility (time);

      double x = deltaCom_ (0);
      double y = deltaCom_ (1);
      double z = deltaCom (2);

      // z-component of center of mass deviation in global frame
      flexZobs_ (0) = stateFlex_ (2, 0) * x + stateFlex_ (2, 1) * y +
	stateFlex_ (2, 2) * z;
      flexZobsSOUT_.setConstant (flexZobs_);

      double theta0, dtheta0, norm, u2x, u2y, u1x, u1y;
      double theta1, dtheta1, delta_x, delta_y, ddxi;

      // compute component of angle orthogonal to the line joining the feet
      delta_x = leftFootPosition (0, 3) - rightFootPosition (0, 3);
      delta_y = leftFootPosition (1, 3) - rightFootPosition (1, 3);
      norm = sqrt (delta_x*delta_x+delta_y*delta_y);
      u2x = delta_x/norm;
      u2y = delta_y/norm;
      u1x = u2y;
      u1y = -u2x;

      cosineLeftFootX_ = leftFootPosition (0, 0)*u1x +
	leftFootPosition (1, 0)*u1y;
      cosineLeftFootY_ = leftFootPosition (0, 1)*u1x +
	leftFootPosition (1, 1)*u1y;;
      cosineRightFootX_ = rightFootPosition (0, 0)*u1x +
	rightFootPosition (1, 0)*u1y;
      cosineRightFootY_ = rightFootPosition (0, 1)*u1x +
	rightFootPosition (1, 1)*u1y;;

      switch (nbSupport_) {
      case 0:
	dcom_ (0) = -gain * x;
	dcom_ (1) = -gain * y;
	dcom_ (2) = -gain * z;
	break;
      case 1: //single support
	//along x
	theta0 = flexAngle_ (0);
	dtheta0 = flexDeriv_ (0);
	d2com_ (0)= -(gain1_ (0)*x + gain1_ (1)*theta0 +
		      gain1_ (2)*dcom_ (0) + gain1_ (3)*dtheta0);
	dcom_ (0) += dt_ * d2com_ (0);
	// along y
	theta1 = flexAngle_ (1);
	dtheta1 = flexDeriv_ (1);
	d2com_ (1) = - (gain1_ (0)*y + gain1_ (1)*theta1 +
			gain1_ (2)*dcom_ (1) + gain1_ (3)*dtheta1);
	dcom_ (1) += dt_ * d2com_ (1);
	d2com_ (2) = -2*gain*dcom_ (2) - gain*gain*z;
	dcom_ (2) += dt_ * d2com_ (2);
	break;
      case 2: //double support
	//along x
	theta0 = flexAngle_ (0);
	dtheta0 = flexDeriv_ (0);
	d2com_ (0)= -(gain2_ (0)*x + gain2_ (1)*theta0 +
		      gain2_ (2)*dcom_ (0) + gain2_ (3)*dtheta0);
	dcom_ (0) += dt_ * d2com_ (0);
	// along y
	theta1 = flexAngle_ (1);
	dtheta1 = flexDeriv_ (1);
	d2com_ (1) = - (gain2_ (0)*y + gain2_ (1)*theta1 +
			gain2_ (2)*dcom_ (1) + gain2_ (3)*dtheta1);
	dcom_ (1) += dt_ * d2com_ (1);
	d2com_ (2) = -2*gain*dcom_ (2) - gain*gain*z;
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
