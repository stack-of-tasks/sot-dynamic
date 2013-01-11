//
// Copyright (c) 2013,
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

#include "flexibility-compensation.hh"

namespace sot {
  namespace dynamic {
    namespace flexibility {
      Compensation::Compensation (const std::string name) :
	Entity (name),
	flexibilityPositionSIN_ (0, "flexibility::Compensation("+name+
				 ")::input(MatrixHomo)::flexPosition"),
	globalReferencePositionSIN_ (0, "flexibility::Compensation("+name+
				     ")::input(MatrixHomo)::globalPosition"),
	flexibilityVelocitySIN_ (0, "flexibility::Compensation("+name+
				 ")::input(vector)::flexVelocity"),
	globalReferenceVelocitySIN_ (0, "flexibility::Compensation("+name+
				     ")::input(vector)::globalVelocity"),
	localReferencePositionSOUT_ ("flexibility::Compensation("+name+
				     ")::output(MatrixHomo)::localPosition"),
	localReferenceVelocitySOUT_ ("flexibility::Compensation("+name+
				     "output(vector)::localVelocity"),
	tFlex_ (3), vFlex_ (3), omegaFlex_ (3),
	tGlob_ (3), vGlob_ (3), omegaGlob_ (3),
	vLoc_ (3), omegaLoc_ (3), Rflex_ (), RflexT_ ()
      {
	signalRegistration (flexibilityPositionSIN_
			    << globalReferencePositionSIN_
			    << flexibilityVelocitySIN_
			    << globalReferenceVelocitySIN_
			    << localReferencePositionSOUT_
			    << localReferenceVelocitySOUT_);
	localReferencePositionSOUT_.addDependency (flexibilityPositionSIN_);
	localReferencePositionSOUT_.addDependency (globalReferencePositionSIN_);

	localReferenceVelocitySOUT_.addDependency (flexibilityPositionSIN_);
	localReferenceVelocitySOUT_.addDependency (globalReferencePositionSIN_);
	localReferenceVelocitySOUT_.addDependency (flexibilityVelocitySIN_);
	localReferenceVelocitySOUT_.addDependency (globalReferenceVelocitySIN_);

	localReferencePositionSOUT_.setFunction
	  (boost::bind (&Compensation::computeLocalPosition, this, _1, _2));
	localReferenceVelocitySOUT_.setFunction
	  (boost::bind (&Compensation::computeLocalVelocity, this, _1, _2));
      }
      MatrixHomogeneous&
      Compensation::computeLocalPosition (MatrixHomogeneous& Mloc,
					  const int& t)
      {
	const MatrixHomogeneous& Mglob = globalReferencePositionSIN_.access (t);
	const MatrixHomogeneous& Mflex = flexibilityPositionSIN_.access (t);
	Mloc = Mflex.inverse ()*Mglob;
	return Mloc;
      }
      Vector& Compensation::computeLocalVelocity (Vector& vloc,
						  const int& t)
      {
	const MatrixHomogeneous& Mglob = globalReferencePositionSIN_.access (t);
	const MatrixHomogeneous& Mflex = flexibilityPositionSIN_.access (t);
	const Vector& vOmegaGlob = globalReferenceVelocitySIN_.access (t);
	const Vector& vOmegaFlex = flexibilityVelocitySIN_.access (t);
	Mglob.extract (tGlob_);
	Mflex.extract (Rflex_);
	Mflex.extract (tFlex_);
	Rflex_.transpose (RflexT_);
	vOmegaGlob.extract (0, 3, vGlob_);
	vOmegaGlob.extract (3, 3, omegaGlob_);
	vOmegaFlex.extract (0, 3, vFlex_);
	vOmegaFlex.extract (3, 3, omegaFlex_);

	vLoc_ =
	  RflexT_*(omegaFlex_.crossProduct (tFlex_ - tGlob_) + vGlob_ - vFlex_);
	omegaLoc_ = RflexT_*(omegaGlob_ - omegaFlex_);
	vloc.resize (6);
	vloc (0) = vLoc_ (0);
	vloc (1) = vLoc_ (1);
	vloc (2) = vLoc_ (2);
	vloc (3) = omegaLoc_ (0);
	vloc (4) = omegaLoc_ (1);
	vloc (5) = omegaLoc_ (2);
	return vloc;
      }

      DYNAMICGRAPH_FACTORY_ENTITY_PLUGIN (Compensation,
					  "FlexibilityCompensation");
    } // namespace flexibility
  } // namespace dynamic
} // namespace sot
