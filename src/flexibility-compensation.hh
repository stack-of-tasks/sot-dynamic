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

#ifndef SOT_DYNAMIC_FLEXBILITY_COMPENSATION_HH
# define SOT_DYNAMIC_FLEXBILITY_COMPENSATION_HH

# include <dynamic-graph/linear-algebra.h>
# include <dynamic-graph/entity.h>
# include <dynamic-graph/factory.h>
# include <dynamic-graph/signal-time-dependent.h>
# include <dynamic-graph/signal-ptr.h>
# include <sot/core/matrix-homogeneous.hh>
# include <sot/core/matrix-rotation.hh>

namespace sot {
  namespace dynamic {
    namespace flexibility {

      using dynamicgraph::Entity;
      using dynamicgraph::SignalTimeDependent;
      using dynamicgraph::SignalPtr;
      using dynamicgraph::Vector;
      using dynamicgraph::sot::MatrixHomogeneous;
      using dynamicgraph::sot::MatrixRotation;

      class Compensation : public Entity
      {
	SignalPtr <MatrixHomogeneous, int> flexibilityPositionSIN_;
	SignalPtr <MatrixHomogeneous, int> globalReferencePositionSIN_;
	SignalPtr <Vector, int> flexibilityVelocitySIN_;
	SignalPtr <Vector, int> globalReferenceVelocitySIN_;
	SignalTimeDependent <MatrixHomogeneous, int>
	localReferencePositionSOUT_;
	SignalTimeDependent <Vector, int> localReferenceVelocitySOUT_;

	DYNAMIC_GRAPH_ENTITY_DECL ();
	Compensation (const std::string name);
	virtual std::string getDocString () const
	{
	  std::string doc =
	    "Ankle flexibility compensation\n"
	    "\n"
	    "This entity computes the reference position and velocity of a trajectory in the\n"
	    "local frame of the flexibility from\n"
	    "  - the reference position and velocity in the global reference frame and\n"
	    "  - the estimated position and velocity of the flexibility.\n"
	    "\n"
	    "The entity has 4 input and 2 output signals:\n"
	    "  - flexPosition: position of the flexibility (matrix homogeneous)\n"
	    "  - flexVelocity: velocity of the flexibility (derivative of the former)\n"
	    "  - globalPosition: reference position of the feature in the global reference\n"
	    "    frame,\n"
	    "  - globalVelocity: reference velocity of the feature in the global reference\n"
	    "    frame\n"
	    "  - localPosition: reference position of the feature in the flexibility\n"
	    "    reference frame,\n"
	    "  - localVelocity: reference velocity of the feature in the flexibility\n"
	    "    reference frame.\n";
	  return doc;
	}
	MatrixHomogeneous&
	computeLocalPosition (MatrixHomogeneous&, const int&);
	Vector&	computeLocalVelocity (Vector&, const int&);

	Vector tFlex_;
	Vector vFlex_;
	Vector omegaFlex_;
	Vector tGlob_;
	Vector vGlob_;
	Vector omegaGlob_;
	Vector vLoc_;
	Vector omegaLoc_;
	MatrixRotation Rflex_, RflexT_;
      }; // class Compensation
    } // namespace flexibility
  } // namespace dynamic
} // namespace sot
#endif // SOT_DYNAMIC_FLEXBILITY_COMPENSATION_HH
