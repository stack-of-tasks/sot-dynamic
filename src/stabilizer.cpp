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
#include <dynamic-graph/signal-ptr.h>
#include <dynamic-graph/command-setter.h>
#include <dynamic-graph/command-getter.h>

#include <sot/core/task-abstract.hh>
#include <sot/core/multi-bound.hh>

namespace sot {
  namespace dynamic {
    using dynamicgraph::sot::TaskAbstract;
    using dynamicgraph::SignalPtr;
    using dynamicgraph::SignalTimeDependent;
    using dynamicgraph::Vector;
    using dynamicgraph::Matrix;
    using dynamicgraph::sot::VectorMultiBound;

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
      /// Constructor by name
      Stabilizer(const std::string& inName) :
	TaskAbstract(inName),
	comSIN_ (NULL, "Stabilizer("+inName+")::input(vector)::com"),
	jacobianSIN_ (NULL, "Stabilizer("+inName+")::input(vector)::Jcom"),
	comDesSIN_ (NULL, "Stabilizer("+inName+")::input(vector)::comDes"),
	zmpSIN_ (NULL, "Stabilizer("+inName+")::input(vector)::zmp"),
	zmpDesSIN_ (NULL, "Stabilizer("+inName+")::input(vector)::zmpDes"),
	comdotSIN_ (NULL, "Stabilizer("+inName+")::input(vector)::comdot"),
	comGain_(0.), zmpGain_ (0.)
      {
	// Register signals into the entity.
	signalRegistration (comSIN_);
	signalRegistration (jacobianSIN_);
	signalRegistration (comDesSIN_);
	signalRegistration (zmpSIN_);
	signalRegistration (zmpDesSIN_);
	signalRegistration (comdotSIN_);

	taskSOUT.addDependency (comSIN_);
	taskSOUT.addDependency (zmpSIN_);
	taskSOUT.addDependency (comDesSIN_);
	taskSOUT.addDependency (zmpDesSIN_);
	taskSOUT.addDependency (comdotSIN_);

	jacobianSOUT.addDependency (jacobianSIN_);

	taskSOUT.setFunction (boost::bind(&Stabilizer::computeControlFeedback,
					  this,_1,_2));
	jacobianSOUT.setFunction (boost::bind(&Stabilizer::computeJacobianCom,
					      this,_1,_2));

	std::string docstring;
	// setGain
	docstring =
	  "\n"
	  "    Set com gain\n"
	  "      - input\n"
	  "        a floating point number\n"
	  "\n";
	addCommand(std::string("setComGain"),
		   new ::dynamicgraph::command::Setter<Stabilizer, double>
		   (*this, &Stabilizer::setComGain, docstring));
	
	// getGain
	docstring =
	  "\n"
	  "    Get gain of controller\n"
	  "      - return a floating point number\n"
	  "\n";
	addCommand(std::string("getComGain"),
		   new ::dynamicgraph::command::Getter<Stabilizer, double>
		   (*this, &Stabilizer::getComGain, docstring));
	// setGain
	docstring =
	  "\n"
	  "    Set zmp gain\n"
	  "      - input\n"
	  "        a floating point number\n"
	  "\n";
	addCommand(std::string("setZmpGain"),
		   new ::dynamicgraph::command::Setter<Stabilizer, double>
		   (*this, &Stabilizer::setZmpGain, docstring));
	
	// getGain
	docstring =
	  "\n"
	  "    Get gain of controller\n"
	  "      - return a floating point number\n"
	  "\n";
	addCommand(std::string("getZmpGain"),
		   new ::dynamicgraph::command::Getter<Stabilizer, double>
		   (*this, &Stabilizer::getZmpGain, docstring));

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
	  "  - com: the position of the center of mass,\n"
	  "  - Jcom: the Jacobian of the center of mass wrt the robot configuration,\n"
	  "  - comDes: the desired position of the center of mass,\n"
	  "  - comdot: the reference velocity of the center of mass \n"
	  "  - zmp: the measure of the Zero Momentum Point,\n"
	  "  - zmpDes: the reference value of the Zero Momentum Point.\n"
	  "  \n"
	  "As any task, the entity provide two output signals:\n"
	  "  - task: the velocity of the center of mass so as to cope with\n"
	  "          perturbations,\n"
	  "  - jacobian: the Jacobian of the center of mass with respect to robot\n"
	  "              configuration.\n";
	return doc;
      }

      /// \name Gains
      /// @{

      /// Set gain relative to center of mass
      void setComGain (const double& inGain) {
	comGain_ = inGain;
      }

      /// Get gain relative to center of mass
      double getComGain () const {
	return comGain_;
      }

      /// Set gain relative to ZMP
      void setZmpGain (const double& inGain) {
	zmpGain_ = inGain;
      }

      /// Get gain relative to ZMP
      double getZmpGain () const {
	return zmpGain_;
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
      /// Compute the control law
      VectorMultiBound& computeControlFeedback(VectorMultiBound& comdot,
					       const int& inTime)
      {
	const Vector& com = comSIN_ (inTime);
	const Vector& zmp = zmpSIN_ (inTime);
	const Vector& comDes = comDesSIN_ (inTime);
	const Vector& zmpDes = zmpDesSIN_ (inTime);
	const Vector& comdotRef = comdotSIN_ (inTime);

	double errorXzmp = zmp (0) - zmpDes (0);
	double errorYzmp = zmp (1) - zmpDes (1);
	double errorXcom = com (0) - comDes (0);
	double errorYcom = com (1) - comDes (1);
	double errorZcom = com (2) - comDes (2);
	double correctionX = zmpGain_ * errorXzmp - comGain_ * errorXcom;
	double correctionY = zmpGain_ * errorYzmp - comGain_ * errorYcom;
	double correctionZ =                      -comGain_ * errorZcom;

	comdot.resize (3);
	comdot [0].setSingleBound (comdotRef (0) + correctionX);
	comdot [1].setSingleBound (comdotRef (1) + correctionY);
	comdot [2].setSingleBound (comdotRef (2) + correctionZ);

	return comdot;
      }

      Matrix& computeJacobianCom(Matrix& jacobian, const int& inTime)
      {
	typedef unsigned int size_t;
	jacobian = jacobianSIN_ (inTime);
	return jacobian;
      }

      /// Position of center of mass
      SignalPtr < ::dynamicgraph::Vector, int> comSIN_;
      /// Position of center of mass
      SignalPtr < ::dynamicgraph::Matrix, int> jacobianSIN_;
      /// Desired position of center of mass
      SignalPtr < ::dynamicgraph::Vector, int> comDesSIN_;
      /// Position of ZMP
      SignalPtr < ::dynamicgraph::Vector, int> zmpSIN_;
      /// Desired position of ZMP
      SignalPtr < ::dynamicgraph::Vector, int> zmpDesSIN_;
      /// Reference velocity of the center of mass
      SignalPtr < ::dynamicgraph::Vector, int> comdotSIN_;

      /// Gain of the controller relative to center of mass
      double comGain_;
      /// Gain of the controller relative to center of pressure
      double zmpGain_;
      // Time sampling period
      double timePeriod_;
    }; // class Stabilizer

    DYNAMICGRAPH_FACTORY_ENTITY_PLUGIN (Stabilizer, "Stabilizer");
  } // namespace dynamic
} // namespace sot
