/*
 * Copyright 2010,
 * François Bleibel,
 * Olivier Stasse,
 *
 * CNRS/AIST
 *
 * This file is part of sot-dynamic.
 * sot-dynamic is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * sot-dynamic is distributed in the hope that it will be
 * useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.  You should
 * have received a copy of the GNU Lesser General Public License along
 * with sot-dynamic.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __SOT_SOTFORCECOMPENSATION_H__
#define __SOT_SOTFORCECOMPENSATION_H__

/* --------------------------------------------------------------------- */
/* --- INCLUDE --------------------------------------------------------- */
/* --------------------------------------------------------------------- */

/* SOT */
#include <dynamic-graph/entity.h>
#include <dynamic-graph/signal-ptr.h>
#include <dynamic-graph/signal-time-dependent.h>
#include <sot/core/matrix-rotation.hh>
#include <sot/core/matrix-force.hh>
#include <sot/core/matrix-homogeneous.hh>

/* STD */
#include <string>

/* --------------------------------------------------------------------- */
/* --- API ------------------------------------------------------------- */
/* --------------------------------------------------------------------- */

#if defined (WIN32) 
#  if defined (force_compensation_EXPORTS)
#    define SOTFORCECOMPENSATION_EXPORT __declspec(dllexport)
#  else  
#    define SOTFORCECOMPENSATION_EXPORT __declspec(dllimport)
#  endif 
#else
#  define SOTFORCECOMPENSATION_EXPORT
#endif


namespace dynamicgraph { namespace sot {
    namespace dg = dynamicgraph;

    /* --------------------------------------------------------------------- */
    /* --- CLASS ----------------------------------------------------------- */
    /* --------------------------------------------------------------------- */

    class SOTFORCECOMPENSATION_EXPORT ForceCompensation
    {
    private:
      static MatrixRotation I3;
    protected:
      bool usingPrecompensation;

    public:
      ForceCompensation( void );
      static MatrixForce& computeHandXworld( 
					    const MatrixRotation & worldRhand,
					    const dg::Vector & transSensorCom,
					    MatrixForce& res );

  
      static MatrixForce& computeHandVsensor( const MatrixRotation & sensorRhand,
					      MatrixForce& res );
      static MatrixForce& computeSensorXhand( const MatrixRotation & sensorRhand,
					      const dg::Vector & transSensorCom,
					      MatrixForce& res );
      /*   static dg::Matrix& computeInertiaSensor( const dg::Matrix& inertiaJoint, */
      /* 					   const MatrixForce& sensorXhand, */
      /* 					   dg::Matrix& res ); */

      static dg::Vector& computeTorsorCompensated( const dg::Vector& torqueInput,
						   const dg::Vector& torquePrecompensation,
						   const dg::Vector& gravity,
						   const MatrixForce& handXworld,
						   const MatrixForce& handVsensor,
						   const dg::Matrix& gainSensor,
						   const dg::Vector& momentum,
						   dg::Vector& res );

      static dg::Vector& crossProduct_V_F( const dg::Vector& velocity,
					   const dg::Vector& force,
					   dg::Vector& res );
      static dg::Vector& computeMomentum( const dg::Vector& velocity,
					  const dg::Vector& acceleration,
					  const MatrixForce& sensorXhand,
					  const dg::Matrix& inertiaJoint,
					  dg::Vector& res );

      static dg::Vector& computeDeadZone( const dg::Vector& torqueInput,
					  const dg::Vector& deadZoneLimit,
					  dg::Vector& res );
  
    public: // CALIBRATION

      std::list<dg::Vector> torsorList;
      std::list<MatrixRotation> rotationList;

      void clearCalibration( void );
      void addCalibrationValue( const dg::Vector& torsor,
				const MatrixRotation & worldRhand );
  
      dg::Vector calibrateTransSensorCom( const dg::Vector& gravity,
					  const MatrixRotation& handRsensor );
      dg::Vector calibrateGravity( const MatrixRotation& handRsensor,
				   bool precompensationCalibration = false,
				   const MatrixRotation& hand0Rsensor = I3 );

    
  
    };

    /* --------------------------------------------------------------------- */
    /* --- PLUGIN ---------------------------------------------------------- */
    /* --------------------------------------------------------------------- */

    class SOTFORCECOMPENSATION_EXPORT ForceCompensationPlugin
      :public dg::Entity, public ForceCompensation
    {
    public:
      static const std::string CLASS_NAME;
      virtual const std::string& getClassName( void ) const { return CLASS_NAME; }
      bool calibrationStarted;


    public: /* --- CONSTRUCTION --- */

      ForceCompensationPlugin( const std::string& name );
      virtual ~ForceCompensationPlugin( void );

    public: /* --- SIGNAL --- */

      /* --- INPUTS --- */
      dg::SignalPtr<dg::Vector,int> torsorSIN; 
      dg::SignalPtr<MatrixRotation,int> worldRhandSIN; 

      /* --- CONSTANTS --- */
      dg::SignalPtr<MatrixRotation,int> handRsensorSIN; 
      dg::SignalPtr<dg::Vector,int> translationSensorComSIN; 
      dg::SignalPtr<dg::Vector,int> gravitySIN; 
      dg::SignalPtr<dg::Vector,int> precompensationSIN; 
      dg::SignalPtr<dg::Matrix,int> gainSensorSIN; 
      dg::SignalPtr<dg::Vector,int> deadZoneLimitSIN; 
      dg::SignalPtr<dg::Vector,int> transSensorJointSIN; 
      dg::SignalPtr<dg::Matrix,int> inertiaJointSIN; 

      dg::SignalPtr<dg::Vector,int> velocitySIN; 
      dg::SignalPtr<dg::Vector,int> accelerationSIN; 

      /* --- INTERMEDIATE OUTPUTS --- */
      dg::SignalTimeDependent<MatrixForce,int> handXworldSOUT; 
      dg::SignalTimeDependent<MatrixForce,int> handVsensorSOUT; 
      dg::SignalPtr<dg::Vector,int> torsorDeadZoneSIN; 

      dg::SignalTimeDependent<MatrixForce,int> sensorXhandSOUT;
      //dg::SignalTimeDependent<dg::Matrix,int> inertiaSensorSOUT;
      dg::SignalTimeDependent<dg::Vector,int> momentumSOUT; 
      dg::SignalPtr<dg::Vector,int> momentumSIN; 

      /* --- OUTPUTS --- */
      dg::SignalTimeDependent<dg::Vector,int> torsorCompensatedSOUT; 
      dg::SignalTimeDependent<dg::Vector,int> torsorDeadZoneSOUT;

      typedef int sotDummyType;
      dg::SignalTimeDependent<sotDummyType,int> calibrationTrigerSOUT; 

    public: /* --- COMMANDLINE --- */

      sotDummyType& calibrationTriger( sotDummyType& dummy,int time );


      virtual void commandLine( const std::string& cmdLine,
				std::istringstream& cmdArgs,
				std::ostream& os );



    };


  } // namaspace sot
} // namespace dynamicgraph

#endif // #ifndef __SOT_SOTFORCECOMPENSATION_H__
