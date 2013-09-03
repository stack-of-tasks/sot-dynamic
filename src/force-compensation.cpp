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

#include <sot-dynamic/force-compensation.h>
#include <sot/core/debug.hh>
#include <dynamic-graph/factory.h>
#include <sot/core/macros-signal.hh>

using namespace dynamicgraph::sot;
using namespace dynamicgraph;
DYNAMICGRAPH_FACTORY_ENTITY_PLUGIN(ForceCompensationPlugin,"ForceCompensation");

/* --- PLUGIN --------------------------------------------------------------- */
/* --- PLUGIN --------------------------------------------------------------- */
/* --- PLUGIN --------------------------------------------------------------- */
ForceCompensation::
ForceCompensation(void)
  :usingPrecompensation(false)
{}


ForceCompensationPlugin::
ForceCompensationPlugin( const std::string & name )
  :Entity(name)
  ,calibrationStarted(false)

  ,torsorSIN(NULL,"sotForceCompensation("+name+")::input(vector6)::torsorIN")
  ,worldRhandSIN(NULL,"sotForceCompensation("+name+")::input(MatrixRotation)::worldRhand")

  ,handRsensorSIN(NULL,"sotForceCompensation("+name+")::input(MatrixRotation)::handRsensor")
  ,translationSensorComSIN(NULL,"sotForceCompensation("+name+")::input(vector3)::sensorCom")
  ,gravitySIN(NULL,"sotForceCompensation("+name+")::input(vector6)::gravity")
  ,precompensationSIN(NULL,"sotForceCompensation("+name+")::input(vector6)::precompensation")
  ,gainSensorSIN(NULL,"sotForceCompensation("+name+")::input(matrix6)::gain")
  ,deadZoneLimitSIN(NULL,"sotForceCompensation("+name+")::input(vector6)::deadZoneLimit")
  ,transSensorJointSIN(NULL,"sotForceCompensation("+name+")::input(vector6)::sensorJoint")  ,inertiaJointSIN(NULL,"sotForceCompensation("+name+")::input(matrix6)::inertiaJoint")
  ,velocitySIN(NULL,"sotForceCompensation("+name+")::input(vector6)::velocity")
  ,accelerationSIN(NULL,"sotForceCompensation("+name+")::input(vector6)::acceleration")

  ,handXworldSOUT( SOT_INIT_SIGNAL_2( ForceCompensation::computeHandXworld,
				      worldRhandSIN,MatrixRotation,
				      translationSensorComSIN,Vector ),
		   "sotForceCompensation("+name+")::output(MatrixForce)::handXworld" )
   ,handVsensorSOUT( SOT_INIT_SIGNAL_1( ForceCompensation::computeHandVsensor,
					handRsensorSIN,MatrixRotation),
		     "sotForceCompensation("+name+")::output(MatrixForce)::handVsensor" )
   ,torsorDeadZoneSIN(NULL,"sotForceCompensation("+name+")::input(vector6)::torsorNullifiedIN")

  ,sensorXhandSOUT( SOT_INIT_SIGNAL_2( ForceCompensation::computeSensorXhand,
				       handRsensorSIN,MatrixRotation,
				       transSensorJointSIN,Vector ),
		   "sotForceCompensation("+name+")::output(MatrixForce)::sensorXhand" )
//   ,inertiaSensorSOUT( SOT_INIT_SIGNAL_2( ForceCompensation::computeInertiaSensor,
// 					  inertiaJointSIN,Matrix,
// 					  sensorXhandSOUT,MatrixForce ),
// 		       "ForceCompensation("+name+")::output(MatrixForce)::inertiaSensor" )
   ,momentumSOUT( SOT_INIT_SIGNAL_4(ForceCompensation::computeMomentum,
				    velocitySIN,Vector,
				    accelerationSIN,Vector,
				    sensorXhandSOUT,MatrixForce,
				    inertiaJointSIN,Matrix),
		  "sotForceCompensation("+name+")::output(Vector6)::momentum" )
  ,momentumSIN(NULL,"sotForceCompensation("+name+")::input(vector6)::momentumIN")

  ,torsorCompensatedSOUT( SOT_INIT_SIGNAL_7(ForceCompensation::computeTorsorCompensated,
					    torsorSIN,Vector,
					    precompensationSIN,Vector,
					    gravitySIN,Vector,
					    handXworldSOUT,MatrixForce,
					    handVsensorSOUT,MatrixForce,
					    gainSensorSIN,Matrix,
					    momentumSIN,Vector),
			  "sotForceCompensation("+name+")::output(Vector6)::torsor" )
  ,torsorDeadZoneSOUT( SOT_INIT_SIGNAL_2(ForceCompensation::computeDeadZone,
					 torsorDeadZoneSIN,Vector,
					 deadZoneLimitSIN,Vector),
		       "sotForceCompensation("+name+")::output(Vector6)::torsorNullified" )
   ,calibrationTrigerSOUT( boost::bind(&ForceCompensationPlugin::calibrationTriger,
				       this,_1,_2),
			   torsorSIN << worldRhandSIN,
			   "sotForceCompensation("+name+")::output(Dummy)::calibrationTriger")
{
  sotDEBUGIN(5);

  signalRegistration(torsorSIN);
  signalRegistration(worldRhandSIN);
  signalRegistration(handRsensorSIN);
  signalRegistration(translationSensorComSIN);
  signalRegistration(gravitySIN);
  signalRegistration(precompensationSIN);
  signalRegistration(gainSensorSIN);
  signalRegistration(deadZoneLimitSIN);
  signalRegistration(transSensorJointSIN);
  signalRegistration(inertiaJointSIN);
  signalRegistration(velocitySIN );
  signalRegistration(accelerationSIN);
  signalRegistration(handXworldSOUT);
  signalRegistration(handVsensorSOUT);
  signalRegistration(torsorDeadZoneSIN);
  signalRegistration(sensorXhandSOUT);
  signalRegistration(momentumSOUT);
  signalRegistration(momentumSIN);
  signalRegistration(torsorCompensatedSOUT);
  signalRegistration(torsorDeadZoneSOUT);
  signalRegistration(calibrationTrigerSOUT);
  torsorDeadZoneSIN.plug(&torsorCompensatedSOUT);

  // By default, I choose: momentum is not compensated.
  //  momentumSIN.plug( &momentumSOUT );
  Vector v(6); v.fill(0); momentumSIN = v;

  sotDEBUGOUT(5);
}


ForceCompensationPlugin::
~ForceCompensationPlugin( void )
{
  return;
}

/* --- CALIB --------------------------------------------------------------- */
/* --- CALIB --------------------------------------------------------------- */
/* --- CALIB --------------------------------------------------------------- */

MatrixRotation ForceCompensation::I3;

void ForceCompensation::
clearCalibration( void )
{
  torsorList.clear();
  rotationList.clear();
}


void ForceCompensation::
addCalibrationValue( const Vector& /*torsor*/,
		     const MatrixRotation & /*worldRhand*/ )
{
  sotDEBUGIN(45);

  //   sotDEBUG(25) << "Add torsor: t"<<torsorList.size() <<" = " << torsor;
  //   sotDEBUG(25) << "Add Rotation: wRh"<<torsorList.size() <<" = " << worldRhand;

  //   torsorList.push_back(torsor);
  //   rotationList.push_back(worldRhand);

  sotDEBUGOUT(45);
}

Vector ForceCompensation::
calibrateTransSensorCom( const Vector& gravity,
			 const MatrixRotation& /*handRsensor*/ )
{
  //   sotDEBUGIN(25);

  //   Vector grav3(3);
  //   Vector Rgrav3(3),tau(3),Rtau(3);
  //   for( unsigned int j=0;j<3;++j ) { grav3(j)=gravity(j); }

  //   std::list< Vector >::iterator iterTorsor = torsorList.begin();
  //   std::list< MatrixRotation >::const_iterator iterRotation
  //     = rotationList.begin();


  //   const unsigned int NVAL = torsorList.size();
  //   if( 0==NVAL )
  //     {
  //       Vector zero(3); zero.fill(0);
  //       return zero;
  //     }

  //   if(NVAL!=rotationList.size() )
  //     {
  // 	  // TODO: ERROR THROW
  //     }
  //   Matrix torsors( 3,NVAL );
  //   Matrix gravitys( 3,NVAL );

  //   for( unsigned int i=0;i<NVAL;++i )
  //     {
  //       if( (torsorList.end()==iterTorsor)||(rotationList.end()==iterRotation) )
  // 	{
  // 	  // TODO: ERROR THROW
  // 	  break;
  // 	}
  //       const Vector & torsor = *iterTorsor;
  //       const MatrixRotation & worldRhand = *iterRotation;

  //       for( unsigned int j=0;j<3;++j ) { tau(j)=torsor(j+3); }
  //       handRsensor.multiply(tau,Rtau);
  //       worldRhand.transpose().multiply( grav3,Rgrav3 );
  //       for( unsigned int j=0;j<3;++j )
  // 	{
  // 	  torsors( j,i ) = -Rtau(j);
  // 	  gravitys( j,i ) = Rgrav3(j);
  // 	}
  //       sotDEBUG(35) << "R" << i << " = " << worldRhand;
  //       sotDEBUG(35) << "Rtau" << i << " = " << Rtau;
  //       sotDEBUG(35) << "Rg" << i << " = " << Rgrav3;

  //       iterTorsor++; iterRotation++;
  //     }

  //   sotDEBUG(35) << "Rgs = " << gravitys;
  //   sotDEBUG(35) << "Rtaus = " << torsors;

  //   Matrix gravsInv( gravitys.nbCols(),gravitys.nbRows() );
  //   sotDEBUG(25) << "Compute the pinv..." << std::endl;
  //   gravitys.pseudoInverse(gravsInv);
  //   sotDEBUG(25) << "Compute the pinv... [DONE]" << std::endl;
  //   sotDEBUG(25) << "gravsInv = " << gravsInv << std::endl;

  //   Matrix Skew(3,3);
  //   torsors.multiply( gravsInv,Skew );
  //   sotDEBUG(25) << "Skew = " << Skew << std::endl;

  //   Vector sc(3);
  //   sc(0)=(Skew(2,1)-Skew(1,2))*.5 ;
  //   sc(1)=(Skew(0,2)-Skew(2,0))*.5 ;
  //   sc(2)=(Skew(1,0)-Skew(0,1))*.5 ;
  //   sotDEBUG(15) << "SC = " << sc << std::endl;
  //   /* TAKE the antisym constraint into account inside the minimization pbm. */

  //   sotDEBUGOUT(25);
  //   return sc;
  return gravity;
}

Vector ForceCompensation::
calibrateGravity( const MatrixRotation& /*handRsensor*/,
		  bool /*precompensationCalibration*/,
		  const MatrixRotation& /*hand0RsensorArg*/ )
{
  sotDEBUGIN(25);

  //   MatrixRotation hand0Rsensor;
  //   if( &hand0Rsensor==&I3 ) hand0Rsensor.setIdentity();
  //   else hand0Rsensor=hand0RsensorArg;

  //   std::list< Vector >::const_iterator iterTorsor = torsorList.begin();
  //   std::list< MatrixRotation >::const_iterator iterRotation
  //     = rotationList.begin();

  //   const unsigned int NVAL = torsorList.size();
  //   if(NVAL!=rotationList.size() )
  //     {
  // 	  // TODO: ERROR THROW
  //     }
  //   if( 0==NVAL )
  //     {
  //       Vector zero(6); zero.fill(0);
  //       return zero;
  //     }

  //   Vector force(3),forceHand(3),forceWorld(3);
  //   Vector sumForce(3); sumForce.fill(0);

  //   for( unsigned int i=0;i<NVAL;++i )
  //     {
  //       if( (torsorList.end()==iterTorsor)||(rotationList.end()==iterRotation) )
  // 	{
  // 	  // TODO: ERROR THROW
  // 	  break;
  // 	}
  //       const Vector & torsor = *iterTorsor;
  //       const MatrixRotation & R = *iterRotation;

  //       /* The sensor read [-] the value, and the grav is [-] the sensor force.
  //        * [-]*[-] = [+] -> force = + torsor(1:3). */
  //       for( unsigned int j=0;j<3;++j ) { force(j)=-torsor(j); }
  //       handRsensor.multiply(force,forceHand);
  //       if( usingPrecompensation )
  // 	{
  // 	  Matrix R_I(3,3); R_I = R.transpose();
  // 	  R_I -= hand0Rsensor;
  // 	  R_I.pseudoInverse(.01).multiply( forceHand,forceWorld );
  // 	}
  //       else
  // 	{ R.multiply( forceHand,forceWorld ); }

  //       sotDEBUG(35) << "R(" << i << "*3+1:" << i << "*3+3,:)  = " << R << "';";
  //       sotDEBUG(35) << "rhFs(" << i << "*3+1:" << i << "*3+3) = " << forceHand;
  //       sotDEBUG(45) << "fworld(" << i << "*3+1:" << i << "*3+3) = " << forceWorld;

  //       sumForce+= forceWorld;

  //       iterTorsor++; iterRotation++;
  //     }

  //   sumForce*= (1./NVAL);
  //   sotDEBUG(35) << "Fmean = " << sumForce;
  //   sumForce.resize(6,false);
  //   sumForce(3)=sumForce(4)=sumForce(5)=0.;

  //   sotDEBUG(25)<<"mg = " << sumForce<<std::endl;

  sotDEBUGOUT(25);
  Vector sumForce(3); sumForce.fill(0);
  return sumForce;
}


/* --- SIGNALS -------------------------------------------------------------- */
/* --- SIGNALS -------------------------------------------------------------- */
/* --- SIGNALS -------------------------------------------------------------- */
MatrixForce& ForceCompensation::
computeHandXworld( const MatrixRotation & worldRhand,
		   const Vector & transSensorCom,
		   MatrixForce& res )
{
  sotDEBUGIN(35);

  sotDEBUG(25) << "wRrh = " << worldRhand <<std::endl;
  sotDEBUG(25) << "SC = " << transSensorCom <<std::endl;

  MatrixRotation R; R = worldRhand.transpose();
  MatrixHomogeneous scRw; scRw.buildFrom( R,transSensorCom);
  sotDEBUG(25) << "scMw = " << scRw <<std::endl;

  res.buildFrom( scRw );
  sotDEBUG(15) << "scXw = " << res <<std::endl;

  sotDEBUGOUT(35);
  return res;
}

MatrixForce& ForceCompensation::
computeHandVsensor( const MatrixRotation & handRsensor,
		    MatrixForce& res )
{
  sotDEBUGIN(35);

  for( unsigned int i=0;i<3;++i )
    for( unsigned int j=0;j<3;++j )
      {
	res(i,j)=handRsensor(i,j);
	res(i+3,j+3)=handRsensor(i,j);
	res(i+3,j)=0.;
	res(i,j+3)=0.;
      }

  sotDEBUG(25) << "hVs" << res << std::endl;

  sotDEBUGOUT(35);
  return res;
}

MatrixForce& ForceCompensation::
computeSensorXhand( const MatrixRotation & /*handRsensor*/,
		    const Vector & transJointSensor,
		    MatrixForce& res )
{
  sotDEBUGIN(35);

  /* Force Matrix to pass from the joint frame (ie frame located
   * at the position of the motor, in which the acc is computed by Spong)
   * to the frame SensorHand where all the forces are expressed (ie
   * frame located at the sensor position bu oriented like the hand). */

  MatrixRotation sensorRhand;  sensorRhand.setIdentity();
  //handRsensor.transpose(sensorRhand);
  MatrixHomogeneous sensorMhand;
  sensorMhand.buildFrom( sensorRhand,transJointSensor );

  res.buildFrom( sensorMhand );
  sotDEBUG(25) << "shXJ" << res << std::endl;

  sotDEBUGOUT(35);
  return res;
}

// Matrix& ForceCompensation::
// computeInertiaSensor( const Matrix& inertiaJoint,
// 		      const MatrixForce& sensorXhand,
// 		      Matrix& res )
// {
//   sotDEBUGIN(35);

//   /* Inertia felt at the sensor position, expressed in the orientation
//    * of the hand. */

//   res.resize(6,6);
//   sensorXhand.multiply( inertiaJoint,res );

//   sotDEBUGOUT(35);
//   return res;
// }


Vector& ForceCompensation::
computeTorsorCompensated( const Vector& torqueInput,
			  const Vector& torquePrecompensation,
			  const Vector& gravity,
			  const MatrixForce& handXworld,
			  const MatrixForce& handVsensor,
			  const Matrix& gainSensor,
			  const Vector& momentum,
			  Vector& res )

{
  sotDEBUGIN(35);
  /* Torsor in the sensor frame: K*(-torsred-gamma)+sVh*hXw*mg  */
  /* Torsor in the hand frame: hVs*K*(-torsred-gamma)+hXw*mg  */
  /* With gamma expressed in the sensor frame  (gamma_s = sVh*gamma_h) */

  sotDEBUG(25) << "t_nc = " << torqueInput;
  Vector torquePrecompensated(6);
  //if( usingPrecompensation )
  { torquePrecompensated = torqueInput + torquePrecompensation; }
  //else { torquePrecompensated = torqueInput; }
  sotDEBUG(25) << "t_pre = " << torquePrecompensated;

  Vector torqueS(6), torqueRH(6);
  torqueS = gainSensor * torquePrecompensated;
  res = handVsensor * torqueS;
  sotDEBUG(25) << "t_rh = " << res;

  Vector grh(6);
  grh = handXworld * gravity;
  grh *= -1;
  sotDEBUG(25) << "g_rh = " << grh;

  res+=grh;
  sotDEBUG(25) << "fcomp = " << res;

  res+=momentum;
  sotDEBUG(25) << "facc = " << res;


  /* TODO res += m xddot */

  sotDEBUGOUT(35);
  return res;
}
void crossProduct(const Vector v1, const Vector v2, Vector& res)
{
  res(0) = v1(1)*v2(2) - v1(2)*v2(1);
  res(1) = v1(2)*v2(0) - v1(0)*v2(2);
  res(2) = v1(0)*v2(1) - v1(1)*v2(0);
}

Vector& ForceCompensation::
crossProduct_V_F( const Vector& velocity,
		  const Vector& force,
		  Vector& res )
{
  /* [ v;w] x [ f;tau ] = [ w x f; v x f + w x tau ] */
  Vector v(3),w(3),f(3),tau(3);
  for( unsigned int i=0;i<3;++i )
    {
      v(i)=velocity(i); w(i) = velocity(i+3);
      f(i) = force(i); tau(i) = force(i+3);
    }
  Vector res1(3),res2a(3),res2b;
  crossProduct(w, f,res1 );
  crossProduct(v, f,res2a );
  crossProduct(w, tau,res2b );
  res2a+= res2b;

  res.resize(6);
  for( unsigned int i=0;i<3;++i )
    {
      res(i)=res1(i); res(i+3)=res2a(i);
    }
  return res;
}
				

Vector& ForceCompensation::
computeMomentum( const Vector& velocity,
		 const Vector& acceleration,
		 const MatrixForce& sensorXhand,
		 const Matrix& inertiaJoint,
		 Vector& res )
{
  sotDEBUGIN(35);

  /* Fs + Fext = I acc + V x Iv */
  Vector Iacc(6); Iacc = inertiaJoint * acceleration;
  res.resize(6); res = sensorXhand * Iacc;

  Vector Iv(6); Iv = inertiaJoint * velocity;
  Vector vIv(6); crossProduct_V_F( velocity,Iv,vIv );
  Vector XvIv(6); XvIv = sensorXhand * vIv;
  res+= XvIv;

  sotDEBUGOUT(35);
  return res;
}

Vector& ForceCompensation::
computeDeadZone( const Vector& torqueInput,
		 const Vector& deadZone,
		 Vector& res )
{
  sotDEBUGIN(35);

  if( torqueInput.size()>deadZone.size() ) return res;
  res.resize(torqueInput.size());
  for( int i=0;i<torqueInput.size();++i )
    {
      const double th = fabs(deadZone(i));
      if( (torqueInput(i)<th)&&(torqueInput(i)>-th) )
	{ res(i)=0; }
      else if(torqueInput(i)<0) res(i)=torqueInput(i)+th;
      else res(i)=torqueInput(i)-th;
    }

  sotDEBUGOUT(35);
  return res;
}


ForceCompensationPlugin::sotDummyType& ForceCompensationPlugin::
calibrationTriger( ForceCompensationPlugin::sotDummyType& dummy,int /*time*/ )
{
  //   sotDEBUGIN(45);
  //   if(! calibrationStarted ) { sotDEBUGOUT(45); return dummy=0; }

  //   addCalibrationValue( torsorSIN(time),worldRhandSIN(time) );
  //   sotDEBUGOUT(45);
  return dummy=1;
}

/* --- COMMANDLINE ---------------------------------------------------------- */
/* --- COMMANDLINE ---------------------------------------------------------- */
/* --- COMMANDLINE ---------------------------------------------------------- */

void ForceCompensationPlugin::
commandLine( const std::string& cmdLine,
	     std::istringstream& cmdArgs,
	     std::ostream& os )
{
  if( "help"==cmdLine )
    {
      os << "ForceCompensation: "
	 << "  - clearCalibration" << std::endl
	 << "  - {start|stop}Calibration [wait <time_sec>]" << std::endl
	 << "  - calibrateGravity\t[only {x|y|z}]" << std::endl
	 << "  - calibratePosition" << std::endl
	 << "  - precomp [{true|false}]:  get/set the "
	 << "precompensation due to sensor calib." << std::endl;
    }
  //   else if( "clearCalibration" == cmdLine )
  //     {
  //       clearCalibration();
  //     }
  //   else if( "startCalibration" == cmdLine )
  //     {
  //       calibrationStarted = true;
  //       cmdArgs >> std::ws;
  //       if( cmdArgs.good() )
  // 	{
  // 	  std::string cmdWait; cmdArgs>>cmdWait>>std::ws;
  // 	  if( (cmdWait == "wait")&&(cmdArgs.good()) )
  // 	    {
  // 	      double timeSec; cmdArgs >> timeSec;
  // 	      unsigned int timeMSec= (unsigned int)(round(timeSec*1000*1000));
  // 	      sotDEBUG(15) << "Calibration: wait for " << timeMSec << "us."<<std::endl;
  // 	      usleep( timeMSec );
  // 	      calibrationStarted = false;
  // 	    }
  // 	}
  //     }
  //   else if( "stopCalibration" == cmdLine )
  //     {
  //       calibrationStarted = false;
  //     }
  //   else if( "calibrateGravity" == cmdLine )
  //     {
  //       if( calibrationStarted )
  // 	{
  // 	  os<< "Calibration phase is on, stop it first."<<std::endl;
  // 	  return;
  // 	}
  //       Vector grav = calibrateGravity( handRsensorSIN.accessCopy(),
  // 					  usingPrecompensation );

  //       cmdArgs >> std::ws;
  //       if( cmdArgs.good() )
  // 	{
  // 	  std::string cmdOnly; cmdArgs>>cmdOnly>>std::ws;
  // 	  if( (cmdOnly == "only")&&(cmdArgs.good()) )
  // 	    {
  // 	      std::string xyz; cmdArgs >> xyz;
  // 	      if( "x"==xyz ) { grav(1)=grav(2)=0.; }
  // 	      else if( "y"==xyz ) { grav(0)=grav(2)=0.; }
  // 	      else if( "z"==xyz ) { grav(0)=grav(1)=0.; }
  // 	    }
  // 	}

  //       gravitySIN = grav;
  //     }
  //   else if( "calibratePosition" == cmdLine )
  //     {
  //       if( calibrationStarted )
  // 	{
  // 	  return;
  //       	  os<< "Calibration phase is on, stop it first."<<std::endl;
  // 	}
  //       Vector position(3);
  //       position = calibrateTransSensorCom( gravitySIN.accessCopy(),
  // 					  handRsensorSIN.accessCopy() );
  //       transSensorComSIN = position;
  //     }
  else if( "precomp" == cmdLine )
    {
      cmdArgs>>std::ws;
      if( cmdArgs.good() )
	{  cmdArgs >>  usingPrecompensation; }
      else { os << "precompensation = " << usingPrecompensation <<std::endl; }
    }
  else if( "compensateMomentum" == cmdLine )
    {
      cmdArgs>>std::ws;
      if( cmdArgs.good() )
	{
	  bool use;  cmdArgs >> use;
	  if( use ) momentumSIN.plug( &momentumSOUT );
	  else
	    {
	      Vector m(6); m.resize(0); momentumSIN = m;
	    }
	}
      else
	{
	  os << "compensateMomentum = " << (momentumSIN.getPtr()!=&momentumSIN)
	     <<std::endl;
	}
    }
  else{ Entity::commandLine( cmdLine,cmdArgs,os ); }


}

