# -*- coding: utf-8 -*-
# Copyright 2012, Florent Lamiraux CNRS
#
# This file is part of sot-dynamic.
# sot-dynamic is free software: you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# sot-dynamic is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Lesser Public License for more details.  You should have
# received a copy of the GNU Lesser General Public License along with
# sot-dynamic. If not, see <http://www.gnu.org/licenses/>.

from dynamic_graph.sot.core import Substract_of_vector, Multiply_double_vector,\
    Multiply_of_matrixHomo, Selec_of_vector, Stack_of_vector, \
    Inverse_of_matrixrotation, Multiply_matrix_vector, Kalman
from dynamic_graph.sot.dynamics import Stabilizer, flexibility_f, \
    flexibility_h, VarianceDoubleSupport, \
    MatrixHomoToYawOrientation
from dynamic_graph import plug

x0 = (0., 0, 0., 0., 425.,)
P0 = ((0.02, 0., 0., 0., 0.,),
      (0., 0.02, 0., 0., 0.,),
      (0., 0., 0.03, 0., 0.,),
      (0., 0., 0., 0.03, 0.,),
      (0., 0., 0., 0., 100.,))
Q = ((2e-6, 0., 0., 0., 0.,),
     (0., 1e-8, 0., 0., 0.,),
     (0., 0., 5e-5, 0., 0.,),
     (0., 0., 0., 1e-4, 0.,),
     (0., 0., 0., 0., 0.,),)

R = ((2.5e-5, 0.),(0., 1.,),)

def createKalmanOneFoot (robot, forceSensor, kalmanState, localDeltaCom,
                         locald2Com, cosineFoot, suffix):
    # Kalman filter along x axis of right foot
    # Entity creation
    ekf = Kalman (robot.name + '_EKF_' + suffix)
    f = flexibility_f (robot.name + '_f_' + suffix)
    h = flexibility_h (robot.name + '_h_' + suffix)
    control = Selec_of_vector (robot.name + '_control_' + suffix)
    obs = Stack_of_vector (robot.name + '_obs_' + suffix)
    vds = VarianceDoubleSupport (robot.name + '_variance_ds_' + suffix)
    setattr (robot, 'ekf_' + suffix, ekf)
    setattr (robot, 'f_' + suffix, f)
    setattr (robot, 'h_' + suffix, h)
    setattr (robot, 'control_' + suffix, control)
    setattr (robot, 'obs_' + suffix, obs)
    setattr (robot, 'variance_ds_' + suffix, vds)
    # Commands
    vds.varianceIn.value = Q
    vds.sigma.value = 1e-2
    ekf.R.value = R
    ekf.setInitialState (x0)
    ekf.setInitialVariance (P0)
    # Plug
    plug (cosineFoot, vds.cosineFoot)
    plug (vds.varianceOut, ekf.Q)
    plug (robot.stabilizer.nbSupport, vds.nbSupport)
    plug (locald2Com.sout, control.sin)
    plug (control.sout, f.control)
    plug (ekf.x_est, f.state)
    plug (ekf.x_est, kalmanState)
    plug (f.newState, h.state)
    plug (f.newState, ekf.x_pred)
    plug (h.observation, ekf.y_pred)
    plug (f.jacobian, ekf.F)
    plug (h.jacobian, ekf.H)
    plug (localDeltaCom.sout, obs.sin1)
    plug (robot.stabilizer.nbSupport, f.nbSupport)
    plug (cosineFoot, f.cosineFoot)
    if suffix [-1] == 'x':
        plug (forceSensor, obs.sin2)
    elif suffix [-1] == 'y':
        minusForce = Multiply_double_vector (robot.name + '_minusForce' +
                                             suffix)
        setattr (robot, 'minusForce_' + suffix, minusForce)
        minusForce.sin1.value = -1.
        plug (forceSensor, minusForce.sin2)
        plug (minusForce.sout, obs.sin2)
    else:
        raise RuntimeError ("Suffix should end by 'x' or 'y'.")
    plug (obs.sout, ekf.y)

def createKalmanFilterFeet (robot):
    createKalmanOneFoot (robot = robot, forceSensor = robot.device.forceRLEG,
                         kalmanState = robot.stabilizer.stateFlex_rfx,
                         localDeltaCom = robot.localDeltaCom_rf,
                         locald2Com = robot.locald2Com_rf,
                         cosineFoot = robot.stabilizer.cosine_rfx,
                         suffix = 'rfx')
    robot.control_rfx.selec (0, 1)
    robot.obs_rfx.selec1 (0, 1)
    robot.obs_rfx.selec2 (4, 5)
    robot.f_rfx.control.recompute (0)

    createKalmanOneFoot (robot = robot, forceSensor = robot.device.forceRLEG,
                         kalmanState = robot.stabilizer.stateFlex_rfy,
                         localDeltaCom = robot.localDeltaCom_rf,
                         locald2Com = robot.locald2Com_rf,
                         cosineFoot = robot.stabilizer.cosine_rfy,
                         suffix = 'rfy')
    robot.control_rfy.selec (1, 2)
    robot.obs_rfy.selec1 (1, 2)
    robot.obs_rfy.selec2 (3, 4)
    robot.f_rfy.control.recompute (0)

    createKalmanOneFoot (robot = robot, forceSensor = robot.device.forceLLEG,
                         kalmanState = robot.stabilizer.stateFlex_lfx,
                         localDeltaCom = robot.localDeltaCom_lf,
                         locald2Com = robot.locald2Com_lf,
                         cosineFoot = robot.stabilizer.cosine_lfx,
                         suffix = 'lfx')
    robot.control_lfx.selec (0, 1)
    robot.obs_lfx.selec1 (0, 1)
    robot.obs_lfx.selec2 (4, 5)
    robot.f_lfx.control.recompute (0)

    createKalmanOneFoot (robot = robot, forceSensor = robot.device.forceLLEG,
                         kalmanState = robot.stabilizer.stateFlex_lfy,
                         localDeltaCom = robot.localDeltaCom_lf,
                         locald2Com = robot.locald2Com_lf,
                         cosineFoot = robot.stabilizer.cosine_lfy,
                         suffix = 'lfy')
    robot.control_lfy.selec (1, 2)
    robot.obs_lfy.selec1 (1, 2)
    robot.obs_lfy.selec2 (3, 4)
    robot.f_lfy.control.recompute (0)

def createStabilizer (robot):
    robot.dynamic.com.recompute(0)
    robot.dynamic.Jcom.recompute(0)
    robot.stabilizer = Stabilizer (robot.name + '_stabilizer')
    robot.stabilizer.comdot.value = (0.,0.,0.,)
    plug (robot.dynamic.Jcom, robot.stabilizer.Jcom)
    robot.deltaCom = Substract_of_vector (robot.name + '_deltaCom')
    robot.comRef = robot.deltaCom.sin2
    robot.comRef.value = robot.dynamic.com.value
    plug (robot.dynamic.com, robot.deltaCom.sin1)
    plug (robot.deltaCom.sout, robot.stabilizer.deltaCom)
    plug (robot.device.forceRLEG, robot.stabilizer.forceRLEG)
    plug (robot.device.forceLLEG, robot.stabilizer.forceLLEG)
    # Position of left foot
    robot.leftFootPosition = Multiply_of_matrixHomo (robot.name +
                                                     '_position_lf')
    plug (robot.leftAnkle.position, robot.leftFootPosition.sin1)
    robot.leftFootPosition.sin2.value = robot.forceSensorInLeftAnkle
    plug (robot.leftFootPosition.sout, robot.stabilizer.leftFootPosition)
    robot.yawOrientation_lf = MatrixHomoToYawOrientation (robot.name +
                                                          '_yawOrientation_lf')
    plug (robot.leftFootPosition.sout, robot.yawOrientation_lf.position)
    robot.localDeltaCom_lf = Multiply_matrix_vector (robot.name +
                                                     '_localDeltaCom_lf')
    plug (robot.yawOrientation_lf.inverse, robot.localDeltaCom_lf.sin1)
    plug (robot.deltaCom.sout, robot.localDeltaCom_lf.sin2)
    robot.locald2Com_lf = Multiply_matrix_vector (robot.name + '_locald2Com_lf')
    plug (robot.yawOrientation_lf.inverse, robot.locald2Com_lf.sin1)
    plug (robot.stabilizer.d2com, robot.locald2Com_lf.sin2)
    # position of right foot
    robot.rightFootPosition = Multiply_of_matrixHomo (robot.name +
                                                      '_position_rf')
    plug (robot.rightAnkle.position, robot.rightFootPosition.sin1)
    robot.rightFootPosition.sin2.value = robot.forceSensorInRightAnkle
    plug (robot.rightFootPosition.sout, robot.stabilizer.rightFootPosition)
    robot.yawOrientation_rf = MatrixHomoToYawOrientation (robot.name +
                                                          '_yawOrientation_rf')
    plug (robot.rightFootPosition.sout, robot.yawOrientation_rf.position)
    robot.localDeltaCom_rf = Multiply_matrix_vector (robot.name +
                                                     '_localDeltaCom_rf')
    plug (robot.yawOrientation_rf.inverse, robot.localDeltaCom_rf.sin1)
    plug (robot.deltaCom.sout, robot.localDeltaCom_rf.sin2)
    robot.locald2Com_rf = Multiply_matrix_vector (robot.name + '_locald2Com_rf')
    plug (robot.yawOrientation_rf.inverse, robot.locald2Com_rf.sin1)
    plug (robot.stabilizer.d2com, robot.locald2Com_rf.sin2)

    createKalmanFilterFeet (robot)

    return robot.stabilizer
