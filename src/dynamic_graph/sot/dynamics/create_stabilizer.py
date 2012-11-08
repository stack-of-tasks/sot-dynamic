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
    Multiply_of_matrixHomo, Selec_of_vector, Kalman
from dynamic_graph.sot.dynamics import Stabilizer, flexibility_f, flexibility_h
from dynamic_graph import plug

x0 = (0., 0, 0., 0., 425.,)
P0 = ((0.02, 0., 0., 0., 0.,),
      (0., 0.02, 0., 0., 0.,),
      (0., 0., 0.03, 0., 0.,),
      (0., 0., 0., 0.03, 0.,),
      (0., 0., 0., 0., 100.,))
Q = ((.00001, 0., 0., 0., 0.,),
     (0., 0., 0., 0., 0.,),
     (0., 0., 0., 0., 0.,),
     (0., 0., 0., 0., 0.,),
     (0., 0., 0., 0., 0.,),)

R = ((1.,),)

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
    # Position of left foot
    prodLeft = Multiply_of_matrixHomo (robot.name + '_prodLeft')
    plug (robot.leftAnkle.position, prodLeft.sin1)
    prodLeft.sin2.value = robot.forceSensorInLeftAnkle
    plug (prodLeft.sout, robot.stabilizer.leftFootPosition)
    # position of right foot
    prodRight = Multiply_of_matrixHomo (robot.name + '_prodRight')
    plug (robot.rightAnkle.position, prodRight.sin1)
    prodRight.sin2.value = robot.forceSensorInRightAnkle
    plug (prodRight.sout, robot.stabilizer.rightFootPosition)
    # Kalman filter along x axis
    robot.ekf_x = Kalman (robot.name + '_EKF_x')
    robot.fx = flexibility_f (robot.name + '_fx')
    robot.hx = flexibility_h (robot.name + '_hx')
    plug (robot.stabilizer.ddx, robot.fx.control)
    plug (robot.ekf_x.x_est, robot.fx.state)
    plug (robot.ekf_x.x_est, robot.stabilizer.stateFlex_x)
    plug (robot.fx.newState, robot.hx.state)
    plug (robot.fx.newState, robot.ekf_x.x_pred)
    plug (robot.hx.observation, robot.ekf_x.y_pred)
    plug (robot.fx.jacobian, robot.ekf_x.F)
    plug (robot.hx.jacobian, robot.ekf_x.H)
    robot.ekf_x.Q.value = Q
    robot.ekf_x.R.value = R
    robot.ekf_x.setInitialState (x0)
    robot.ekf_x.setInitialVariance (P0)
    robot.My = Selec_of_vector(robot.name + '_My')
    robot.My.selec (4, 5)
    plug (robot.device.forceRLEG, robot.My.sin)
    plug (robot.My.sout, robot.ekf_x.y)
    # Kalman filter along y axis
    robot.ekf_y = Kalman (robot.name + '_EKF_y')
    robot.fy = flexibility_f (robot.name + '_fy')
    robot.hy = flexibility_h (robot.name + '_hy')
    plug (robot.stabilizer.ddy, robot.fy.control)
    plug (robot.ekf_y.x_est, robot.fy.state)
    plug (robot.ekf_y.x_est, robot.stabilizer.stateFlex_y)
    plug (robot.fy.newState, robot.hy.state)
    plug (robot.fy.newState, robot.ekf_y.x_pred)
    plug (robot.hy.observation, robot.ekf_y.y_pred)
    plug (robot.fy.jacobian, robot.ekf_y.F)
    plug (robot.hy.jacobian, robot.ekf_y.H)
    robot.ekf_y.Q.value = Q
    robot.ekf_y.R.value = R
    robot.ekf_y.setInitialState (x0)
    robot.ekf_y.setInitialVariance (P0)
    robot.My = Selec_of_vector(robot.name + '_Mx')
    robot.My.selec (3, 4)
    robot.mMy = Multiply_double_vector (robot.name + '_mMy')
    robot.mMy.sin1.value = -1.
    plug (robot.device.forceRLEG, robot.My.sin)
    plug (robot.My.sout, robot.mMy.sin2)
    plug (robot.mMy.sout, robot.ekf_y.y)
    
    return robot.stabilizer


