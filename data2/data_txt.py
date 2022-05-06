#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import rospy
import numpy as np
from geometry_msgs.msg import TwistStamped
from geometry_msgs.msg import PoseStamped
from sensor_msgs.msg import Imu
from nlink_parser.msg import LinktrackNodeframe3
from nav_msgs.msg import Odometry

bagnumber = 6

class data_txt:
    def __init__(self):

        self.imu_uav = rospy.Subscriber("/mavros/imu/data", Imu, self.callback1)
        
        self.vel_uav = rospy.Subscriber("/mavros/local_position/velocity_body", TwistStamped, self.callback2)
        
        self.uwb_uav = rospy.Subscriber("/nlink_linktrack_nodeframe3", LinktrackNodeframe3, self.callback3)
        
        self.vi_pos_base = rospy.Subscriber("/yyjbase/viconros/mocap/pos", PoseStamped, self.callback4)
        self.vi_pos_fly = rospy.Subscriber("/yyjfly/viconros/mocap/pos", PoseStamped, self.callback5)
        
        self.vi_vel_base = rospy.Subscriber("/yyjbase/viconros/mocap/vel", TwistStamped, self.callback6)
        self.vi_vel_fly = rospy.Subscriber("/yyjfly/viconros/mocap/vel", TwistStamped, self.callback7)
        
        self.vi_acc_base = rospy.Subscriber("/yyjbase/viconros/mocap/acc", TwistStamped, self.callback8)
        self.vi_acc_fly = rospy.Subscriber("/yyjfly/viconros/mocap/acc", TwistStamped, self.callback9)

        self.imu_acce_uav = np.zeros((3,1))
        self.imu_ang_uav = np.zeros((3,1))
        self.imu_qua_uav = np.zeros((4,1))

        self.v_uav = np.zeros((3,1))
        
        self.dis_uav = 0.0

        self.base_pos = np.zeros((3,1))
        self.fly_pos = np.zeros((3,1))

        self.base_vel = np.zeros((3,1))
        self.fly_vel = np.zeros((3,1))

        self.base_acc = np.zeros((3,1))
        self.fly_acc = np.zeros((3,1))


        self.t_fly = 0.0
        self.t_vicon = 0.0

    def callback1(self, imu):
        # the linear acceleration from IMU
        self.imu_acce_uav[0, 0] = imu.linear_acceleration.x
        self.imu_acce_uav[1, 0] = imu.linear_acceleration.y
        self.imu_acce_uav[2, 0] = imu.linear_acceleration.z

        # the angular velocity from IMU
        # self.imu_ang_uav[0, 0] = imu.angular_velocity.x
        # self.imu_ang_uav[1, 0] = imu.angular_velocity.y
        # self.imu_ang_uav[2, 0] = imu.angular_velocity.z

        # the quaternion from IMU
        self.imu_qua_uav[0, 0] = imu.orientation.x
        self.imu_qua_uav[1, 0] = imu.orientation.y
        self.imu_qua_uav[2, 0] = imu.orientation.z
        self.imu_qua_uav[3, 0] = imu.orientation.w

    def callback2(self, vel):
        # the linear velocity from optical flow method
        self.v_uav[0, 0] = vel.twist.linear.x
        self.v_uav[1, 0] = vel.twist.linear.y
        self.v_uav[2, 0] = vel.twist.linear.z

        self.t_fly = vel.header.stamp.nsecs

        uav_t = str(self.t_fly)
        ax_uav = str(self.imu_acce_uav[0, 0])
        ay_uav = str(self.imu_acce_uav[1, 0])
        az_uav = str(self.imu_acce_uav[2, 0])

        qx_uav = str(self.imu_qua_uav[0, 0])
        qy_uav = str(self.imu_qua_uav[1, 0])
        qz_uav = str(self.imu_qua_uav[2, 0])
        qw_uav = str(self.imu_qua_uav[3, 0])

        vx_uav = str(self.v_uav[0, 0])
        vy_uav = str(self.v_uav[1, 0])
        vz_uav = str(self.v_uav[2, 0])

        d_uav = str(self.dis_uav)

        vicon_t = str(self.t_vicon)
        vicon_x = str(self.fly_pos[0, 0]-self.base_pos[0, 0])
        vicon_y = str(self.fly_pos[1, 0]-self.base_pos[1, 0])
        vicon_z = str(self.fly_pos[2, 0]-self.base_pos[2, 0])
        vicon_vx = str(self.fly_vel[0, 0]-self.base_vel[0, 0])
        vicon_vy = str(self.fly_vel[1, 0]-self.base_vel[1, 0])
        vicon_vz = str(self.fly_vel[2, 0]-self.base_vel[2, 0])
        vicon_ax = str(self.fly_acc[0, 0]-self.base_acc[0, 0])
        vicon_ay = str(self.fly_acc[1, 0]-self.base_acc[1, 0])
        vicon_az = str(self.fly_acc[2, 0]-self.base_acc[2, 0])

        file_handle = open("/home/jojo/Desktop/yyj_fly/uav%d.txt" %(bagnumber), "a+")
        file_handle.write(uav_t)
        file_handle.write('    ')
        
        file_handle.write(ax_uav)
        file_handle.write('    ')
        file_handle.write(ay_uav)
        file_handle.write('    ')
        file_handle.write(az_uav)
        file_handle.write('    ')

        file_handle.write(qx_uav)
        file_handle.write('    ')
        file_handle.write(qy_uav)
        file_handle.write('    ')
        file_handle.write(qz_uav)
        file_handle.write('    ')        
        file_handle.write(qw_uav)
        file_handle.write('    ')

        file_handle.write(vx_uav)
        file_handle.write('    ')
        file_handle.write(vy_uav)
        file_handle.write('    ')
        file_handle.write(vz_uav)
        file_handle.write('    ')

        file_handle.write(d_uav)
        file_handle.write('    ')        

        file_handle.write('\n')
        file_handle.close()         

        file_handle = open("/home/jojo/Desktop/yyj_fly/vicon%d.txt" %(bagnumber), "a+")
        file_handle.write(vicon_t)
        file_handle.write('    ')

        file_handle.write(vicon_x)
        file_handle.write('    ')
        file_handle.write(vicon_y)
        file_handle.write('    ')
        file_handle.write(vicon_z)
        file_handle.write('    ')

        file_handle.write(vicon_vx)
        file_handle.write('    ')
        file_handle.write(vicon_vy)
        file_handle.write('    ')
        file_handle.write(vicon_vz)
        file_handle.write('    ')

        file_handle.write(vicon_ax)
        file_handle.write('    ')
        file_handle.write(vicon_ay)
        file_handle.write('    ')
        file_handle.write(vicon_az)
        file_handle.write('    ')

        file_handle.write('\n')
        file_handle.close() 


        

    def callback3(self, uwb):
        # the distance from UWB
        i = len(uwb.nodes)
        if i == 1:
            self.dis_uav = uwb.nodes[0].dis
        else:
            return
    
    def callback4(self, mocap):
        # the position from VICON
        self.base_pos[0, 0] = mocap.pose.position.x
        self.base_pos[1, 0] = mocap.pose.position.y
        self.base_pos[2, 0] = mocap.pose.position.z

        self.t_vicon = mocap.header.stamp.nsecs
    
    def callback5(self, mocap):
        # the position from VICON
        self.fly_pos[0, 0] = mocap.pose.position.x
        self.fly_pos[1, 0] = mocap.pose.position.y
        self.fly_pos[2, 0] = mocap.pose.position.z


    def callback6(self, mocap):
        # the velocity from VICON
        self.base_vel[0, 0] = mocap.twist.linear.x
        self.base_vel[1, 0] = mocap.twist.linear.y
        self.base_vel[2, 0] = mocap.twist.linear.z

    def callback7(self, mocap):
        # the velocity from VICON
        self.fly_vel[0, 0] = mocap.twist.linear.x
        self.fly_vel[1, 0] = mocap.twist.linear.y
        self.fly_vel[2, 0] = mocap.twist.linear.z

    def callback8(self, mocap):
        # the acceleration from VICON
        self.base_acc[0, 0] = mocap.twist.linear.x
        self.base_acc[1, 0] = mocap.twist.linear.y
        self.base_acc[2, 0] = mocap.twist.linear.z

    def callback9(self, mocap):
        # the acceleration from VICON
        self.fly_acc[0, 0] = mocap.twist.linear.x
        self.fly_acc[1, 0] = mocap.twist.linear.y
        self.fly_acc[2, 0] = mocap.twist.linear.z


if __name__ == '__main__':
    file1 = open("/home/jojo/Desktop/yyj_fly/uav%d.txt" %(bagnumber), "r+")
    file1.truncate()

    file2 = open("/home/jojo/Desktop/yyj_fly/vicon%d.txt" %(bagnumber), "r+")
    file2.truncate()

    try:
        rospy.init_node("data_txt")
        rospy.loginfo("starting the node")
        data_txt()
        rospy.spin()
    
    except KeyboardInterrupt:
        print("shutting down the node")
       