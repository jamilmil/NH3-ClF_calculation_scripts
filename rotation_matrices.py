import numpy as np
from math import cos, sin, radians
R=np.zeros((3,3))
def x_axis_rot(i):
    angle=radians(i)
    R[0][0]=1
    R[0][1]=0
    R[0][2]=0
    R[1][0]=0
    R[1][1]=cos(angle)
    R[1][2]=-sin(angle)
    R[2][0]=0
    R[2][1]=sin(angle)
    R[2][2]=cos(angle)
    return R
def y_axis_rot(i):
    angle=radians(i)
    R[0][0]=cos(angle)
    R[0][1]=0
    R[0][2]=sin(angle)
    R[1][0]=0
    R[1][1]=1
    R[1][2]=0
    R[2][0]=-sin(angle)
    R[2][1]=0
    R[2][2]=cos(angle)
    return R
def z_axis_rot(i):
    angle=radians(i)
    R[0][0]=cos(angle)
    R[0][1]=-sin(angle)
    R[0][2]=0
    R[1][0]=sin(angle)
    R[1][1]=cos(angle)
    R[1][2]=0
    R[2][0]=0
    R[2][1]=0
    R[2][2]=1
    return R