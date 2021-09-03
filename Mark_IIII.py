# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import time
t=time.time()
import numpy as np
import scipy.optimize
import math
import matplotlib.pyplot as plt

def cosd(theta):
    return np.cos(theta*(math.pi/180))
def sind(theta):
    return np.sin(theta*(math.pi/180))
def fun_I(tet2,tet):
    R = np.array([203.2 , 57.15 , 184.15 , 177.8 , 50.8 , 127])
    R44=127
    R11=101.6   
    tet=np.array([np.NAN ,tet2 ,  tet[0] , tet[1] , tet[2] , tet[3]])
    [F0 , F1]=R[1]*np.array([cosd(tet[1]) , sind(tet[1])]) + R[2]*np.array([cosd(tet[2]) , sind(tet[2])]) - R[3]*np.array([cosd(tet[3]) , sind(tet[3])])-R[0]*np.array([1 , 0]) 
    [F2 ,F3]=R[5]*np.array([cosd(tet[5]) ,sind(tet[5]) ]) + R[4]*np.array([cosd(tet[4]) , sind(tet[4])]) - R44*np.array([cosd(tet[3]) , sind(tet[3])])-R11*np.array([1 , 0])
    return F0 , F1 , F2 , F3
tet2=np.arange(0 , 360 , 0.1)
n=len(tet2)
init_guess=np.zeros([n+1,4]) #hads avaliye dar har marhale v javbe marhaleye ghabl
init_guess[0 , :]=np.array([64 , 110 , 136 ,42 ])
#optimizing and finding theta
for i in range(0 ,n):
    fun=lambda tet: fun_I( tet2[i] ,tet)
    s=scipy.optimize.fsolve(fun , init_guess[i , :] );
    init_guess[i+1 , :]=s
post= init_guess[1:,:]
del init_guess, i, s
plt.figure(1)
plt.plot(tet2 , post[:,0],label = 'theta3')
plt.plot(tet2 , post[:,1] , label = 'theta4')
plt.plot(tet2 , post[:,2] , label = 'theta5')
plt.plot(tet2 , post[:,3] , label = 'theta6')
plt.show()
plt.legend()
a = np.transpose

# velocity analysis
R = np.array([203.2 , 57.15 , 184.15 , 177.8 , 50.8 , 127])
R44=127
R11=101.6
KC1 = np.zeros([n,4]) #1st order kinematic coefficient
KC2 = np.zeros([n,4])
AV = np.zeros([n,6])
def AV2(tet): #angular velocity2
    if tet>=0 and tet<50:
        a = np.sqrt((-5*math.pi*AA2(tet)/9)*(1-tet/500))
    elif tet >= 50 and tet<150:
        a = np.sqrt((10*np.pi*AA2(tet)/9)*(tet - 50)/100)
    elif tet>=150 and tet <150:
        a=np.sqrt(-10*np.pi*AA2(tet)/9*(1-(tet-150)/100))
    elif tet>=250 and tet<330:
        a=np.sqrt(8*np.pi*AA2(tet)/9*((tet-250)/80))
    else:
        a=np.sqrt(-8*np.pi*AA2(tet)/9*(1-(tet-330)/80))
    return a

def AA2(tet):
    if tet>=0 and tet<50:
        a=-5*np.pi
    elif tet>=50 and tet<150:
        a=6*np.pi
    elif tet>=150 and tet<250:
        a=-6*np.pi;
    elif tet>=250 and tet<330:
        a=5*np.pi
    else:
        a=-5*np.pi
    return a
DET1 = np.zeros([n,1])
A = np.zeros([4 ,4])
for i in range(len(KC1)):
    A[0 , :] = np.array([R[2]*-sind(post[i,0]),R[3]*sind(post[i,1]), 0 ,0 ])
    A[1 , :] = np.array([R[2]*cosd(post[i,0]), -R[3]*cosd(post[i,1]) , 0 , 0])
    A[2 , :] = np.array([0 , R44*sind(post[i , 1]),-R[4]*sind(post[i , 2]), R[5]*-sind(post[i,3])])
    A[3 , :] = np.array([0 , -R44*cosd(post[i,1]) , R[4]*cosd(post[i,2]) , R[5]*cosd(post[i , 3]) ])
    B = np.array([[-R[1]*-sind(tet2[i])],[-R[1]*cosd(tet2[i])] , [0] , [0]])
    x = np.linalg.solve(A,B)
    KC1[i] = [x[0] , x[1] , x[2] , x[3]]
    AV[i , 1] = AV2(tet2[i])
    DET1[i] = np.linalg.det(A)
AV[: , 2]=KC1[: , 0]*AV[: , 1] #%AV of link 3
AV[: , 3]=KC1[: , 1]*AV[: , 1]
AV[: , 4]=KC1[: , 2]*AV[: , 1]
AV[: , 5]=KC1[: , 3]*AV[: , 1]
##Accelartion Analysis (Ax=B)
AA=np.zeros([n , 6])
DET2=np.zeros([n , 1])
for i in range(n):
    A[0 , :] = np.array([R[2]*-sind(post[i,0]),R[3]*sind(post[i,1]), 0 ,0 ])
    A[1 , :] = np.array([R[2]*cosd(post[i,0]), -R[3]*cosd(post[i,1]) , 0 , 0])
    A[2 , :] = np.array([0 , R44*sind(post[i , 1]),-R[4]*sind(post[i , 2]), R[5]*-sind(post[i,3])])
    A[3 , :] = np.array([0 , -R44*cosd(post[i,1]) , R[4]*cosd(post[i,2]) , R[5]*cosd(post[i , 3]) ])
    B=[[-R[1]*-cosd(tet2[i])-R[2]*(KC1[i , 0])**2*-cosd(post[i , 0])+R[3]*(KC1[i , 1])**2*-cosd(post[i , 1]) ],\
       [-R[1]*-sind(tet2[i])-R[2]*(KC1[i , 0])**2*-sind(post[i , 0])+R[3]*(KC1[i , 1])**2*-sind(post[i , 1])],\
       [-R[5]*(KC1[i , 3])**2*-cosd(post[i , 3])+R[4]*(KC1[i , 2])**2*-cosd(post[i , 2]) + R44*(KC1[i , 1])**2*-cosd(post[i , 1])],\
       [-R[5]*(KC1[i , 3])**2*-sind(post[i , 3])+R[4]*(KC1[i , 2])**2*-sind(post[i , 2]) + R44*(KC1[i , 1])**2*-sind(post[i , 1])]];
    x = np.linalg.solve(A,B)
    KC2[i] = [x[0] , x[1] , x[2] , x[3]]
    AA[i , 1] = AA2(tet2[i])
    DET2[i] = np.linalg.det(A)
AA[: , 2]=KC2[: , 0]*(AV[: , 1]**2)+KC1[: , 0]*AA[: , 1]; #AA of link 3
AA[: , 3]=KC2[: , 1]*(AV[: , 1]**2)+KC1[: , 1]*AA[: , 1];
AA[: , 4]=KC2[: , 2]*(AV[: , 1]**2)+KC1[: , 2]*AA[: , 1];
AA[: , 5]=KC2[: , 3]*(AV[: , 1]**2)+KC1[: , 3]*AA[: , 1];
## calculations of points P , Z
PP=np.zeros([3 , n]); #Position of point P
PZ=PP #position of point z
dp=0 #distance gone by P
dz=0 #distance gone by z
VP=PP # V of P
VZ=PP
AP=PP #acceleration of P
AZ=PP #  .. . . .. . . Z
AP_n=PP
AP_t=PP
AZ_n=PP
AZ_t=PP
KC1p=np.zeros([3 , n]) # 1st order KC for P
KC2p=KC1p
rho=np.zeros([1 , n]) # R of curv. for z
C=np.zeros([3 , n]) # center of curv. for Z
PP[: , 0]= (R[0]-R11)*np.array([1 , 0 , 0])+R[5]*np.array([cosd(post[0 , 3]) , sind(post[0 , 3]) , 0])+R[4]/2*np.array([cosd(post[0 , 2]) ,sind(post[0 , 2])  , 0]);
PZ[: , 0]= (R[0]-R11)*np.array([1 , 0 , 0])+R[5]*np.array([cosd(post[0 , 3]), sind(post[0 , 3]) , 0])+R[4]/4*np.array([cosd(post[0 , 2]) ,sind(post[0 , 2])  , 0 ]); 
for i in range(1 , n):
    PP[: , i]= (R[0]-R11)*np.array([1 , 0 , 0])+R[5]*np.array([cosd(post[i , 3]) , sind(post[i , 3]) , 0])+R[4]/2*np.array([cosd(post[i , 2]) ,sind(post[i , 2])  , 0]);
    PZ[: , i]= (R[0]-R11)*np.array([1 , 0 , 0])+R[5]*np.array([cosd(post[i , 3]), sind(post[i , 3]) , 0])+R[4]/4*np.array([cosd(post[i , 2]) ,sind(post[i , 2])  , 0 ]); 
    dp=dp + np.linalg.norm((PP[: , i]- PP[: , i-1]));
    dz=dz + np.linalg.norm((PZ[: , i] - PZ[: , i-1]));
elap=time.time()-t