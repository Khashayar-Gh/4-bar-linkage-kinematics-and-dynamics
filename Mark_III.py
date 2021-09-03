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
#import sage 
R = np.array([203.2 , 57.15 , 184.15 , 177.8 , 50.8 , 127])
R44=127
R11=101.6
''
def cosd(theta):
    return np.cos(theta*(math.pi)/180)
def sind(theta):
    return np.sin(theta*(math.pi)/180)
def fun_I(tet2,tet):
    R = np.array([203.2 , 57.15 , 184.15 , 177.8 , 50.8 , 127])
    R44=127
    R11=101.6   
    tet=np.array([np.NAN ,tet2 ,  tet[0] , tet[1] , tet[2] , tet[3]])
    [F0 , F1]=R[1]*np.array([cosd(tet[1]) , sind(tet[1])]) + R[2]*np.array([cosd(tet[2]) , sind(tet[2])]) - R[3]*np.array([cosd(tet[3]) , sind(tet[3])])-R[0]*np.array([1 , 0]) 
    [F2 ,F3]=R[5]*np.array([cosd(tet[5]) ,sind(tet[5]) ]) + R[4]*np.array([cosd(tet[4]) , sind(tet[4])]) - R44*np.array([cosd(tet[3]) , sind(tet[3])])-R11*np.array([1 , 0])
    return F0 , F1 , F2 , F3
tet2=np.arange(0 , 360 , .1)
n=len(tet2)
init_guess=np.zeros([n+1,4]); #hads avaliye dar har marhale v javbe marhaleye ghabl
init_guess[0 , :]=np.array([64 , 110 , 136 ,42 ])
#options = optimset('Display','off');
for i in range(n):
    fun=lambda tet: fun_I( tet2[i] ,tet);
    #s=scipy.optimize.minimize(fun , init_guess[i , :] );
    #s=scipy.optimize.root(fun , init_guess[i , :] , tol=1e-12);
    s=scipy.optimize.fsolve(fun , init_guess[i , :] );
    #s=sage.numerical.optimize.minimize(fun , init_guess[i , :] )
    init_guess[i+1 , :]=s;
    #init_guess[i+1 , :]=s;
'''    
init_guess(1 , :)=[];
post=[zeros(n , 1) , tet2' ,init_guess]; #postures (tet1 , tet2 , tet3 , tet4 , tet5 , tet6)
% post=mod(post , 360); %tabdile zavaye be adade beyne 0 , 360
plot(tet2 , post(: , 3) , tet2 , post(: , 4) , tet2 , post(: , 5) , tet2 , post(: , 6))
legend('tet3' , 'tet4' , 'tet5' , 'tet6')
xlabel('tet2')
'''
elap=time.time()-t