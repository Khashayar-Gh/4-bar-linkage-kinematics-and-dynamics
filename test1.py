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
R = np.array([203.2 , 57.15 , 184.15 , 177.8 , 50.8 , 127])
R44=127
R11=101.6
# import useful libs
 
 
def newton_raphson(f, x_guess=None, max_num_iter= 100, tolerance=1e-4, alpha=1.0, print_info=True ):
    '''
    Author: Christian Howard
    Function for representing a Newton-Raphson iteration for multidimensional systems of equations
    :param f: function class that must define the following methods:
        - numDims(): Method that returns an integer number of variables in the system of equations		- __call__(np.ndarray): Method to make this class act like a function operating on some input x
    :param x_guess: an initial guess for the Newton-Raphson iteration
    :param max_num_iter: a maximum number of iterations that will be taken
    :param tolerance: a tolerance that will stop the sequence once the error drops below it
    :param alpha: A coefficient that can tune the Newton-Raphson stepsize. Recommend setting alpha <= 1.
    :return: A tuple with the root estimate, final error for the root, and the number of iterations it took
    '''
    # set the initial guess
    if x_guess is None:
        x_guess = np.random.rand(f.numDims())
    x = x_guess
 
    # compute function value at initial guess6
    fx = f(x)
 
    # define the initial value for the error and the starting iteration count
    err = np.linalg.norm(fx)
    iter = 0
 
    if print_info:
        print("Iteration {0}: Error of {1} with an estimate of {2}".format(iter,err,x))
 
    # perform the Newton-Raphson iteration algo
    while err > tolerance and iter < max_num_iter:
 
        # perform newton step
        x = x - alpha*np.linalg.solve(f.getJacobian(x),fx)
 
        # update the function value at the new root estimate
        fx = f(x)
 
        # compute the current root error
        err = np.linalg.norm(fx)
 
        # update the iteration counter
        iter = iter + 1
 
        # print useful message
        if print_info:
            print("Iteration {0}: Error of {1} with an estimate of {2}".format(iter, err, x))
 
    # return the root estimate, 2-norm error of the estimate, and iteration count we ended at
    return (x, err, iter)
def cosd(theta):
    return np.cos(theta*(math.pi)/180)
def sind(theta):
    return np.sin(theta*(math.pi)/180)
def fun_I(tet2,tet):
    R = np.array([203.2 , 57.15 , 184.15 , 177.8 , 50.8 , 127])
    R44=127
    R11=101.6   
    tet=np.array([np.NaN ,tet2 ,  tet[0] , tet[1] , tet[2] , tet[3]])
    F=np.array([0, 0 ,0, 0])
    [F[0] , F[1]]=R[1]*np.array([cosd(tet[1]) , sind(tet[1])]) + R[2]*np.array([cosd(tet[2]) , sind(tet[2])]) - R[3]*np.array([cosd(tet[3]) , sind(tet[3])])-R[0]*np.array([1 , 0]) 
    [F[2] ,F[3]]=R[5]*np.array([cosd(tet[5]) ,sind(tet[5]) ]) + R[4]*np.array([cosd(tet[4]) , sind(tet[4])]) - R44*np.array([cosd(tet[3]) , sind(tet[3])])-R11*np.array([1 , 0])
    return F
tet2=np.arange(0 , 360 , .01)
n=len(tet2)
init_guess=np.zeros([n+1,4]); #hads avaliye dar har marhale v javbe marhaleye ghabl
init_guess[0 , :]=np.array([64 , 110 , 136 ,42 ])
#options = optimset('Display','off');
for i in range(n+1):
    fun=lambda tet: fun_I( tet2[i] ,tet);
    #s=scipy.optimize.root(fun , init_guess[i , :] , method='df-sane');
    [s , err , itr]=newton_raphson(fun, x_guess=init_guess[i , :])
    init_guess[i+1 , :]=s;
'''    
init_guess(1 , :)=[];
post=[zeros(n , 1) , tet2' ,init_guess]; #postures (tet1 , tet2 , tet3 , tet4 , tet5 , tet6)
% post=mod(post , 360); %tabdile zavaye be adade beyne 0 , 360
plot(tet2 , post(: , 3) , tet2 , post(: , 4) , tet2 , post(: , 5) , tet2 , post(: , 6))
legend('tet3' , 'tet4' , 'tet5' , 'tet6')
xlabel('tet2')
'''
elap=time.time()-t