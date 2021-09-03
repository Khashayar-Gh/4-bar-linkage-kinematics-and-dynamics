import numpy as np
def cosd(theta):
    return np.cos(theta*(np.pi/180))
def sind(theta):
    return np.sin(theta*(np.pi/180))
n=2
PP=np.zeros([3 , n])
PZ=np.zeros([3 , n])
dp , dz=0 ,  0
post=np.array([[63.942	,111.497	,139.935	,42.2969],[63.9029,	111.458,	139.843	,42.2781]])
R = np.array([203.2 , 57.15 , 184.15 , 177.8 , 50.8 , 127])
R44=127
R11=101.6
PP[: , 0]= (R[0]-R11)*np.array([1 , 0 , 0])+R[5]*np.array([cosd(post[0 , 3]) , sind(post[0 , 3]) , 0])+R[4]/2*np.array([cosd(post[0 , 2]) ,sind(post[0 , 2])  , 0])
PZ[: , 0]= (R[0]-R11)*np.array([1 , 0 , 0])+R[5]*np.array([cosd(post[0 , 3]), sind(post[0 , 3]) , 0])+R[4]/4*np.array([cosd(post[0 , 2]) ,sind(post[0 , 2])  , 0 ])
d=0
for i in range(1 , n):

    PP[: , i]= (R[0]-R11)*np.array([1 , 0 , 0]) + R[5]*np.array([cosd(post[i , 3]) , sind(post[i , 3]) , 0]) + (R[4]/2)*np.array([cosd(post[i , 2]) ,sind(post[i , 2])  , 0])
    PZ[: , i]= (R[0]-R11)*np.array([1 , 0 , 0]) + R[5]*np.array([cosd(post[i , 3]) , sind(post[i , 3]) , 0]) + (R[4]/4)*np.array([cosd(post[i , 2]) ,sind(post[i , 2])  , 0])
    '''
    PP[: , i]= (R[0]-R11)*np.array([1 , 0 , 0])+R[5]*np.array([cosd(post[i , 3]) , sind(post[i , 3]) , 0])+R[4]/2*np.array([cosd(post[i , 2]) ,sind(post[i , 2])  , 0]);
    PZ[: , i]= (R[0]-R11)*np.array([1 , 0 , 0])+R[5]*np.array([cosd(post[i , 3]), sind(post[i , 3]) , 0])+R[4]/4*np.array([cosd(post[i, 2]) ,sind(post[i , 2])  , 0 ]); 
    '''
    dp=dp + np.linalg.norm((PP[: , i]- PP[: , i-1]));
    dz=dz + np.linalg.norm((PZ[: , i] - PZ[: , i-1]));
