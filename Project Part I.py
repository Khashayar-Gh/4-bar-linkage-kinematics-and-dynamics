#%% importing requiring moduls
import time
t=time.time()
import numpy as np
import scipy.optimize
from tkinter import ttk
import tkinter as tk

from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.figure import Figure
import matplotlib.animation as animation 
#%% functions
def cosd(theta):
    return np.cos(theta*(np.pi/180))
def sind(theta):
    return np.sin(theta*(np.pi/180))
def fun_I(tet2,tet):
    R = np.array([203.2 , 57.15 , 184.15 , 177.8 , 50.8 , 127])
    R44=127
    R11=101.6   
    tet=np.array([np.NAN ,tet2 ,  tet[0] , tet[1] , tet[2] , tet[3]])
    [F0 , F1]=R[1]*np.array([cosd(tet[1]) , sind(tet[1])]) + R[2]*np.array([cosd(tet[2]) , sind(tet[2])]) - R[3]*np.array([cosd(tet[3]) , sind(tet[3])])-R[0]*np.array([1 , 0]) 
    [F2 ,F3]=R[5]*np.array([cosd(tet[5]) ,sind(tet[5]) ]) + R[4]*np.array([cosd(tet[4]) , sind(tet[4])]) - R44*np.array([cosd(tet[3]) , sind(tet[3])])-R11*np.array([1 , 0])
    return F0 , F1 , F2 , F3
def AV2(tet): #angular velocity2
    if tet>=0 and tet<50:
        a = np.sqrt(-5*np.pi*AA2(tet)/9*(1-tet/50))
    elif tet >= 50 and tet<150:
        a = np.sqrt(10*np.pi*AA2(tet)/9*(tet - 50)/100)
    elif tet>=150 and tet <250:
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

tet2=np.arange(0 , 360 , 0.1)
n=len(tet2)
init_guess=np.zeros([n+1,4]) #hads avaliye dar har marhale v javbe marhaleye ghabl
init_guess[0 , :]=np.array([64 , 110 , 136 ,42 ])
#%%optimizing and finding theta
for i in range(0 ,n):
    fun=lambda tet: fun_I( tet2[i] ,tet)
    s=scipy.optimize.fsolve(fun , init_guess[i , :] );
    init_guess[i+1 , :]=s
post= init_guess[1:,:]
del init_guess, i, s
#%% velocity analysis
R = np.array([203.2 , 57.15 , 184.15 , 177.8 , 50.8 , 127])
R44=127
R11=101.6
KC1 = np.zeros([n,4]) #1st order kinematic coefficient
KC2 = np.zeros([n,4])
AV = np.zeros([n,6])

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
#%%Accelartion Analysis (Ax=B)
AA=np.zeros([n , 6])
DET2=np.zeros([n , 1])
for i in range(n):
    A[0 , :] = np.array([R[2]*-sind(post[i,0]),R[3]*sind(post[i,1]), 0 ,0 ])
    A[1 , :] = np.array([R[2]*cosd(post[i,0]), -R[3]*cosd(post[i,1]) , 0 , 0])
    A[2 , :] = np.array([0 , R44*sind(post[i , 1]),-R[4]*sind(post[i , 2]), R[5]*-sind(post[i,3])])
    A[3 , :] = np.array([0 , -R44*cosd(post[i,1]) , R[4]*cosd(post[i,2]) , R[5]*cosd(post[i , 3]) ])
    B=[[-R[1]*-cosd(tet2[i])-R[2]*(KC1[i , 0])**2*-cosd(post[i , 0])+R[3]*(KC1[i , 1])**2*-cosd(post[i , 1]) ],\
       [-R[1]*-sind(tet2[i])-R[2]*(KC1[i , 0])**2*-sind(post[i , 0])+R[3]*(KC1[i , 1])**2*-sind(post[i , 1])],\
       [-R[5]*(KC1[i , 3])**2*-cosd(post[i , 3])-R[4]*(KC1[i , 2])**2*-cosd(post[i , 2]) + R44*(KC1[i , 1])**2*-cosd(post[i , 1])],\
       [-R[5]*(KC1[i , 3])**2*-sind(post[i , 3])-R[4]*(KC1[i , 2])**2*-sind(post[i , 2]) + R44*(KC1[i , 1])**2*-sind(post[i , 1])]];
    x = np.linalg.solve(A,B)
    KC2[i] = [x[0] , x[1] , x[2] , x[3]]
    AA[i , 1] = AA2(tet2[i])
    DET2[i] = np.linalg.det(A)
AA[: , 2]=KC2[: , 0]*(AV[: , 1]**2)+KC1[: , 0]*AA[: , 1]; #AA of link 3
AA[: , 3]=KC2[: , 1]*(AV[: , 1]**2)+KC1[: , 1]*AA[: , 1];
AA[: , 4]=KC2[: , 2]*(AV[: , 1]**2)+KC1[: , 2]*AA[: , 1];
AA[: , 5]=KC2[: , 3]*(AV[: , 1]**2)+KC1[: , 3]*AA[: , 1];
#%% calculations of points P , Z
PP=np.zeros([3 , n]) #Position of point P
PZ=np.zeros([3 , n]) #position of point z
dp=0 #distance gone by P
dz=0 #distance gone by z
VP=np.zeros([3 , n]) # V of P
VZ=np.zeros([3 , n])
AP=np.zeros([3 , n]) #acceleration of P
AZ=np.zeros([3 , n]) #  .. . . .. . . Z
AP_n=np.zeros([3 , n])
AP_t=np.zeros([3 , n])
AZ_n=np.zeros([3 , n]) 
AZ_t=np.zeros([3 , n])
P_t=np.zeros([3 , n])
P_n=np.zeros([3 , n])
Z_t=np.zeros([3 , n])
Z_n=np.zeros([3 , n])
KC1p=np.zeros([3 , n]) # 1st order KC for P
KC2p=np.zeros([3 , n])
rho=np.zeros([1 , n]) # R of curv. for z
O=np.zeros([3 , n]) # center of curv. for Z
PP[: , 0]= (R[0]-R11)*np.array([1 , 0 , 0])+R[5]*np.array([cosd(post[0 , 3]) , sind(post[0 , 3]) , 0])+R[4]/2*np.array([cosd(post[0 , 2]) ,sind(post[0 , 2])  , 0])
PZ[: , 0]= (R[0]-R11)*np.array([1 , 0 , 0])+R[5]*np.array([cosd(post[0 , 3]), sind(post[0 , 3]) , 0])+R[4]/4*np.array([cosd(post[0 , 2]) ,sind(post[0 , 2])  , 0 ])
for i in range(1 , n):

    PP[: , i]= (R[0]-R11)*np.array([1 , 0 , 0]) + R[5]*np.array([cosd(post[i , 3]) , sind(post[i , 3]) , 0]) + (R[4]/2)*np.array([cosd(post[i , 2]) ,sind(post[i , 2])  , 0])
    PZ[: , i]= (R[0]-R11)*np.array([1 , 0 , 0]) + R[5]*np.array([cosd(post[i , 3]) , sind(post[i , 3]) , 0]) + (R[4]/4)*np.array([cosd(post[i , 2]) ,sind(post[i , 2])  , 0])
    dp=dp + np.linalg.norm((PP[: , i]- PP[: , i-1]))
    dz=dz + np.linalg.norm((PZ[: , i] - PZ[: , i-1]))
for i in range(n):
    V_D= np.cross(AV[i , 5]*np.array([0, 0, 1]) , R[5]*np.array([cosd(post[i , 3]) , sind(post[i , 3]) , 0]))
    V_PD=np.cross(AV[i , 4]*np.array([0, 0, 1]) , (R[4]/2)*np.array([cosd(post[i , 2]) , sind(post[i , 2]) , 0])) # V of P|D
    V_ZD=np.cross(AV[i , 4]*np.array([0, 0, 1]) , R[4]/4*np.array([cosd(post[i , 2]) , sind(post[i , 2]) , 0]))# V of Z|D
    VP[:,i] = V_D + V_PD
    VZ[: , i] = V_D + V_ZD
    A_D=np.cross(AA[i , 5]*np.array([0, 0, 1]) , R[5]*np.array([cosd(post[i , 3]) , sind(post[i , 3]) , 0])) + np.cross(AV[i , 5]*np.array([0, 0, 1]) , np.cross(AV[i , 5]*np.array([0, 0, 1]) , R[5]*np.array([cosd(post[i , 3]) , sind(post[i , 3]) , 0])))
    A_PD=np.cross(AA[i , 4]*np.array([0, 0, 1]) , R[4]/2*np.array([cosd(post[i , 2]) , sind(post[i , 2]) , 0])) + np.cross(AV[i , 4]*np.array([0, 0, 1]) , np.cross(AV[i , 4]*np.array([0, 0, 1]) , R[4]/2*np.array([cosd(post[i , 2]) , sind(post[i , 2]) , 0])))
    A_ZD=np.cross(AA[i , 4]*np.array([0, 0, 1]) , R[4]/4*np.array([cosd(post[i , 2]) , sind(post[i , 2]) , 0])) + np.cross(AV[i , 4]*np.array([0, 0, 1]) , np.cross(AV[i , 4]*np.array([0, 0, 1]) , R[4]/4*np.array([cosd(post[i , 2]) , sind(post[i , 2]) , 0])))
    AP[:,i] = A_D + A_PD
    AZ[:,i] = A_D + A_ZD
    AP_t[:,i] = VP[: , i]*np.dot(AP[: , i] , VP[: , i])/(np.linalg.norm(VP[: , i]))**2
    AP_n[: , i]=AP[: , i]-AP_t[: , i]
    AZ_t[: , i]=np.dot(AZ[: , i] , VZ[: , i])*VZ[: , i]/(np.linalg.norm(VZ[: , i]))**2
    AZ_n[: , i]=AZ[: , i]-AZ_t[: , i]
    P_t[: , i]=AP_t[: , i]/np.linalg.norm(AP_t[: , i])
    P_n[: , i]=AP_n[: , i]/np.linalg.norm(AP_n[: , i])
    Z_t[: , i]=AZ_t[: , i]/np.linalg.norm(AZ_t[: , i])
    Z_n[: , i]=AZ_n[: , i]/np.linalg.norm(AZ_n[: , i])
    rho[0,i]=(np.linalg.norm(VZ[: , i]))**3/np.linalg.norm(np.cross(VZ[: , i] , AZ[: , i] ))
    O[: , i]=PZ[: , i]+rho[0,i]*AZ_n[: , i]/np.linalg.norm(AZ_n[: , i])
    KC1p[: ,i] = R[5]*KC1[i , 3]*np.array([-sind(post[i , 3]) , cosd(post[i , 3]) , 0]) + R[4]/2*KC1[i , 2]*np.array([-sind(post[i , 2]) , cosd(post[i , 2]) , 0])
    KC2p[:,i] = R[5]*KC2[i ,3]*np.array([-sind(post[i , 3]) , cosd(post[i , 3]) , 0]) + R[4]/2*KC2[i , 2]*np.array([-sind(post[i , 2]) , cosd(post[i , 2]) , 0])+R[5]*((KC1[i , 3])**2)*np.array([-cosd(post[i , 3]) , -sind(post[i , 3]) , 0])+R[4]/2*((KC1[i , 2])**2)*np.array([-cosd(post[i , 2]) , -sind(post[i , 2]) , 0])

#%% preallocatings related to finding Is
I=np.zeros([6 , 6])
Ix=np.zeros([n ,6 , 6])
Iy=np.zeros([n , 6 , 6])
o6=(R[0]-R11)*np.array([1, 0])
o4=R[0]*np.array([1, 0])

xx=np.zeros([2 , 2])
yy=np.zeros([2 , 2])
t=np.array([0 , 0])
AAA=np.zeros([2 , 2])
BBB=np.zeros([2 , 1])
#%% GUI
def plot_tet():
    listBox.place(relx=.5 , rely=.5,relheight=0 , relwidth=0)


    try: 
        ax1.clear()
    except AttributeError: 
        pass     
    ax1.plot(tet2 , post[:,0],label = '\u03B8\u2083')
    ax1.plot(tet2 , post[:,1] , label = '\u03B8\u2084')
    ax1.plot(tet2 , post[:,2] , label = '\u03B8\u2085')
    ax1.plot(tet2 , post[:,3] , label = '\u03B8\u2086')
    ax1.set_xlabel('\u03B8\u2082(degree)')
    ax1.set_ylabel('\u03B8(degree)')
    ax1.legend(fontsize=15)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    toolbar.update()
    
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
def plot_AV():
    listBox.place(relx=.5 , rely=.5,relheight=0 , relwidth=0)


    try: 

        ax1.clear()

    except AttributeError: 
        pass     
    ax1.plot(tet2 , AV[: , 2],label = '\u03C9\u2083')
    ax1.plot(tet2 , AV[: , 3] , label = '\u03C9\u2084')
    ax1.plot(tet2 , AV[: , 4], label = '\u03C9\u2085')
    ax1.plot(tet2 , AV[: , 5] , label = '\u03C9\u2086')
    ax1.set_xlabel('\u03B8\u2082(degree)')
    ax1.set_ylabel('\u03C9(rad.s\u207B\u00B9)')
    ax1.legend(fontsize=15)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    toolbar.update()
    
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
def plot_KC1():
    listBox.place(relx=.5 , rely=.5,relheight=0 , relwidth=0)
    try: 
        ax1.clear()
        
    except:
        pass   
    ax1.plot(tet2 , KC1[: , 0],label = '\u03B8\u2032\u2083')
    ax1.plot(tet2 , KC1[: , 1] , label = '\u03B8\u2032\u2084')
    ax1.plot(tet2 , KC1[: , 2], label = '\u03B8\u2032\u2085')
    ax1.plot(tet2 , KC1[: , 3] , label = '\u03B8\u2032\u2086')
    ax1.set_xlabel('\u03B8\u2082(degree)')
    ax1.set_ylabel('\u03B8\u2032')
    ax1.legend(fontsize=15)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    

    toolbar.update()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
def plot_KC2():
    listBox.place(relx=.5 , rely=.5,relheight=0 , relwidth=0)
    try: 
        ax1.clear()
        
    except:
        pass   
    ax1.plot(tet2 , KC2[: , 0],label = '\u03B8\u2032\u2032\u2083')
    ax1.plot(tet2 , KC2[: , 1] , label = '\u03B8\u2032\u2032\u2084')
    ax1.plot(tet2 , KC2[: , 2], label = '\u03B8\u2032\u2032\u2085')
    ax1.plot(tet2 , KC2[: , 3] , label = '\u03B8\u2032\u2032\u2086')
    ax1.set_xlabel('\u03B8\u2082(degree)')
    ax1.set_ylabel('\u03B8\u2032\u2032')
    ax1.legend(fontsize=15)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    toolbar.update()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
def plot_AA():
    listBox.place(relx=.5 , rely=.5,relheight=0 , relwidth=0)
    try: 
        ax1.clear()
        
    except:
        pass   
    ax1.plot(tet2 , AA[: , 2],label = '\u03B1\u2083')
    ax1.plot(tet2 , AA[: , 3] , label = '\u03B1\u2084')
    ax1.plot(tet2 , AA[: , 4], label = '\u03B1\u2085')
    ax1.plot(tet2 , AA[: , 5] , label = '\u03B1\u2086')
    ax1.set_xlabel('\u03B8\u2082(degree)')
    ax1.set_ylabel('\u03B1(rad.s\u207B\u00B2)')
    ax1.legend(fontsize=15)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    toolbar.update()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
def plot_PathP():
    listBox.place(relx=.5 , rely=.5,relheight=0 , relwidth=0)
    try: 
        ax1.clear()
        
    except:
        pass   
    ax1.plot(PP[0 , :] , PP[1 , :],label = 'Path of Point P')
    ax1.set_xlabel('X(mm)')
    ax1.set_ylabel('Y(mm)')
    ax1.legend(fontsize=15)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    
    toolbar.update()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
def plot_PathZ():
    listBox.place(relx=.5 , rely=.5,relheight=0 , relwidth=0)
    try: 
        ax1.clear()
        
    except:
        pass   
    ax1.plot(PZ[0 , :] , PZ[1 , :],label = 'Path of Point Z')
    ax1.set_xlabel('X(mm)')
    ax1.set_ylabel('Y(mm)')
    ax1.legend(fontsize=15)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    
    toolbar.update()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
def plot_VP():
    listBox.place(relx=.5 , rely=.5,relheight=0 , relwidth=0)
    try: 
        ax1.clear()
        
    except:
        pass   
    ax1.plot(tet2 , VP[0 , :],label = 'V\u2093')
    ax1.plot(tet2 , VP[1 , :] , label = 'Vᵧ')
    ax1.set_xlabel('\u03B8\u2082(degree)')
    ax1.set_ylabel('V(mm.s\u207B\u00B9)')
    ax1.legend(fontsize=15)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    
    toolbar.update()
def plot_VZ():
    listBox.place(relx=.5 , rely=.5,relheight=0 , relwidth=0)
    try: 
        ax1.clear()
        fig.get_children().forget()
        
    except:
        pass   
    ax1.plot(tet2 , VZ[0 , :],label = 'V\u2093')
    ax1.plot(tet2 , VZ[1 , :] , label = 'Vᵧ')
    ax1.set_xlabel('\u03B8\u2082(degree)')
    ax1.set_ylabel('V(mm.s\u207B\u00B9)')
    ax1.legend(fontsize=15)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    
    toolbar.update()
def plot_AZ():
    listBox.place(relx=.5, rely=.5,relheight=0 , relwidth=0)
    try: 
        ax1.clear()
        fig.get_children().forget()
        
    except:
        pass   
    ax1.plot(tet2 , AZ[0 , :],label = 'a\u2093')
    ax1.plot(tet2 , AZ[1 , :] , label = 'aᵧ')
    ax1.set_xlabel('\u03B8\u2082(degree)')
    ax1.set_ylabel('a(mm.s\u207B\u00B2)')

    ax1.legend(fontsize=15)
    canvas.draw()

    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    
    #toolbar = NavigationToolbar2Tk(canvas, results)
    toolbar.update()
def plot_AP():
    listBox.place(relx=.5, rely=.5,relheight=0 , relwidth=0)

    try: 

        ax1.clear()
        fig.get_children().forget()
        
    except:
        pass   
    ax1.plot(tet2 , AP[0 , :],label = 'a\u2093')
    ax1.plot(tet2 , AP[1 , :] , label = 'aᵧ')
    ax1.set_xlabel('\u03B8\u2082(degree)')
    ax1.set_ylabel('a(mm.s\u207B\u00B2)')
    ax1.legend(fontsize=15)
    canvas.draw()

    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    toolbar.update()
def plot_mechanism():
    listBox.place(relx=.5 , rely=.5,relheight=0 , relwidth=0)
    try: 
        

        ax1.clear()
    
    except:
        pass  
    ax1.set_xlim([-60, 240])
    ax1.set_ylim([-70, 200])
    ax1.set_xlabel('X(mm)')
    ax1.set_ylabel('Y(mm)')
    line, = ax1.plot([], [], lw=2.5 , color = '#ffcc00' ,label = 'Position of point P' , linestyle = ':') 
    line2, = ax1.plot([], [], lw=2.5 ,color = '#ff0055' , label = 'Position of point Z' , linestyle = ':')
    line3, = ax1.plot([], [], lw=3 ,color = 'black')
    line4, = ax1.plot([], [], lw=3 ,color = 'black')
    line5, = ax1.plot([], [], lw=3 ,color = 'black')
    line6, = ax1.plot([], [], lw=3 ,color = 'black')
    line7, = ax1.plot([], [], lw=3 ,color = 'black')
    
    # initialization function 
    
    
    # lists to store x and y axis points 
    xdata, ydata = [], [] 
    x1data, y1data = [], []
    # animation function 
    o6=(R[0]-R11)*np.array([1, 0])
    o4=R[0]*np.array([1, 0])

    def init(): 
    	# creating an empty plot/frame 
    	line.set_data([], []) 
    	return line, 
    
    def animate(i): 
        x = PP[0,i] 
        y = PP[1,i]
        x1 = PZ[0,i] 
        y1= PZ[1,i]
        Ax= [0,R[1]*cosd(tet2[i])]
        Ay= [0 ,R[1]*sind(tet2[i])]
        Bx= [R[1]*cosd(tet2[i]) , R[1]*cosd(tet2[i])+R[2]*cosd(post[i,0])]
        By= [R[1]*sind(tet2[i]) , R[1]*sind(tet2[i]) + R[2]*sind(post[i,0])]
        Cx= [o4[0] ,Bx[1]]
        Cy= [o4[1] ,By[1]]
        Dx= [o6[0] ,o6[0] + R[5]*cosd(post[i,3])] 
        Dy= [o6[1] ,o6[1] + R[5]*sind(post[i,3])]
        Ex= [o4[0] + R44*cosd(post[i,1]) , Dx[1] ]
        Ey= [o4[1] + R44*sind(post[i,1]) ,Dy[1]]
    	# appending new points to x, y axes points list 
        xdata.append(x) 
        ydata.append(y) 
        x1data.append(x1) 
        y1data.append(y1)
        line.set_data(xdata, ydata) 
        line2.set_data(x1data, y1data)
        line3.set_data(Ax, Ay)
        line4.set_data(Bx, By)
        line5.set_data(Cx, Cy)
        line6.set_data(Dx, Dy)
        line7.set_data(Ex, Ey) 
        #time.sleep(t[i])
        return line, line2,line3, line4, line5,line6, line7
    

    animation.FuncAnimation(fig, animate, init_func=init, frames=3600, interval=1, blit=True,repeat = False)
    ax1.legend(fontsize=15)

    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    toolbar.update()
        
def plot_Cent():
    listBox.place(relx=.5 , rely=.5,relheight=0 , relwidth=0)
    try: 
        

        ax1.clear()
    
    except:
        pass  
    ax1.set_xlim([-60, 350])
    ax1.set_ylim([-70, 500])
    ax1.set_xlabel('X(mm)')
    ax1.set_ylabel('Y(mm)')
    #ax1.axes(xlim=(-60, 240), ylim=(-70, 200))

    line, = ax1.plot([], [], lw=2.5 , color = '#ffcc00' ,label = 'I₁₃' , linestyle = ':') 
    line2, = ax1.plot([], [], lw=2.5 ,color = '#ff0055' , label = 'I₁₅' , linestyle = ':')
    line3, = ax1.plot([], [], lw=3 ,color = 'black')
    line4, = ax1.plot([], [], lw=3 ,color = 'black')
    line5, = ax1.plot([], [], lw=3 ,color = 'black')
    line6, = ax1.plot([], [], lw=3 ,color = 'black')
    line7, = ax1.plot([], [], lw=3 ,color = 'black')
    
    # initialization function 
    
    
    # lists to store x and y axis points 
    xdata, ydata = [], [] 
    x1data, y1data = [], []
    # animation function 
    o6=(R[0]-R11)*np.array([1, 0])
    o4=R[0]*np.array([1, 0])

    def init(): 
    	# creating an empty plot/frame 
    	line.set_data([], []) 
    	return line, 
    
    def animate(i): 
        x = Ix[i ,0,2] 
        y = Iy[i ,0 , 2]
        x1 = Ix[i ,0,4]
        y1= Iy[i ,0 , 4]
        Ax= [0,R[1]*cosd(tet2[i])]
        Ay= [0 ,R[1]*sind(tet2[i])]
        Bx= [R[1]*cosd(tet2[i]) , R[1]*cosd(tet2[i])+R[2]*cosd(post[i,0])]
        By= [R[1]*sind(tet2[i]) , R[1]*sind(tet2[i]) + R[2]*sind(post[i,0])]
        Cx= [o4[0] ,Bx[1]]
        Cy= [o4[1] ,By[1]]
        Dx= [o6[0] ,o6[0] + R[5]*cosd(post[i,3])] 
        Dy= [o6[1] ,o6[1] + R[5]*sind(post[i,3])]
        Ex= [o4[0] + R44*cosd(post[i,1]) , Dx[1] ]
        Ey= [o4[1] + R44*sind(post[i,1]) ,Dy[1]]
    	# appending new points to x, y axes points list 
        xdata.append(x) 
        ydata.append(y) 
        x1data.append(x1) 
        y1data.append(y1)
        line.set_data(xdata, ydata) 
        line2.set_data(x1data, y1data)
        line3.set_data(Ax, Ay)
        line4.set_data(Bx, By)
        line5.set_data(Cx, Cy)
        line6.set_data(Dx, Dy)
        line7.set_data(Ex, Ey) 
        #time.sleep(t[i])
        return line, line2,line3, line4, line5,line6, line7
    

    animation.FuncAnimation(fig, animate, init_func=init, frames=3600, interval=1, blit=True,repeat = False)
    ax1.legend(fontsize=15)

    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    toolbar.update()
            
def plot_RZx():
    listBox.place(relx=.5, rely=.5,relheight=0 , relwidth=0)

    try: 

        ax1.clear()
        fig.get_children().forget()
        
    except:
        pass   
    ax1.plot(tet2 , O[0 , :] , label = 'Oₓ')
    ax1.set_xlabel('\u03B8\u2082(degree)')
    ax1.set_ylabel('X(mm)')

    ax1.legend(fontsize=15)
    canvas.draw()

    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    toolbar.update()
def plot_RZy():
    listBox.place(relx=.5, rely=.5,relheight=0 , relwidth=0)

    try: 

        ax1.clear()
        fig.get_children().forget()
        
    except:
        pass   
    ax1.plot(tet2 , O[1 , :] , label = 'Oᵧ')
    ax1.set_xlabel('\u03B8\u2082(degree)')
    ax1.set_ylabel('Y(mm)')

    ax1.legend(fontsize=15)
    canvas.draw()

    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    toolbar.update()
def plot_RZ():
    listBox.place(relx=.5, rely=.5,relheight=0 , relwidth=0)

    try: 
        ax1.clear()
        fig.get_children().forget()
        
    except:
        pass   
    ax1.plot(tet2 , rho[0 , :],label = '\u03C1')
    ax1.set_xlabel('\u03B8\u2082(degree)')
    ax1.set_ylabel('\u03C1(mm)')

    ax1.legend(fontsize=15)
    canvas.draw()

    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    toolbar.update()
def plot_DET():
    listBox.place(relx=.5, rely=.5,relheight=0 , relwidth=0)

    try: 

        ax1.clear()
        fig.get_children().forget()
        
    except:
        pass   
    ax1.plot(tet2 , DET1 ,label = 'det')
    ax1.set_ylabel('det(mm\u00B2)')
    ax1.set_xlabel('\u03B8\u2082(degree)')

    ax1.legend(fontsize=15)
    canvas.draw()

    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    toolbar.update()
def plot_KC1px():
    listBox.place(relx=.5, rely=.5,relheight=0 , relwidth=0)

    try: 

        ax1.clear()
        fig.get_children().forget()
        
    except:
        pass   
    ax1.plot(tet2 , KC1p[0 , :] ,label = 'x\u2032\u209a')
    ax1.set_ylabel('x\u2032\u209a((mm.rad\u207B\u00B9)')
    ax1.set_xlabel('\u03B8\u2082(degree)')

    ax1.legend(fontsize=15)
    canvas.draw()

    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    toolbar.update()
def plot_KC1py():
    listBox.place(relx=.5, rely=.5,relheight=0 , relwidth=0)

    try: 

        ax1.clear()
        fig.get_children().forget()
        
    except:
        pass   
    ax1.plot(tet2 , KC1p[1 , :] ,label = 'y\u2032\u209a')
    ax1.set_ylabel('y\u2032\u209a((mm.rad\u207B\u00B9)')
    ax1.set_xlabel('\u03B8\u2082(degree)')

    ax1.legend(fontsize=15)
    canvas.draw()

    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    toolbar.update()
def plot_KC2px():
    listBox.place(relx=.5, rely=.5,relheight=0 , relwidth=0)

    try: 

        ax1.clear()
        fig.get_children().forget()
        
    except:
        pass   
    ax1.plot(tet2 , KC2p[0 , :] ,label = 'x\u2032\u209a')
    ax1.set_ylabel('x\u2032\u2032\u209a(mm.rad\u207B\u00B2)')
    ax1.set_xlabel('\u03B8\u2082(degree)')

    ax1.legend(fontsize=15)
    canvas.draw()

    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    toolbar.update()
def plot_KC2py():
    listBox.place(relx=.5, rely=.5,relheight=0 , relwidth=0)

    try: 

        ax1.clear()
        fig.get_children().forget()
        
    except:
        pass   
    ax1.plot(tet2 , KC2p[1 , :] ,label = 'y\u2032\u209a')
    ax1.set_ylabel('y\u2032\u2032\u209a(mm.rad\u207B\u00B2)')
    ax1.set_xlabel('\u03B8\u2082(degree)')
    #plt.tick_params

    ax1.legend(fontsize=15)
    canvas.draw()

    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    toolbar.update()
def plot_I():
    listBox.place(relx=.5, rely=.5,relheight=0 , relwidth=0)

    try: 

        ax1.clear()
        fig.get_children().forget()
        
    except:
        pass   
    kj=ent.get()
    k=int(kj[0])-1
    j=int(kj[1])-1
    label1='I\u2093 (mm)'
    label2='Iᵧ (mm)'
    ax1.plot(tet2 , Ix[: , k , j] ,label = label1)
    ax1.plot(tet2 , Iy[: , k , j] ,label = label2)
    ax1.set_ylabel('I(mm)')
    ax1.set_xlabel('\u03B8\u2082(degree)')
    ax1.legend(fontsize=15)
    canvas.draw()

    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    toolbar.update()
def plot_PI():
    listBox.place(relx=.5, rely=.5,relheight=0 , relwidth=0)

    try: 

        ax1.clear()
        fig.get_children().forget()
        
    except:
        pass   
    kj=ent.get()
    k=int(kj[0])-1
    j=int(kj[1])-1
    ax1.plot(Ix[: , k , j] , Iy[: , k , j],label ='Path of I')
    ax1.set_ylabel('Y(mm)')
    ax1.set_xlabel('X(mm)')
    #plt.tick_params

    ax1.legend(fontsize=15)
    canvas.draw()

    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    toolbar.update()
def table_tet():
    try: 
        canvas.get_tk_widget().pack_forget()
        ax1.clear()
        
    except:
        pass   
    col = ( '\u03B8\u2082(degree)', '\u03B8\u2083(degree)','\u03B8\u2084(degree)','\u03B8\u2085(degree)','\u03B8\u2086(degree)')
    for i in range(5):
        listBox.heading(i, text=col[i])  
    #listBox.column(0 , width=130)
    listBox.delete(*listBox.get_children())
    listBox.place(relx=.255,rely=.035,relheight=.9, relwidth=.735)

    for i in range(0 , n , 50):
        listBox.insert("", "end", values=(tet2[i],post[i , 0], post[i , 1], post[i , 2] ,post[i , 3] ))
def table_AV():
    try: 
        canvas.get_tk_widget().pack_forget()
        ax1.clear()
        
    except:
        pass   
    col = ( '\u03B8\u2082(degree)', '\u03C9\u2083(rad.s\u207B\u00B9)','\u03C9\u2084(rad.s\u207B\u00B9)','\u03C9\u2085(rad.s\u207B\u00B9)','\u03C9\u2086(rad.s\u207B\u00B9)')
    for i in range(5):
        listBox.heading(i, text=col[i])  
    listBox.delete(*listBox.get_children())
    listBox.place(relx=.255,rely=.035,relheight=.9, relwidth=.735)

    for i in range(0 , n , 50):
        listBox.insert("", "end", values=(tet2[i],AV[i , 2], AV[i , 3], AV[i , 4] ,AV[i , 5] ))

def table_AA():
    
    try: 
        canvas.get_tk_widget().pack_forget()
        ax1.clear()
        
    except:
        pass   
    col = ( '\u03B8\u2082(degree)', '\u03B1\u2083(rad.s\u207B\u00B2)','\u03B1\u2084(rad.s\u207B\u00B2)','\u03B1\u2085(rad.s\u207B\u00B2)','\u03B1\u2086(rad.s\u207B\u00B2)')
    for i in range(5):
        listBox.heading(i, text=col[i])  
    listBox.delete(*listBox.get_children())
    listBox.place(relx=.255,rely=.035,relheight=.9, relwidth=.735)

    for i in range(0 , n , 50):
        listBox.insert("", "end", values=(tet2[i],AA[i , 2], AA[i , 3], AA[i , 4] ,AA[i , 5] ))
def table_KC1():
    
    try: 
        canvas.get_tk_widget().pack_forget()
        ax1.clear()
        
    except:
        pass 
    col = ( '\u03B8\u2082(degree)', '\u03B8\u2032\u2083','\u03B8\u2032\u2084','\u03B8\u2032\u2085','\u03B8\u2032\u2086')
    for i in range(5):
        listBox.heading(i, text=col[i])  

    listBox.delete(*listBox.get_children())
    listBox.place(relx=.255,rely=.035,relheight=.9, relwidth=.735)
    
    for i in range(0 , n , 50):
        listBox.insert("", "end", values=(tet2[i],KC1[i , 0], KC1[i , 1], KC1[i , 2] ,KC1[i , 3] ))

def table_KC2():
    
    try: 
        canvas.get_tk_widget().pack_forget()
        ax1.clear()
        
    except:
        pass   
    col = ( '\u03B8\u2082(degree)', '\u03B8\u2032\u2032\u2083','\u03B8\u2032\u2032\u2084','\u03B8\u2032\u2032\u2085','\u03B8\u2032\u2032\u2086')
    for i in range(5):
        listBox.heading(i, text=col[i])  
    listBox.delete(*listBox.get_children())
    listBox.place(relx=.255,rely=.035,relheight=.9, relwidth=.735)



    for i in range(0 , n , 50):
        listBox.insert("", "end", values=(tet2[i] , KC2[i , 0], KC2[i , 1], KC2[i , 2] ,KC2[i , 3] ))
def table_TN():
    
    try: 
        canvas.get_tk_widget().pack_forget()
        ax1.clear()
        
    except:
        pass   
    col = ( '\u03B8\u2082(degree)', 'T of Z','N of Z', 'T of P','N of P' )
    for i in range(5):
        listBox.heading(i, text=col[i])  
    listBox.delete(*listBox.get_children())
    listBox.place(relx=.255,rely=.035,relheight=.9, relwidth=.735)


    for i in range(0 , n , 50):
        listBox.insert("", "end", values=(tet2[i] , str(round(Z_t[0 , i] , 4)) + ' i'+' + '+str(round(Z_t[1 , i] , 4))+' j',\
                                          str(round(Z_n[0 , i] , 4))+ ' i'+' + '+str(round(Z_n[1 , i], 4))+' j',\
                                          str(round(P_t[0 , i ] , 4))+ ' i'+' + '+str(round(P_t[1 , i], 4))+' j',\
                                          str(round(P_n[0 , i] , 4))+ ' i'+' + '+str(round(P_n[1 , i], 4))+' j'))
def table_dPZ():
    
    try: 
        canvas.get_tk_widget().pack_forget()
        ax1.clear()
        
    except:
        pass   
    col = ( 'dP(mm)','dZ(mm)', '','' , '' )
    for i in range(5):
        listBox.heading(i, text=col[i])  
    listBox.delete(*listBox.get_children())
    listBox.place(relx=.255,rely=.035,relheight=.9, relwidth=.735)


    listBox.insert("", "end", values=(dp , dz))
def table_RZ():
    
    try: 
        canvas.get_tk_widget().pack_forget()
        ax1.clear()
        
    except:
        pass   
    col = ( '\u03B8\u2082(degree)', '\u03C1 of Z (mm)','O of Z (mm)', '','' )
    for i in range(5):
        listBox.heading(i, text=col[i])  
    listBox.delete(*listBox.get_children())
    listBox.place(relx=.255,rely=.035,relheight=.9, relwidth=.735)

    for i in range(0 , n , 50):
        listBox.insert("", "end", values=(tet2[i] , rho[0 , i] , str(round(O[0 , i] , 5)) +' i ' + '+' + str(round(O[1 , i] , 5))+ ' j' ))

def table_KC1px():
    
    try: 
        canvas.get_tk_widget().pack_forget()
        ax1.clear()
        
    except:
        pass 
    col = ( '\u03B8\u2082(degree)', 'x\u2032\u209a(mm.rad\u207B\u00B9)','','','')
    for i in range(5):
        listBox.heading(i, text=col[i])  

    listBox.delete(*listBox.get_children())
    listBox.place(relx=.255,rely=.035,relheight=.9, relwidth=.735)
    
    for i in range(0 , n , 50):
        listBox.insert("", "end", values=(tet2[i],KC1p[0 , i] ))
def table_KC1py():
    
    try: 
        canvas.get_tk_widget().pack_forget()
        ax1.clear()
        
    except:
        pass 
    col = ( '\u03B8\u2082(degree)', 'y\u2032\u209a(mm.rad\u207B\u00B9)','','','')
    for i in range(5):
        listBox.heading(i, text=col[i])  

    listBox.delete(*listBox.get_children())
    listBox.place(relx=.255,rely=.035,relheight=.9, relwidth=.735)
    
    for i in range(0 , n , 50):
        listBox.insert("", "end", values=(tet2[i],KC1p[1 , i] ))
def table_KC2px():
    
    try: 
        canvas.get_tk_widget().pack_forget()
        ax1.clear()
        
    except:
        pass 
    col = ( '\u03B8\u2082(degree)', 'x\u2032\u2032\u209a(mm.rad\u207B\u00B2)','','','')
    for i in range(5):
        listBox.heading(i, text=col[i])  

    listBox.delete(*listBox.get_children())
    listBox.place(relx=.255,rely=.035,relheight=.9, relwidth=.735)
    
    for i in range(0 , n , 50):
        listBox.insert("", "end", values=(tet2[i],KC2p[0 , i] ))
def table_KC2py():
    
    try: 
        canvas.get_tk_widget().pack_forget()
        ax1.clear()
        
    except:
        pass 
    col = ( '\u03B8\u2082(degree)', 'y\u2032\u2032\u209a(mm.rad\u207B\u00B2)','','','')
    for i in range(5):
        listBox.heading(i, text=col[i])  

    listBox.delete(*listBox.get_children())
    listBox.place(relx=.255,rely=.035,relheight=.9, relwidth=.735)
    
    for i in range(0 , n , 50):
        listBox.insert("", "end", values=(tet2[i],KC2p[1 , i] ))

def Cal_I():
    for i in range(n):
        I=np.zeros([6 , 6])
        I[0 ,  1] , I[1 ,  0] , Ix[i ,0 ,  1] , Ix[i ,1 ,  0] ,Iy[i ,0 ,  1],Iy[i ,1 ,  0]=1 , 1 , 0 , 0, 0 , 0
        I[1 ,  2] , I[2 ,  1]=1 , 1
        I[2 ,  3] , I[3 ,  2]=1 , 1
        I[3 ,  4] , I[4 ,  3]=1 , 1
        I[4 ,  5] , I[5 ,  4]=1 , 1
        I[5 ,  0] , I[0 ,  5], Ix[i ,5 ,  0] ,Ix[i ,0 ,  5]  ,Iy[i ,5 ,  0],Iy[i ,0 ,  5]=1 , 1, o6[0 ,] ,  o6[0 ], o6[ 1], o6[ 1]
        I[3 ,  0] , I[0 ,  3] ,Ix[i ,3 ,  0] ,Ix[i ,0 ,  3]  ,Iy[i ,3 ,  0],Iy[i ,0 ,  3]=1 , 1 , o4[0 ] ,  o4[0 ], o4[ 1], o4[ 1]
        for Thing in range(6):
            I[Thing , Thing]=1
        A=R[1]*np.array([cosd(tet2[i]) , sind(tet2[i])]);
        B=A+R[2]*np.array([cosd(post[i  , 0]) , sind(post[i ,0]) ]);
        C=o4+R44*np.array([cosd(post[i , 1]) , sind(post[i , 1])]);
        D=o6+R[5]*np.array([cosd(post[i , 3]) , sind(post[i , 3])]);
        Ix[i , 1 , 2] , Ix[i ,2 , 1]=A[0] , A[0]
        Iy[ i ,1 , 2] , Iy[i ,2 , 1]=A[1] , A[1]
        Ix[ i ,3 , 2] , Ix[i ,2 , 3]=B[0] , B[0]
        Iy[ i ,3 , 2] , Iy[i ,2 , 3]=B[1] , B[1]
        Ix[i , 3 , 4] , Ix[i ,4 , 3]=C[0] , C[0]
        Iy[i , 3 , 4] , Iy[i ,4 , 3]=C[1] , C[1]
        Ix[ i ,4 , 5] , Ix[i ,5 , 4]=D[0] , D[0]
        Iy[ i ,4 , 5] , Iy[i ,5 , 4]=D[1] , D[1]
    
        j=0
        MaxIter=0
        while np.any(I==0):
            MaxIter+=1
            if MaxIter==30:
                break
            
            if np.count_nonzero(I[: , j]==1)>2 and  np.count_nonzero(I[: , j]==1)<6:
                y=np.where(I[: , j]==1)
                y=y[0][:]
                y=np.delete(y ,np.where(y==j))
                for k in range (6):
                    if k in y:
                        continue
                    if np.count_nonzero(I[k , :]==1)>2 :
                        x=np.where(I[k , :]==1)
                        x=x[0][:]
                        x=np.delete(x , np.where(x==k))
                        l=0
                        for ii in range(len(y)):
                            if y[ii] in x:
                                t[l]=y[ii]
                                l=l+1
                                if l==2:
                                     break
                    try:
                        if l==2:
                            xx[0 , 0]=Ix[i ,k , t[0]]
                            xx[0 , 1]=Ix[i ,t[0] , j]
                            xx[1 , 0]=Ix[i ,k , t[1]]
                            xx[1 , 1]=Ix[i ,t[1] , j]
                            yy[0 , 0]=Iy[i ,k , t[0]]
                            yy[0 , 1]=Iy[i ,t[0] , j]
                            yy[1 , 0]=Iy[i ,k , t[1]]
                            yy[1 , 1]=Iy[i ,t[1] , j]
                            AAA=np.array([[(yy[0 , 1]-yy[0 , 0]) , -1*(xx[0 , 1]-xx[0 , 0])] ,[ yy[1 , 1]-yy[1 , 0], -1*(xx[1 , 1]-xx[1 , 0]) ]])
                            BBB[0]=xx[0 ,0]*(yy[0 , 1]-yy[0 , 0])-yy[0 , 0]*(xx[0 , 1]-xx[0 , 0])
                            BBB[1]=xx[1 ,0]*(yy[1 , 1]-yy[1 , 0])-yy[1 , 0]*(xx[1 , 1]-xx[1 , 0])
                            m=np.linalg.solve(AAA , BBB)
                            I[k , j] , I[j , k]= 1 , 1
                            Ix[i ,k , j] , Ix[i ,j , k]=m[0] , m[0]
                            Iy[i ,k , j] , Iy[i ,j , k]=m[1] , m[1]
    
                    except:
                        pass
                                
            if j==5:
                j=0
            else:
                j+=1
def table_I():
    
    try: 
        canvas.get_tk_widget().pack_forget()
        ax1.clear()
        
    except:
        pass   
    kj=ent.get()
    k=int(kj[0])-1
    j=int(kj[1])-1
    label1='I(mm)'
    col = ( '\u03B8\u2082(degree)',label1,'', '','' )
    for i in range(5):
        listBox.heading(i, text=col[i])  
    listBox.delete(*listBox.get_children())
    listBox.place(relx=.255,rely=.035,relheight=.9, relwidth=.735)


    for i in range(0 , n , 50):
        listBox.insert("", "end", values=(tet2[i] , str(round(Ix[i , k , j] , 4)) + ' i' + ' + ' + str(round(Iy[i , k ,j] , 4))+' j'))


elap=time.time()-t


root=tk.Tk()
root.title('Project Part I')
root.wm_iconbitmap('logo.ico')
root.geometry('1920x950')
root.configure(background='navy')

results=tk.LabelFrame(root, text="Results", padx=5, pady=5 , bg='navy' , fg='white')
results.place(relx=.25 , rely=0 , relheight=.95 , relwidth=.75)
fig = Figure(figsize=(5, 4), dpi=100)
ax1=fig.add_subplot(111)
ax1.set_xlabel('',fontsize=15)
ax1.set_ylabel('',fontsize=15)
canvas = FigureCanvasTkAgg(fig, master=results)
toolbar = NavigationToolbar2Tk(canvas, root)
toolbar.config(background='navy')
cols = ( '\u03B8\u2082', 'bar3','bar4','bar5','bar6')
listBox = ttk.Treeview(root, columns=cols, show='headings')
for col in cols:
    listBox.heading(col, text=col)   





PA_Label = tk.LabelFrame(root, text="Posture Analisys", padx=5, pady=5 , bg='slategray')
PA_Label.place(relx=0.035 , rely=0.0005)
VA_Label=tk.LabelFrame(root, text="Velocity Analisys", padx=5, pady=5 , bg='slategray')
VA_Label.place(relx=0.005 , rely=0.2)
AA_Label=tk.LabelFrame(root, text="Acceleration Analisys", padx=5, pady=5 , bg='slategray')
AA_Label.place(relx=0.12 , rely=0.2)
I_Label=tk.LabelFrame(root, text="Instant Center of Rotation", padx=5, pady=5 , bg='slategray')
I_Label.place(relx=0.05 , rely=0.45)
others_Label=tk.LabelFrame(root, text="Others", padx=5, pady=5 , bg='slategray')
others_Label.place(relx=0.035 , rely=0.65)

but_Ptet=tk.Button(PA_Label ,text='plot \u03B8s vs \u03B8\u2082', command=plot_tet  , padx=10, pady=5 )
but_Ptet.grid(column=0 , row=0)
but_Ttet=tk.Button(PA_Label ,text='table \u03B8s vs \u03B8\u2082', command=table_tet , padx=5, pady=5)
but_Ttet.grid(column=1 , row=0)
but_PathP=tk.Button(PA_Label ,text='draw path of Point P', command=plot_PathP  , padx=5, pady=5 )
but_PathP.grid(column=0 , row=1 , columnspan=1)
but_PathZ=tk.Button(PA_Label ,text='draw path of Point Z', command=plot_PathZ  , padx=5, pady=5 )
but_PathZ.grid(column=1 , row=1 , columnspan=1)
but_PathZ=tk.Button(PA_Label ,text='distance gone by points P , Z ', command=table_dPZ  , padx=5, pady=5 )
but_PathZ.grid(column=0 , row=2 , columnspan=2)


but_Pw=tk.Button(VA_Label ,text='plot \u03C9s vs \u03B8\u2082', command=plot_AV )
but_Pw.grid(column=0 , row=0)
but_Tw=tk.Button(VA_Label ,text='table \u03C9s vs \u03B8\u2082', command=table_AV )
but_Tw.grid(column=1 , row=0)
but_PKC1=tk.Button(VA_Label ,text='plot \u03B8\u2032s vs \u03B8\u2082', command=plot_KC1 )
but_PKC1.grid(column=0 , row=1)
but_TKC1=tk.Button(VA_Label ,text='table \u03B8\u2032s vs \u03B8\u2082', command=table_KC1 )
but_TKC1.grid(column=1 , row=1)
but_PVP=tk.Button(VA_Label ,text='plot V\u209a vs \u03B8\u2082', command=plot_VP )
but_PVP.grid(column=0 , row=2)
but_PVZ=tk.Button(VA_Label ,text='plot Vz vs \u03B8\u2082', command=plot_VZ )
but_PVZ.grid(column=1 , row=2)
but_PKC1px=tk.Button(VA_Label ,text='plot x\u2032\u209a vs \u03B8\u2082', command=plot_KC1px )
but_PKC1px.grid(column=0 , row=3)
but_PKC1py=tk.Button(VA_Label ,text='plot y\u2032\u209a vs \u03B8\u2082', command=plot_KC1py )
but_PKC1py.grid(column=0 , row=4)
but_TKC1px=tk.Button(VA_Label ,text='table x\u2032\u209a vs \u03B8\u2082', command=table_KC1px )
but_TKC1px.grid(column=1 , row=3)
but_TKC1py=tk.Button(VA_Label ,text='table y\u2032\u209a vs \u03B8\u2082', command=table_KC1py )
but_TKC1py.grid(column=1 , row=4)



but_Palph=tk.Button(AA_Label ,text='plot \u03B1s vs \u03B8\u2082', command=plot_AA )
but_Palph.grid(column=0 , row=0)
but_Talph=tk.Button(AA_Label ,text='table \u03B1s vs \u03B8\u2082', command=table_AA )
but_Talph.grid(column=1 , row=0)
but_PKC2=tk.Button(AA_Label ,text='plot \u03B8\u2032\u2032s vs \u03B8\u2082', command=plot_KC2 )
but_PKC2.grid(column=0 , row=1)
but_TKC2=tk.Button(AA_Label ,text='table \u03B8\u2032\u2032s vs \u03B8\u2082', command=table_KC2 )
but_TKC2.grid(column=1 , row=1)
but_PAP=tk.Button(AA_Label ,text='plot a\u209a vs \u03B8\u2082', command=plot_AP )
but_PAP.grid(column=0 , row=2)
but_PAZ=tk.Button(AA_Label ,text='plot az vs \u03B8\u2082', command=plot_AZ )
but_PAZ.grid(column=1 , row=2)
but_PKC2px=tk.Button(AA_Label ,text='plot x\u2032\u2032\u209a vs \u03B8\u2082', command=plot_KC2px )
but_PKC2px.grid(column=0 , row=3)
but_PKC2py=tk.Button(AA_Label ,text='plot y\u2032\u2032\u209a vs \u03B8\u2082', command=plot_KC2py )
but_PKC2py.grid(column=0 , row=4)
but_TKC2px=tk.Button(AA_Label ,text='table x\u2032\u2032\u209a vs \u03B8\u2082', command=table_KC2px )
but_TKC2px.grid(column=1 , row=3)
but_TKC2py=tk.Button(AA_Label ,text='table y\u2032\u2032\u209a vs \u03B8\u2082', command=table_KC2py )
but_TKC2py.grid(column=1 , row=4)

but_I=tk.Button(I_Label ,text='Calculate Instant Center of Rotations', command=Cal_I )
but_I.grid(column=0 , row=0 , columnspan=2)
but_Cent=tk.Button(I_Label ,text='Draw Centrodes', command=plot_Cent )
but_Cent.grid(column=0 , row=1, columnspan=2)
l1=tk.Label(I_Label , text='for I')
l1.grid(column=0 , row=2 , sticky=tk.W)
ent=tk.Entry(I_Label  ,width=7)
ent.grid( column=0 , row=2 ,sticky=tk.W , padx=35)
l2=tk.Label(I_Label , text=':')
l2.grid(column=0 , row=2  , sticky=tk.W , padx=98)
but_PI=tk.Button(I_Label , text='Plot I vs \u03B8\u2082' , command=plot_I)
but_PI.grid(column=0 , row=2 , sticky=tk.E)
but_TI=tk.Button(I_Label , text='Table I vs \u03B8\u2082' , command=table_I)
but_TI.grid(column=1 , row=2)
but_PathI=tk.Button(I_Label , text='Draw Path of I' , command=plot_PI)
but_PathI.grid(column=0 , row=3 , columnspan=2)



but_mech=tk.Button(others_Label ,text='Run the mechanism', command=plot_mechanism )
but_mech.grid(column=0 , row=0)
but_RZx=tk.Button(others_Label ,text='plot X of center of curvature of point Z vs \u03B8\u2082', command=plot_RZx )
but_RZx.grid(column=0 , row=1)
but_RZy=tk.Button(others_Label ,text='plot Y of center of curvature of point Z vs \u03B8\u2082', command=plot_RZy )
but_RZy.grid(column=0 , row=2)
but_RZ=tk.Button(others_Label ,text='plot radii of curvature of point Z vs \u03B8\u2082', command=plot_RZ )
but_RZ.grid(column=0 , row=3)
but_DET=tk.Button(others_Label ,text='plot Determinant of loops vs \u03B8\u2082 ', command=plot_DET ) #table T , N Of P , Z
but_DET.grid(column=0 , row=5)
but_TTN=tk.Button(others_Label ,text='table T & N of points P & Z vs \u03B8\u2082 ', command=table_TN )
but_TTN.grid(column=0 , row=6)
but_TTN=tk.Button(others_Label ,text='table radii and center of curvature of point Z vs \u03B8\u2082 ', command=table_RZ )
but_TTN.grid(column=0 , row=7)

root.mainloop()

    
    
    
    
    
    
    
    
