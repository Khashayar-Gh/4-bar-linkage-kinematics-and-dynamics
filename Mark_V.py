#%% importing requiring moduls
import time
t=time.time()
import numpy as np
import scipy.optimize
from tkinter import ttk

#import math
#import matplotlib.pyplot as plt
import tkinter as tk

from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.
#from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
#import matplotlib
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
'''
plt.figure(1)
plt.plot(tet2 , post[:,0],label = 'theta3')
plt.plot(tet2 , post[:,1] , label = 'theta4')
plt.plot(tet2 , post[:,2] , label = 'theta5')
plt.plot(tet2 , post[:,3] , label = 'theta6')
plt.show()
plt.legend()
'''
#a = np.transpose

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
KC1p=np.zeros([3 , n]) # 1st order KC for P
KC2p=KC1p
rho=np.zeros([1 , n]) # R of curv. for z
C=np.zeros([3 , n]) # center of curv. for Z
PP[: , 0]= (R[0]-R11)*np.array([1 , 0 , 0])+R[5]*np.array([cosd(post[0 , 3]) , sind(post[0 , 3]) , 0])+R[4]/2*np.array([cosd(post[0 , 2]) ,sind(post[0 , 2])  , 0])
PZ[: , 0]= (R[0]-R11)*np.array([1 , 0 , 0])+R[5]*np.array([cosd(post[0 , 3]), sind(post[0 , 3]) , 0])+R[4]/4*np.array([cosd(post[0 , 2]) ,sind(post[0 , 2])  , 0 ])
for i in range(1 , n):

    PP[: , i]= (R[0]-R11)*np.array([1 , 0 , 0]) + R[5]*np.array([cosd(post[i , 3]) , sind(post[i , 3]) , 0]) + (R[4]/2)*np.array([cosd(post[i , 2]) ,sind(post[i , 2])  , 0])
    PZ[: , i]= (R[0]-R11)*np.array([1 , 0 , 0]) + R[5]*np.array([cosd(post[i , 3]) , sind(post[i , 3]) , 0]) + (R[4]/4)*np.array([cosd(post[i , 2]) ,sind(post[i , 2])  , 0])
    dp=dp + np.linalg.norm((PP[: , i]- PP[: , i-1]))
    dz=dz + np.linalg.norm((PZ[: , i] - PZ[: , i-1]))

#%% GUI
def plot_tet():
    listBox.place(relx=.5 , rely=.5,relheight=0 , relwidth=0)

    #root.geometry("{0}x{1}+0+0".format(root.winfo_screenwidth(), root.winfo_screenheight()))
    try: 
        #canvas.get_tk_widget().pack_forget()
        ax1.clear()
    except AttributeError: 
        pass     
    #fig = Figure(figsize=(5, 4), dpi=100)
    #t = np.arange(0, 3, .01)
    #fig.add_subplot(111).plot(t, 2 * np.sin(2 * np.pi * t))
    #plt.figure(1)
    #ax1=fig.add_subplot(111)
    ax1.plot(tet2 , post[:,0],label = 'theta3')
    ax1.plot(tet2 , post[:,1] , label = 'theta4')
    ax1.plot(tet2 , post[:,2] , label = 'theta5')
    ax1.plot(tet2 , post[:,3] , label = 'theta6')
    # A tk.DrawingArea.
    ax1.legend()
    canvas.draw()
    #canvas.geometry('100x100')
    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    ''' 
    toolbar = NavigationToolbar2Tk(canvas, results)
    '''
    toolbar.update()
    
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
def plot_AV():
    listBox.place(relx=.5 , rely=.5,relheight=0 , relwidth=0)

    #root.geometry("{0}x{1}+0+0".format(root.winfo_screenwidth(), root.winfo_screenheight()))
    try: 
        #canvas.get_tk_widget().pack_forget()
        ax1.clear()

    except AttributeError: 
        pass     
    #fig = Figure(figsize=(5, 4), dpi=100)
    #t = np.arange(0, 3, .01)
    #fig.add_subplot(111).plot(t, 2 * np.sin(2 * np.pi * t))
    #plt.figure(1)
    #ax1=fig.add_subplot(111)
    ax1.plot(tet2 , AV[: , 2],label = 'w3')
    ax1.plot(tet2 , AV[: , 3] , label = 'w4')
    ax1.plot(tet2 , AV[: , 4], label = 'w5')
    ax1.plot(tet2 , AV[: , 5] , label = 'w6')
    #canvas = FigureCanvasTkAgg(fig, master=results)  # A tk.DrawingArea.
    ax1.legend()
    canvas.draw()
    #canvas.geometry('100x100')
    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    '''
    toolbar = NavigationToolbar2Tk(canvas, results)
    '''
    toolbar.update()
    
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
def plot_KC1():
    listBox.place(relx=.5 , rely=.5,relheight=0 , relwidth=0)
    #root.geometry("{0}x{1}+0+0".format(root.winfo_screenwidth(), root.winfo_screenheight()))
    try: 
        #canvas.get_tk_widget().pack_forget()
        ax1.clear()
        
    except:
        pass   
    #fig = Figure(figsize=(5, 4), dpi=100)
    #t = np.arange(0, 3, .01)
    #fig.add_subplot(111).plot(t, 2 * np.sin(2 * np.pi * t))
    #plt.figure(1)
    #ax1=fig.add_subplot(111)
    ax1.plot(tet2 , KC1[: , 0],label = 'bar3')
    ax1.plot(tet2 , KC1[: , 1] , label = 'bar4')
    ax1.plot(tet2 , KC1[: , 2], label = 'bar5')
    ax1.plot(tet2 , KC1[: , 3] , label = 'bar6')
    #canvas = FigureCanvasTkAgg(fig, master=results)  # A tk.DrawingArea.
    ax1.legend()
    canvas.draw()
    #canvas.geometry('100x100')
    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    
    #toolbar = NavigationToolbar2Tk(canvas, results)
    toolbar.update()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
def plot_KC2():
    listBox.place(relx=.5 , rely=.5,relheight=0 , relwidth=0)
    #root.geometry("{0}x{1}+0+0".format(root.winfo_screenwidth(), root.winfo_screenheight()))
    try: 
        #canvas.get_tk_widget().pack_forget()
        ax1.clear()
        
    except:
        pass   
    #fig = Figure(figsize=(5, 4), dpi=100)
    #t = np.arange(0, 3, .01)
    #fig.add_subplot(111).plot(t, 2 * np.sin(2 * np.pi * t))
    #plt.figure(1)
    #ax1=fig.add_subplot(111)
    ax1.plot(tet2 , KC2[: , 0],label = 'bar3')
    ax1.plot(tet2 , KC2[: , 1] , label = 'bar4')
    ax1.plot(tet2 , KC2[: , 2], label = 'bar5')
    ax1.plot(tet2 , KC2[: , 3] , label = 'bar6')
    #canvas = FigureCanvasTkAgg(fig, master=results)  # A tk.DrawingArea.
    ax1.legend()
    canvas.draw()
    #canvas.geometry('100x100')
    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    
    #toolbar = NavigationToolbar2Tk(canvas, results)
    toolbar.update()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
def plot_AA():
    listBox.place(relx=.5 , rely=.5,relheight=0 , relwidth=0)
    #root.geometry("{0}x{1}+0+0".format(root.winfo_screenwidth(), root.winfo_screenheight()))
    try: 
        #canvas.get_tk_widget().pack_forget()
        ax1.clear()
        
    except:
        pass   
    #fig = Figure(figsize=(5, 4), dpi=100)
    #t = np.arange(0, 3, .01)
    #fig.add_subplot(111).plot(t, 2 * np.sin(2 * np.pi * t))
    #plt.figure(1)
    #ax1=fig.add_subplot(111)
    ax1.plot(tet2 , AA[: , 2],label = '\u03B1\u2083')
    ax1.plot(tet2 , AA[: , 3] , label = '\u03B1\u2084')
    ax1.plot(tet2 , AA[: , 4], label = '\u03B1\u2085')
    ax1.plot(tet2 , AA[: , 5] , label = '\u03B1\u2086')
    #canvas = FigureCanvasTkAgg(fig, master=results)  # A tk.DrawingArea.
    ax1.legend()
    canvas.draw()
    #canvas.geometry('100x100')
    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    
    #toolbar = NavigationToolbar2Tk(canvas, results)
    toolbar.update()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
def plot_PathP():
    listBox.place(relx=.5 , rely=.5,relheight=0 , relwidth=0)
    #root.geometry("{0}x{1}+0+0".format(root.winfo_screenwidth(), root.winfo_screenheight()))
    try: 
        #canvas.get_tk_widget().pack_forget()
        ax1.clear()
        
    except:
        pass   
    #fig = Figure(figsize=(5, 4), dpi=100)
    #t = np.arange(0, 3, .01)
    #fig.add_subplot(111).plot(t, 2 * np.sin(2 * np.pi * t))
    #plt.figure(1)
    #ax1=fig.add_subplot(111)
    ax1.plot(PP[0 , :] , PP[1 , :],label = 'Path of Point P')
    #canvas = FigureCanvasTkAgg(fig, master=results)  # A tk.DrawingArea.
    ax1.legend()
    canvas.draw()
    #canvas.geometry('100x100')
    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    
    #toolbar = NavigationToolbar2Tk(canvas, results)
    toolbar.update()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
def plot_PathZ():
    listBox.place(relx=.5 , rely=.5,relheight=0 , relwidth=0)
    #root.geometry("{0}x{1}+0+0".format(root.winfo_screenwidth(), root.winfo_screenheight()))
    try: 
        #canvas.get_tk_widget().pack_forget()
        ax1.clear()
        
    except:
        pass   
    #fig = Figure(figsize=(5, 4), dpi=100)
    #t = np.arange(0, 3, .01)
    #fig.add_subplot(111).plot(t, 2 * np.sin(2 * np.pi * t))
    #plt.figure(1)
    #ax1=fig.add_subplot(111)
    ax1.plot(PZ[0 , :] , PZ[1 , :],label = 'Path of Point Z')
    #canvas = FigureCanvasTkAgg(fig, master=results)  # A tk.DrawingArea.
    ax1.legend()
    canvas.draw()
    #canvas.geometry('100x100')
    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    
    #toolbar = NavigationToolbar2Tk(canvas, results)
    toolbar.update()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    
    
def table_tet():
    try: 
        canvas.get_tk_widget().pack_forget()
        ax1.clear()
        
    except:
        pass   
    
    listBox.delete(*listBox.get_children())
    listBox.place(relx=.255,rely=.035,relheight=.9, relwidth=.735)

    for i in range(0 , n , 50):
        listBox.insert("", "end", values=(tet2[i],post[i , 0], post[i , 1], post[i , 2] ,post[i , 3] ))

    #toolbar.update()
def table_AV():
    try: 
        canvas.get_tk_widget().pack_forget()
        ax1.clear()
        
    except:
        pass   
    listBox.delete(*listBox.get_children())
    listBox.place(relx=.255,rely=.035,relheight=.9, relwidth=.735)

    for i in range(0 , n , 50):
        listBox.insert("", "end", values=(tet2[i],AV[i , 2], AV[i , 3], AV[i , 4] ,AV[i , 5] ))

    #toolbar.update()
def table_AA():
    
    try: 
        canvas.get_tk_widget().pack_forget()
        ax1.clear()
        
    except:
        pass   
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
    listBox.delete(*listBox.get_children())
    listBox.place(relx=.255,rely=.035,relheight=.9, relwidth=.735)

    for i in range(0 , n , 50):
        listBox.insert("", "end", values=(tet2[i] , KC2[i , 0], KC2[i , 1], KC2[i , 2] ,KC2[i , 3] ))


elap=time.time()-t

root=tk.Tk()
root.geometry('1280x720')
root.configure(background='cyan')

results=tk.LabelFrame(root, text="Results", padx=5, pady=5 , bg='green')
results.place(relx=.25 , rely=0 , relheight=.95 , relwidth=.75)
fig = Figure(figsize=(5, 4), dpi=100)
ax1=fig.add_subplot(111)
canvas = FigureCanvasTkAgg(fig, master=results)
toolbar = NavigationToolbar2Tk(canvas, root)
toolbar.config(background='cyan')
cols = ( '\u03B8\u2082', 'bar3','bar4','bar5','bar6')
listBox = ttk.Treeview(root, columns=cols, show='headings')
for col in cols:
    listBox.heading(col, text=col)   

#listBox.place(rely=.035,relheight=.001 , relwidth=.75)
# set column headings
#buttest=tk.Button(results , text=';fdghelkrjgh' , bg='green')
#buttest.pack()

PA_Label = tk.LabelFrame(root, text="Posture Analisys", padx=5, pady=5 , bg='pink')
PA_Label.place(relx=0.0005 , rely=0.0005)
VA_Label=tk.LabelFrame(root, text="Velocity Analisys", padx=5, pady=5 , bg='pink')
VA_Label.place(relx=0.0005 , rely=0.3)
AA_Label=tk.LabelFrame(root, text="Acceleration Analisys", padx=5, pady=5 , bg='pink')
AA_Label.place(relx=0.0005 , rely=0.5)

#PA_Label.pack(padx=10, pady=10)

but_Ptet=tk.Button(PA_Label ,text='plot \u03B8s vs \u03B8\u2082', command=plot_tet  , padx=10, pady=5 )
#but2=tk.Button(root ,text='sk,hfgqwkdgfhwquli', command=plot_tet)
but_Ptet.grid(column=0 , row=0)
but_Ttet=tk.Button(PA_Label ,text='table \u03B8s vs \u03B8\u2082', command=table_tet , padx=5, pady=5)
but_Ttet.grid(column=1 , row=0)
but_PathP=tk.Button(PA_Label ,text='draw path of Point P', command=plot_PathP  , padx=5, pady=5 )
but_PathP.grid(column=0 , row=1 , columnspan=2)
but_PathZ=tk.Button(PA_Label ,text='draw path of Point Z', command=plot_PathZ  , padx=5, pady=5 )
but_PathZ.grid(column=0 , row=2 , columnspan=2)



but_Pw=tk.Button(VA_Label ,text='plot \u03C9s vs \u03B8\u2082', command=plot_AV )
but_Pw.grid(column=0 , row=0)
but_Tw=tk.Button(VA_Label ,text='table \u03C9s vs \u03B8\u2082', command=table_AV )
but_Tw.grid(column=1 , row=0)
but_PKC1=tk.Button(VA_Label ,text='plot KC1 vs \u03B8\u2082', command=plot_KC1 )
but_PKC1.grid(column=0 , row=1)
but_TKC1=tk.Button(VA_Label ,text='table KC1 vs \u03B8\u2082', command=table_KC1 )
but_TKC1.grid(column=1 , row=1)



but_Palph=tk.Button(AA_Label ,text='plot \u03B1s vs \u03B8\u2082', command=plot_AA )
but_Palph.grid(column=0 , row=0)
but_Talph=tk.Button(AA_Label ,text='table \u03B1s vs \u03B8\u2082', command=table_AA )
but_Talph.grid(column=1 , row=0)
but_PKC2=tk.Button(AA_Label ,text='plot KC2 vs \u03B8\u2082', command=plot_KC2 )
but_PKC2.grid(column=0 , row=1)
but_TKC2=tk.Button(AA_Label ,text='table KC2 vs \u03B8\u2082', command=table_KC2 )
but_TKC2.grid(column=1 , row=1)
#but2.pack(side=tk.RIGHT)

root.mainloop()

    
    
    
    
    
    
    
    
