# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 23:50:07 2019

@author: Khashigh
"""

#%%plot 
import matplotlib.pyplot as plt
fig = plt.figure() 
ax = plt.axes(xlim=(-60, 240), ylim=(-70, 200)) 
line, = ax.plot([], [], lw=2.5 , color = '#ffcc00' ,label = 'Position of point P' , linestyle = ':') 
line2, = ax.plot([], [], lw=2.5 ,color = '#ff0055' , label = 'Position of point Z' , linestyle = ':')
line3, = ax.plot([], [], lw=3 ,color = 'black')
line4, = ax.plot([], [], lw=3 ,color = 'black')
line5, = ax.plot([], [], lw=3 ,color = 'black')
line6, = ax.plot([], [], lw=3 ,color = 'black')
line7, = ax.plot([], [], lw=3 ,color = 'black')
# initialization function 
def init(): 
	# creating an empty plot/frame 
	line.set_data([], []) 
	return line, 

# lists to store x and y axis points 
xdata, ydata = [], [] 
x1data, y1data = [], []
# animation function 
o6=(R[0]-R11)*np.array([1, 0])
o4=R[0]*np.array([1, 0])
def animate(i): 
	x = Ix[i , 0,2] 
	y = Iy[i, 0 , 2]
	x1 = Ix[i , 0,4]  
	y1= Iy[i, 0 , 4]
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
	return line, line2,line3, line4, line5,line6, line7
	
# hiding the axis details 
plt.axis('off') 

# call the animator	 
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=3600, interval=1, blit=True , repeat=False)
plt.legend()
plt.show()