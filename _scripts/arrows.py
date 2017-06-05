# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 10:06:57 2016

@author: JanKM
"""

import matplotlib.pyplot as plt
import numpy as np
from numpy.random import seed, randint

# Build function that connects two points with a curved line, 
# and an arrow on the middle of it

#seed(1679)

narrow = 3
rad_one = 50
numpoints = 3

random_points = list(randint(1,20,[numpoints,4]))
rpoints = [[(a,b),(c,d)] for a,b,c,d in random_points]

def curvline(start,end,rad,t=100,arrows=1,push=0.8):
    #Compute midpoint
    rad = rad/100.    
    x1, y1 = start
    x2, y2 = end
    y12 = (y1 + y2) / 2
    dy = (y2 - y1)
    cy = y12 + (rad) * dy
    #Prepare line
    tau = np.linspace(0,1,t)
    xsupport = np.linspace(x1,x2,t)
    ysupport = [(1-i)**2 * y1 + 2*(1-i)*i*cy + (i**2)*y2 for i in tau]
    #Create arrow data    
    arset = list(np.linspace(0,1,arrows+2))
    c = zip([xsupport[int(t*a*push)] for a in arset[1:-1]],
                      [ysupport[int(t*a*push)] for a in arset[1:-1]])
    dt = zip([xsupport[int(t*a*push)+1]-xsupport[int(t*a*push)] for a in arset[1:-1]],
                      [ysupport[int(t*a*push)+1]-ysupport[int(t*a*push)] for a in arset[1:-1]])
    arrowpath = zip(c,dt)
    return xsupport, ysupport, arrowpath

def plotcurv(start,end,rad,t=100,arrows=1,arwidth=.25):
    x, y, c = curvline(start,end,rad,t,arrows)
    plt.plot(x,y,'k-')
    for d,dt in c:
        plt.arrow(d[0],d[1],dt[0],dt[1], shape='full', lw=0, 
                  length_includes_head=False, head_width=arwidth)
    return c

#Create figure
figure = plt.figure()
ax = plt.subplot(111)
for n1,n2 in rpoints:    
    #First line
    plotcurv(n1,n2,rad_one,200,narrow,0.5)
    #Second line
    plotcurv(n2,n1,rad_one,200,narrow,0.5)
ax.set_xlim(0,20)
ax.set_ylim(0,20)
plt.show