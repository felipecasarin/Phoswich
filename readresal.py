# -*- coding: cp1252 -*-
from __future__ import division 
from math import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
tabname="resal.txt"
# reads table with experimental centroids
intable=np.loadtxt(tabname)
#print ("ixtab, iytab, xtab, ytab, x, y, rms, dL, dC")
#print(intable)
t=np.transpose(intable)
#print(t)
rs=np.split(t,t.shape[0])
ixtab=rs[0]
iytab=rs[1]
xtab=rs[2]
ytab=rs[3]
x=rs[4]
y=rs[5]
rms=rs[6]
dL=rs[7:11]
dC=rs[11:15]
scal=100

print(xtab)
print(ytab)

#fig = plt.figure()



#ax.scatter(xtab,ytab,rms,s=scal)




# Plot the surface.
#surf = ax.plot_wireframe(x, y, ((x-xtab)**2 + (y-ytab)**2)**(0.5), cmap=cm.coolwarm, linewidth=0, antialiased=True)
#ax.set_title('Surface plot')
'''
#3d plot of points and error
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
graph = ax.scatter(xtab,ytab,((x-xtab)**2 + (y-ytab)**2)**(0.5),s=100)
ax.hold(True)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('Módulo do erro')
plt.show()

'''



'''
#2d figure with reference circle
fig, ax = plt.subplots()
ax.quiver(0, 0,x-xtab, y-ytab, scale=1, angles='xy', scale_units = 'xy')
circle1 = plt.Circle((0, 0), 3, color='r', fill=False)
ax.add_patch(circle1)
plt.show()
'''


#2d figure of every hole and its measure
fig, ax = plt.subplots()
ept = plt.scatter(x,y, label="Pontos experimentais")
tpt = plt.scatter(xtab, ytab, label="Localização dos buracos")
vec= plt.quiver(xtab, ytab, x-xtab, y-ytab, label="Erro", scale=1, angles='xy', scale_units = 'xy')
plt.legend(bbox_to_anchor=(1, 1), loc=1, borderaxespad=0.)


plt.grid()

major_ticks = np.arange(-8.4, 8.4, 4.2)
minor_ticks = np.arange(-8.4, 8.4, 0.01)
ax.set_xticks(major_ticks)
ax.set_xticks(minor_ticks, minor=True)
ax.set_yticks(major_ticks)
ax.set_yticks(minor_ticks, minor=True)

'''
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Posição dos furos')
plt.show()
'''
'''
fig, ax = plt.subplots()
tpt = plt.scatter(xtab, ytab, label="")
plt.legend(bbox_to_anchor=(1, 1), loc=1, borderaxespad=0.)

plt.grid()

major_ticks = np.arange(-8.4, 8.4, 4.2)
minor_ticks = np.arange(-8.4, 8.4, 0.01)
ax.set_xticks(major_ticks)
ax.set_xticks(minor_ticks, minor=True)
ax.set_yticks(major_ticks)
ax.set_yticks(minor_ticks, minor=True)
plt.show()


'''
#ax.scatter(xtab,ytab,dL[0],s=scal)
#ax.scatter(xtab,ytab,dL[1],s=scal)
#ax.scatter(xtab,ytab,dL[2],s=scal)
#ax.scatter(xtab,ytab,dC[3],s=scal)
#'''
