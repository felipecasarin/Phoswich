# -*- coding: cp1252 -*-
from __future__ import division 
from math import *
import numpy as np
import matplotlib.pyplot as plt
from numpy import cross, eye, dot
import ImageModelClass as IMC
import ExpPat as EX
import ExpReCal as ERC
import Fits as Ft
import scipy as sp
from scipy.optimize import least_squares

#list of experimental data point positions
xtab=EX.xtab
ytab=EX.ytab
Nx=xtab.size
Ny=ytab.size

ExDataTab=EX.ExDataTab #data: ExDataTab=intable.reshape(Nx,Ny,2,4)
RelErrorDat=EX.RelErrorDat

# hack ExDataTab=ExDataTab/ExDataTab

#Initial parameters
eps=np.loadtxt("par_eps.txt")
Ybg=np.loadtxt("par_Ybg0.txt")
a=np.loadtxt("par_a.txt")
a=np.reshape(a,(4,4))
b=np.loadtxt("par_b.txt")
b=np.reshape(b,(2,4))
c=np.loadtxt("par_c.txt")
c=np.reshape(c,(2,4))
xoff=0.
yoff=0.

IM=IMC.ImageModel(eps,Ybg)
ExN=ERC.normdata(ERC.recaldata(ExDataTab,b,c))
errorscalefactor=1.
SigDat=ExN*RelErrorDat*errorscalefactor


YpixModel=IM.YpixModel #model: YpixPattern(x,y,ep=eps, a=np.ones((4,4)),istonormalize=1):
YPattern=IM.YPattern	#model: YPattern(Ypm,direction,a=np.ones((4,4)),b=np.ones((2,4)),istonormalize=1)

Ytable=Ft.Ytable
errorfunc=Ft.errorfunc
PatternError=Ft.PatternError

fitpat=Ft.fitpat

Yfinal=Ytable(xoff,yoff,1,a)

#hack 2 ExDataTab=Yfinal




def PlotPatt(ixtab,iytab):
	X=IM.ruler
	fig, axs = plt.subplots(1, 2)
	#axs[0].clear()
	#axs[1].clear()
	axs[0].errorbar(X,ExN[ixtab][iytab][0],SigDat[ixtab][iytab][0])
	axs[0].plot(X,Yfinal[ixtab][iytab][0])
	axs[1].errorbar(X,ExN[ixtab][iytab][1],SigDat[ixtab][iytab][1])
	axs[1].plot(X,Yfinal[ixtab][iytab][1])
	#axs[0].set_ylim(bottom=0.)
	#axs[1].set_ylim(bottom=0.)
	axs[0].set(xlabel='xl (mm)', ylabel='Y',
       title='Lines '+'x='+str(xtab[ixtab]-xoff))
	axs[1].set(xlabel='yc (mm)', title='Columns '+'y='+str(ytab[iytab]-yoff))
	#print(ixtab,iytab,"rmsdev=",sqrt(np.average(np.power(Yfinal[ixtab][iytab]-ExN[ixtab][iytab],2))))
	PE=PatternError(xtab[ixtab]-xoff,ytab[iytab]-yoff,eps,a,b,ExN[ixtab][iytab][0:2])
	axs[0].scatter(X,PE[0:4],marker='s')
	axs[1].scatter(X,PE[4:8],marker='s')
	axs[0].grid()
	axs[1].grid()
	plt.show()
	#plt.draw()
	#plt.pause(0.01)

def PlotPattSeq():
	for ixtab in range(0,Nx):
		for iytab in range(0,Ny):
			PlotPatt(ixtab,iytab)
			#plt.show()

vYpx=np.vectorize(YpixModel)
vYpt=np.vectorize(YPattern)

def PlotFunX(iytab,ilc):
	clc=np.array(['Lines','Columns'])
	fig, ax = plt.subplots()
	for ifil in range(0,4):	
		plt.scatter(xtab-xoff,ExN[:,iytab,ilc,ifil],marker="s")	
	plt.gca().set_prop_cycle(None)
	for ifil in range(0,4):
		plt.plot(xtab-xoff,Yfinal[:,iytab,ilc,ifil],label=str(ifil))
	ax.set(xlabel='x (mm)', ylabel='Y',
       title=clc[ilc]+' y='+str(ytab[iytab]))
	ax.legend()	
	ax.set_ylim(bottom=0.)
	plt.show()

def PlotFunY(ixtab,ilc):
	clc=np.array(['Lines','Columns'])
	fig, ax = plt.subplots()
	for ifil in range(0,4):	
		plt.scatter(ytab-yoff,ExN[ixtab,:,ilc,ifil],marker="s")	
	plt.gca().set_prop_cycle(None)
	for ifil in range(0,4):
		plt.plot(ytab-yoff,Yfinal[ixtab,:,ilc,ifil],label=str(ifil))
	ax.set(xlabel='y (mm)', ylabel='Y',
       title=clc[ilc]+' x='+str(xtab[ixtab]-xoff))
	ax.legend()	
	ax.set_ylim(bottom=0.)
	plt.show()
	
def PlotFunsX(ilc):
	for iytab in range(0,Ny):
		PlotFunX(iytab,ilc)
		print ("Y=",ytab[iytab]-yoff)

def PlotFunsY(ilc):
	for ixtab in range(0,Nx):
		PlotFunY(ixtab,ilc)
		print ("X=",xtab[ixtab]-xoff)

# --------------------------------------------
# --------------------------------------------

def PlotPattFit(x,y,ixtab,iytab):
	Ypm=YpixModel(x,y)
	X=IM.ruler
	fig, axs = plt.subplots(1, 2)
	axs[0].clear()
	axs[1].clear()
	axs[0].scatter(X,ExN[ixtab][iytab][0])
	axs[0].plot(X,YPattern(Ypm,0),label='{:6.2f}'.format(x))
	axs[1].scatter(X,ExN[ixtab][iytab][1])
	axs[1].plot(X,YPattern(Ypm,1),label='{:6.2f}'.format(y))
	#axs[0].set_ylim(bottom=0.)
	#axs[1].set_ylim(bottom=0.)
	axs[0].set(xlabel='xl (mm)', ylabel='Y',
       title='Lines '+'x='+str(xtab[ixtab]-xoff))
	axs[1].set(xlabel='yc (mm)', ylabel='Y',
       title='Columns '+'y='+str(ytab[iytab])-yoff)
	YLC=np.stack((YPattern(Ypm,0),YPattern(Ypm,1)))
	#print(ixtab,iytab,"rmsdev=",sqrt(np.average(np.power(YLC-ExN[ixtab][iytab],2))))
	PE=PatternError(x,y,eps,a,b,c,ExN[ixtab][iytab][0:2])
	axs[0].scatter(X,PE[0:4],marker='s')
	axs[1].scatter(X,PE[4:8],marker='s')
	axs[0].grid()
	axs[1].grid()
	plt.show()
	#plt.draw()
	#plt.pause(0.01)
	
def PlotBoth(x,y,ixtab,iytab,istoshow=1):
	Ypm=YpixModel(x,y)
	X=IM.ruler
	fig, axs = plt.subplots(1, 2)
	#axs[0].clear()
	#axs[1].clear()
	axs[0].errorbar(X,ExN[ixtab][iytab][0],SigDat[ixtab][iytab][0],xerr=0.5,fmt='o',color='c')
	axs[0].plot(X,Yfinal[ixtab][iytab][0],color='c')	
	axs[0].plot(X,YPattern(Ypm,0),linestyle='dotted',color='c',label='{:5.2f}'.format(x))
	axs[1].errorbar(X,ExN[ixtab][iytab][1],SigDat[ixtab][iytab][0],xerr=0.5,fmt='o',color='c')
	axs[1].plot(X,Yfinal[ixtab][iytab][1],color='c')
	axs[1].plot(X,YPattern(Ypm,1),linestyle='dotted',color='c',label='{:5.2f}'.format(y))
	#axs[0].set_ylim(bottom=0.)
	#axs[1].set_ylim(bottom=0.)
	axs[0].set(xlabel='xl (mm)', ylabel='Y',
       title='Lines '+'x='+"{:5.2f}".format(xtab[ixtab]-xoff))
	axs[1].set(xlabel='yc (mm)', title='Columns '+'y='+"{:5.2f}".format(ytab[iytab]-yoff))
	desv=Yfinal[ixtab][iytab]-ExN[ixtab][iytab]
	print(ixtab,iytab,"{:6.3f}".format(xtab[ixtab]-xoff),"{:6.3f}".format(ytab[iytab]-yoff),"{:6.3f}".format(x),"{:6.3f}".format(y),sqrt(np.average(np.power(desv,2))),desv[0],desv[1])
	PE=PatternError(xtab[ixtab]-xoff,ytab[iytab]-yoff,ExN[ixtab][iytab][0:2])
	axs[0].scatter(X,PE[0:4],marker='s',color='r')
	axs[1].scatter(X,PE[4:8],marker='s',color='r')
	PE=PatternError(x,y,ExN[ixtab][iytab][0:2])
	axs[0].scatter(X,PE[0:4],marker='x',color='m')
	axs[1].scatter(X,PE[4:8],marker='x',color='m')
	axs[0].grid()
	axs[1].grid()
	axs[0].legend()	
	axs[1].legend()	
	plt.savefig("Patt"+str(ixtab)+str(iytab)+".png", dpi=150)
	if istoshow:
		plt.show()
	plt.close()
	
def PlotAllBoth(istoshow=1):
	print("eps=",eps)
	print("Ybg",Ybg)
	print("b=\n",b)
	print("c=\n",c)
	print("ixtab, iytab, xtab, ytab, x, y, rms")
	for ixtab in range(0,Nx):
		for iytab in range(0,Ny):
			xy=fitpat(ixtab,iytab)
			PlotBoth(xy[0],xy[1],ixtab,iytab,istoshow)

def PlotYSumLminusC(b=np.ones((2,4)),c=np.zeros((2,4))):
	scal=100
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')		
	Y0=Ytable(0.,0.,0)
	scal=30
	for ixtab in range(0,Nx):
		for iytab in range(0,Ny):
			if not ( (ixtab==6) and (iytab==5) ):
				L=np.sum(Y0[ixtab][iytab][0][0:4]*b[0]+c[0])
				C=np.sum(Y0[ixtab][iytab][1][0:4]*b[1]+c[1])
				ax.scatter(ixtab,iytab,2*(L-C)/(L+C),c='b')				
				#ax.scatter(ixtab,iytab,np.sum(ExDataTab[ixtab][iytab][1][0:4]),c='r')
	plt.show()
	
def PlotYSumLandC(b=np.ones((2,4)),c=np.zeros((2,4))):
	scal=100
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')		
	Y0=Ytable(0.,0.,0)	
	for ixtab in range(0,Nx):
		for iytab in range(0,Ny):
			if not ( (ixtab==6) and (iytab==5) ):
				L=np.sum(Y0[ixtab][iytab][0][0:4]*b[0]+c[0])
				C=np.sum(Y0[ixtab][iytab][1][0:4]*b[1]+c[1])
				ax.scatter(ixtab,iytab,L,c='b')				
				ax.scatter(ixtab,iytab,C,c='r')
	plt.show()

#PlotYSumLminusC(btest,ctest)


def PlotSumLminusC(ExData):
	scal=100
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')	
	for ixtab in range(0,Nx):
		for iytab in range(0,Ny):
			if not ( (ixtab==6) and (iytab==5) ):
				L=np.sum(ExData[ixtab][iytab][0][0:4])
				C=np.sum(ExData[ixtab][iytab][1][0:4])
				#ax.scatter(ixtab,iytab,2*(L-C)/(L+C),c='b')				
				ax.scatter(ixtab,iytab,(L-C),c='b')				
				#ax.scatter(ixtab,iytab,np.sum(ExDataTab[ixtab][iytab][1][0:4]),c='r')
	plt.show()
	
def PlotSumLandC(ExData):
	scal=100
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')	
	for ixtab in range(0,Nx):
		for iytab in range(0,Ny):
			if not ( (ixtab==6) and (iytab==5) ):
				L=np.sum(ExData[ixtab][iytab][0][0:4])
				C=np.sum(ExData[ixtab][iytab][1][0:4])
				ax.scatter(ixtab,iytab,L,c='b')				
				ax.scatter(ixtab,iytab,C,c='r')
	plt.show()
	
	
def PlotSumLoverC(ExData):
	scal=100
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')	
	ax.set(xlabel="x",ylabel="y",zlabel="sL/sC")
	for ixtab in range(0,Nx):
		for iytab in range(0,Ny):
			if not ( (ixtab==6) and (iytab==5) ):
				L=np.sum(ExData[ixtab][iytab][0][0:4])
				C=np.sum(ExData[ixtab][iytab][1][0:4])
				ax.scatter(xtab[ixtab],ytab[iytab],L/C,c='b')				
	plt.show()
    
    
    
def teste():
    b=np.array([])
    for ix in range (0,7):
        b=np.append(b, xtab)
    print(b)


#PlotSumLoverC(ExDataTab)

#PlotAllBoth(0)

#b=np.ones((2,4))
#c=np.zeros((2,4))
#b[0][3]=10.
#

print("eps, Ybg:",eps,Ybg)
print("a:\n",a)
print("b:\n",b)
print("c:\n",c)
#PlotSumLminusC(ERC.recaldata(ExDataTab,b,c))
#print(ERC.recaldata(ExDataTab,b,c))
#PlotSumLminusC(ExDataTab)
#PlotSumLoverC(ERC.recaldata(ExDataTab,b,c))
#PlotYSumLminusC(b,c)

