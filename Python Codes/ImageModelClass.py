# -*- coding: cp1252 -*-
from __future__ import division 
from math import *
import numpy as np
import math as m
import matplotlib.pyplot as plt
from numpy import cross, eye, dot

'''
Program to calculate pixel yields due to interaction point at surface
Based on Zhi Li et al 2010 Phys. Med. Biol. 55 6515, essentially, 
particularized to phoswich where interaction is at surface at height h from
4 x 4 SiPM array device.
'''
class ImageModel:
    def __init__(self,eps=0.1,Ybg=0.,Lpitch=4.2,h=10.,n1=1.6,n2=1.): # pixel array pitch and model parameters		
        self.eps=eps
        self.Ybg=Ybg
        self.Lpitch=Lpitch
        self.h=h
        self.thcrit=m.asin(n2/n1)
        self.ruler=Lpitch*np.array([-1.5,-0.5,0.5,1.5])


	#images:
    def iplus(self,x):
        return self.Lpitch*4.-x
    def iminus(self,x):
        return -self.Lpitch*4-x

    def image(self,x,y,face):
		#face = 0 (+x), 1(-y), 2(-x) or 3(+y)
        if face==0:
            return (self.iplus(x), y)
        if face==1:
            return (x, self.iminus(y))
        if face==2:
            return (self.iminus(x), y)
        if face==3:
            return (x, self.iplus(y))	 

    def imagelist0(self,x,y):
        return [(x,y)]

    def imagelist1(self,x,y):
        il=[self.image(x,y,0)]
        for f in range(1,4):
            il.append(self.image(x,y,f))
        return il		#########
	
    def imagelist2(self,x,y):
        il=[]
        f1=0
        for xy in self.imagelist1(x,y):
            xi1,yi1=xy
            for n in range(1,4):
                f2=f1+n
                if f2>3:
                    f2=f2-4
				#print("f1=",f1," f2=",f2)
                il.append(self.image(xi1,yi1,f2))
            f1=f1+1
        return il
		
    def Threfl(self,xi,yi,xl,yl):
        return m.atan(self.h/m.sqrt(pow(xi-xl,2)+pow(yi-yl,2)))

    def fRefl(self,xi,yi,xl,yl):
		#if(self.Threfl(xi,yi,xl,yl)>self.thcrit):
			##else:
                return self.eps

    def YpixModel(self,x,y):
        Ypat=np.array([])
        for xp in self.ruler:
            for yp in self.ruler:
                Ypat=np.append(Ypat,self.Yield(x,y,xp,yp)+self.Ybg)
        Ypat=Ypat.reshape(4,4)
        return Ypat
    '''
    def Path(self,x,y,xl,yl):
        r = ((x-xl)**2 + (y-yl)**2 + self.h**2)**(1/2)
        return r

    def abs_exp(self,r,absl):
        absfactor = exp(-r/absl)
        return absfactor
    '''
    def Yield(self,x,y,xl,yl): # Y at pixel xl, yl
        #r = self.Path(x,y,xl,yl)
        #teste = self.abs_exp(r,self.absl)
        #print(teste)
        Y=W(x-xl,y-yl,1.,1.,self.h)# image 0 (object)
        for xy in self.imagelist1(x,y): # first reflections
            (xi,yi)=xy	
            Y=Y+self.fRefl(xi,yi,xl,yl)*W(xi-xl,yi-yl,1.,1.,self.h)      #Colocar uma função de absorção com fator de absorção minimizado também
        for xy in self.imagelist2(x,y): # second reflections
            (xi,yi)=xy
            Y=Y+pow(self.fRefl(xi,yi,xl,yl),2)*W(xi-xl,yi-yl,1.,1.,self.h)
        return Y
    '''
    def Yield(self,x,y,xl,yl): # Y at pixel xl, yl
        r = self.Path(x,y,xl,yl)
        #teste = self.abs_exp(r,self.absl)
        #print(teste)
        Y=W(x-xl,y-yl,1.,1.,self.h)*self.abs_exp(r,self.absl)# image 0 (object)
        for xy in self.imagelist1(x,y): # first reflections
            (xi,yi)=xy	
            Y=Y+self.fRefl(xi,yi,xl,yl)*W(xi-xl,yi-yl,1.,1.,self.h)*self.abs_exp(r,self.absl)      #Colocar uma função de absorção com fator de absorção minimizado também
        for xy in self.imagelist2(x,y): # second reflections
            (xi,yi)=xy
            Y=Y+pow(self.fRefl(xi,yi,xl,yl),2)*W(xi-xl,yi-yl,1.,1.,self.h)*self.abs_exp(r,self.absl)
        return Y
   '''
    def YPattern(self,Ypm,direction,istonormalize=1): # Lines pattern: direction=0, Columns: direction=1
        Ypa=Ypm
        Ypat = np.sum(Ypa,axis=(1-direction))     #Can be altered to use pixels alone
        S=1.
        if istonormalize:
            S=np.sum(Ypat)
        return Ypat/S


#print (YpixPattern(-6,0))


def Wel(x,y,D): # Solid angle of "elementary" rectangle with center at x,y and corner at 0,0
	alfa=2.*x/D
	beta=2.*y/D
	return m.atan(alfa*beta/m.sqrt(1+alfa*alfa+beta*beta)) # negative, if x*y<0

def W(x,y,a,b,D): # Solid angle of rect. with center at x,y, side a (x) and  height b (y) 
	sx=abs(x)/2.
	sy=abs(y)/2.
	qa=a/4
	qb=b/4
	X=np.array([sx+qa,sx+qa,sx-qa,sx-qa])
	Y=np.array([sy+qb,sy-qb,sy+qb,sy-qb])
	W=Wel(X[0],Y[0],D)-Wel(X[1],Y[1],D)-Wel(X[2],Y[2],D)+Wel(X[3],Y[3],D)
	return W
#test
'''
IM=ImageModel(0.1,0.,4.2,10.,1.6,1)
print (IM.YpixModel(0.,0.))
'''

