from __future__ import division 
import ROOT
from ROOT import TNtuple
import sys
import numpy as np
from math import *
import matplotlib.pyplot as plt
from numpy import cross, eye, dot
import ImageModelClass as IMC
import ExpPat as EX
import ExpReCal as ERC
import scipy as sp
import Fits as ft
from scipy.optimize import least_squares
from ROOT import TFile, TTree, gROOT, addressof
from array import array
from scipy import linalg



#eps0=0.5
#Ybg0=1.E-4


tab_eps0 = "./par_eps.txt"
eps0 = np.loadtxt(tab_eps0)

tab_Ybg0 = "./par_eps.txt"
Ybg0 = np.loadtxt(tab_Ybg0)

#Import of calibration parameters b calculated by fitb() from Fits.py
tab_b="./par_b.txt"

b=np.loadtxt(tab_b)

tab_c="./par_c.txt"

c=np.loadtxt(tab_c)

#c = np.array([0,0,0,0,0,0,0,0])

IM=IMC.ImageModel(eps0,Ybg0)

Lpitch = 4.2

ruler=Lpitch*np.array([-1.5,-0.5,0.5,1.5])


####################################################################

#Opens the root file with the data to be processed
#data = ROOT.TFile.Open('/home/casarin/Desktop/rootdata/nicedata.root', 'read')
#data = ROOT.TFile.Open('/home/casarin/Desktop/rootdata/testepequeno.root', 'read') 
#data = ROOT.TFile.Open('/home/casarin/Desktop/rootdata/2000_entries_per_point.root', 'read')
#data = ROOT.TFile.Open('/home/casarin/Desktop/rootdata/2000_entries_per_point_new.root', 'read')
data = ROOT.TFile.Open('/home/casarin/Desktop/rootdata/300_entries.root', 'read')

#Gets the ntuple from the given file
ntuple = data.Get('ntuple')


#Creates the file in which the processed data will be saved afterwards
#outfile =  ROOT.TFile("/home/casarin/Desktop/rootdata/fit_3000_entries_per_point.root",'recreate')
#outfile =  ROOT.TFile("/home/casarin/Desktop/rootdata/fit_2000_entries_per_point.root",'recreate')
#outfile =  ROOT.TFile("/home/casarin/Desktop/rootdata/novotesteconsb.root",'recreate')
outfile =  ROOT.TFile("/home/casarin/Desktop/rootdata/Fitting_test/fit_300_newoffset2.root",'recreate')


print('Programa rodando')

#Calculates the Pattern Error Function using xguess, yguess and the Yield Data
def fpatxy(X,Y_dat): 
	x=X[0]
	y=X[1]
	f=ft.PatternErrorOffset(x,y,Y_dat)   #<====PatternErrorOffset??
	return f



def Ytable(xhole,yhole,xoff=-0.4565,yoff=-0.0122,istonorm=1,a=np.ones((4,4))): # theoretical yields calculated using the hole position (xhole,yhole)
	YT=np.zeros((2,4))  
	x= xhole-xoff
	y= yhole-yoff
	print(x,y)        
	Ypm=IM.YpixModel(x,y)*a    
	YT[0][0:4]=IM.YPattern(Ypm,0,istonorm) #lines
	YT[1][0:4]=IM.YPattern(Ypm,1,istonorm) #columns  
	#print('break') 
	print(YT)       
	return YT

def fitpatxy(): # pattern fit of x,y
	# loads model and calib. fit parameters
    IM.eps=np.loadtxt("par_eps.txt")
    IM.Ybg=np.loadtxt("par_Ybg0.txt")
    xguess=np.sum(ruler*newExN[0][0:4])*3.
    #print('old xguess=', xguess)
    
    #This part changes the value of xguess and yguess in case they lie outside the detector, as this would lead to miscalculations of the x and y coordinate
    if xguess < -10.:
        xguess = -8.
    if xguess > 10.:
        xguess = 8.
        
    #print('new xguess=', xguess)
    
    yguess=np.sum(ruler*newExN[1][0:4])*3.
    #print('old yguess=', yguess)
    
    if yguess < -10.:
        yguess=-8.
    if yguess > 10.:
        yguess = 8.
    
    #print('new yguess=', yguess)
    
	#print("x,y guess:",xguess,yguess)
    X=np.array([xguess,yguess]) # initial x,y guess
    Lmax=25. #Determines the maximum value for x and y
    bdi=np.array([-Lmax,-Lmax]) # lower limits
    bds=np.array([Lmax,Lmax]) # upper limits	
    
    global result
    #Calculates the x and y coordinates of the interaction point of the alpha particle on the phoswich detector
    result = sp.optimize.least_squares(fpatxy,X,bounds=(bdi,bds),args = ([newExN[0:2]]), ftol=1e-08, xtol=1e-08, gtol=1e-08, loss='cauchy', tr_solver='lsmr')
    #result = sp.optimize.leastsq(fpatxy,X,bounds=(bdi,bds),args = ([newExN[0:2]]), ftol=1e-08, xtol=1e-08, gtol=1e-08, loss='cauchy', tr_solver='lsmr')
    print(result)
    print('--------------------')
    '''
    chi2dof = np.sum(result.fun**2)/(result.fun.size - result.x.size)
    cov = chi2dof
    perr = np.sqrt(cov)
    print("sigma")
    print(perr)

    '''
    J1 = result.jac
    print("Jacobiano.transpose")
    print(J1)
    #print('*------------------*')
    #print(J[1])
    J1t = J1.transpose()
    hess = np.dot(J1t,J1)
    cov = linalg.inv(hess)
    var = np.sqrt(np.diagonal(cov))
    scalecov = cov*(result.cost)
    print("var")
    print(var)
    
    global unc
    
    unc = np.sqrt(np.diagonal(scalecov))
    
    #print("uncertainty")
    #print(unc)

    return result.x

def calib(l1,l2,l3,l4,c1,c2,c3,c4):
    
   
    #Calibrated lines
    calL1 = (l1 + c[0])*b[0]
    calL2 = (l2 + c[1])*b[1]
    calL3 = (l3 + c[2])*b[2]
    calL4 = (l4 + c[3])*b[3]
    sumL = calL1 + calL2 + calL3 + calL4
    
    
    #Calibrated columns
    
    calC1 = (c1 + c[4])*b[4]
    calC2 = (c2 + c[5])*b[5]
    calC3 = (c3 + c[6])*b[6]
    calC4 = (c4 + c[7])*b[7]
    sumC = calC1 + calC2 + calC3 + calC4
    
    
    #Normalized and calibrated lines
    cL1 = calL1/(sumL)
    cL2 = calL2/(sumL)
    cL3 = calL3/(sumL)
    cL4 = calL4/(sumL)
    #print('sumcL =',cL1+cL2+cL3+cL4)
    
    #Normalized and calibrated columns
    cC1 = calC1/(sumC)
    cC2 = calC2/(sumC)
    cC3 = calC3/(sumC)
    cC4 = calC4/(sumC)
   # print('sumcC =',cC1+cC2+cC3+cC4)
    

    global calb
    calb = np.array([cL1,cL2,cL3,cL4,cC1,cC2,cC3,cC4])
   # print(l1,l2,l3,l4,c1,c2,c3,c4)
    #print(calb)

    return calb


def pos0():
    #Creates the Tree in which the new data will be saved
    tree  = ROOT.TTree( 'tree','tree')
    
    
    nu = ntuple.GetEntries()
    
    #Defines the arrays, filled with float type variables, with an initial value of zero
    aL1 = array('f',[0])
    aL2 = array('f',[0])
    aL3 = array('f',[0])
    aL4 = array('f',[0])
    aC1 = array('f',[0])
    aC2 = array('f',[0])
    aC3 = array('f',[0])
    aC4 = array('f',[0])
    aL1c = array('f',[0])
    aL2c = array('f',[0])
    aL3c = array('f',[0])
    aL4c = array('f',[0])
    aC1c = array('f',[0])
    aC2c = array('f',[0])
    aC3c = array('f',[0])
    aC4c = array('f',[0])
    aEt = array('f',[0])
    aEvent = array('f',[0])
    aXh = array('f',[0])
    aYh = array('f',[0])
    aXp = array('f',[0])
    aYp = array('f',[0])
    aCost = array('f',[0])
    anfev = array('f',[0])
    anjev = array('f',[0])
    aXt = array('f',[0])
    aYt = array('f',[0])
    aYTL1 = array('f',[0])
    aYTL2 = array('f',[0])
    aYTL3 = array('f',[0])
    aYTL4 = array('f',[0])
    aYTC1 = array('f',[0])
    aYTC2 = array('f',[0])
    aYTC3 = array('f',[0])
    aYTC4 = array('f',[0])
    aYTL1_fit = array('f',[0])
    aYTL2_fit = array('f',[0])
    aYTL3_fit = array('f',[0])
    aYTL4_fit = array('f',[0])
    aYTC1_fit = array('f',[0])
    aYTC2_fit = array('f',[0])
    aYTC3_fit = array('f',[0])
    aYTC4_fit = array('f',[0])
    auncertx = array('f',[0])
    auncerty = array('f',[0])
    
    
    
    #Creates the Branches, which lay inside the Tree and divides the data on subcategories. The /F is to indicate we're storing float variables
    bL1 = tree.Branch('L1',aL1,'L1/F')
    bL2 = tree.Branch('L2',aL2,'L2/F')
    bL3 = tree.Branch('L3',aL3,'L3/F')
    bL4 = tree.Branch('L4',aL4,'L4/F')
    bC1 = tree.Branch('C1',aC1,'C1/F')
    bC2 = tree.Branch('C2',aC2,'C2/F')
    bC3 = tree.Branch('C3',aC3,'C3/F')
    bC4 = tree.Branch('C4',aC4,'C4/F')
    bL1c = tree.Branch('L1c',aL1c,'L1c/F')
    bL2c = tree.Branch('L2c',aL2c,'L2c/F')
    bL3c = tree.Branch('L3c',aL3c,'L3c/F')
    bL4c = tree.Branch('L4c',aL4c,'L4c/F')
    bC1c = tree.Branch('C1c',aC1c,'C1c/F')
    bC2c = tree.Branch('C2c',aC2c,'C2c/F')
    bC3c = tree.Branch('C3c',aC3c,'C3c/F')
    bC4c = tree.Branch('C4c',aC4c,'C4c/F')
    bEt = tree.Branch('Et',aEt,'Et/F')
    bEvent = tree.Branch('Event',aEvent,'Event/F')
    bXh = tree.Branch('Xh',aXh,'Xh/F')
    bYh = tree.Branch('Yh',aYh,'Yh/F')
    bXp = tree.Branch('Xp',aXp,'Xp/F')
    bYp = tree.Branch('Yp',aYp,'Yp/F')
    bCost = tree.Branch('Cost',aCost,'Cost/F')
    bnfev = tree.Branch('nfev',anfev,'nfev/F')
    bnjev = tree.Branch('njev',anjev,'njev/F')
    bXt = tree.Branch('Xt',aXt,'Xt/F')
    bYt = tree.Branch('Yt',aYt,'Yt/F')
    bYTL1 = tree.Branch('YTL1',aYTL1,'YTL1/F')
    bYTL2 = tree.Branch('YTL2',aYTL2,'YTL2/F')
    bYTL3 = tree.Branch('YTL3',aYTL3,'YTL3/F')
    bYTL4 = tree.Branch('YTL4',aYTL4,'YTL4/F')
    bYTC1 = tree.Branch('YTC1',aYTC1,'YTC1/F')
    bYTC2 = tree.Branch('YTC2',aYTC2,'YTC2/F')
    bYTC3 = tree.Branch('YTC3',aYTC3,'YTC3/F')
    bYTC4 = tree.Branch('YTC4',aYTC4,'YTC4/F')
    bYTL1_fit = tree.Branch('YTL1_fit',aYTL1_fit,'YTL1_fit/F')
    bYTL2_fit = tree.Branch('YTL2_fit',aYTL2_fit,'YTL2_fit/F')
    bYTL3_fit = tree.Branch('YTL3_fit',aYTL3_fit,'YTL3_fit/F')
    bYTL4_fit = tree.Branch('YTL4_fit',aYTL4_fit,'YTL4_fit/F')
    bYTC1_fit = tree.Branch('YTC1_fit',aYTC1_fit,'YTC1_fit/F')
    bYTC2_fit = tree.Branch('YTC2_fit',aYTC2_fit,'YTC2_fit/F')
    bYTC3_fit = tree.Branch('YTC3_fit',aYTC3_fit,'YTC3_fit/F')
    bYTC4_fit = tree.Branch('YTC4_fit',aYTC4_fit,'YTC4_fit/F')
    buncertx = tree.Branch('uncertx',auncertx,'uncertx/F')
    buncerty = tree.Branch('uncerty',auncerty,'uncerty/F')
     
    
    
    #Gets the variables event from event
    for entryNum in range(0, ntuple.GetEntries()):    
        ntuple.GetEntry(entryNum)
        
        
        print('-------------------------------')
        print(entryNum)
        
        Xh = getattr(ntuple, 'Xh')
        Yh = getattr(ntuple, 'Yh')
        
        #For some reason, the point (5,5) has the order of lines and columns swipped, so it's taken care of in this part of the program
        #In future experiments, this part may not be needed
        if Xh == 5 and Yh == 5:
            L1 = getattr(ntuple, 'C1')
            #print(L1)
            L2 = getattr(ntuple, 'C2')
            #print(L2)
            L3 = getattr(ntuple, 'C3')
            #print(L3)
            L4 = getattr(ntuple, 'C4')
            #print(L4)
            C1 = getattr(ntuple, 'L4')
            #print(C1)
            C2 = getattr(ntuple, 'L3')
            #print(C2)
            C3 = getattr(ntuple, 'L2')
            #print(C3)
            C4 = getattr(ntuple, 'L1')
            #print(C4)
            Et = getattr(ntuple, 'Et')
            Event = getattr(ntuple, 'Event')
            print(Xh,Yh)
            

        else:
            L1 = getattr(ntuple, 'C4')
            #print(L1)
            L2 = getattr(ntuple, 'C3')
            #print(L2)
            L3 = getattr(ntuple, 'C2')
            #print(L3)
            L4 = getattr(ntuple, 'C1')
            #print(L4)
            C1 = getattr(ntuple, 'L1')
            #print(C1)
            C2 = getattr(ntuple, 'L2')
            #print(C2)
            C3 = getattr(ntuple, 'L3')
            #print(C3)
            C4 = getattr(ntuple, 'L4')
            #print(C4)
            Et = getattr(ntuple, 'Et')
            Event = getattr(ntuple, 'Event')
            print(Xh,Yh)
           
            
        calib(L1,L2,L3,L4,C1,C2,C3,C4)
        #################################### 
        #teste = Ytable(7.,-3.)
        #print("TESTE")
        #################################### 
       
        global newExN
        
        #Organizes the calibrated values from each line and column on an array holding two vectors (Lines and Columns)
        newExN = np.array([[calb[0],calb[1],calb[2],calb[3]],[calb[4],calb[5],calb[6],calb[7]]])
        
        #Associates the values on the right to the pointers on the left
        aL1[0] = L1
        aL2[0] = L2
        aL3[0] = L3
        aL4[0] = L4
        aC1[0] = C1
        aC2[0] = C2
        aC3[0] = C3
        aC4[0] = C4
        aL1c[0] = newExN[0][0]
        aL2c[0] = newExN[0][1]
        aL3c[0] = newExN[0][2]
        aL4c[0] = newExN[0][3]
        aC1c[0] = newExN[1][0]
        aC2c[0] = newExN[1][1]
        aC3c[0] = newExN[1][2]
        aC4c[0] = newExN[1][3]
        aEvent[0] = Event
        aXh[0] = Xh
        aYh[0] = Yh
        aXt[0] = Xh - 8.3
        aYt[0] = Yh - 8.3
        aEt[0] = Et
        
        #teste = Ytable(7.,-3.)
        #print("TESTE")
        
        #Executes the optimize least squared routine and stores the useful information on given pointers
        xypat = fitpatxy()
        
        
        YT = Ytable(aXt[0],aYt[0])
        aYTL1[0] = YT[0][0]
        aYTL2[0] = YT[0][1]
        aYTL3[0] = YT[0][2]
        aYTL4[0] = YT[0][3]
        aYTC1[0] = YT[1][0]
        aYTC2[0] = YT[1][1]
        aYTC3[0] = YT[1][2]
        aYTC4[0] = YT[1][3]
        auncertx[0] = unc[0]
        auncerty[0] = unc[1]
        print(aXt[0],aYt[0])
        print(YT)
        
        
        aXp[0] = result.x[0]
        aYp[0] = result.x[1]
        aCost[0] = result.cost
        anfev[0] = result.nfev
        anjev[0] = result.njev
        
        YT_fit = Ytable(aXp[0],aYp[0])
        aYTL1_fit[0] = YT_fit[0][0]
        aYTL2_fit[0] = YT_fit[0][1]
        aYTL3_fit[0] = YT_fit[0][2]
        aYTL4_fit[0] = YT_fit[0][3]
        aYTC1_fit[0] = YT_fit[1][0]
        aYTC2_fit[0] = YT_fit[1][1]
        aYTC3_fit[0] = YT_fit[1][2]
        aYTC4_fit[0] = YT_fit[1][3]
        print("YT_fit")
        print(YT_fit)
        
    
       
        #Fills the Tree with the data of the given event
        tree.Fill()
       
    
    
    tree.SetDirectory(0)
    tree.Fill()
    data.Close()
    print('Acabou')
    
    #Writes all data from every event on the Tree
    tree.Write()
    
    #Writes everything on the exit file and closes it
    outfile.Write()
    outfile.Close()
    
    #Congratulations! It's the end of the program
