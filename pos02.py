#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import ROOT
from ROOT import TNtuple
import sys

print('Programa rodando')
Nmax = 10


def pos0():
    
    data = ROOT.TFile.Open('/home/casarin/Desktop/rootdata/lessdata.root', 'read')
    
    ntuple = data.Get('ntuple')
    
    test3 = ROOT.TFile.Open( '/home/casarin/Desktop/rootdata/testepequeno.root', 'RECREATE', 'ROOT file with an NTuple' )
    ntupleout  = TNtuple( 'ntuple','ntuple','L1:L2:L3:L4:C1:C2:C3:C4:Xh:Yh:Et:Event')

 
    for entryNum in range(0, ntuple.GetEntries()):
        ntuple.GetEntry(entryNum)
        
        print(entryNum)
        L1 = getattr(ntuple, 'L1')
        L2 = getattr(ntuple, 'L2')
        L3 = getattr(ntuple, 'L3')
        L4 = getattr(ntuple, 'L4')
        C1 = getattr(ntuple, 'C1')
        C2 = getattr(ntuple, 'C2')
        C3 = getattr(ntuple, 'C3')
        C4 = getattr(ntuple, 'C4')
        Xh = getattr(ntuple, 'Xh')
        Yh = getattr(ntuple, 'Yh')
        Et = getattr(ntuple, 'Et')
        Event = getattr(ntuple, 'Event')
    
        #Manual cut
        if Event < Nmax and (C1+C2+C3+C4 < 1500) and (C1+C2+C3+C4 > 650) and Et > 550 and Et < 1800:
          print('ok')
          ntupleout.Fill(L1,L2,L3,L4,C1,C2,C3,C4,Xh,Yh,Et,Event)
    
    
    ntupleout.SetDirectory(0)
    
    data.Close()
    

    print('Acabou')
    ntupleout.Write()
    test3.Write()
    test3.Close()
