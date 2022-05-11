import ROOT
from array import array

def plot(entryNum):

    file = ROOT.TFile.Open('/home/casarin/Desktop/rootdata/fit_2000_entries_per_point_newB.root','read')
    #file = ROOT.TFile.Open('/home/casarin/Desktop/rootdata/novotesteconsb.root','read')
    
    t = file.tree
    t.GetEntry(entryNum)
    L1 = getattr(t,'L1')
    L2 = getattr(t,'L2')
    L3 = getattr(t,'L3')
    L4 = getattr(t,'L4')
    
    C1 = getattr(t,'C1')
    C2 = getattr(t,'C2')
    C3 = getattr(t,'C3')
    C4 = getattr(t,'C4')
    
    L1c = getattr(t,'L1c')
    L2c = getattr(t,'L2c')
    L3c = getattr(t,'L3c')
    L4c = getattr(t,'L4c')
    
    C1c = getattr(t,'C1c')
    C2c = getattr(t,'C2c')
    C3c = getattr(t,'C3c')
    C4c = getattr(t,'C4c')
    
    YTL1 = getattr(t,'YTL1')
    YTL2 = getattr(t,'YTL2')
    YTL3 = getattr(t,'YTL3')
    YTL4 = getattr(t,'YTL4')
    
    YTC1 = getattr(t,'YTC1')
    YTC2 = getattr(t,'YTC2')
    YTC3 = getattr(t,'YTC3')
    YTC4 = getattr(t,'YTC4')
    
    YTL1_fit = getattr(t,'YTL1_fit')
    YTL2_fit = getattr(t,'YTL2_fit')
    YTL3_fit = getattr(t,'YTL3_fit')
    YTL4_fit = getattr(t,'YTL4_fit')
    
    YTC1_fit = getattr(t,'YTC1_fit')
    YTC2_fit = getattr(t,'YTC2_fit')
    YTC3_fit = getattr(t,'YTC3_fit')
    YTC4_fit = getattr(t,'YTC4_fit')
    
    Xh = getattr(t,'Xh')
    Yh = getattr(t,'Yh')
    
    tot1 = L1+L2+L3+L4
    tot2 = C1+C2+C3+C4
    title = str("(Xh,Yh)="+"("+str(Xh)+","+str(Yh)+")")
    en = str("entry"+str(entryNum))
    c1 = ROOT.TCanvas(en,title,1000,500)
    #c1.SeTitle(Xh,',',Yh)
    yl = array('f',[L1/tot1,L2/tot1,L3/tot1,L4/tot1])
    yc = array('f',[C1/tot2,C2/tot2,C3/tot2,C4/tot2])
    ytl = array('f',[YTL1,YTL2,YTL3,YTL4])
    ytc = array('f',[YTC1,YTC2,YTC3,YTC4])
    ytlf = array('f',[YTL1_fit,YTL2_fit,YTL3_fit,YTL4_fit])
    ytcf = array('f',[YTC1_fit,YTC2_fit,YTC3_fit,YTC4_fit])
    ylc = array('f',[L1c,L2c,L3c,L4c])
    ycc = array('f',[C1c,C2c,C3c,C4c])
    pitch = array('f',[-6.3,-2.1,2.1,6.3])
	
	
    pad1 = ROOT.TPad('pad1','This is pad1',0,0,0.48,1,0)
    pad2 = ROOT.TPad('pad2','This is pad2',0.52,0,1,1,0)
	
    pad1.SetRightMargin(0.09);  
    pad1.SetLeftMargin(0.15);
    pad1.SetBottomMargin(0.15);
	
    pad2.SetRightMargin(0.09);
    pad2.SetLeftMargin(0.15);
    pad2.SetBottomMargin(0.15);
    	
    pad1.Draw()
    pad2.Draw()
	
	###################################################
	#Dados Experimentais das Linhas
	
    pad1.cd()
    
    gr11 = ROOT.TGraph(4,pitch,ytl)
    ax11 = gr11.GetYaxis()
    ax11.SetRangeUser(0.,0.5)
    ax11.SetTitle("Yield normalizado")
    ax12 = gr11.GetXaxis()
    ax12.SetRangeUser(-7.,7.)
    ax12.SetTitle("Coordenadas y dos centros das linhas (mm)")
    gr11.SetMarkerStyle(20)
    gr11.SetMarkerSize(1)
    gr11.SetMarkerColor(2)
    gr11.SetLineColor(2)
    gr11.SetTitle('Pattern das linhas')
    
    gr11.Draw('APL')
    
    gr1 = ROOT.TGraph(4,pitch,yl)
    gr1.SetMarkerStyle(20)
    gr1.SetMarkerSize(1)
    gr1.Draw('PL')
    
    gr111 = ROOT.TGraph(4,pitch,ytlf)
    gr111.SetMarkerStyle(20)
    gr111.SetMarkerSize(1)
    gr111.SetMarkerColor(3)
    gr111.SetLineColor(3)
    gr111.Draw('PL')
    
    gr1111 = ROOT.TGraph(4,pitch,ylc)
    gr1111.SetMarkerStyle(20)
    gr1111.SetMarkerSize(1)
    gr1111.SetMarkerColor(4)
    gr1111.SetLineColor(4)
    gr1111.Draw('PL')
	
    
    c1.Update()
	
	
	
	###################################################
	#Dados Experimentais das Colunas
	
    pad2.cd()
    
    gr21 = ROOT.TGraph(4,pitch,ytc)
    ax21 = gr21.GetYaxis()
    ax21.SetRangeUser(0.,0.5)
    ax21.SetTitle("Yield normalizado")
    ax22 = gr21.GetXaxis()
    ax22.SetRangeUser(-7.,7.)
    ax22.SetTitle("Coordenadas x dos centros das colunas (mm)")
    
    ax21.Draw()
    gr21.SetMarkerStyle(20)
    gr21.SetMarkerSize(1)
    gr21.SetMarkerColor(2)
    gr21.SetLineColor(2)
    gr21.SetTitle('Pattern das colunas')
    
    gr21.Draw('APL')
    
    gr2 = ROOT.TGraph(4,pitch,yc)
    gr2.SetMarkerStyle(20)
    gr2.SetMarkerSize(1)
    gr2.Draw("PL")
    
    
    gr122 = ROOT.TGraph(4,pitch,ytcf)
    gr122.SetMarkerStyle(20)
    gr122.SetMarkerSize(1)
    gr122.SetMarkerColor(3)
    gr122.SetLineColor(3)
    gr122.Draw('PL')
    
    gr1222 = ROOT.TGraph(4,pitch,ycc)
    gr1222.SetMarkerStyle(20)
    gr1222.SetMarkerSize(1)
    gr1222.SetMarkerColor(4)
    gr1222.SetLineColor(4)
    gr1222.Draw('PL')
	
    c1.Update()
	
    c1.Draw()
    c1.Print()
    print('Imagem salva no diretório')
