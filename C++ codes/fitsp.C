#include "Riostream.h"

double pi=acos(-1.);

double ECompton(double E0, double theta){
	double M=511.;
	double E = E0/(1.+E0/M*(1-cos(theta)));
	return E;
}

double WS(double R0,double A0, double x){
	double w=1./(1.+exp((x-R0)/A0));
	//if((x-R0)/A0>30.)w=0.;
	return w;
}

double WS3(double R0,double A0, double B0, double x){
	double w=(1.+pow(x/B0,2))/(1.+exp((x-R0)/A0));
	//if((x-R0)/A0>30.)w=0.;
	return w;
}

double bgy1(Double_t *x,Double_t *p) {
	//p[0]=yi, p[1]=yf (ajustáveis), p[2]=xi, p[3]=xf (fixos)
	return p[0]+(p[1]-p[0])/(p[3]-p[2])*(x[0]-p[2]);
   }  

double bgf(Double_t *x,Double_t *par) {//6 parameters
	//[0]:b0, [1]:b1, [2]: b2, [3]:fwhm1, [4]:Area1, [5]:Centr
	double b0=par[0];
	double b1=par[1];
	double b2=par[2];
	
	double bg=b0+b1*x[0]+b2*x[0]*x[0];
	return bg;
   }  
   
double f1p(Double_t *x,Double_t *par) {//7 parameters
	//[0]:yi, [1]:yf, [2]:xi, [3]:xf, [4]:fwhm1, [5]:Area1, [6]:Centr


	double fwhm1=par[4];
	double Area1=par[5];
	double sig1=fwhm1/2.355;// standard deviations
	double Centr=par[6];

	double arg1=1.;
	
	if (sig1!=0) {
		arg1 = -0.5*pow((x[0] - Centr)/sig1,2);
	}
	double f = Area1/sqrt(2.*pi)/sig1*exp(arg1)+bgy1(x,par);
	return f;
   }


TH1D * GetLastHist(){
   TList *l=gPad->GetListOfPrimitives();
   Int_t ns=l->GetSize();
   int found=0;
   TH1D *hist=0;
   for(Int_t n=1;n<ns;n++){
   TObject *tob=l->At(n);
   if(tob->InheritsFrom("TH1")){
    hist =(TH1D*)tob;
    found++;
        }
   }
   if(found==0)cout<<"f1p - no Histogram found"<<endl;
   return hist;
} 


TFitResultPtr Fit1p(double xi, double xf){

   TH1D *hist=GetLastHist();

	//double xf=gPad->GetUxmax();
	TF1 *func=new TF1("fit",f1p,xi,xf,7);
	TF1 *bg=new TF1("bgf",bgy1,xi,xf,4);
	//TH1D *hist = (TH1D*)gDirectory->Get(histname.c_str());

	Int_t bini=hist->FindBin(xi);
	Int_t binf=hist->FindBin(xf);
	hist->GetXaxis()->SetRange(bini,binf);
	hist->Sumw2();
	TH1D *histcopy=(TH1D*)hist->Clone();
	histcopy->SetName("histcopy");
	Int_t b=hist->GetMaximumBin();
	double dbin=binf-bini;
	//double xm=(xf+xi)/2.;
	double xm=hist->GetBinCenter(b);
	double dx=xf-xi;
	double yi=hist->GetBinContent(bini);//hist->GetBinWidth(bini);
	double yf=hist->GetBinContent(binf);//hist->GetBinWidth(binf);
	double ym=(yi+yf)/2.;
	double b1=(yf-yi)/dx;
	double b0=yi-b1*xi;	
	double b2=0.;
	double Ab=dx*ym;
	double Atot=hist->Integral(bini,binf);
	double A=0;//(Atot-Ab);
	double Centr=xm;
	func->SetParNames("yi","yf","xi","xf","fwhm","Area1*binw", "Centr");
	func->SetParameter(0,yi);
	func->SetParLimits(0,0.,5.*yi+0.1);
	func->SetParameter(1,yf);
	func->SetParLimits(1,0.,5.*yf+0.1);
	func->SetParameter(2,xi);
	func->SetParLimits(2,xi,xi);
	func->FixParameter(2,xi);
	func->SetParameter(3,xf);
	func->SetParLimits(3,xf,xf);
	func->FixParameter(3,xf);

	bg->SetParameter(0,yi);
	bg->SetParLimits(0,0.,10.*yi);
	bg->SetParameter(1,yf);
	bg->SetParLimits(1,0.,10.*yf);
	bg->SetParameter(2,xi);
	bg->SetParLimits(2,xi,xi);
	bg->FixParameter(2,xi);
	bg->SetParameter(3,xf);
	bg->SetParLimits(3,xf,xf);
	bg->FixParameter(3,xf);
	histcopy->GetXaxis()->SetRange(bini,binf);
	histcopy->Add(bg,-1.);
	double west = 2.355*histcopy->GetStdDev()+hist->GetBinWidth(b)/2.;

	func->SetParameter(4,west);
	func->SetParLimits(4,0.05*west,2.*west);
	func->SetParameter(5,A);
	func->SetParLimits(5,A/10.,10.*A);
	func->SetParameter(6,xm);
	func->SetParLimits(6,xi,xf);
	func->SetLineColor(3);
	
	TFitResultPtr r=hist->Fit("fit","RSB");
	TMatrixDSym cov = r->GetCovarianceMatrix();   //  to access the covariance matrix
	Double_t chi2   = r->Chi2();                  // to retrieve the fit chi2
	/*
	Double_t e0   = r->Parameter(0);            // retrieve the value for the parameter 0
	Double_t e1   = r->Parameter(1);            // retrieve the value for the parameter 0
	Double_t err0   = r->ParError(0);             // retrieve the error for the parameter 0
	Double_t err1   = r->ParError(1);             // retrieve the error for the parameter 0
	*/
	//r->Print("V");           // print full information of fit including covariance matrix
	r->Print();    
	double areacounts=r->Parameter(5)/hist->GetBinWidth(hist->FindBin(xm))  ;          
	cout<<"Area="<<areacounts<<endl; 
	//TFile *fres = new TFile("fitresults.root","recreate");
	//fres->cd();
	//r->Write();
	return r;                                   // store the result in a file
	
}

TFitResultPtr Fit1p(double xp){
    TH1D *hist=GetLastHist();
    Int_t peakbin=hist->FindBin(xp);
    Double_t binw=hist->GetBinWidth(peakbin);// assuming constant
    
    double minder=0.;
    double b=peakbin;
    
    while(hist->GetBinContent(b+2)-hist->GetBinContent(b)<minder){
		if(hist->GetBinContent(b+1)>hist->GetBinContent(b)+sqrt(hist->GetBinContent(b)))peakbin=b+1;
		minder=hist->GetBinContent(b+2)-hist->GetBinContent(b);
		b++;
		cout<<b<<" minder="<<minder<<" ";
	}
	double hwr=2.355*(b-peakbin+0.5)*binw/2.;
	cout<<"peakbin="<<peakbin<<" bmind="<<b+1<<" hwr="<<hwr<<endl;

    double maxder=0.;
    b=peakbin;
    while(hist->GetBinContent(b)-hist->GetBinContent(b-2)>maxder){
		if(hist->GetBinContent(b-1)>hist->GetBinContent(b)+sqrt(hist->GetBinContent(b)))peakbin=b-1;
		maxder=hist->GetBinContent(b)-hist->GetBinContent(b-2);
		b--;
		cout<<b<<" maxder="<<maxder<<" ";
	}
	double hwl=2.355*(peakbin-b+0.5)*binw/2.;
	cout<<"peakbin="<<peakbin<<" bmaxd="<<b-1<<" hwl="<<hwl<<endl;
	xp=hist->GetBinCenter(peakbin);
	double fwest=(hwl+hwr);
	double xi=xp-3.*fwest;
	double xf=xp+3.*fwest;
	cout<<"xi,xf="<<xi<<","<<xf<<endl;
	TFitResultPtr r=Fit1p(xi,xf);
	double fwhm=r->Parameter(4);
	xp=r->Parameter(6);
	xi=xp-3.*fwhm;
	xf=xp+3.*fwhm;
	cout<<"xi,xf="<<xi<<","<<xf<<endl;
	r=Fit1p(xi,xf);
	while( abs(r->Parameter(4)-fwhm)/r->Parameter(4) > r->Parameter(4)/4.){	
		xp=r->Parameter(6);
		fwhm=r->Parameter(4);
		xi=xp-3.*fwhm;
		xf=xp+3.*fwhm;
		cout<<"xi,xf="<<xi<<","<<xf<<endl;
		r=Fit1p(xi,xf);
	}
	
	/*
			xp=r->Parameter(5);
			xi=xp-2.*fwhm;
			xf=xp+2.*fwhm;
			r=Fit1p(xi,xf);
			cout<<"xi,xf="<<xi<<","<<xf<<endl;
	*/	
	return r;
}

double Fit1p(){
	
   TList *l=gPad->GetListOfPrimitives();
	double xi=gPad->GetUxmin();
	double xf=gPad->GetUxmax();
	TFitResultPtr r=Fit1p(xi,xf);
	return r;
}

double fCo60ws(Double_t *x,Double_t *par) {//9 parameters
	//[0]:e0, [1]:e1, [2]:fwhm1, [3]:fwhm2, [4]:ws1, [5]:ws2, [6]:Area1, [7]:Area2, [8]:bg0	

	double e0=par[0];
	double e1=par[1]; //keV/channel
	double fwhm1=par[2];
	double fwhm2=par[3];
	double ws1=par[4];
	double ws2=par[5];
	double Area1=par[6];
	double Area2=par[7];
	double bg0=par[8];//counts/keV

	double sig1=fwhm1/2.355;// standard deviations
	double sig2=fwhm2/2.355;
	double E=e0+e1*x[0]; // linear calibration
	double E1=1173.2;// energies in keV
	double E2=1332.5;
	double Elim1=E1-ECompton(E1,pi);
	double Elim2=E2-ECompton(E2,pi);
	
	//cout<<"Elims=" <<Elim1<<" "<<Elim2<<endl;
	
	double bg=(ws1*WS(Elim1,fwhm1/2.355,E)+ws2*WS(Elim2,fwhm2/2.355,E)+bg0);
	
	
	double arg1=1.;
	double arg2=1.;
	

	
	if ((sig1*sig2)!=0) {
		arg1 = -0.5*pow((E - E1)/sig1,2);
		arg2 = -0.5*pow((E - E2)/sig2,2);
	}
	double fitval = e1*(Area1/sqrt(2.*pi)/sig1*exp(arg1)+Area2/sqrt(2.*pi)/sig2*exp(arg2)+bg);//counts/channel
	return fitval;
   }


double FitCo60ws(double xi, double xf){
	TH1D *hist = GetLastHist();	
	TF1 *func=new TF1("fit",fCo60ws,xi,xf,9);
	double dx=xf-xi;
	double slest=2505./(xi+xf);
	double west = slest*dx/4.;
	double yi=hist->GetBinContent(xi);
	double yf=hist->GetBinContent(xf);
	double bs=(yf-yi)/(xf-xi);
	double b0=yf-bs*xf;
	double Atot=hist->Integral(xi,xf);
	double Ab=(xf-xi)*(yf+yi)/2.;
	double A=(Atot-Ab)/2.;
	func->SetParNames("e0(keV)","e1(keV/ch)","fwhm1(keV)","fwhm2(keV)","ws1","ws2","Area1","Area2","bg0");
	func->SetParameter(0,0.);
	func->SetParLimits(0,-2.*slest*dx,2.*slest*dx);
	func->SetParameter(1,slest);
	func->SetParLimits(1,0.5*slest,1.5*slest);
	func->SetParameter(2,west);
	func->SetParLimits(2,west/100.,west*3.5);
	func->SetParameter(3,west);
	func->SetParLimits(3,west/100.,west*3.5);
	func->SetParameter(4,A/200.);
	func->SetParLimits(4,A/2000.,10.*A);
	func->SetParameter(5,A/200.);
	func->SetParLimits(5,A/2000.,10.*A);
	func->SetParameter(6,A);
	func->SetParLimits(6,A/4.,4.*A);
	func->SetParameter(7,A);
	func->SetParLimits(7,A/4.,4.*A);
	func->SetParameter(8,yf);
	func->SetParLimits(8,0.,5.*yf);
	
	TFitResultPtr r=hist->Fit("fit","RS");
	TMatrixDSym cov = r->GetCovarianceMatrix();   //  to access the covariance matrix
	Double_t chi2   = r->Chi2();                  // to retrieve the fit chi2
	Double_t e0   = r->Parameter(0);            // retrieve the value for the parameter 0
	Double_t e1   = r->Parameter(1);            // retrieve the value for the parameter 0
	Double_t err0   = r->ParError(0);             // retrieve the error for the parameter 0
	Double_t err1   = r->ParError(1);             // retrieve the error for the parameter 0
	//r->Print("V");           // print full information of fit including covariance matrix
	r->Print();                                
	double E1=1173.2;// energies in keV
	double E2=1332.5;
	double ch1=(E1-e0)/e1;
	double ch2=(E2-e0)/e1;
	double s1=sqrt( pow(err0/e1,2)+pow(err1*ch1/e1,2) );
	double s2=sqrt( pow(err0/e1,2)+pow(err1*ch2/e1,2) );
	cout<<"chan1="<<ch1<<"("<<s1<<") chan2="<<ch2<<"("<<s2<<")"<<endl;
	//TFile *fres = new TFile("fitresults.root","recreate");
	//fres->cd();
	//r->Write();                                   // store the result in a file
	return 1.;
}


double fCs137ws(Double_t *x,Double_t *par) {//8 parameters
	//[0]:e0, [1]:e1, [2]:fwhm1, [3]:aws, [4]:ws1, [5]:Area1, [6]:bg0, [7]:bws

	double e0=par[0];
	double e1=par[1];
	double fwhm1=par[2];
	double aws=par[3];
	double ws1=par[4];
	double Area1=par[5];
	double bg0=par[6];
	double bws=par[7];
	double sig1=fwhm1/2.355;// standard deviations
	double sig2=aws;
	double E=e0+e1*x[0]; // linear calibration
	double E1=661.657;// energies in keV
	double Elim1=E1-ECompton(E1,pi);// Compton edge of E1
	
	//cout<<"Elims=" <<Elim1<<" "<<Elim2<<endl;
	
	double bg=ws1*WS3(Elim1,aws,bws,E)+bg0;
	
	double arg1=1.;
	
	if (sig1!=0) {
		arg1 = -0.5*pow((E - E1)/sig1,2);
	}
	double fitval = e1*(Area1/sqrt(2.*pi)/sig1*exp(arg1)+bg);
	return fitval;
   }

   
double FitCo60ws(){
	double xi=gPad->GetUxmin();
	double xf=gPad->GetUxmax();
	return FitCo60ws(xi,xf);
}


  double FitCs137ws(double xi, double xf){
	TF1 *func=new TF1("fit",fCs137ws,xi,xf,8);
	TH1D *hist = GetLastHist();	
	double E1=661.657;// energies in keV
	double Elim1=E1-ECompton(E1,pi);
	double dx=xf-xi;
	double slest=1.2*(E1-Elim1)/dx;
	double west = slest*dx/4.;
	double yi=hist->GetBinContent(xi);
	double yf=hist->GetBinContent(xf);
	double b0=yf*slest;
	double Atot=hist->Integral(xi,xf);
	double Ab=(xf-xi)*(yf+yi)/2.5;
	double A=(Atot-Ab)/2.;
	double bws=10.*slest*xf;
	func->SetParNames("e0(keV)","e1(keV/ch)","fwhm1(keV)","aws(keV)","ws1","Area1","bg0","bws");
	func->SetParameter(0,0.);
	func->SetParLimits(0,-slest*dx/4.,slest*dx/4.);
	func->SetParameter(1,slest);
	func->SetParLimits(1,slest/4.,4.*slest);
	func->SetParameter(2,west);
	func->SetParLimits(2,west/100.,west*3.);
	func->SetParameter(3,west);
	func->SetParLimits(3,west/100.,west*3.);
	func->SetParameter(4,A/50.);
	func->SetParLimits(4,A/200.,A/2.);
	func->SetParameter(5,A);
	func->SetParLimits(5,A/5.,5.*A);
	func->SetParameter(6,b0);
	func->SetParLimits(6,0.,2.*b0);
	func->SetParameter(7,bws);
	func->SetParLimits(7,0.01*bws,1000.*bws);
	
	TFitResultPtr r=hist->Fit("fit","RS");
	TMatrixDSym cov = r->GetCovarianceMatrix();   //  to access the covariance matrix
	Double_t chi2   = r->Chi2();                  // to retrieve the fit chi2
	Double_t e0   = r->Parameter(0);            // retrieve the value for the parameter 0
	Double_t e1   = r->Parameter(1);            // retrieve the value for the parameter 0
	Double_t err0   = r->ParError(0);             // retrieve the error for the parameter 0
	Double_t err1   = r->ParError(1);             // retrieve the error for the parameter 0
	//r->Print("V");           // print full information of fit including covariance matrix
	r->Print();                                

	double ch1=(E1-e0)/e1;
	double ch2=(Elim1-e0)/e1;
	double s1=sqrt(err0*err0)+pow(err1*ch1/e1,2);
	double s2=sqrt(err0*err0)+pow(err1*ch2/e1,2);
	cout<<"chan1(peak)="<<ch1<<"("<<s1<<") chan2(edge)="<<ch2<<"("<<s2<<")"<<endl;
	//TFile *fres = new TFile("fitresults.root","recreate");
	//fres->cd();
	//r->Write();                                   // store the result in a file
	return 1.;
}

  double FitCs137ws(){ 
	double xi=gPad->GetUxmin();
	double xf=gPad->GetUxmax();
	return FitCs137ws(xi,xf);
}
