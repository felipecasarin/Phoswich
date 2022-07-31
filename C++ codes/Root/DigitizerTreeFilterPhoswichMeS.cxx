#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TVector3.h"
#include <TFile.h>
#include <TNtuple.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <vector>

void FilterCoinc(const char* nameFile0="SDataR_run_7.root", const char* nameFileOut="PhoswichMeS.root") {

	
	TFile *myFile0 = TFile::Open(nameFile0);
	if(myFile0->IsZombie()){cout<<"Problem opening file 0"<<endl; return;}
	TFile outputFile (nameFileOut,"RECREATE");
	if(outputFile.IsZombie()){cout<<"Problem opening output file"<<endl; return;}
	TTree *tree= new TTree("tree","tree"); // output tree
	
	Double_t El0;
	Double_t Es0;
	Double_t Elp[16];
	Double_t Esp[16];
	Int_t Ch[16];
	Int_t NCh;
	unsigned long long int T0;
	Long_t Tdif[16];
	
	tree->Branch("NCh",&NCh,"NCh/I");	
	tree->Branch("El0",&El0,"El0/D");
	tree->Branch("Es0",&Es0,"Es0/D");
	tree->Branch("Elp",Elp,"Elp/D");
	tree->Branch("Esp",Esp,"Esp/D");
	tree->Branch("T0",&T0,"T0/L");
	tree->Branch("Tdif",Tdif,"Tdif/L");
	tree->Branch("Ch",Ch,"Ch/I");


   TTreeReader myReader0("Data_R", myFile0);
   TTreeReaderValue<UShort_t> myCh0(myReader0, "Channel");
   TTreeReaderValue<unsigned long long int> myTst0(myReader0, "Timestamp");
   TTreeReaderValue<UShort_t> myBoard0(myReader0, "Board");
   TTreeReaderValue<UShort_t> myEn0(myReader0, "Energy");
   TTreeReaderValue<UShort_t> myEns0(myReader0, "EnergyShort");
   TTreeReaderValue<UInt_t> myFlag0(myReader0, "Flags");


   myReader0.Restart();

   // Loop over all entries of the trees

   NCh=0;
   myReader0.Next();
        if (*myBoard0==0 && *myCh0==0)
     {
      T0=(unsigned long long int)*myTst0;
      El0=(Double_t)*myEn0;
      Es0=(Double_t)*myEns0;
      cout<<"Começou!!!"<<endl;
     }
     else
      cout<<"Começou errado!!!"<<endl;
     


    while(myReader0.Next())
    {
     if(*myBoard0==1){
      Elp[NCh]=(Double_t)*myEn0;
      Esp[NCh]=(Double_t)*myEns0;
      Tdif[NCh]=(unsigned long long int)*myTst0-T0;
      Ch[NCh]=(Int_t)*myCh0;
      NCh++;
     }
     if(*myBoard0==0 && *myCh0==0){
     tree->Fill();
      T0=(unsigned long long int)*myTst0;
      El0=(Double_t)*myEn0;
      Es0=(Double_t)*myEns0;
      NCh=0;
      }
    }
      tree->Fill();   
	tree->Write();

	outputFile.Close();
	myFile0->Close();

	return;
} 
