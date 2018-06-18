#ifndef Axion_hh
#define Axion_hh 1


#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <ctime>
#include <algorithm>
#include <cmath>


#include "TMath.h"
#include "TGraph.h"
#include "TGraph2DErrors.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"
#include "TRandom3.h"

#include "AxionFunctions.hh"


using namespace std;

struct Axion{

  TGraph *gCS;
  TGraph *gFlux;
  TGraph *gCSall;
  TGraph *gM;
  TGraph *gNeE ;
  TGraph *gENe;

  TH1F *hGalactic;
  TH1F *hGalactic_smeared;
  TH1F *hSolar;
  TH1F *hSolar_smeared;

  TRandom3 *rnd;
  
  double OccCore;
  double OccCentralPMT;


};




void initAxion(Axion& ax){

  ax.OccCore       = 0.4310;
  ax.OccCentralPMT = 0.0620;
  
  TFile *fne = new TFile("escale.root");
  ax.gNeE = (TGraph*) fne->Get("gNeE");
  fne->Close();


  ax.rnd   = new TRandom3(1234);

  ax.gCS   = new TGraph("data/cs_photoelectric.dat");
  ax.gFlux = new TGraph("data/gae_flux.dat");
  
  ax.hGalactic         = new TH1F("hGalactic","hGalactic",1001,-0.5,200.5);
  ax.hGalactic_smeared = new TH1F("hGalactic_smeared","hGalactic (with Ne resolution)",1001,-0.5,200.5);

  ax.hSolar          = new TH1F("hSolar","hSolar",501,-0.5,100.5);
  //ax.hSolar          = new TH1F("hSolar","hSolar",200,0,200); 
  ax.hSolar_smeared  = new TH1F("hSolar_smeared","hSolar (with Ne resolution)",501,-0.5,100.5);

}

void resetAxion(Axion& ax){

  ax.hGalactic->Reset();
  ax.hGalactic_smeared->Reset();

  ax.hSolar->Reset();
  ax.hSolar_smeared->Reset();
}
  

double GetNeRMS(bool IsCentralPMT, int Ne )
{
  double rms;
  if (IsCentralPMT )
    {
      if ( Ne==1 )     rms=0.22;
      else if (Ne==2 ) rms=0.35;
      else if (Ne==3 ) rms=0.44;
      else             rms=0.44+(double(Ne)-3)*0.0423;
    }
  else                 rms=0.29+(double(Ne)-1.)*0.17;
  return rms;
}



void GenerateSolarAxionSpec(Axion &ax, double gAe){


  TF1 *f1 = new TF1("f1",fSolar,0,100,2);
  f1->SetParameters(gAe,0);
  f1->SetNpx(1000);

  cout << f1->Integral(0.5,100.5) <<endl;
  //f1->Draw();

  bool IsCentral = false;

  for(int i=0; i< ax.hSolar->GetNbinsX(); i++){

    double Ne = ax.hSolar->GetBinCenter(i);
    double nsig = f1->Eval(Ne);
   
    ax.hSolar->SetBinContent(i, nsig);
  }

  ax.hSolar->Rebin(5);

  for(int i=0; i<100000;++i){

    double Ne  = ax.hSolar->GetRandom();
    double E   = ax.gNeE->Eval(Ne);

    if(ax.rnd->Rndm() < ax.OccCentralPMT/ax.OccCore) IsCentral = true;
   
    double Ne_bin = ax.rnd->Binomial((int)(E/0.0195), Ne/(E/0.0195));

    double rms  = GetNeRMS(IsCentral,Ne_bin);
    
    double Ne_f = ax.rnd->Gaus(Ne_bin, 0.2*sqrt(Ne_bin));
    
    ax.hSolar_smeared->Fill(Ne_f);

  }
  
  ax.hSolar_smeared->Rebin(5);


  float n1 = ax.hSolar->Integral();
  float n2 = ax.hSolar_smeared->Integral();

  float scale = n1/n2;

  for(int i=0; i<=ax.hSolar_smeared->GetNbinsX(); ++i){

    float A   = ax.hSolar_smeared->GetBinContent(i);
    float err = ax.hSolar_smeared->GetBinError(i);

    ax.hSolar_smeared->SetBinContent(i, A * scale);
    ax.hSolar_smeared->SetBinError(i, err * scale);
  }

}


void GenerateGalacticAxionSpec(Axion &ax, double mA){

  TF1 *f1 = new TF1("f1",fGalactic,0,200,2);
  f1->SetParameters(1e-13,mA);


  for(int i=0; i< ax.hGalactic->GetNbinsX(); i++){

    double ene = ax.hGalactic->GetBinCenter(i);
    double nsig = f1->Eval(ene);
  
    ax.hGalactic->SetBinContent(i, nsig);

  }
  
  bool IsCentral = false;

  for(int i=0; i<100000;++i){

    double Ne  = ax.hGalactic->GetRandom();
    double E   = ax.gNeE->Eval(Ne);

    if(ax.rnd->Rndm() < ax.OccCentralPMT/ax.OccCore) IsCentral = true;

    double Ne_bin = ax.rnd->Binomial((int)(E/0.0195), Ne/(E/0.0195));
    
    double rms  = GetNeRMS(IsCentral,Ne_bin);

    double Ne_f = ax.rnd->Gaus(Ne_bin,0.2*Ne_bin);
    ax.hGalactic_smeared->Fill(Ne_f);

  }

  float n1 = ax.hGalactic->Integral();
  float n2 = ax.hGalactic_smeared->Integral();

  float scale = n1/n2;

  for(int i=0; i<=ax.hGalactic_smeared->GetNbinsX(); ++i){

    float A   = ax.hGalactic_smeared->GetBinContent(i);
    float err = ax.hGalactic_smeared->GetBinError(i);
    
    ax.hGalactic_smeared->SetBinContent(i, A * scale);
    ax.hGalactic_smeared->SetBinError(i, err * scale);
  }

}

#endif
