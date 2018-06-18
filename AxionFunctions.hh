#ifndef AxionFunctions_hh
#define AxionFunctions_hh 1

#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <ctime>
#include <algorithm>
#include <cmath>


#include "TMath.h"
#include "TGraph2DErrors.h"
#include "TF1.h"
#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"


#include "Axion.hh"

using namespace std;


double alpha   = 1./137.035 ;
double me      = 510.9989461; // keV
double density = 1.3954 ; // g/cm^3 
double day     = 24*3600; 
double barn    = 1e-24; //cm2
double kg      = 6.022e23/40.*1000; // # atomi / kg

TGraph *gCS = new TGraph("/Users/anavrera/DarkSide/Axions/data/cs_photoelectric.dat");
TGraph *gFlux = new TGraph("/Users/anavrera/DarkSide/Axions/data/gae_flux.dat");
TGraph *gENe = new TGraph();

double getNe(double Eee) {
  return TMath::Log(3.765*Eee + 1)*9.33093e-01*(1.69420e+01+Eee*1.43766e+00);

}

double getNe_new(double Eee) {
  return TMath::Log(3.76531*Eee + 1)*9.33093e-01*(1.69420e+01+Eee*1.43766e+00);

}


void init_gENe(){

  int nbins = 2000;
  double step = 40./nbins ;
  for(int i=0;i<nbins;++i) {
    double E = i*step + step/3. ;
    double Ne  = getNe(E);
    gENe->SetPoint(i,Ne,E);
  } 


}


double fGalacticFlux(double *x, double *p) {
  double E   = x[0];
  double mA  = p[0];

  double     beta = 1 ; 
  if(mA > 0) beta = sqrt(2*E/mA);
  
  return 9e15*beta / mA;

}


double fSolar(double *x, double *p) {
  
 
  double E   = gENe->Eval(x[0]);
  double gAe = p[0];
  double mA  = p[1];
  
  double     beta = 1 ; 
  if(mA > 0) beta = sqrt(2*E/mA);
     
  double cs =   pow(gAe,2)/beta
    * 3*pow(E,2)/(16.*TMath::Pi()*alpha*pow(me,2))
    * (1-pow(beta,2./3)/3.);
  
  //cout << E << " " << cs << endl;
  
  cs *= gCS->Eval(E);//*barn;
  cs *= gFlux->Eval(E)*pow(gAe,2) ; // cm2 / (cm2 keV s)
  cs *= day ; 
  cs *= kg ; 
  
  //cout << E << " " << cs << endl;

  if(cs < 0) return 0 ;
  return cs;
 
}




double fGalactic(double *x, double *p) {
  
  init_gENe();
  
  double E   = gENe->Eval(x[0]);
  double gAe = p[0];
  double mA  = p[1];
  
  if(E < mA) return 0 ;
  double     beta = 1 ; 
  if(mA > 0) beta = sqrt(2*E/mA);
     
  double cs =   pow(gAe,2)/beta
    * 3*pow(E,2)/(16.*TMath::Pi()*alpha*pow(me,2))
    * (1-pow(beta,2./3)/3.);
   
  cs *= fGalacticFlux(&E,&mA);
  cs *= gCS->Eval(E) ; // cm2 / (cm2 keV s)
  cs *= day ; 
  cs *= kg ; 
  if(cs < 0) return 0 ;
  return cs;
 
}

#endif
