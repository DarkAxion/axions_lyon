#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>


#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TMath.h"
#include "TProfile.h"
#include "TMinuit.h"
#include "TGaxis.h"
#include "TKey.h"
#include "TROOT.h"

bool fix = 0;
double scalefit = 1.;
using namespace std;
using namespace TMath;

TH2F *hdata;
TH2F *hdata2;
TH2F *hcryo;
TH2F *hpmt;
TH2F *hint;

TH1D *hmodel ;




TGraph *gCS ;
TGraph *gFlux;
TGraph *gENe;

double fSolar(double *x, double *p) {
  double alpha   = 1./137.035 ;
  double me      = 510.9989461; // keV
  double density = 1.3954 ; // g/cm^3 
  double day     = 24*3600; 
  double barn    = 1e-24; //cm2
  double kg      = 6.022e23/40.*1000; // # atomi / kg
  double eff     = 1.;

  double exposure = 6786;
  
  double E   = gENe->Eval(x[0]);
  double gAe = p[0];
  double mA  = p[1];
  
  double     beta = 1 ; 
  if(mA > 0) beta = sqrt(2*E/mA);
     
  double cs =   pow(gAe,2)/beta
    * 3*pow(E,2)/(16.*TMath::Pi()*alpha*pow(me,2))
    * (1-pow(beta,2./3)/3.);
  
   
  cs *= gCS->Eval(E);//*barn;
  
  
  cs *= gFlux->Eval(E)*pow(gAe,2) ; // cm2 / (cm2 keV s)
  cs *= day ; 
  cs *= kg ; 
  cs *= eff ;
  cs *= exposure ;

  if(cs < 0) return 0 ;
  return cs;
 
}

//-------------------------------------------------------------------------
//    Build Model
//-------------------------------------------------------------------------

void buildModel() {
  gCS   = new TGraph("input/cs_photoelectric.dat");
  gFlux = new TGraph("input/gae_flux.dat");  

  TFile *fNe = TFile::Open("../MC/ERscale_Mar30.root");
  TGraph *gNeE  = (TGraph*) fNe->Get("gScale");
  gENe = new TGraph ;
  int N = gNeE->GetN();
  for(int i=0;i<N;++i) {
    double x,y ;
    gNeE->GetPoint(i,x,y);
    gENe->SetPoint(i,y,x);
  }
 
  int nbinsX = hdata->GetXaxis()->GetNbins();
  hmodel = new TH1D("hmodel","",
     hdata->GetXaxis()->GetNbins(),
     hdata->GetXaxis()->GetBinLowEdge(1),
		    hdata->GetXaxis()->GetBinLowEdge(nbinsX-1)+hdata->GetXaxis()->GetBinWidth(3)+1);

    
  TF1 *fun = new TF1("fun",fSolar,0,200,2);
  fun->SetParameters(3.5e-12,0);
  fun->Draw();
 
  
  for(int i=1;i<hmodel->GetXaxis()->GetNbins();++i) {
    float   ne   = hmodel->GetBinCenter(i);
    double  val = fun->Eval(ne);
    hmodel->SetBinContent(i,val);
  }
  
}
//-------------------------------------------------------------------------
//    Chi2 minimization
//-------------------------------------------------------------------------
double getChi2(double *par) {
  double logL = 0;

  //cout<<"hdata "<<hdata->GetXaxis()->GetNbins()<<endl;
  
  for(int ix=1;ix<=hdata->GetXaxis()->GetNbins();++ix) {
  
    double Ne    = hdata->GetXaxis()->GetBinCenter(ix);    
    if(Ne < 20)     continue ;
  
    double axion = hmodel->GetBinContent(ix);
    //cout<<" Ne "<<Ne<<endl;
    if(Ne > 100) {
      for(int iy=1;iy<hdata->GetYaxis()->GetNbins();++iy) {
	double td = hdata->GetYaxis()->GetBinCenter(iy);
	if(td < 100)      continue ;
	if(td > 250)       continue ;
	double data     = hdata->GetBinContent(ix,iy);
	double bg_int   = hint->GetBinContent(ix,iy);
	double bg_cryo  = hcryo->GetBinContent(ix,iy);
	double bg_pmt   = hpmt->GetBinContent(ix,iy);
	double width    = hdata->GetYaxis()->GetBinWidth(3);
        double signal   = axion/375.*width ;
	double model    = (par[0]*bg_cryo + par[1]*bg_pmt + par[2]*bg_int + par[3]*signal) ;
	
	if(model == 0)    continue ;
	//if(data  == 0)    continue ;

	//cout<<"2D fit "<<Ne<<" "<<td<<" "<<data<<" "<<model<<endl;
	
	//logL += pow(data-model,2)/(data) ;
	logL -= 2*TMath::Log(TMath::Poisson(data,model));
      }
    } else {
      double data     = 0 ;
      double bg_int   = 0;
      double bg_cryo  = 0;
      double bg_pmt   = 0;
      for(int iy=1;iy<hdata->GetYaxis()->GetNbins();++iy) {
        data     += hdata->GetBinContent(ix,iy);
        bg_int   += hint->GetBinContent(ix,iy);
        bg_cryo  += hcryo->GetBinContent(ix,iy);
        bg_pmt   += hpmt->GetBinContent(ix,iy);
      }
      double width    = hdata->GetYaxis()->GetBinWidth(3);
      double model    = (par[0]*bg_cryo + par[1]*bg_pmt + par[2]*bg_int + par[3]*axion) ;

      if(model <= 0)     continue ;
      //if(data  <= 0)     continue ;
      //logL += pow(data-model,2)/data ;
      //cout << Ne<< " 2   " <<  data << " " << model << " " <<  par[0]*bg_cryo<<" "<<par[1]*bg_pmt<<" "<<par[2]*bg_int<<" "<<par[2]<<endl;//2*TMath::Log(TMath::Poisson(data,model)) << endl ;
      logL -= 2*TMath::Log(TMath::Poisson(data,model));

      
      
    }
  }
  double err0 = 0.1;
  double err1 = 0.1;
  double err2 = 0.02;
  
  double scale = 1;

  logL += pow(par[0]-1*scale,2)/pow(err0,2);
  logL += pow(par[1]-1*scale,2)/pow(err1,2);
  logL += pow(par[2]-1*scale,2)/pow(err2,2);
  
  return logL ;
}

//-------------------------------------------------------------------------
//    FCN
//-------------------------------------------------------------------------
void minFunc(int& Dim, double* out, double& result, double par[], int flg) {
  result = getChi2(&par[0]); 
}




//-------------------------------------------------------------------------
//    Main
//-------------------------------------------------------------------------

void fitter_claudio() {
  gStyle->SetOptStat(0);
  //  gStyle->SetOptLogy(1);
  
  // data
  TFile *_fdata = TFile::Open("data_spectra.root");
  hdata = (TH2F*) _fdata->Get("hdata_tdrift_ne");
  //TH1D *hdx2   = (TH1D*) hdata2->ProjectionX("hdx2");
   
 

  //cout<<"data before "<<hdata->Integral()<<endl;
  
  
  // MC
  //TFile *_fin = TFile::Open("bg_histos.root");
  //hcryo = (TH2F*) _fin->Get("hcryo");
  //hpmt  = (TH2F*) _fin->Get("hpmt");
  //hint  = (TH2F*) _fin->Get("hint");


  //Paolo's file
  TFile *_fin = TFile::Open("lowEne2D.root");
  //hdata = (TH2F*) _fin->Get("hdata_tdrift_ne");
  hcryo = (TH2F*) _fin->Get("hsumTDcry");
  hpmt  = (TH2F*) _fin->Get("hsumTDpmt");
  hint  = (TH2F*) _fin->Get("hsumTDint");

  //hdata->RebinY(10);


  //TH1D *hdx1   = (TH1D*) hdata->ProjectionX("hdx1");

  //hdx2->Draw();
  //hdx1->Draw("same");
  //return;

  
  cout<<hdata->GetNbinsX()<<" "<<hdata->GetNbinsY()<<" X min "<<hdata->GetXaxis()->GetXmin()<<" Xmax "<<hdata->GetXaxis()->GetXmax()<<" Y "<<hdata->GetYaxis()->GetXmin()<<" "<<hdata->GetYaxis()->GetXmax()<<endl;
  cout<<hint->GetNbinsX()<<" "<<hint->GetNbinsY()<<" "<<hint->GetXaxis()->GetXmin()<<" "<<hint->GetXaxis()->GetXmax()<<" Y "<<hint->GetYaxis()->GetXmin()<<" "<<hint->GetYaxis()->GetXmax()<<endl;

  hdata->Draw("colz");
  hint->Draw("box same");
  //return;

   // model 
  buildModel();
  hmodel->Draw();
  //return ;

  
  int rebinX = 1;
  int rebinY = 1;
  
  //hdata->Sumw2();//
  //hcryo->Sumw2();//RebinX(rebinX); 
  //hpmt->Sumw2();//RebinX(rebinX); 
  //hint->Sumw2();//RebinX(rebinX); 
  //hmodel->Sumw2();//RebinX(rebinX);

  hdata->RebinX(rebinX); 
  hcryo->RebinX(rebinX); 
  hpmt->RebinX(rebinX); 
  hint->RebinX(rebinX); 
  hmodel->RebinX(rebinX); 
  

  hdata->RebinY(rebinY); 
  hcryo->RebinY(rebinY); 
  hpmt->RebinY(rebinY); 
  hint->RebinY(rebinY); 

  cout<<"data after "<<hdata->Integral()<<endl;


  const int NVAR = 4;
  double arglist[30];
  
  TMinuit *minLogL  = minLogL = new TMinuit(NVAR);
  int ierflg=0;  //error code for Minuit

  minLogL->SetFCN(minFunc);//Definition of the function which is minimized
  minLogL->mnexcm("SET PRINT", arglist, 1, ierflg);
  
  
  minLogL->mnparm(0,"cryo",     1,       0.5 ,     0.,     20,      ierflg);
  minLogL->mnparm(1,"pmt",      1,       0.5,      0.,     20,      ierflg);
  minLogL->mnparm(2,"int",      1,       0.5,      0.,     20,      ierflg);
  //minLogL->mnparm(0,"cryo",     1,       0.5 ,     0,     20,      ierflg);
  //minLogL->mnparm(1,"pmt",      1,       0.5,      0,     20,      ierflg);
  //minLogL->mnparm(2,"int",      1,       0.5,      0,     20,      ierflg);
  minLogL->mnparm(3,"axion",    0,       0.1,      0,     1e6,     ierflg);

  
  if(fix){
    minLogL->FixParameter(0);  
    minLogL->FixParameter(1);  
    minLogL->FixParameter(2);
  }
  //minLogL->FixParameter(3);  
 // minLogL->FixParameter(3);
  //for(int i=0;i<int(vdata.size());++i) minLogL->FixParameter(5+i*4+2);
  
  minLogL->mnexcm("MIGRAD",arglist,0,ierflg);//minimization with Migrad
  minLogL->mnexcm("MIGRAD",arglist,0,ierflg);//minimization with Migrad
  minLogL->mnexcm("MIGRAD",arglist,0,ierflg);//minimization with Migrad
  minLogL->mnexcm("MIGRAD",arglist,0,ierflg);//minimization with Migrad
  minLogL->mnexcm("MIGRAD",arglist,0,ierflg);//minimization with Migrad
  minLogL->mnexcm("MIGRAD",arglist,0,ierflg);//minimization with Migrad

  double best[NVAR], err[NVAR];
  for(int i=0;i<NVAR;++i) {
    minLogL->GetParameter(i,best[i],err[i]);
    cout << best[i] << " +- " << err[i] << endl ;
  } 
  double ga  = pow(best[3],1/4.)*1e-12;
  double ega = pow(err[3],1/4.)*1e-12;
  cout << ga << " +- " << ega <<  endl ;
  
  
  
  hcryo->Scale(best[0]);
  hpmt->Scale(best[1]);
  hint->Scale(best[2]);
  hmodel->Scale(best[3]);

  
  TH1D *hdx   = (TH1D*) hdata->ProjectionX("hdx");
  TH1D *hcx   = (TH1D*) hcryo->ProjectionX("hcx");
  TH1D *hpx   = (TH1D*) hpmt->ProjectionX("hpx");
  TH1D *hix   = (TH1D*) hint->ProjectionX("hix");
  TH1D *haxx = (TH1D*) hix->Clone("hallx");
 TH1D *hallx = (TH1D*) hix->Clone("hallx");
  hallx->Reset();
  //hdx->Reset();
  hcx->Reset();
  hpx->Reset();
  hix->Reset();
  haxx->Reset();
  hdx->Reset();

  int minX = hdata->GetXaxis()->FindBin(100.);
  int maxX = hdata->GetXaxis()->FindBin(199.);
  int minY = hdata->GetXaxis()->FindBin(0.);
  int maxY= hdata->GetXaxis()->FindBin(400.);
  TH1D *hdy = (TH1D*) hdata->ProjectionY("hdy");//,minX,maxX);
  TH1D *hcy = (TH1D*) hcryo->ProjectionY("hcy");
  TH1D *hpy = (TH1D*) hpmt->ProjectionY("hpy");
  TH1D *hiy = (TH1D*) hint->ProjectionY("hiy");
  TH1D *hally = (TH1D*) hiy->Clone("hally");
  hally->Reset();
  hcy->Reset();
  hpy->Reset();
  hiy->Reset();
  hdy->Reset();

  double dtX=0;
  double mcallX=0;
  double cryoX=0;
  double pmtX=0;
  double intX=0;
  double axX=0;
  for(int i=1;i<=hdata->GetNbinsX();i++){
    for(int k=0;k<=hdata->GetNbinsY();k++){//maxY;k++){
      dtX += hdata->GetBinContent(i,k);
      mcallX += hcryo->GetBinContent(i,k)+hpmt->GetBinContent(i,k)+hint->GetBinContent(i,k)+hmodel->GetBinContent(i,k);
      cryoX += hcryo->GetBinContent(i,k);
      pmtX += hpmt->GetBinContent(i,k);
      intX += hint->GetBinContent(i,k);
      axX += hmodel->GetBinContent(i,k);
    }
    hallx->SetBinContent(i,mcallX);
    hcx->SetBinContent(i,cryoX);
    hpx->SetBinContent(i,pmtX);
    hix->SetBinContent(i,intX);
    hdx->SetBinContent(i,dtX);
    haxx->SetBinContent(i,axX);
    mcallX=0;
    cryoX=0;
    pmtX=0;
    intX=0;
    dtX=0;
    axX=0;
  }

  double dtY=0;
  double mcallY=0;
  double cryoY=0;
  double pmtY=0;
  double intY=0;
  for(int i=1;i<=hdata->GetNbinsY();i++){
    for(int k=minX;k<=maxX;k++){//hdata->GetNbinsX();k++){
      //for(int k=1;k<=hdata->GetNbinsX();k++){
      dtY += hdata->GetBinContent(k,i);
      mcallY += hcryo->GetBinContent(k,i)+hpmt->GetBinContent(k,i)+hint->GetBinContent(k,i);
      cryoY += hcryo->GetBinContent(k,i);
      pmtY += hpmt->GetBinContent(k,i);
      intY += hint->GetBinContent(k,i);
    }
    hally->SetBinContent(i,mcallY);
    hcy->SetBinContent(i,cryoY);
    hpy->SetBinContent(i,pmtY);
    hiy->SetBinContent(i,intY);
    hdy->SetBinContent(i,dtY);
    dtY=0;
    mcallY=0;
    cryoY=0;
    pmtY=0;
    intY=0;
  }

  TH1D *hall = (TH1D*) hint->Clone("hall");
  hall->Reset();
  hall->Add(hcryo,1);
  hall->Add(hpmt,1);
  hall->Add(hint,1);

  
  //hallx->Add(hcx,1);
  //hallx->Add(hpx,1);
  //hallx->Add(hix,1);
  //hallx->Add(hmodel,1);
  TCanvas *c3 = new TCanvas("c3","all",700,500);
  hdata->Draw("colz");
  hall->Draw("box,same");

  TCanvas *c1 =new TCanvas("c1","c1");
  hdx->GetXaxis()->SetRangeUser(6,200);
  hdx->GetYaxis()->SetRangeUser(0,2000*rebinX);
  hdx->Draw("e");
  hdx->SetLineColor(1);
  hallx->SetLineColor(2);
  hcx->SetLineColor(kOrange);
  hpx->SetLineColor(kGreen);
  hix->SetLineColor(kBlue);
  haxx->SetLineColor(64);

  double dt = hdx->Integral(minX,maxX);
  double mc = hallx->Integral(minX,maxX);

  hallx->Scale(dt/mc);
  hcx->Scale(dt/mc);
  hpx->Scale(dt/mc);
  hix->Scale(dt/mc);

    
  cout<<dt<<" "<<mc<<" "<<dt/mc<<endl;
  
  hallx->Draw("same,histo");
  hcx->Draw("same,histo");
  hpx->Draw("same,histo");
  hix->Draw("same,histo");
  haxx->Draw("same,histo");
  
  TLegend *leg = new TLegend(0.12,0.7,0.2,0.88);
  leg->AddEntry(hdy,"Data","L");
  leg->AddEntry(hally,"MC tot","L");
  leg->AddEntry(hcy,"MC cryo","L");
  leg->AddEntry(hpy,"MC pmts","L");
  leg->AddEntry(hiy,"MC Kr+Ar","L");
  leg->SetLineColor(0);
  
  leg->Draw();
 


  
  TCanvas *c2 =new TCanvas("c2","c2");
  
  double nax = hmodel->Integral(minX,maxX); 
  for(int i=1;i<hally->GetXaxis()->GetNbins();++i) {
    int nbins = hally->GetXaxis()->GetNbins();
    double v = hally->GetBinContent(i);
    hally->SetBinContent(i,v+nax/nbins);
  }
  
  //hally->Add(hmodel,best[3]);
  
  hdy->Draw("e");
  //hally->SetLineColor(2);

  hdy->SetLineColor(1);
  hally->SetLineColor(2);
  hcy->SetLineColor(kOrange);
  hpy->SetLineColor(kGreen);
  hiy->SetLineColor(kBlue);


  
  
  
  hally->Draw("same,histo");
  hcy->Draw("same,histo");
  hpy->Draw("same,histo");
  hiy->Draw("same,histo");

  leg->Draw();
  
  cout<<"all "<<hdata->Integral()<<" MC "<<hall->Integral()<<endl;
  cout<<"pX "<<hdy->Integral()<<" MC "<<hally->Integral()<<endl;
  cout<<"pY "<<hdx->Integral()<<" MC "<<hallx->Integral()<<endl;
  

  cout<<getChi2(best)<<endl;
  
  return ;
}
