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


using namespace std;
using namespace TMath;

bool fix = true;
bool testFake = false;
bool testShapes = true;
bool separate = true;
bool ShowChi2 = false;

TH2F *hdata;
TH2F *hcryo;
TH2F *hpmt;
TH2F *hint;
TH2F *har39;
TH2F *hkr85;

TH1F *hint_new;
TH1D *hmodel ;

double final_chi2;

int rebinX = 2;
int rebinY = 8;

TGraph *gr_chi2 = new TGraph();
int idx = 0;

TGraph *gCS ;
TGraph *gFlux;
TGraph *gENe;

double fSolar(double *x, double *p) {
  double alpha   = 1./137.035 ;
  double me      = 510.9989461; // keV
  double density = 1.3954 ; // g/cm^3
  double day     = 432.*24*3600;
  double barn    = 1e-24; //cm2
  double kg      = 46.7*6.022e23/40.*1000; // # atomi / kg


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


  if(cs < 0) return 0 ;
  return cs;

}

//-------------------------------------------------------------------------
//    Build Model
//-------------------------------------------------------------------------

void buildModel(double gAe=1e-12) {
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

  int nbinsX = hcryo->GetXaxis()->GetNbins();
  hmodel = new TH1D("hmodel","",
     hcryo->GetXaxis()->GetNbins(),
     hcryo->GetXaxis()->GetBinLowEdge(1),
     hcryo->GetXaxis()->GetBinLowEdge(nbinsX-1)+hcryo->GetXaxis()->GetBinWidth(3)+1);


  TF1 *fun = new TF1("fun",fSolar,0,200,2);
  fun->SetParameters(gAe,0);
  //fun->Draw();


  for(int i=1;i<hmodel->GetXaxis()->GetNbins()+1;++i) {
    float   ne   = hmodel->GetBinCenter(i);
    double  val = fun->Eval(ne);
    hmodel->SetBinContent(i,val);
  }

    hmodel->Draw();

}

//-------------------------------------------------------------------------
//    Creation of a fake sample to check fit validity
//-------------------------------------------------------------------------

TH2F* createFakeSample(double A1, double A2, double A3){


  int nbinsX = hcryo->GetXaxis()->GetNbins();
  int nbinsY = hcryo->GetYaxis()->GetNbins();

  TH2F *hFake = new TH2F("hFake","",
      hcryo->GetNbinsX(),
      hcryo->GetXaxis()->GetBinLowEdge(1),
      hcryo->GetXaxis()->GetBinLowEdge(nbinsX-1)+hcryo->GetYaxis()->GetBinWidth(3)+1,
      hcryo->GetNbinsY(),
      hcryo->GetYaxis()->GetBinLowEdge(1),
      hcryo->GetYaxis()->GetBinLowEdge(nbinsY-1)+hcryo->GetYaxis()->GetBinWidth(3)+1);

cout << "Creation of fake sample" << endl;
    for(int ix=1; ix<hcryo->GetNbinsX()+1; ++ix){

      for(int iy=1; iy<hcryo->GetNbinsY()+1;++iy){
        double bg_int   = hint->GetBinContent(ix,iy);
        double bg_cryo  = hcryo->GetBinContent(ix,iy);
        double bg_pmt   = hpmt->GetBinContent(ix,iy);
        double sum = A1 * bg_int + A2 * bg_cryo + A3 * bg_pmt;

        hFake->SetBinContent(ix, iy, sum);
      }

    }
    cout << "Fake sample created" << endl;

    return hFake;

}

//-------------------------------------------------------------------------
//    Fill 2D histo from 1D
//-------------------------------------------------------------------------
void fill_2D(TH1F *hh, TH2F *h2){

  for(int ix=1; ix<hcryo->GetNbinsX()+1; ++ix){
    float val = hh->GetBinContent(ix);
    for(int iy=1; iy<hcryo->GetNbinsY()+1;++iy){
        if(iy == 1) h2->SetBinContent(ix, iy, val);
        else h2->SetBinContent(ix, iy, 0);
    }
  }

}


//-------------------------------------------------------------------------
//    2D histogram for internal bkg to test shapes
//-------------------------------------------------------------------------
void makeInt(TH1F *s_hkr85, TH1F* s_har39){

  double norm = s_har39->Integral();
  s_har39->Scale(0.806*1e-3*46.7*432.*24*3600*1*0.0206169/norm);

  norm = s_hkr85->Integral();
  s_hkr85->Scale(1.961*1e-3*46.7*432.*24*3600*0.998949*0.0224819/norm);

  fill_2D(s_har39, har39);
  fill_2D(s_hkr85, hkr85);

  for(int ix=1; ix<hcryo->GetNbinsX()+1; ++ix){

    float val = s_har39->GetBinContent(ix) + s_hkr85->GetBinContent(ix);
    for(int iy=1; iy<hcryo->GetNbinsY()+1;++iy){
      if(iy == 1) hint->SetBinContent(ix, iy, val);
       else hint->SetBinContent(ix, iy, 0);
      }
    }
  cout << hint->Integral() << endl;
  //hint->Draw("colz");
  cout << "Histograms created" << endl;
}

//-------------------------------------------------------------------------
//    Chi2 minimization
//-------------------------------------------------------------------------
double getChi2(double *par) {
  double logL = 0;
  for(int ix=1;ix<hdata->GetXaxis()->GetNbins();++ix) {

    double Ne    = hdata->GetXaxis()->GetBinCenter(ix);
    if(Ne < 10)     continue ;

    double axion = hmodel->GetBinContent(ix);

    if(Ne > 200) {
      for(int iy=1;iy<hdata->GetYaxis()->GetNbins();++iy) {
	double td = hdata->GetYaxis()->GetBinCenter(iy);
	if(td > 250)      continue ;
	if(td < 30)       continue ;
	double data     = hdata->GetBinContent(ix,iy);
	double bg_cryo  = hcryo->GetBinContent(ix,iy);
	double bg_pmt   = hpmt->GetBinContent(ix,iy);
	double width    = hdata->GetYaxis()->GetBinWidth(3);
  double bg_ar39 = har39->GetBinContent(ix,iy);
  double bg_kr85 = hkr85->GetBinContent(ix,iy);
  double bg_int   = hint->GetBinContent(ix,iy);
  double signal   = axion/375.*width ;
  double model;
	if(separate) {
    model = (par[1]*bg_cryo + par[2]*bg_pmt + par[3]*bg_ar39 + par[4]*bg_kr85 + par[0]*signal) ;
  }
  else {
    model = (par[1]*bg_cryo + par[2]*bg_pmt + par[3]*bg_int + par[0]*signal) ;
  }
	if(model == 0)    continue ;
	if(data  == 0)    continue ;
    if(TMath::Poisson(data,model) == 0) continue ;
	//logL += pow(data-model,2)/data 	//if(hdata->GetXaxis()->GetNbins())
    logL -= 2*TMath::Log(TMath::Poisson(data,model));
    //cout << "1   " <<  data << " " << model << " " << logL << " " << " " << TMath::Poisson(data,model) << " " << TMath::Log(TMath::Poisson(data,model)) <<endl ;
      }
    } else {
      double data     = 0 ;
      double bg_int   = 0;
      double bg_cryo  = 0;
      double bg_pmt   = 0;
      double bg_ar39  = 0;
      double bg_kr85  = 0;
      for(int iy=1;iy<hdata->GetYaxis()->GetNbins();++iy) {
        data     += hdata->GetBinContent(ix,iy);
        bg_int   += hint->GetBinContent(ix,iy);
        bg_cryo  += hcryo->GetBinContent(ix,iy);
        bg_pmt   += hpmt->GetBinContent(ix,iy);
        bg_ar39  += har39->GetBinContent(ix,iy);
        bg_kr85  += hkr85->GetBinContent(ix,iy);

      }
      //double width    = hdata->GetYaxis()->GetBinWidth(3);
      double model = 0;
      if(separate){
        model = (par[1]*bg_cryo + par[2]*bg_pmt + par[3]*bg_ar39 + par[4]*bg_kr85 + par[0]*axion);
      }
      else {
        model = (par[1]*bg_cryo + par[2]*bg_pmt + par[3]*bg_int + par[0]*axion) ;
      }

      if(model <= 0)     continue ;
      //if(data  <= 0)     continue ;
      if(TMath::Poisson(data,model) == 0) continue ;
      //logL += pow(data-model,2)/data ;
      logL -= 2*TMath::Log(TMath::Poisson(data,model));
      //cout << "2   " <<  data << " " << model << " " << logL << " " <<  TMath::Poisson(data,model) << " " << TMath::Log(TMath::Poisson(data,model)) << endl;
    }
  }
  double err0 = 0.1;
  double err1 = 0.1;
  double err2 = 0.01;

  //logL += pow(par[1]-1,2)/pow(err0,2);
  //logL += pow(par[2]-1,2)/pow(err2,2);
  logL += pow(par[3]-1,2)/pow(err2,2);
  //logL += pow(par[1]-1,2)/pow(err2,2);
  //cout << logL << endl;
  gr_chi2->SetPoint(idx, par[0], logL);
  //cout << idx << " " << par[3] << endl;
  idx++;
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

void fitter() {
  //gStyle->SetOptLogx(1);
  //gStyle->SetOptLogy(1);

 TFile *_fin = TFile::Open("lowEne2D.root");
 //hdata = (TH2F*) _fin->Get("hdata_tdrift_ne");
 hcryo = (TH2F*) _fin->Get("hsumTDcry");
 hpmt  = (TH2F*) _fin->Get("hsumTDpmt");
 hint  = (TH2F*) _fin->Get("hsumTDint");

  har39 = (TH2F *) hint->Clone("har39");
  hkr85 = (TH2F *) hint->Clone("hkr85");

  if(testShapes){
    TFile *_fin_int = TFile::Open("shapes.root");
    TH1F *s_har39 = (TH1F*) _fin_int->Get("har39");
    TH1F *s_hkr85 = (TH1F*) _fin_int->Get("hkr85");

    makeInt(s_har39, s_hkr85);
  }


  // data
  TFile *_fdata = TFile::Open("data_spectra.root");
  if(testFake) hdata = createFakeSample(0.6,0.7,0.9);
  else hdata = (TH2F*) _fdata->Get("hdata_tdrift_ne");
  TH1F* hdata_ne = (TH1F *) _fdata->Get("hdata_ne");


  // model
  buildModel();

  // MC
  //TFile *_fin = TFile::Open("../MC/bg_histos.root");

  hdata->RebinX(rebinX);
  hcryo->RebinX(rebinX);
  hpmt->RebinX(rebinX);
  hint->RebinX(rebinX);
  har39->RebinX(rebinX);
  hkr85->RebinX(rebinX);
  //hint_new->RebinX(rebinX);
  hmodel->RebinX(rebinX);

  hdata->RebinY(rebinY);
  hcryo->RebinY(rebinY);
  hpmt->RebinY(rebinY);
  hint->RebinY(rebinY);
  har39->RebinY(rebinY);
  hkr85->RebinY(rebinY);
  cout <<hdata->GetNbinsX()<< " " << hint->GetNbinsX() << endl;

  int N;

  if(separate) N = 5;
  else N = 4;

  const int NVAR = N;

  double arglist[30];

  TMinuit *minLogL  = minLogL = new TMinuit(NVAR);
  int ierflg=0;  //error code for Minuit

  minLogL->SetFCN(minFunc);//Definition of the function which is minimized
  minLogL->mnexcm("SET PRINT", arglist, 1, ierflg);

  minLogL->mnparm(0,"axion",    0,       0.1,      0,     1e6,     ierflg);
  minLogL->mnparm(1,"cryo",     1,       0.5,      0,     20,      ierflg);
  minLogL->mnparm(2,"pmt",      1,       0.5,      0,     20,      ierflg);
  if(separate){
      minLogL->mnparm(3,"ar39",      1,       0.5,      0,     20,      ierflg);
      minLogL->mnparm(4,"kr85",      1,       0.5,      0,     20,      ierflg);
  }
  else minLogL->mnparm(3,"int",      1,       0.5,      0,     20,      ierflg);




  minLogL->FixParameter(0);
  //minLogL->FixParameter(2);
  //minLogL->FixParameter(1);
  //minLogL->FixParameter(3);


  minLogL->mnexcm("MIGRAD",arglist,0,ierflg);
  minLogL->mnexcm("MIGRAD",arglist,0,ierflg);
  minLogL->mnexcm("MIGRAD",arglist,0,ierflg);
  minLogL->mnexcm("MIGRAD",arglist,0,ierflg);
  minLogL->mnexcm("MIGRAD",arglist,0,ierflg);
  minLogL->mnexcm("MIGRAD",arglist,0,ierflg);



  double best[NVAR], err[NVAR];
  for(int i=0;i<NVAR;++i) {
    minLogL->GetParameter(i,best[i],err[i]);
    cout << best[i] << " +- " << err[i] << endl ;
  }
  double ga  = pow(best[0],1/4.)*1e-12;
  double ega = pow(err[0],1/4.)*1e-12;
  cout << ga << " +- " << ega <<  endl ;


  hmodel->Scale(best[0]);
  hcryo->Scale(best[1]);
  hpmt->Scale(best[2]);
  if(separate){
      har39->Scale(best[3]);
      hkr85->Scale(best[4]);
  }
  else hint->Scale(best[3]);

  TH1D *hdx   = (TH1D*) hdata->ProjectionX("hdx");
  TH1D *hcx   = (TH1D*) hcryo->ProjectionX("hcx");
  TH1D *hpx   = (TH1D*) hpmt->ProjectionX("hpx");
  TH1D *hkrx  = (TH1D*) hkr85->ProjectionX("harx");
  TH1D *harx  = (TH1D*) har39->ProjectionX("hkrx");
  TH1D *hix   = (TH1D*) hint->ProjectionX("hix");
  TH1D *haxx  = (TH1D*) hix->Clone("haxx");
  TH1D *hallx = (TH1D*) hix->Clone("hallx");
  hix->Reset();
  hkrx->Reset();
  harx->Reset();
  hallx->Reset();
  hcx->Reset();
  hpx->Reset();
  haxx->Reset();
  hdx->Reset();

  int minX = hdata->GetXaxis()->FindBin(100.);
  int maxX = hdata->GetXaxis()->FindBin(199.);
  int minY = hdata->GetXaxis()->FindBin(0.);
  int maxY= hdata->GetXaxis()->FindBin(400.);
  TH1D *hdy   = (TH1D*) hdata->ProjectionY("hdy");//,minX,maxX);
  TH1D *hcy   = (TH1D*) hcryo->ProjectionY("hcy");
  TH1D *hpy   = (TH1D*) hpmt->ProjectionY("hpy");
  TH1D *hkry  = (TH1D*) hkr85->ProjectionY("hary");
  TH1D *hary  = (TH1D*) har39->ProjectionY("hkry");
  TH1D *hiy   = (TH1D*) hint->ProjectionY("hiy");
  TH1D *hally = (TH1D*) hiy->Clone("hally");
  hkry->Reset();
  hary->Reset();
  hiy->Reset();
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
  double arX = 0;
  double krX = 0;
  double axX=0;

  for(int i=1;i<=hdata->GetNbinsX();i++){
    for(int k=0;k<=hdata->GetNbinsY();k++){//maxY;k++){
      dtX += hdata->GetBinContent(i,k);
      if(separate) mcallX += hcryo->GetBinContent(i,k)+hpmt->GetBinContent(i,k)+ har39->GetBinContent(i, k) + hkr85->GetBinContent(i, k);
      else mcallX += hcryo->GetBinContent(i,k)+hpmt->GetBinContent(i,k)+hint->GetBinContent(i,k);
      cryoX += hcryo->GetBinContent(i,k);
      pmtX += hpmt->GetBinContent(i,k);
      arX += har39->GetBinContent(i,k);
      krX += hkr85->GetBinContent(i,k);
      intX += hint->GetBinContent(i,k);
      //axX += hmodel->GetBinContent(i,k);
    }
    mcallX += hmodel->GetBinContent(i);
    axX = hmodel->GetBinContent(i);
    hallx->SetBinContent(i,mcallX);
    hcx->SetBinContent(i,cryoX);
    hpx->SetBinContent(i,pmtX);
    harx->SetBinContent(i,arX);
    hkrx->SetBinContent(i,krX);
    hix->SetBinContent(i,intX);
    hdx->SetBinContent(i,dtX);
    haxx->SetBinContent(i,axX);
    mcallX=0;
    cryoX=0;
    pmtX=0;
    intX=0;
    krX = 0;
    arX = 0;
    dtX=0;
    axX=0;
  }

  double dtY=0;
  double mcallY=0;
  double cryoY=0;
  double pmtY=0;
  double intY=0;
  double arY = 0;
  double krY = 0;

  for(int i=1;i<=hdata->GetNbinsY();i++){
    for(int k=minX;k<=maxX;k++){//hdata->GetNbinsX();k++){
      //for(int k=1;k<=hdata->GetNbinsX();k++){
      dtY += hdata->GetBinContent(k,i);
      if(separate) mcallY += hcryo->GetBinContent(k,i)+hpmt->GetBinContent(k,i)+ har39->GetBinContent(k, i) + hkr85->GetBinContent(k, i);
      else mcallY += hcryo->GetBinContent(k,i)+hpmt->GetBinContent(k,i)+hint->GetBinContent(k,i);
      cryoY += hcryo->GetBinContent(k,i);
      pmtY += hpmt->GetBinContent(k,i);
      if(separate){
        arY += har39->GetBinContent(k,i);
        krY += hkr85->GetBinContent(k,i);
      }
      else intY += hint->GetBinContent(k,i);
    }
    hally->SetBinContent(i,mcallY);
    hcy->SetBinContent(i,cryoY);
    hpy->SetBinContent(i,pmtY);
    if(separate){
      hary->SetBinContent(i,arY);
      hkry->SetBinContent(i,krY);
    }
    else hiy->SetBinContent(i,intY);
    hdy->SetBinContent(i,dtY);
    dtY=0;
    mcallY=0;
    cryoY=0;
    pmtY=0;
    intY=0;
    arY = 0;
    krY = 0;
  }

  TH1D *hall = (TH1D*) hint->Clone("hall");
  hall->Reset();
  hall->Add(hcryo,1);
  hall->Add(hpmt,1);
  hall->Add(hint,1);

  gStyle->SetOptStat(0);
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
  harx->SetLineColor(kViolet);
  hkrx->SetLineColor(49);
  haxx->SetLineColor(64);
  hmodel->SetLineColor(64);

  double dt = hdx->Integral(minX,maxX);
  double mc = hallx->Integral(minX,maxX);

  hallx->Scale(dt/mc);
  hcx->Scale(dt/mc);
  hpx->Scale(dt/mc);
  if(separate){
    harx->Scale(dt/mc);
    hkrx->Scale(dt/mc);
  }
  else hix->Scale(dt/mc);
  cout<<dt<<" "<<mc<<" "<<dt/mc<<endl;

  hallx->Draw("same,histo");
  hcx->Draw("same,histo");
  hpx->Draw("same,histo");
  if(separate){
    harx->Draw("same, histo");
    hkrx->Draw("same, histo");
  }
  else hix->Draw("same,histo");
  //haxx->Draw("same,histo");
  hmodel->Draw("same, histo");
  TLegend *leg = new TLegend(0.12,0.7,0.2,0.88);
  leg->AddEntry(hdy,"Data","L");
  leg->AddEntry(hally,"MC tot","L");
  leg->AddEntry(hcy,"MC cryo","L");
  leg->AddEntry(hpy,"MC pmts","L");
  leg->AddEntry(hiy,"MC Kr+Ar","L");
  leg->AddEntry(harx,"MC Ar","L");
  leg->AddEntry(hkrx,"MC Kr","L");
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
  if(separate){
    hary->Draw("same, histo");
    hkry->Draw("same, histo");
  }
  else hiy->Draw("same,histo");

  leg->Draw();

  cout<<"all "<<hdata->Integral()<<" MC "<<hall->Integral()<<endl;
  cout<<"pX "<<hdy->Integral()<<" MC "<<hally->Integral()<<endl;
  cout<<"pY "<<hdx->Integral()<<" MC "<<hallx->Integral()<<endl;

  if(ShowChi2){
    new TCanvas;
    gr_chi2->Draw("a*.");
    cout<<getChi2(best)<<endl;
  }

  return ;
}
