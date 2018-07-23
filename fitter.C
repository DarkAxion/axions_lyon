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
#include "TRandom3.h"


using namespace std;
using namespace TMath;

bool fix        = false;
bool testFake   = false;
bool testShapes = true;
bool separate   = false;
bool ShowChi2   = false;
bool newStat    = true;
bool doPlots    = true;

TH2F *hdata;
TH2F *hcryo;
TH2F *hpmt;
TH2F *hint;
TH2F *har39;
TH2F *hkr85;
TH2F *hmodel_2D;

TH1F *hint_new;
TH1F *hmodel ;


double final_chi2;

int rebinX = 1;
int rebinY = 8;

int fit_Nemin = 15;
int fit_Nemax = 200;

int Ne_tdrift = 200;

int Ne_min = 0;
int Ne_max = 200;

TGraph *gr_chi2 = new TGraph();
int idx = 0;


//-------------------------------------------------------------------------
//    Fill 2D histo from 1D
//-------------------------------------------------------------------------
void fill_2D(TH1F *hh, TH2F *h2){

  for(int ix=1; ix<hcryo->GetNbinsX()+1; ++ix){
    float n_ev = hh->GetBinContent(ix);
    float val = n_ev/hcryo->GetNbinsY();
    for(int iy=1; iy<hcryo->GetNbinsY()+1;++iy){
      h2->SetBinContent(ix, iy, val);
    }
  }

}

//-------------------------------------------------------------------------
//    Solar axion spectrum function
//-------------------------------------------------------------------------

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

void buildModel(double gAe=1e-12) {
  gCS   = new TGraph("input/cs_photoelectric.dat");
  gFlux = new TGraph("input/gae_flux.dat");

  TFile *fNe = TFile::Open("../MC/ERscale_Mar30.root");
  //TFile *fNe = TFile::Open("scale_30.root");
  TGraph *gNeE  = (TGraph*) fNe->Get("gScale");
  gENe = new TGraph ;
  int N = gNeE->GetN();
  for(int i=0;i<N;++i) {
    double x,y ;
    gNeE->GetPoint(i,x,y);
    gENe->SetPoint(i,y,x);
  }

  int nbinsX = hcryo->GetXaxis()->GetNbins();

  hmodel = new TH1F("hmodel","",
     hcryo->GetXaxis()->GetNbins(),
     hcryo->GetXaxis()->GetBinLowEdge(1),
     hcryo->GetXaxis()->GetBinLowEdge(nbinsX-1)+hcryo->GetXaxis()->GetBinWidth(3)+1);

  TH1F *htemp = (TH1F *) hmodel->Clone("htemp");

  TF1 *fun = new TF1("fun",fSolar,0,200,2);
  fun->SetParameters(gAe,0);
  //fun->Draw();


  for(int i=1;i<htemp->GetXaxis()->GetNbins()+1;++i) {
    float   ne   = htemp->GetBinCenter(i);
    double  val = fun->Eval(ne);
    htemp->SetBinContent(i,val);
  }

  TRandom3 *rnd = new TRandom3();
  bool IsCentral = false;
  for(int i=0; i<100000;++i){

    double Ne  = htemp->GetRandom();
    double E   = gNeE->Eval(Ne);

    if(rnd->Rndm() < 0.0620/0.4310) IsCentral = true;

    double Ne_bin = rnd->Binomial((int)(E/0.0195), Ne/(E/0.0195));

    double rms  = GetNeRMS(IsCentral,Ne_bin);

    double Ne_f = rnd->Gaus(Ne_bin, 0.2*sqrt(Ne_bin));

    hmodel->Fill(Ne_f);

  }

  float n1 = htemp->Integral();
  float n2 = hmodel->Integral();

  float scale = n1/n2;
  hmodel->Scale(scale);
  hmodel->Draw();

  fill_2D(hmodel, hmodel_2D);

  hmodel_2D->Draw("colz");

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
//   Correction of energy scale difference
//-------------------------------------------------------------------------
TH1F* scale_correction(TH1F *hh){

  TH1F *hnew = (TH1F *) hh->Clone("hnew");
  hnew->Reset();

  TFile *f = new TFile("scale_correction.root");
  TGraph *gr = (TGraph *) f->Get("gr_corr");

  for(int i=1; i < hh->GetNbinsX(); ++i){

    double Ne = hh->GetBinCenter(i);
    double val = hh->GetBinContent(i);
    double Ne_new = gr->Eval(Ne);
    int bin = hnew->FindBin(Ne_new);
    hnew->SetBinContent(bin, val);

  }

  hnew->SetLineColor(2);

  //new TCanvas;
  //hh->Draw("histo");
  //hnew->Draw("same, histo");

}

//-------------------------------------------------------------------------
//    2D histogram for internal bkg to test shapes
//-------------------------------------------------------------------------
void makeInt(TH1F *s_hkr85, TH1F* s_har39){

  double norm = s_har39->Integral();
  s_har39->Scale(0.806*1e-3*46.7*0.46*432.*24*3600*1*0.0206169/norm);

  norm = s_hkr85->Integral();
  s_hkr85->Scale(1.961*1e-3*46.7*0.46*432.*24*3600*0.998949*0.0224819/norm);

  //TH1F *s_har39_corr = scale_correction(s_har39);
  //return;

  fill_2D(s_har39, har39);
  fill_2D(s_hkr85, hkr85);

  for(int ix=1; ix<hcryo->GetNbinsX()+1; ++ix){

    float n_ev = s_har39->GetBinContent(ix) + s_hkr85->GetBinContent(ix);
    float val = n_ev/hcryo->GetNbinsY();
    for(int iy=1; iy<hcryo->GetNbinsY()+1;++iy){
      hint->SetBinContent(ix, iy, val);
      }
    }
  cout << hint->Integral() << endl;
  //hint->Draw("colz");
  cout << "Histograms created" << endl;
}

//-------------------------------------------------------------------------
//    Plot all shapes normalized to 1
//-------------------------------------------------------------------------
void plotShapes(TH1D *hdata, TH1D *hint, TH1D *hpmt, TH1D *hcryo, TH1D *hkr, TH1D *har){

    double norm = hdata->Integral();
    hdata->Scale(1/norm);

    norm = hint->Integral();
    hint->Scale(1/norm);

    norm = hpmt->Integral();
    hpmt->Scale(1/norm);

    norm = hcryo->Integral();
    hcryo->Scale(1/norm);

    norm = hkr->Integral();
    hkr->Scale(1/norm);

    norm = har->Integral();
    har->Scale(1/norm);

    new TCanvas;
    hdata->Draw();
    //hint->Draw("same,histo");
    hpmt->Draw("same,histo");
    hcryo->Draw("same,histo");
    har->Draw("same, histo");
    hkr->Draw("same, histo");

}

//-------------------------------------------------------------------------
//    Chi2 minimization
//-------------------------------------------------------------------------
double getChi2(double *par) {
  double logL = 0;
  for(int ix=1;ix<hdata->GetNbinsX();++ix) {

    double Ne    = hdata->GetXaxis()->GetBinCenter(ix);
    if(Ne < fit_Nemin) continue ;
    if(Ne > fit_Nemax) continue;
    //double axion = hmodel->GetBinContent(ix);
    //if(Ne > 50 && Ne < 55) continue ;

    if(Ne > Ne_tdrift) {
      for(int iy=1;iy<hdata->GetNbinsY();++iy) {
	       double td = hdata->GetYaxis()->GetBinCenter(iy);
	       if(td > 300)      continue ;
	       if(td < 50)      continue ;
	       double data     = hdata->GetBinContent(ix,iy);
	       double bg_cryo  = hcryo->GetBinContent(ix,iy);
	       double bg_pmt   = hpmt->GetBinContent(ix,iy);
	       //double width    = hdata->GetYaxis()->GetBinWidth(3);
         double bg_ar39 = har39->GetBinContent(ix,iy);
         double bg_kr85 = hkr85->GetBinContent(ix,iy);
         double bg_int   = hint->GetBinContent(ix,iy);
         //double signal   = axion/375.*width ;
         double signal = hmodel_2D->GetBinContent(ix,iy);
         double model;
	       if(separate) {
           model = (par[1]*bg_cryo + par[2]*bg_pmt + par[3]*bg_ar39 + par[4]*bg_kr85 + par[0]*signal) ;
         }
         else {
           model = (par[1]*bg_cryo + par[2]*bg_pmt + par[3]*bg_int + par[0]*signal) ;
         }

         if(Ne>65 && Ne<66)
            cout<<iy<<" "<<Ne<<" "<<par[0]<< " "<<signal<<" data "<<data<<" "<<model<<endl;

	       if(model == 0)    continue ;
	       if(data  == 0)    continue ;
         if(TMath::Poisson(data,model) == 0) continue ;
	       //logL += pow(data-model,2)/data; 	//if(hdata->GetXaxis()->GetNbins())
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
      double axions = 0;
      for(int iy=1;iy<=hdata->GetNbinsY();++iy) {
        data     += hdata->GetBinContent(ix,iy);
        bg_int   += hint->GetBinContent(ix,iy);
        bg_cryo  += hcryo->GetBinContent(ix,iy);
        bg_pmt   += hpmt->GetBinContent(ix,iy);
        bg_ar39  += har39->GetBinContent(ix,iy);
        bg_kr85  += hkr85->GetBinContent(ix,iy);
        axions   += hmodel_2D->GetBinContent(ix,iy);
        //if(Ne>55 && Ne<56)cout<<iy<<" "<<axions<< " " <<hmodel_2D->GetBinContent(ix,iy)<< endl;

      }
      //double width    = hdata->GetYaxis()->GetBinWidth(3);


      double model = 0;
      if(separate){
        model = (par[1]*bg_cryo + par[2]*bg_pmt + par[3]*bg_ar39 + par[4]*bg_kr85 + par[0]*axions);
      }
      else {
        model = (par[1]*bg_cryo + par[2]*bg_pmt + par[3]*bg_int + par[0]*axions) ;
      }
      if(Ne == 55.5)cout<<"no tdrift "<< par[1]<< " " << data << " "<<model<<endl;

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
  double err2 = 0.1;

  logL += pow(par[1]-1,2)/pow(err0,2);
  logL += pow(par[2]-1,2)/pow(err0,2);
  logL += pow(par[3]-1,2)/pow(err0,2);
  //if(separate) logL += pow(par[4]-1,2)/pow(err0,2);
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

void fitter(int fixAxion) {
  //gStyle->SetOptLogx(1);
  //gStyle->SetOptLogy(1);

  if(newStat){

    TFile *_fin = TFile::Open("../MC/bg_histos_new.root");
    hcryo = (TH2F*) _fin->Get("hcryo");
    hpmt  = (TH2F*) _fin->Get("hpmt");
    hint  = (TH2F*) _fin->Get("hint");
    if(!testShapes){
      har39 = (TH2F *) _fin->Get("h_int_ar39");
      hkr85 = (TH2F *) _fin->Get("h_int_kr85");
    }

  }else{

    TFile *_fin = TFile::Open("lowEne2D.root");
    //hdata = (TH2F*) _fin->Get("hdata_tdrift_ne");
    hcryo = (TH2F*) _fin->Get("hsumTDcry");
    hpmt  = (TH2F*) _fin->Get("hsumTDpmt");
    hint  = (TH2F*) _fin->Get("hsumTDint");
    har39 = (TH2F *) hint->Clone("har39");
    hkr85 = (TH2F *) hint->Clone("hkr85");
  }

  if(testShapes){

    har39 = (TH2F *) hint->Clone("har39");
    hkr85 = (TH2F *) hint->Clone("hkr85");

    TFile *_fin_int = TFile::Open("shapes.root");
    TH1F *s_har39 = (TH1F*) _fin_int->Get("har39");
    TH1F *s_hkr85 = (TH1F*) _fin_int->Get("hkr85");

    makeInt(s_har39, s_hkr85);
    //return;
  }



  // data
  //if(newStat){
    TFile *_fdata = TFile::Open("data_spectra.root");
    if(testFake) hdata = createFakeSample(0.6,0.7,0.9);
    else hdata = (TH2F*) _fdata->Get("hdata_tdrift_ne");
    TH1F* hdata_ne = (TH1F *) _fdata->Get("hdata_ne");
  //}

  // model
  hmodel_2D = (TH2F *) hint->Clone("hmodel_2D");
  hmodel_2D->Reset();
  buildModel();
  //return;


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
  hmodel_2D->RebinY(rebinY);

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
  minLogL->mnparm(1,"cryo",     1,       0.05,      0,     5,      ierflg);
  minLogL->mnparm(2,"pmt",      1,       0.05,      0,     5,      ierflg);
  if(separate){
      minLogL->mnparm(3,"ar39",      1,       0.05,      0,     5,      ierflg);
      minLogL->mnparm(4,"kr85",      1,       0.05,      0,     5,      ierflg);
  }
  else minLogL->mnparm(3,"int",      1,       0.05,      0,     5,      ierflg);



  if(fixAxion) minLogL->FixParameter(0);
  if(fix){
    minLogL->FixParameter(1);
    minLogL->FixParameter(2);
    minLogL->FixParameter(3);
  }


  minLogL->mnexcm("MIGRAD",arglist,0,ierflg);
  minLogL->mnexcm("MIGRAD",arglist,0,ierflg);
  minLogL->mnexcm("MIGRAD",arglist,0,ierflg);
  minLogL->mnexcm("MIGRAD",arglist,0,ierflg);
  minLogL->mnexcm("MIGRAD",arglist,0,ierflg);
  minLogL->mnexcm("MIGRAD",arglist,0,ierflg);



  double best[NVAR], err[NVAR];
  const int Nnames = NVAR;
  //string names[4] = {"axions","cryo","pmt","int"};
  string names[5] = {"axions","cryo","pmt","ar39","kr85"};

  for(int i=0;i<NVAR;++i) {
    minLogL->GetParameter(i,best[i],err[i]);
    cout << names[i] << ": " << best[i] << " +- " << err[i] << endl ;
  }
  double ga  = pow(best[0],1/4.)*1e-12;
  double ega = pow(err[0],1/4.)*1e-12;
  cout << "gAe: " << ga << " +- " << ega <<  endl ;

  if(doPlots){


    hmodel_2D->Scale(best[0]);
    hcryo->Scale(best[1]);
    hpmt->Scale(best[2]);
    if(separate){
        har39->Scale(best[3]);
        hkr85->Scale(best[4]);
    } else {
      hint->Scale(best[3]);
    }

/*
    //hmodel->Scale(best[0]);
    hmodel_2D->Scale(best[0]);
    hcryo->Scale(best[1]); hmodel_2D->Add(hcryo,1);
    hpmt->Scale(best[2]);  hmodel_2D->Add(hpmt,1);
    if(separate){
        har39->Scale(best[3]); hmodel_2D->Add(har39,1);
        hkr85->Scale(best[4]); hmodel_2D->Add(hkr85,1);
    } else {
      hint->Scale(best[3]); hmodel_2D->Add(hint,1);
    }

  hmodel_2D->Draw("colz");
  int tbin1 = hmodel_2D->GetYaxis()->FindBin(50);
  int tbin2 = hmodel_2D->GetYaxis()->FindBin(300);

  int ebin1 = hmodel_2D->GetXaxis()->FindBin(Ne_tdrift);
  int ebin2 = hmodel_2D->GetXaxis()->FindBin(199);

  TH1D *hne = (TH1D*) hmodel_2D->ProjectionX("hne",tbin1,tbin2);
  TH1D *hne2 = (TH1D*) hmodel_2D->ProjectionX("hne2",1,399);
  TH1D *hnt = (TH1D*) hmodel_2D->ProjectionY("hnt",ebin1,ebin2);


  TCanvas *ce = new TCanvas("ce","ce");
  hdata->ProjectionX("hde",tbin1,tbin2)->Draw("e");
  hne->SetLineColor(2);
  hne->Draw("same,histo");

  TCanvas *ce2 = new TCanvas("ce2","ce2");
  hdata->ProjectionX("hde2",1,1)->Draw("e");
  hne2->SetLineColor(2);
  hne2->Draw("same,histo");


  TCanvas *ct = new TCanvas("ct","ct");
  hdata->ProjectionY("hdt",ebin1,ebin2)->Draw("e");
  hnt->SetLineColor(2);
  hnt->Draw("same,histo");




  return;
*/

  TH1D *hdx   = (TH1D*) hdata->ProjectionX("hdx");
  TH1D *hcx   = (TH1D*) hcryo->ProjectionX("hcx");
  TH1D *hpx   = (TH1D*) hpmt->ProjectionX("hpx");
  TH1D *hkrx  = (TH1D*) hkr85->ProjectionX("harx");
  TH1D *harx  = (TH1D*) har39->ProjectionX("hkrx");
  TH1D *hix   = (TH1D*) hint->ProjectionX("hix");
  //TH1D *haxx  = (TH1D*) hix->Clone("haxx");
  TH1D *haxx  = (TH1D*) hmodel_2D->ProjectionX("haxx");
  TH1D *hallx = (TH1D*) hix->Clone("hallx");
  hix->Reset();
  hkrx->Reset();
  harx->Reset();
  hallx->Reset();
  hcx->Reset();
  hpx->Reset();
  haxx->Reset();
  hdx->Reset();

  int minX = hdata->GetXaxis()->FindBin(Ne_tdrift);
  int maxX = hdata->GetXaxis()->FindBin(199.);
  int minY = hdata->GetXaxis()->FindBin(0.);
  int maxY= hdata->GetXaxis()->FindBin(400.);
  TH1D *hdy   = (TH1D*) hdata->ProjectionY("hdy");//,minX,maxX);
  TH1D *hcy   = (TH1D*) hcryo->ProjectionY("hcy");
  TH1D *hpy   = (TH1D*) hpmt->ProjectionY("hpy");
  TH1D *hkry  = (TH1D*) hkr85->ProjectionY("hary");
  TH1D *hary  = (TH1D*) har39->ProjectionY("hkry");
  TH1D *hiy   = (TH1D*) hint->ProjectionY("hiy");
  TH1D *haxy  = (TH1D*) hmodel_2D->ProjectionY("haxy");
  TH1D *hally = (TH1D*) hiy->Clone("hally");
  hkry->Reset();
  hary->Reset();
  hiy->Reset();
  hally->Reset();
  hcy->Reset();
  hpy->Reset();
  hiy->Reset();
  haxy->Reset();
  hdy->Reset();

  double dtX=0;
  double mcallX=0;
  double cryoX=0;
  double pmtX=0;
  double intX=0;
  double arX = 0;
  double krX = 0;
  double axX=0;

  int bin = hdata->GetXaxis()->FindBin(55.5);

  for(int i=1;i<=hdata->GetNbinsX();i++){
    for(int k=1;k<=hdata->GetNbinsY();k++){//maxY;k++){
      dtX += hdata->GetBinContent(i,k);
      if(separate) mcallX += hcryo->GetBinContent(i,k)+hpmt->GetBinContent(i,k)+ har39->GetBinContent(i, k) + hkr85->GetBinContent(i, k) + hmodel_2D->GetBinContent(i, k);
      else mcallX += hcryo->GetBinContent(i,k)+hpmt->GetBinContent(i,k)+hint->GetBinContent(i,k) + hmodel_2D->GetBinContent(i, k);
      cryoX += hcryo->GetBinContent(i,k);
      pmtX += hpmt->GetBinContent(i,k);
      arX += har39->GetBinContent(i,k);
      krX += hkr85->GetBinContent(i,k);
      intX += hint->GetBinContent(i,k);
      axX += hmodel_2D->GetBinContent(i,k);
    }

    hallx->SetBinContent(i,mcallX);
    if(i == bin) cout << "After fit " << dtX << " " << mcallX << " " << hallx->GetBinContent(i) << " " << hallx->GetBinCenter(i) << endl;
    //mcallX += hmodel->GetBinContent(i);
    //axX = hmodel->GetBinContent(i);


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

  //hkr85->Draw("colz");
  //return;

  double dtY=0;
  double mcallY=0;
  double cryoY=0;
  double pmtY=0;
  double intY=0;
  double arY = 0;
  double krY = 0;
  double axy = 0;

  for(int i=1;i<=hdata->GetNbinsY();i++){
    for(int k=minX;k<=maxX;k++){//hdata->GetNbinsX();k++){
      //for(int k=1;k<=hdata->GetNbinsX();k++){
      dtY += hdata->GetBinContent(k,i);
      if(separate) mcallY += hcryo->GetBinContent(k,i)+hpmt->GetBinContent(k,i)+ har39->GetBinContent(k, i) + hkr85->GetBinContent(k, i) + hmodel_2D->GetBinContent(k, i);
      else mcallY += hcryo->GetBinContent(k,i)+hpmt->GetBinContent(k,i)+hint->GetBinContent(k,i) + hmodel_2D->GetBinContent(k, i);
      cryoY += hcryo->GetBinContent(k,i);
      pmtY += hpmt->GetBinContent(k,i);
      if(separate){
        arY += har39->GetBinContent(k,i);
        krY += hkr85->GetBinContent(k,i);
      }
      else intY += hint->GetBinContent(k,i);
      axy += hmodel_2D->GetBinContent(k,i);
    }
    hally->SetBinContent(i,mcallY);
    hcy->SetBinContent(i,cryoY);
    hpy->SetBinContent(i,pmtY);
    if(separate){
      hary->SetBinContent(i,arY);
      hkry->SetBinContent(i,krY);
    }
    else hiy->SetBinContent(i,intY);
    haxy->SetBinContent(i,axy);
    hdy->SetBinContent(i,dtY);
    dtY=0;
    mcallY=0;
    cryoY=0;
    pmtY=0;
    intY=0;
    arY = 0;
    krY = 0;
    axy = 0;
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
  hdx->GetXaxis()->SetRangeUser(Ne_min,Ne_max);
  hdx->GetYaxis()->SetRangeUser(2,2000*rebinX);
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

  //plotShapes(hdx, hix, hpx, hcx, hkrx, harx);
  //return;

  double dt = hdx->Integral(minX,maxX);
  double mc = hallx->Integral(minX,maxX);
/*
  hallx->Scale(dt/mc);
  hcx->Scale(dt/mc);
  hpx->Scale(dt/mc);
  if(separate){
    harx->Scale(dt/mc);
    hkrx->Scale(dt/mc);
  }
  else hix->Scale(dt/mc);
  */

  cout<<dt<<" "<<mc<<" "<<dt/mc<<endl;

  hallx->Draw("same,histo");
  hcx->Draw("same,histo");
  hpx->Draw("same,histo");
  if(separate){
    harx->Draw("same, histo");
    hkrx->Draw("same, histo");
  }
  else hix->Draw("same,histo");
  haxx->Draw("same,histo");
  //hmodel->Draw("same, histo");
  TLegend *leg = new TLegend(0.12,0.7,0.2,0.88);
  leg->AddEntry(hdy,"Data","L");
  leg->AddEntry(hally,"MC tot","L");
  leg->AddEntry(hcy,"MC cryo","L");
  leg->AddEntry(hpy,"MC pmts","L");
  leg->AddEntry(hiy,"MC Kr+Ar","L");
  leg->AddEntry(harx,"MC Ar","L");
  leg->AddEntry(hkrx,"MC Kr","L");
  leg->AddEntry(haxx,"axions","L");
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
  haxy->SetLineColor(64);


  hally->Draw("same,histo");
  hcy->Draw("same,histo");
  hpy->Draw("same,histo");
  if(separate){
    hary->Draw("same, histo");
    hkry->Draw("same, histo");
  }
  else hiy->Draw("same,histo");
  haxy->Draw("same, histo");

  //plotShapes(hdy, hiy, hpy, hcy, hkry, hary);
  //return;

  leg->Draw();

  cout<<"all "<<hdata->Integral()<<" MC "<<hall->Integral()<<endl;
  cout<<"pX "<<hdy->Integral()<<" MC "<<hally->Integral()<<endl;
  cout<<"pY "<<hdx->Integral()<<" MC "<<hallx->Integral()<<endl;

  if(ShowChi2){
    new TCanvas;
    gr_chi2->Draw("a*.");
    cout<<getChi2(best)<<endl;
  }
}
  return ;
}
