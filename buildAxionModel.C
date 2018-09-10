#include <fstream>
#include <iostream>
#include <string>

#include "TFile.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMath.h"
#include "math.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TLegend.h"

#include "RooWorkspace.h"
#include "RooNLLVar.h"
#include "RooStats/ModelConfig.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"


using namespace std;
using namespace RooFit;
using namespace RooStats;


TString InputData   = "./data_spectra.root";
TString InputBkg    = "../MC/bg_histos_newScale_Au14_err.root";
//TString InputBkg    = "../MC/bg_histos_new.root";
TString InputSignal = "./axion_spectra.root";

bool plotData   = true;
bool plotSignal = true;
bool plotBkg    = true;

void importHistos(RooWorkspace w){


}

void buildAxionModel(){

  //Create the workspace
  RooWorkspace w("w");

  // Access the file.
  TFile* fileData = new TFile(InputData);
  // Load the histogram.
  TH1* hdata = (TH1*) fileData->Get("hdata_ne");
  // Declare an observable x.
  RooRealVar ne("ne", "ne", 0, 200);
  // Create a binned dataset that imports the contents of TH1 and associates its contents to observable 'x'.
  RooDataHist data("obsData","obsData",RooArgList(ne), hdata);
  w.import(data);
  cout <<"Data imported" << endl;
  // Plot the imported dataset.
  if(plotData){
    new TCanvas;
    RooPlot* frameData = ne.frame();
    data.plotOn(frameData);
    frameData->Draw();
  }


  //Import signal histogram
  TFile* fileSignal = new TFile(InputSignal);
  // Load the histogram.
  TH1* haxion = (TH1*) fileSignal->Get("haxion_0");
  //haxion->Scale(1./haxion->Integral());
  cout << "Integral " << haxion->Integral() << endl;
  // Create a binned dataset that imports the contents of TH1 and associates its contents to observable 'x'.
  RooDataHist hsignal("hsignal","hsignal", ne, haxion);



  //Import background Histograms
  //Import signal histogram
  TFile* fileBkg = new TFile(InputBkg);
  // Load the histogram.
  TH1* har = (TH1*) fileBkg->Get("hNe_ar39_0");
  //har->Scale(1./har->Integral());
  TH1* hkr = (TH1*) fileBkg->Get("hNe_kr85_0");
  //hkr->Scale(1./hkr->Integral());
  TH1* hpmt  = (TH1*) fileBkg->Get("hNe_pmt_0");
  //hpmt->Scale(1./hpmt->Integral());
  TH1* hcryo = (TH1*) fileBkg->Get("hNe_cryo_0");
  //hcryo->Scale(1./hcryo->Integral());

  // Create a binned dataset that imports the contents of TH1 and associates its contents to observable 'x'.
  RooDataHist har39("har39","har39",ne, har);
  RooDataHist hkr85("hkr85","hkr85",ne, hkr);
  RooDataHist hPMT("hpmt","hpmt",ne, hpmt);
  RooDataHist hCryo("hcryo","hcryo",ne, hcryo);

  cout << "Number of evts: " <<  har39.sum(kFALSE) << " " << hkr85.sum(kFALSE) << " " << hPMT.sum(kFALSE) << " " << hCryo.sum(kFALSE) << endl;

  //Transform histograms into pdfs
  RooHistPdf ar39("ar39", "ar39", ne, har39);
  RooHistPdf kr85("kr85", "kr85", ne, hkr85);
  RooHistPdf PMT("PMT",  "PMT",  ne, hPMT);
  RooHistPdf cryo("cryo", "cryo", ne, hCryo);
  RooHistPdf sig("signal", "signal", ne, hsignal);

  if(plotSignal){
    new TCanvas;
    RooPlot* frameSignal = ne.frame();
    sig.plotOn(frameSignal);
    frameSignal->Draw();
  }


  if(plotBkg){
    new TCanvas;
    RooPlot* frameBkg = ne.frame();
    ar39.plotOn(frameBkg);
    kr85.plotOn(frameBkg);
    PMT.plotOn(frameBkg);
    cryo.plotOn(frameBkg);
    frameBkg->Draw();
  }

  //Create normalization constants
  RooRealVar nsig("nsig","nnsig", 0., 0., 1e6);
  RooRealVar nbkg_ar("nbkg_ar","nbkg_ar", har39.sumEntries("ne >= 0 && ne <=200."), 0., 1e5);
  RooRealVar nbkg_kr("nbkg_kr","nbkg_kr", hkr85.sumEntries("ne >= 0 && ne <=200."), 0., 1e5);
  RooRealVar nbkg_pmt("nbkg_pmt","nbkg_pmt", hPMT.sumEntries("ne >= 0 && ne <=200."), 0., 1e5);
  RooRealVar nbkg_cryo("nbkg_cryo","nbkg_cryo", hCryo.sumEntries("ne >= 0 && ne <=200."), 0., 1e5);

  nsig.setConstant();

  cout << "Model creation" << endl;

  //w.factory("SUM::bmodel(nbkg_ar*ar39, nbkg_kr*kr85, nbkg_pmt*PMT, nbkg_cryo*cryo)");
  RooAddPdf *model = new RooAddPdf("model","Total model",RooArgList(sig, ar39,kr85,PMT,cryo),RooArgList(nsig,nbkg_ar,nbkg_kr,nbkg_pmt, nbkg_cryo));
  //w.factory("SUM::bmodel(RooArgList(ar39,kr85,PMT,cryo),RooArgList(nbkg_ar,nbkg_kr,nbkg_pmt, nbkg_cryo))");

  //model->fitTo(data,Extended(kTRUE), SumW2Error(kTRUE), Save(kTRUE));

  new TCanvas;
  RooPlot * plot = ne.frame();
  data.plotOn(plot);
  model->plotOn(plot);
  //draw the two separate pdf's
  model->plotOn(plot, RooFit::Components("ar39"), RooFit::LineColor(kViolet), RooFit::LineStyle(kDashed) );
  model->plotOn(plot, RooFit::Components("kr85"), RooFit::LineColor(49), RooFit::LineStyle(kDashed) );
  model->plotOn(plot, RooFit::Components("PMT"), RooFit::LineColor(kGreen), RooFit::LineStyle(kDashed) );
  model->plotOn(plot, RooFit::Components("cryo"), RooFit::LineColor(kOrange), RooFit::LineStyle(kDashed) );
  model->plotOn(plot, RooFit::Components("sig"), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed) );
  plot->Draw();   // to show the RooPlot in the current ROOT Canvas

  cout << "Value " << nbkg_ar.getValV() << endl;

  RooStats::ModelConfig mc("ModelConfig",&w);
  mc.SetPdf(*model);
  mc.SetParametersOfInterest(nsig);
  mc.SetObservables(ne);

  // define set of nuisance parameters
  w.defineSet("nuisParams","nbkg_ar,nbkg_kr,nbkg_pmt,nbkg_cryo");
  mc.SetNuisanceParameters(*w.set("nuisParams"));


  w.import(mc);


  w.Print();



  w.writeToFile("AxionModel.root");

}
