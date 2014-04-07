#include <TH2F.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <iostream>
#include "TStyle.h"
#include <ostream>
#include <fstream> 
#include <string.h>
#include <TTree.h>
#include <TH1.h>
#include <TH3.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <string.h>
#include <TMath.h>
#include <TF1.h>
#include <stdio.h>
#include "TROOT.h"
#include "TMinuit.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "Math/FitMethodFunction.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TList.h"
#include "TROOT.h"
//#include "functions.h"

#ifndef _mal2_hpp_
#define _mal2_hpp_
// variables
bool debug_;
char buffer_[100];
int n_;
TString TTemp_("");
vector<double> lowPuBinEdeges_, truthPuMeans_;
vector <vector<TF1*> > TF1DPUInclusive_;
vector<TString> inputFolderNames_, inputTH2Names_; 
double quantilesLow_, quantilesHigh_;
double lowPTConstFit_, highPTConstFit_;
vector<TH1D*> cbValues_;
vector<double> cbFixValues_;
vector<double> etaDeptendentFixValues_;
// input output file folders
TDirectory *InputD_, *OutPutD_;
TFile *outF_, *inF_;
// methods
bool exampleFunctions(TDirectory *OutPutF);
vector<TH1D*> th1creator ( TH2D *inputTH2D, TDirectory *outPutFolder);
TF1* cbFitting(TH1D *th1D, std::vector<double> cbFixValues, bool useFixValues);
void writeTH1D(TH1D* inputTH1D, TDirectory *outPutFolder, std::vector<TString> FunctionNames);
// functions
Double_t cb(Double_t *x, Double_t *par);
Double_t gauss (Double_t *x, Double_t *par);

#endif
