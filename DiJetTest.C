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

#include "GlobalVars.h"
#include "DiJetSelector.h"


void DiJetTest()
{
	etaBins_.push_back(std::make_pair(0.0,0.1) );
	etaBins_.push_back(std::make_pair(0.1,0.2) );
	etaBins_.push_back(std::make_pair(0.2,0.3) );
	etaBins_.push_back(std::make_pair(0.3,0.4) );
	etaBins_.push_back(std::make_pair(0.4,0.5) );
	etaBins_.push_back(std::make_pair(0.5,0.6) );
	etaBins_.push_back(std::make_pair(0.6,0.7) );
	etaBins_.push_back(std::make_pair(0.7,0.8) );
	etaBins_.push_back(std::make_pair(0.8,0.9) );
	etaBins_.push_back(std::make_pair(0.9,1.0) );
	etaBins_.push_back(std::make_pair(1.0,1.1) );
	etaBins_.push_back(std::make_pair(1.1,1.2) );
	etaBins_.push_back(std::make_pair(1.2,1.3) );
	etaBins_.push_back(std::make_pair(1.3,1.4) );
	etaBins_.push_back(std::make_pair(1.4,1.5) );
	etaBins_.push_back(std::make_pair(1.5,1.6) );
	etaBins_.push_back(std::make_pair(1.6,1.7) );
	etaBins_.push_back(std::make_pair(1.7,1.8) );
	etaBins_.push_back(std::make_pair(1.8,1.9) );
	etaBins_.push_back(std::make_pair(1.9,2.0) );
	etaBins_.push_back(std::make_pair(2.0,2.1) );
	etaBins_.push_back(std::make_pair(2.1,2.2) );
	etaBins_.push_back(std::make_pair(2.2,2.3) );
	etaBins_.push_back(std::make_pair(2.3,2.4) );
	etaBins_.push_back(std::make_pair(2.4,2.5) );
	etaBins_.push_back(std::make_pair(2.5,2.8) );
	etaBins_.push_back(std::make_pair(2.8,3.0) );
	etaBins_.push_back(std::make_pair(3.0,3.2) );
	etaBins_.push_back(std::make_pair(3.2,5.0) );
	lowPuBinEdeges_.push_back(0); lowPuBinEdeges_.push_back(10); lowPuBinEdeges_.push_back(20); lowPuBinEdeges_.push_back(30); lowPuBinEdeges_.push_back(50);
	
	
	outPutFileName_ = "DiJetTest.root";
	outPutFile_ = new TFile(outPutFileName_,"RECREATE");
	// starting recomputing eff and expectation
	std::cout<<"Started"<<std::endl;
	TChain analyseDataChain2("DiJetTree");
//	analyseDataChain2.Add("/nfs/dust/cms/user/adraeger/kalibri/nTuples/MC/Z2star_pythia_v3/Summer12_DR53X_QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6_v3_ak5PFCHS.root");
	analyseDataChain2.Add("/nfs/dust/cms/user/adraeger/kalibri/nTuples/MC/Z2star_pythia_v3/Summer12_DR53X_QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6_v3_ak5FastPF.root");
	analyseDataChain2.Process("DiJetSelector.C+");
}