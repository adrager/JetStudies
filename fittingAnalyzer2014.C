#include <TH2F.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <iostream>
#include "TStyle.h"
#include <ostream>
#include <fstream> 
#include <iomanip>
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
#include "TVirtualFitter.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
//#include "functions.h"

#ifndef _mal2_hpp_
#define _mal2_hpp_
// variables
double ptRangeLow_, ptRangeHigh_;
bool debug_;
vector<std::pair <double,double> > lowPuBinEdeges_;
vector<std::pair <double,double> > etaBins_;
std::vector<TString> ptBinTH1Names_;
char buffer_[100];
int n_;
TString TTemp_("");
vector<TH1D*> cbValues_;
vector<double> meanPUs_;
TH1D *npuVSRho_;
vector<TH1D*> cbForFourDMinimizer_;
vector<TString> inputFolderNames_, inputTH2Names_;
ofstream *FittingAnalyzer2014ErrorReport, *FitResults;
// input output file folders
TDirectory *InputD_, *OutPutD_;
TFile *outF_, *inF_;
// methods
//bool exampleFunctions(TDirectory *OutPutF);
//vector<TH1D*> th1creator ( TH2D *inputTH2D, TDirectory *outPutFolder);
//TF1* cbFitting(TH1D *th1D);
// functions
Double_t exp0 (Double_t *x, Double_t *par);
Double_t fConst (Double_t *x, Double_t *par);
Double_t fLinear (Double_t *x, Double_t *par);
Double_t fWurzelN(Double_t *x, Double_t *par);
Double_t fSigmaOld (Double_t *x, Double_t *par);
Double_t fSigmaNoPU (Double_t *x, Double_t *par);
Double_t fSigmaTwoTerms (Double_t *x, Double_t *par);
Double_t fSigma3Terms (Double_t *x, Double_t *par);
Double_t fSigma3TermsMinus (Double_t *x, Double_t *par);
Double_t fSigma4Terms (Double_t *x, Double_t *par);
Double_t fSigma4TermsMinus (Double_t *x, Double_t *par);
double fourDminimizer(const double *par);
Double_t TwoDSigmaFitFunction(Double_t *x,Double_t *par);
TGraphErrors* TGraphCreator(TH1D* inputTH1D, unsigned int etabin, unsigned int puBin, TFile* inputFile);
void combinedPlots(vector<TH1D*> th1s, TDirectory *outPutFolderSpecial, TString name, bool rename);
void differencePlots (vector<TH1D*> th1s, TDirectory *outPutFolderSpecial, TString name, vector<double> meanPUs_);
void PUInclusivePlots (unsigned int parameter, TString parName);

//Double_t cb(Double_t *x, Double_t *par);
//Double_t gauss (Double_t *x, Double_t *par);

#endif



void fittingAnalyzer2014()
{
	gROOT->SetBatch(true);
	ptRangeLow_=20;
	ptRangeHigh_=1200;
	cout<<"Analyzer started..."<<endl;

	inF_ = TFile::Open("Fitting2014.root","UPDATE");
	TFile * ptMeanFile = TFile::Open("PTBinsMeans.root","UPDATE");
	npuVSRho_ = (TH1D*) ptMeanFile->Get("npuVsRho");
	int nn = sprintf(buffer_,"Results2014_ptRange_%.0f_%.0f.root",ptRangeLow_,ptRangeHigh_); // mean1
	FittingAnalyzer2014ErrorReport = new ofstream("FittingAnalyzer2014ErrorReport.txt");
	if(FittingAnalyzer2014ErrorReport->is_open() ) std::cout<<"FittingAnalyzer2014ErrorReport"<<std::endl;
	*FittingAnalyzer2014ErrorReport << "----------------- fittingAnalyzer2014 started--------------\n";
	*FittingAnalyzer2014ErrorReport << "This file containes error messages and warnings.\n";
	*FittingAnalyzer2014ErrorReport << "In particular warnings are saved if a free fit did not converge these entries have an error of 0 they will be reported.\n";
	
	FitResults = new ofstream("FitResults.txt");
	if(FitResults->is_open() ) std::cout<<"FitResults opend"<<std::endl;
	*FitResults << "[resolution]\n";
	*FitResults << "\n";
	*FitResults << "[1 |JetEta| 1 JetPt 1 NPU  CBGaussSigma Resolution sigma]\n";
	*FitResults << "[sigma]\n";
	*FitResults << "{1  |JetEta|  1 JetPt    NPU    (sqrt([0] *[0] + [1]*[1]/x[0] + TMath::Sign(1.,[2]) * ([2])*([2])/(x[0]*x[0]) + ([3]*sqrt(x[1])) * ([3]*sqrt(x[1]))/(x[0]*x[0])) ) PAR0 \\sigma } \n";
	
	TTemp_ = buffer_;
	outF_ = new TFile(TTemp_,"RECREATE");
	outF_->mkdir(TTemp_);
	OutPutD_ = (TDirectory*)outF_->Get(TTemp_);
	TH1D *TruthPUMeans = (TH1D*) inF_->Get("TruthPUMeans")->Clone();
	for(int i=0; i < 4; i++)
	{
		cout<<"Truth mean Pileups"<<TruthPUMeans->GetBinContent(i+1)<<endl;
		meanPUs_.push_back(TruthPUMeans->GetBinContent(i+1));
	}
	inputFolderNames_.push_back("MCTruthResolPUEta00"); 	etaBins_.push_back(std::make_pair(0.0,0.1) );
	inputFolderNames_.push_back("MCTruthResolPUEta01"); 	etaBins_.push_back(std::make_pair(0.1,0.2) );
	inputFolderNames_.push_back("MCTruthResolPUEta02"); 	etaBins_.push_back(std::make_pair(0.2,0.3) );
	inputFolderNames_.push_back("MCTruthResolPUEta03"); 	etaBins_.push_back(std::make_pair(0.3,0.4) );
	inputFolderNames_.push_back("MCTruthResolPUEta04"); 	etaBins_.push_back(std::make_pair(0.4,0.5) );
	inputFolderNames_.push_back("MCTruthResolPUEta05"); 	etaBins_.push_back(std::make_pair(0.5,0.6) );
	inputFolderNames_.push_back("MCTruthResolPUEta06"); 	etaBins_.push_back(std::make_pair(0.6,0.7) );
	inputFolderNames_.push_back("MCTruthResolPUEta07"); 	etaBins_.push_back(std::make_pair(0.7,0.8) );
	inputFolderNames_.push_back("MCTruthResolPUEta08"); 	etaBins_.push_back(std::make_pair(0.8,0.9) );
	inputFolderNames_.push_back("MCTruthResolPUEta09"); 	etaBins_.push_back(std::make_pair(0.9,1.0) );
	inputFolderNames_.push_back("MCTruthResolPUEta10"); 	etaBins_.push_back(std::make_pair(1.0,1.1) );
	inputFolderNames_.push_back("MCTruthResolPUEta11"); 	etaBins_.push_back(std::make_pair(1.1,1.2) );
	inputFolderNames_.push_back("MCTruthResolPUEta12"); 	etaBins_.push_back(std::make_pair(1.2,1.3) );
	inputFolderNames_.push_back("MCTruthResolPUEta13"); 	etaBins_.push_back(std::make_pair(1.3,1.4) );
	inputFolderNames_.push_back("MCTruthResolPUEta14"); 	etaBins_.push_back(std::make_pair(1.4,1.5) );
	inputFolderNames_.push_back("MCTruthResolPUEta15"); 	etaBins_.push_back(std::make_pair(1.5,1.6) );
	inputFolderNames_.push_back("MCTruthResolPUEta16"); 	etaBins_.push_back(std::make_pair(1.6,1.7) );
	inputFolderNames_.push_back("MCTruthResolPUEta17"); 	etaBins_.push_back(std::make_pair(1.7,1.8) );
	inputFolderNames_.push_back("MCTruthResolPUEta18"); 	etaBins_.push_back(std::make_pair(1.8,1.9) );
	inputFolderNames_.push_back("MCTruthResolPUEta19"); 	etaBins_.push_back(std::make_pair(1.9,2.0) );
	inputFolderNames_.push_back("MCTruthResolPUEta20"); 	etaBins_.push_back(std::make_pair(2.0,2.1) );
	inputFolderNames_.push_back("MCTruthResolPUEta21"); 	etaBins_.push_back(std::make_pair(2.1,2.2) );
	inputFolderNames_.push_back("MCTruthResolPUEta22"); 	etaBins_.push_back(std::make_pair(2.2,2.3) );
	inputFolderNames_.push_back("MCTruthResolPUEta23"); 	etaBins_.push_back(std::make_pair(2.3,2.4) );
	inputFolderNames_.push_back("MCTruthResolPUEta24"); 	etaBins_.push_back(std::make_pair(2.4,2.5) );
	inputFolderNames_.push_back("MCTruthResolPUEta25"); 	etaBins_.push_back(std::make_pair(2.5,2.8) );
	inputFolderNames_.push_back("MCTruthResolPUEta26"); 	etaBins_.push_back(std::make_pair(2.8,3.0) );
	inputFolderNames_.push_back("MCTruthResolPUEta27"); 	etaBins_.push_back(std::make_pair(3.0,3.2) );
	inputFolderNames_.push_back("MCTruthResolPUEta28"); 	etaBins_.push_back(std::make_pair(3.2,5.0) );
	inputTH2Names_.push_back("GenJetResponseVsGenJetPt_Z2star_L2L3_NPU0"); inputTH2Names_.push_back("GenJetResponseVsGenJetPt_Z2star_L2L3_NPU1"); inputTH2Names_.push_back("GenJetResponseVsGenJetPt_Z2star_L2L3_NPU2"); inputTH2Names_.push_back("GenJetResponseVsGenJetPt_Z2star_L2L3_NPU3"); 
	lowPuBinEdeges_.push_back(std::make_pair(0,10) );
	lowPuBinEdeges_.push_back(std::make_pair(10,20) );
	lowPuBinEdeges_.push_back(std::make_pair(20,30) );
	lowPuBinEdeges_.push_back(std::make_pair(30,50) );
	lowPuBinEdeges_.push_back(std::make_pair(50,1000) );
	TString cbFittingProcedure("NoFixFitting");
	TH1D *inputTH1D = new TH1D();
	TGraphErrors *inputGraph = new TGraphErrors();
	TH1D *inputTH1DError = new TH1D();
	TH1D *inputTH1DErrorG = new TH1D();
	TGraphErrors *inputGraphError = new TGraphErrors();
		// create result th1 to hold the result fit parameters for each eta bin
	double etaBins[etaBins_.size()+1];
	cout<<"Eta Bin edges for result comparison plots: ";
	for (unsigned int i=0; i< etaBins_.size(); i++)
	{
		etaBins[i]=etaBins_[i].first;
		cout<<etaBins[i]<<", ";
	}	
	etaBins[etaBins_.size()]=etaBins_[etaBins_.size()-1].second;
	cout<<etaBins[etaBins_.size()];
	cout<<endl;
	TH1D *resultParC = new TH1D("resultParC","resultParC; |#eta|;C",etaBins_.size(),etaBins);
	TH1D *resultParS = new TH1D("resultParS","resultParS; |#eta|;S",etaBins_.size(),etaBins);
	TH1D *resultParN = new TH1D("resultParN","resultParN; |#eta|;N",etaBins_.size(),etaBins);
	TH1D *resultParPU = new TH1D("resultParPU","resultParPU; |#eta|;PU",etaBins_.size(),etaBins);
	TH1D *resultParChi2 = new TH1D("resultParChi2","resultParChi2; |#eta|;chi2/NDF",etaBins_.size(),etaBins);
	TCanvas *resultCanvasComparePU10 = new TCanvas("PU10Compare","PU10Compare",200,10,700,500);
	TLegend *resultLegendComparePU10 = new TLegend(0.3,0.7,0.9,0.9);
	
	TCanvas *resultCanvasComparePU30 = new TCanvas("PU30Compare","PU30Compare",200,10,700,500);
	TLegend *resultLegendComparePU30 = new TLegend(0.3,0.7,0.9,0.9);
	cout<<"1";
	PUInclusivePlots(0,"Normalization");
	PUInclusivePlots(1,"Gaus Mean");
	PUInclusivePlots(2,"Gaus Sigma");
	PUInclusivePlots(3,"n left");
	PUInclusivePlots(4,"#alpha left");
	PUInclusivePlots(5,"n right");
	PUInclusivePlots(6,"#alpha right");
	cout<<endl;
	cout<<"2"<<std::endl;

	// Plotting of all fit values for puinclusive each inputFolderNames_ within one histogramm
	for (unsigned int i=0; i < inputFolderNames_.size(); i++)
	{
		TString EtaBin("ErrorNoEtaBinNameSelected");
//		if(i==0) EtaBin="[0,0.5]";
//		if(i==1) EtaBin="[0.5,1.1]";
//		if(i==2) EtaBin="[1.1,1.7]";
//		if(i==3) EtaBin="[1.7,2.3]";
//		if(i==4) EtaBin="[2.3,5.0]";
		TString etaTString("");
		TString npuTString("");
		TString tempTstring("");
		etaTString=Form (" %.1f  %.1f       ", etaBins_[i].first, etaBins_[i].second);
		char buffer[100];
		int n;
		TString TTemp("");
		n = sprintf(buffer,"[%.1f,%.1f]",etaBins_[i].first,etaBins_[i].second);
		TTemp=buffer;
		EtaBin=TTemp;
		vector<TH1D*> th1ds;
		double noPuResult[3] ={0,0,0};
		double noPuResultG[3] ={0,0,0};
		double PuResultG[4] ={0,0,0,0};
		n = sprintf(buffer,"Eta: %d, sigma of gauss core",i);
		TTemp=buffer;
		TCanvas *cCanvas2 = new TCanvas(TTemp_,TTemp_);
		TLegend *lLeg2 = new TLegend(0.3,0.7,0.9,0.9);
		lLeg2->SetFillColor(0);
		TCanvas *cCanvas2G = new TCanvas(TTemp_+"G",TTemp_+"G");
		TLegend *lLeg2G = new TLegend(0.3,0.7,0.9,0.9);
		lLeg2G->SetFillColor(0);

		for(unsigned int ii=0; ii < inputTH2Names_.size();ii++)
		{
			npuTString=Form ("%.0f  %.0f       ", lowPuBinEdeges_[ii].first, lowPuBinEdeges_[ii].second);
			//*FitResults <<setw(10)<<etaTString;
			n = sprintf(buffer,"NPU[%.0f]",meanPUs_[ii]); // mean1
			TTemp = buffer;
			cout<<"Round:"<<inputFolderNames_[i]<<", "<<inputTH2Names_[ii]<<endl;
			TTemp_ = "CB_FitValue2_"+inputFolderNames_[i] + "_" + inputTH2Names_[ii];
			cout<<"inputTh1D"<<TTemp_<<endl;
			inputTH1D =      (TH1D*) ((TDirectory*)((TDirectory*) ((TDirectory*) inF_->Get(inputFolderNames_[i]) )->Get(inputTH2Names_[ii]))->Get(cbFittingProcedure))->Get(TTemp_)->Clone();
			unsigned int nbinsX = inputTH1D->GetNbinsX();
			cbForFourDMinimizer_.push_back(inputTH1D);
			//std::cout<<"Number of nBinsX"<<nbinsX<<std::endl;
			Double_t *xbins = new Double_t[nbinsX];
			inputTH1D->GetXaxis()->GetLowEdge(xbins);
			//for(int iii=0; iii< nbinsX;iii++)std::cout<<"bin["<<iii<<"]: "<<xbins[iii];
			//std::cout<<std::endl;
			inputTH1DError = new TH1D("ErrorBand"+TTemp, ".95 conf.band",nbinsX-1, xbins );
			inputTH1DErrorG = new TH1D("ErrorBand"+TTemp, ".95 conf.band",nbinsX-1, xbins );
			th1ds.push_back( (TH1D*) ((TDirectory*)((TDirectory*) ((TDirectory*) inF_->Get(inputFolderNames_[i]) )->Get(inputTH2Names_[ii]))->Get(cbFittingProcedure))->Get(TTemp_)->Clone() );
			TCanvas *cCanvas = new TCanvas(TTemp_,TTemp_);
			TLegend *lLeg = new TLegend(0.3,0.7,0.9,0.9);
			lLeg->SetFillColor(0);
			cCanvas->cd();
			cCanvas->SetName(TTemp_);
			cCanvas->SetTitle("Crystal Ball: Gaus Sigma");
			inputTH1D->SetMarkerColor(1);
			inputTH1D->SetMarkerSize(2);
			inputTH1D->SetLineColor(ii+1);
			inputTH1D->SetLineWidth(0);
			inputTH1D->SetMarkerStyle(20);
			inputTH1D->SetMarkerSize(1);
			inputTH1D->SetMarkerColor(ii+1);
			lLeg->AddEntry(inputTH1D,TTemp,"f");
			if(ii==0)
			{
				TF1 *fNoPUMinus = new TF1 ("fNoPUMinus",fSigma3TermsMinus, ptRangeLow_,ptRangeHigh_,3);
				fNoPUMinus->SetParameters(0.03,1.8,0.7);
				fNoPUMinus->SetLineColor(12);
				fNoPUMinus->SetParLimits(0,0,100);
				fNoPUMinus->SetParLimits(1,0,100);
				fNoPUMinus->SetParLimits(2,-100,100);
				inputTH1D->GetXaxis()->SetRange(inputTH1D->GetXaxis()->FindBin(ptRangeLow_),inputTH1D->GetXaxis()->FindBin(ptRangeHigh_));
				inputTH1D->SetTitle("Gauss #sigma Eta:"+EtaBin);
				inputTH1D->SetMaximum(inputTH1D->GetMaximum()*1.6);
				inputTH1D->Fit(fNoPUMinus,"QBRI+");
				(TVirtualFitter::GetFitter())->GetConfidenceIntervals(inputTH1DError);
				inputTH1DError->GetXaxis()->SetRange(inputTH1D->GetXaxis()->FindBin(ptRangeLow_),inputTH1D->GetXaxis()->FindBin(ptRangeHigh_));
				n = sprintf(buffer,"y=sqrt([%.3f]^2 + [%.3f]^2/x + [%.4f]^2/x^2 ",fNoPUMinus->GetParameter(0),fNoPUMinus->GetParameter(1),fNoPUMinus->GetParameter(2)); // mean1
				TTemp = TTemp+","+buffer;
				lLeg->AddEntry(fNoPUMinus,TTemp,"f");
				lLeg2->AddEntry(fNoPUMinus,TTemp,"f");
				n = sprintf(buffer,"Error=%.3f, %.3f, %.4f",fNoPUMinus->GetParError(0),fNoPUMinus->GetParError(1),fNoPUMinus->GetParError(2)); // mean1
				TTemp = buffer;
				lLeg->AddEntry(fNoPUMinus,TTemp,"f");
				noPuResult[0]= fNoPUMinus->GetParameter(0);
				noPuResult[1]= fNoPUMinus->GetParameter(1);
				noPuResult[2]=fNoPUMinus->GetParameter(2);

			}
			if(ii>0)
			{
				TF1 *fPileup = new TF1 ("PileUpIncluded",fSigma4TermsMinus,ptRangeLow_,ptRangeHigh_,4);
				fPileup->SetParameters(noPuResult[0],noPuResult[1],noPuResult[2],noPuResult[2]);
				fPileup->FixParameter(0,noPuResult[0]);
				fPileup->FixParameter(1,noPuResult[1]);
				fPileup->FixParameter(2,noPuResult[2]);
				fPileup->SetParLimits(3,0,100);
				fPileup->SetLineColor(ii+1);
				inputTH1D->GetXaxis()->SetRange(inputTH1D->GetXaxis()->FindBin(ptRangeLow_),inputTH1D->GetXaxis()->FindBin(ptRangeHigh_));
				inputTH1D->SetMaximum(inputTH1D->GetMaximum()*1.6);
				inputTH1D->Fit(fPileup,"QRBI+");
				(TVirtualFitter::GetFitter())->GetConfidenceIntervals(inputTH1DError);
				inputTH1DError->GetXaxis()->SetRange(inputTH1D->GetXaxis()->FindBin(ptRangeLow_),inputTH1D->GetXaxis()->FindBin(ptRangeHigh_));
				n = sprintf(buffer,"y= sqrt([%.3f]^2 + [%.3f]^2/x + [%.4f]^2/x^2 + [%.3f]^2/x^2",fPileup->GetParameter(0),fPileup->GetParameter(1),fPileup->GetParameter(2),fPileup->GetParameter(3)); // mean1
				TTemp = TTemp +","+ buffer;
				lLeg->AddEntry(fPileup,TTemp,"f");
				lLeg2->AddEntry(fPileup,TTemp,"f");
				n = sprintf(buffer,"Error=%.3f, %.3f, %.4f, %.3f,",fPileup->GetParError(0),fPileup->GetParError(1),fPileup->GetParError(2),fPileup->GetParError(3)); // mean1
				TTemp = buffer;
				lLeg->AddEntry(fPileup,TTemp,"f");
			}
			cCanvas->cd();
			gStyle->SetErrorX(0);
			inputTH1D->SetStats(kFALSE);
			inputTH1D->Draw("P");
			inputTH1DError->SetStats(kFALSE);
			inputTH1DError->SetFillStyle(0);
		//	inputTH1DError->SetLineStyle(2);
			inputTH1DError->SetFillColor(ii+1);
		//	inputTH1DError->SetLineColor(ii+1);
			if(ii==0) inputTH1DError->SetFillColor(12);
			inputTH1DError->Draw("e3 same");
			cCanvas->SetLogy();
			cCanvas->SetLogx();
			lLeg->Draw();
			cCanvas->Update();
			outF_->cd();
			cCanvas->Write();
			cCanvas2->cd();
			if(ii==0)
			{
				inputTH1D->Draw("P");
				inputTH1DError->Draw("e3 same");
			}
			if(ii>0)
			{
				inputTH1D->Draw("PSame");
				inputTH1DError->Draw("e3 same");
			}
			cCanvas2->Update();
			delete cCanvas;
			delete lLeg;
			// tgraph stuff
			n = sprintf(buffer,"NPU[%.0f]",meanPUs_[ii]); // mean1
			TTemp = buffer;
			TTemp_ = "CB_FitValue2_"+inputFolderNames_[i] + "_" + inputTH2Names_[ii];
			inputGraph = TGraphCreator(inputTH1D, i,ii,ptMeanFile);
			inputGraph->SetTitle("Gauss #sigma Eta:"+EtaBin);
			inputGraph->SetMarkerColor(ii+1);
			std::cout<<std::endl;
			TCanvas *cCanvasG = new TCanvas(TTemp_+"G",TTemp_+"G",200,10,700,500);
			TLegend *lLegG = new TLegend(0.3,0.7,0.9,0.9);
			lLegG->SetFillColor(0);
			cCanvasG->cd();
			cCanvasG->SetName(TTemp_+"G");
			cCanvasG->SetTitle("Crystal Ball: Gaus Sigma(Graph)");
			lLegG->AddEntry(inputGraph,TTemp+"G","f");
			if(ii==0)
			{
				TF1 *fNoPUMinus = new TF1 ("fNoPUMinus",fSigma3TermsMinus, ptRangeLow_,ptRangeHigh_,3);
				fNoPUMinus->SetParameters(noPuResult[0],noPuResult[1],noPuResult[2]);
				fNoPUMinus->SetLineColor(12);
			//	fNoPUMinus->SetParLimits(noPuResult[0],noPuResult[0]*0.3,noPuResult[0]*2);
			//	fNoPUMinus->SetParLimits(noPuResult[1],noPuResult[1]*0.3,noPuResult[1]*2);
			//	fNoPUMinus->SetParLimits(noPuResult[2],noPuResult[2]*0.3,noPuResult[2]*2);
				// do the advanced fitting
				ROOT::Math::WrappedMultiTF1 fSigma3TermsMinus(*fNoPUMinus,1);
				ROOT::Fit::DataOptions opt; 
				ROOT::Fit::DataRange rangeTGraph;
				rangeTGraph.SetRange(ptRangeLow_,ptRangeHigh_);
				ROOT::Fit::BinData dataTGraph(opt,rangeTGraph);
				ROOT::Fit::FillData(dataTGraph, inputGraph);
				ROOT::Fit::Chi2Function chi2_CB(dataTGraph, fSigma3TermsMinus);
				ROOT::Fit::Fitter fitter;
				fitter.Config().SetParamsSettings(3,noPuResult);
				//fitter.Config().MinimizerOptions().SetMinimizerType("Genetic");
				//fitter.Config().SetMinimizer("Minuit2","Migrad");
				//fitter.Config().SetMinimizer("GSLMultiMin","BFGS");// bad results
			//	fitter.Config().MinimizerOptions().SetMinimizerType("GSLSimAn"); // bad result
				//fitter.Config().MinimizerOptions().SetMinimizerType("Fumili"); // bad result
				fitter.Config().SetMinimizer("Minuit","Migrad");
				//fitter.Config().MinimizerOptions().SetMinimizerType("Minuit2","Migrad");
				if(noPuResult[0]>0)fitter.Config().ParSettings(0).SetLimits(noPuResult[0]* 0.4,noPuResult[0]*1.8);
				else fitter.Config().ParSettings(0).SetLimits(noPuResult[0]* 1.8,noPuResult[0]*0.4);
				if(noPuResult[1]>0)fitter.Config().ParSettings(1).SetLimits(noPuResult[1]* 0.4,noPuResult[1]*1.8);
				else fitter.Config().ParSettings(1).SetLimits(noPuResult[1]* 1.8,noPuResult[1]*0.4);
				if(noPuResult[2]>0)fitter.Config().ParSettings(2).SetLimits(noPuResult[2]* 0.4,noPuResult[2]*1.8);
				else fitter.Config().ParSettings(2).SetLimits(noPuResult[2]* 1.8,noPuResult[2]*0.4);
			
				//fitter.Config().ParSettings(0).SetLimits(0,10);
				fitter.Config().ParSettings(0).SetStepSize(0.01);
				//fitter.Config().ParSettings(1).SetLimits(0,10);
				fitter.Config().ParSettings(1).SetStepSize(0.01);
				//fitter.Config().ParSettings(2).SetLimits(-10,10);
				fitter.Config().ParSettings(2).SetStepSize(0.01);
				fitter.FitFCN(chi2_CB,0,dataTGraph.Size(),true);
				ROOT::Fit::FitResult resultprint = fitter.Result();
				resultprint.Print(std::cout);
				fNoPUMinus->SetFitResult( resultprint, 0);
				fNoPUMinus->SetNDF(resultprint.Ndf());
				fNoPUMinus->SetChisquare(resultprint.Chi2());
				inputGraph->GetListOfFunctions()->Add(fNoPUMinus);
				
				//inputGraph->Fit(fNoPUMinus,"QMBR+");
				inputGraphError = new TGraphErrors(inputGraph->GetN());
				inputGraphError->SetTitle("Gauss #sigma Eta:"+EtaBin);
				for (int iii=0; iii<inputGraph->GetN();iii++)
				{
					inputGraphError->SetPoint(iii, inputGraph->GetX()[iii], 0);
				//std::cout<<"GraphValues ["<<ii<<"]: "<<inputGraph->GetX()[ii]<<", ";
				}
				(TVirtualFitter::GetFitter())->GetConfidenceIntervals(inputTH1DErrorG);
				inputTH1DErrorG->GetXaxis()->SetRange(inputTH1DErrorG->GetXaxis()->FindBin(ptRangeLow_)-1,inputTH1DErrorG->GetXaxis()->FindBin(ptRangeHigh_)+1);
				//inputTH1DErrorG->GetXaxis()->SetRange(inputTH1D->GetXaxis()->FindBin(ptRangeLow_),inputTH1D->GetXaxis()->FindBin(ptRangeHigh_));
				n = sprintf(buffer,"y=sqrt([%.3f]^2 + [%.3f]^2/x + [%.4f]^2/x^2 ) ",fNoPUMinus->GetParameter(0),fNoPUMinus->GetParameter(1),fNoPUMinus->GetParameter(2)); // mean1
				TTemp = TTemp+","+buffer;
				lLegG->AddEntry(fNoPUMinus,TTemp,"f");
				lLeg2G->AddEntry(fNoPUMinus,TTemp,"f");
				n = sprintf(buffer,"Error=%.3f, %.3f, %.4f",fNoPUMinus->GetParError(0),fNoPUMinus->GetParError(1),fNoPUMinus->GetParError(2)); // mean1
				TTemp = buffer;
				lLegG->AddEntry(fNoPUMinus,TTemp,"f");
				noPuResultG[0]= fNoPUMinus->GetParameter(0);
				noPuResultG[1]= fNoPUMinus->GetParameter(1);
				noPuResultG[2]=fNoPUMinus->GetParameter(2);
				PuResultG[0]=fNoPUMinus->GetParameter(0);
				PuResultG[1]=fNoPUMinus->GetParameter(1);
				PuResultG[2]=fNoPUMinus->GetParameter(2);
				PuResultG[3]=fNoPUMinus->GetParameter(2);
				tempTstring=Form ("C=%.5f   S=%.5f   N1=%.5f   N2=0.00000 \n", fNoPUMinus->GetParameter(0), fNoPUMinus->GetParameter(1),fNoPUMinus->GetParameter(2));
			}
			if(ii>0)
			{
				TF1 *fPileup = new TF1 ("PileUpIncluded",fSigma4TermsMinus,ptRangeLow_,ptRangeHigh_,4);
				fPileup->SetParameters(noPuResultG[0],noPuResultG[1],noPuResultG[2],noPuResultG[2]);
				fPileup->FixParameter(0,noPuResultG[0]);
				fPileup->FixParameter(1,noPuResultG[1]);
				fPileup->FixParameter(2,noPuResultG[2]);
				fPileup->SetParLimits(3,0,100);
				fPileup->SetLineColor(ii+1);
				ROOT::Math::WrappedMultiTF1 fSigma4TermsMinus(*fPileup,1);
				ROOT::Fit::DataOptions opt; 
				ROOT::Fit::DataRange rangeTGraph;
				rangeTGraph.SetRange(ptRangeLow_,ptRangeHigh_);
				ROOT::Fit::BinData dataTGraph(opt,rangeTGraph);
				ROOT::Fit::FillData(dataTGraph, inputGraph);
				ROOT::Fit::Chi2Function chi2_CB(dataTGraph, fSigma4TermsMinus);
				ROOT::Fit::Fitter fitter;
				fitter.Config().SetParamsSettings(4,PuResultG);
				//fitter.Config().MinimizerOptions().SetMinimizerType("Genetic");
				//fitter.Config().SetMinimizer("Minuit2","Migrad");
				//fitter.Config().SetMinimizer("GSLMultiMin","BFGS"); // bad results
				//fitter.Config().MinimizerOptions().SetMinimizerType("GSLSimAn");// bad result
				//fitter.Config().MinimizerOptions().SetMinimizerType("Fumili");  //bad result
				fitter.Config().SetMinimizer("Minuit","Migrad");
				//fitter.Config().MinimizerOptions().SetMinimizerType("Minuit2","Migrad");
			/*	if(noPuResult[0]>0)fitter.Config().ParSettings(0).SetLimits(noPuResult[0]* 0.4,noPuResult[0]*1.8);
				else fitter.Config().ParSettings(0).SetLimits(noPuResult[0]* 1.8,noPuResult[0]*0.4);
				if(noPuResult[1]>0)fitter.Config().ParSettings(1).SetLimits(noPuResult[1]* 0.4,noPuResult[1]*1.8);
				else fitter.Config().ParSettings(1).SetLimits(noPuResult[1]* 1.8,noPuResult[1]*0.4);
				if(noPuResult[2]>0)fitter.Config().ParSettings(2).SetLimits(noPuResult[2]* 0.4,noPuResult[2]*1.8);
				else fitter.Config().ParSettings(2).SetLimits(noPuResult[2]* 1.8,noPuResult[2]*0.4);
			*/
				fitter.Config().ParSettings(0).Fix();
				fitter.Config().ParSettings(1).Fix();
				fitter.Config().ParSettings(2).Fix();
				fitter.Config().ParSettings(3).SetLimits(0,20);
				fitter.Config().ParSettings(3).SetStepSize(0.01);
				fitter.FitFCN(chi2_CB,0,dataTGraph.Size(),true);
				ROOT::Fit::FitResult resultprint = fitter.Result();
				resultprint.Print(std::cout);
				fPileup->SetFitResult( resultprint, 0);
				fPileup->SetNDF(resultprint.Ndf());
				fPileup->SetChisquare(resultprint.Chi2());
				inputGraph->GetListOfFunctions()->Add(fPileup);
				inputGraphError = new TGraphErrors(inputGraph->GetN());
				inputGraphError->SetTitle("Gauss #sigma Eta:"+EtaBin);
				for (int iii=0; iii<inputGraph->GetN();iii++)
				{
					inputGraphError->SetPoint(iii, inputGraph->GetX()[iii], 0);
				//std::cout<<"GraphValues ["<<ii<<"]: "<<inputGraph->GetX()[ii]<<", ";
				}
				(TVirtualFitter::GetFitter())->GetConfidenceIntervals(inputTH1DErrorG);
				inputTH1DErrorG->GetXaxis()->SetRange(inputTH1DErrorG->GetXaxis()->FindBin(ptRangeLow_)-1,inputTH1DErrorG->GetXaxis()->FindBin(ptRangeHigh_)+1);
				//inputTH1DErrorG->GetXaxis()->SetRange(inputTH1D->GetXaxis()->FindBin(ptRangeLow_),inputTH1D->GetXaxis()->FindBin(ptRangeHigh_));
				n = sprintf(buffer,"y= sqrt([%.3f]^2 + [%.3f]^2/x + [%.4f]^2/x^2 + [%.3f]^2/x^2)",fPileup->GetParameter(0),fPileup->GetParameter(1),fPileup->GetParameter(2),fPileup->GetParameter(3)); // mean1
				TTemp = TTemp +","+ buffer;
				lLegG->AddEntry(fPileup,TTemp,"f");
				lLeg2G->AddEntry(fPileup,TTemp,"f");
				n = sprintf(buffer,"Error=%.3f, %.3f, %.4f, %.3f,",fPileup->GetParError(0),fPileup->GetParError(1),fPileup->GetParError(2),fPileup->GetParError(3)); // mean1
				TTemp = buffer;
				lLegG->AddEntry(fPileup,TTemp,"f");
				tempTstring=Form ("C=%.5f   S=%.5f   N1=%.5f   N2=%.5f  \n", fPileup->GetParameter(0), fPileup->GetParameter(1),fPileup->GetParameter(2),fPileup->GetParameter(3));
				
			}
			cCanvasG->cd();
//			inputGraphError->SetLineColor(kRed);
//			inputGraphError->Draw("e3 AP");
			inputGraph->GetXaxis()->SetLimits(ptRangeLow_*0.9, ptRangeHigh_*1.1);
			inputGraph->SetMarkerColor(ii+1);
			inputGraph->SetMarkerSize(1);
			inputGraph->SetMarkerStyle(20);
			inputGraph->SetLineColor(ii+1);
			inputGraph->Draw("APSame");
		//	inputGraph->SetStats(kFALSE);
			inputGraphError->SetFillStyle(0);
		//	inputTH1DError->SetLineStyle(2);
			inputTH1DError->SetFillColor(ii+1);
		//	inputTH1DError->SetLineColor(ii+1);
			if(ii==0) inputGraphError->SetFillColor(12);
			inputGraphError->GetXaxis()->SetLimits(ptRangeLow_*0.9, ptRangeHigh_*1.1);
			//inputGraphError->Draw("3Same");
			inputTH1DErrorG->SetStats(kFALSE);
			inputTH1DErrorG->SetFillStyle(0);
		//	inputTH1DError->SetLineStyle(2);
			inputTH1DErrorG->SetFillColor(ii+1);
		//	inputTH1DError->SetLineColor(ii+1);
			if(ii==0) inputTH1DErrorG->SetFillColor(12);
			inputTH1DErrorG->Draw("e3 same");
			cCanvasG->SetLogy();
			cCanvasG->SetLogx();
			lLegG->Draw();
			cCanvasG->Update();
			outF_->cd();
			cCanvasG->Write();
			cCanvas2G->cd();
			if(ii==0)
			{
				inputGraph->Draw("AP");
				//inputGraphError->Draw("3Same");
				inputTH1DErrorG->Draw("e3 same");
			}
			if(ii>0)
			{
				inputGraph->Draw("PSame");
				//inputGraphError->Draw("3Same");
				inputTH1DErrorG->Draw("e3 same");
			}
			cCanvas2G->Update();
			//*FitResults<<tempTstring; ///////////////////////////////// use to show results in result file
			delete cCanvasG;
			delete lLegG;
		}
		// starting combined minimization

		if(cbForFourDMinimizer_.size()==0) 
		{
			std::cout<<"Error cbForFourDMinimizer_ is empty!"<<std::endl;
			*FitResults<<"Error cbForFourDMinimizer_ is empty!\n";
		}
		if(meanPUs_.size()==0) 
		{
			std::cout<<"Error meanPUs_ is empty!"<<std::endl;
			*FitResults<<"Error meanPUs_ is empty!\n";
		}
		//ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
		ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Genetic");
		//min->SetMaxFunctionCalls(1000000);
//		min->SetTolerance(0.001);
		min->SetPrintLevel(1);
		ROOT::Math::Functor f(&fourDminimizer,4);
		min->SetFunction(f);
		min->SetFixedVariable(0,"x",noPuResult[0]);
		min->SetFixedVariable(1,"y",noPuResult[1]);
		min->SetLimitedVariable(2,"z",-1,0.00001,-4.5,4.5);
		min->SetLimitedVariable(3,"w",4,0.0001,0.001,2.5);
		min->Minimize();
		const double *pars2 = min->X();
		std::cout << "Fixed:" << pars2[0] << "," << pars2[1] << "," << pars2[2] << ","<< pars2[3]<< std::endl;
		//*FitResults<<"\nCombined fixed fit Result: C="<<pars2[0]<<" S:"<<pars2[1]<<" N:"<<pars2[2]<<" PU"<<pars2[3]<<"\n";

		//ROOT::Math::Minimizer* min2 = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
		ROOT::Math::Minimizer* min2 = ROOT::Math::Factory::CreateMinimizer("Genetic");
		//min->SetMaxFunctionCalls(1000000);
//		min->SetTolerance(0.001);
		min2->SetPrintLevel(1);
		ROOT::Math::Functor f2(&fourDminimizer,4);
		min2->SetFunction(f);
		min2->SetLimitedVariable(0,"x",pars2[0],0.00001,pars2[0]*0.8,pars2[0]*1.2);
		min2->SetLimitedVariable(1,"y",pars2[1],0.00001,pars2[1]*0.8,pars2[1]*1.2);
		if(pars2[2]>0)min2->SetLimitedVariable(2,"z",pars2[2],0.00001,pars2[2]*0.8,pars2[2]*1.2);
		if(pars2[2]<0)min2->SetLimitedVariable(2,"z",pars2[2],0.00001,pars2[2]*1.2,pars2[2]*0.8);
		if(pars2[3]>0)min2->SetLimitedVariable(3,"w",pars2[3],0.00001,pars2[3]*0.8,pars2[3]*1.2);
		if(pars2[3]<0)min2->SetLimitedVariable(3,"w",pars2[3],0.00001,pars2[3]*1.2,pars2[3]*0.8);
		min2->Minimize();
		const double *pars3 = min2->X();
		std::cout << "Free Genetic:" << pars3[0] << "," << pars3[1] << "," << pars3[2] << ","<< pars3[3] << std::endl;
		//*FitResults<<"Combined free Genetic fit Result: C="<<pars3[0]<<" S:"<<pars3[1]<<" N:"<<pars3[2]<<" PU"<<pars3[3]<<"\n";
		double pars4[4];
		for(int ii=0; ii<4;ii++)
		{
			if(abs((pars3[ii]-pars2[ii])/pars2[ii])>0.13)pars4[ii]=pars2[ii];
			else pars4[ii]=pars3[ii];
		}
		ROOT::Math::Minimizer* min3 = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
		//ROOT::Math::Minimizer* min2 = ROOT::Math::Factory::CreateMinimizer("Genetic");
		min3->SetMaxFunctionCalls(1000000);
//		min->SetTolerance(0.001);
		min3->SetPrintLevel(1);
		ROOT::Math::Functor f3(&fourDminimizer,4);
		min3->SetFunction(f3);
		min3->SetLimitedVariable(0,"x",pars4[0],0.00001,pars4[0]*0.8,pars4[0]*1.2);
		min3->SetLimitedVariable(1,"y",pars4[1],0.00001,pars4[1]*0.8,pars4[1]*1.2);
		if(pars4[2]>0)min3->SetLimitedVariable(2,"z",pars4[2],0.00001,pars4[2]*0.8,pars4[2]*1.2);
		if(pars4[2]<0)min3->SetLimitedVariable(2,"z",pars4[2],0.00001,pars4[2]*1.2,pars4[2]*0.8);
		if(pars4[3]>0)min3->SetLimitedVariable(3,"w",pars4[3],0.00001,pars4[3]*0.8,pars4[3]*1.2);
		if(pars4[3]<0)min3->SetLimitedVariable(3,"w",pars4[3],0.00001,pars4[3]*1.2,pars4[3]*0.8);
		min3->Minimize();
		const double *pars = min3->X();
		const double *parsError = min3->Errors();
		std::cout << "Free Minuit2:" << pars[0] << "," << pars[1] << "," << pars[2] << ","<< pars[3] << std::endl;
		//*FitResults<<"Combined free Minuit2 fit Result: C="<<pars[0]<<" S:"<<pars[1]<<" N:"<<pars[2]<<" PU"<<pars[3]<<"\n";
		*FitResults <<setw(7)<<etaTString;
		*FitResults <<setw(3)<<"8"<<setw(6)<<"10"<<setw(6)<<"9999"<<setw(6)<<"0"<<setw(6)<<"120";
		*FitResults<<setw(14)<<pars[0]<<setw(10)<<pars[1]<<setw(10)<<pars[2]<<setw(10)<<pars[3]<<"\n";
		// set results for the reult plot
		resultParC->SetBinContent(i+1,pars[0]);
		resultParC->SetBinError(i+1,parsError[0]);
		resultParS->SetBinContent(i+1,pars[1]);
		resultParS->SetBinError(i+1,parsError[1]);
		resultParN->SetBinContent(i+1,pars[2]);
		resultParN->SetBinError(i+1,parsError[2]);
		resultParPU->SetBinContent(i+1,pars[3]);
		resultParPU->SetBinError(i+1,parsError[3]);
		if(min3->NFree()>0)resultParChi2->SetBinContent(i+1,min3->MinValue()/min3->NFree());
		if(min3->NFree()>0)resultParChi2->SetBinError(i+1,1/sqrt(min3->NFree()));
		TF1 *twoSigmaPU10 = new TF1 ("twoSigmaPU10",TwoDSigmaFitFunction, ptRangeLow_,ptRangeHigh_,5);
		twoSigmaPU10->FixParameter(0,pars[0]);
		twoSigmaPU10->FixParameter(1,pars[1]);
		twoSigmaPU10->FixParameter(2,pars[2]);
		twoSigmaPU10->FixParameter(3,pars[3]);
		twoSigmaPU10->FixParameter(4,10);
		twoSigmaPU10->SetLineColor(i+1);
		twoSigmaPU10->SetLineWidth(3);
		resultLegendComparePU10->AddEntry(twoSigmaPU10,etaTString);
		resultCanvasComparePU10->cd();
		if(i==0)twoSigmaPU10->DrawCopy();
		if(i>0) twoSigmaPU10->DrawCopy("Same");
		
		TF1 *twoSigmaPU30 = new TF1 ("twoSigmaPU30",TwoDSigmaFitFunction, ptRangeLow_,ptRangeHigh_,5);
		twoSigmaPU30->FixParameter(0,pars[0]);
		twoSigmaPU30->FixParameter(1,pars[1]);
		twoSigmaPU30->FixParameter(2,pars[2]);
		twoSigmaPU30->FixParameter(3,pars[3]);
		twoSigmaPU30->FixParameter(4,30);
		twoSigmaPU30->SetLineColor(i+1);
		twoSigmaPU30->SetLineWidth(3);
		resultLegendComparePU30->AddEntry(twoSigmaPU30,etaTString);
		resultCanvasComparePU30->cd();
		if(i==0)twoSigmaPU30->DrawCopy();
		if(i>0) twoSigmaPU30->DrawCopy("Same");

		// draw the results
		TCanvas *cCanvasCombined = new TCanvas("Eta_"+etaTString,"Eta_"+etaTString,200,10,700,500);
		TLegend *lLegCombined = new TLegend(0.3,0.7,0.9,0.9);
		cCanvasCombined->SetLogx();
		cCanvasCombined->SetLogy();
		for(unsigned int ii=0; ii <cbForFourDMinimizer_.size(); ii++)
		{
			npuTString=Form ("%.0f  %.0f       ", lowPuBinEdeges_[ii].first, lowPuBinEdeges_[ii].second);
			TCanvas *cCanvas3 = new TCanvas("Eta_"+etaTString+"PU_"+npuTString,"Eta_"+etaTString+"PU_"+npuTString,200,10,700,500);
			TLegend *lLeg3 = new TLegend(0.3,0.7,0.9,0.9);
			lLeg3->SetFillColor(0);
			TF1 *twoSigma = new TF1 ("twoSigma",TwoDSigmaFitFunction, ptRangeLow_,ptRangeHigh_,5);
			twoSigma->FixParameter(0,pars[0]);
			twoSigma->FixParameter(1,pars[1]);
			twoSigma->FixParameter(2,pars[2]);
			twoSigma->FixParameter(3,pars[3]);
			twoSigma->FixParameter(4,meanPUs_[ii]);
			//cbForFourDMinimizer_[i]->GetListOfFunctions()->Add(twoSigma);
			cCanvas3->cd();
			if(cbForFourDMinimizer_[ii]->GetFunction("fNoPUMinus"))cbForFourDMinimizer_[ii]->GetFunction("fNoPUMinus")->SetBit(TF1::kNotDraw);
			if(cbForFourDMinimizer_[ii]->GetFunction("PileUpIncluded"))cbForFourDMinimizer_[ii]->GetFunction("PileUpIncluded")->SetBit(TF1::kNotDraw);
			cbForFourDMinimizer_[ii]->DrawCopy();
			twoSigma->SetLineColor(ii+1);
			twoSigma->Draw("Same");
			TString tttemp = Form("C=%.4f, S=%.4f, N=%.3f, PUTerm=%.3f, MeanPU=%.0f",pars[0],pars[1],pars[2],pars[3],meanPUs_[ii]);
			lLeg3->AddEntry(twoSigma,tttemp);
			cCanvas3->SetLogx();
			cCanvas3->SetLogy();
			lLeg3->Draw();
			outF_->cd();
			cCanvas3->Write();
			cCanvasCombined->cd();
			if(ii==0) cbForFourDMinimizer_[ii]->DrawCopy();
			else cbForFourDMinimizer_[ii]->DrawCopy("Same");
			twoSigma->Draw("Same");
			if(ii>0)tttemp = Form("MeanPU=%.0f",meanPUs_[ii]);
			lLegCombined->AddEntry(twoSigma,tttemp);
			delete cCanvas3;
			delete lLeg3;
			
		}
		lLegCombined->Draw();
		outF_->cd();
		cCanvasCombined->Write();
		delete cCanvasCombined;
		delete lLegCombined;

		// clearing the input th1collection
		cbForFourDMinimizer_.clear();
		
		cout<<"1";
		cCanvas2->cd();
		cCanvas2->SetLogx();
		lLeg2->Draw();
		cCanvas2->Update();
		cCanvas2G->cd();
		cCanvas2G->SetLogx();
		lLeg2G->Draw();
		cCanvas2G->Update();
		outF_->cd();
		cCanvas2->Write();
		cCanvas2G->Write();
		delete inputTH1D;
		delete inputGraph;
		delete inputTH1DError;
		delete inputGraphError;
		delete cCanvas2;
		delete cCanvas2G;
		delete lLeg2;
		delete lLeg2G;
		
		cout<<"2";
		TTemp_ = "EtaBin_"+EtaBin + "-pT_Gaus_Sigma_CB";
		combinedPlots(th1ds, OutPutD_, TTemp_, true);
		OutPutD_->mkdir("PUStudies");
		differencePlots (th1ds, (TDirectory*)OutPutD_->Get("PUStudies"), TTemp_, meanPUs_);
		th1ds.clear();
		
	}
	cout<<"Almost finished"<<endl;
	outF_->cd();
	resultParC->Write();
	cout<<"1";
	resultParS->Write();	cout<<"2";
	resultParN->Write();
	cout<<"3";
	resultParPU->Write();	cout<<"4";
	resultParChi2->Write();
	cout<<"5";
	resultCanvasComparePU10->cd();
	cout<<"5.1";
	resultLegendComparePU10->SetFillColor(0);
	cout<<"5.2";
	resultLegendComparePU10->Draw();
	cout<<"5.3";
	resultCanvasComparePU10->Update();
	cout<<"5.4";
	resultCanvasComparePU10->Write();
	cout<<"6";
	resultCanvasComparePU30->cd();
	resultLegendComparePU30->SetFillColor(0);
	resultLegendComparePU30->Draw();
	resultCanvasComparePU30->Update();
	resultCanvasComparePU30->Write();
	cout<<"7";
	outF_->Write();
	outF_->Close();
	FittingAnalyzer2014ErrorReport->close();
	FitResults->close();
}

TGraphErrors* TGraphCreator(TH1D* inputTH1D, unsigned int etabin, unsigned int puBin, TFile* inputFile)
{
	if(ptBinTH1Names_.size()==0)
	{
		for (int i=0; i<inputTH1D->GetNbinsX();i++)
		{
			ptBinTH1Names_.push_back(Form("MeanPTForBin=%d_%d",(int)inputTH1D->GetBinLowEdge(i),(int)inputTH1D->GetBinLowEdge(i+1)));
		}
	}
	// check if the fit converged if fit status was -2 (failed) the errors will be zero. these points are removed from the interpretation...
	unsigned int notConvergedBins=0;
	for(int i=0; i< inputTH1D->GetNbinsX();i++) if(inputTH1D->GetBinError(i) < 0.00000001) notConvergedBins++;
	TGraphErrors *gr = new TGraphErrors(inputTH1D->GetNbinsX()-notConvergedBins);
	gr->SetName(inputTH1D->GetName());
	TString etaFolder("");
	TString PuFolder("");
	TString ptBin("");
	for (int i=0; i<inputTH1D->GetNbinsX();i++)
	{
		ptBin=ptBinTH1Names_[i];
		if((int)inputTH1D->GetBinLowEdge(i)<0 ) ptBin=Form("MeanPTForBin=0_%d",(int)inputTH1D->GetBinLowEdge(i+1));
		etaFolder=Form ("%.1f_%.1f", etaBins_[etabin].first, etaBins_[etabin].second);
		PuFolder=Form ("%.0f_%.0f", lowPuBinEdeges_[puBin].first, lowPuBinEdeges_[puBin].second);
	//	std::cout<<"Getting ptBin(TH1Name): "<<ptBin<<", from etaFolder: "<<etaFolder<<", from PUFolder: "<<PuFolder<<std::endl;
		
		TH1D *tempTH1D = (TH1D*) ((TDirectory*)((TDirectory*) ((TDirectory*) inputFile->Get("PTBinmeans") )->Get(etaFolder))->Get(PuFolder))->Get(ptBin)->Clone();
		//cout<<"SetPoint i="<<i<<", toX:"<<tempTH1D->GetMean(1)<<std::endl;
		//cout<<"SetPoint i="<<i<<", toXerror:"<<tempTH1D->GetMeanError(1)<<std::endl;
		//cout<<"SetPoint i="<<i<<", inY:"<<inputTH1D->GetBinContent(1)<<std::endl;
		//cout<<"SetPoint i="<<i<<", inYerror:"<<inputTH1D->GetBinError(1)<<std::endl;
		if(inputTH1D->GetBinError(i) < 0.00000001)continue;
		gr->SetPoint(i,tempTH1D->GetMean(1),inputTH1D->GetBinContent(i));
		gr->SetPointError(i,tempTH1D->GetMeanError(1),inputTH1D->GetBinError(i));
		delete tempTH1D;
	}
	return gr;
	
}

void combinedPlots(vector<TH1D*> th1s, TDirectory *outPutFolderSpecial, TString name, bool rename)
{
	
	TCanvas *cCanvas = new TCanvas(name,name);
	TLegend *lLeg = new TLegend(0.6,0.7,0.9,0.9);
	lLeg->SetFillColor(0);
	TH1D *first = new TH1D();
	for (unsigned int iv=0; iv<th1s.size();iv++)
	{
		cCanvas->cd();
		cCanvas->SetName(name);
		cCanvas->SetTitle(name);
		if(iv==0)
		{
			first = th1s[0];
			if (rename) first->SetName("name");
			if (rename) first->SetTitle(name+";p_{T}");
			first->SetMarkerColor(1);
			first->SetMarkerSize(1);
			first->SetLineColor(1);
			first->SetLineWidth(2);
			first->SetMarkerStyle(21);
			lLeg->AddEntry(first,first->GetTitle(),"f");
			first->Draw("P");
		}
		else
		{
			first = th1s[iv];
			first->SetLineColor(iv+4);
			first->SetLineWidth(2);
			first->SetMarkerColor(iv+4);
			first->SetMarkerSize(1);
			first->SetMarkerStyle(21);
			lLeg->AddEntry(first,first->GetTitle(),"f");
			first->Draw("SameP");
		}
	}
	outPutFolderSpecial->cd();
	cCanvas->cd();
	cCanvas->SetLogy();
	cCanvas->SetLogx();
	lLeg->Draw();
	cCanvas->Update();
	cCanvas->Write();
	delete cCanvas;
	delete lLeg;
}

void differencePlots (vector<TH1D*> th1s, TDirectory *outPutFolderSpecial, TString name, vector<double> meanPU_)
{
	TCanvas *cCanvas = new TCanvas(name,name);
	TLegend *lLeg = new TLegend(0.6,0.7,0.9,0.9);
	lLeg->SetFillColor(0);
	vector<TF1*> fitts;
	
	for (unsigned int i=1; i < th1s.size(); i++)
	{
		TF1 *constant =  new TF1("Constant",fConst, ptRangeLow_,th1s[i]->GetXaxis()->GetBinCenter(th1s[i]->GetNbinsX()-1),1);
		fitts.push_back(constant);
		for (int ii=0; ii<th1s[i]->GetNbinsX();ii++)
		{
			double lowPU = th1s[0]->GetBinContent(ii);
			double lowPuError = th1s[0]->GetBinError(ii);
			double highPU = th1s[i]->GetBinContent(ii);
			double highPuError = th1s[i]->GetBinError(ii);
			double result=-1;
			double fehler =-1;
			if ((highPU * highPU - (lowPU * lowPU) )>0)
			{
				result = TMath::Power(std::abs(highPU * highPU - (lowPU * lowPU) ),0.5 );
				fehler =  sqrt(TMath::Power(lowPuError * ( lowPU / TMath::Power(std::abs(highPU * highPU - (lowPU * lowPU) ),0.5 ) ),2) + TMath::Power(highPuError * (highPU / TMath::Power(std::abs(highPU * highPU - (lowPU * lowPU) ),0.5 ) ),2) );
			}
			else 
			{
				result = -1 * TMath::Power(std::abs(highPU * highPU - (lowPU * lowPU) ),0.5 );
				fehler =  sqrt(TMath::Power(lowPuError * ( lowPU / TMath::Power(std::abs(highPU * highPU - (lowPU * lowPU) ),0.5 ) ),2) + TMath::Power(highPuError * (highPU / TMath::Power(std::abs(highPU * highPU - (lowPU * lowPU) ),0.5 ) ),2));
			}
			fehler = fehler * th1s[i]->GetXaxis()->GetBinCenter(ii);
			result = result * th1s[i]->GetXaxis()->GetBinCenter(ii);
			th1s[i]->SetBinContent(ii,result);
			th1s[i]->SetBinError(ii,fehler);
			th1s[i]->SetLineColor(i+4);
			th1s[i]->SetLineWidth(2);
			th1s[i]->SetMarkerColor(i+4);
			th1s[i]->SetMarkerSize(1);
			th1s[i]->SetMarkerStyle(21);
			
		}
		lLeg->AddEntry(th1s[i],th1s[i]->GetTitle(),"f");
		cCanvas->cd();

		fitts[i-1]->SetLineColor(i+4);
		fitts[i-1]->SetRange(15,14000); // range for fits in pt
		th1s[i]->Fit(fitts[i-1],"Rq0N");
		if(i==1) th1s[i]->Draw("PE");
		else th1s[i]->Draw("SamePE");
	}
	outPutFolderSpecial->cd();
	cCanvas->SetLogx();
	lLeg->Draw();
	cCanvas->Update();
	TH1D *SigmaGauss = new TH1D("SigmaGauss"+name,"SigmaGauss"+name+";N(PileUP);p_{T}GeV",101,0,100);
	for (unsigned int i=0; i<fitts.size();i++)
	{
		cCanvas->cd();
		fitts[i]->Draw("Same");
		double mean= meanPU_[i+1]-meanPU_[0]; 
		SigmaGauss->SetBinContent(SigmaGauss->GetXaxis()->FindBin(mean),fitts[i]->GetParameter(0));
		SigmaGauss->SetBinError(SigmaGauss->GetXaxis()->FindBin(mean),fitts[i]->GetParameter(0));
	}
	cCanvas->Write();
	cCanvas->Clear();
/*	SigmaGauss->SetBinContent(SigmaGauss->GetXaxis()->FindBin(8),fitts[0]->GetParameter(0)); // use here average of each pileup bin minus the avarge of the lowest bin (bin0) on which the rest is realtivily calculated to
	SigmaGauss->SetBinError(SigmaGauss->GetXaxis()->FindBin(8),fitts[0]->GetParError(0));
	SigmaGauss->SetBinContent(SigmaGauss->GetXaxis()->FindBin(16.5),fitts[1]->GetParameter(0));
	SigmaGauss->SetBinError(SigmaGauss->GetXaxis()->FindBin(16.5),fitts[1]->GetParError(0));	
	SigmaGauss->SetBinContent(SigmaGauss->GetXaxis()->FindBin(26),fitts[2]->GetParameter(0));
	SigmaGauss->SetBinError(SigmaGauss->GetXaxis()->FindBin(26),fitts[2]->GetParError(0)); */
	TF1 *linear =  new TF1("Linear",fLinear, 1,100,2);
	TLegend *lLeg2 = new TLegend(0.6,0.7,0.9,0.9);
	lLeg2->SetFillColor(0);
	linear->SetParameters(1,1);
	linear->FixParameter(0,0);
	SigmaGauss->Fit(linear,"qNR");
	char buffer[50];
	int n;
	TString TTemp("");
	n = sprintf(buffer,"Linear_f=_%.2f+_%.2f*x",linear->GetParameter(0),linear->GetParameter(1)); // mean1
	TTemp = buffer;
//	lLeg2->AddEntry(linear,TTemp,"f");
	n = sprintf(buffer,"Error=_%.2f, _%.2f",linear->GetParError(0),linear->GetParError(1)); // mean1
	TTemp = buffer;
//	lLeg2->AddEntry(linear,TTemp,"f");
	cout<<"WurzelN fuer"<<name<<endl;
	TF1 *wurzelN = new TF1("WurzelN", fWurzelN,0,100,2);
	wurzelN->SetParameters(1,1);
	wurzelN->SetLineColor(4);
	SigmaGauss->Fit(wurzelN,"qNR+");

	cout<<"WurzelN fix par[0] fuer"<<name<<endl;
	TF1 *wurzelN2 = new TF1("WurzelNFix", fWurzelN,0,100,2);
	wurzelN2->SetParameters(1,1);
	wurzelN2->FixParameter(0,0);
	wurzelN2->SetLineColor(4);
	wurzelN2->SetLineWidth(3);
	SigmaGauss->Fit(wurzelN2,"qR+");
	n = sprintf(buffer,"y=_%.3f*sqrt(x)",wurzelN2->GetParameter(1)); // mean1
	TTemp = buffer;
	lLeg2->AddEntry(wurzelN2,TTemp,"f");
	n = sprintf(buffer,"Error= _%.3f",wurzelN2->GetParError(1)); // mean1
	TTemp = buffer;
	lLeg2->AddEntry(wurzelN2,TTemp,"f");
	cCanvas->cd();
	SigmaGauss->SetMarkerStyle(2);
	SigmaGauss->SetMarkerSize(2);
	SigmaGauss->Draw();
	TF1 *wurzelError = new TF1("Error band",fWurzelN,0,100,2);
	wurzelError->SetParameters(0,wurzelN2->GetParameter(1)+wurzelN2->GetParError(1));
	wurzelError->SetLineColor(2);
	wurzelError->SetLineStyle(3);
	wurzelError->Draw("Same");
	TF1 *wurzelError2 = new TF1("Error band2",fWurzelN,0,100,2);
	wurzelError2->SetParameters(0,wurzelN2->GetParameter(1)-wurzelN2->GetParError(1));
	wurzelError2->SetLineColor(2);
	wurzelError2->SetLineStyle(3);
	wurzelError2->Draw("Same");
	lLeg2->Draw("Same");
	cCanvas->SetName("SigmaGauss"+name);
	cCanvas->Update();
	cCanvas->Write();
	SigmaGauss->Write();
	
	
}

void PUInclusivePlots (unsigned int parameter, TString parName)
{	
	TCanvas *cCanvas = new TCanvas(parName,parName);
	TLegend *lLeg = new TLegend(0.6,0.7,0.9,0.9);
	lLeg->SetFillColor(0);
	outF_->cd();
	TH1D *TH1DTemp;
	vector<TF1*> fitts;
	for( unsigned int i=0; i < inputFolderNames_.size(); i++)
	{
		TF1 *constant =  new TF1("Constant",fConst,30,1000,1);
		fitts.push_back(constant);
		fitts[i]->SetRange(30,1000);
		char buffer[50];
		int n;
		TString TTemp("");
		TString TTemp2("");
		n = sprintf(buffer,"CB_FitValue%d",parameter); // mean1
		TTemp = buffer;
		TTemp = TTemp+"_"+inputFolderNames_[i]+"_PUInclusiveFits";
		TTemp2 = inputFolderNames_[i] + "/PUInclusiveFits";
		//std::cout<<"TTemp currently"<< TTemp<<std::endl;
		cCanvas->cd();
		//(TH1D*) ((TDirectory*) inF_->Get(inputFolderNames_[i]) )->Get(TTemp_)->Clone()
		TH1DTemp = (TH1D*) ((TDirectory*) inF_->Get(TTemp2))->Get(TTemp)->Clone();
		if( parameter==4 || parameter==6) TH1DTemp-> Fit(fitts[i],"0NR");

		if(i==0)
		{
			TH1DTemp->SetMarkerColor(1);
			TH1DTemp->SetName(parName);
			TH1DTemp->SetMarkerSize(1);
			TH1DTemp->SetLineColor(1);
			TH1DTemp->SetLineWidth(2);
			TH1DTemp->SetMarkerStyle(21);
			TH1DTemp->Draw("P");
			if(parameter==4 || parameter==6) 
			{
				((TF1*)TH1DTemp->GetFunction("Constant"))->SetLineColor(1);
				//((TF1*)TH1DTemp->GetFunction("Constant"))->Draw("Same");
			}
		}
		if(i>0)
		{
			TH1DTemp->SetName(parName);
			TH1DTemp->SetLineColor(i+4);
			TH1DTemp->SetLineWidth(2);
			TH1DTemp->SetMarkerColor(i+4);
			TH1DTemp->SetMarkerSize(1);
			TH1DTemp->SetMarkerStyle(21);
			TH1DTemp->Draw("SameP");
			if(parameter==4 || parameter==6) 
			{
				((TF1*)TH1DTemp->GetFunction("Constant"))->SetLineColor(i+4);
				//((TF1*)TH1DTemp->GetFunction("Constant"))->Draw("Same");
			}
		}
		n = sprintf(buffer,"Eta:%d",i); // mean1
		TTemp = buffer + parName;
		if(parameter==4 || parameter==6) 
		{
			n = sprintf(buffer,"Fit, Const:%.2f+-%.2f",((TF1*)TH1DTemp->GetFunction("Constant"))->GetParameter(0),((TF1*)TH1DTemp->GetFunction("Constant"))->GetParError(0));
			TTemp += buffer;
		}
		lLeg->AddEntry(TH1DTemp,TTemp,"f");
		
	}
//	fitts.clear();
	cCanvas->cd();
	if(parameter==0)cCanvas->SetLogy();
	cCanvas->SetLogx();
	cCanvas->SetTitle(parName);
	lLeg->Draw();
	OutPutD_->cd();
	cCanvas->Update();
	cCanvas->Write();
	delete cCanvas;
	delete lLeg;
}
Double_t fSigmaOld (Double_t *x, Double_t *par)
{
	//(sqrt(par[0]/x[0])+(sqrt(par[1])/x[0]))+sqrt(par[2]))+(par[3]/x[0]))+((par[4]/x[0])/sqrt(x[0])
	//(sqrt((sq([0]/x)+(sq([1])/x))+sq([2]))+([3]/x))+(([4]/x)/sqrt(x))
	Double_t result = sqrt( sqrt(par[0]/x[0] + sqrt(par[1])/x[0] ) + sqrt(par[2])) +par[3]/x[0] + (par[4]/x[0])/sqrt(x[0]);
	return result;
}
Double_t fSigmaNoPU (Double_t *x, Double_t *par)
{

	Double_t result = par[1]/TMath::Power(x[0],par[2]) + par[0];
	return result;
}
Double_t fSigmaTwoTerms (Double_t *x, Double_t *par)
{

	Double_t result = par[0] + par[1]/TMath::Power(x[0],par[2]) + par[3]/TMath::Power(x[0],par[4]);
	return result;
}
Double_t fSigma3Terms (Double_t *x, Double_t *par) // sqrt([0]^2 + [1]^2/x + [2]^2/x^2)
{
	Double_t result = sqrt(par[0] *par[0] + par[2]*par[2]/(x[0]*x[0])+ par[1]*par[1]/x[0]);
	return result;
}
Double_t fSigma3TermsMinus (Double_t *x, Double_t *par) // sqrt([0]^2 + [1]^2/x + [2]^2/x^2)
{
	Double_t result = sqrt(par[0] *par[0] + TMath::Sign(1.,par[2]) * par[2]*par[2]/(x[0]*x[0])+ par[1]*par[1]/x[0]);
	return result;
}
Double_t TwoDSigmaFitFunction(Double_t *x,Double_t *par) // par 0=Const, 1=Stoch, 2=Noise, 3=PU, 4=MeanNumber of PileUp in that bin
		//sqrt([0] *[0] + [1]*[1]/x[0] + TMath::Sign(1.,[2]) * ([2])*([2])/(x[0]*x[0]) + ([3]*sqrt([4])) * ([3]*sqrt([4]))/(x[0]*x[0]))
{
	return sqrt(par[0] *par[0] + par[1]*par[1]/x[0] + TMath::Sign(1.,par[2]) * (par[2])*(par[2])/(x[0]*x[0]) + (par[3]*sqrt(par[4])) * (par[3]*sqrt(par[4]))/(x[0]*x[0]));
}
Double_t fSigma4Terms (Double_t *x, Double_t *par) // sqrt([0]^2 + [1]^2/x + [2]^2/x^2 + [3]^2/x^2
{
	Double_t result = sqrt(par[0] *par[0] + par[2]*par[2]/(x[0]*x[0])+ par[1]*par[1]/x[0] + par[3] * par[3]/(x[0] * x[0]));
	return result;
}
Double_t fSigma4TermsMinus (Double_t *x, Double_t *par) // sqrt([0]^2 + [1]^2/x + [2]^2/x^2 + [3]^2/x^2
{
	Double_t result = sqrt(par[0] *par[0] + TMath::Sign(1.,par[2]) * par[2]*par[2]/(x[0]*x[0])+ par[1]*par[1]/x[0] + par[3] * par[3]/(x[0] * x[0]));
	return result;
}
Double_t exp0 (Double_t *x, Double_t *par)
{
	Double_t exp0 = par[0]*TMath::Exp(-par[1] *x[0]) + par[2];
	return exp0;	
}
Double_t fLinear (Double_t *x, Double_t *par)
{
	Double_t result = par[0] + par[1] * x[0];
	return result;
}
Double_t fWurzelN(Double_t *x, Double_t *par)
{
	Double_t result = par[0] + par[1] * sqrt(x[0]);
	return result;
}
Double_t fConst (Double_t *x, Double_t *par)
{
	
	Double_t para = par[0];
	return para;
}

double fourDminimizer(const double *par)  // par[0]=Const, par[1]=Stochastic, par[2]=NoiseTerm, par[3]= term to descripe PU contribution
{
	double result=0;
	double bincontent=0;
	double bincenter=0;
	double binerror=0;
	double puMean=0;
	double funcValue=0;
	if(cbForFourDMinimizer_.size()==0) 
	{
		std::cout<<"fourDminimizer::Error cbForFourDMinimizer_ is empty!"<<std::endl;
		*FitResults<<"fourDminimizer::Error cbForFourDMinimizer_ is empty!\n";
	}
	if(meanPUs_.size()==0) 
	{
		std::cout<<"fourDminimizer::Error meanPUs_ is empty!"<<std::endl;
		*FitResults<<"fourDminimizer::Error meanPUs_ is empty!\n";
	}
	for(unsigned int i=0; i< cbForFourDMinimizer_.size(); i++)
	{
		puMean = meanPUs_[i];
//		*FitResults<<"fourDminimizer::puMean"<<puMean<<"\n";
		for (int ii=0; ii< cbForFourDMinimizer_[i]->GetNbinsX();ii++)
		{
			if(cbForFourDMinimizer_[i]->GetBinError(ii)< 0.0000001) continue;
			if(cbForFourDMinimizer_[i]->GetBinCenter(ii)<ptRangeLow_ || cbForFourDMinimizer_[i]->GetBinCenter(ii)>ptRangeHigh_)continue;
			bincontent = cbForFourDMinimizer_[i]->GetBinContent(ii);
			bincenter = cbForFourDMinimizer_[i]->GetBinCenter(ii);
			binerror = cbForFourDMinimizer_[i]->GetBinError(ii);
			funcValue = sqrt(par[0] *par[0] + par[1]*par[1]/bincenter + TMath::Sign(1.,par[2]) * (par[2])*(par[2])/(bincenter*bincenter) + (par[3]*sqrt(puMean)) * (par[3]*sqrt(puMean))/(bincenter*bincenter));
			result+= (funcValue - bincontent) * (funcValue - bincontent) / ( binerror * binerror );
		}
	}
	return result;
}

