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
double ptRangeLow_, ptRangeHigh_;
bool debug_;
char buffer_[100];
int n_;
TString TTemp_("");
vector<TH1D*> cbValues_;
vector<TString> inputFolderNames_, inputTH2Names_;
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
Double_t fSigma4Terms (Double_t *x, Double_t *par);

void combinedPlots(vector<TH1D*> th1s, TDirectory *outPutFolderSpecial, TString name, bool rename);
void differencePlots (vector<TH1D*> th1s, TDirectory *outPutFolderSpecial, TString name, vector<double> meanPUs);
void PUInclusivePlots (unsigned int parameter, TString parName);
//Double_t cb(Double_t *x, Double_t *par);
//Double_t gauss (Double_t *x, Double_t *par);

#endif



void fittingAnalyzer2014()
{
	gROOT->SetBatch(true);
	ptRangeLow_=30;
	ptRangeHigh_=1200;
	cout<<"Analyzer started..."<<endl;

	inF_ = TFile::Open("Fitting2014.root","UPDATE");
	outF_ = new TFile("Results2014.root","RECREATE");
	outF_->mkdir("Results");
	OutPutD_ = (TDirectory*)outF_->Get("Results");
	TH1D *TruthPUMeans = (TH1D*) inF_->Get("TruthPUMeans")->Clone();
	vector<double> meanPUs;
	for(int i=0; i < 4; i++)
	{
		cout<<"Truth mean Pileups"<<TruthPUMeans->GetBinContent(i+1)<<endl;
		meanPUs.push_back(TruthPUMeans->GetBinContent(i+1));
	}
	inputFolderNames_.push_back("MCTruthResolPUEta0"); inputFolderNames_.push_back("MCTruthResolPUEta1"); inputFolderNames_.push_back("MCTruthResolPUEta2"); inputFolderNames_.push_back("MCTruthResolPUEta3"); inputFolderNames_.push_back("MCTruthResolPUEta4"); 
	inputTH2Names_.push_back("GenJetResponseVsGenJetPt_Z2star_L2L3_NPU0"); inputTH2Names_.push_back("GenJetResponseVsGenJetPt_Z2star_L2L3_NPU1"); inputTH2Names_.push_back("GenJetResponseVsGenJetPt_Z2star_L2L3_NPU2"); inputTH2Names_.push_back("GenJetResponseVsGenJetPt_Z2star_L2L3_NPU3"); 
	TString cbFittingProcedure("NoFixFitting");
	TH1D *inputTH1D = new TH1D();
	PUInclusivePlots(0,"Normalization");
	PUInclusivePlots(1,"Gaus Mean");
	PUInclusivePlots(2,"Gaus Sigma");
	PUInclusivePlots(3,"n left");
	PUInclusivePlots(4,"#alpha left");
	PUInclusivePlots(5,"n right");
	PUInclusivePlots(6,"#alpha right");
	cout<<endl;
	// Plotting of all fit values for puinclusive each inputFolderNames_ within one histogramm
	for (unsigned int i=0; i < inputFolderNames_.size(); i++)
	{
		char buffer[100];
		int n;
		TString TTemp("");
		vector<TH1D*> th1ds;
		double noPuResult[3] ={0,0,0};
		n = sprintf(buffer,"Eta: %d, sigma of gauss core",i);
		TTemp=buffer;
		TCanvas *cCanvas2 = new TCanvas(TTemp_,TTemp_);
		TLegend *lLeg2 = new TLegend(0.3,0.7,0.9,0.9);
		for(unsigned int ii=0; ii < inputTH2Names_.size();ii++)
		{

			
			n = sprintf(buffer,"Eta: %d, NPU: %d",i,ii); // mean1
			TTemp = buffer;
			cout<<"Round:"<<inputFolderNames_[i]<<", "<<inputTH2Names_[ii]<<endl;
			TTemp_ = "CB_FitValue2_"+inputFolderNames_[i] + "_" + inputTH2Names_[ii];
			cout<<"inputTh1D"<<TTemp_<<endl;
			inputTH1D =      (TH1D*) ((TDirectory*)((TDirectory*) ((TDirectory*) inF_->Get(inputFolderNames_[i]) )->Get(inputTH2Names_[ii]))->Get(cbFittingProcedure))->Get(TTemp_)->Clone();
			th1ds.push_back( (TH1D*) ((TDirectory*)((TDirectory*) ((TDirectory*) inF_->Get(inputFolderNames_[i]) )->Get(inputTH2Names_[ii]))->Get(cbFittingProcedure))->Get(TTemp_)->Clone() );
			TCanvas *cCanvas = new TCanvas(TTemp_,TTemp_);
			TLegend *lLeg = new TLegend(0.3,0.7,0.9,0.9);
			cCanvas->cd();
			cCanvas->SetName(TTemp_);
			cCanvas->SetTitle("Crystal Ball: Gaus Sigma");
			inputTH1D->SetMarkerColor(1);
			inputTH1D->SetMarkerSize(2);
			inputTH1D->SetLineColor(ii+1);
			inputTH1D->SetLineWidth(2);
			inputTH1D->SetMarkerStyle(21);
			inputTH1D->SetMarkerSize(1);
			inputTH1D->SetMarkerColor(ii+1);
			lLeg->AddEntry(inputTH1D,TTemp,"f");
			if(ii==0)
			{
				TF1 *fNoPU = new TF1 ("fNoPU",fSigma3Terms, ptRangeLow_,ptRangeHigh_,3);
				fNoPU->SetParameters(0.03,1.8,0.7);
				fNoPU->SetLineColor(ii+1);
				fNoPU->SetParLimits(0,0,100);
				fNoPU->SetParLimits(1,0,100);
				fNoPU->SetParLimits(2,0,100);
				inputTH1D->Fit(fNoPU,"RB+");
				n = sprintf(buffer,"PU[%d],y=sqrt([%.3f]^2 + [%.3f]^2/x + [%.4f]^2/x^2 ",ii,fNoPU->GetParameter(0),fNoPU->GetParameter(1),fNoPU->GetParameter(2)); // mean1
				TTemp = buffer;
				lLeg->AddEntry(fNoPU,TTemp,"f");
				lLeg2->AddEntry(fNoPU,TTemp,"f");
				n = sprintf(buffer,"Error=%.3f, %.3f, %.4f",fNoPU->GetParError(0),fNoPU->GetParError(1),fNoPU->GetParError(2)); // mean1
				TTemp = buffer;
				lLeg->AddEntry(fNoPU,TTemp,"f");
				noPuResult[0]= fNoPU->GetParameter(0);
				noPuResult[1]= fNoPU->GetParameter(1);
				noPuResult[2]=fNoPU->GetParameter(2);
			}
			if(ii>0)
			{
				TF1 *fPileup = new TF1 ("PileUpIncluded",fSigma4Terms,ptRangeLow_,ptRangeHigh_,4);
				fPileup->SetParameters(noPuResult[0],noPuResult[1],noPuResult[2],noPuResult[2]);
				fPileup->FixParameter(0,noPuResult[0]);
				fPileup->FixParameter(1,noPuResult[1]);
				fPileup->FixParameter(2,noPuResult[2]);
				fPileup->SetParLimits(3,0,100);
				fPileup->SetLineColor(ii+1);
				inputTH1D->Fit(fPileup,"R+");
				n = sprintf(buffer,"PU[%d],y= sqrt([%.3f]^2 + [%.3f]^2/x + [%.4f]^2/x^2 + [%.3f]^2/x^2",ii,fPileup->GetParameter(0),fPileup->GetParameter(1),fPileup->GetParameter(2),fPileup->GetParameter(3)); // mean1
				TTemp = buffer;
				lLeg->AddEntry(fPileup,TTemp,"f");
				lLeg2->AddEntry(fPileup,TTemp,"f");
				n = sprintf(buffer,"Error=%.3f, %.3f, %.4f, %.3f,",fPileup->GetParError(0),fPileup->GetParError(1),fPileup->GetParError(2),fPileup->GetParError(3)); // mean1
				TTemp = buffer;
				lLeg->AddEntry(fPileup,TTemp,"f");
			}
			cCanvas->cd();
			inputTH1D->Draw("P");
			//cCanvas->SetLogy();
			cCanvas->SetLogx();
			lLeg->Draw();
			cCanvas->Update();
			outF_->cd();
			cCanvas->Write();
			cCanvas2->cd();
			if(ii==0)inputTH1D->Draw("P");
			if(ii>0)inputTH1D->Draw("PSame");
			cCanvas2->Update();
			delete cCanvas;
			delete lLeg;
			//delete inputTH1D;
		}
		cout<<"1";
		cCanvas2->cd();
		cCanvas2->SetLogx();
		lLeg2->Draw();
		cCanvas2->Update();
		outF_->cd();
		cCanvas2->Write();
		delete inputTH1D;
		delete cCanvas2;
		delete lLeg2;
		
		cout<<"2";
		TTemp_ = inputFolderNames_[i] + "-pT_Gaus_Sigma_CB";
		combinedPlots(th1ds, OutPutD_, TTemp_, true);
		OutPutD_->mkdir("PUStudies");
		differencePlots (th1ds, (TDirectory*)OutPutD_->Get("PUStudies"), TTemp_, meanPUs);
		th1ds.clear();
		
	}
}

void combinedPlots(vector<TH1D*> th1s, TDirectory *outPutFolderSpecial, TString name, bool rename)
{
	
	TCanvas *cCanvas = new TCanvas(name,name);
	TLegend *lLeg = new TLegend(0.6,0.7,0.9,0.9);
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
			first->SetMarkerSize(2);
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
			first->SetMarkerSize(2);
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

void differencePlots (vector<TH1D*> th1s, TDirectory *outPutFolderSpecial, TString name, vector<double> meanPUs)
{
	TCanvas *cCanvas = new TCanvas(name,name);
	TLegend *lLeg = new TLegend(0.6,0.7,0.9,0.9);
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
		double mean= meanPUs[i+1]-meanPUs[0]; 
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
		if( parameter==4 || parameter==6) TH1DTemp-> Fit(fitts[i],"0R");

		if(i==0)
		{
			TH1DTemp->SetMarkerColor(1);
			TH1DTemp->SetName(parName);
			TH1DTemp->SetMarkerSize(2);
			TH1DTemp->SetLineColor(1);
			TH1DTemp->SetLineWidth(2);
			TH1DTemp->SetMarkerStyle(21);
			TH1DTemp->Draw("P");
			if(parameter==4 || parameter==6) 
			{
				((TF1*)TH1DTemp->GetFunction("Constant"))->SetLineColor(1);
				((TF1*)TH1DTemp->GetFunction("Constant"))->Draw("Same");
			}
		}
		if(i>0)
		{
			TH1DTemp->SetName(parName);
			TH1DTemp->SetLineColor(i+4);
			TH1DTemp->SetLineWidth(2);
			TH1DTemp->SetMarkerColor(i+4);
			TH1DTemp->SetMarkerSize(2);
			TH1DTemp->SetMarkerStyle(21);
			TH1DTemp->Draw("SameP");
			if(parameter==4 || parameter==6) 
			{
				((TF1*)TH1DTemp->GetFunction("Constant"))->SetLineColor(i+4);
				((TF1*)TH1DTemp->GetFunction("Constant"))->Draw("Same");
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
Double_t fSigma4Terms (Double_t *x, Double_t *par) // sqrt([0]^2 + [1]^2/x + [2]^2/x^2 + [3]^2/x^2
{
	Double_t result = sqrt(par[0] *par[0] + par[2]*par[2]/(x[0]*x[0])+ par[1]*par[1]/x[0] + par[3] * par[3]/(x[0] * x[0]));
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
