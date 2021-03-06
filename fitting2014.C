#include "GlobalVars.h"
Double_t fConst (Double_t *x, Double_t *par);
void fitting2014()
{
// define constants and names
	debug_=false; // use if you want addional information
	gROOT->SetBatch(true);
	inF_ = TFile::Open("KalibriPlots_PU_Eta28Bins.root","UPDATE"); // input file 
	outF_ = new TFile("Fitting2014.root","RECREATE"); // file to store output in
	outF_->mkdir("FitFunctions"); // stores example plots of the used functions
	std::vector<TString> cbNames;
	cbNames.push_back("CrystalBallFix");
	cbNames.push_back("CrystalBall");
	// ouput text file for error messages
	Fitting2014errorReport = new ofstream("ErrorsFromfitting2014.txt");
	if(Fitting2014errorReport->is_open() ) std::cout<<"FileOpendForErrorOutPut"<<std::endl;
	*Fitting2014errorReport << "----------------- fitting2014 started--------------\n";
	*Fitting2014errorReport << "This file containes error messages and warnings.\n";
	*Fitting2014errorReport << "In particular warnings are saved if a free fit did not converge and the previous fix fit is being stored instead.\n";
	//
	TString truthPileUpFolderAndName ("MCTruthResponseVsNPUTruth/MCTruthResponseVsNPUTruth_Z2star_NPUTruthSpectrum");
	//set parameters
	//fixed parameters: etaDeptendentFixValues_ containes 10 values first 2 being the fourth and sixst parameter to be fixed for the cb obtained from PU and PT inclusive. 5 pairs each for the corresponding eta bin.

	
	lowPuBinEdeges_.push_back(0); lowPuBinEdeges_.push_back(10); lowPuBinEdeges_.push_back(20); lowPuBinEdeges_.push_back(30); lowPuBinEdeges_.push_back(50);
	quantilesLow_=0.0005; quantilesHigh_=0.9999995;
	lowPTConstFit_=30; highPTConstFit_=1200;
	// input folder structure
	/*inputFolderNames_.push_back("MCTruthResolPUEta0"); inputFolderNames_.push_back("MCTruthResolPUEta1"); inputFolderNames_.push_back("MCTruthResolPUEta2"); inputFolderNames_.push_back("MCTruthResolPUEta3"); inputFolderNames_.push_back("MCTruthResolPUEta4"); 
	*/
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
	cbFixValues_.push_back(-10000); cbFixValues_.push_back(-10000); cbFixValues_.push_back(-10000); cbFixValues_.push_back(-10000); cbFixValues_.push_back(-10000); cbFixValues_.push_back(-10000); cbFixValues_.push_back(-10000); // same order as value definiton in cb funciton Normierung N gaus mean, gaus sigma, n1,alpha1,n2,alph2
	
	// create example fit funcitons
	TDirectory *functionF = (TDirectory*)outF_->Get("FitFunctions");
	if(exampleFunctions(functionF))cout<<"Example functions created."<<endl;
	else cout<<"There was a problem with creating the example functions!!!!"<<endl;
	
	// calculated the mean of each pileup bin for normation purpose
	cout<<"truthPileUpFolderAndName:"<<truthPileUpFolderAndName<<endl;
	TH1F *TH1Temp = (TH1F*) inF_->Get(truthPileUpFolderAndName)->Clone();
	TH1D *TruthPUMeans = new TH1D ("TruthPUMeans","TruthPUMeans",lowPuBinEdeges_.size(),0,lowPuBinEdeges_.size());
	for(unsigned int i=0; i < (lowPuBinEdeges_.size()-1);i++)
	{
		TH1Temp->GetXaxis()->SetRange(TH1Temp->GetXaxis()->FindBin(lowPuBinEdeges_[i]),TH1Temp->GetXaxis()->FindBin(lowPuBinEdeges_[i+1]));
		cout<<"Mean for PU bin["<<lowPuBinEdeges_[i]<<","<<lowPuBinEdeges_[i+1]<<"] is"<<TH1Temp->GetMean()<<std::endl;
		truthPuMeans_.push_back(TH1Temp->GetMean());
		TruthPUMeans->SetBinContent(i+1,TH1Temp->GetMean());
	}
	outF_->cd();
	TruthPUMeans->Write();
	TDirectory *EtaFolder = new TDirectory();
	TDirectory *PuInclusiveFolder = new TDirectory();
	TDirectory *PUFolder = new TDirectory();
	// Use sum over all pileups in the differnt eta bins to get an inclusvie cb fit
	std::vector<std::vector<TF1*> > CBPileupInclusiveFitForEachEta; // use this double vector to get fixed parameters for the cb fitting later on
	quantilesLow_=0.00005; quantilesHigh_=0.9999995;
	for (unsigned int i=0; i < inputFolderNames_.size(); i++)
	{
		outF_->mkdir(inputFolderNames_[i]);
		EtaFolder = (TDirectory*) outF_->Get(inputFolderNames_[i]);
		EtaFolder->mkdir("PUInclusiveFits");
		PuInclusiveFolder = (TDirectory*) EtaFolder->Get("PUInclusiveFits");
		PuInclusiveFolder->cd();
		TTemp_ = inputFolderNames_[i] +"_" + inputTH2Names_[0];
		TH2D *inputTH2D = (TH2D*) ((TDirectory*) inF_->Get(inputFolderNames_[i]) )->Get(TTemp_)->Clone();
		inputTH2D->Sumw2();
		for(unsigned int ii=1; ii < inputTH2Names_.size();ii++)
		{
			TTemp_ = inputFolderNames_[i] +"_" + inputTH2Names_[ii];
			TH2D *inputTH2DTemp = (TH2D*) ((TDirectory*) inF_->Get(inputFolderNames_[i]) )->Get(TTemp_)->Clone();
			inputTH2DTemp->Sumw2();
			inputTH2D->Add(inputTH2DTemp,1);
		}
		vector<TH1D*> th1ds = th1creator (inputTH2D, PuInclusiveFolder );
		vector<TF1*> th1fCBPar46EtaDependentConstantFitFix;
		for (unsigned int iii=0; iii< th1ds.size(); iii++)
		{
			th1fCBPar46EtaDependentConstantFitFix.push_back( cbFitting(th1ds[iii], cbFixValues_, true));
			writeTH1D(th1ds[iii], PuInclusiveFolder, cbNames);
			//th1ds[iii]->Write();
			TF1 *constant;
			for(int iv=0; iv <7;iv++)
			{
				cbValues_[iv]->SetBinContent(iii,th1fCBPar46EtaDependentConstantFitFix[iii]->GetParameter(iv));
				cbValues_[iv]->SetBinError(iii,th1fCBPar46EtaDependentConstantFitFix[iii]->GetParError(iv));
					//cout<<"resultCBValues::bin iv="<<iv<<":"<<th1fCBPar46EtaDependentConstantFitFix[iii]->GetParameter(iv)<<"+-"<<th1fCBPar46EtaDependentConstantFitFix[iii]->GetParError(iv)<<endl;
			}
		}
		CBPileupInclusiveFitForEachEta.push_back(th1fCBPar46EtaDependentConstantFitFix);
		
		TF1DPUInclusive_.push_back(th1fCBPar46EtaDependentConstantFitFix);
		PuInclusiveFolder->cd();
		for(int iii=0; iii <7;iii++) 
		{
			if(iii==4 || iii==6)
			{
				if(debug_) cout<<"Used lowPTConstFit_"<<lowPTConstFit_<<", highPTConstFit_"<<highPTConstFit_<<endl;
				TF1 *constant = new TF1("Constant",fConst,lowPTConstFit_,highPTConstFit_,1);
				cbValues_[iii]-> Fit(constant,"0R");
				if(debug_) cout<<"Const Fit done for"<<iii <<" "<<EtaFolder->GetName()<<", with"<<constant->GetParameter(0)<<"+-"<<constant->GetParError(0)<<endl;
				etaDeptendentFixValues_.push_back(constant->GetParameter(0));
				if(debug_) cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			}
		/*	if(iii=4 || iii==6)
			{
				TF1 *constant = new TF1("Constant",fConst,lowPTConstFit_,highPTConstFit_,1);
				constant->SetRange(lowPTConstFit_,highPTConstFit_);
				cbValues_[iii]-> Fit(constant,"0R");
				cout<<"Const Fit done for"<<iii <<" "<<EtaFolder->GetName()<<", with"<<constant->GetParameter(0)<<"+-"<<constant->GetParError(0)<<endl;
		}*/
			TString temp1 ( Form ("%d", iii));
			TTemp_ = "CB_FitValue"+ temp1 +"_"+inputFolderNames_[i] + "_PUInclusiveFits";
			cbValues_[iii]->SetName(TTemp_);
			cbValues_[iii]->SetTitle(TTemp_);
			cbValues_[iii]->Write();
		}
		TTemp_ = "CB_Chi2OverNDF_"+inputFolderNames_[i] + "_PUInclusiveFits";
		cbValues_[7]->SetName(TTemp_);
		cbValues_[7]->SetTitle(TTemp_);
		cbValues_[7]->Write();
		cbValues_.clear();
	}

	
	// separately for each Eta Pu and jet PT
	// temp folders etc 
	quantilesLow_=0.0005; quantilesHigh_=0.9999995;
	for (unsigned int i=0; i < inputFolderNames_.size(); i++)
	{
		cout<<"Starting folder:"<<inputFolderNames_[i]<<"!!!!"<<endl;
		*Fitting2014errorReport <<"inputFolderNames_["<<i<<"]="<<inputFolderNames_[i]<<" opend...";
		// create output folder for each Eta bin
		//outF_->mkdir(inputFolderNames_[i]); laready done in inclusive pu operation see above
		EtaFolder = (TDirectory*) outF_->Get(inputFolderNames_[i]);
		*Fitting2014errorReport<<" EtaFolder->Name()"<<EtaFolder->GetName()<<"\n";
		for(unsigned int ii=0; ii < inputTH2Names_.size();ii++)
		{
			*Fitting2014errorReport<<"inputTH2Names_["<<ii<<"]="<<inputTH2Names_[ii]<<"Opend...";
			cout<<"Starting folder:"<<inputTH2Names_[ii]<<"!!!!"<<endl;
			EtaFolder->mkdir(inputTH2Names_[ii]);
			PUFolder = (TDirectory*) EtaFolder->Get(inputTH2Names_[ii]);
			TTemp_ = inputFolderNames_[i] +"_" + inputTH2Names_[ii];
			CurrentFolder_ = TTemp_;
			TH2D *inputTH2DTemp = (TH2D*) ((TDirectory*) inF_->Get(inputFolderNames_[i]) )->Get(TTemp_)->Clone();
			// do the fitting for the different fixed parameters
			// first for const fitted 4 and 6 parameter for pileup inclusive input 
			PUFolder->cd();
			PUFolder->mkdir("FixToLinearFitPUInclusiveEtaBined46");
			TDirectory *FixToLinearFitPUInclusiveEtaBined46 = (TDirectory*) PUFolder->Get("FixToLinearFitPUInclusiveEtaBined46");
			vector<TH1D*> th1ds = th1creator (inputTH2DTemp, FixToLinearFitPUInclusiveEtaBined46 );
			vector<TF1*> th1fCBPar46EtaDependentConstantFitFix;
			FixToLinearFitPUInclusiveEtaBined46->cd();
			for (unsigned int iii=0; iii< th1ds.size(); iii++)
			{
				quantilesLow_=0.0005; quantilesHigh_=0.9999995;
				cbFixValues_[4]=etaDeptendentFixValues_[i*2];
				cbFixValues_[6]=etaDeptendentFixValues_[i*2+1];
				th1fCBPar46EtaDependentConstantFitFix.push_back( cbFitting(th1ds[iii], cbFixValues_, true) );
				writeTH1D(th1ds[iii], FixToLinearFitPUInclusiveEtaBined46, cbNames);
				th1ds[iii]->Write();
				for(int iv=0; iv <7;iv++)
				{
					cbValues_[iv]->SetBinContent(iii,th1fCBPar46EtaDependentConstantFitFix[iii]->GetParameter(iv));
					cbValues_[iv]->SetBinError(iii,th1fCBPar46EtaDependentConstantFitFix[iii]->GetParError(iv));
					//cout<<"resultCBValues::bin iv="<<iv<<":"<<th1fCBPar46EtaDependentConstantFitFix[iii]->GetParameter(iv)<<"+-"<<th1fCBPar46EtaDependentConstantFitFix[iii]->GetParError(iv)<<endl;
				}
				if(th1fCBPar46EtaDependentConstantFitFix[iii]->GetNDF()>0.1)cbValues_[7]->SetBinContent(iii,th1fCBPar46EtaDependentConstantFitFix[iii]->GetChisquare()/th1fCBPar46EtaDependentConstantFitFix[iii]->GetNDF());
				else cbValues_[7]->SetBinContent(iii,0.00001);
			}
			FixToLinearFitPUInclusiveEtaBined46->cd();
			for(int iii=0; iii <7;iii++) 
			{
				TString temp1 ( Form ("%d", iii));
				TTemp_ = "CB_FitValue"+ temp1 +"_"+inputFolderNames_[i] + "_" + inputTH2Names_[ii];
				cbValues_[iii]->SetName(TTemp_);
				cbValues_[iii]->SetTitle(TTemp_);
				cbValues_[iii]->Write();
			}
			TTemp_ = "CB_Chi2OverNDF_"+inputFolderNames_[i] + "_" + inputTH2Names_[ii];
			cbValues_[7]->SetName(TTemp_);
			cbValues_[7]->SetTitle(TTemp_);
			cbValues_[7]->Write();
			cbValues_.clear();
			
			// do the fitting for the different fixed parameters
			// second for each pt bin for 4 and 6 parameter for pileup inclusive input in eta separately
			PUFolder->cd();
			PUFolder->mkdir("FixToPTPUInclusiveEtaBined46");
			TDirectory *FixToPTPUInclusiveEtaBined46 = (TDirectory*) PUFolder->Get("FixToPTPUInclusiveEtaBined46");
			FixToPTPUInclusiveEtaBined46->cd();
			vector<TH1D*> th1ds2 = th1creator (inputTH2DTemp, FixToPTPUInclusiveEtaBined46 );
			vector<TF1*> th1fCBPar46EtaPTDependent;
			FixToPTPUInclusiveEtaBined46->cd();
			for (unsigned int iii=0; iii< th1ds.size(); iii++)
			{
				quantilesLow_=0.0005; quantilesHigh_=0.9999995;
				cbFixValues_[4]=CBPileupInclusiveFitForEachEta[i][iii]->GetParameter(4);
				cbFixValues_[6]=CBPileupInclusiveFitForEachEta[i][iii]->GetParameter(6);
				th1fCBPar46EtaPTDependent.push_back( cbFitting(th1ds2[iii], cbFixValues_, true) );
				writeTH1D(th1ds2[iii], FixToPTPUInclusiveEtaBined46, cbNames);
				th1ds2[iii]->Write();
				for(int iv=0; iv <7;iv++)
				{
					cbValues_[iv]->SetBinContent(iii,th1fCBPar46EtaPTDependent[iii]->GetParameter(iv));
					cbValues_[iv]->SetBinError(iii,th1fCBPar46EtaPTDependent[iii]->GetParError(iv));
					//cout<<"resultCBValues::bin iv="<<iv<<":"<<th1fCBPar46EtaDependentConstantFitFix[iii]->GetParameter(iv)<<"+-"<<th1fCBPar46EtaDependentConstantFitFix[iii]->GetParError(iv)<<endl;
				}
				if(th1fCBPar46EtaDependentConstantFitFix[iii]->GetNDF()>0.1)cbValues_[7]->SetBinContent(iii,th1fCBPar46EtaDependentConstantFitFix[iii]->GetChisquare()/th1fCBPar46EtaDependentConstantFitFix[iii]->GetNDF());
				else cbValues_[7]->SetBinContent(iii,0.00001);
			}
			FixToPTPUInclusiveEtaBined46->cd();
			for(int iii=0; iii <7;iii++) 
			{
				TString temp1 ( Form ("%d", iii));
				TTemp_ = "CB_FitValue"+ temp1 +"_"+inputFolderNames_[i] + "_" + inputTH2Names_[ii];
				cbValues_[iii]->SetName(TTemp_);
				cbValues_[iii]->SetTitle(TTemp_);
				cbValues_[iii]->Write();
			}
			TTemp_ = "CB_Chi2OverNDF_"+inputFolderNames_[i] + "_" + inputTH2Names_[ii];
			cbValues_[7]->SetName(TTemp_);
			cbValues_[7]->SetTitle(TTemp_);
			cbValues_[7]->Write();
			cbValues_.clear();
			
			// do the fitting for non fixed fitting
			// no fixed parameters
			PUFolder->cd();
			PUFolder->mkdir("NoFixFitting");
			TDirectory *NoFixFitting = (TDirectory*) PUFolder->Get("NoFixFitting");
			NoFixFitting->cd();
			vector<TH1D*> th1ds3 = th1creator (inputTH2DTemp, NoFixFitting );
			vector<TF1*> th1fCBNoFix;
			NoFixFitting->cd();
			for (unsigned int iii=0; iii< th1ds3.size(); iii++)
			{
				quantilesLow_=0.0005; quantilesHigh_=0.9999995;
				cbFixValues_[4]=CBPileupInclusiveFitForEachEta[i][iii]->GetParameter(4);
				cbFixValues_[6]=CBPileupInclusiveFitForEachEta[i][iii]->GetParameter(6);
				th1fCBNoFix.push_back( cbFitting(th1ds3[iii], cbFixValues_, false) );
				writeTH1D(th1ds3[iii], NoFixFitting, cbNames);
				th1ds3[iii]->Write();
				for(int iv=0; iv <7;iv++)
				{
					cbValues_[iv]->SetBinContent(iii,th1fCBNoFix[iii]->GetParameter(iv));
					cbValues_[iv]->SetBinError(iii,th1fCBNoFix[iii]->GetParError(iv));
					//cout<<"resultCBValues::bin iv="<<iv<<":"<<th1fCBPar46EtaDependentConstantFitFix[iii]->GetParameter(iv)<<"+-"<<th1fCBPar46EtaDependentConstantFitFix[iii]->GetParError(iv)<<endl;
				}
				if(th1fCBPar46EtaDependentConstantFitFix[iii]->GetNDF()>0.1)cbValues_[7]->SetBinContent(iii,th1fCBPar46EtaDependentConstantFitFix[iii]->GetChisquare()/th1fCBPar46EtaDependentConstantFitFix[iii]->GetNDF());
				else cbValues_[7]->SetBinContent(iii,0.00001);
			}
			NoFixFitting->cd();
			for(int iii=0; iii <7;iii++) 
			{
				TString temp1 ( Form ("%d", iii));
				TTemp_ = "CB_FitValue"+ temp1 +"_"+inputFolderNames_[i] + "_" + inputTH2Names_[ii];
				cbValues_[iii]->SetName(TTemp_);
				cbValues_[iii]->SetTitle(TTemp_);
				cbValues_[iii]->Write();
			}
			TTemp_ = "CB_Chi2OverNDF_"+inputFolderNames_[i] + "_" + inputTH2Names_[ii];
			cbValues_[7]->SetName(TTemp_);
			cbValues_[7]->SetTitle(TTemp_);
			cbValues_[7]->Write();
			cbValues_.clear();
			
			
			
		}
		*Fitting2014errorReport<<"\n";
	}

}


// methods
bool exampleFunctions(TDirectory *OutPutF)
{
	bool result=true;
	OutPutF->cd();
	TCanvas *temp = new TCanvas("CrystalBall","crystal ball");
	temp->cd();
	TF1 *fCB = new TF1 ("fCB",cb,-4,9,7);
	double cbNorm = 10;
	double cbGaussMean = 2;
	double cbGaussSigma =1;
	double cbNLeftCB = 1;
	double cbAlphaLeft = 1.5; // uebergang zwischen gauss und cb on the left side
	double cbNRightCB = 2;
	double cbAlphaRight = 2;
	fCB->SetParameters(cbNorm,cbGaussMean,cbGaussSigma,cbNLeftCB,cbAlphaLeft,cbNRightCB,cbAlphaRight);
	fCB->SetLineColor(1);
	fCB->SetLineWidth(3);
	fCB->Draw();
	n_ = sprintf(buffer_,"Crystal ball,Norm%.0f,Mean%.0f,Sigma%.0f,NLeftCB%.0f,AlphaLeftCB%.0f,NRightCB%.0f,AlphaRightCB%.0f",cbNorm,cbGaussMean,cbGaussSigma,cbNLeftCB,cbAlphaLeft,cbNRightCB,cbAlphaRight); // mean1
	TTemp_ = buffer_;
	TLegend *legCB = new TLegend(0.6,0.7,0.9,0.9);
	legCB->AddEntry(fCB,TTemp_,"f");
	legCB->SetFillColor(0);
	legCB->Draw();
	temp->Update();
	temp->Write();
	return result;
}
vector<TH1D*> th1creator ( TH2D *inputTH2D, TDirectory *outPutFolder)
{
	vector<TH1D*> result;
	TCanvas *c2 = new TCanvas();
	outPutFolder->cd();
	// create output histogram borders
	unsigned int nbinsX = inputTH2D->GetNbinsX();
	unsigned int nbinsY = inputTH2D->GetNbinsY();
	double lowY = inputTH2D->GetYaxis()->GetXmin();
	double highY = inputTH2D->GetYaxis()->GetXmax();
	Double_t *xbins = new Double_t[nbinsX+1];
	inputTH2D->GetXaxis()->GetLowEdge(xbins);
	xbins[nbinsX]=inputTH2D->GetXaxis()->GetBinLowEdge(nbinsX+1);
	for (unsigned int i=0; i< nbinsX;i++)
	{
		double xlowtemp = inputTH2D->GetXaxis()->GetBinLowEdge(i);
		int xlow= (int)xlowtemp;
		if(xlow<0) xlow=0;
		double xhightemp = inputTH2D->GetXaxis()->GetBinLowEdge(i+1);
		int xhigh= (int)xhightemp;
		if (xlowtemp<0)xlowtemp=0;
		TString temp1 ( Form ("%d", xlow));
		TString temp2 ( Form ("%d", xhigh));
		TString temp3 ();
		if (i+1==nbinsX)
		{
			TString temp3 ("Inf");
		}
		else 
		{
			TString temp3 ( (int) inputTH2D->GetXaxis()->GetBinLowEdge(i+1));
		}
		TH1D *output1 = new TH1D("PT_"+temp1+"_"+temp2,"PT=["+temp1+","+temp2+"]",nbinsY,lowY,highY);
		c2 = new TCanvas("PT_"+temp1+"_"+temp2,"PT=["+temp1+","+temp2+"]",0,0,700,900);
		c2->cd();
		TLegend *leg2 = new TLegend(0.6,0.7,0.9,0.9);
		leg2->SetFillColor(0);
		leg2->AddEntry(output1, output1->GetTitle(),"f");
		for (unsigned int ii=0;ii<nbinsY;ii++)
		{
			output1->SetBinContent(ii,inputTH2D->GetBinContent(i,ii));
			//cout<<"TH1D creator::bin ii="<<ii<<":"<<inputTH2D->GetBinContent(i,ii)<<"+-"<<inputTH2D->GetBinError(i,ii)<<endl;
			output1->SetBinError(ii,inputTH2D->GetBinError(i,ii));
		}
		delete c2;
		result.push_back(output1);
		
	}
	TH1D *cbValues = new TH1D("CBValue0","CBValue0",nbinsX,xbins);
	TH1D *cbValues1 = new TH1D("CBValue1","CBValue1",nbinsX,xbins);
	TH1D *cbValues2 = new TH1D("CBValue2","CBValue2",nbinsX,xbins);
	TH1D *cbValues3 = new TH1D("CBValue3","CBValue3",nbinsX,xbins);
	TH1D *cbValues4 = new TH1D("CBValue4","CBValue4",nbinsX,xbins);
	TH1D *cbValues5 = new TH1D("CBValue5","CBValue5",nbinsX,xbins);
	TH1D *cbValues6 = new TH1D("CBValue6","CBValue6",nbinsX,xbins);
	TH1D *cbValues7 = new TH1D("CBChiSqOverNDF","CBChiSqOverNDF",nbinsX,xbins);
	//TH1D *cbValues8 = new TH1D("CBChiSqOverNDF","CBChiSqOverNDF",nbinsX,xbins);
	cbValues_.push_back(cbValues);
	cbValues_.push_back(cbValues1);
	cbValues_.push_back(cbValues2);
	cbValues_.push_back(cbValues3);
	cbValues_.push_back(cbValues4);
	cbValues_.push_back(cbValues5);
	cbValues_.push_back(cbValues6);
	cbValues_.push_back(cbValues7);
	//cbValues_.push_back(cbValues8);
	return result;
}
TF1* cbFitting(TH1D *th1D, std::vector<double> cbFixValues, bool useFixValues)
{
	//if(debug_)cout<<"cbFitting:: started for th1d name:"<<th1D->GetName()<<endl;
	// calc qunatils for TH1ds
	const Int_t nq = 2;
	Double_t xq[nq];
	xq[0]=quantilesLow_;
	xq[1]=quantilesHigh_;
	Double_t yq[nq];
	if (th1D->ComputeIntegral() <0.0000000001 )
	{
		cout<<"Warning integral of TH1D"<<th1D->GetName()<<", is 0. Skipping th1 for fitting..."<<endl;
		TF1 *CB = new TF1 ("CrystalBallFix",cb,0.5,1.5,7);
		for(int i=0; i<7; i++)CB->FixParameter(i,0.1);
		return CB;
	}
	th1D->GetQuantiles(nq,yq,xq);
	// hard coded set the borders
	//yq[1]=1.6;!!!!!!!!!!!!!!!!!!!!! beginning!!
	TF1 *fFullGauss = new TF1 ("fFullGauss",gauss, 0.2,1.7,3);
	fFullGauss->SetParameters(th1D->GetMaximum(),th1D->GetMean(),th1D->GetRMS());
	th1D->Fit(fFullGauss,"qNr");
	TF1 *fRestrictedGauss = new TF1 ("fRestrictedGauss",gauss, fFullGauss->GetParameter(1)-fFullGauss->GetParameter(2)*0.9, fFullGauss->GetParameter(1)+fFullGauss->GetParameter(2)*0.9,3);
	fRestrictedGauss->SetParameters(th1D->GetMaximum(),fFullGauss->GetParameter(1),th1D->GetRMS()*1.1);
	th1D->Fit(fRestrictedGauss,"NRq+");
	const int Npar = 7;
	double parametersd[Npar] = { (double)th1D->GetMaximum(),(double)fRestrictedGauss->GetParameter(1),(double)fRestrictedGauss->GetParameter(2),3,1.5,3,1.5};
	// some checks to make sure the starting values are reasonable, they could be very small or large if the fRestrictedGauss did not converge 
	if( std::abs((fRestrictedGauss->GetParameter(1)-th1D->GetMean())/th1D->GetMean() )>.7)
	{
		cout<<"!!!!!!!!!!!!!!WARNING::the starting parameter of the mean of the cb fit["<<fRestrictedGauss->GetParameter(1)<<"] is very differnt from the mean of the used th1d["<<th1D->GetMean()<<". starting parameter is set to the mean of the th1d"<<endl;
		parametersd[1]=th1D->GetMean(); 
	}
	if(std::abs( (fRestrictedGauss->GetParameter(2) - th1D->GetRMS())/th1D->GetRMS()) > 0.8 ) 
	{
		cout<<"!!!!!!!!!!!!!!WARNING::the starting parameter of the sigma of the cb fit["<<fRestrictedGauss->GetParameter(2)<<"] is very differnt from the RMS of the used th1d["<<th1D->GetRMS()<<". starting parameter is set to the RMS of the th1d"<<endl;
		parametersd[2]=th1D->GetRMS(); 
	}
	//if(debug_)cout<<"cbFitting::parameters set:[";
	for(int i=0; i<Npar;i++)
	{
		if(parametersd[i]<0.000000000001 )
		{
			cout<<"!!!!!!!!!!!!!!!!!!!!!WARNING STARTING VALUE["<<i<<"]="<<parametersd[i]<<" FOR FIT IS SMALLER OR EQUAL TO 0, setting it to 0.000001!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			parametersd[i]=0.000001;
		}
		//if(debug_)cout<<parametersd[i]<<",";
	}
	for (unsigned int i=0; i< cbFixValues_.size();i++)
	{
		if(cbFixValues_[i]>-9999) {
			
			parametersd[i]=cbFixValues_[i];
		}
	}
	//if(debug_)cout<<"]"<<endl;
	TF1 *fcbFix =  new TF1 ("CrystalBallFix",cb,yq[0],yq[1],7);
	// set up[ the wrappedMultiTF1 etc for fitting
	ROOT::Math::WrappedMultiTF1 wCB(*fcbFix,1);
	ROOT::Fit::DataOptions opt; 
	ROOT::Fit::DataRange rangeTH1d;
	rangeTH1d.SetRange(yq[0],yq[1]);
	ROOT::Fit::BinData dataTH1d(opt,rangeTH1d);
	ROOT::Fit::FillData(dataTH1d, th1D);
	std::cout<<"dataTH1d.size()"<<dataTH1d.Size()<<endl;
	if(dataTH1d.Size() < 10) 
	{
		cout<<"Warning dataTH1d.Size()<10 true: "<<dataTH1d.Size()<<". Skipping th1 for fitting..."<<endl;
		TF1 *CB = new TF1 ("CrystalBallFix",cb,0.5,1.5,7);
		for(int i=0; i<7; i++)CB->FixParameter(i,0.1);
		return CB;
	}
	ROOT::Fit::Chi2Function chi2_CB(dataTH1d, wCB);
	ROOT::Fit::Fitter fitter;
	for (int i=0;i <7;i++)std::cout<<"parametersd["<<i<<"]="<<parametersd[i]<<" ";
	std::cout<<std::endl;
	fitter.Config().SetParamsSettings(7,parametersd);
	std::cout<<"1";
	//fitter2.Config().MinimizerOptions().SetPrintLevel(1);
	fitter.Config().MinimizerOptions().SetMinimizerType("Genetic"); // define which method should be used
	std::cout<<"2";
	fitter.Config().ParSettings(1).Fix();
	//if(debug_)cout<<"Parameter ranges"<<endl;
	for (int i=0; i < Npar;i++)
	{
		fitter.Config().ParSettings(i).SetLimits(parametersd[i]/4,parametersd[i]*2);
		fitter.Config().ParSettings(i).SetStepSize(parametersd[i]*0.5);
	//	if(debug_)cout<<"Parameter:["<<i<<"] limits:"<<parametersd[i]/2<<","<<parametersd[i]*2<<"Step size:"<<parametersd[i]*0.005<<endl;
	}
	fitter.Config().ParSettings(0).SetLimits(parametersd[0]* 0.8,parametersd[0]*1.2);
	fitter.Config().ParSettings(0).SetStepSize(parametersd[0]*0.1);
	fitter.Config().ParSettings(1).SetLimits(parametersd[1] * 0.8,parametersd[1]*1.2);
	fitter.Config().ParSettings(1).SetStepSize(parametersd[1]*0.05);
	fitter.Config().ParSettings(2).SetLimits(parametersd[2] * 0.8,parametersd[2]*1.2);
	fitter.Config().ParSettings(2).SetStepSize(parametersd[2]*0.05);
	fitter.Config().ParSettings(3).SetLimits(parametersd[3]/4,40);
	fitter.Config().ParSettings(3).SetStepSize(parametersd[3]*1);
	std::cout<<"3 ";
	fitter.Config().ParSettings(4).SetLimits(.9, 15);
	fitter.Config().ParSettings(4).SetStepSize(parametersd[2] * 0.05);
	std::cout<<"3.1 ";
	fitter.Config().ParSettings(5).SetLimits(parametersd[5]/4,40);
	fitter.Config().ParSettings(5).SetStepSize(parametersd[5]*1);
	fitter.Config().ParSettings(6).SetLimits(.9, 15);
	std::cout<<"3.2 ";
	fitter.Config().ParSettings(6).SetStepSize(parametersd[2] * 0.05);
	std::cout<<"4";
	for (unsigned int i=0; i< cbFixValues_.size();i++)
	{
		if(cbFixValues_[i]>-9999) {

			fitter.Config().ParSettings(i).Fix();
		}
	}
	std::cout<<"5";
	fitter.FitFCN(chi2_CB,0,dataTH1d.Size(),true);
	ROOT::Fit::FitResult resultprint = fitter.Result();
	std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! beginning!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
	resultprint.Print(std::cout);
	std::cout<<"---------------------------------"<<std::endl;
	fcbFix->SetFitResult( resultprint, 0);
	for (int i=0;i <7;i++)
	{
		parametersd[i]=resultprint.Parameter(i);
		fcbFix->FixParameter(i,resultprint.Parameter(i));
		if(debug_) cout<<"Parameter:"<<i<<":value="<<resultprint.Parameter(i)<<"+-"<<resultprint.ParError(i)<<" ";
		fcbFix->SetParError(i,resultprint.ParError(i));
	}
	if( useFixValues) // set values to fixed parameters if wanted
	{
		for (unsigned int i=0; i< cbFixValues_.size();i++)
		{
			if(cbFixValues_[i]>-9999) 
			{
			
				parametersd[i]=cbFixValues_[i];
			}
		}
		
	}
	if(debug_) std::cout<<std::endl;
	fcbFix->SetNDF(resultprint.Ndf());
	fcbFix->SetChisquare(resultprint.Chi2());
	// second iteration with other algroithm
	ROOT::Math::WrappedMultiTF1 wCB2(*fcbFix,1);
	ROOT::Fit::DataOptions opt2; 
	ROOT::Fit::DataRange rangeTH1d2;
	rangeTH1d2.SetRange(yq[0],yq[1]);
	ROOT::Fit::BinData dataTH1d2(opt2,rangeTH1d2);
	ROOT::Fit::FillData(dataTH1d2, th1D);
	ROOT::Fit::Chi2Function chi2_CB2(dataTH1d2, wCB2);
	ROOT::Fit::Fitter fitter2;
	fitter2.Config().SetMinimizer("GSLMultiFit", "BFGS2");
	fitter2.Config().SetParamsSettings(7,parametersd);
	ROOT::Math::Minimizer* min = fitter2.Config().CreateMinimizer();
	min->SetMaxIterations(10000);
	min->SetTolerance(0.001);
	min->SetValidError(true);
	for(unsigned int i=0; i < 7;i++)
	{
		if(parametersd[i]<0) fitter2.Config().ParSettings(i).SetLimits(0.000001,10);
		if(i>0 && parametersd[i]>10000 ) fitter2.Config().ParSettings(i).SetLimits(parametersd[i] * 0.1,100);
		if(parametersd[i] >0 && parametersd[i]<10000)fitter2.Config().ParSettings(i).SetLimits(parametersd[i] * 0.1,parametersd[i]*1.5);

	}
	if(useFixValues)
	{
		for (unsigned int i=0; i< cbFixValues_.size();i++)
		{
			if(cbFixValues_[i]>-9999) {
				std::cout<<"Vlaue being fixed in second fitting process"<<endl;
				fitter2.Config().ParSettings(i).Fix();
			}
		}
	}
	fitter2.FitFCN(chi2_CB2,parametersd,dataTH1d2.Size(),true);
	ROOT::Fit::FitResult resultprint2 = fitter2.Result();
	resultprint2.Print(std::cout);
	TF1 *fcb2 =  new TF1 ("CrystalBall",cb,yq[0],yq[1],7);
	fcb2->SetFitResult( resultprint2, 0);
	for (int i=0;i <7;i++)
	{
		parametersd[i]=resultprint2.Parameter(i);
		fcb2->FixParameter(i,resultprint2.Parameter(i));
		if(debug_) cout<<"Parameter:"<<i<<":value="<<resultprint2.Parameter(i)<<"+-"<<resultprint2.ParError(i)<<" ";
		fcb2->SetParError(i,resultprint2.ParError(i));
	}
	fcb2->SetNDF(resultprint2.Ndf());
	fcb2->SetChisquare(resultprint2.Chi2());
	// end of fitting
	std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! end of fitting!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
	fcbFix->SetLineColor(3);
	fcb2->SetLineColor(2);
	if( resultprint2.Status()<0)
	{
		ROOT::Math::WrappedMultiTF1 wCB3(*fcbFix,1);
		ROOT::Fit::DataOptions opt3; 
		ROOT::Fit::DataRange rangeTH1d3;
		rangeTH1d3.SetRange(yq[0],yq[1]);
		ROOT::Fit::BinData dataTH1d3(opt3,rangeTH1d3);
		ROOT::Fit::FillData(dataTH1d3, th1D);
		ROOT::Fit::Chi2Function chi2_CB3(dataTH1d3, wCB3);
		ROOT::Fit::Fitter fitter3;
//		fitter3.Config().MinimizerOptions().SetMinimizerType("GSLSimAn");
		fitter3.Config().SetMinimizer("Minuit2", "Fumili2");
		fitter3.Config().SetParamsSettings(7,parametersd);
		ROOT::Math::Minimizer* min3 = fitter3.Config().CreateMinimizer();
		min3->SetMaxIterations(10000);
		min3->SetTolerance(0.001);
		min3->SetValidError(true);
		for(unsigned int i=0; i < 7;i++)
		{

			if(resultprint.Parameter(i)>0 )
			{
				fitter3.Config().ParSettings(i).SetLimits(resultprint.Parameter(i)*0.1,resultprint2.Parameter(i)*10);
				std::cout<<"Limits par["<<i<<"]: "<<resultprint.Parameter(i)*0.1<<" "<<resultprint.Parameter(i)*10<<std::endl;
			}
			if(resultprint.Parameter(i)<0 )
			{
				fitter3.Config().ParSettings(i).SetLimits(resultprint.Parameter(i)*10,resultprint2.Parameter(i)*0.1);
				std::cout<<"Limits par["<<i<<"]: "<<resultprint.Parameter(i)*10<<" "<<resultprint.Parameter(i)*0.1<<std::endl;
			}
		}
		fitter3.Config().ParSettings(3).SetLimits(0.1,40);
		fitter3.Config().ParSettings(5).SetLimits(0.1,40);
		for(unsigned int i=0; i < 7;i++) std::cout<<"fitter3 limits par["<<i<<"] is bound: "<<fitter3.Config().ParSettings(i).IsBound()<<", is fixed: "<<fitter3.Config().ParSettings(i).IsFixed()<<" lowerLimit: "<<fitter3.Config().ParSettings(i).LowerLimit()<<" upperLimit: "<<fitter3.Config().ParSettings(i).UpperLimit()<<std::endl;
		fitter3.FitFCN(chi2_CB3,parametersd,dataTH1d3.Size(),true);
		ROOT::Fit::FitResult resultprint3 = fitter3.Result();
		resultprint3.Print(std::cout);
		fcb2->SetFitResult( resultprint3, 0);
		for (int i=0;i <7;i++)
		{
			parametersd[i]=resultprint3.Parameter(i);
			fcb2->FixParameter(i,resultprint3.Parameter(i));
			if(debug_) cout<<"Parameter:"<<i<<":value="<<resultprint3.Parameter(i)<<"+-"<<resultprint3.ParError(i)<<" ";
			fcb2->SetParError(i,resultprint3.ParError(i));
		}
		fcb2->SetNDF(resultprint3.Ndf());
		fcb2->SetChisquare(resultprint3.Chi2());
		
		
		*Fitting2014errorReport <<"Warning free fit did not converge for:\n";
		*Fitting2014errorReport <<"CurrentFolder:"<<CurrentFolder_<<"\n";
		*Fitting2014errorReport <<"th1D->GetName():"<<th1D->GetName()<<"\n";
		*Fitting2014errorReport <<"Using now Minut2 with Migrad: (status:"<<resultprint2.Status()<<"):\n";
		for(int i=0; i<7; i++)*Fitting2014errorReport <<"Par["<<i<<"]: "<<resultprint2.Parameter(i)<<" +- "<<resultprint2.ParError(i)<<"\n";
		*Fitting2014errorReport <<"Genetic(status:"<<resultprint.Status()<<")\n";
		for(int i=0; i<7; i++)*Fitting2014errorReport <<"Par["<<i<<"]: "<<resultprint.Parameter(i)<<" +- "<<resultprint.ParError(i)<<"\n";
		*Fitting2014errorReport <<"Minuit2(status:"<<resultprint3.Status()<<")\n";
		for(int i=0; i<7; i++)*Fitting2014errorReport <<"Par["<<i<<"]: "<<resultprint3.Parameter(i)<<" +- "<<resultprint3.ParError(i)<<"\n";
		for(unsigned int i=0; i < 7;i++) *Fitting2014errorReport<<"fitter3 limits par["<<i<<"] is bound: "<<fitter3.Config().ParSettings(i).IsBound()<<", is fixed: "<<fitter3.Config().ParSettings(i).IsFixed()<<" lowerLimit: "<<fitter3.Config().ParSettings(i).LowerLimit()<<" upperLimit: "<<fitter3.Config().ParSettings(i).UpperLimit()<<"\n";
		*Fitting2014errorReport <<"Updated!\n"<<std::flush;
	}
	th1D->GetListOfFunctions()->Add(fcbFix);
	th1D->GetListOfFunctions()->Add(fcb2);
	std::cout<<"Resulting parameters: ";
	for(int i=0; i<7;i++)
	{
		std::cout<<i<<":"<<fcb2->GetParameter(i)<<"+-"<<fcb2->GetParError(i)<<" ";
	}
	std::cout<<std::endl;
	std::cout<<"Should be: ";
	for(int i=0; i<7;i++)
	{
		std::cout<<i<<":"<<resultprint2.Parameter(i)<<"+-"<<resultprint2.ParError(i)<<" ";
	}
	std::cout<<std::endl;
	return fcb2;
	
}

void writeTH1D(TH1D* inputTH1D, TDirectory *outPutFolder, std::vector<TString> FunctionNames)
{
	char buffer[50];
	int n;
	TString TTemp("");
	TCanvas *cCanvas = new TCanvas(inputTH1D->GetName(),inputTH1D->GetTitle());
	TLegend *lLeg = new TLegend(0.6,0.7,0.9,0.9);
	lLeg->SetFillColor(0);
	cCanvas->cd();
	inputTH1D->Draw("P");
	lLeg->AddEntry(inputTH1D,"Entries","f");
	TF1 *fcb;
	for (unsigned int i=0; i<FunctionNames.size();i++)
	{
		fcb= (TF1*)inputTH1D->GetListOfFunctions()->FindObject(FunctionNames[i]);
		if( !fcb) 
		{
			std::cout<<"Error Function name :"<<FunctionNames[i]<<", not found!"<<std::endl;
			continue;
		}
		TTemp= FunctionNames[i];
		double chi2Temp=fcb->GetChisquare()/fcb->GetNDF();
		n = sprintf(buffer,"CB:Weighted Chi2=%.2f", chi2Temp); // mean1
		TTemp = TTemp + "_" + buffer;
		lLeg->AddEntry(fcb,TTemp,"f");	
	}
	outPutFolder->cd();
	lLeg->Draw();
	cCanvas->Update();
	cCanvas->Write();
	delete cCanvas;
	delete lLeg;
}

// functions
Double_t cb(Double_t *x, Double_t *par)
{
		/*
	par[0] Normierung N of total function
	par[1]  gaus mean
	par[2]  gaus sigma
	par[3]  n cb 1  free parameter of the potenz law
	par[4] alpha units of sigma at which the gauss turns into the potenz law
	par[5]  n cb 2 free parameter of the potenz law
	par[6] alpha proportional to the cutoff at which cb2 takes over
		*/
	Double_t N = par[0];
	Double_t p  = (x[0]-par[1])/par[2];
	Double_t n1 = par[3];
	Double_t alpha1 = par[4];
	Double_t n2 = par[5];
	Double_t alpha2 = par[6];
	Double_t fitval=0;
	Double_t gauss = N * TMath::Exp(-0.5*p*p);
	Double_t cb1 = N * TMath::Power(n1/std::abs(alpha1),n1) * TMath::Exp(-0.5*alpha1*alpha1) * TMath::Power( (n1/std::abs(alpha1) - alpha1 - p),-n1);
	Double_t cb2 = N * TMath::Power(n2/std::abs(alpha2),n2) * TMath::Exp(-0.5*alpha2*alpha2) * TMath::Power( (n2/std::abs(alpha2) - alpha2 + p),-n2);
	if ( p <= -alpha1)fitval=cb1;
	else if(p >= alpha2) fitval=cb2;
	else if(p> -alpha1 && p <alpha2) fitval=gauss;
	if(fitval!=fitval)
	{
		std::cout<<"par0:"<<par[0]<<", par1:"<<par[1]<<", par2:"<<par[2]<<", par3:"<<par[3]<<", par4:"<<par[4]<<", par5:"<<par[5]<<", par6:"<<par[6]<<", x:"<<x[0]<<", gauss"<<gauss<<", cb1"<<cb1<<", cb2"<<cb2<<std::endl;
	}
	assert(fitval==fitval);
	return fitval;
}

Double_t gauss (Double_t *x, Double_t *par)
{
	Double_t gaus = par[0]*TMath::Exp(-0.5*((x[0]-par[1])/par[2]) * ((x[0]-par[1])/par[2]));
	return gaus;	
}

Double_t fConst (Double_t *x, Double_t *par)
{
	Double_t para = par[0];
	return para;
}
