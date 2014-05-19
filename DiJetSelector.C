#define DiJetSelector_cxx
// The class definition in DiJetSelector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("DiJetSelector.C")
// Root > T->Process("DiJetSelector.C","some options")
// Root > T->Process("DiJetSelector.C+")
//

#include "DiJetSelector.h"
#include <TH2.h>
#include <TStyle.h>


void DiJetSelector::Begin(TTree * /*tree*/)
{
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
	inputFolderNames_.push_back("MCTruthResolPUEta25"); 	etaBins_.push_back(std::make_pair(2.5,3.0) );
	inputFolderNames_.push_back("MCTruthResolPUEta26"); 	etaBins_.push_back(std::make_pair(3.0,5.0) );
	inputTH2Names_.push_back("GenJetResponseVsGenJetPt_Z2star_L2L3_NPU0"); inputTH2Names_.push_back("GenJetResponseVsGenJetPt_Z2star_L2L3_NPU1"); inputTH2Names_.push_back("GenJetResponseVsGenJetPt_Z2star_L2L3_NPU2"); inputTH2Names_.push_back("GenJetResponseVsGenJetPt_Z2star_L2L3_NPU3"); 
	lowPuBinEdeges_.push_back(0); lowPuBinEdeges_.push_back(10); lowPuBinEdeges_.push_back(20); lowPuBinEdeges_.push_back(30); lowPuBinEdeges_.push_back(50);
   std::cout<<"Starting DiJetSelector to determine the composition of events in each PT bin"<<std::endl;
   if(etaBins_.size() != inputFolderNames_.size()) std::cout<<"Error amount of inputFolderNames_ (should be for each eta bin one folder) is not the same as the amount of input etabin pairs etaBins_!"<<std::endl;
   TDirectory *EtaFolder = new TDirectory();
   TDirectory *PuInclusiveFolder = new TDirectory();
   TDirectory *PUFolder = new TDirectory();
   inF_ = TFile::Open("KalibriPlots_PU_Eta26Bins.root","UPDATE");
   
   TTemp_ = inputFolderNames_[0] +"_" + inputTH2Names_[0];
   std::cout<<"Opening TH2D: inF_->Get("<<inputFolderNames_[0]<<") )->Get("<<TTemp_<<")->Clone()"<<std::endl;
   TH2D *inputTH2Dtemp = (TH2D*) ((TDirectory*) inF_->Get(inputFolderNames_[0]) )->Get(TTemp_)->Clone();
   th1Results_ = new TH1StoringClass(etaBins_,lowPuBinEdeges_,inputTH2Dtemp);
   th1Results_->Fill(0.01,25,60,1);
   th1Results_->Fill(0.14,60,200,1);
   th1Results_->Fill(1.14,60,200,1);
   th1Results_->Fill(1.000000001,60,200,1);

   for (unsigned int i=0; i < inputFolderNames_.size(); i++) // loop over all eta bin folders
   {
	   EtaFolder = (TDirectory*) inF_->Get(inputFolderNames_[i]);
	   for(unsigned int ii=0; ii < inputTH2Names_.size();ii++) // loop over all PU bins
	   {
		   PUFolder = (TDirectory*) EtaFolder->Get(inputTH2Names_[ii]);
		   TTemp_ = inputFolderNames_[i] +"_" + inputTH2Names_[ii];
		   TH2D *inputTH2D = (TH2D*) ((TDirectory*) inF_->Get(inputFolderNames_[i]) )->Get(TTemp_)->Clone();
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
			   TH1D *output1 = new TH1D("MeanPTForBin="+temp1+"_"+temp2,"MeanPTForBin=["+temp1+","+temp2+"]",100,xlow,xhigh);
			   MeanPTTH1D_.push_back(output1);

		   }
	   }
	   
   }

   std::map<unsigned int,std::map<unsigned int,TH1D*> > mapPUOfEta_;
   for (unsigned int i=0; i<inputTH2Names_.size(); i++)
   {
	   std::map<unsigned int,TH1D*> mapOfEta;
	   for (unsigned int ii=0; ii < inputFolderNames_.size(); ii++)
	   {
		   int xlow= (int)1;
		   int xhigh= (int)2;
		   TString temp1 ( Form ("%d", xlow));
		   TString temp2 ( Form ("%d", xhigh));
		   mapOfEta[ii]= new TH1D("MeanPTForBin="+temp1+"_"+temp2,"MeanPTForBin=["+temp1+","+temp2+"]",100,xlow,xhigh);
	   }
	   mapPUOfEta_[i] = mapOfEta;
   }

   TString option = GetOption();
   std::cout<<"Done: starting prediction"<<std::endl;
count=0;
}

void DiJetSelector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t DiJetSelector::Process(Long64_t entry)
{
   	resetValues();
	fChain->GetTree()->GetEntry(entry);
  	count++;
	if(count==1000) {std::cout<<"|";count=0;}
	if(NobjGenJet<3)return kTRUE;
	double NPU = PUMCNumVtx+0.001;
	for (unsigned in=0; in<2; in++)
	{
		double genJetPt=GenJetPt[in];
		double genJetEta = GenJetPhi[in];
   		th1Results_->Fill(genJetEta,NPU,genJetPt,1);

		
		
	}

   return kTRUE;
}

void DiJetSelector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void DiJetSelector::Terminate()
{
std::cout<<"DiJetSelector done startintg termiate process..."<<std::endl;
   TFile *outF = new TFile("PTBinsMeans.root","RECREATE");
   outF->mkdir("PTBinmeans");
   TDirectory *ptBinOut = (TDirectory*) outF->Get("PTBinmeans");
   th1Results_->SaveResultToFile(ptBinOut);
   std::cout<<"Done exiting DiJetSelector"<<std::endl;

}
void DiJetSelector::resetValues()
{
}
