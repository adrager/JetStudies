#include "GlobalVars.h"

class TH1StoringClass : public TObject{
	public:
		TH1StoringClass();
		TH1StoringClass(vector<std::pair <double,double> > ,vector<double> );
		TH1StoringClass(vector<std::pair <double,double> > ,vector<double>, TH2D* inputTH2D);
		void InitSetOfTH1D(std::pair <double,double>, double lowPUBinEdge, TH2D* inputTH2D);
		void Fill(double eta, double NPU, double pT,double weight);
		void SaveResultToFile(TDirectory *ResultFolder);
	private:
		void ResetVariables();
		std::vector<std::vector<std::vector< TH1D* > > > resultTH1Ds_;
		vector<std::pair <double,double> > EtaBins_;
		vector<std::pair <double,double> > LowPuBinEdeges_;
		std::vector<std::vector<std::vector< std::pair <double,double> > > > resultTH1DEdges_;

};


TH1StoringClass::TH1StoringClass() : TObject()
{
	std::cout<<"Warning TauTemplateClass initialized without setting variables might lead to weared results!"<<std::endl;
	
}
TH1StoringClass::TH1StoringClass(vector<std::pair <double,double> > EtaBins,vector<double> LowPuBinEdeges) : TObject() // use this constructor if variable pt bining in different eta regions is needed else use next constructor which initializes all vector entries already with the right (same) PT bined TH1Ds
{
	std::cout<<"Warning! Use this constructor if variable pt bining in different eta regions is needed else use next constructor which initializes all vector entries already with the right (same) PT bined TH1Ds!!!!!"<<std::endl;
	std::cout<<"Warning need to initialize each bin now with the function InitSetOfTH1D, else exception will be throughn!!!!!!!"<<std::endl;
	EtaBins_=EtaBins;
	for(unsigned int i=0; i < LowPuBinEdeges.size();i++)
	{
		if((i+1)<LowPuBinEdeges.size() ) LowPuBinEdeges_.push_back(std::make_pair(LowPuBinEdeges[i],LowPuBinEdeges[i+1] ) );
		else LowPuBinEdeges_.push_back(std::make_pair(LowPuBinEdeges[i],1000) );
	}
	std::cout<<"New object TH1StoringClass created"<<std::endl;
	std::cout<<"EtaBins["<<EtaBins.size()<<"] with PuBins["<<LowPuBinEdeges.size()<<"]"<<std::endl;
	for (unsigned int i=0; i< EtaBins_.size();i++)
	{
		std::cout<<"Eta["<<EtaBins_[i].first<<","<<EtaBins_[i].second<<"] with PU: ";
		std::vector<std::vector< TH1D* > > v1;
		resultTH1Ds_.push_back(v1);
		for (unsigned int ii=0; ii < LowPuBinEdeges_.size();ii++)
		{
			std::cout<<"["<<LowPuBinEdeges_[ii].first<<","<<LowPuBinEdeges_[ii].second<<"]";

			std::vector<TH1D*> v2;
			resultTH1Ds_[i].push_back(v2);
		}
		std::cout<<std::endl;
	}
}
TH1StoringClass::TH1StoringClass(vector<std::pair <double,double> > EtaBins,vector<double> LowPuBinEdeges, TH2D* inputTH2D) : TObject()
{
	std::vector<TH1D*> th1dsVector;
	std::vector< std::pair <double,double> > th1dsEdges;
	unsigned int nbinsX = inputTH2D->GetNbinsX();
	Double_t *xbins = new Double_t[nbinsX+1];
	inputTH2D->GetXaxis()->GetLowEdge(xbins);
	xbins[nbinsX]=inputTH2D->GetXaxis()->GetBinLowEdge(nbinsX+1);
	std::cout<<"Used PT bins["<<nbinsX<<"]: ";
	for (unsigned int i=0; i< nbinsX;i++)
	{
		double xlowtemp = inputTH2D->GetXaxis()->GetBinLowEdge(i);
		if (xlowtemp<0) xlowtemp=0;
		int xlow= (int)xlowtemp;
		double xhightemp = inputTH2D->GetXaxis()->GetBinLowEdge(i+1);
		int xhigh= (int)xhightemp;
		TString temp1 ( Form ("%d", xlow));
		TString temp2 ( Form ("%d", xhigh));
		th1dsVector.push_back(new TH1D("MeanPTForBin="+temp1+"_"+temp2,"MeanPTForBin=["+temp1+","+temp2+"]",100,xlow,xhigh));
		th1dsEdges.push_back(std::make_pair(xlowtemp,xhightemp));
		std::cout<<"["<<th1dsEdges[i].first<<","<<th1dsEdges[i].second<<"], ";

	}
	std::cout<<std::endl;
	EtaBins_=EtaBins;
	for(unsigned int i=0; i < LowPuBinEdeges.size();i++)
	{
		if((i+1)<LowPuBinEdeges.size() ) LowPuBinEdeges_.push_back(std::make_pair(LowPuBinEdeges[i],LowPuBinEdeges[i+1] ) );
		else LowPuBinEdeges_.push_back(std::make_pair(LowPuBinEdeges[i],1000) );
	}
	std::cout<<"New object TH1StoringClass created"<<std::endl;
	std::cout<<"EtaBins["<<EtaBins.size()<<"] with PuBins["<<LowPuBinEdeges.size()<<"]"<<std::endl;
	for (unsigned int i=0; i< EtaBins_.size();i++)
	{
		std::cout<<"Eta["<<EtaBins_[i].first<<","<<EtaBins_[i].second<<"] with PU: ";
		std::vector<std::vector< TH1D* > > v1;
		resultTH1Ds_.push_back(v1);
		std::vector<std::vector<std::pair <double, double> > > vv1;
		resultTH1DEdges_.push_back(vv1);
		for (unsigned int ii=0; ii < LowPuBinEdeges_.size();ii++)
		{
			std::cout<<"["<<LowPuBinEdeges_[ii].first<<","<<LowPuBinEdeges_[ii].second<<"]";
			resultTH1Ds_[i].push_back(th1dsVector);
			resultTH1DEdges_[i].push_back(th1dsEdges);
		}
		std::cout<<std::endl;
	}
}

void TH1StoringClass::InitSetOfTH1D(std::pair <double,double> etabin, double lowPUBinEdge, TH2D* inputTH2D)
{
	unsigned int nbinsX = inputTH2D->GetNbinsX();
	Double_t *xbins = new Double_t[nbinsX+1];
	inputTH2D->GetXaxis()->GetLowEdge(xbins);
	xbins[nbinsX]=inputTH2D->GetXaxis()->GetBinLowEdge(nbinsX+1);
	for (unsigned int i=0; i< nbinsX;i++)
	{
		double xlowtemp = inputTH2D->GetXaxis()->GetBinLowEdge(i);
		int xlow= (int)xlowtemp;
		double xhightemp = inputTH2D->GetXaxis()->GetBinLowEdge(i+1);
		int xhigh= (int)xhightemp;
		TString temp1 ( Form ("%d", xlow));
		TString temp2 ( Form ("%d", xhigh));
		TH1D *output1 = new TH1D("MeanPTForBin="+temp1+"_"+temp2,"MeanPTForBin=["+temp1+","+temp2+"]",100,xlow,xhigh);
	}
	
}

void TH1StoringClass::Fill(double eta, double NPU, double pT,double weight)
{
	unsigned int etaBin=10000;
	unsigned int puBin=10000;
	unsigned int ptBin=10000;
	eta=std::abs(eta);
	for (unsigned int i=0; i<EtaBins_.size(); i++) if(eta > EtaBins_[i].first && eta< EtaBins_[i].second) etaBin=i;
	for (unsigned int i=0; i<LowPuBinEdeges_.size(); i++) if(NPU > LowPuBinEdeges_[i].first && NPU< LowPuBinEdeges_[i].second) puBin=i;
	if(etaBin==10000)
	{
		std::cout<<"Error etabin not found!!! for eta: "<<eta<<std::endl;
		eta+=0.0001;
		for (unsigned int i=0; i<EtaBins_.size(); i++) if(eta > EtaBins_[i].first && eta< EtaBins_[i].second) etaBin=i;
	}
	if(puBin==10000)std::cout<<"Error puBin not found!!! for NPU: "<<NPU<<std::endl;
	for (unsigned int i=0; i < resultTH1DEdges_[etaBin][puBin].size();i++) if(pT > resultTH1DEdges_[etaBin][puBin][i].first && pT < resultTH1DEdges_[etaBin][puBin][i].second) ptBin=i;
	if(ptBin==10000)
	{
		//std::cout<<"Error ptBin not found!!! for pt: "<<pT<<std::endl;	
		//std::cout<<"Setting ptBin to max bin!["<<resultTH1DEdges_[etaBin][puBin].size()-1<<"]"<<std::endl;
		ptBin=resultTH1DEdges_[etaBin][puBin].size()-1;
	}
	
//	std::cout<<"EtaBin found["<<etaBin<< "] for eta: "<<eta<<" Binedges ["<<EtaBins_[etaBin].first<<","<<EtaBins_[etaBin].second<<"]"<<std::endl;
//	std::cout<<"PUBin found["<<puBin<< "] for NPU: "<<NPU<<" Binedges ["<<LowPuBinEdeges_[puBin].first<<","<<LowPuBinEdeges_[puBin].second<<"]"<<std::endl;
//	std::cout<<"PTBin found["<<ptBin<< "] for pT: "<<pT<<" Binedges ["<<resultTH1DEdges_[etaBin][puBin][ptBin].first<<","<<resultTH1DEdges_[etaBin][puBin][ptBin].second<<"]"<<std::endl;
	resultTH1Ds_[etaBin][puBin][ptBin]->Fill(pT,weight);

		
}
void TH1StoringClass::SaveResultToFile(TDirectory *ResultFolder)
{
	ResultFolder->cd();
	TDirectory *EtaFolder = new TDirectory();
	TDirectory *PUFolder = new TDirectory();
	TString temp1 ("");
	for (unsigned int i=0; i< EtaBins_.size();i++)
	{
		temp1=Form ("%.1f_%.1f", EtaBins_[i].first, EtaBins_[i].second);
		ResultFolder->mkdir(temp1);
		EtaFolder = (TDirectory*) ResultFolder->Get(temp1);
		for (unsigned int ii=0; ii < LowPuBinEdeges_.size();ii++)
		{
			temp1=Form ("%.0f_%.0f", LowPuBinEdeges_[ii].first, LowPuBinEdeges_[ii].second);
			EtaFolder->mkdir(temp1);
			PUFolder = (TDirectory*) EtaFolder->Get(temp1);
			for (unsigned int iii=0; iii < resultTH1DEdges_[i][ii].size();iii++)
			{
				PUFolder->cd();
				resultTH1Ds_[i][ii][iii]->Write();
			}
		}
	}
}
void TH1StoringClass::ResetVariables()
{
}
