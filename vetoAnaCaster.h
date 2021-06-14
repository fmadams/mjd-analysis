//
// vetoAnaCaster.h
//
// used by vetoAnaCaster.C
//
//-------------------------------------------------------------------------------------------------------------------------
// histogram operations

// define histograms
TH1F *hrqdc[32];    //raw
TH1F *hcqdc[32];    //cut on thresh, overflow
TH1F *hrun;
TH1F *hMultip0;
TH1F *hMultip1;
TH1F *hMultip2;
TH1F *hMultip3;
TH1F *hMultip4;
TH1F *hMultip5;
TH1D *hiDet;

TH1F *ht1;
#include <iostream>

using namespace std;

void plotQDCs(string savePath)
{
  	bool HistView = true;
  	//if (!(gROOT->IsBatch()) && HistView)
	//{
  		TCanvas *vcan0 = new TCanvas("vcan0","raw veto QDC",0,0,500,1000);
		vcan0->Divide(4,8,0,0);
    // Draw them in order of pannel
		// for (Int_t i=0; i<32; i++)
		// {
		// 	vcan0->cd(i+1);
		// 	hrqdc[i]->Draw();
		// }

    //Draw them in with 4 top centered on Top

    for(Int_t i=0; i<2; i++)  // 2 centerd in row 1
    {
      vcan0->cd(i+2);
      hrqdc[i]->Draw();  // need panel # to add to i
    }
    for(Int_t i=0; i<2; i++)  // 2 centerd in row 2
    {
      vcan0->cd(i+6);
      hrqdc[i]->Draw(); // need panel # to add to i
    }
    for(Int_t i=0; i<16; i++)  // bottom 4 rows
    {
      vcan0->cd(i+17);
      hrqdc[i]->Draw(); // need panel # to add to i
    }
		TCanvas *vcan1 = new TCanvas("vcan1","thresh cut veto QDC",50,0,500,1000);
		gStyle->SetOptStat("e");
		vcan1->UseCurrentStyle();
		vcan1->Divide(4,8,0,0);
    // Draw them in order of pannel
		// for (Int_t i=0; i<32; i++)
		// {
		// 	vcan1->cd(i+1);
		// 	hcqdc[i]->Draw();
		// }
		// std::cout << "Call to save/print qdc" << std::endl;
		// vcan1->Print(savePath.append(".pdf").c_str(),"pdf");
		// vcan0->~TCanvas();
	 //}
}

void plotMultip(string savePath)
{
  	bool HistView = true;
  	//if (!(gROOT->IsBatch()) && HistView)
	//{
  		TCanvas *mcan0 = new TCanvas("mcan0","multiplicity",0,0,800,800);
		//mcan0->Divide(3,1);
		//mcan0->cd(1);
    		gStyle->SetOptStat(0);

    		hMultip0->SetXTitle("Veto Panel Multiplicity");
    		hMultip0->SetTitle("");
    		hMultip0->GetXaxis()->SetTitleOffset(1.2);
    		hMultip0->GetYaxis()->SetTitleOffset(1.5);
    		hMultip0->GetXaxis()->CenterTitle();
    		hMultip0->GetYaxis()->CenterTitle();

   		hMultip0->GetXaxis()->SetRange(2,20);
		//mcan0->cd();
   		//hMultip0->GetXaxis()->SetRange(2,32);
   		//hMultip0->GetXaxis()->SetRangeUser(2,32);
		hMultip0->Draw();

		//mcan0->cd(2);
		//hMultip1->Draw();
		//mcan0->cd(3);
		//hMultip4->Draw();
		//mcan0->cd(4);
		//hMultip4->Draw();
		//mcan0->cd(5);
		//hMultip0->Draw();
		//mcan0->cd(6);
		//hMultip5->Draw();

		// create output file
		std::cout << "Call to save/print multips" << std::endl;
		mcan0->Print(savePath.append(".pdf").c_str(),"pdf");
		//mcan0->Print(savePath + ".pdf","pdf");
		//}
}

void plotTimeHists()
{
  	bool HistView = true;
  	if (!(gROOT->IsBatch()) && HistView)
	{
  		TCanvas *tcan0 = new TCanvas("tcan0","time sequence",100,0,800,800);
		ht1->GetXaxis()->SetLabelSize(0.02);
    		ht1->SetXTitle("Date");
	 	ht1->SetYTitle("Muon Event Count");
	    	ht1->GetXaxis()->SetTitleOffset(1.2);
	    	ht1->GetYaxis()->SetTitleOffset(1.5);
  	 	ht1->GetXaxis()->CenterTitle();
    		ht1->GetYaxis()->CenterTitle();
		ht1->Draw();
	}
}

//-------------------------------------------------------------------------------------------------------------------------
int PanelMap(int qdcChan, int runNum)
{
	// For Dave and Bradley.
	// Using the 32-panel configuration as "standard", map the QDC index to a physical panel location for all of the run ranges.
	// The figure "Veto Panels, View From The Top" in Veto System Change Log is the map I use for ALL configurations.
	// key: "qdc channel"  value: "panel location"

	map<int,int> panels;

	// 32-panel config (default) - began 7/10/15
	if (runNum < 45000000 && runNum >= 3057) {
		map<int,int> tempMap =
		{{0,1},  {1,2},  {2,3},  {3,4}, {4,5}, {5,6},
		 {6,7},  {7,8},  {8,9},  {9,10}, {10,11}, {11,12},
		 {12,13}, {13,14}, {14,15}, {15,16},
		 {16,17}, {17,18}, {18,19}, {19,20},
		 {20,21}, {21,22}, {22,23}, {23,24},
		 {24,25}, {25,26}, {26,27}, {27,28},
		 {28,29}, {29,30}, {30,31}, {31,32}};
		panels = tempMap;
	}
	// 1st prototype config (24 panels)
	else if (runNum >= 45000509 && runNum <= 45004116) {
		map<int,int> tempMap =
		{{0,1},  {1,2},  {2,3},  {3,4}, {4,5}, {5,6},
		 {6,7},  {7,8},  {8,9},  {9,10}, {10,11}, {11,12},
		 {12,21}, {13,22}, {14,15}, {15,16},
		 {16,17}, {17,18}, {18,19}, {19,20},
		 {20,13}, {21,14}, {22,23}, {23,24},
		 {24,-1}, {25,-1}, {26,-1}, {27,-1},
		 {28,-1}, {29,-1}, {30,-1}, {31,-1}};
		 panels = tempMap;
	}
	// 2nd prototype config (24 panels)
	else if (runNum >= 45004117 && runNum <= 45008659) {
		map<int,int> tempMap =
		{{0,1},  {1,2},  {2,3},  {3,4}, {4,5}, {5,6},
		 {6,7},  {7,8},  {8,9},  {9,10}, {10,11}, {11,12},
		 {12,13}, {13,14}, {14,15}, {15,16},
		 {16,17}, {17,18}, {18,19}, {19,20},
		 {20,21}, {21,22}, {22,23}, {23,24},
		 {24,-1}, {25,-1}, {26,-1}, {27,-1},
		 {28,-1}, {29,-1}, {30,-1}, {31,-1}};
		 panels = tempMap;
	}
	// 1st module 1 (P3JDY) config, 6/24/15 - 7/7/15
	else if (runNum > 0 && runNum <=3056){
		map<int,int> tempMap =
		{{0,1},  {1,2},  {2,3},  {3,4}, {4,5}, {5,6},
		 {6,7},  {7,8},  {8,9},  {9,10}, {10,11}, {11,12},
		 {12,13}, {13,14}, {14,15}, {15,16},
		 {16,17}, {17,18}, {18,19}, {19,20},
		 {20,21}, {21,22}, {22,23}, {23,24},
		 {24,-1}, {25,-1}, {26,-1}, {27,-1},
		 {28,-1}, {29,-1}, {30,-1}, {31,-1}};
		 panels = tempMap;
	}
	else {
		cout << "Panel map not known for this run number!\n";
		return -1;
	}
	int panel = -1;
	auto search = panels.find(qdcChan);
	if(search != panels.end()) panel=search->second;
	return panel;
}
//-------------------------------------------------------------------------------------------------------------------------

//This function takes in 2 top, 2 bottom panels and returns the
//index (1-144) of the combination of those panels
int iDetIndex(int t1,int t2,int b1,int b2)
{
	int theIdet = 0;
	int hi = 0;
	int lo = 0;

	if(t1 == 19 || t2 == 19)
		theIdet += 144/2;
	if(t1 == 22 || t2 == 22)
		theIdet += 144/4;
	if(b1 > b2)
	{
		hi = b1;
		lo = b2;
	}
	else
	{
		lo = b1;
		hi = b2;
	}
	theIdet += lo;
	theIdet += ((hi-7)*6);

	return theIdet;
}

//-------------------------------------------------------------------------------------------------------------------------
