//
//
// A simple TTreeReader adapted by DJT for MJD
//
// run as:
// root vetoAnaCasterNew.C++
//
// ------------------------------------------------------------------
//
//Modifying Author: Benjamin Ranson
//Spring and Summer of 2019
//Modified to detect certain patterns in panel firings, as well as act as a detector for multiplicity firings above a certain threshold
//Has the ability to "cast" onto multiple data files and use that information to output for each file
//
//Modifying Author: Franklin Adams
//Summer of 2021
//

#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include <TCanvas.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include "TBenchmark.h"
#include "TStyle.h"

#include "MJVetoEvent.hh"

#include <vector>
#include <iostream>
#include <ctime>
#include <utility>

#include "vetoAnaCasterNew.h"

#define HIGH_MULTIP_OUTPUT_LIST_NAME "high-multip-list.txt"
#define HIGH_MULTIP_THRESHOLD 16
#define MULTIP_TABLE_OUTPUT_NAME "multips-table.txt"
#define SKIM_CUT_MODIFIER "-sc"

using namespace std;

struct SetData
{
  string name;
  int fourPanelEvents;
  double totalTime;
  double maxRunDuration; //Longest run of the set
  vector<int> zeroDurationRuns;
  vector<int> zeroDurationFourPanelRuns;

  double getFourPanelRate() {return (double)(this->fourPanelEvents / this->totalTime);}
  double getFourPanelRateErr() {return (double)sqrt(this->fourPanelEvents)/(this->totalTime);}
};

struct RunStats
{
  double runDuration;
  int runNumber;
  bool zeroScalerDuration; //If this should be a zero-duration event based on data
  int numEvents; //ANY event
  int fourPanelEvents; //4 panel events

  bool operator< (const RunStats& other) const {return runDuration < other.runDuration;}
};

struct RunSet
{
  string baseName; //Basic name, like P3LTP, P3JDY, etc
  string extName; //Extended name, like P3LTPNz
  string path; //some/thing/to/where/file.root

  RunSet(string baseName, string extName, string path)
  {
    this->baseName = baseName; this->extName = extName; this->path = path;
  }
};

SetData ana(RunSet runSet);

const string rootFileFolder = "/Users/franklinadams/muon_cross_section/veto-skim/"; //Folder where many .root files can be found
const string mDataFolder = "/Users/franklinadams/muon-reanalysis/muon-data";
const string QDC_AGGLOM_FILE_NAME = "qdc-agglom";
const string RUN_TIMING_FILE_NAME = "run-timing.pdf";
const string RUN_COMPARISON_FILE_NAME = "run-comparison.pdf";
const string ZERO_RUN_OUTPUT_FILE = "zero-runs.txt";
//const string SIM_COMP_OUTPUT_FILE = "sim-comp.root";

const bool DO_HI_MULTIP_CUT = false; //Set to true if you want the program to cut high multiplicites into a high-multip-list
const bool DO_MULTIP_TABLE = false; //Set to true if you want the program to output "run-number multip" into a table text file
const bool DO_QDC_AGGLOM = false; //Set to true if you want the program to agglomerate all QDC data from each set into a saved .root file
const bool DO_RUN_TIMING = true; //Set to true to graph run duration vs run number
const bool DO_ZERO_RUN = false; //Set to true if you want to output the run numbers of 0 duration runs

void vetoAnaCasterNew()
{
        //Clear the QDD agglom file
        TFile agglomFile((QDC_AGGLOM_FILE_NAME + ".root").c_str(), "RECREATE"); //Creates the file or clears the existing one
	agglomFile.Close();

        //Begin clearing the multip table
        ofstream fileClearer;
        fileClearer.open(MULTIP_TABLE_OUTPUT_NAME, ios::trunc);
        fileClearer.close(); //Clear the multipTable output file of all contents


	const RunSet targets[] = //All the files to execute this script on, including their directory and indentifier
	{
	 // RunSet("P3JDY", "P3JDYNz", mDataFolder + "/P3JDY/skimVeto_P3JDYNz-skim-cut.root"),
	 // RunSet("P3KJR", "P3JKJRNz", mDataFolder + "/P3KJR/skimVeto_P3KJRNz-skim-cut.root"),
	 // RunSet("P3LQKa", "P3LQKaNz", mDataFolder + "/P3LQKa/skimVeto_P3LQKaNz-skim-cut.root"),
	 // RunSet("P3LQK2", "P3LQK2Nz", mDataFolder + "/P3LQK2/skimVeto_P3LQK2Nz-skim-cut.root"),
	 // RunSet("P3LTP", "P3LTPNz", mDataFolder + "/P3LTP/skimVeto_P3LTPNz-skim-cut.root"),
	 // RunSet("P3LTP2", "P3LTP2Nz", mDataFolder + "/P3LTP2/skimVeto_P3LTP2Nz-skim-cut.root"),
	 // RunSet("P3LTP3", "P3LTP3Nz", mDataFolder + "/P3LTP3/skimVeto_P3LTP3Nz-skim-cut.root"),
   // New Runs
   // RunSet("P3LTP4", "P3LTP4Nz", mDataFolder + "/P3LTP4/skimVeto_P3LTP4.root"),
   // RunSet("P3LTP5", "P3LTP5Nz", mDataFolder + "/P3LTP5/skimVeto_P3LTP5.root"),
   // RunSet("P3N991", "P3N991Nz", mDataFolder + "/P3N991/skimVeto_P3N991.root"),
   RunSet("P3N992", "P3N992Nz", mDataFolder + "/P3N992/skimVeto_P3N992.root"),
	};

	int numTargets = sizeof(targets) / sizeof(targets[0]);

	//Aggregate data about the sets
	SetData allData[numTargets];

	for (int i = 0; i < numTargets; i++)
	{
	  cout << "Casting onto the data file " << targets[i].extName << endl;
	  cout << "\tlocated at " << targets[i].path << endl;
	  allData[i] = ana(targets[i]);
	  cout << "allData[" << i << "].fourPanelEvents = " << allData[i].fourPanelEvents << endl;
	}

	if (DO_RUN_TIMING)
	{
	  TCanvas* comparisonCanvas = new TCanvas("comparisonCanvas", "Comparison of Muons/time", 500, 500);
	  //comparisonCanvas->SetFillColor(42);
	  //TGraphErrors* runComparison = new TGraphErrors(numTargets, x, y, xe, ye);
	  TH1D* runComparison = new TH1D("runComparison", "Comparison of Muons/time", numTargets, 0, numTargets);
	  runComparison->GetYaxis()->SetTitle("Counts per second");
	  runComparison->GetXaxis()->SetTitle("Configuration");
	  for (int i = 0; i < numTargets; i++)
	  {
	    runComparison->SetBinContent(i + 1, allData[i].getFourPanelRate());
	    cout << "Bin " << i + 1 << " value: " << allData[i].getFourPanelRate() << endl;
	    cout << "Bin " << i + 1 << " 4P#: " << allData[i].fourPanelEvents << endl;
	    runComparison->SetBinError(i + 1, allData[i].getFourPanelRateErr());
	    cout << "Bin " << i + 1 << " error: " << allData[i].getFourPanelRateErr() << endl;
	    runComparison->GetXaxis()->SetBinLabel(i + 1, allData[i].name.c_str());
	  }
	  runComparison->Draw("E");
	  comparisonCanvas->Print(RUN_COMPARISON_FILE_NAME.c_str(), "pdf");
	  delete comparisonCanvas;
	  delete runComparison;
	}

	//Exporting the QDC agglom to be identical to any other QDC
	if (DO_QDC_AGGLOM)
	{
	  TFile file((QDC_AGGLOM_FILE_NAME + ".root").c_str(), "READ"); //Read all the qdc agglom data
	  for (int j = 0; j < 32; j++)
	  {
	    hcqdc[j] = (TH1F*)file.Get(("hcqdc" + to_string(j)).c_str()); //Let the agglomerated qdc be copied into the hcqdc
	  }
	  plotQDCs((mDataFolder + "/" + QDC_AGGLOM_FILE_NAME).c_str()); //Plot the agglom data just like a regular qdc
	  file.Close();
	}

	if (DO_ZERO_RUN)
	{
	  ofstream zeroRunOutput(ZERO_RUN_OUTPUT_FILE.c_str(), ios::out | ios::trunc);
	  for (int i = 0; i < numTargets; i++)
	  {
	    cout << "In " << targets[i].extName << ": zero duration runs with four panel events:" << allData[i].zeroDurationFourPanelRuns.size() << endl;
	    cout << "\tzero duration runs (of any time): " << allData[i].zeroDurationRuns.size() << endl;
	    for (int a = 0; a < allData[i].zeroDurationRuns.size(); a++)
	    {
	      zeroRunOutput << allData[i].zeroDurationRuns[a] << endl;
	    }
	  }
	  zeroRunOutput.close();
	}
}

SetData ana(RunSet runSet) {

        SetData setData;
	setData.name = runSet.extName;

        cout << "Start of ana() on " << runSet.extName << endl;
        vector<pair<int, int>> multipTable; //First valule is run #, second value is multiplicity

	//Set up collection for run duration vs run number
	vector<RunStats> runStats = vector<RunStats>();

	// set up execution timer
	gBenchmark->Reset();
	gBenchmark->Start("VetoAna");

	// Create histograms
	static Int_t nqdc_bins=100;
	static Float_t ll_qdc=0.;
	static Float_t ul_qdc=4200.;
	Char_t hname[50];
	for (Int_t i=0; i<32; i++)
	{
		sprintf(hname,"hrqdc%d",i);
		hrqdc[i]=new TH1F(hname,hname,nqdc_bins,ll_qdc,ul_qdc);
		hrqdc[i]->SetMarkerColor(3);
		sprintf(hname,"hcqdc%d",i);
		hcqdc[i]=new TH1F(hname,hname,nqdc_bins,ll_qdc,ul_qdc);
		hcqdc[i]->Sumw2();
	}

   	hrun = new TH1F("hrun", "run numbers", 1000, 0, 30000);
   	hMultip0 = new TH1F("hMultip0", "Multiplicity", 33, -0.5, 32.5);
   	hMultip1 = new TH1F("hMultip1", "Multiplicity", 33, -0.5, 32.5);
   	hMultip2 = new TH1F("hMultip2", "Multiplicity", 33, -0.5, 32.5);
   	hMultip3 = new TH1F("hMultip3", "Multiplicity", 33, -0.5, 32.5);
   	hMultip4 = new TH1F("hMultip4", "Multiplicity", 33, -0.5, 32.5);
   	hMultip5 = new TH1F("hMultip5", "Multiplicity", 33, -0.5, 32.5);

	hiDet = new TH1D("hiDet","hiDet",145,0.,144);

	// counts vs time hist
	// limits are in seconds relative to 1/1/95
	// one hour bins = 16667
	// one day bins = 695
	ht1 = new TH1F("ht1","ht1",695,640000000,750000000);
	ht1->GetXaxis()->SetTimeDisplay(1);

   	// Open the file containing the tree.
	//RC: Changed changed the path to work from my own directory, but it is using the main file
	//TODO: Ask if I should copy the main data file OR if there's a read-only way to access it
   	//TFile *myFile = TFile::Open("/Users/Shared/muon_cross_section/veto-skim/skimVeto_DS5.root");
   	TFile *myFile = TFile::Open(runSet.path.c_str());
   	//TFile *myFile = TFile::Open("/Users/Shared/muon-reanalysis/muon-data/P3JDY/skimVeto_P3JDY-skim-cut.root");
	//^ Change input file here

	// The simulation-comparison file helps us compare a run with the simulation data
	// A number of histograms will be stored in this file in order to be later graphable
	//const string simCompFilePath = mDataFolder + "/" + runSet.baseName + "/" + runSet.extName + "-" + SIM_COMP_OUTPUT_FILE;
	//TFile *simCompFile = new TFile(simCompFilePath.c_str(), "RECREATE");

   	// Create a TTreeReader for the tree, for instance by passing the TTree's name and the TDirectory / TFile it is in.
	TTreeReader reader("vetoTree",myFile);

   	// Address some branches
		TTreeReaderValue<Int_t> run(reader, "run");
	TTreeReaderValue<MJVetoEvent> vetoEvent(reader,"vetoEvent");
   		TTreeReaderValue<Int_t> fCard1(reader, "fCard1");
   		TTreeReaderValue<Int_t>  fCard2(reader, "fCard2");
   		TTreeReaderValue<Int_t>  fRun(reader, "fRun");
   		TTreeReaderValue<Int_t>  fEntry(reader, "fEntry");
		TTreeReaderValue<Bool_t> fBadScaler(reader,"fBadScaler");
		TTreeReaderArray<int> fQDC(reader,"fQDC[32]");
		TTreeReaderArray<int> fSWThresh(reader,"fSWThresh[32]");
   		TTreeReaderValue<Int_t>  fMultip(reader, "fMultip");

	TTreeReaderValue<double> xTime(reader,"xTime");
	TTreeReaderValue<Long64_t> start(reader,"start");
	TTreeReaderValue<Long64_t> stop(reader,"stop");
	TTreeReaderValue<double> unixDuration(reader,"unixDuration");
	TTreeReaderValue<double> scalerDuration(reader,"scalerDuration");

	//CoinType[0]: "weak" muon candidate.  2 or more panels fired over 500QDC
	//CoinType[1]: "strong" muon candidate.  Require 2 top + 2 bottom coincidence
	//CoinType[2]: both planes of a side + both bottom planes (edited)
	//CoinType[3]: both top planes + both planes of a side
	TTreeReaderArray<int> CoinType(reader,"CoinType");	//[32]


	//some useful vars
	Int_t iDet=0;
	Int_t FourPanelHits[145]={0};  //why 145?
	Int_t FourPanelHitsTOT=0.;

	Double_t total_run_time=0;
	Int_t last_run = 0;
	Int_t total_runs = 0;

	//ctime structures
	tm *last_run_tm;
	tm *run_tm;
	int last_run_year=0;
	int last_run_mon=0;
	int last_run_mday=0;

	Int_t iday=-1;      // get incremented to zero on first day
	Int_t ndays=0;
	Double_t daily_Duration[1000];
	Double_t daily_count[1000];
	Double_t daily_day[1000];
	Double_t daily_mon[1000];
	Double_t daily_year[1000];


	//initialize arrays
	for (Int_t i=0;i<1000;i++){
		daily_Duration[i]=0.;
		daily_count[i]=0;
	}

	// Loop over all entries of the TTree
	// loop over panels for the event, save 4 hit info
	// save run time (use scalerDuration) - clint says perhaps use unixDuration...

	int totalThreePanelEvents = 0; //RC: Increases with each hit
       	int classBCount = 0; //RC
	int classCCount = 0; //RC
	int classDCount = 0; //RC
	int totalCutCount = 0; //RC
	Int_t firstRun = -1; //RC
	//RC:
	//Class B = Any event
	//Class C = Coin Type 1 events
	//Class D = Coin Type 0 events

	vector<int> highMultipRuns; //Contains all run numbers considered to be high multiplicity

	//const std::string DISPLAY_INPUT_PATH = "/Users/ranson/Documents/Veto_Display_Input/3_Panel_Special_Events_2.txt";
	//ofstream writer(DISPLAY_INPUT_PATH.c_str(), ios::out | ios::trunc); //Used to output collected data

	Int_t finalRun = 0; //Holds the final run number looped over
	tm finalTimeValue; //Holds the last run time value for the data set
	tm firstTimeValue; //Holds the first run time value for the data set

	while (reader.Next()) {
		//cout << *fEntry << " " << *fRun << " " << *fCard1 << " " << *fCard2 << endl;
		//cout << vetoEvent->fQDC[32]->size() << endl;
		//if (CoinType[1] || CoinType[2] || CoinType[3] || CoinType[4]){
		//	cout << CoinType[1] << " "<< CoinType[2] << " "<< CoinType[3] << " "<< CoinType[4] << endl;
		//	cout << *run << " "<< *start << " " << *stop << " " << *unixDuration << " " << *scalerDuration << endl;
		//}


		// keep track of the time
		if (*run != last_run){

		        if (DO_ZERO_RUN && *scalerDuration == 0) // minimum run time to consider
			{
			  setData.zeroDurationRuns.push_back(*run);
			  //cout << "Event inside a run with duration 0! Run #" << *run << endl;
			}

		        if (DO_RUN_TIMING)
			{
			  //Add this run's duration and number to the graph data
		          RunStats currStats;

			  currStats.runDuration = *scalerDuration;
			  currStats.zeroScalerDuration = (*scalerDuration == 0);
			  currStats.runNumber = *run;
			  runStats.push_back(currStats);
			}

			total_run_time += *scalerDuration;
			total_runs++;
			last_run = *run;

                        finalRun = *run; //Updates every loop, so the final value is the final run

			// check if current run is on same day as last_run
			// have to save last_run to ints to avoid struct tm issues
			time_t run_start = *start; //The start long int value is converted to a time_t object
			run_tm = localtime(&run_start);         //convert unix time to tm structure //The time_t is converted to a C time struct
			//^ localtime returns a tm* data

			if (last_run_year ==  run_tm->tm_year){
				if (last_run_mon ==  run_tm->tm_mon){
					if (last_run_mday ==  run_tm->tm_mday){             // same day
						daily_Duration[iday]+= *scalerDuration;
					}
					else{                                               //different day
						ndays++;
						iday++;
						daily_Duration[iday]+= *scalerDuration;
						daily_day[iday]=run_tm->tm_mday;
						daily_mon[iday]=run_tm->tm_mon+1;
						daily_year[iday]=run_tm->tm_year+1900;
					}
				}
				else{                                                     //different month
					ndays++;
					iday++;
					daily_Duration[iday]+= *scalerDuration;
					daily_day[iday]=run_tm->tm_mday;
					daily_mon[iday]=run_tm->tm_mon+1;
					daily_year[iday]=run_tm->tm_year+1900;
				}
			}
			else{                                                          //different year
				ndays++;
				iday++;
				daily_Duration[iday]+= *scalerDuration;
				daily_day[iday]=run_tm->tm_mday;
				daily_mon[iday]=run_tm->tm_mon+1;
				daily_year[iday]=run_tm->tm_year+1900;
			}

			last_run_year=run_tm->tm_year;
			last_run_mon=run_tm->tm_mon;
			last_run_mday=run_tm->tm_mday;

		}
		//^End of Time-If's



	        if (firstRun == -1) //RC: Recording first run
		{
		  firstRun = *run; //Save the first run
		  firstTimeValue = *run_tm; //Save the first time value
		}

		finalTimeValue = *run_tm; //Updates every loop so that on the final loop it has the final value

		hrun->Fill(*run);

		// check different multiplicities
		if (CoinType[0] && *fMultip >= 3) hMultip0->Fill(*fMultip);
		if (CoinType[1]) hMultip1->Fill(*fMultip);
		if (CoinType[2]) hMultip2->Fill(*fMultip);
		if (CoinType[3]) hMultip3->Fill(*fMultip);
		if (CoinType[1] || CoinType[2] || CoinType[3]) hMultip4->Fill(*fMultip);

		if (DO_MULTIP_TABLE)
		{
		  //Add the run# and multiplicity value onto the table
		  pair<int, int> multipPair(*run, *fMultip);
		  multipTable.push_back(multipPair);
		}

		//RC
		//Working area for 3-panel events.
		//Intended to require  1-top, 2-bottom panels
		//TODO: move these constants to the header or a more appropriate place
		const int numTopPanels = 4;
		const int numBotXPanels = 6;
		const int numBotYPanels = 6;
		const int topPanelNumbers [numTopPanels] = {18, 19, 21, 22}; //Panel numbers that are on top
		const int botXPanelNumbers [numBotXPanels] = {1, 2, 3, 4, 5, 6}; //Panel numbers that are on bottom in x-direction
		const int botYPanelNumbers [numBotYPanels] = {7, 8, 9, 10, 11, 12}; //Panel numbers that are on bottom in y-direction
		classBCount++;
		bool atLeastOne = false;
		for (int h = 0; h < 32; h++)
		{
		        if (fQDC[h] >= fSWThresh[h])
			{
			      //cout << "fSWThresh[" << h << "]: " << fSWThresh[h] << endl;
			      atLeastOne = true;
			}
		}
		if (atLeastOne) {totalCutCount++;}
		if ((CoinType[1])) {classCCount++;} //RC
		if ((CoinType[0])) //Filter down to 2+ panel firings
		{
		        classDCount++;
		        //cout << "CoinType[0] event found." << endl;
		        if (*fMultip >= 3) //Events where at least 3 panels file
			{
			      //cout << "\t That event had multiplicity >= 3." << endl;
			      //Begin checking to see if 1-top and 2-bottom panels were hit
			      vector<int> hitTracker; //Fills up with the panel numbers that are hit in this event
			      for (int o = 0; o < 32; o++) //For each panel
			      {
				    //cout << "\t For panel o = " << o << endl;
				    //cout << "\t\t fQDC = " << fQDC[o] << endl;
				    //cout << "\t\t fSWThresh = " << fSWThresh[o] << endl;
				    if (fQDC[o] >= fSWThresh[o])
				    {
				          //Assumed: fQDC is some sort of raw data for this even organized by panel index
				          //Assumed: fSWThresh is some sort of threshold of current for an event (possibly specific to muons)
				          int eventPanelNum = PanelMap(o, *run);
					  //^Returns the panel NUMBER for the panel where the threshhold was exceeded
					  hitTracker.push_back(eventPanelNum); //Add this panel number to the list of hit panels
				    }
			      }
			      int numBotXHits = 0; //Ups for each bottom hit along x
			      int numBotYHits = 0; //Ups for each bottom hit along y
			      int numTopHits = 0; //Ups for each top hit
			      unsigned int hitIndex = 0;
			      for (unsigned int hitIndex = 0; hitIndex < hitTracker.size(); hitIndex++)
			      {
				     //Check each of the assumed 3 hits
				     for (int botIX = 0; botIX < numBotXPanels; botIX++) //Check each bottom panel to see if this event fired there
				     {
				           numBotXHits += (botXPanelNumbers[botIX] == hitTracker[hitIndex]) ? 1 : 0;
				           //^Add 1 to bot hits x if that hit was in a bottom panel along the x
				     }
				     for (int botIY = 0; botIY < numBotYPanels; botIY++) //Check each bottom panel to see if this event fired there
				     {
				           numBotYHits += (botYPanelNumbers[botIY] == hitTracker[hitIndex]) ? 1 : 0;
				           //^Add 1 to bot hits y if that hit was in a bottom panel along the y
				     }
				     for (int topI = 0; topI < numTopPanels; topI++) //Check each top panel to see if this event fired there
				     {
				           numTopHits += (topPanelNumbers[topI] == hitTracker[hitIndex]) ? 1 : 0;
				           //^Add 1 to top hits if that hit was in a top panel
				     }
			      }
			      // if (numBotXHits == 1 && numBotYHits == 1  && numTopHits == 1)
			      // {
				    //  //This is a satisfactory event!
				    //  totalThreePanelEvents++;
				    //  //cout << "Panel count: " << hitTracker.size() << endl;
				    //  //cout << "\tEvent run?: " << total_runs << endl;
				    //  writer << *run << " "; //Add the relevant variable data to the file
				    //  writer << *fEntry << " "; //Q: Really need to confirm this
				    //  writer << 0  << " "; //Q: Is this the write order and variable?
				    //  writer << *scalerDuration << " "; //Q: Scaler time vs scaler duration?
				    //  for (int p = 0; p < 32; p++) //Go through each panel
				    //  {
				    //        if (fQDC[p] > fSWThresh[p]) //Cut ones below threshhold
					  //  {
					  //        writer << fQDC[p] << " "; //Add in QDC value
					  //  }
					  //  else
					  //  {
					  //        writer << "0 "; //Add 0 value
					  //  }
				    //  }
				    //  writer << endl;
     			  //     }
			}
			if ((DO_HI_MULTIP_CUT) && (*fMultip >= HIGH_MULTIP_THRESHOLD)) //Events where at least N panels file
			{
			      cout << "Detected a high multiplicity (>=" << HIGH_MULTIP_THRESHOLD << ") event in run: #" << *run << endl;
			      //if (highMultipRuns.size() > 0 && highMultipRuns.back() != *run)
			      //{
			      //       highMultipRuns.push_back(*run); //Append that run number onto the vector
			      //}

			      bool alreadyListed = false; //Assume this is the first time this run has been registered
			      for (vector<int>::iterator k = highMultipRuns.begin(); k != highMultipRuns.end(); k++)
				{ //Iterate through all the current high-multip runs
				    if (*k == *run) //If ANY of those runs have this run number
				    {
				            alreadyListed = true; //This run has already been registered
				    }
			      }
			      if (alreadyListed == false) //Only add the run to the list if it does not exist yet
			      {
				    highMultipRuns.push_back(*run);
				    cout << "\tAdded run " << *run << " since it was not added before" << endl;
			      }
			}
		} //End of coinType 0's

		//End of RC:3-panel modifications
		//RC: Determining the range of run numbers
		//cout << "RC: Run Being Examined: " << *run << "; ";

		if ((CoinType[1]) && (*fMultip == 4)){      //two top and two bottom panels fired

		        if (DO_ZERO_RUN && *scalerDuration == 0)
	                {
			  setData.zeroDurationFourPanelRuns.push_back(*run);
			  cout << "Event inside a 4-panel run with duration 0! Run #" << *run << endl;
		        }


		        hMultip5->Fill(*fMultip);
			Int_t nPanel=0;
			Int_t hitPanels[4];    // will eventually be larger...

			for (int j=0; j<32; j++) {
				hrqdc[j]->Fill(fQDC[j]);
				if (fQDC[j] >= fSWThresh[j]) {
					hcqdc[j]->Fill(fQDC[j]);
					hitPanels[nPanel]=PanelMap(j,*run);
					nPanel++;
				}
			}
			iDet = iDetIndex(hitPanels[2],hitPanels[3],hitPanels[0],hitPanels[1]);
			hiDet->Fill(iDet);
			FourPanelHits[iDet]++;
			FourPanelHitsTOT++;
			daily_count[iday]++;

			// start is unix time - seconds since 1,1,1970
			// root time axis start is 1/1/95
			// difference is 788918400 s
			//cout << *start-788918400 << endl;
			ht1->Fill(*start-788918400);

		}
	} //END OF RUN LOOP

	cout << endl;
	cout << "RC: Total 3-panel (1 top + 1 bot_x + 1 bot_y) events...(Class A): " << totalThreePanelEvents << endl; //RC
        cout << "RC: Total events.......................................(Class B): " << classBCount << endl; //RC
	cout << "RC: Coin Type 1 (2 top + 2 bottom) events..............(Class C): " << classCCount << endl; //RC
	cout << "RC: Coin Type 0 (Multiplicity 2) events................(Class D): " << classDCount << endl; //RC
	cout << "------------------------------------------------------------------" << endl; //RC
	cout << "RC: Efficiency_1 (N_A/N_B): " << static_cast<double>(totalThreePanelEvents)/static_cast<double>(classBCount) << endl;
	cout << "RC: Efficiency_2 (N_A/N_C): " << static_cast<double>(totalThreePanelEvents)/static_cast<double>(classCCount) << endl;
	cout << "RC: Efficiency_3 (N_A/N_D): " << static_cast<double>(totalThreePanelEvents)/static_cast<double>(classDCount) << endl;
	cout << "RC: First run number: " << firstRun << endl;
	cout << "RC: Total Events where at least one panel meets the cut(Class E): " << totalCutCount << endl; //RC
	//writer.close(); //RC

	if (DO_RUN_TIMING)
	{
	  //Construct TGraph with the collected data
	  TCanvas* runCanvas = new TCanvas("runCanvas", "A Run Number vs Run Duration Graph", 500, 1000);
	  runCanvas->SetLogy();
	  runCanvas->Divide(1, 2);
	  //Top graph: duration v run #
	  runCanvas->cd(1);
	  gPad->SetLogy();
	  int runNumbers[runStats.size()];
	  int runDurations[runStats.size()];
	  for (int i = 0; i < runStats.size(); i++)
	  {
	    runNumbers[i] = runStats[i].runNumber;
	    runDurations[i] = runStats[i].runDuration;
	  }
	  TGraph* runGraph = new TGraph(runStats.size(), runNumbers, runDurations);
	  runGraph->SetTitle("Run Number vs. Run Duration");
	  runGraph->Draw("AC");

	  //Sorting
	  sort(runStats.rbegin(), runStats.rend());

	  cout << ">>> Longest-Duration Runs <<<" << endl;
	  setData.maxRunDuration = runStats[0].runDuration; //Save the longest run
	  for (int i = 0; i < 10; i++)
	    { cout << "Rank " << i+1 << ": " << (double)runStats[i].runDuration << "s; Run #" << runStats[i].runNumber << endl;}
	  cout << "<<< --------------------- >>>" << endl;

	  sort(runStats.begin(), runStats.end());

	  cout << ">>> Shortest-Duration Runs <<<" << endl;
	  for (int i = 0; i < 10; i++)
	    { cout << "Rank " << i+1 << ": " << (double)runStats[i].runDuration << "s; Run #" << runStats[i].runNumber << "; zeroScalerDuration? " << runStats[i].zeroScalerDuration << endl;}
	  cout << "<<< --------------------- >>>" << endl;

	  //Bottom graph: distribution of durations
	  runCanvas->cd(2);
	  gPad->SetLogy();
	  cout << "Maximum Run Duration in this Set: " << setData.maxRunDuration << endl;
	  TH1I* runDistrib = new TH1I("runDistrib", "Frequency of Run Durations", 30, 0, setData.maxRunDuration);
	  for (int a = 0; a < runStats.size(); a++)
	  {
	    runDistrib->Fill(runDurations[a]);
	  }
	  runDistrib->Draw();
	  runCanvas->Print(("muon-data/" + runSet.baseName + "/" + runSet.extName + "-" + RUN_TIMING_FILE_NAME).c_str(), "pdf");
	  delete runGraph;
	  delete runCanvas;


	}

	if (DO_HI_MULTIP_CUT)
	{
	  //Write high multiplicities list to file
	  ofstream highMultipOutput;
	  highMultipOutput.open(HIGH_MULTIP_OUTPUT_LIST_NAME, ios::trunc);
	  for (vector<int>::iterator i = highMultipRuns.begin(); i != highMultipRuns.end(); i++)
	  {
	      highMultipOutput << *i << endl;
	  }
	  highMultipOutput.close();
	}

	if (DO_MULTIP_TABLE)
	{
	  //Append run# and multiplicities pair to file
	  ofstream multipTableOutput;
	  multipTableOutput.open(MULTIP_TABLE_OUTPUT_NAME, ios::app);
	  for (vector<pair<int, int>>::iterator p = multipTable.begin(); p != multipTable.end(); p++)
	  {
	        multipTableOutput << p->first << " " << p->second << endl;
	  }
	  multipTableOutput.close();
	}

	//See what the numeric values are
	for (int g = 0; g < hMultip0->GetSize(); g++)
	{
	   // cout << "Bin " << g << ": " << hMultip0->GetBinContent(g) << endl;
	}


	//
	// calculate 4-track efficiency
	//
	Int_t MuonHits_0 = 0;
	Int_t MuonHits_1 = 0;
	Int_t MuonHits_4 = 0;
	Double_t MuonEff_4t_0=0;
	Double_t MuonEff_4t_1=0;
	Double_t MuonEff_4t_4=0;
	for (Int_t i=4; i<=15;i++ ){
		MuonHits_0=MuonHits_0+hMultip0->GetBinContent(i);
		MuonHits_1=MuonHits_1+hMultip1->GetBinContent(i);
		MuonHits_4=MuonHits_4+hMultip4->GetBinContent(i);
	}
	MuonEff_4t_0=hMultip0->GetBinContent(5)/MuonHits_0;
	MuonEff_4t_1=hMultip1->GetBinContent(5)/MuonHits_1;
	MuonEff_4t_4=hMultip4->GetBinContent(5)/MuonHits_4;
	cout << "MuonEff_4t_0= " << MuonEff_4t_0 << endl;
	cout << "MuonEff_4t_1= " << MuonEff_4t_1 << endl;
	cout << "MuonEff_4t_4= " << MuonEff_4t_4 << endl;


	//Agglomerate QDC graphs
  //Check that this isn't adding to the agglom file every time the macro is run
	if (DO_QDC_AGGLOM) //If agglomerating
	{
	  TFile file((QDC_AGGLOM_FILE_NAME + ".root").c_str(), "UPDATE"); //Open the accumulator file. UPDATE mode allows writing to file
	  TH1F* qdcAggloms[32]; //Declare memory for all 32 qdc plots
	  for (int w = 0; w < 32; w++) //Loop over each plot
	  {
	    qdcAggloms[w] = (TH1F*)file.Get(("hcqdc" + to_string(w)).c_str()); //Pull reference to each object

	    if (qdcAggloms[w] == nullptr) //
	    {
	      cout << "qdcAgglom[w=" << w << "] was returned with nullptr when trying to read." << endl;
	      cout << "qdcAgglom[w=" << w << "] is being written to for the first time." << endl;
	      qdcAggloms[w] = new TH1F(*hcqdc[w]); //Make the agglomerated qdc have the same data as the existing one
	    }
	    else
	    {
	      qdcAggloms[w]->Add(hcqdc[w]); //Critical step: add the histograms together
	    }

	    //Independent of previous state...
	    cout << "QDC test output from agglomerator: qdcAggloms[w=" << w << "] mean: " << qdcAggloms[w]->GetMean() << endl;
	    qdcAggloms[w]->Write(); //Write back to file
	    delete qdcAggloms[w];
	  }
	  file.Close();
	}

// verification
//Double_t time_sum_check=0.;
//for (Int_t i=0; i < ndays; i++){
//	time_sum_check += daily_Duration[i];
//	cout << i << " " << daily_Duration[i] << " " << daily_count[i] << endl;
//}
//cout << time_sum_check << endl;
	setData.totalTime = total_run_time;
	setData.fourPanelEvents = FourPanelHitsTOT;

	cout << "=========================================" << endl;
	cout << "Data Set: " << runSet.extName << endl;
	cout << "=========================================" << endl;
	cout << "total_run_time (sec) = "  << total_run_time  << endl;
	cout << "total nruns = " << total_runs << endl;
	cout << "total ndays = " << ndays << endl;
	cout << "fourPanelHitsTOT = "  << FourPanelHitsTOT  << endl;
	cout << "first run: " << firstRun  << endl;
	cout << "final run: " << finalRun << endl;
	cout << "start date: " << asctime(&firstTimeValue);
	cout << "end date: " << asctime(&finalTimeValue);
	cout << "total events: " << classBCount << endl;
	cout << "total hours: " << (total_run_time / 3600) << endl;
	cout << "4-panel muons/second: " << FourPanelHitsTOT/total_run_time << endl;
	cout << "4-panel %: " << ((double)FourPanelHitsTOT)/((double)classBCount) << endl;
	cout << "=========================================" << endl;

//-----------------------------------------------------------------------------------------
// save counts for each detector to file

	ofstream outputfile;
	outputfile.open(mDataFolder + "/" + runSet.baseName + "/" + runSet.extName + "-vetoAna_det" + SKIM_CUT_MODIFIER + ".txt");
	cout << "Saving ..._det.txt" << endl;

	for (int j=1; j<145; j++) {
		//cout << j << " " << FourPanelHits[j] << endl;
		outputfile << j << " " << FourPanelHits[j] << " " << total_run_time << endl;
	}
	outputfile.close();

//-----------------------------------------------------------------------------------------
// save counts for each DAY to file

	ofstream outputfile2;
	outputfile2.open(mDataFolder + "/" + runSet.baseName + "/" + runSet.extName + "-vetoAna_day" + SKIM_CUT_MODIFIER + ".txt");
	cout << "Saving ..._day.txt" << endl;

	for (int j=0; j<ndays; j++) {
		//cout << daily_day[j] << " " << daily_mon[j] << " " << daily_year[j] << " " << daily_count[j] << " " << daily_Duration[j] << endl;
		outputfile2 << daily_day[j] << " " << daily_mon[j] << " " << daily_year[j] << " " << daily_count[j] << " " << daily_Duration[j] << endl;
	}
	outputfile2.close();


//-----------------------------------------------------------------------------------------
// create plots


	plotMultip((mDataFolder + "/" + runSet.baseName + "/" + runSet.extName + "-multip" + SKIM_CUT_MODIFIER).c_str());
	plotQDCs((mDataFolder + "/" + runSet.baseName + "/" + runSet.extName + "-qdc" + SKIM_CUT_MODIFIER).c_str());
 	plotTimeHists();


//-----------------------------------------------------------------------------------------
// Save hiDet to file
	//simCompFile->Write();
	//simCompFile->Close();


	////////

    gBenchmark->Show("VetoAna");

    return setData;
}

//-----------------------------------------------------------------------------------------
