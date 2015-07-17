#include "TStopwatch.h"
#include "TRandom3.h"
#include "TChain.h"
#include "TProfile.h"
#include "TCut.h"
#include "TLegend.h"
#include <cmath>
#include <cstdlib>
#include "TROOT.h"
#include "TTree.h"
#include "TBranch.h" 
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TLegend.h"
#include "TMemberInspector.h"
#include "TObjArray.h"
#include "TSeqCollection.h"
#include "TNamed.h"
#include <iostream>
#include <fstream>
#include <TRandom.h>
#include <vector>
#include <stdio.h>

TStopwatch timer;

TFile *myFile;
TTree *myTree;
TBranch *myBranch;
TLeaf *pfId;
TLeaf *sumpt;

int nEvents;
int nParticles;
int nSumPt;

TH1D histoPFID("histoPFID", "Particle Flow Candidate ID", 8, -0.5, 7.5);
TH1D histoSumPt("histoSumPt", "Sum Transverse Momentum", 100, -1000, 1000);

void macro(const int job,const int nJobs) {

	// Start the timer
	timer.Start();

	// Setting tree and branch locations
	myFile = TFile::Open("root://xrootd.cmsaf.mit.edu//store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat220_Track9_Jet30_matchEqR_merged_forest_0.root");
	pfcandAnalyzer->cd();
	myTree = pfTree;
	myBranch = pfTree->GetBranch("pfId");
//	pfTree->Print();
//	gDirectory->ls();

	// Get number of entries
	nEvents = myTree->GetEntries();

	// Loop over the entries and perform actions
	for (Int_t ii = 100/nJobs*job; ii <= 100/nJobs*(job+1)-1/*nEvents*/; ii++) {
		myTree->GetEntry(ii);

		explicit vector muonPosition(const allocator_type& alloc = allocator_type());

		// Code to print out how far through I am
		cout << "Processing event : " << ii << "/" << nEvents << endl;

		pfId = myTree->GetLeaf("pfId");
		sumpt = myTree->GetLeaf("sumpt");

		int muonCounter = 0;
		
		nParticles = pfId->GetLen();
		nSumPt = sumpt->GetLen();

		// Fill histograms
		for (Int_t jj = 0; jj < nParticles; jj++) {
			histoPFID.Fill(pfId->GetValue(jj));

			// Find number of muons
			if (pfId->GetValue(jj) == 3) {
				muonPosition[muonCounter] = jj;
				muonCounter++;
			}
		}

		for (Int_t jj = 0; jj < nSumPt; jj++) {
			histoSumPt.Fill(sumpt->GetValue(0));

		}

//		if (muonCounter >= 2) {
			
//		}
	}

	// Close the file
	myFile->Close();

	// Open root file, store histograms and close it
	TFile out_file(Form("myhisto%d.root", job), "RECREATE");

	histoPFID.Write();
	histoSumPt.Write();

	out_file.Close();

	// Stop timer and print time taken
	timer.Stop();
	float rTime = timer.RealTime();
	float cTime = timer.CpuTime();

	cout << "\t" << endl;
	cout << Form("RealTime=%f seconds, CpuTime=%f seconds", rTime, cTime) << endl;
	cout << "\t" << endl;
	cout << "Good bye : " << "\t" << endl;
}
