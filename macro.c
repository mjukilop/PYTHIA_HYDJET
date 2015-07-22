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
TLeaf *pfId;
TLeaf *pfPt;
TLeaf *eta;
TLeaf *phi;

int nEvents;
int nParticles;

TH1D histoPFID("histoPFID", "Particle Flow Candidate ID", 8, -0.5, 7.5);
TH1D histoSumPt("histoSumPt", "Sum Transverse Momentum", 12, 60, 120);
TH1D histoExcludedMuons("histoExcludedMuons", "How many events were excluded (0 = included, 1 = excluded)", 2, -0.5, 1.5);
TH1D histoMuonPt("histoMuonPt", "Individual muon transverse momentum", 40, 0, 200);
TH1D histoMuonPerEvent("histoMuonPerEvent", "Number of Muons per Event", 10, -0.5, 9.5);

void macro(const int job,const int nJobs) {

	// Start the timer
	timer.Start();

	// Setting tree and branch locations
	myFile = TFile::Open("root://xrootd.cmsaf.mit.edu//store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat220_Track9_Jet30_matchEqR_merged_forest_0.root");
	pfcandAnalyzer->cd();
	myTree = pfTree;
//	pfTree->Print();
//	gDirectory->ls();

	// Get number of entries
	nEvents = myTree->GetEntries();

	// Loop over the entries and perform actions
	for (Int_t ii = 100/nJobs*job; ii <= 100/nJobs*(job+1)-1/*nEvents*/; ii++) {
		myTree->GetEntry(ii);

		std::vector<int> muonPosition;
		std::vector< vector<int> > muonPairPosition;
		std::vector<int> muonExclusion;

		// Code to print out how far through I am
		cout << "Processing event : " << ii << "/" << nEvents << endl;

		pfId = myTree->GetLeaf("pfId");
		pfPt = myTree->GetLeaf("pfPt");
		eta = myTree->GetLeaf("pfEta");
		phi = myTree->GetLeaf("pfPhi");

		int muonCounter = 0;
		int muonPairExcludedCount = 0;
		
		nParticles = pfId->GetLen();

		// Fill histograms
		for (Int_t jj = 0; jj < nParticles; jj++) {
			histoPFID.Fill(pfId->GetValue(jj));

			// Find number of muons, add positions to vector
			if (pfId->GetValue(jj) == 3 /*&& pfPt->GetValue(jj) > 20*/) {
				muonPosition.push_back(jj);
				histoMuonPt.Fill(pfPt->GetValue(jj));
				muonCounter++;
			}
		}

		// Fill histogram showing number of muons per event
		histoMuonPerEvent.Fill(muonCounter);

		// Comparing muons to find close pairs.
		if (muonCounter >= 2) {
			outerLoop:
			for (Int_t jj = 0; jj < muonPosition.size(); jj++) {
				for (Int_t kk = 0; kk < muonPosition.size(); kk++) {
					if (jj != kk) {

						// Set variables to see if muon has paired more than once
						int muonPairNum = 0;
						double combinedPT = pfPt->GetValue(muonPosition[jj]) + pfPt->GetValue(muonPosition[kk]);

						// Check muons for pairs within 60-120 GeV, exclude muons which pair multiple times
						if (combinedPT > 60 && combinedPT < 120) {
							vector<int> temp(2);
							temp[0] = jj;
							temp[1] = kk;
							muonPairNum++;

							if (muonPairNum > 1) {
								muonPairPosition.pop_back();
								muonExclusion.push_back(jj);
								muonPairExcludedCount++;
								goto outerLoop;
							} else {
								muonPairPosition.push_back(temp);
								histoSumPt.Fill(combinedPT);
							}
						}
					}
				}
			}
		}

		// Fill number of excluded muon pair count
		histoExcludedMuons.Fill(muonPairExcludedCount);
	}

	// Close the file
	myFile->Close();

	// Open root file, store histograms and close it
	TFile out_file(Form("myhisto%d.root", job), "RECREATE");

	histoPFID.Write();
	histoSumPt.Write();
	histoExcludedMuons.Write();
	histoMuonPt.Write();
	histoMuonPerEvent.Write();

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
