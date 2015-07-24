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



void macro(const int job,const int nJobs) {

	// Declare variables
	TStopwatch timer;

	TFile *myFile;
	TTree *myTree;
	TLeaf *pfId;
	TLeaf *pfPt;
	TLeaf *eta;
	TLeaf *phi;

	int nEvents;
	int nParticles;
	int muonCounter;
	int muonPairExcludedCount;

	std::vector<int> muonPosition;
	std::vector< vector<int> > muonPairPosition;
	std::vector<int> muonExclusion;

	TH1D histoPFID("histoPFID", "Particle Flow Candidate ID", 8, -0.5, 7.5);
	TH1D histoSumPt("histoSumPt", "Sum Transverse Momentum", 12, 60, 120);
	TH1D histoExcludedMuons("histoExcludedMuons", "How many events were excluded (0 = included, 1 = excluded)", 2, -0.5, 1.5);
	TH1D histoMuonPt("histoMuonPt", "Individual muon transverse momentum (after cuts)", 40, 0, 200);
	TH1D histoMuonPtAll("histoMuonPtAll", "Individual muon transverse momentum (all)", 40, 0, 200);
	TH1D histoMuonPerEvent("histoMuonPerEvent", "Number of Muons per Event", 10, -0.5, 9.5);
	TH1D histoMuonsPassedEvent("histoMuonsPassedEvent", "0 = muon, 1 = didn't pass cut, 2 = passed cut", 3, -0.5, 2.5);

	// Start the timer
	timer.Start();

	// Setting tree and branch locations
	myFile = TFile::Open("root://xrootd.cmsaf.mit.edu//store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat170_Track9_Jet30_matchEqR_merged_forest_0.root");
	pfcandAnalyzer->cd();
	myTree = pfTree;
//	pfTree->Print();
//	gDirectory->ls();

	// Get number of entries
	nEvents = myTree->GetEntries();

	// Loop over the entries and perform actions
	for (Int_t ii = 100/nJobs*job; ii <= 100/nJobs*(job+1)-1/*nEvents*/; ii++) {
		
		// Access entry data
		myTree->GetEntry(ii);

		// Declare or reset event specific variables
		muonPosition.clear();
		muonPairPosition.clear();
		muonExclusion.clear();

		pfId = myTree->GetLeaf("pfId");
		pfPt = myTree->GetLeaf("pfPt");
		eta = myTree->GetLeaf("pfEta");
		phi = myTree->GetLeaf("pfPhi");

		muonCounter = 0;
		muonPairExcludedCount = 0;
		
		nParticles = pfId->GetLen();

		// Code to print out how far through I am
		cout << "Processing event : " << ii << "/" << nEvents << endl;

		// Loop through particles and find muons
		for (Int_t jj = 0; jj < nParticles; jj++) {
			histoPFID.Fill(pfId->GetValue(jj));

			// Find number of muons, add positions to vector, APPLY CUTS
			if (pfId->GetValue(jj) == 3) {
				if (pfPt->GetValue(jj) > 7 /*&& abs(eta->GetValue(jj)) < 2.4*/) {
					muonPosition.push_back(jj);
					// Debug
					histoMuonsPassedEvent.Fill(2);

					muonCounter++;
				}
				else {
					// Debug
					histoMuonsPassedEvent.Fill(1);
				}
				
				histoMuonsPassedEvent.Fill(0);
				histoMuonPtAll.Fill(pfPt->GetValue(jj));
			}
		}

		// Fill histogram showing number of muons per event
		histoMuonPerEvent.Fill(muonCounter);

		// Comparing muons to find close pairs.
		if (muonCounter >= 2) {
			// Debug
			TFile debug("debug.txt","RECREATE");
			debug << "We are in the main if with at least 2 muons. \n";

			Int_t jj = -1;
			goto outerLoop;
			outerLoop:

			while (jj < muonPosition.size()-1) {
				// Debug
				debug << "WHILE LOOP";

				// Increment loop at start to make up for goto breaking before end of while loop
				jj++;
				
				// Fill histogram and set variables for data storage and keeping track of times paired
				histoMuonPt.Fill(pfPt->GetValue(muonPosition[jj]));
				int sumPtStorage= 0;
				int muonPairNum = 0;

				for (Int_t kk = 0; kk < muonPosition.size(); kk++) {
					if (jj != kk) {

						// Find SumPt for later use
						double combinedPt = pfPt->GetValue(muonPosition[jj]) + pfPt->GetValue(muonPosition[kk]);

						// Check muons for pairs within 60-120 GeV, exclude muons which pair multiple times
						if (combinedPT > 60 && combinedPT < 120) {
							// Debug
							debug << "---------------------------COMBINED PT----------------------------";

							vector<int> temp(2);
							temp[0] = jj;
							temp[1] = kk;
							sumPtStorage = combinedPt;
							muonPairNum++;

							if (muonPairNum == 2) {
								// Debug
								debug << ">>>> EXCLUDED <<<<";

								muonPairPosition.pop_back();
								muonPairExcludedCount++;
								muonExclusion.push_back(jj);
								goto outerLoop;
							}
							else {
								muonPairPosition.push_back(temp);
							}
						}
					}
				}
				// Debug
				debug << "Filling sumPt";
				debug.Close();

				// Fill histogram if muon had only one pairing
				if (muonPairNum == 1) {
					histoSumPt.Fill(sumPtStorage);
				}
			}
		}
			
		// Fill number of excluded muon pair count
		histoExcludedMuons.Fill(muonPairExcludedCount);
			
/*			for (Int_t jj = 0; jj < muonPosition.size(); jj++) {
				histoMuonPt.Fill(pfPt->GetValue(jj));
				
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
*/	}

	// Close the file
	myFile->Close();

	// Open root file, store histograms and close it
	TFile out_file(Form("myhisto%d.root", job), "RECREATE");

	histoPFID.Write();
	histoSumPt.Write();
	histoExcludedMuons.Write();
	histoMuonPt.Write();
	histoMuonPtAll.Write();
	histoMuonPerEvent.Write();
	histoMuonsPassedEvent.Write();

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
