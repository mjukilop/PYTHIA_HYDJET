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
#include "TLeaf.h"
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
#include <sstream>
#include <iostream>
#include <fstream>
#include <TRandom.h>
#include <vector>
#include <stdio.h>

// Declare structs
struct jetData{
	string filePath;
	int job;
	double nJobs;
	std::vector<double>* motherPhi;
	TH1D* histoJetPt;
} jet;

struct pfData{
	string filePath;
	int job;
	double nJobs;
//	std::vector<int>* muonPosition;
//	std::vector< vector<int> >* muontPairPosition;
//	std::vector<int>* muonExclusion;
	std::vector<double>* motherMass;
	std::vector<double>* motherEta;
	std::vector<double>* motherPhi;
	TH1D* histoMotherInvMass;
	TH1D* histoMuonPt;
} pf;

void analysePF(struct pfData pf);
void analyseJet(struct jetData jet);

void macro(const int job,const double nJobs, const int file) {

	// Declare variables
	TStopwatch timer;

//	std::vector<int> muonPosition;
//	std::vector< vector<int> > muonPairPosition;
//	std::vector<int> muonExclusion;
	std::vector<double> motherMass;
	std::vector<double> motherPhi;
	std::vector<double> motherEta;
	std::vector<double> motherP;
	std::vector<double> motherE;
	
/*	// Debug
	std::vector<string> debug;
	stringstream stringStream;
	string str;
	int variables[] = { 10521, 12716, 15477, 18384, 1994, 21839, 23405, 2376, 25056, 25775, 26891, 27260, 27818, 29880, 30252, 30830, 3335, 3476, 36127, 36571, 36965, 3769, 39407, 39490, 3979, 39937, 40672, 40917, 42671, 44443, 46620, 47017, 51866, 54781, 54957, 55273, 55405, 56243, 57267, 57413, 60193, 6407, 65470, 67504, 67669, 68357, 69141, 69212, 69481, 69849, 71256, 71964, 72153, 72461, 7549, 9338 };
	std::vector<int> testVariables(56);
	for (int i = 0; i < testVariables.size(); i++) {
		testVariables[i] = variables[i];
	}
*/

	TH1D histoPFID("histoPFID", "Particle Flow Candidate ID", 8, -0.5, 7.5);
	TH1D histoSumPt("histoSumPt", "Sum Transverse Momentum", 12, 60, 120);
	TH1D histoExcludedMuons("histoExcludedMuons", "How many events were excluded (0 = included, 1 = excluded)", 2, -0.5, 1.5);
	TH1D histoMuonPt("histoMuonPt", "Individual muon transverse momentum (after cuts)", 200, 0, 200);
	TH1D histoMuonPtAll("histoMuonPtAll", "Individual muon transverse momentum (all)", 200, 0, 200);
	TH1D histoMuonPerEvent("histoMuonPerEvent", "Number of Muons per Event", 10, -0.5, 9.5);
	TH1D histoMuonsPassedEvent("histoMuonsPassedEvent", "0 = muon, 1 = didn't pass cut, 2 = passed cut", 3, -0.5, 2.5);
	TH1D histoMotherInvMass("histoMotherInvMass", "Invariant mass of dimuon mother particle", 150, 0, 150);

	TH1D histoJetPt("histoJetPt", "Jet transverse momentum for valid jets", 200, 0, 200);

	// Start the timer
	timer.Start();

	// Opening file and setting tree
	string inFile;
	inFile = "rootfiles.txt";
	ifstream instr(inFile.c_str(), std::ifstream::in);
	string fileName;
	// Reading up to and stop on our file
	for (int lineNum = 1; lineNum <= file; lineNum++) {
		instr >> fileName;
	}

	// analysePF
	pf.filePath = fileName;
	pf.job = job;
	pf.nJobs = nJobs;
//	pf.muonPosition = &muonPosition;
//	pf.muonPairPosition = &muonPairPosition;
//	pf.muonExclusion = &muonExclusion;
	pf.motherMass = &motherMass;
	pf.motherEta = &motherEta;
	pf.motherPhi = &motherPhi;
	pf.histoMotherInvMass = &histoMotherInvMass;
	pf.histoMuonPt = &histoMuonPt;

	analysePF(pf);

	// analyseJet
	jet.filePath = fileName;
	jet.job = job;
	jet.nJobs = nJobs;
	jet.motherPhi = &motherPhi;
	jet.histoJetPt = &histoJetPt;

	analyseJet(jet);

/*	cout << "File: " << fileName << endl;
//	cout << file << endl;
	myFile = TFile::Open(fileName.c_str());
//	myFile = TFile::Open("root://xrootd.cmsaf.mit.edu//store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat170_Track9_Jet30_matchEqR_merged_forest_0.root");

//	jetTree = (TTree*)myFile->Get("akPu3PFJetAnalyzer/t");
	pfTree = (TTree*)myFile->Get("pfcandAnalyzer/pfTree");
//	pfTree->Print();
//	gDirectory->ls();

	// Get number of entries
	nEvents = pfTree->GetEntries();


	// Loop over the entries and perform actions
	for (Int_t ii = (int)(100/nJobs*job); ii <= (int)(100/nJobs*(job+1)-1); ii++) {
	//for (Int_t pp = 0; pp < testVariables.size(); pp++) {
		//Int_t ii = testVariables[pp];

		// Access entry data
		pfTree->GetEntry(ii);

		// Declare or reset event specific variables
		muonPosition.clear();
		muonPairPosition.clear();
		muonExclusion.clear();

		pfId = pfTree->GetLeaf("pfId");
		pfPt = pfTree->GetLeaf("pfPt");
		pfEta = pfTree->GetLeaf("pfEta");
		pfPhi = pfTree->GetLeaf("pfPhi");

		muonCounter = 0;
		muonPairExcludedCount = 0;
		
		nParticles = pfId->GetLen();

		// Loop through particles and find muons
		for (Int_t jj = 0; jj < nParticles; jj++) {
			histoPFID.Fill(pfId->GetValue(jj));

			// Find number of muons, add positions to vector, APPLY CUTS
			if (pfId->GetValue(jj) == 3) {
				if (pfPt->GetValue(jj) > 7) {
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

		// Debug
		if (muonCounter >= 2) {
			stringStream.str("");
			stringStream << "We are looking at event " << ii << "\n";
			str = stringStream.str();
			debug.push_back(str);

			stringStream.str("");
			stringStream << muonCounter << " belongs to event " << ii << "\n";
			str = stringStream.str();
			debug.push_back(str);

			debug.push_back("/////////////////////////////////////END EVENT\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\n");
		}

		// Check muon combined pT and calculate mother particle invariant mass and plot data
		if (muonCounter >= 2) {
			
			int jj = -1;
			outerLoop:
			int whileLoopSize = muonPosition.size() - 1;
			
			while (jj < whileLoopSize) {
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
						if (combinedPt > 60 && combinedPt < 120) {
							// Debug
							debug.push_back("---------------------------COMBINED PT----------------------------\n");

							vector<int> temp(2);
							temp[0] = jj;
							temp[1] = kk;
							sumPtStorage = combinedPt;
							muonPairNum++;

							if (muonPairNum == 2) {
								// Debug
								debug.push_back(">>>> EXCLUDED <<<<\n");

								muonPairPosition.pop_back();
								muonPairExcludedCount++;
								muonExclusion.push_back(jj);
								goto outerLoop;
							}
							else {
								muonPairPosition.push_back(temp);
							}
						}

						// Calculating invariant mass of mother particle
						double eta1 = pfEta->GetValue(muonPosition[jj]);
						double eta2 = pfEta->GetValue(muonPosition[kk]);
						double phi1 = pfPhi->GetValue(muonPosition[jj]);
						double phi2 = pfPhi->GetValue(muonPosition[kk]);
						double muonInvMass = 0.1056583715; //GeV/c^2
						double pt1 = pfPt->GetValue(muonPosition[jj]);
						double pt2 = pfPt->GetValue(muonPosition[kk]);
						double e1;
						double e2;
						double px1;
						double px2;
						double py1;
						double py2;
						double pz1;
						double pz2;
						double e1;
						double e2;
						double p1;
						double p2;
						double motherInvMass;
						double motherPTemp;
						double motherETemp;
						double motherPhiTemp;
						double motherThetaTemp;
						double motherEtaTemp;

						double theta1 = 2 * atan(exp(-1 * eta1));
						double theta2 = 2 * atan(exp(-1 * eta2));

						px1 = pt1 * cos(phi1);
						px2 = pt2 * cos(phi2);
						py1 = pt1 * sin(phi1);
						py2 = pt2 * sin(phi2);
						pz1 = pt1 / tan(theta1);
						pz2 = pt2 / tan(theta2);
						motherP = sqrt((pow((px1 + px2), 2) + pow((py1 + py2), 2) + pow((pz1 + pz2), 2)));

						p1 = sqrt((pow(px1, 2) + pow(py1, 2) + pow(pz1, 2)));
						p2 = sqrt((pow(px2, 2) + pow(py2, 2) + pow(pz2, 2)));
						e1 = sqrt((pow(muonInvMass, 2) + pow(p1, 2)));
						e2 = sqrt((pow(muonInvMass, 2) + pow(p2, 2)));
						motherETemp = e1 + e2;

						motherInvMass = sqrt((pow(motherETemp, 2) - pow(motherPTemp, 2)));

						motherPhiTemp = tan((py1 + py2) / (px1 + px2));
						motherThetaTemp = tan((sqrt((px1 + px2) + (py1 + py2))) / (pz1 + pz2));
						motherEtaTemp = -1 * log(tan(motherThetaTemp / 2));

						// Store info about mother particle
						if (motherInvMass > 60 && motherInvMass < 120) {
							motherMass.push_back(motherInvMass);
							motherEta.push_back(motherEtaTemp);
							motherPhi.push_back(motherPhiTemp);

						}

						histoMotherInvMass.Fill(motherInvMass);
					}
				}

				// Fill histogram if muon had only one pairing
				if (muonPairNum == 1) {
					histoSumPt.Fill(sumPtStorage);
				}
			}
		}

		// Fill number of excluded muon pair count
		histoExcludedMuons.Fill(muonPairExcludedCount);
	}

	//	analysePF(pf);
	jet.motherPhi = &motherPhi;
	jet.histoJetPt = &histoJetPt;

	analyseJet(jet);

	// Close the file
	myFile->Close();

	ofstream debugtxt(Form("debug%d_file%d.txt", job, file), ios::out | ios::trunc);
	if (debugtxt.is_open()) {
		for (Int_t ii = 0; ii < debug.size(); ii++) {
			debugtxt << debug[ii];
		}
	}
*/

	// Open root file, store histograms and close it
	TFile out_file(Form("myhisto%d_file%d.root", job, file), "RECREATE");

	// Write histograms filled in analysePF
	histoPFID.Write();
	histoSumPt.Write();
	histoExcludedMuons.Write();
	histoMuonPt.Write();
	histoMuonPtAll.Write();
	histoMuonPerEvent.Write();
	histoMuonsPassedEvent.Write();
	histoMotherInvMass.Write();

	// Write histograms filled in analyseJet
	histoJetPt.Write();

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

void analyseJet(struct jetData arg) {
	cout << "analyseJet\n";
	// ------------------------------------- Declaring variables -------------------------------------
	// Open file
	TFile *myFile = TFile::Open(arg.filePath.c_str());

	// Set tree and leaves
	TTree *jetTree = (TTree*)myFile->Get("akPu3PFJetAnalyzer/t");
	TLeaf *jtpt = jetTree->GetLeaf("jtpt");
	TLeaf *jtPhi = jetTree->GetLeaf("jtphi");
//	TLeaf *jtEta = jetTree->GetLeaf("jteta");

	// Find number of events in tree
//	int nEvents = jetTree->GetEntries();
	
	// ------------------------------------- Doing calculations and filling histograms -------------------------------------
	// Loop over all events and perform actions
	for (Int_t ii = (int)(100 / arg.nJobs*arg.job); ii <= (int)(100 / arg.nJobs*(arg.job + 1) - 1); ii++) {
		
		// Access the ii'th event
		jetTree->GetEntry(ii);

		// Find motherPhi.length() for use in the loop
		int motherLength = arg.motherPhi->size();

		// Loop over all mother particles and compare them with each jet
		for (Int_t jj = 0; jj < motherLength; jj++) {
			if (abs(jtPhi->GetValue(jj) - (*(arg.motherPhi))[jj]) > 2 * M_PI / 3){
				arg.histoJetPt->Fill(jtpt->GetValue(ii));
			}
		}
	}
}

void analysePF(struct pfData arg) {
	
	// ------------------------------------- Declaring variables -------------------------------------
	// Open file
	TFile *myFile = TFile::Open(arg.filePath.c_str());

	// Set tree and leaves
	TTree *pfTree = (TTree*)myFile->Get("pfcandAnalyzer/pfTree");
	TLeaf *pfId = pfTree->GetLeaf("pfId");
	TLeaf *pfPt = pfTree->GetLeaf("pfPt");
	TLeaf *pfEta = pfTree->GetLeaf("pfEta");
	TLeaf *pfPhi = pfTree->GetLeaf("pfPhi");

	// Declare vectors to hold muon locations etc
	std::vector<int> muonPosition;
	std::vector< vector<int> > muonPairPosition;
	std::vector<int> muonExclusion;

	// Find number of events in tree
	int nEvents = pfTree->GetEntries();

	// ------------------------------------- Doing calculations and filling histograms -------------------------------------
	// Loop over all events and perform actions
	for (Int_t ii = (int)(100 / arg.nJobs*arg.job); ii <= (int)(100 / arg.nJobs*(arg.job + 1) - 1); ii++) {

		// Access the ii'th event
		pfTree->GetEntry(ii);

		// Set counters
		int muonCounter = 0;
		int muonPairExcludedCount = 0;
		cout << "counters set\n";

		// Set vector size for use later
		int position = muonPosition.size();
		cout << position << endl;

		// Clear the vectors used in the last event
		if (ii != 0) {
//			std::vector<int, allocator<int> >::iterator begin = pf.muonPosition->begin(), end = pf.muonPosition->end();
//			pf.muonPosition->erase(begin, end);
			cout << "erased\n";
			muonPosition.clear();
			muonPairPosition.clear();
			muonExclusion.clear();
		}
		cout << "cleared\n";
		// The number of particles in the event
		int nParticles = pfId->GetLen();

		// Loop through each particle in the event
		for (Int_t jj = 0; jj < nParticles; jj++) {

			// Find number of muons, add particle position to vector, apply cuts.
			if (pfId->GetValue(jj) == 3 && pfPt->GetValue(jj) > 7) {
				muonPosition.push_back(jj);
				muonCounter++;
				cout << "muon\n";

			}
		}

		// Find muons pairs mother particle invariant mass and plot data
		if (muonCounter >= 2) {
			 int jj = -1;
			 outperLoop:
			 int whileLoopSize = (muonPosition.size()) -1;

			 while (jj < whileLoopSize) {
				// Increment loop at the start to make up for the goto breaking before the end of the while loop
				jj++;

				// Fill histogram and set variables for data storage and keeping track of how many times a muon has been paired
				arg.histoMuonPt->Fill(pfPt->GetValue(muonPosition[jj]));

				// Loop over each muon pair
				for (unsigned Int_t kk = 0; kk < muonPosition.size(); kk++) {
					if (jj != kk) {

						// Calculating invariant mass
						double eta1 = pfEta->GetValue(muonPosition.at(jj));
						double eta2 = pfEta->GetValue(muonPosition.at(jj));
						double phi1 = pfPhi->GetValue((muonPosition.at(jj));
						double phi2 = pfPhi->GetValue((muonPosition.at(jj));
						double muonInvMass = 0.1056583715; //GeV/c^2
						double pt1 = pfPt->GetValue((muonPosition.at(jj));
						double pt2 = pfPt->GetValue((muonPosition.at(jj));
						double e1;
						double e2;
						double px1;
						double px2;
						double py1;
						double py2;
						double pz1;
						double pz2;
						double e1;
						double e2;
						double p1;
						double p2;
						double motherInvMass;
						double motherPTemp;
						double motherETemp;
						double motherPhiTemp;
						double motherThetaTemp;
						double motherEtaTemp;

						double theta1 = 2 * atan(exp(-1 * eta1));
						double theta2 = 2 * atan(exp(-1 * eta2));

						px1 = pt1 * cos(phi1);
						px2 = pt2 * cos(phi2);
						py1 = pt1 * sin(phi1);
						py2 = pt2 * sin(phi2);
						pz1 = pt1 / tan(theta1);
						pz2 = pt2 / tan(theta2);
						motherPTemp = sqrt((pow((px1 + px2), 2) + pow((py1 + py2), 2) + pow((pz1 + pz2), 2)));

						p1 = sqrt((pow(px1, 2) + pow(py1, 2) + pow(pz1, 2)));
						p2 = sqrt((pow(px2, 2) + pow(py2, 2) + pow(pz2, 2)));
						e1 = sqrt((pow(muonInvMass, 2) + pow(p1, 2)));
						e2 = sqrt((pow(muonInvMass, 2) + pow(p2, 2)));
						motherETemp = e1 + e2;

						motherInvMass = sqrt((pow(motherETemp, 2) - pow(motherPTemp, 2)));

						motherPhiTemp = tan((py1 + py2) / (px1 + px2));
						motherThetaTemp = tan((sqrt((px1 + px2) + (py1 + py2))) / (pz1 + pz2));
						motherEtaTemp = -1 * log(tan(motherThetaTemp / 2));

						// Store info about mother particle
						if (motherInvMass > 60 && motherInvMass < 120) {
							pf.motherMass.push_back(motherInvMass);
							pf.motherEta.push_back(motherEtaTemp);
							pf.motherPhi.push_back(motherPhiTemp);

						}

						arg.histoMotherInvMass.Fill(motherInvMass);
					}
				}
			}
		}
	}
}