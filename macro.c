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
#include <string.h>

using namespace std;

// Declare structs
struct data{
	TFile* myFile;
	int job;
	double nJobs;
	Int_t eventNum;
	std::vector<double>* motherMass;
	std::vector<double>* motherEta;
	std::vector<double>* motherPhi;
	std::vector<double>* motherPt;
	TH1D* histoMotherInvMass;
	TH1D* histoMuonPtCut;
	TH1D* histoSumPt;
	TH1D* histoSumPtCut;
	TH1D* histoMuonPt;
	TH1D* histoMuonPerEvent;
	TH1D* histoMotherPhi;
	TH1D* histoMotherEta;
	TH1D* histoMotherTheta;
	TH1D* histoMotherPt;
	TH1D* histoJetPt3Small;
	TH1D* histoJetPt4Small;
	TH1D* histoJetPt5Small;
	TH1D* histoJetB3Small;
	TH1D* histoJetSumE3Small;
	TH2D* histoPtCompare3Small;
	TH2D* histoPtCompare4Small;
	TH2D* histoPtCompare5Small;
	TH1D* histoJetPt3Large;
	TH1D* histoJetPt4Large;
	TH1D* histoJetPt5Large;
	TH1D* histoJetB3Large;
	TH1D* histoJetSumE3Large;
	TH2D* histoPtCompare3Large;
	TH2D* histoPtCompare4Large;
	TH2D* histoPtCompare5Large;
} args;

void analysePF(struct data args);
void analyseJet(struct data args);

void macro(const int job,const double nJobs, const int file) {

	// Declare variables
	TStopwatch timer;

	int eventNum = 0;

	std::vector<double> motherMass;
	std::vector<double> motherPhi;
	std::vector<double> motherEta;
	std::vector<double> motherPt;
	std::vector<double> motherE;
	
	TH1D histoSumPt("histoSumPt", "Sum Transverse Momentum", 200, 0, 200);
	TH1D histoSumPtCut("histoSumPtCut", "Sum Muon Transverse Momentum for Z Events", 200, 0, 200);
	TH1D histoMuonPtCut("histoMuonPtCut", "Individual muon transverse momentum (after cuts)", 200, 0, 200);
	TH1D histoMuonPt("histoMuonPt", "Individual muon transverse momentum (all)", 200, 0, 200);
	TH1D histoMuonPerEvent("histoMuonPerEvent", "Number of Muons per Event", 10, -0.5, 9.5);
	TH1D histoMotherInvMass("histoMotherInvMass", "Invariant mass of dimuon mother particle", 200, 0, 200);
	TH1D histoMotherPhi("histoMotherPhi", "Phi of dimuon mother particle", 628, 0, 2 * M_PI);
	TH1D histoMotherEta("histoMotherEta", "Eta of dimuon mother particle", 150, 0, 150);
	TH1D histoMotherTheta("histoMotherTheta", "Theta of dimuon mother particle", 1256, -2 * M_PI, 2 * M_PI);
	TH1D histoMotherPt("histoMotherPt", "Transverse momentum of dimuon mother particle (calculated component-wise from muons)", 200, 0, 200);

	TH1D histoJetPt3Small("histoJetPt3Small", "Jet transverse momentum for valid jets, kt=0.3", 200, 0, 200);
	TH1D histoJetPt4Small("histoJetPt4Small", "Jet transverse momentum for valid jets, kt=0.4", 200, 0, 200);
	TH1D histoJetPt5Small("histoJetPt5Small", "Jet transverse momentum for valid jets, kt=0.5", 200, 0, 200);
	TH1D histoJetB3Small("histoJetB3Small", "Jet 'b' parameter, explore", 500, -250, 250);
	TH1D histoJetSumE3Small("histoSumE3Small", "Checking what eSum is", 500, -250, 250);
	TH1D histoJetPt3Large("histoJetPt3Large", "Jet transverse momentum for valid jets, wide angle, kt=0.3", 200, 0, 200);
	TH1D histoJetPt4Large("histoJetPt4Large", "Jet transverse momentum for valid jets, wide angle, kt=0.4", 200, 0, 200);
	TH1D histoJetPt5Large("histoJetPt5Large", "Jet transverse momentum for valid jets, wide angle, kt=0.5", 200, 0, 200);
	TH1D histoJetB3Large("histoJetB3Large", "Jet 'b' parameter, explore", 500, -250, 250);
	TH1D histoJetSumE3Large("histoSumE3Large", "Checking what eSum is", 500, -250, 250);

	TH2D histoPtCompare3Small("histoPtCompare3Small", "Comparing x=MotherPt, y=JetPt", 200, 0, 200, 200, 0, 200);
	TH2D histoPtCompare4Small("histoPtCompare4Small", "Comparing x=MotherPt, y=JetPt", 200, 0, 200, 200, 0, 200);
	TH2D histoPtCompare5Small("histoPtCompare5Small", "Comparing x=MotherPt, y=JetPt", 200, 0, 200, 200, 0, 200);
	TH2D histoPtCompare3Large("histoPtCompare3Large", "Comparing x=MotherPt, y=JetPt", 200, 0, 200, 200, 0, 200);
	TH2D histoPtCompare4Large("histoPtCompare4Large", "Comparing x=MotherPt, y=JetPt", 200, 0, 200, 200, 0, 200);
	TH2D histoPtCompare5Large("histoPtCompare5Large", "Comparing x=MotherPt, y=JetPt", 200, 0, 200, 200, 0, 200);

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
	TFile *myFile = TFile::Open(fileName.c_str());
	cout << "Opening file: " << fileName.c_str() << endl;

	// Set struct variables for use in functions
	args.myFile = myFile;
	args.eventNum = eventNum;
	args.job = job;
	args.nJobs = nJobs;
	args.motherMass = &motherMass;
	args.motherEta = &motherEta;
	args.motherPhi = &motherPhi;
	args.motherPt = &motherPt;
	args.histoMotherInvMass = &histoMotherInvMass;
	args.histoMuonPtCut = &histoMuonPtCut;
	args.histoSumPt = &histoSumPt;
	args.histoSumPtCut = &histoSumPtCut;
	args.histoMuonPt = &histoMuonPt;
	args.histoMuonPerEvent = &histoMuonPerEvent;
	args.histoMotherPhi = &histoMotherPhi;
	args.histoMotherEta = &histoMotherEta;
	args.histoMotherTheta = &histoMotherTheta;
	args.histoMotherPt = &histoMotherPt;
	args.histoJetPt3Small = &histoJetPt3Small;
	args.histoJetPt4Small = &histoJetPt4Small;
	args.histoJetPt5Small = &histoJetPt5Small;
	args.histoJetB3Small = &histoJetB3Small;
	args.histoJetSumE3Small = &histoJetSumE3Small;
	args.histoPtCompare3Small = &histoPtCompare3Small;
	args.histoPtCompare4Small = &histoPtCompare4Small;
	args.histoPtCompare5Small = &histoPtCompare5Small;
	args.histoJetPt3Large = &histoJetPt3Large;
	args.histoJetPt4Large = &histoJetPt4Large;
	args.histoJetPt5Large = &histoJetPt5Large;
	args.histoJetB3Large = &histoJetB3Large;
	args.histoJetSumE3Large = &histoJetSumE3Large;
	args.histoPtCompare3Large = &histoPtCompare3Large;
	args.histoPtCompare4Large = &histoPtCompare4Large;
	args.histoPtCompare5Large = &histoPtCompare5Large;

	// analysePF
	analysePF(args);

	// Open root file, store histograms and close it
	TFile out_file(Form("myhisto%d_file%d.root", job, file), "RECREATE");

	// Write histograms filled in analysePF
	histoSumPt.Write();
	histoSumPtCut.Write();
	histoMuonPtCut.Write();
	histoMuonPt.Write();
	histoMuonPerEvent.Write();
	histoMotherInvMass.Write();
	histoMotherPt.Write();
	histoMuonPt.Write();
	histoMotherPhi.Write();
	histoMotherEta.Write();
	histoMotherTheta.Write();

	// Write histograms filled in analyseJet
	histoJetPt3Small.Write();
	histoJetPt4Small.Write();
	histoJetPt5Small.Write();
	histoJetB3Small.Write();
	histoJetSumE3Small.Write();
	histoJetPt3Large.Write();
	histoJetPt4Large.Write();
	histoJetPt5Large.Write();
	histoJetB3Large.Write();
	histoJetSumE3Large.Write();

	// Write histograms both
	histoPtCompare3Small.Write();
	histoPtCompare4Small.Write();
	histoPtCompare5Small.Write();
	histoPtCompare3Large.Write();
	histoPtCompare4Large.Write();
	histoPtCompare5Large.Write();

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

void analyseJet(struct data arg) {

	// ------------------------------------- Declaring variables -------------------------------------
	// Set tree and leaves
	TTree *jetTree3 = (TTree*)arg.myFile->Get("akPu3PFJetAnalyzer/t");
	TLeaf *jtPt3 = jetTree3->GetLeaf("jtpt");
	TLeaf *jtPhi3 = jetTree3->GetLeaf("jtphi");
//	TLeaf *jtEta3 = jetTree3->GetLeaf("jteta");
	TLeaf *jtB3 = jetTree3->GetLeaf("b");
	TLeaf *jtSumE3 = jetTree3->GetLeaf("eSum");

	TTree *jetTree4 = (TTree*)arg.myFile->Get("akPu4PFJetAnalyzer/t");
	TLeaf *jtPt4 = jetTree4->GetLeaf("jtpt");
	TLeaf *jtPhi4 = jetTree4->GetLeaf("jtphi");

	TTree *jetTree5 = (TTree*)arg.myFile->Get("akPu5PFJetAnalyzer/t");
	TLeaf *jtPt5 = jetTree5->GetLeaf("jtpt");
	TLeaf *jtPhi5 = jetTree5->GetLeaf("jtphi");

	cout << "analysing Jet with eventNum = " << arg.eventNum << endl;

	// Find number of events in TLeaf jtPhi
	int nEvents3 = jtPhi3->GetLen();
	int nEvents4 = jtPhi4->GetLen();
	int nEvents5 = jtPhi5->GetLen();

	// Find motherPhi.length() for use in the loop
	int motherLength = arg.motherPhi->size();

	// ------------------------------------- Doing calculations and filling histograms -------------------------------------
	// Access the ii'th event, 3
	jetTree3->GetEntry(arg.eventNum);
	
	// Loop over all mother particles and compare them with each jet
	for (int ii = 0; ii < motherLength; ii++) {
		for (int jj = 0; jj < nEvents3; jj++) {
			if (abs(jtPhi3->GetValue(jj) - arg.motherPhi->at(ii)) > 2 * M_PI / 3) {
				arg.histoJetPt3Small->Fill(jtPt3->GetValue(jj));
				arg.histoJetB3Small->Fill(jtB3->GetValue(jj));
				arg.histoJetSumE3Small->Fill(jtSumE3->GetValue(jj));
				arg.histoPtCompare3Small->Fill(arg.motherPt->at(ii), jtPt3->GetValue(jj));
			}
		}
	}

	// 4
	jetTree4->GetEntry(arg.eventNum);

	for (int ii = 0; ii < motherLength; ii++) {
		for (int jj = 0; jj < nEvents4; jj++) {
			if (abs(jtPhi4->GetValue(jj) - arg.motherPhi->at(ii)) > 2 * M_PI / 3) {
				arg.histoJetPt4Small->Fill(jtPt4->GetValue(jj));
				arg.histoPtCompare4Small->Fill(arg.motherPt->at(ii), jtPt4->GetValue(jj));
			}
		}
	}

	// 5
	jetTree5->GetEntry(arg.eventNum);

	for (int ii = 0; ii < motherLength; ii++) {
		for (int jj = 0; jj < nEvents5; jj++) {
			if (abs(jtPhi5->GetValue(jj) - arg.motherPhi->at(ii)) > 2 * M_PI / 3) {
				arg.histoJetPt5Small->Fill(jtPt5->GetValue(jj));
				arg.histoPtCompare5Small->Fill(arg.motherPt->at(ii), jtPt5->GetValue(jj));
			}
		}
	}

	// Do the same calculations and filling with a larger angle range for jets
	// Access the ii'th event, 3
	jetTree3->GetEntry(arg.eventNum);

	// Loop over all mother particles and compare them with each jet
	for (int ii = 0; ii < motherLength; ii++) {
		for (int jj = 0; jj < nEvents3; jj++) {
			if (abs(jtPhi3->GetValue(jj) - arg.motherPhi->at(ii)) > M_PI / 2) {
				arg.histoJetPt3Large->Fill(jtPt3->GetValue(jj));
				arg.histoJetB3Large->Fill(jtB3->GetValue(jj));
				arg.histoJetSumE3Large->Fill(jtSumE3->GetValue(jj));
				arg.histoPtCompare3Large->Fill(arg.motherPt->at(ii), jtPt3->GetValue(jj));
			}
		}
	}

	// 4
	jetTree4->GetEntry(arg.eventNum);

	for (int ii = 0; ii < motherLength; ii++) {
		for (int jj = 0; jj < nEvents4; jj++) {
			if (abs(jtPhi4->GetValue(jj) - arg.motherPhi->at(ii)) > M_PI / 2) {
				arg.histoJetPt4Large->Fill(jtPt4->GetValue(jj));
				arg.histoPtCompare4Large->Fill(arg.motherPt->at(ii), jtPt4->GetValue(jj));
			}
		}
	}

	// 5
	jetTree5->GetEntry(arg.eventNum);

	for (int ii = 0; ii < motherLength; ii++) {
		for (int jj = 0; jj < nEvents5; jj++) {
			if (abs(jtPhi5->GetValue(jj) - arg.motherPhi->at(ii)) > M_PI / 2) {
				arg.histoJetPt5Large->Fill(jtPt5->GetValue(jj));
				arg.histoPtCompare5Large->Fill(arg.motherPt->at(ii), jtPt5->GetValue(jj));
			}
		}
	}
}

void analysePF(struct data arg) {
	
	// ------------------------------------- Declaring variables -------------------------------------
	// Set tree and leaves
	TTree *pfTree = (TTree*)arg.myFile->Get("pfcandAnalyzer/pfTree");
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
	for (Int_t ii = (int)(nEvents / arg.nJobs*arg.job); ii <= (int)(nEvents / arg.nJobs*(arg.job + 1) - 1); ii++) {

		// Access the ii'th event
		pfTree->GetEntry(ii);
		arg.eventNum = ii;
		if (ii % 1000 == 0) {
			cout << ii << "/" << nEvents << endl;
		}

		// Set counters
		int muonCounter = 0;

		// Clear the vectors used in the last event
		if (ii != 0) {
			muonPosition.clear();
			muonPairPosition.clear();
			muonExclusion.clear();
			arg.motherPt->clear();
		}

		// The number of particles in the event
		int nParticles = pfId->GetLen();

		// Loop through each particle in the event
		for (Int_t jj = 0; jj < nParticles; jj++) {

			// Find number of muons, add particle position to vector, apply cuts.
			if (pfId->GetValue(jj) == 3 && pfPt->GetValue(jj) > 7) {
				muonPosition.push_back(jj);
				muonCounter++;

				arg.histoMuonPt->Fill(pfPt->GetValue(jj));
			}
		}

		// Fill number of muon in event
		arg.histoMuonPerEvent->Fill(muonCounter);

		// Find muons pairs mother particle invariant mass and plot data
		if (muonCounter >= 2) {
			int jj = -1;
			int whileLoopSize = (muonPosition.size()) -1;

			while (jj < whileLoopSize) {
				// Increment loop at the start to make up for the goto breaking before the end of the while loop
				jj++;

				// Fill histogram
				arg.histoMuonPtCut->Fill(pfPt->GetValue(muonPosition[jj]));

				// Loop over each muon pair
				for (int kk = (jj+1); kk < (int)muonPosition.size(); kk++) {

					// Calculating invariant mass
					double eta1 = pfEta->GetValue(muonPosition.at(jj));
					double eta2 = pfEta->GetValue(muonPosition.at(kk));
					double phi1 = pfPhi->GetValue(muonPosition.at(jj));
					double phi2 = pfPhi->GetValue(muonPosition.at(kk));
					double muonInvMass = 0.1056583715; //GeV/c^2
					double pt1 = pfPt->GetValue(muonPosition.at(jj));
					double pt2 = pfPt->GetValue(muonPosition.at(kk));
					double e1;
					double e2;
					double px1;
					double px2;
					double py1;
					double py2;
					double pz1;
					double pz2;
					double p1;
					double p2;
					double motherInvMass;
					double motherPTemp;
					double motherETemp;
					double motherPtTemp;
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
					motherPtTemp = sqrt(pow((px1 + px2), 2) + pow((py1 + py2), 2));

					p1 = sqrt((pow(px1, 2) + pow(py1, 2) + pow(pz1, 2)));
					p2 = sqrt((pow(px2, 2) + pow(py2, 2) + pow(pz2, 2)));
					e1 = sqrt((pow(muonInvMass, 2) + pow(p1, 2)));
					e2 = sqrt((pow(muonInvMass, 2) + pow(p2, 2)));
					motherETemp = e1 + e2;

					motherInvMass = sqrt((pow(motherETemp, 2) - pow(motherPTemp, 2)));
						
					// debug
					cout << "MotherP = " << motherPTemp << "| MotherE = " << motherETemp << "| MotherInvMass = " << motherInvMass << endl;
					cout << "p1 pT = " << pt1 << "| p2 pT = " << pt2 << endl;
						
					motherPhiTemp = tan((py1 + py2) / (px1 + px2));
					motherThetaTemp = tan((sqrt((px1 + px2) + (py1 + py2))) / (pz1 + pz2));
					motherEtaTemp = -1 * log(tan(motherThetaTemp / 2));

					// Store info about mother particle
					if (motherInvMass > 60 && motherInvMass < 120) {
						arg.motherMass->push_back(motherInvMass);
						arg.motherEta->push_back(motherEtaTemp);
						arg.motherPhi->push_back(motherPhiTemp);
						arg.motherPt->push_back(motherPtTemp);
						arg.histoSumPtCut(pt1 + pt2);
						cout << "We have a big muon" << endl;

					}
					arg.histoMotherInvMass->Fill(motherInvMass);
					arg.histoMotherPhi->Fill(motherPhiTemp);
					arg.histoMotherEta->Fill(motherEtaTemp);
					arg.histoMotherTheta->Fill(motherThetaTemp);
					arg.histoMotherPt->Fill(motherPtTemp);
					arg.histoSumPt->Fill(pt1 + pt2);
				}
			}
			analyseJet(arg);
		}
	}
}