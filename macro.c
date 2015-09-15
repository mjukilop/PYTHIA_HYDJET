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
#include <exception>
#include <stdio.h>
#include <string.h>

using namespace std;

// Declare structs
struct data{
	TFile* myFile;
	int job;
	double nJobs;
	Int_t eventNum;
	bool hasPlottedContour;
	std::vector<double>* motherMass;
	std::vector<double>* motherEta;
	std::vector<double>* motherPhi;
	std::vector<double>* motherTheta;
	std::vector<double>* motherPt;
	std::vector<double>* motherY;
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
	TH1D* histoPtDiff3Small;
	TH2D* histoPtCompare3Large;
	TH2D* histoPtCompare4Large;
	TH2D* histoPtCompare5Large;
	TH1D* histoZhaveJet3Small;
	TH1D* histoZhaveJet3Large;
	TH1D* histoPtDiff3Large;
	TH1D* histoEventPlanes;
	TH2D* histoJetContourSingle;
	TH2D* histoZContour;
	TH2D* histoJetContour;
	TH2D* histoZMuMuContourSingle;
	TH1D* histoXVtx;
	TH1D* histoPtDiff4Small;
	TH1D* histoPtDiff4Large;
	TH1D* histoPtDiff5Small;
	TH1D* histoPtDiff5Large;
	TH1D* histoZhaveJet4Small;
	TH1D* histoZhaveJet4Large;
	TH1D* histoZhaveJet5Small;
	TH1D* histoZhaveJet5Large;
	TH1D* histoPtRatio3Small;
	TH1D* histoPtRatio3Large;
	TH1D* histoPtRatio4Small;
	TH1D* histoPtRatio4Large;
	TH1D* histoPtRatio5Small;
	TH1D* histoPtRatio5Large;
	TH2D* histoPtRatioEvtPlane3Small;
	TH2D* histoPtRatioEvtPlane3Large;
	TH2D* histoPtRatioEvtPlane4Small;
	TH2D* histoPtRatioEvtPlane4Large;
	TH2D* histoPtRatioEvtPlane5Small;
	TH2D* histoPtRatioEvtPlane5Large;

} args;

void analysePF(struct data args);
void analyseJet(struct data args);
void analyseEventPlane(struct data args);
void analyseVertex(struct data args);
void jetContour(struct data args);
double evtPlane(struct data args);
double thetaFromEta(double eta);

void macro(const int job,const double nJobs, const int file) {
	try {
		// Declare variables
		TStopwatch timer;

		int eventNum = 0;

		std::vector<double> motherMass;
		std::vector<double> motherPhi;
		std::vector<double> motherTheta;
		std::vector<double> motherEta;
		std::vector<double> motherPt;
		std::vector<double> motherE;
		std::vector<double> motherY;

		args.hasPlottedContour = false;

		// 1-D histograms for muons and their parent particle
		TH1D histoSumPt("histoSumPt", "Sum Transverse Momentum", 200, 0, 200);
		histoSumPt.SetTitle("Transverse Momentum Sum of Muons;Transverse Momentum;Number of Events;");

		TH1D histoSumPtCut("histoSumPtCut", "Sum Muon Transverse Momentum for Z Events", 200, 0, 200);
		histoSumPtCut.SetTitle("Transverse Momentum Sum of Muon with a Z Parent;Transverse Momentum Sum;Number of Events;");

		TH1D histoMuonPtCut("histoMuonPtCut", "Individual muon transverse momentum (after cuts)", 200, 0, 200);
		histoMuonPtCut.SetTitle("Individual Muon Transverse Momentum, after cuts;Transverse Momentum;Number of Events;");

		TH1D histoMuonPt("histoMuonPt", "Individual muon transverse momentum (all)", 200, 0, 200);
		histoMuonPt.SetTitle("Individual Muon Transverse Momentum, All Muons;Transverse Momentum;Number of Events;");

		TH1D histoMuonPerEvent("histoMuonPerEvent", "Number of Muons per Event", 21, -0.5, 20.5);
		histoMuonPerEvent.SetTitle("Number of Muons per Event;Number of Muons;Number of Events;");

		TH1D histoMotherInvMass("histoMotherInvMass", "Invariant mass of dimuon mother particle", 200, 0, 200);
		histoMotherInvMass.SetTitle("Invariant Mass of Dimuon MotherParticle;Invariant Mass;Number of Particles;");

		TH1D histoMotherPhi("histoMotherPhi", "Phi of dimuon mother particle", 628, -1 * M_PI, M_PI);
		histoMotherPhi.SetTitle("Azimuthal Angle of Dimuon Mother Particle;Azimuthal Angle (Phi);Number of Particles;");

		TH1D histoMotherEta("histoMotherEta", "Eta of dimuon mother particle", 200, -10, 10);
		histoMotherEta.SetTitle("Pseudorapidity of Dimuon Mother Particle;Pseudorapidity;Number of Particles;");

		TH1D histoMotherTheta("histoMotherTheta", "Theta of dimuon mother particle", 314, 0, M_PI);
		histoMotherTheta.SetTitle("Polar Angle of Dimuon Mother Particle;Polar Angle (Theta);Number of Particles;");

		TH1D histoMotherPt("histoMotherPt", "Transverse momentum of dimuon mother particle (calculated component-wise from muons)", 200, 0, 200);
		histoMotherPt.SetTitle("Transverse Momentum of Dimuon Parent Particle;Transverse Momentum;Number of Particles;");

		// 1-D histograms for jets
		TH1D histoJetPt3Small("histoJetPt3Small", "Jet transverse momentum for valid jets, R=0.3", 200, 0, 200);
		histoJetPt3Small.SetTitle("Jet Transverse Momentum for Valid Jets, R=0.3;Transverse Momentum;Number of Jets;");

		TH1D histoJetPt4Small("histoJetPt4Small", "Jet transverse momentum for valid jets, R=0.4", 200, 0, 200);
		histoJetPt4Small.SetTitle("Jet Transverse Momentum for Valid Jets, R=0.4;Transverse Momentum;Number of Jets;");

		TH1D histoJetPt5Small("histoJetPt5Small", "Jet transverse momentum for valid jets, R=0.5", 200, 0, 200);
		histoJetPt5Small.SetTitle("Jet Transverse Momentum for Valid Jets, R=0.5;Transverse Momentum;Number of Jets;");

		TH1D histoJetB3Small("histoJetB3Small", "Jet 'b' parameter, explore", 500, -250, 250);
		histoJetB3Small.SetTitle("Plotting Jet 'b' Parameter;b;Number of Jets;");

		TH1D histoJetSumE3Small("histoSumE3Small", "Checking what eSum is", 500, -250, 250);
		histoJetSumE3Small.SetTitle("Plotting Jet 'eSum';eSum;Number of Jets;");

		TH1D histoJetPt3Large("histoJetPt3Large", "Jet transverse momentum for valid jets, wide angle, R=0.3", 200, 0, 200);
		histoJetPt3Large.SetTitle("Jet Transverse Momentum for Valid Jets with Wide Angle, R=0.3;Transverse Momentum;Number of Jets;");

		TH1D histoJetPt4Large("histoJetPt4Large", "Jet transverse momentum for valid jets, wide angle, R=0.4", 200, 0, 200);
		histoJetPt4Large.SetTitle("Jet Transverse Momentum for Valid Jets with Wide Angle, R=0.4;Transverse Momentum;Number of Jets;");

		TH1D histoJetPt5Large("histoJetPt5Large", "Jet transverse momentum for valid jets, wide angle, R=0.5", 200, 0, 200);
		histoJetPt5Large.SetTitle("Jet Transverse Momentum for Valid Jets with Wide Angle, R=0.5;Transverse Momentum;Number of Jets;");

		TH1D histoJetB3Large("histoJetB3Large", "Jet 'b' parameter, explore", 500, -250, 250);
		histoJetB3Large.SetTitle("Plotting Jet 'b' Parameter with Wide Angle;b;Number of Jets;");

		TH1D histoJetSumE3Large("histoSumE3Large", "Checking what eSum is", 500, -250, 250);
		histoJetSumE3Large.SetTitle("Plotting Jet 'eSum' Parameter with Wide Angle;eSum;Number of Jets;");

		// 1-D histograms from comparisons
		TH1D histoPtDiff3Small("histoPtDiff3Small", "The difference in transverse momentum between the Z and the highest pT jet, R=0.3", 200, -100, 100);
		histoPtDiff3Small.SetTitle("Transverse Momentum Difference, Jet - Z, R=0.3;Transverse Momentum Difference;Number of Events;");

		TH1D histoPtDiff3Large("histoPtDiff3Large", "The difference in transverse momentum between the Z and the highest pT jet, R=0.3", 200, -100, 100);
		histoPtDiff3Large.SetTitle("Transverse Momentum Difference, Jet - Z with Wide Angle, R=0.3;Transverse Momentum Difference;Number of Events;");

		TH1D histoZhaveJet3Small("histoZhaveJet3Small", "How many Z's have a jet in the opposite direction, 0=no, 1=yes", 3, -0.5, 2.5);
		histoZhaveJet3Small.SetTitle("Number of Z's with a back-to-back Jet. 0=No, 1=Yes;Back-to-back Z-Jet;Number of Events;");

		TH1D histoZhaveJet3Large("histoZhaveJet3Large", "How many Z's have a jet in the opposite direction, 0=no, 1=yes", 3, -0.5, 2.5);
		histoZhaveJet3Large.SetTitle("Number of Z's with a back-to-back Jet with Wide Angle. 0=No, 1=Yes;Back-to-back Z-Jet;Number of Events;");

		TH1D histoPtDiff4Small("histoPtDiff4Small", "The difference in transverse momentum between the Z and the highest pT jet, R=0.4", 200, -100, 100);
		histoPtDiff4Small.SetTitle("Transverse Momentum Difference, Jet - Z, R=0.4;Transverse Momentum Difference;Number of Events;");

		TH1D histoPtDiff4Large("histoPtDiff4Large", "The difference in transverse momentum between the Z and the highest pT jet, R=0.4", 200, -100, 100);
		histoPtDiff4Large.SetTitle("Transverse Momentum Difference, Jet - Z with Wide Angle, R=0.4;Transverse Momentum Difference;Number of Events;");

		TH1D histoZhaveJet4Small("histoZhaveJet4Small", "How many Z's have a jet in the opposite direction, 0=no, 1=yes", 3, -0.5, 2.5);
		histoZhaveJet4Small.SetTitle("Number of Z's with a back-to-back Jet. 0=No, 1=Yes;Back-to-back Z-Jet;Number of Events;");

		TH1D histoZhaveJet4Large("histoZhaveJet4Large", "How many Z's have a jet in the opposite direction, 0=no, 1=yes", 3, -0.5, 2.5);
		histoZhaveJet4Large.SetTitle("Number of Z's with a back-to-back Jet with Wide Angle. 0=No, 1=Yes;Back-to-back Z-Jet;Number of Events;");

		TH1D histoPtDiff5Small("histoPtDiff5Small", "The difference in transverse momentum between the Z and the highest pT jet, R=0.5", 200, -100, 100);
		histoPtDiff5Small.SetTitle("Transverse Momentum Difference, Jet - Z, R=0.5;Transverse Momentum Difference;Number of Events;");

		TH1D histoPtDiff5Large("histoPtDiff5Large", "The difference in transverse momentum between the Z and the highest pT jet, R=0.5", 200, -100, 100);
		histoPtDiff5Large.SetTitle("Transverse Momentum Difference, Jet - Z with Wide Angle, R=0.5;Transverse Momentum Difference;Number of Events;");

		TH1D histoZhaveJet5Small("histoZhaveJet5Small", "How many Z's have a jet in the opposite direction, 0=no, 1=yes", 3, -0.5, 2.5);
		histoZhaveJet5Small.SetTitle("Number of Z's with a back-to-back Jet. 0=No, 1=Yes;Back-to-back Z-Jet;Number of Events;");

		TH1D histoZhaveJet5Large("histoZhaveJet5Large", "How many Z's have a jet in the opposite direction, 0=no, 1=yes", 3, -0.5, 2.5);
		histoZhaveJet5Large.SetTitle("Number of Z's with a back-to-back Jet with Wide Angle. 0=No, 1=Yes;Back-to-back Z-Jet;Number of Events;");

		TH1D histoPtRatio3Small("histoPtRatio3Small", "", 1000, -5, 5);
		histoPtRatio3Small.SetTitle("Ratio of transverse momentum Jet/Z;Ratio;Number of Events");

		TH1D histoPtRatio3Large("histoPtRatio3Large", "", 1000, -5, 5);
		histoPtRatio3Large.SetTitle("Ratio of transverse momentum Jet/Z;Ratio;Number of Events");

		TH1D histoPtRatio4Small("histoPtRatio4Small", "", 1000, -5, 5);
		histoPtRatio4Small.SetTitle("Ratio of transverse momentum Jet/Z;Ratio;Number of Events");

		TH1D histoPtRatio4Large("histoPtRatio4Large", "", 1000, -5, 5);
		histoPtRatio4Large.SetTitle("Ratio of transverse momentum Jet/Z;Ratio;Number of Events");

		TH1D histoPtRatio5Small("histoPtRatio5Small", "", 1000, -5, 5);
		histoPtRatio5Small.SetTitle("Ratio of transverse momentum Jet/Z;Ratio;Number of Events");

		TH1D histoPtRatio5Large("histoPtRatio5Large", "", 1000, -5, 5);
		histoPtRatio5Large.SetTitle("Ratio of transverse momentum Jet/Z;Ratio;Number of Events");

		// 2-D histograms from comparisons
		TH2D histoPtCompare3Small("histoPtCompare3Small", "Comparing x=MotherPt, y=JetPt", 200, 0, 200, 200, 0, 200);
		histoPtCompare3Small.SetTitle("Comparing x=Z Transverse Momentum, y=Jet Transverse Momentum, R=0.3;Z Transverse Momentum;Jet Transverse Momentum;Number of Pairs");
		histoPtCompare3Small.SetOption("colz");

		TH2D histoPtCompare4Small("histoPtCompare4Small", "Comparing x=MotherPt, y=JetPt", 200, 0, 200, 200, 0, 200);
		histoPtCompare4Small.SetTitle("Comparing x=Z Transverse Momentum, y=Jet Transverse Momentum, R=0.4;Z Transverse Momentum;Jet Transverse Momentum;Number of Pairs");
		histoPtCompare4Small.SetOption("colz");

		TH2D histoPtCompare5Small("histoPtCompare5Small", "Comparing x=MotherPt, y=JetPt", 200, 0, 200, 200, 0, 200);
		histoPtCompare5Small.SetTitle("Comparing x=Z Transverse Momentum, y=Jet Transverse Momentum, R=0.5;Z Transverse Momentum;Jet Transverse Momentum;Number of Pairs");
		histoPtCompare5Small.SetOption("colz");

		TH2D histoPtCompare3Large("histoPtCompare3Large", "Comparing x=MotherPt, y=JetPt", 200, 0, 200, 200, 0, 200);
		histoPtCompare3Large.SetTitle("Comparing x=Z Transverse Momentum, y=Jet Transverse Momentum, Wide Angle, R=0.3;Z Transverse Momentum;Jet Transverse Momentum;Number of Pairs");
		histoPtCompare3Large.SetOption("colz");

		TH2D histoPtCompare4Large("histoPtCompare4Large", "Comparing x=MotherPt, y=JetPt", 200, 0, 200, 200, 0, 200);
		histoPtCompare4Large.SetTitle("Comparing x=Z Transverse Momentum, y=Jet Transverse Momentum, Wide Angle, R=0.4;Z Transverse Momentum;Jet Transverse Momentum;Number of Pairs");
		histoPtCompare4Large.SetOption("colz");

		TH2D histoPtCompare5Large("histoPtCompare5Large", "Comparing x=MotherPt, y=JetPt", 200, 0, 200, 200, 0, 200);
		histoPtCompare5Large.SetTitle("Comparing x=Z Transverse Momentum, y=Jet Transverse Momentum, Wide Angle, R=0.5;Z Transverse Momentum;Jet Transverse Momentum;Number of Pairs");
		histoPtCompare5Large.SetOption("colz");

		// 1-D histograms from event plane
		TH1D histoEventPlanes("histoEventPlanes", "Event Planes", 400, -20, 20);
		histoEventPlanes.SetTitle("Event Planes;Event Plane Angle? Phi?;Number of Planes;");

		// 2-D Contour plots
		TH2D histoJetContourSingle("histoJetContourSingle", "Jet pT weighted in x=eta and y=phi", 100, -5, 5, 100, -1 * M_PI, M_PI);
		histoJetContourSingle.SetTitle("Jet Transverse Momentum for one Z in x=Pseudorapidity, y=Azimuthal Angle;Pseudorapidity (Eta);Azimuthal Angle (Phi);");
		histoJetContourSingle.SetOption("colz");

		TH2D histoZContour("histoZContour", "", 400, -10, 10, 100, -1 * M_PI, M_PI);
		histoZContour.SetTitle("Z Transverse Momentum in x=Pseudorapidity, y=Azimuthal Angle;Pseudorapidity (Eta);Azimuthal Angle (Phi);");
		histoZContour.SetOption("colz");

		TH2D histoJetContour("histoJetContour", "", 100, -5, 5, 100, -1 * M_PI, M_PI);
		histoJetContour.SetTitle("Jet Transverse Momentum in x=Pseudorapidity, y=Azimuthal Angle;Pseudorapidity (Eta);Azimuthal Angle (Phi)");
		histoJetContour.SetOption("colz");

		TH2D histoZMuMuContourSingle("histoZMuMuContourSingle", "", 400, -10, 10, 100, -1 * M_PI, M_PI);
		histoZMuMuContourSingle.SetTitle("Z and Muon Transverse Momentum in x=Pseudorapidity, y=Azimuthal Angle;Pseudorapidity (Eta);Azimuthal Angle (Phi)");
		histoZMuMuContourSingle.SetOption("colz");

		// 1-D histograms from vertex
		TH1D histoXVtx("histoXVtx", "", 2000, -10, 10);
		histoXVtx.SetTitle("X Vertex Position of Particles;X Vertex Position;Number of Particles;");

		// 2-D contour plots of ptRatio vs evtPlane
		TH2D histoPtRatioEvtPlane3Small("histoPtRatioEvtPlane3Small", "", 100, 0, 2 * M_PI, 200, -10, 10);
		histoPtRatioEvtPlane3Small.SetTitle("Transverse Momentum Ratio, Jet/Z, w.r.t Phi, the Azimuthal Angle of the Z from the Event Plane, Clockwise;Phi, Angle Between the Jet and Event Plane in the Azimuth;Transverse Momentum Ratio of Jet/Z");
		histoPtRatioEvtPlane3Small.SetOption("colz");

		TH2D histoPtRatioEvtPlane3Large("histoPtRatioEvtPlane3Large", "", 100, 0, 2 * M_PI, 200, -10, 10);
		histoPtRatioEvtPlane3Large.SetTitle("Transverse Momentum Ratio, Jet/Z, w.r.t Phi, the Azimuthal Angle of the Z from the Event Plane, Clockwise;Phi, Angle Between the Jet and Event Plane in the Azimuth;Transverse Momentum Ratio of Jet/Z");
		histoPtRatioEvtPlane3Large.SetOption("colz");

		TH2D histoPtRatioEvtPlane4Small("histoPtRatioEvtPlane4Small", "", 100, 0, 2 * M_PI, 200, -10, 10);
		histoPtRatioEvtPlane4Small.SetTitle("Transverse Momentum Ratio, Jet/Z, w.r.t Phi, the Azimuthal Angle of the Z from the Event Plane, Clockwise;Phi, Angle Between the Jet and Event Plane in the Azimuth;Transverse Momentum Ratio of Jet/Z");
		histoPtRatioEvtPlane4Small.SetOption("colz");

		TH2D histoPtRatioEvtPlane4Large("histoPtRatioEvtPlane4Large", "", 100, 0, 2 * M_PI, 200, -10, 10);
		histoPtRatioEvtPlane4Large.SetTitle("Transverse Momentum Ratio, Jet/Z, w.r.t Phi, the Azimuthal Angle of the Z from the Event Plane, Clockwise;Phi, Angle Between the Jet and Event Plane in the Azimuth;Transverse Momentum Ratio of Jet/Z");
		histoPtRatioEvtPlane4Large.SetOption("colz");

		TH2D histoPtRatioEvtPlane5Small("histoPtRatioEvtPlane5Small", "", 100, 0, 2 * M_PI, 200, -10, 10);
		histoPtRatioEvtPlane5Small.SetTitle("Transverse Momentum Ratio, Jet/Z, w.r.t Phi, the Azimuthal Angle of the Z from the Event Plane, Clockwise;Phi, Angle Between the Jet and Event Plane in the Azimuth;Transverse Momentum Ratio of Jet/Z");
		histoPtRatioEvtPlane5Small.SetOption("colz");

		TH2D histoPtRatioEvtPlane5Large("histoPtRatioEvtPlane5Large", "", 100, 0, 2 * M_PI, 200, -10, 10);
		histoPtRatioEvtPlane5Large.SetTitle("Transverse Momentum Ratio, Jet/Z, w.r.t Phi, the Azimuthal Angle of the Z from the Event Plane, Clockwise;Phi, Angle Between the Jet and Event Plane in the Azimuth;Transverse Momentum Ratio of Jet/Z");
		histoPtRatioEvtPlane5Large.SetOption("colz");


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
			//std::getline(instr, fileName);
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
		args.motherTheta = &motherTheta;
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
		args.histoPtDiff3Small = &histoPtDiff3Small;
		args.histoPtCompare3Large = &histoPtCompare3Large;
		args.histoPtCompare4Large = &histoPtCompare4Large;
		args.histoPtCompare5Large = &histoPtCompare5Large;
		args.histoZhaveJet3Small = &histoZhaveJet3Small;
		args.histoZhaveJet3Large = &histoZhaveJet3Large;
		args.histoPtDiff3Large = &histoPtDiff3Large;
		args.histoEventPlanes = &histoEventPlanes;
		args.histoJetContourSingle = &histoJetContourSingle;
		args.histoZContour = &histoZContour;
		args.histoJetContour = &histoJetContour;
		args.histoZMuMuContourSingle = &histoZMuMuContourSingle;
		args.histoXVtx = &histoXVtx;
		args.motherY = &motherY;
		args.histoPtDiff4Small = &histoPtDiff4Small;
		args.histoPtDiff4Large = &histoPtDiff4Large;
		args.histoPtDiff5Small = &histoPtDiff5Small;
		args.histoPtDiff5Large = &histoPtDiff5Large;
		args.histoZhaveJet4Small = &histoZhaveJet4Small;
		args.histoZhaveJet4Large = &histoZhaveJet4Large;
		args.histoZhaveJet5Small = &histoZhaveJet5Small;
		args.histoZhaveJet5Large = &histoZhaveJet5Large;
		args.histoPtRatio3Small = &histoPtRatio3Small;
		args.histoPtRatio3Large = &histoPtRatio3Large;
		args.histoPtRatio4Small = &histoPtRatio4Small;
		args.histoPtRatio4Large = &histoPtRatio4Large;
		args.histoPtRatio5Small = &histoPtRatio5Small;
		args.histoPtRatio5Large = &histoPtRatio5Large;
		args.histoPtRatioEvtPlane3Small = &histoPtRatioEvtPlane3Small;
		args.histoPtRatioEvtPlane3Large = &histoPtRatioEvtPlane3Large;
		args.histoPtRatioEvtPlane4Small = &histoPtRatioEvtPlane4Small;
		args.histoPtRatioEvtPlane4Large = &histoPtRatioEvtPlane4Large;
		args.histoPtRatioEvtPlane5Small = &histoPtRatioEvtPlane5Small;
		args.histoPtRatioEvtPlane5Large = &histoPtRatioEvtPlane5Large;

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
		histoZContour.Write();
		histoZMuMuContourSingle.Write();

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
		histoPtDiff3Small.Write();
		histoZhaveJet3Small.Write();
		histoZhaveJet3Large.Write();
		histoPtDiff3Large.Write();
		histoJetContourSingle.Write();
		histoJetContour.Write();
		histoPtDiff4Small.Write();
		histoPtDiff4Large.Write();
		histoPtDiff5Small.Write();
		histoPtDiff5Large.Write();
		histoZhaveJet4Small.Write();
		histoZhaveJet4Large.Write();
		histoZhaveJet5Small.Write();
		histoZhaveJet5Large.Write();
		histoPtRatio3Small.Write();
		histoPtRatio3Large.Write();
		histoPtRatio4Small.Write();
		histoPtRatio4Large.Write();
		histoPtRatio5Small.Write();
		histoPtRatio5Large.Write();
		histoPtRatioEvtPlane3Small.Write();
		histoPtRatioEvtPlane3Large.Write();
		histoPtRatioEvtPlane4Small.Write();
		histoPtRatioEvtPlane4Large.Write();
		histoPtRatioEvtPlane5Small.Write();
		histoPtRatioEvtPlane5Large.Write();
		
		// Write histograms both
		histoPtCompare3Small.Write();
		histoPtCompare4Small.Write();
		histoPtCompare5Small.Write();
		histoPtCompare3Large.Write();
		histoPtCompare4Large.Write();
		histoPtCompare5Large.Write();
		
		// Write histograms filled in analyseEventPlane
		histoEventPlanes.Write();
		
		// Write histograms filled in analyseVertex
		histoXVtx.Write();
		
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
		catch (std::exception& e) {
		cout << "EXCEPTION CAUGHT! It says: " << e.what() << std::endl;
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
		int validMuonCounter = 0;
		int muonCounter = 0;

		// Clear the vectors used in the last event
		if (ii != 0) {
			muonPosition.clear();
			muonPairPosition.clear();
			muonExclusion.clear();
			arg.motherPt->clear();
			arg.motherPhi->clear();
			arg.motherMass->clear();
			arg.motherEta->clear();
			arg.motherY->clear();
			arg.motherTheta->clear();
		}

		// The number of particles in the event
		int nParticles = pfId->GetLen();

		// Loop through each particle in the event
		for (Int_t jj = 0; jj < nParticles; jj++) {

			// Find number of muons, add particle position to vector, apply cuts.
			if (pfId->GetValue(jj) == 3) {
				muonCounter++;
				arg.histoMuonPt->Fill(pfPt->GetValue(jj));
				if (pfPt->GetValue(jj) > 7) {
					muonPosition.push_back(jj);
					validMuonCounter++;
				}
			}
		}

		for (int jj = 0; jj < (int)(muonPosition.size() - 2); jj++) {
			int lowest = jj;
			for (int kk = 0; kk < (int)muonPosition.size(); kk++) {
				if (pfPt->GetValue(muonPosition.at(lowest)) > pfPt->GetValue(muonPosition.at(kk))) {
					lowest = kk;
				}
			}
			muonPosition.erase(muonPosition.begin() + lowest);
		}

		// Fill number of muon in event
		arg.histoMuonPerEvent->Fill(muonCounter);

		// Find muons pairs mother particle invariant mass and plot data
		if (validMuonCounter >= 2) {
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
					double motherYTemp;
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
						
					if ((px1 + px2) > 0) {
						if ((py1 + py2) > 0) {
							motherPhiTemp = atan((py1 + py2) / (px1 + px2));
						}
						else if ((py1 + py2) < 0) {
							motherPhiTemp = atan((py1 + py2) / (px1 + px2));
						}
					}
					else if ((px1 + px2) < 0) {
						if ((py1 + py2) > 0) {
							motherPhiTemp = (M_PI) + atan((py1 + py2) / (px1 + px2));
						}
						else if ((py1 + py2) < 0) {
							motherPhiTemp = atan((py1 + py2) / (px1 + px2)) - (M_PI);
						}
					}

					if ((pz1 + pz2) > 0) {
						motherThetaTemp = atan((sqrt(pow((px1 + px2), 2) + pow((py1 + py2), 2))) / (pz1 + pz2));
					}
					else if ((pz1 + pz2) < 0) {
						motherThetaTemp = (M_PI) + atan((sqrt(pow((px1 + px2), 2) + pow((py1 + py2), 2))) / (pz1 + pz2));
					}
					
					motherEtaTemp = -1 * log(tan(motherThetaTemp / 2));

					motherYTemp = 0.5 * log((motherETemp + (pz1 + pz2)) / (motherETemp - (pz1 + pz2)));

					// Store info about mother particle
					if (motherInvMass > 60 && motherInvMass < 120) {
						arg.motherMass->push_back(motherInvMass);
						arg.motherEta->push_back(motherEtaTemp);
						arg.motherPhi->push_back(motherPhiTemp);
						arg.motherTheta->push_back(motherThetaTemp);
						arg.motherPt->push_back(motherPtTemp);
						arg.motherY->push_back(motherYTemp);
						arg.histoSumPtCut->Fill(pt1 + pt2);
						arg.histoZContour->Fill(motherEtaTemp, motherPhiTemp, motherPtTemp);

						if (arg.hasPlottedContour == false) {
							arg.histoZMuMuContourSingle->Fill(eta1, phi1, pt1);
							arg.histoZMuMuContourSingle->Fill(eta2, phi2, pt2);
							arg.histoZMuMuContourSingle->Fill(motherEtaTemp, motherPhiTemp, motherPtTemp);
							jetContour(args);
						}
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
			analyseEventPlane(arg);
			analyseVertex(arg);
		}
	}
}

void analyseJet(struct data arg) {

	// ------------------------------------- Declaring variables -------------------------------------
	// Set tree and leaves
	TTree *jetTree3 = (TTree*)arg.myFile->Get("akPu3PFJetAnalyzer/t");
	TLeaf *jtPt3 = jetTree3->GetLeaf("jtpt");
	TLeaf *jtPhi3 = jetTree3->GetLeaf("jtphi");
	TLeaf *jtEta3 = jetTree3->GetLeaf("jteta");
	TLeaf *jtB3 = jetTree3->GetLeaf("b");
	TLeaf *jtSumE3 = jetTree3->GetLeaf("eSum");
	TLeaf *jtY3 = jetTree3->GetLeaf("jty");

	TTree *jetTree4 = (TTree*)arg.myFile->Get("akPu4PFJetAnalyzer/t");
	TLeaf *jtPt4 = jetTree4->GetLeaf("jtpt");
	TLeaf *jtPhi4 = jetTree4->GetLeaf("jtphi");
	TLeaf *jtEta4 = jetTree3->GetLeaf("jteta");
	TLeaf *jtY4 = jetTree4->GetLeaf("jty");

	TTree *jetTree5 = (TTree*)arg.myFile->Get("akPu5PFJetAnalyzer/t");
	TLeaf *jtPt5 = jetTree5->GetLeaf("jtpt");
	TLeaf *jtPhi5 = jetTree5->GetLeaf("jtphi");
	TLeaf *jtEta5 = jetTree3->GetLeaf("jteta");
	TLeaf *jtY5 = jetTree5->GetLeaf("jty");

	cout << "analysing Jet with eventNum = " << arg.eventNum << endl;

	// Set variable to hold number of jets
	int nEvents3;
	int nEvents4;
	int nEvents5;

	// Find motherPhi length for use in the loop
	int motherLength = arg.motherPhi->size();

	// Set int storing position of highest transverse momentum jet
	double smallestR = 9999;
	int smallestRPos = 0;

	// Bool to see if each Z has a corresonding jet
	bool zHaveJet = false;

	// ------------------------------------- Doing calculations and filling histograms -------------------------------------
	// Access the ii'th event, 3
	jetTree3->GetEntry(arg.eventNum);

	// Find number of events in TLeaf jtPhi
	nEvents3 = jtPhi3->GetLen();

	// Loop over all mother particles and compare them with each jet
	for (int ii = 0; ii < motherLength; ii++) {

		double mPhi = arg.motherPhi->at(ii);
		double mPt = arg.motherPt->at(ii);
		double mY = arg.motherY->at(ii);
		double mTheta = arg.motherTheta->at(ii);

		for (int jj = 0; jj < nEvents3; jj++) {

			double jetTheta = thetaFromEta((jtEta3->GetValue(jj)));
			double jetPhi = jtPhi3->GetValue(jj);
			double jetY = jtY3->GetValue(jj);
			double jetPt = jtPt3->GetValue(jj);
			double R = sqrt((jetTheta - mTheta) + (jetPhi - mPhi));

			//if ((abs(jtPhi3->GetValue(jj) - arg.motherPhi->at(ii))) >((2 * M_PI) / 3)) {
			if (sqrt((jetPhi - mPhi) + (jetY - mY)) > 0.4 && R > (2 * M_PI / 3)) {
				arg.histoJetPt3Small->Fill(jetPt);
				arg.histoJetB3Small->Fill(jtB3->GetValue(jj));
				arg.histoJetSumE3Small->Fill(jtSumE3->GetValue(jj));
				arg.histoPtCompare3Small->Fill(mPt, jetPt);
				zHaveJet = true;
				if (R < smallestR) {
					smallestR = R;
					smallestRPos = jj;
				}
			}
		}

		if (zHaveJet == true) {
			double ptDiff = mPt - jtPt3->GetValue(smallestRPos);
			arg.histoPtDiff3Small->Fill(ptDiff);

			double ptRatio = jtPt3->GetValue(smallestRPos) / mPt;
			arg.histoPtRatio3Small->Fill(ptRatio);

			double eventPlane = evtPlane(arg);

			double jetPhi = jtPhi3->GetValue(smallestRPos);
			if (jetPhi < 0) {
				jetPhi = jetPhi + (2 * M_PI);
			}

			double phiDiff = eventPlane - jetPhi;
			if (phiDiff < 0) {
				phiDiff = 2 * M_PI + phiDiff;
			}
			arg.histoPtRatioEvtPlane3Small->Fill(phiDiff, ptRatio);
		}

		arg.histoZhaveJet3Small->Fill(((int)zHaveJet));
		zHaveJet = false;

		smallestR = 9999;
		smallestRPos = 0;
	}



	// 4
	jetTree4->GetEntry(arg.eventNum);

	// Find number of events in TLeaf jtPhi
	nEvents4 = jtPhi4->GetLen();

	for (int ii = 0; ii < motherLength; ii++) {

		double mPhi = arg.motherPhi->at(ii);
		double mPt = arg.motherPt->at(ii);
		double mY = arg.motherY->at(ii);
		double mTheta = arg.motherTheta->at(ii);

		for (int jj = 0; jj < nEvents4; jj++) {

			double jetTheta = thetaFromEta((jtEta4->GetValue(jj)));
			double jetPhi = jtPhi4->GetValue(jj);
			double jetY = jtY4->GetValue(jj);
			double jetPt = jtPt4->GetValue(jj);
			double R = sqrt((jetTheta - mTheta) + (jetPhi - mPhi));

			//if ((abs(jtPhi4->GetValue(jj) - arg.motherPhi->at(ii))) >((2 * M_PI) / 3)) {
			if (sqrt((jetPhi - mPhi) + (jetY - mY)) > 0.4 && R > (2 * M_PI / 3)) {
				arg.histoJetPt4Small->Fill(jetPt);
				arg.histoPtCompare4Small->Fill(mPt, jetPt);
				zHaveJet = true;
				if (R < smallestR) {
					smallestR = R;
					smallestRPos = jj;
				}
			}
		}

		if (zHaveJet == true) {
			double ptDiff = mPt - jtPt4->GetValue(smallestRPos);
			arg.histoPtDiff4Small->Fill(ptDiff);

			double ptRatio = jtPt4->GetValue(smallestRPos) / mPt;
			arg.histoPtRatio4Small->Fill(ptRatio);

			double eventPlane = evtPlane(arg);

			double jetPhi = jtPhi4->GetValue(smallestRPos);
			if (jetPhi < 0) {
				jetPhi = jetPhi + (2 * M_PI);
			}

			double phiDiff = eventPlane - jetPhi;
			if (phiDiff < 0) {
				phiDiff = 2 * M_PI + phiDiff;
			}
			arg.histoPtRatioEvtPlane4Small->Fill(phiDiff, ptRatio);
		}

		arg.histoZhaveJet4Small->Fill(((int)zHaveJet));
		zHaveJet = false;

		smallestR = 9999;
		smallestRPos = 0;
	}

	// 5
	jetTree5->GetEntry(arg.eventNum);

	// Find number of events in TLeaf jtPhi
	nEvents5 = jtPhi5->GetLen();

	for (int ii = 0; ii < motherLength; ii++) {

		double mPhi = arg.motherPhi->at(ii);
		double mPt = arg.motherPt->at(ii);
		double mY = arg.motherY->at(ii);
		double mTheta = arg.motherTheta->at(ii);

		for (int jj = 0; jj < nEvents5; jj++) {

			double jetTheta = thetaFromEta((jtEta5->GetValue(jj)));
			double jetPhi = jtPhi5->GetValue(jj);
			double jetY = jtY5->GetValue(jj);
			double jetPt = jtPt5->GetValue(jj);
			double R = sqrt((jetTheta - mTheta) + (jetPhi - mPhi));

			//if ((abs(jtPhi5->GetValue(jj) - arg.motherPhi->at(ii))) >((2 * M_PI) / 3)) {
			if (sqrt((jetPhi - mPhi) + (jetY - mY)) > 0.4 && R > (2 * M_PI / 3)) {
				arg.histoJetPt5Small->Fill(jetPt);
				arg.histoPtCompare5Small->Fill(mPt, jetPt);
				zHaveJet = true;
				if (R < smallestR) {
					smallestR = R;
					smallestRPos = jj;
				}
			}
		}

		if (zHaveJet == true) {
			double ptDiff = mPt - jtPt5->GetValue(smallestRPos);
			arg.histoPtDiff5Small->Fill(ptDiff);

			double ptRatio = jtPt5->GetValue(smallestRPos) / mPt;
			arg.histoPtRatio5Small->Fill(ptRatio);

			double eventPlane = evtPlane(arg);

			double jetPhi = jtPhi5->GetValue(smallestRPos);
			if (jetPhi < 0) {
				jetPhi = jetPhi + (2 * M_PI);
			}

			double phiDiff = eventPlane - jetPhi;
			if (phiDiff < 0) {
				phiDiff = 2 * M_PI + phiDiff;
			}
			arg.histoPtRatioEvtPlane5Small->Fill(phiDiff, ptRatio);
		}

		arg.histoZhaveJet5Small->Fill(((int)zHaveJet));
		zHaveJet = false;

		smallestR = 9999;
		smallestRPos = 0;
	}

	// Do the same calculations and filling with a larger angle range for jets
	// Access the ii'th event, 3
	jetTree3->GetEntry(arg.eventNum);

	// Find number of events in TLeaf jtPhi
	nEvents3 = jtPhi3->GetLen();

	// Loop over all mother particles and compare them with each jet
	for (int ii = 0; ii < motherLength; ii++) {

		double mPhi = arg.motherPhi->at(ii);
		double mPt = arg.motherPt->at(ii);
		double mY = arg.motherY->at(ii);
		double mTheta = arg.motherTheta->at(ii);

		for (int jj = 0; jj < nEvents3; jj++) {

			double jetTheta = thetaFromEta((jtEta3->GetValue(jj)));
			double jetPhi = jtPhi3->GetValue(jj);
			double jetY = jtY3->GetValue(jj);
			double jetPt = jtPt3->GetValue(jj);
			double R = sqrt((jetTheta - mTheta) + (jetPhi - mPhi));

			//if ((abs(jtPhi3->GetValue(jj) - arg.motherPhi->at(ii))) >(M_PI / 2)) {
			if (sqrt((jetPhi - mPhi) + (jetY - mY)) > 0.4 && R > (M_PI / 2)) {
				arg.histoJetPt3Large->Fill(jetPt);
				arg.histoJetB3Large->Fill(jtB3->GetValue(jj));
				arg.histoJetSumE3Large->Fill(jtSumE3->GetValue(jj));
				arg.histoPtCompare3Large->Fill(mPt, jetPt);
				zHaveJet = true;
				if (R < smallestR) {
					smallestR = R;
					smallestRPos = jj;
				}
			}
		}

		if (zHaveJet == true) {
			double ptDiff = mPt - jtPt3->GetValue(smallestRPos);
			arg.histoPtDiff3Large->Fill(ptDiff);

			double ptRatio = jtPt3->GetValue(smallestRPos) / mPt;
			arg.histoPtRatio3Large->Fill(ptRatio);

			double eventPlane = evtPlane(arg);

			double jetPhi = jtPhi3->GetValue(smallestRPos);
			if (jetPhi < 0) {
				jetPhi = jetPhi + (2 * M_PI);
			}

			double phiDiff = eventPlane - jetPhi;
			if (phiDiff < 0) {
				phiDiff = 2 * M_PI + phiDiff;
			}
			arg.histoPtRatioEvtPlane3Large->Fill(phiDiff, ptRatio);
		}

		arg.histoZhaveJet3Large->Fill(((int)zHaveJet));
		zHaveJet = false;

		smallestR = 9999;
		smallestRPos = 0;
	}

	// 4
	jetTree4->GetEntry(arg.eventNum);

	// Find number of events in TLeaf jtPhi
	nEvents4 = jtPhi4->GetLen();

	for (int ii = 0; ii < motherLength; ii++) {

		double mPhi = arg.motherPhi->at(ii);
		double mPt = arg.motherPt->at(ii);
		double mY = arg.motherY->at(ii);
		double mTheta = arg.motherTheta->at(ii);

		for (int jj = 0; jj < nEvents4; jj++) {

			double jetTheta = thetaFromEta((jtEta4->GetValue(jj)));
			double jetPhi = jtPhi4->GetValue(jj);
			double jetY = jtY4->GetValue(jj);
			double jetPt = jtPt4->GetValue(jj);
			double R = sqrt((jetTheta - mTheta) + (jetPhi - mPhi));

			//if ((abs(jtPhi4->GetValue(jj) - arg.motherPhi->at(ii))) >(M_PI / 2)) {
			if (sqrt((jetPhi - mPhi) + (jetY - mY)) > 0.4 && R > (M_PI / 2)) {
				arg.histoJetPt4Large->Fill(jetPt);
				arg.histoPtCompare4Large->Fill(mPt, jetPt);
				zHaveJet = true;
				if (R < smallestR) {
					smallestR = R;
					smallestRPos = jj;
				}
			}
		}

		if (zHaveJet == true) {
			double ptDiff = mPt - jtPt4->GetValue(smallestRPos);
			arg.histoPtDiff4Large->Fill(ptDiff);

			double ptRatio = jtPt4->GetValue(smallestRPos) / mPt;
			arg.histoPtRatio4Large->Fill(ptRatio);

			double eventPlane = evtPlane(arg);

			double jetPhi = jtPhi4->GetValue(smallestRPos);
			if (jetPhi < 0) {
				jetPhi = jetPhi + (2 * M_PI);
			}

			double phiDiff = eventPlane - jetPhi;
			if (phiDiff < 0) {
				phiDiff = 2 * M_PI + phiDiff;
			}
			arg.histoPtRatioEvtPlane4Large->Fill(phiDiff, ptRatio);
		}

		arg.histoZhaveJet4Large->Fill(((int)zHaveJet));
		zHaveJet = false;

		smallestR = 9999;
		smallestRPos = 0;
	}

	// 5
	jetTree5->GetEntry(arg.eventNum);

	// Find number of events in TLeaf jtPhi
	nEvents5 = jtPhi5->GetLen();

	for (int ii = 0; ii < motherLength; ii++) {

		double mPhi = arg.motherPhi->at(ii);
		double mPt = arg.motherPt->at(ii);
		double mY = arg.motherY->at(ii);
		double mTheta = arg.motherTheta->at(ii);

		for (int jj = 0; jj < nEvents5; jj++) {

			double jetTheta = thetaFromEta((jtEta5->GetValue(jj)));
			double jetPhi = jtPhi5->GetValue(jj);
			double jetY = jtY5->GetValue(jj);
			double jetPt = jtPt5->GetValue(jj);
			double R = sqrt((jetTheta - mTheta) + (jetPhi - mPhi));

			//if ((abs(jtPhi5->GetValue(jj) - arg.motherPhi->at(ii))) >(M_PI / 2)) {
			if (sqrt((jetPhi - mPhi) + (jetY - mY)) > 0.4 && R > (M_PI / 2)) {
				arg.histoJetPt5Large->Fill(jetPt);
				arg.histoPtCompare5Large->Fill(mPt, jetPt);
				zHaveJet = true;
				if (R < smallestR) {
					smallestR = R;
					smallestRPos = jj;
				}
			}
		}

		if (zHaveJet == true) {
			double ptDiff = mPt - jtPt5->GetValue(smallestRPos);
			arg.histoPtDiff5Large->Fill(ptDiff);

			double ptRatio = jtPt5->GetValue(smallestRPos) / mPt;
			arg.histoPtRatio5Large->Fill(ptRatio);

			double eventPlane = evtPlane(arg);

			double jetPhi = jtPhi5->GetValue(smallestRPos);
			if (jetPhi < 0) {
				jetPhi = jetPhi + (2 * M_PI);
			}

			double phiDiff = eventPlane - jetPhi;
			if (phiDiff < 0) {
				phiDiff = 2 * M_PI + phiDiff;
			}
			arg.histoPtRatioEvtPlane5Large->Fill(phiDiff, ptRatio);
		}

		arg.histoZhaveJet5Large->Fill(((int)zHaveJet));
		zHaveJet = false;

		smallestR = 9999;
		smallestRPos = 0;
	}

	// Contour Plot for ak 3
	jetTree3->GetEntry(arg.eventNum);

	// Find number of events in TLeaf jtPhi
	nEvents3 = jtPhi3->GetLen();

	for (int ii = 0; ii < nEvents3; ii++) {
		if (jtPt3->GetValue(ii) < 400) {
			arg.histoJetContour->Fill(jtEta3->GetValue(ii), jtPhi3->GetValue(ii), jtPt3->GetValue(ii));
		}
	}
}

void analyseEventPlane(struct data arg) {

	// ------------------------------------- Declaring variables -------------------------------------
	// Set tree and leaves
	TTree *evtTree = (TTree*)arg.myFile->Get("hiEvtAnalyzer/HiTree");
	TLeaf *evtPlanes = evtTree->GetLeaf("hiEvtPlanes");

	// Filling histograms

	evtTree->GetEntry(arg.eventNum);

	int numEvtPlanes = evtPlanes->GetLen();

	cout << "eventPlanes has " << numEvtPlanes << " entries in event " << arg.eventNum << endl;

	for (int ii = 0; ii < numEvtPlanes; ii++) {
		arg.histoEventPlanes->Fill(evtPlanes->GetValue(ii));
	}
}

void analyseVertex(struct data arg) {

	// ------------------------------------- Declaring variables -------------------------------------
	// Set tree and leaves
	TTree *trackTree = (TTree*)arg.myFile->Get("pixelTrack/trackTree");
	TLeaf *xVtx = trackTree->GetLeaf("xVtx");
//	TLeaf *yVtx = trackTree->GetLeaf("yVtx");
//	TLeaf *zVtx = trackTree->GetLeaf("zVtx");

	trackTree->GetEntry(arg.eventNum);

	int numVtx = xVtx->GetLen();

	for (int ii = 0; ii < numVtx; ii++) {
		arg.histoXVtx->Fill(xVtx->GetValue(ii));
	}
}

void jetContour(struct data arg) {

	// ------------------------------------- Declaring variables -------------------------------------
	// Set tree and leaves
	TTree *jetTree3 = (TTree*)arg.myFile->Get("akPu3PFJetAnalyzer/t");
	TLeaf *jtPt3 = jetTree3->GetLeaf("jtpt");
	TLeaf *jtPhi3 = jetTree3->GetLeaf("jtphi");
	TLeaf *jtEta3 = jetTree3->GetLeaf("jteta");

	int nEvents3;

	// Contour Plot for ak 3

	jetTree3->GetEntry(arg.eventNum);

	// Find number of events in TLeaf jtPhi
	nEvents3 = jtPhi3->GetLen();

	if (arg.hasPlottedContour == false) {
		for (int ii = 0; ii < nEvents3; ii++) {
			arg.histoJetContourSingle->Fill(jtEta3->GetValue(ii), jtPhi3->GetValue(ii), jtPt3->GetValue(ii));
		}
		arg.hasPlottedContour = true;
	}
}

double evtPlane(struct data arg) {

	// ------------------------------------- Declaring variables -------------------------------------
	// Set tree and leaves
	TTree *evtTree = (TTree*)arg.myFile->Get("hiEvtAnalyzer/HiTree");
	TLeaf *evtPlane = evtTree->GetLeaf("NPhi0");

	evtTree->GetEntry(arg.eventNum);

	double temp = evtPlane->GetValue(0);
	if (temp < 0) {
		temp = temp + (2 * M_PI);
	}

	return temp;
}

double thetaFromEta(double eta) {

	double theta = 2 * atan(exp(-1 * eta));
	return theta;
}