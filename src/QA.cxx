//Isaac Mooney, WSU - June 2019
//This simple code produces some QA histograms for the pp, pA, and (eventually) AA data from Y12, Y15, and Y14 respectively

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TFile.h>

#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TString.h>
#include <TChain.h>
#include <TBranch.h>
#include <TMath.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <chrono>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequencePassiveArea.hh"
#include "fastjet/ClusterSequenceActiveArea.hh"
#include "fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh"
#include "fastjet/Selector.hh"

#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/FunctionOfPseudoJet.hh"

#include "fastjet/contrib/SoftDrop.hh"

#include "TStarJetPicoReader.h"
#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventHeader.h"
#include "TStarJetPicoEventCuts.h"

#include "TStarJetPicoPrimaryTrack.h"
#include "TStarJetPicoTower.h"
#include "TStarJetPicoTrackCuts.h"
#include "TStarJetPicoTowerCuts.h"

#include "TStarJetVectorContainer.h"
#include "TStarJetVector.h"
#include "TStarJetPicoTriggerInfo.h"
#include "TStarJetPicoUtils.h"

#include "funcs.hh"
#include "params.hh"

using namespace std;
using namespace fastjet;
using namespace Analysis;
typedef fastjet::contrib::SoftDrop SD;

// -------------------------
// Command line arguments: ( Defaults
// Defined for debugging in main )
// [0]: output directory
// [1]: name for the output root file containing histograms of observables
// [2]: string for selecting events which pass a trigger
// [3]: a dummy for now - if I add more features I may use it for something
// [3]: input data: can be a single .root or a .txt or .list of root files - should always be last argument

// DEF MAIN()
int main ( int argc, const char** argv) {
  
  //Start a timer
  TStopwatch TimeKeeper;
  TimeKeeper.Start( );
  
  //starting timing for function duration
  typedef std::chrono::high_resolution_clock clock;
  
  // Read in command line arguments
  // ------------------------------
  // Defaults
  std::string executable = "./bin/ppjetmass"; // placeholder
  std::string outputDir = "out/"; // directory where everything will be saved
  std::string outFileName = "test.root"; // histograms will be saved here
  std::string chainList = "list.txt"; // input file: can be .root, .txt, .list
  std::string chainName = "JetTree"; // Tree name in input file
  std::string trigger = "ppJP2";
  std::string dummy = " "; //This keeps the number of parameters the same between analysis segments for easy compatibility. If I later have a flag I want to pass in, it can go here.

  // Now check to see if we were given modifying arguments
  switch ( argc ) {
  case 1: // Default case
    __OUT("Using Default Settings");
    break;
  case 6: { // Custom case
    __OUT("Using Custom Settings");
    std::vector<std::string> arguments( argv+1, argv+argc );
    
    // Set non-default values
    // ----------------------
    
    // output and file names
    outputDir         = arguments[0];
    outFileName       = arguments[1];
    trigger           = arguments[2];
    dummy             = arguments[3]; //can be replaced by another flag if I want to add it in later
    chainList         = arguments[4];
    
    std::cout << "Running QA on the " << trigger << "-triggered data located in " << outputDir << "." << std::endl;
    std::cout << "The input file is " << chainList << " and the output file is " << outFileName << "." << std::endl;
    
    break;
  }
  default: { // Error: invalid custom settings
    __ERR("Invalid number of command line arguments");
    return -1;
    break;
  }
  }
  
  //in place for now; will encapsulate in a function if it gets much more involved
  int tID1 = -9999, tID2 = -9999;
  if (trigger == "ppJP2") {tID1 = tppJP2; tID2 = -8888;} //-8888 just ensures it won't accidentally match a trigger
  if (trigger == "ppVPDMB") {tID1 = tppVPDMB_nobsmd; tID2 = -8888;}
  if (trigger == "pAuJP2") {tID1 = tpAuJP2a; tID2 = tpAuJP2b;}
  if (trigger == "pAuBBCMB") {tID1 = tpAuBBCMBa; tID2 = tpAuBBCMBb;}

  // Build our input now
  // --------------------
  TChain* chain = new TChain( chainName.c_str() );
  
  // Check to see if the input is a .root file or a .txt
  bool inputIsRoot = Analysis::HasEnding( chainList.c_str(), ".root" );
  bool inputIsTxt  = Analysis::HasEnding( chainList.c_str(), ".txt"  );
  bool inputIsList = Analysis::HasEnding( chainList.c_str(), ".list" );
  
  // If its a recognized file type, build the chain
  // If its not recognized, exit
  if ( inputIsRoot ) { chain->Add( chainList.c_str() ); }
  else if ( inputIsTxt )  { chain = TStarJetPicoUtils::BuildChainFromFileList( chainList.c_str() ); }
  else if ( inputIsList)  { chain = TStarJetPicoUtils::BuildChainFromFileList( chainList.c_str() ); }
  else { __ERR("data file is not recognized type: .root or .txt only.") return -1; }
	
  // Build the event structure w/ cuts
  // ---------------------------------
  TStarJetPicoReader* reader = new TStarJetPicoReader();
  //NOTE: I am implementing the trigger in this file, rather than pulling it from src/params.hh. So in this case, what is set in that file is irrelevant!
  InitReader(reader, chain,100000 /*nEvents*/, "All" /*det_triggerString*/, det_absMaxVz, det_vZDiff, det_evPtMax, det_evEtMax, det_evEtMin, det_DCA, det_NFitPts, det_FitOverMaxPts, dat_maxEtTow, 0.9999, false, det_badTowers, dat_bad_run_list);
  
  // Data classes
  // ------------
  TStarJetVectorContainer<TStarJetVector>* container;
  TStarJetPicoEventHeader* header;
  TStarJetPicoEvent* event;
  
  // Histograms
  // ----------
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();

  //1D event observables
  TH1D* hn_trks = new TH1D("hn_trks","hn_trks",150,0,150);
  TH1D* hn_tows = new TH1D("hn_tows","hn_tows",250,0,500);
  TH1D* hbbc_coinc = new TH1D("hbbc_coinc","hbbc_coinc",200,0,2600000);
  TH1D* hevt_vtx = new TH1D("hevt_vtx","hevt_vtx",100,-35,35);
  TH1D* hvpdvz = new TH1D("hvpdvz","hvpdvz",100,-35,35);
  TH1D* hvzdiff = new TH1D("hvzdiff","hvzdiff",100,0,5);
  TH1D* hrunId = new TH1D("hrunId","hrunId",1000,13000000,14000000);//this is actually a dummy for now
  TH1D* hn_vertices = new TH1D("hn_vertices","hn_vertices",40,0,40);
  TH1D* hn_globals = new TH1D("hn_globals","hn_globals",200,0,4000);

  //1D raw track info (before QA cuts):
  TH1D* htrackNhits = new TH1D("htrackNhits","htrackNhits",50,0,50);
  TH1D* htrackNhitsposs = new TH1D("htrackNhitsposs","htrackNhitsposs",50,0,50);
  TH1D* htrackNhitsratio = new TH1D("htrackNhitsratio","htrackNhitsratio",100,0,1);

  //1D track observables
  TH1D* htrackPt = new TH1D("htrackPt","htrackPt",200,0,40);
  TH1D* htrackEta = new TH1D("htrackEta","htrackEta",40,-1,1);
  TH1D* htrackPhi = new TH1D("htrackPhi","htrackPhi",120,-M_PI,M_PI);
  TH1D* htrackDCA = new TH1D("htrackDCA","htrackDCA",200,0,1);

  //1D tower observables
  TH1D* htowerEt = new TH1D("htowerEt","htowerEt",200,0,40);
  TH1D* htowerEta = new TH1D("htowerEta","htowerEta",40,-1,1);
  TH1D* htowerPhi = new TH1D("htowerPhi","htowerPhi",120,-M_PI,M_PI);
  TH1D* htowerId = new TH1D("htowerId","htowerId",4800,0.5,4800.5);

  //2D event observables
  TH2D* hn_globals_n_primaries = new TH2D("hn_globals_n_primaries",";N_{glob};N_{prim}",400,0,4000,400,0,400);
  
  //2D track observables
  TH2D* htrackEta_Phi = new TH2D("htrackEta_Phi",";#eta_{trk};#phi_{trk}",40,-1,1,120,-M_PI,M_PI);
  TH2D* htrackDCA_n_trks = new TH2D("htrackDCA_n_trks",";track DCA;N_{trk}",200,0,1,150,0,150); //this one's empty for now - just annoying to produce
  
  //2D tower observables
  TH2D* htowerEta_Phi = new TH2D("htowerEta_Phi",";#eta_{tow};#phi_{tow}",40,-1,1,120,-M_PI,M_PI);
  TH2D* htowerId_Et = new TH2D("htowerId_Et",";#tower ID;E^{tow}_{T} [GeV]",4800,0.5,4800.5,200,0,40); //NEW!~~~

  //3D event, track, tower observables
  TH3D* htrackPt_Eta_Phi = new TH3D("htrackPt_Eta_Phi",";#p^{trk}_{T};#eta_{trk};#phi_{trk}",200,0,40,40,-1,1,120,-M_PI,M_PI);
  TH3D* htowerEt_Eta_Phi = new TH3D("htowerEt_Eta_Phi",";#E^{tow}_{T};#eta_{tow};#phi_{tow}",200,0,40,40,-1,1,120,-M_PI,M_PI);

  //2D BBC coincidence rate dependence of event observables
  TH2D* hbbc_coinc_evt_vtx = new TH2D("hbbc_coinc_evt_vtx",";BBC coincidence;v_{z} [cm]",200,0,2600000,100,-35,35);
  TH2D* hbbc_coinc_n_trks = new TH2D("hbbc_coinc_n_trks",";BBC coincidence;N_{trk}",200,0,2600000,150,0,150);
  TH2D* hbbc_coinc_n_tows = new TH2D("hbbc_coinc_n_tows",";BBC coincidence;N_{tow}",200,0,2600000,250,0,500);
  TH2D* hbbc_coinc_n_vertices = new TH2D("hbbc_coinc_n_vertices",";BBC coincidence;N_{vertices}",200,0,2600000,40,0,40);
  TH2D* hbbc_coinc_n_globals = new TH2D("hbbc_coinc_n_globals",";BBC coincidence;N_{globals}",200,0,2600000,200,0,4000);
  TH2D* hbbc_coinc_trackDCA = new TH2D("hbbc_coinc_trackDCA",";BBC coincidence;track DCA",200,0,2600000,200,0,1);
  TH2D* hbbc_coinc_trackPt = new TH2D("hbbc_coinc_trackPt",";BBC coincidence;p^{trk}_{T} [GeV/c]",200,0,2600000,200,0,40);
  TH2D* hbbc_coinc_towerEt = new TH2D("hbbc_coinc_towerEt",";BBC coincidence;E^{tow}_{T} [GeV]",200,0,2600000,200,0,40);
  
  //all runID histograms are currently in a very wide range (the whole calendar year), will narrow once I've seen the actual range in the plots
  TH2D* hrunId_trackPt = new TH2D("hrunId_trackPt",";run ID; p^{trk}_{T} [GeV/c]",1000,13015000,13075000,200,0,40);
  TH2D* hrunId_towerEt = new TH2D("hrunId_towerEt",";run ID; E^{tow}_{T} [GeV]",1000,13015000,13075000,200,0,40);
  TH2D* hrunId_towerId = new TH2D("hrunId_towerId",";run ID; tower ID",1000,13015000,13075000,4800,0.5,4800.5);
  
  //2D run ID dependence of events and tracks
  TH2D* hrunId_evt_vtx = new TH2D("hrunId_evt_vtx","; runID; v_{z} [cm]",1000,13015000,13075000,100,-35,35);
  TH2D* hrunId_bbc_coinc = new TH2D("hrunId_bbc_coinc","; runID; BBC coincidence [kHz]",1000,13015000,13075000,200,0,2600000);
  TH2D* hrunId_trackDCA = new TH2D("hrunId_trackDCA","; runID; track DCA",1000,13015000,13075000,200,0,1);
  
  //3D run ID dependence of tracks and towers
  TH3D* hrunId_trackEta_Phi = new TH3D("hrunId_trackEta_Phi",";runID;#eta_{trk};#phi_{trk}",1000,13015000,13075000,40,-1,1,120,-M_PI,M_PI);
  TH3D* hrunId_towerEta_Phi = new TH3D("hrunId_towerEta_Phi",";runID;#eta_{tow};#phi_{tow}",1000,13015000,13075000,40,-1,1,120,-M_PI,M_PI);

  
  //vector of hists for easy writing
  vector<TH1D*> events1D = {hn_trks,hn_tows,hbbc_coinc,hevt_vtx,hvpdvz,hvzdiff,hrunId,hn_vertices,hn_globals};
  vector<TH1D*> tracks1D = {htrackNhits,htrackNhitsposs,htrackNhitsratio,htrackPt,htrackEta,htrackPhi,htrackDCA};
  vector<TH1D*> towers1D = {htowerEt,htowerEta,htowerPhi,htowerId};
  vector<TH2D*> events2D = {hbbc_coinc_evt_vtx,hbbc_coinc_n_trks,hbbc_coinc_n_tows,hbbc_coinc_n_vertices,hbbc_coinc_n_globals,
			    hbbc_coinc_trackDCA,hbbc_coinc_trackPt,hbbc_coinc_towerEt,hrunId_evt_vtx,hrunId_bbc_coinc,hn_globals_n_primaries};
  vector<TH2D*> tracks2D = {htrackEta_Phi,htrackDCA_n_trks,hrunId_trackPt,hrunId_trackDCA};
  vector<TH2D*> towers2D = {htowerEta_Phi,htowerId_Et,hrunId_towerEt,hrunId_towerId};
  vector<TH3D*> events3D = {};
  vector<TH3D*> tracks3D = {htrackPt_Eta_Phi,hrunId_trackEta_Phi};
  vector<TH3D*> towers3D = {htowerEt_Eta_Phi,hrunId_towerEta_Phi};

  // Helpers
  // -------
  vector<PseudoJet> particles;
  
  //in theory there could also be jet QA for pp, so I will leave some jet machinery here just in case
  int nJets = 0;

  // Constituent selectors
  // ---------------------
  Selector select_track_rap = fastjet::SelectorAbsRapMax(max_track_rap);
  Selector select_lopt      = fastjet::SelectorPtMin( partMinPt );
  Selector select_loptmax   = fastjet::SelectorPtMax( partMaxPt );
  Selector slo = select_track_rap * select_lopt * select_loptmax;

  // Jet candidate selectors
  // -----------------------
  Selector select_jet_rap     = fastjet::SelectorAbsRapMax(max_rap);
  Selector select_jet_pt_min  = fastjet::SelectorPtMin( det_jet_ptmin );
  Selector select_jet_pt_max  = fastjet::SelectorPtMax( jet_ptmax );
  Selector select_jet_m_min = fastjet::SelectorMassMin( mass_min );
  Selector sjet = select_jet_rap && select_jet_pt_min && select_jet_pt_max && select_jet_m_min;
  
  // Choose a jet and area definition
  // --------------------------------
  JetDefinition jet_def = fastjet::JetDefinition(fastjet::antikt_algorithm, R);
  
  // create an area definition for the clustering - but not used for pp!
  //----------------------------------------------------------
  // ghosts should go up to the acceptance of the detector or
  // (with infinite acceptance) at least 2R beyond the region
  // where you plan to investigate jets.
  GhostedAreaSpec area_spec = fastjet::GhostedAreaSpec( ghost_maxrap, ghost_repeat, ghost_area );
  AreaDefinition  area_def = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts, area_spec);

  //Creating SoftDrop grooming object
  contrib::SoftDrop sd(Beta,z_cut,R0);
  cout << "SoftDrop groomer is: " << sd.description() << endl;
  
  cout << "Performing analysis." << endl;
  // Cycle through events
  // --------------------  
  //int nEvents = -1;
  int nEventsUsed = 0;

  //for later use looking up PDG masses using particle PID
  TDatabasePDG *pdg = new TDatabasePDG();

  unsigned n_trks = 0, n_tows = 0, n_evts = 0;
  
  try{
    while ( reader->NextEvent() ) {
      
      //clearing vectors if there are any, which there currently aren't
      //so nothing goes on this line for now

      //initialize variables to -9999 if there are any, which there currently aren't
      //so nothing goes on this line for now
      
      reader->PrintStatus(10);
      
      //get the event header
      event = reader->GetEvent();
      header = event->GetHeader();
      
      particles.clear();
      
      // Get the output container from the reader (not useful unless doing jet QA)
      // ----------------------------------------
      container = reader->GetOutputContainer();

      // if it HAS the trigger, do the QA for the triggered events. Should ALWAYS be satisfied for pp, since I selected the JP2 trigger in initialization!
      if (header->HasTriggerId(tID1) || header->HasTriggerId(tID2)) {//see function "SetTriggers()" for assignment of tID1, tID2 (or above until it's written) 
	++ n_evts;
	
	//~~~global observables~~~//
	//1Ds
	double n_tows_raw = header->GetNOfTowers();        //!this would be before cuts - see later for actual assignment [this isn't used currently]
	double n_primaries = header->GetNOfPrimaryTracks(); //!this would be before cuts - see later for actual assignment
	//first, assigning some info to variables for ease of repeated use later
	double bbc = header->GetBbcCoincidenceRate();
	double runId = header->GetRunId();
	hbbc_coinc->Fill(bbc); //BBC coincidence rate
        hn_vertices->Fill(header->GetNumberOfVertices()); //vertex multiplicity
	hn_globals->Fill(header->GetNGlobalTracks()); //global multiplicity
	hrunId->Fill(runId); //run ID
	//~~~event vertex information~~~//
	double evt_vtx = header->GetPrimaryVertexZ();
	double vpdvz = header->GetVpdVz();
	hevt_vtx->Fill(evt_vtx); //TPC event vertex
	hvpdvz->Fill(vpdvz); //VPD event vertex
	hvzdiff->Fill(abs(evt_vtx - vpdvz)); //|TPC - VPD Vz|
	//2Ds [trk and tow multiplicity vs. bbc coincidence to be filled later since cuts must be applied first]
	hbbc_coinc_evt_vtx->Fill(bbc, evt_vtx);
	hbbc_coinc_n_vertices->Fill(bbc, header->GetNumberOfVertices());
	hbbc_coinc_n_globals->Fill(bbc, header->GetNGlobalTracks());
	hrunId_evt_vtx->Fill(runId, evt_vtx);
	hrunId_bbc_coinc->Fill(runId, bbc);
	hn_globals_n_primaries->Fill(header->GetNGlobalTracks(), header->GetNOfPrimaryTracks());
	
	//~~~raw track observables~~~//
	for (int i = 0; i < header->GetNOfPrimaryTracks(); ++ i) {
	  double nhits = event->GetPrimaryTrack(i)->GetNOfFittedHits();
	  double nhitsposs = event->GetPrimaryTrack(i)->GetNOfPossHits();
          htrackNhits->Fill(nhits); //number of used pad hits
          htrackNhitsposs->Fill(nhitsposs); //number of possible pad hits
	  htrackNhitsratio->Fill(nhits/(double)nhitsposs); //ratio between the above two
        }//end raw track loop
	
	//~~~applying selection criteria~~~//
	TList *trks = reader->GetListOfSelectedTracks();
        TIter nxt_trk(trks);

	unsigned n_trks_in_evt = 0;
	//now looping over the tracks which pass selections (including hardcoded eta and pT cuts)
        while (TStarJetPicoPrimaryTrack* trk = (TStarJetPicoPrimaryTrack*) nxt_trk()) {
          if (fabs(trk->GetEta()) > 1.0 || trk->GetPt() < 0.2) {//!my cuts are on particle-level, not primary-level, so have to add these in by hand
            continue;
	  }
	  //1Ds
          htrackDCA->Fill(trk->GetDCA()); //distance of closest approach to the primary vertex
          htrackPt->Fill(trk->GetPt());
          htrackEta->Fill(trk->GetEta());
          htrackPhi->Fill(trk->GetPhi());
	  //2Ds
	  htrackEta_Phi->Fill(trk->GetEta(), trk->GetPhi());
	  hbbc_coinc_trackDCA->Fill(bbc, trk->GetDCA());
	  hbbc_coinc_trackPt->Fill(bbc, trk->GetPt());
	  hrunId_trackPt->Fill(runId, trk->GetPt());
	  hrunId_trackDCA->Fill(runId, trk->GetDCA());
	  //3Ds
	  htrackPt_Eta_Phi->Fill(trk->GetPt(), trk->GetEta(), trk->GetPhi());
	  hrunId_trackEta_Phi->Fill(runId, trk->GetEta(), trk->GetPhi());
	  
	  //  trackNhits.push_back(trk->GetNOfFittedHits()); //this would be after cuts - see earlier for actual assignment
          //trackNhitsposs.push_back(trk->GetNOfPossHits()); //this would be after cuts - see earlier for actual assignment
          ++n_trks; //number of tracks passing all cuts in the event. Will help get an average number / event later, for debugging output
	  ++n_trks_in_evt; //same as above but reset to zero for each event for a per-event count
        }//end track loop
	hn_trks->Fill(n_trks_in_evt); //trk multiplicity

	//~~~applying selection criteria~~~//
        TList *tows = reader->GetListOfSelectedTowers();
        TIter nxt_tow(tows);
	
	unsigned n_tows_in_evt = 0;
	//now looping over the towers which pass selections (including hardcoded eta and Et cuts)
        while (TStarJetPicoTower* tow = (TStarJetPicoTower*) nxt_tow()) {
          if (fabs(tow->GetEta()) > 1.0 || tow->GetEt() < 0.2)
            continue;
	  //1Ds
          htowerEta->Fill(tow->GetEta());
          htowerPhi->Fill(tow->GetPhi());
          htowerEt->Fill(tow->GetEt());
          htowerId->Fill(tow->GetId());
	  //2Ds
	  htowerEta_Phi->Fill(tow->GetEta(), tow->GetPhi());
	  htowerId_Et->Fill(tow->GetId(), tow->GetEt());
	  hbbc_coinc_towerEt->Fill(bbc, tow->GetEt());
	  hrunId_towerEt->Fill(runId, tow->GetEt());
	  hrunId_towerId->Fill(runId, tow->GetId());
	  //3Ds
	  htowerEt_Eta_Phi->Fill(tow->GetEt(), tow->GetEta(), tow->GetPhi());
	  hrunId_towerEta_Phi->Fill(runId, tow->GetEta(), tow->GetPhi());
	  
	  ++n_tows; //number of towers passing all cuts in the event
	  ++n_tows_in_evt;
        }//tower loop
	hn_tows->Fill(n_tows_in_evt); //tow multiplicity
      }// if statement (trigger)
    } // Event loop
  }catch ( std::exception& e) {
    std::cerr << "Caught " << e.what() << std::endl;
    return -1;
  }

  // Output
  // ------                                                                                                                                                                                               
  TFile* fout = new TFile((outputDir + outFileName).c_str(), "RECREATE");

  // Close up shop
  // -------------
  //fout->Write();
  //writing all histograms by looping over containing vectors
  for (int i = 0; i < events1D.size(); ++ i) {events1D[i]->Write();}
  for (int i = 0; i < tracks1D.size(); ++ i) {tracks1D[i]->Write();}
  for (int i = 0; i < towers1D.size(); ++ i) {towers1D[i]->Write();}
  for (int i = 0; i < events2D.size(); ++ i) {events2D[i]->Write();}
  for (int i = 0; i < tracks2D.size(); ++ i) {tracks2D[i]->Write();}
  for (int i = 0; i < towers2D.size(); ++ i) {towers2D[i]->Write();}
  for (int i = 0; i < events3D.size(); ++ i) {events3D[i]->Write();}
  for (int i = 0; i < tracks3D.size(); ++ i) {tracks3D[i]->Write();}
  for (int i = 0; i < towers3D.size(); ++ i) {towers3D[i]->Write();}
  
  fout->Close();

  //quick debugging info
  cout << "There were " << n_evts << " events passing all cuts. " << endl
       << "There were on average " << n_trks/(double)n_evts << " tracks passing all cuts per accepted event" << endl
       << "                  and " << n_tows/(double)n_evts << " towers passing all cuts per accepted event" << endl;

  cout << "Wrote to " << fout->GetName() << endl;
  cout << "Bye :)" << endl;
    
  return 0;
}
  
