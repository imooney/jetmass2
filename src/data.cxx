//Isaac Mooney, WSU - June 2019
//This file will run the initial analysis on data for the jet mass project.
//It takes in the Picos, performs selections, clusters particles, performs selections on the resulting jets,
//applies the Soft Drop grooming procedure to a copy of the jet population, fills jet trees, and writes to files.
//It is currently capable of doing this for pp or pA and should be easily extendable to AA when I start working on it.

#include <iostream>
#include <iomanip>

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
// Command line arguments:                                                                                                                                    
// (Defaults defined for debugging)                                                                                                                           
// [0]: output directory                                                                                                                                      
// [1]: name for the output root file containing histograms of observables                                                                                    
// [2]: flag determining if we run over ch+ne (i.e. "full") jets or just ch jets. If it's "ch", runs over ch jets. If anything else, runs over full jets.
// [3]: flag determining trigger: ppJP2, ppHT2, ppVPDMB, pAuJP2, pAuHT2, pAuBBCMB, or AA[TBD]. Also used to determine collision species
// [4]: input data: can be a single .root or a .txt or .list of root files - should always be last argument 

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
  std::string executable = "./bin/data"; // placeholder                                                                                    
  std::string outputDir = "out/"; // directory where everything will be saved                                                                                  
  std::string outFileName = "test.root"; // histograms will be saved here                                                                                      
  std::string chainList = "list.txt"; // input file: can be .root, .txt, .list                                                                                 
  std::string chainName = "JetTree"; // tree name in input file                                                                                                
  std::string trigger = "ppJP2"; // trigger name: ppJP2, ppHT2, ppVPDMB, pAuJP2, pAuBBCMB, AA[TBD] 
  double radius = 0.4; //jet radius parameter; input value can range from 0.1 to 9.9.
  bool full = 1; //If TRUE, run over full (ch+ne) jets. If FALSE, run over ch jets.
  
  // Now check to see if we were given modifying arguments
  switch ( argc ) {
  case 1: // Default case
    __OUT("Using Default Settings");
      break;
  case 7: { // Custom case
    __OUT("Using Custom Settings");
      std::vector<std::string> arguments( argv+1, argv+argc );
            
    // Set non-default values
    // ----------------------
            
    // output and file names
    outputDir         = arguments[0];
    outFileName       = arguments[1];
    radius            = radius_str_to_double (arguments[2]);
    trigger           = arguments[3]; //ppJP2, ppHT2, ppVPDMB, pAuJP2, pAuHT2, pAuBBCMB, or AA [TBD]
    if (arguments[4] == "ch") {full = 0;} else {full = 1;} //either ch+ne jets (default) or ch jets (if "ch")
    chainList         = arguments[5];
            
    std::cout << "Running analysis of " << arguments[4] << " jets in the " << trigger << "-triggered data. Results will be output to " << outputDir << "." << std::endl;
    std::cout << "The input file is " << chainList << " and the output file is " << outFileName << "." << std::endl;
    
    break;
  }
  default: { // Error: invalid custom settings
    __ERR("Invalid number of command line arguments");
      return -1;
    break;
  }
  }

  //Setting up specifics of analysis based on the flags that were received above!
  string badtows = "", badruns = "";
  double vzdiff = -1;
  if (trigger.find("pp") != string::npos) {badtows = det_badTowers; badruns = dat_bad_run_list; vzdiff = det_vZDiff;}
  //TEMP (COMBINED_BADTOWERS FOR COMPARISON TO EMBEDDING):
  if (trigger.find("pA") != string::npos) {badtows = /*pAu_badTowers*/combined_badTowers; badruns = pAu_bad_run_list; vzdiff = pAu_vZDiff;}
  if (trigger.find("AA") != string::npos) {badtows = ""; badruns = ""; vzdiff = -1;} //TBD

  //in place for now; will encapsulate in a function if it gets much more involved. Hardcodes the trigger IDs.                                                 
  int tID1 = -9999, tID2 = -9999, tID3 = -9999;
  if (trigger == "ppJP2") {tID1 = tppJP2; tID2 = -8888; tID3 = -8888;} //-8888 just ensures it won't accidentally match a trigger
  if (trigger == "ppHT2") {tID1 = tppHT2a; tID2 = tppHT2b; tID3 = tppHT2c;}
  if (trigger == "ppVPDMB") {tID1 = tppVPDMB_nobsmd; tID2 = -8888; tID3 = -8888;}
  if (trigger == "pAuJP2") {tID1 = tpAuJP2a; tID2 = tpAuJP2b; tID3 = -8888;}
  if (trigger == "pAuHT2") {tID1 = tpAuHT2a; tID2 = tpAuHT2b; tID3 = -8888;}
  if (trigger == "pAuBBCMB") {tID1 = tpAuBBCMBa; tID2 = tpAuBBCMBb; tID3 = -8888;}
  
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
  TStarJetPicoReader * reader = new TStarJetPicoReader();
  InitReader(reader, chain, nEvents, "All"/*det_triggerString*/, det_absMaxVz, vzdiff, det_evPtMax, det_evEtMax, det_evEtMin, det_DCA, det_NFitPts, det_FitOverMaxPts, dat_maxEtTow, 0.9999, false, badtows, badruns);

  // Data classes
  // ------------
  TStarJetVectorContainer<TStarJetVector>* container;
  TStarJetVector* sv; // TLorentzVector* would be sufficient
  TStarJetPicoEventHeader* header;
  TStarJetPicoEvent* event;
  TStarJetPicoEventCuts* EventCuts = reader->GetEventCuts();//will use this for hardcoding checks if events passed selections
  
  // Histograms
  // ----------
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();
  
  // Trees
  // -----
  int dummy_int;
  double dummy_double;

  //variables to link to branches in the tree
  double n_jets; double refMult, bbc_east_rate, bbc_east_sum, zdc_coinc_rate, UE_avgpt, UEnpart;
  vector<double> Pt; vector<double> Eta; vector<double> Phi; vector<double> M; vector<double> E;
  vector<double> ch_e_frac;
  vector<double> zg; vector<double> rg; vector<double> mg; vector<double> ptg;
  vector<double> mcd;
  vector<double> tau0; vector<double> tau05; vector<double> tau_05; vector<double> tau_1;
  vector<double> tau0_g; vector<double> tau05_g; vector<double> tau_05_g; vector<double> tau_1_g;
  
  vector<double> ncons_ch; vector<double> ncons_neut;
  vector<vector<double> > conspt_ch; vector<vector<double> > conspt_neut;
  vector<vector<double> > conseta_ch; vector<vector<double> > conseta_neut;

  double eventID; double runID;
  //contains all (important) jet observables for both groomed and ungroomed jets
  TTree *eventTree = new TTree("event","event");
  eventTree->Branch("eventID",&eventID);
  eventTree->Branch("runID",&runID);
  eventTree->Branch("n_jets", &n_jets);
  eventTree->Branch("refMult", &refMult);//reference multiplicity
  eventTree->Branch("UE_avgpt", &UE_avgpt);
  eventTree->Branch("UEnpart", &UEnpart);
  eventTree->Branch("bbc_east_rate", &bbc_east_rate); //~lumi
  eventTree->Branch("bbc_east_sum", &bbc_east_sum); //~centrality
  eventTree->Branch("zdc_coinc_rate", &zdc_coinc_rate); //~lumi
  eventTree->Branch("ncons_ch", &ncons_ch);
  eventTree->Branch("ncons_neut", &ncons_neut);
  eventTree->Branch("conspt_ch", &conspt_ch);
  eventTree->Branch("conspt_neut", &conspt_neut);
  eventTree->Branch("conseta_ch", &conseta_ch);
  eventTree->Branch("conseta_neut", &conseta_neut);
  eventTree->Branch("Pt", &Pt); eventTree->Branch("Eta",&Eta); eventTree->Branch("Phi",&Phi); eventTree->Branch("M",&M); eventTree->Branch("E",&E);
  eventTree->Branch("ch_e_frac",&ch_e_frac);
  eventTree->Branch("zg", &zg); eventTree->Branch("rg", &rg); eventTree->Branch("mg", &mg); eventTree->Branch("ptg",&ptg);
  eventTree->Branch("mcd",&mcd);
  eventTree->Branch("tau0",&tau0); eventTree->Branch("tau05",&tau05); eventTree->Branch("tau_05",&tau_05); eventTree->Branch("tau_1",&tau_1);
  eventTree->Branch("tau0_g",&tau0_g); eventTree->Branch("tau05_g",&tau05_g); eventTree->Branch("tau_05_g",&tau_05_g); eventTree->Branch("tau_1_g",&tau_1_g);

  //TEMP check of the dN/d#eta of particles to compare to embedding for pA!
  vector<double> part_pt; vector<double> part_eta;
  TTree *particle_tree = new TTree("particle_tree","particle_tree");
  particle_tree->Branch("refMult", &refMult);
  particle_tree->Branch("bbc_east_sum", &bbc_east_sum);
  particle_tree->Branch("Pt",&part_pt);
  particle_tree->Branch("Eta",&part_eta);

  // Helpers
  // -------
  vector<PseudoJet> particles;
  
  unsigned nJets = 0, nJetsSD = 0;

  // Constituent selectors
  // ---------------------
  Selector select_track_eta  = fastjet::SelectorAbsEtaMax(max_track_eta);
  Selector select_pt_min     = fastjet::SelectorPtMin( partMinPt );
  Selector select_pt_max     = fastjet::SelectorPtMax( partMaxPt );
  Selector spart = select_track_eta * select_pt_min * select_pt_max;

  // Jet candidate selectors
  // -----------------------
  Selector select_jet_eta     = fastjet::SelectorAbsEtaMax(max_eta);
  Selector select_jet_pt_min  = fastjet::SelectorPtMin( det_jet_ptmin );
  Selector select_jet_pt_max  = fastjet::SelectorPtMax( jet_ptmax );
  Selector select_jet_m_min   = fastjet::SelectorMassMin( mass_min );
  Selector sjet = select_jet_eta && select_jet_pt_min && select_jet_pt_max && select_jet_m_min;
  
  // Choose a jet and area definition
  // --------------------------------
  JetDefinition jet_def = fastjet::JetDefinition(fastjet::antikt_algorithm, radius/*R*/);
  
  // create an area definition for the clustering
  //----------------------------------------------------------
  // ghosts should go up to the acceptance of the detector or
  // (with infinite acceptance) at least 2R beyond the region
  // where you plan to investigate jets.
  GhostedAreaSpec area_spec = fastjet::GhostedAreaSpec( ghost_maxrap, ghost_repeat, ghost_area );
  AreaDefinition  area_def = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts, area_spec);

  //Creating SoftDrop grooming object
  contrib::SoftDrop sd(Beta,z_cut,R0);
  cout << "SoftDrop groomer is: " << sd.description() << endl;
  
  //for later use looking up PDG masses using particle PID            
  TDatabasePDG *pdg = new TDatabasePDG();

  cout << "Performing analysis." << endl;
  // Cycle through events
  // --------------------  
  int nEventsUsed = 0;
  
  try{
    while ( reader->NextEvent() ) {
      
      //clearing vectors
      //temp!
      part_pt.clear(); part_eta.clear();
      
      Pt.clear(); Eta.clear(); Phi.clear(); M.clear(); E.clear();
      zg.clear(); rg.clear(); mg.clear(); ptg.clear();
      ch_e_frac.clear(); mcd.clear();
      tau0.clear(); tau05.clear(); tau_05.clear(); tau_1.clear(); //for now, angularity observables may be wrong. Don't trust until further vetting
      tau0_g.clear(); tau05_g.clear(); tau_05_g.clear(); tau_1_g.clear();
      //initializing variables to -9999
      n_jets = -9999; refMult = -9999; bbc_east_rate = -9999; bbc_east_sum = -9999; zdc_coinc_rate = -9999;
      eventID = -9999; runID = -9999;
      ncons_neut.clear(); ncons_ch.clear();
      conspt_neut.clear(); conspt_ch.clear();
      conseta_neut.clear(); conseta_ch.clear();

      UE_avgpt = -9999;
      UEnpart = -9999;
      
      reader->PrintStatus(10);
      
      //get the event header
      event = reader->GetEvent();
      header = event->GetHeader();

      particles.clear();
      
      // Get the output container from the reader
      // ----------------------------------------
      container = reader->GetOutputContainer();
      
      eventID = reader->GetNOfCurrentEvent();
      runID = header->GetRunId();
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~Skipping undesired events!~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                                                    

      //if the event lacks the desired trigger, skip it                                                                                                        
      //see function "SetTriggers()" for assignment of tID1, tID2 (or above until I write it)                                                                  
      if ( ! (header->HasTriggerId(tID1) || header->HasTriggerId(tID2) || header->HasTriggerId(tID3) ) ) {//cout << "DEBUG: skipping this event because it lacks appropriate triggers. Does it have trigger ID " << tID1 << "? " << header->HasTriggerId(tID1) << endl;
      	continue;
      }
      if (trigger.find("pA") != string::npos) {//removing some runs by hand in pA until we have bad run/tower lists
        //TEMPORARILY SKIPPING THESE RUNS for pA [should define these runIDs somewhere later so they're not magic numbers]                                     
        if (header->GetRunId() >= 16142059 && header->GetRunId() <= 16149001) {cout << "DEBUG: should never see this for pp!" << endl; continue;}
        //something weird happened to the towers in run 16135032 (and it looks like it started at the end of run 16135031), so excluding both                  
        if (header->GetRunId() == 16135031 || header->GetRunId() == 16135032) {cout << "DEBUG: should never see this for pp!" << endl; continue;}
      	//the event cuts don't check if the vzdiff is acceptable, so I have to hardcode it here. UPDATE: I believe Nick updated the eventstructuredAu to check this condition, so this line should be redundant now
      	//	if (!EventCuts->IsVertexZDiffOK(event)) {cout << "DEBUG: shouldn't see this now!" << endl; continue;}
      	//Above 64000 seems like detectors saturate (tower multiplicity explodes).                                                                             
        if (header->GetBbcAdcSumEast() >= pAu_BBCE_ADC_sum_max) {continue;}
      }
      
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~// 

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~STARTING ANALYSIS ON THE EVENT!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      
      // Transform TStarJetVectors into (FastJet) PseudoJets; assign mass to particles
      // ----------------------------------------------------------
      GatherParticles(container, sv, particles, full, 0, pdg); //"pdg" here finds the rest mass for particles with a given PID; 0 means detector-level.


      
      // Analysis
      // --------
      // Apply selector, spart (defined above, under "Constituent selectors"), to the full particle set
      vector<PseudoJet> good_parts = spart( particles );

      //calculating refMult (N_ch in |eta| < 0.5)                                                                                                              
      double refmult = 0;
      for (int i = 0; i < good_parts.size(); ++ i) {
      	//temporarily using a modified definition of refmult that includes neutrals to check something in the unfolding
        if (good_parts[i].user_index() != 0 && good_parts[i].user_index() != -9999 && fabs(good_parts[i].eta()) < 0.5) {
      	  refmult ++;
      	}
      }	
      //refmult = particles.size();
      
      // find corresponding jets with soft constituents
      // ----------------------------------------------
      ClusterSequence/*Area*/ csa ( good_parts, jet_def/*, area_def */); // WITHOUT background subtraction
	
      /*
      // Background initialization - for later use with AA
      // -------------------------
      // Background selector - Exclude two hardest jets for background extermination
      Selector selector_bkgd = fastjet::SelectorAbsRapMax( max_rap ) * (!fastjet::SelectorNHardest(2));
      // Area - same as for jets
      AreaDefinition area_def_bkgd ( area_def );
      // Jet definition - use kT instead of anti-kT algorithm here
      JetDefinition jet_def_bkgd (fastjet::kt_algorithm, R );
      // Energy density estimate from median ( pt_i / area_i )
      JetMedianBackgroundEstimator bkgd_estimator (selector_bkgd, jet_def_bkgd, area_def_bkgd);
      bkgd_estimator.set_particles( pLo );
      // Subtract A*rho from the original pT & m
      Subtractor bkgd_subtractor (&bkgd_estimator);
      bkgd_subtractor.set_use_rho_m();*/
      
      vector<PseudoJet> initial = fastjet::sorted_by_pt(sjet(/*bkgd_subtractor(*/csa.inclusive_jets())/*)*/);//applies jet selector to clustered jets
      vector<PseudoJet> good_jets; //jets in population "initial" then have NEF selection applied (since this is data)

      //Implementing a neutral energy fraction cut of 90% on inclusive jets
      ApplyNEFSelection(initial, good_jets);
      

      //temp: get UE                                                                                                                                           
      vector<fastjet::PseudoJet> UE;
      if (good_jets.size() > 0) {
        GatherUE (good_jets[0], container , UE );
      }
      double avgUEpt = -9999;
      if (UE.size() != 0) {
	avgUEpt = 0; // prevents filling with 0 entries for avg pT when no particles in the away-region.
        for (int i = 0; i < UE.size(); ++ i) {
          avgUEpt += UE[i].pt();
        }
        avgUEpt /= (double) UE.size();
      }
      UE_avgpt = avgUEpt;

      UEnpart = UE.size();


      vector<PseudoJet> GroomedJets;
      //loop over the jets which passed cuts, groom them, and add to a vector (sorted by pt of the original jet)
      for (int i = 0; i < good_jets.size(); ++ i) {
	GroomedJets.push_back(sd(good_jets[i])); //in this step, we Soft Drop groom the jet and add it to the GroomedJet population [whether or not z > 0.1]
	if (sd(good_jets[i]).structure_of<SD>().symmetry() >= 0.1) { //incrementing total N of groomed jets with z > 0.1 for debug later
	  nJetsSD ++;
	}
      }
      
      if (good_jets.size() != 0) {
	//SoftDrop is a groomer not a tagger, so if we have at least one ungroomed jet, we should also have a SoftDrop'd jet.
	nJets += good_jets.size();
	n_jets = good_jets.size();
	refMult = refmult;
	bbc_east_rate = header->GetBbcEastRate();
	bbc_east_sum = header->GetBbcAdcSumEast(); 
	zdc_coinc_rate = header->GetZdcCoincidenceRate();
	if (good_jets.size() >= 3) {
	  if (good_jets[0].pt() > 30 && good_jets[1].pt() > 25 && good_jets[2].pt() > 15) {
	    cout << "FOR RAGHAV! : " << std::setprecision(10) << eventID << " " << runID << endl;
	  }
	}
	//temp! particle loop for good events
	for (int i = 0; i < good_parts.size(); ++ i) {
	  part_pt.push_back(good_parts[i].pt());
	  part_eta.push_back(good_parts[i].eta());
	}

	for (int i = 0; i < n_jets; ++ i) {
	  Pt.push_back(good_jets[i].pt()); Eta.push_back(good_jets[i].eta()); Phi.push_back(good_jets[i].phi());
	  M.push_back(good_jets[i].m()); E.push_back(good_jets[i].e());
	  zg.push_back(GroomedJets[i].structure_of<SD>().symmetry()); rg.push_back(GroomedJets[i].structure_of<SD>().delta_R());
	  mg.push_back(GroomedJets[i].m()); ptg.push_back(GroomedJets[i].pt());
	  
	  //squaring mass and groomed mass for use in calculating M_cd [collinear dropped mass] which represents mass that was groomed away
	  double m2 = (good_jets[i].m())*(good_jets[i].m()); double gm2 = (GroomedJets[i].m())*(GroomedJets[i].m());
	  double m_cd = (double) sqrt(m2 - gm2); if ((m2 - gm2) < 1e-10) {m_cd = 0;}
	  mcd.push_back(m_cd);
	  
	  //for calculating angularities - still haven't vetted these, so use with 15 grains of salt
	  double ch_e = 0; double tot_e = 0;
	  double sum_tau0 = 0; double sum_tau05 = 0; double sum_tau_05 = 0; double sum_tau_1 = 0;
	  double sum_tau0_g = 0; double sum_tau05_g = 0; double sum_tau_05_g = 0; double sum_tau_1_g = 0;
	  vector<PseudoJet> cons = good_jets[i].constituents(); //ungroomed jet's constituents
	  vector<PseudoJet> cons_g = GroomedJets[i].constituents(); //groomed jet's constituents

	  int nconsch = 0, nconsneut = 0;
	  vector<double> consptch; vector<double> consptneut;
	  vector<double> consetach; vector<double> consetaneut;
	  for (int j = 0; j < cons.size(); ++ j) {
	    if (cons[j].user_index() == 0) { nconsneut ++; consptneut.push_back(cons[j].pt()); consetaneut.push_back(cons[j].eta());}
	    if ((cons[j].user_index() != 0) && fabs(cons[j].user_index()) < 5) { nconsch ++; consptch.push_back(cons[j].pt()); consetach.push_back(cons[j].eta());}
	  }
	  ncons_ch.push_back(nconsch); ncons_neut.push_back(nconsneut);
	  conspt_ch.push_back(consptch); conspt_neut.push_back(consptneut);
	  conseta_ch.push_back(consetach); conseta_neut.push_back(consetaneut);

	  //using pT, not energy, so some names are misnomers
	  for (int j = 0; j < cons.size(); ++ j) { //looping over ungroomed jet's constituents
	    if (cons[j].user_index() != 0) {ch_e += cons[j].pt();}
	    tot_e += cons[j].pt();
	    //angularity:
	    sum_tau0 += (cons[j].pt()*pow(good_jets[i].delta_R(cons[j]), 2 - 0));
	    sum_tau05 += (cons[j].pt()*pow(good_jets[i].delta_R(cons[j]), 2 - 0.5));
	    sum_tau_05 += (cons[j].pt()*pow(good_jets[i].delta_R(cons[j]), 2 + 0.5));
	    sum_tau_1 += (cons[j].pt()*pow(good_jets[i].delta_R(cons[j]), 2 + 1));
	  }//ungroomed jet
	  for (int j = 0; j < cons_g.size(); ++j) { //looping over groomed jet's constituents 
	    sum_tau0_g += (cons_g[j].pt()*pow(GroomedJets[i].delta_R(cons_g[j]), 2 - 0));
	    sum_tau05_g += (cons_g[j].pt()*pow(GroomedJets[i].delta_R(cons_g[j]), 2 - 0.5));
	    sum_tau_05_g += (cons_g[j].pt()*pow(GroomedJets[i].delta_R(cons_g[j]), 2 + 0.5));
	    sum_tau_1_g += (cons_g[j].pt()*pow(GroomedJets[i].delta_R(cons_g[j]), 2 + 1));
	  }//groomed jet
	  //angularity is scaled by the jet pT
	  tau0.push_back(sum_tau0 / (double) good_jets[i].pt());
	  tau05.push_back(sum_tau05 / (double) good_jets[i].pt());
	  tau_05.push_back(sum_tau_05 / (double) good_jets[i].pt());
	  tau_1.push_back(sum_tau_1 / (double) good_jets[i].pt());
	  
	  tau0_g.push_back(sum_tau0_g / (double) GroomedJets[i].pt());
	  tau05_g.push_back(sum_tau05_g / (double) GroomedJets[i].pt());
	  tau_05_g.push_back(sum_tau_05_g / (double) GroomedJets[i].pt());
	  tau_1_g.push_back(sum_tau_1_g / (double) GroomedJets[i].pt());

	  ch_e_frac.push_back(ch_e/(double)tot_e); //CEF [charged energy fraction] in the jet
	}//for loop over jets
      }//if we had jets
      
      // And we're done!
      // -----------------------------
      //temp, uncomment the if (good_jets ) later
      if (good_jets.size() != 0) {
	nEventsUsed++; //this event was accepted and had at least one jet passing all criteria
	eventTree->Fill();
	particle_tree->Fill();
      }
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
  eventTree->Write();
  particle_tree->Write();
  //fout->Write();
  fout->Close();

  cout << "In " << nEventsUsed << " events, found " << endl
       << nJets << " jets above 5 GeV, with constituents above 0.2 GeV," << endl
       << "for an average of " << nJets/(double)nEventsUsed << " jets per event" << endl;
  cout << "In those events, also found " << endl
       << nJetsSD << " groomed jets with z > 0.1" << endl
       << "for an average of " << nJetsSD/(double)nEventsUsed << " groomed jets per event" << endl;
  cout << "So we lose about " << (1-(nJetsSD/(double)nJets))*100 << "% of our jets to grooming by using it as a tagger" << endl;
  
  cout << "Wrote to " << fout->GetName() << endl;
  cout << "Bye :)" << endl;
    
  return 0;
}
  
