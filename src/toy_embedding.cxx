//Isaac Mooney, WSU - November 2019
//This file embeds ppJP2 data into pAuMB data

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
// [2]: full jets (Ch+Ne) or Ch jets
// [3]: flag determining if we run over ch+ne (i.e. "full") jets or just ch jets. If it's "ch", runs over ch jets. If anything else, runs over full jets.
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
  std::string chainList_pp = "list.txt"; // input file: can be .root, .txt, .list                                                                           
  std::string chainList_pAu = "list.txt"; // input file: can be .root, .txt, .list                                                                                 
  std::string chainName = "JetTree"; // tree name in input file                                                                                              
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
    if (arguments[3] == "ch") {full = 0;} else {full = 1;} //either ch+ne jets (default) or ch jets (if "ch")   
    chainList_pp      = arguments[4];
    chainList_pAu     = arguments[5];
            
    std::cout << "Running embedding of jets in the JP2-triggered pp data into the MB pAu data. Results will be output to " << outputDir << "." << std::endl;
    std::cout << "The input pp file is " << chainList_pp << ", the input pAu file is " << chainList_pAu << " and the output file is " << outFileName << "." << std::endl;
    
    break;
  }
  default: { // Error: invalid custom settings
    __ERR("Invalid number of command line arguments");
      return -1;
    break;
  }
  }

  //Setting up specifics of analysis based on the flags that were received above!
  string badtows_pp = "", badruns_pp = "", badtows_pAu = "", badruns_pAu = "";
  double vzdiff_pp = -1, vzdiff_pAu = -1;
  badtows_pp = det_badTowers; badruns_pp = dat_bad_run_list; vzdiff_pp = det_vZDiff;
  badtows_pAu = pAu_badTowers; badruns_pAu = pAu_bad_run_list; vzdiff_pAu = pAu_vZDiff;

  //in place for now; will encapsulate in a function if it gets much more involved. Hardcodes the trigger IDs.                                                 
  int tID_pp = -9999, tID1_pAu = -9999, tID2_pAu = -9999;
  tID_pp = tppJP2; //ppJP2
  tID1_pAu = tpAuBBCMBa; tID2_pAu = tpAuBBCMBb; //pAuBBCMB
  
  // Build our input now
  // --------------------
  TChain* chain_pp = new TChain( chainName.c_str() );
  TChain* chain_pAu = new TChain( chainName.c_str() );
  
  
  // Check to see if the input is a .root file or a .txt
  bool inputIsRoot_pp = Analysis::HasEnding( chainList_pp.c_str(), ".root" );
  bool inputIsTxt_pp  = Analysis::HasEnding( chainList_pp.c_str(), ".txt"  );
  bool inputIsList_pp = Analysis::HasEnding( chainList_pp.c_str(), ".list" );
  
  // If its a recognized file type, build the chain
  // If its not recognized, exit
  if ( inputIsRoot_pp ) { chain_pp->Add( chainList_pp.c_str() ); }
  else if ( inputIsTxt_pp )  { chain_pp = TStarJetPicoUtils::BuildChainFromFileList( chainList_pp.c_str() ); }
  else if ( inputIsList_pp )  { chain_pp = TStarJetPicoUtils::BuildChainFromFileList( chainList_pp.c_str() ); }
  else { __ERR("data file is not recognized type: .root or .txt only.") return -1; }

    // Check to see if the input is a .root file or a .txt
  bool inputIsRoot_pAu = Analysis::HasEnding( chainList_pAu.c_str(), ".root" );
  bool inputIsTxt_pAu  = Analysis::HasEnding( chainList_pAu.c_str(), ".txt"  );
  bool inputIsList_pAu = Analysis::HasEnding( chainList_pAu.c_str(), ".list" );
  
  // If its a recognized file type, build the chain
  // If its not recognized, exit
  if ( inputIsRoot_pAu ) { chain_pAu->Add( chainList_pAu.c_str() ); }
  else if ( inputIsTxt_pAu )  { chain_pAu = TStarJetPicoUtils::BuildChainFromFileList( chainList_pAu.c_str() ); }
  else if ( inputIsList_pAu )  { chain_pAu = TStarJetPicoUtils::BuildChainFromFileList( chainList_pAu.c_str() ); }
  else { __ERR("data file is not recognized type: .root or .txt only.") return -1; }


  // Build the event structure w/ cuts
  // ---------------------------------
  TStarJetPicoReader * reader_pp = new TStarJetPicoReader();
  InitReader(reader_pp, chain_pp, nEvents, "All"/*det_triggerString*/, det_absMaxVz, vzdiff_pp, det_evPtMax, det_evEtMax, det_evEtMin, det_DCA, det_NFitPts, det_FitOverMaxPts, dat_maxEtTow, 0.9999, false, badtows_pp, badruns_pp);
  TStarJetPicoReader * reader_pAu = new TStarJetPicoReader();
  InitReader(reader_pAu, chain_pAu, nEvents, "All", det_absMaxVz, vzdiff_pAu, det_evPtMax, det_evEtMax, det_evEtMin, det_DCA, det_NFitPts, det_FitOverMaxPts, dat_maxEtTow, 0.9999, false, badtows_pAu, badruns_pAu);

  // Data classes
  // ------------
  TStarJetVectorContainer<TStarJetVector>* container_pp;
  TStarJetVector* sv_pp; // TLorentzVector* would be sufficient
  TStarJetPicoEventHeader* header_pp;
  TStarJetPicoEvent* event_pp;
  TStarJetPicoEventCuts* EventCuts_pp = reader_pp->GetEventCuts();//will use this for hardcoding checks if events passed selections
  
  TStarJetVectorContainer<TStarJetVector>* container_pAu;
  TStarJetVector* sv_pAu; // TLorentzVector* would be sufficient
  TStarJetPicoEventHeader* header_pAu;
  TStarJetPicoEvent* event_pAu;
  TStarJetPicoEventCuts* EventCuts_pAu = reader_pAu->GetEventCuts();//will use this for hardcoding checks if events passed selections

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
  double n_jets;
  double bbc_east_rate_pp, bbc_east_sum_pp;
  double bbc_east_rate_pAu, bbc_east_sum_pAu;
  
  vector<double> Pt; vector<double> Eta; vector<double> Phi; vector<double> M; vector<double> E;
  vector<double> ch_e_frac;
  vector<double> zg; vector<double> rg; vector<double> mg; vector<double> ptg;
  vector<double> mcd;
  vector<double> tau0; vector<double> tau05; vector<double> tau_05; vector<double> tau_1;
  vector<double> tau0_g; vector<double> tau05_g; vector<double> tau_05_g; vector<double> tau_1_g;
  
  //contains all (important) jet observables for both groomed and ungroomed jets
  TTree *eventTree = new TTree("event","event");
  eventTree->Branch("n_jets", &n_jets);
  eventTree->Branch("bbc_east_rate_pp", &bbc_east_rate_pp); //~lumi
  eventTree->Branch("bbc_east_sum_pp", &bbc_east_sum_pp); //~centrality
  eventTree->Branch("bbc_east_rate_pAu", &bbc_east_rate_pAu); //~lumi
  eventTree->Branch("bbc_east_sum_pAu", &bbc_east_sum_pAu); //~centrality

  eventTree->Branch("Pt", &Pt); eventTree->Branch("Eta",&Eta); eventTree->Branch("Phi",&Phi); eventTree->Branch("M",&M); eventTree->Branch("E",&E);
  eventTree->Branch("ch_e_frac",&ch_e_frac);
  eventTree->Branch("zg", &zg); eventTree->Branch("rg", &rg); eventTree->Branch("mg", &mg); eventTree->Branch("ptg",&ptg);
  eventTree->Branch("mcd",&mcd);
  eventTree->Branch("tau0",&tau0); eventTree->Branch("tau05",&tau05); eventTree->Branch("tau_05",&tau_05); eventTree->Branch("tau_1",&tau_1);
  eventTree->Branch("tau0_g",&tau0_g); eventTree->Branch("tau05_g",&tau05_g); eventTree->Branch("tau_05_g",&tau_05_g); eventTree->Branch("tau_1_g",&tau_1_g);

  // Helpers
  // -------
  vector<PseudoJet> particles, particles_pp, particles_pAu;
  
  unsigned nJets = 0, nJetsSD = 0;

  // Constituent selectors
  // ---------------------
  Selector select_track_rap  = fastjet::SelectorAbsRapMax(max_track_rap);
  Selector select_pt_min     = fastjet::SelectorPtMin( partMinPt );
  Selector select_pt_max     = fastjet::SelectorPtMax( partMaxPt );
  Selector spart = select_track_rap * select_pt_min * select_pt_max;

  // Jet candidate selectors
  // -----------------------
  Selector select_jet_rap     = fastjet::SelectorAbsRapMax(max_rap);
  Selector select_jet_pt_min  = fastjet::SelectorPtMin( det_jet_ptmin );
  Selector select_jet_pt_max  = fastjet::SelectorPtMax( jet_ptmax );
  Selector select_jet_m_min   = fastjet::SelectorMassMin( mass_min );
  Selector sjet = select_jet_rap && select_jet_pt_min && select_jet_pt_max && select_jet_m_min;
  
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
    while ( reader_pp->NextEvent() && reader_pAu->NextEvent()) { //is this ok?
      
      //clearing vectors
      Pt.clear(); Eta.clear(); Phi.clear(); M.clear(); E.clear();
      zg.clear(); rg.clear(); mg.clear(); ptg.clear();
      ch_e_frac.clear(); mcd.clear();
      tau0.clear(); tau05.clear(); tau_05.clear(); tau_1.clear(); //for now, angularity observables may be wrong. Don't trust until further vetting
      tau0_g.clear(); tau05_g.clear(); tau_05_g.clear(); tau_1_g.clear();
      //initializing variables to -9999
      n_jets = -9999;
      bbc_east_rate_pp = -9999; bbc_east_sum_pp = -9999;
      bbc_east_rate_pAu = -9999; bbc_east_sum_pAu = -9999;

      reader_pp->PrintStatus(10);
      reader_pAu->PrintStatus(10);
      
      //get the event header
      event_pp = reader_pp->GetEvent();
      header_pp = event_pp->GetHeader();
      event_pAu = reader_pAu->GetEvent();
      header_pAu = event_pAu->GetHeader();

      particles.clear(); particles_pAu.clear(); particles_pp.clear();
      
      // Get the output container from the reader
      // ----------------------------------------
      container_pp = reader_pp->GetOutputContainer();
      container_pAu = reader_pAu->GetOutputContainer();
      
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~Skipping undesired events!~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                                                    

      //if the event lacks the desired trigger, skip it                                                                                                        
      //see function "SetTriggers()" for assignment of tID1, tID2 (or above until I write it)                                                                  
      if ( ! header_pp->HasTriggerId(tID_pp) ) {
	cout << "DEBUG: skipping this event because the pp lacks appropriate triggers." << endl;
	continue;
      }
      if ( ! (header_pAu->HasTriggerId(tID1_pAu) || header_pAu->HasTriggerId(tID2_pAu) ) ) {
	cout << "DEBUG: skipping this event because the pAu lacks appropriate triggers." << endl;
	continue;
      }
      //removing some runs by hand in pA until we have bad run/tower lists
      //TEMPORARILY SKIPPING THESE RUNS for pA [should define these runIDs somewhere later so they're not magic numbers]                                   
      if (header_pAu->GetRunId() >= 16142059 && header_pAu->GetRunId() <= 16149001) {cout << "skipping bad pAu runs." << endl; continue;}
      //something weird happened to the towers in run 16135032 (and it looks like it started at the end of run 16135031), so excluding both                  
      if (header_pAu->GetRunId() == 16135031 || header_pAu->GetRunId() == 16135032) {cout << "skipping bad pAu runs." << endl; continue;}
      //Above 64000 seems like detectors saturate (tower multiplicity explodes).                                                                             
      if (header_pAu->GetBbcAdcSumEast() >= pAu_BBCE_ADC_sum_max) {continue;}
      
      
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~// 

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~STARTING ANALYSIS ON THE EVENT!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      
      // Transform TStarJetVectors into (FastJet) PseudoJets; assign mass to particles
      // ----------------------------------------------------------
      GatherParticles(container_pp, sv_pp, particles_pp, full, 0, pdg); //"pdg" here finds the rest mass for particles with a given PID; 0 means detector-level.
      GatherParticles(container_pAu, sv_pAu, particles_pAu, full, 0, pdg); //"pdg" here finds the rest mass for particles with a given PID; 0 means detector-level.
      //adding both pp and pA particles to the vector of pseudojets containing the events' particles
      //particles.push_back(particles_pp); particles.push_back(particles_pAu);
      particles.insert(particles.end(), particles_pp.begin(), particles_pp.end());
      particles.insert(particles.end(), particles_pAu.begin(), particles_pAu.end());

      // Analysis
      // --------
      // Apply selector, spart (defined above, under "Constituent selectors"), to the full particle set
      vector<PseudoJet> good_parts = spart( particles );
      
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
	bbc_east_rate_pp = header_pp->GetBbcEastRate();
	bbc_east_sum_pp = header_pp->GetBbcAdcSumEast();
	bbc_east_rate_pAu = header_pAu->GetBbcEastRate();
	bbc_east_sum_pAu = header_pAu->GetBbcAdcSumEast(); 
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
      if (good_jets.size() != 0) {
	nEventsUsed++; //this event was accepted and had at least one jet passing all criteria
	eventTree->Fill();
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
  
