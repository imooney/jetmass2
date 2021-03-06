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
  std::string chainNameMc = "JetTreeMc";
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
  badtows_pp = /*det_badTowers*/combined_badTowers; badruns_pp = dat_bad_run_list; vzdiff_pp = det_vZDiff;
  badtows_pAu = /*pAu_badTowers*/combined_badTowers; badruns_pAu = pAu_bad_run_list; vzdiff_pAu = pAu_vZDiff;

  //TEMP:
  //badtows_pp = sim_badTowers;
  //badtows_pAu = sim_badTowers;

  //badruns_pp = sim_bad_run_list;
  //badruns_pAu = sim_bad_run_list;
  
  //in place for now; will encapsulate in a function if it gets much more involved. Hardcodes the trigger IDs.                                                 
  int tID_pp = -9999, tID1_pAu = -9999, tID2_pAu = -9999;
  tID_pp = tppJP2; //ppJP2
  tID1_pAu = tpAuBBCMBa; tID2_pAu = tpAuBBCMBb; //pAuBBCMB
  
  // Build our input now
  // --------------------
  TChain* chain_pp = new TChain( chainName.c_str() );
  TChain* chain_pAu = new TChain( chainName.c_str() );
  TChain* chain_pyth = new TChain( chainNameMc.c_str() );
  
  // Check to see if the input is a .root file or a .txt
  bool inputIsRoot_pp = Analysis::HasEnding( chainList_pp.c_str(), ".root" );
  bool inputIsTxt_pp  = Analysis::HasEnding( chainList_pp.c_str(), ".txt"  );
  bool inputIsList_pp = Analysis::HasEnding( chainList_pp.c_str(), ".list" );
  
  // If its a recognized file type, build the chain
  // If its not recognized, exit
  if ( inputIsRoot_pp ) { chain_pp->Add( chainList_pp.c_str()); }
  else if ( inputIsTxt_pp )  { chain_pp = TStarJetPicoUtils::BuildChainFromFileList( chainList_pp.c_str(), chainName.c_str() ); }
  else if ( inputIsList_pp )  { chain_pp = TStarJetPicoUtils::BuildChainFromFileList( chainList_pp.c_str(), chainName.c_str() ); }
  else { __ERR("data file is not recognized type: .root or .txt only.") return -1; }

  std::string chainList_pyth = chainList_pp;

  // If its a recognized file type, build the chain                                                                                                            
  // If its not recognized, exit                                                                                                                               
  if ( inputIsRoot_pp ) { chain_pyth->Add( chainList_pyth.c_str() ); }
  else if ( inputIsTxt_pp )  { chain_pyth = TStarJetPicoUtils::BuildChainFromFileList( chainList_pyth.c_str(), chainNameMc.c_str() ); }
  else if ( inputIsList_pp )  { chain_pyth = TStarJetPicoUtils::BuildChainFromFileList( chainList_pyth.c_str(), chainNameMc.c_str() ); }
  else { __ERR("data file is not recognized type: .root or .txt only.") return -1; }

    // Check to see if the input is a .root file or a .txt
  bool inputIsRoot_pAu = Analysis::HasEnding( chainList_pAu.c_str(), ".root" );
  bool inputIsTxt_pAu  = Analysis::HasEnding( chainList_pAu.c_str(), ".txt"  );
  bool inputIsList_pAu = Analysis::HasEnding( chainList_pAu.c_str(), ".list" );
  
  // If its a recognized file type, build the chain
  // If its not recognized, exit
  if ( inputIsRoot_pAu ) { chain_pAu->Add( chainList_pAu.c_str() ); }
  else if ( inputIsTxt_pAu )  { chain_pAu = TStarJetPicoUtils::BuildChainFromFileList( chainList_pAu.c_str(), chainName.c_str() ); }
  else if ( inputIsList_pAu )  { chain_pAu = TStarJetPicoUtils::BuildChainFromFileList( chainList_pAu.c_str(), chainName.c_str() ); }
  else { __ERR("data file is not recognized type: .root or .txt only.") return -1; }


  TFile *frefmultweight = new TFile("~/jetmass2/out/sim/hists/pA_peripheral_refmultbbcweight.root","READ");
  //  TH2D *hrefmultbbcweight = (TH2D*) frefmultbbcweight->Get("refmult_bbce_weight");
  TH1D* hrefmultrat = (TH1D*) frefmultweight->Get("refrat");
  hrefmultrat->SetDirectory(0);
  frefmultweight->Close();
 
  
  TFile *fbbcweight = new TFile("~/jetmass2/out/data/pAJP2_hbbcweight.root","READ");
  TH1D* hbbcrat = (TH1D*) fbbcweight->Get("bbcrat");
  //hrefmultbbcweight->SetDirectory(0);
  hbbcrat->SetDirectory(0);
  fbbcweight->Close();

  // Build the event structure w/ cuts
  // ---------------------------------
  TStarJetPicoReader * reader_pyth = new TStarJetPicoReader();
  InitReader(reader_pyth, chain_pyth, nEvents, "All", truth_absMaxVz, truth_vZDiff, truth_evPtMax, truth_evEtMax, truth_evEtMin, truth_DCA, truth_NFitPts, truth_FitOverMaxPts, sim_maxEtTow, 0.9999, false, sim_badTowers, sim_bad_run_list);
  TStarJetPicoReader * reader_pp = new TStarJetPicoReader();
  InitReader(reader_pp, chain_pp, nEvents, "All"/*det_triggerString*/, det_absMaxVz, vzdiff_pp, det_evPtMax, det_evEtMax, det_evEtMin, det_DCA, det_NFitPts, det_FitOverMaxPts, dat_maxEtTow, 0.9999, false, badtows_pp, badruns_pp);
  TStarJetPicoReader * reader_pAu = new TStarJetPicoReader();
  InitReader(reader_pAu, chain_pAu, nEvents, "All", det_absMaxVz, vzdiff_pAu, det_evPtMax, det_evEtMax, det_evEtMin, det_DCA, det_NFitPts, det_FitOverMaxPts, dat_maxEtTow, 0.9999, false, badtows_pAu, badruns_pAu);

  // Data classes
  // ------------
  TStarJetVectorContainer<TStarJetVector>* container_pyth;
  TStarJetVector* sv_pyth; // TLorentzVector* would be sufficient
  TString filename_pyth;
  TStarJetPicoEventHeader* header_pyth;
  TStarJetPicoEvent* event_pyth;
  TStarJetPicoEventCuts* EventCuts_pyth = reader_pyth->GetEventCuts();//will use this for hardcoding checks if events passed selections
  
  TStarJetVectorContainer<TStarJetVector>* container_pp;
  TStarJetVector* sv_pp; // TLorentzVector* would be sufficient
  TString filename_pp;
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
  double refMult_pp, refMult_pAu;
  double bbc_east_rate_pp, bbc_east_sum_pp;
  double bbc_east_rate_pAu, bbc_east_sum_pAu;
  
  double zdc_coinc_rate_pp, zdc_coinc_rate_pAu;
  
  double bkgrho, bkgrho_m, bkgsigma, bkgsigma_m;
  
  double pp_weight;

  vector<double> ncons_ch; vector<double> ncons_neut;
  vector<vector<double> > conspt_ch; vector<vector<double> > conspt_neut;
  vector<vector<double> > conseta_ch; vector<vector<double> > conseta_neut;

  vector<double> Pt_pyth; vector<double> Eta_pyth; vector<double> Phi_pyth; vector<double> M_pyth; vector<double> E_pyth;
  vector<double> ch_e_frac_pyth;
  vector<double> zg_pyth; vector<double> rg_pyth; vector<double> mg_pyth; vector<double> ptg_pyth;

  vector<double> Pt; vector<double> Eta; vector<double> Phi; vector<double> M; vector<double> E;
  vector<double> ch_e_frac;
  vector<double> zg; vector<double> rg; vector<double> mg; vector<double> ptg;
  vector<double> mcd;
  vector<double> tau0; vector<double> tau05; vector<double> tau_05; vector<double> tau_1;
  vector<double> tau0_g; vector<double> tau05_g; vector<double> tau_05_g; vector<double> tau_1_g;
  
  //contains all (important) jet observables for both groomed and ungroomed jets
  TTree *eventTree = new TTree("event","event");
  eventTree->Branch("n_jets", &n_jets);
  eventTree->Branch("refMult_pp", &refMult_pp);
  eventTree->Branch("refMult_pAu", &refMult_pAu);
  eventTree->Branch("bbc_east_rate_pp", &bbc_east_rate_pp); //~lumi
  eventTree->Branch("bbc_east_sum_pp", &bbc_east_sum_pp); //~centrality
  eventTree->Branch("bbc_east_rate_pAu", &bbc_east_rate_pAu); //~lumi
  eventTree->Branch("bbc_east_sum_pAu", &bbc_east_sum_pAu); //~centrality
  eventTree->Branch("zdc_coinc_rate_pp",&zdc_coinc_rate_pp);
  eventTree->Branch("zdc_coinc_rate_pAu",&zdc_coinc_rate_pAu); //to match pp to pA lumi
  eventTree->Branch("rho",&bkgrho);
  eventTree->Branch("rho_m",&bkgrho_m);
  eventTree->Branch("sigma",&bkgsigma);
  eventTree->Branch("sigma_m",&bkgsigma_m);
  eventTree->Branch("pp_weight",&pp_weight);

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
  eventTree->Branch("Pt_pyth", &Pt_pyth); eventTree->Branch("Eta_pyth",&Eta_pyth); eventTree->Branch("Phi_pyth",&Phi_pyth); eventTree->Branch("M_pyth",&M_pyth); eventTree->Branch("E_pyth",&E_pyth);
  eventTree->Branch("ch_e_frac_pyth",&ch_e_frac_pyth);
  eventTree->Branch("zg_pyth", &zg_pyth); eventTree->Branch("rg_pyth", &rg_pyth); eventTree->Branch("mg_pyth", &mg_pyth); eventTree->Branch("ptg_pyth",&ptg_pyth);
 
  eventTree->Branch("tau0",&tau0); eventTree->Branch("tau05",&tau05); eventTree->Branch("tau_05",&tau_05); eventTree->Branch("tau_1",&tau_1);
  eventTree->Branch("tau0_g",&tau0_g); eventTree->Branch("tau05_g",&tau05_g); eventTree->Branch("tau_05_g",&tau_05_g); eventTree->Branch("tau_1_g",&tau_1_g);

  //hists
  TH2D* respT_v_partpT = new TH2D("respT_v_partpT","",25,0,2,15,5,80);
  TH2D* resM_v_partpT = new TH2D("resM_v_partpT","",25,0,2,15,5,80);
  
  TH1D* truth_unweight = new TH1D("truth_unweight","",75,5,80);
  TH1D* truth_bbcweight = new TH1D("truth_bbcweight","",75,5,80);
  TH1D* truth_refweight = new TH1D("truth_refweight","",75,5,80);
  TH1D* truth_weight = new TH1D("truth_weight","",75,5,80);
 
  TH1D* truth_matches_unweight = new TH1D("truth_matches_unweight","",75,5,80);
  TH1D* truth_matches_bbcweight = new TH1D("truth_matches_bbcweight","",75,5,80);
  TH1D* truth_matches_refweight = new TH1D("truth_matches_refweight","",75,5,80);
  TH1D* truth_matches_weight = new TH1D("truth_matches_weight","",75,5,80);
 
  TH2D* weight_v_pt = new TH2D("weight_v_pt","",100,0,10,75,5,80);
  
  //hists for use in responses - note:                                                                                                                         
  //"hist_measured and hist_truth are used to specify the dimensions of the distributions (the histogram contents are not used here), eg. for 2D or 3D distributions or non-uniform binning." - http://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html                                                                   
  //I.e. these histograms are just used as a template, they're not actually what is being filled when the responses are constructed.                            
  TH2D *pyMvPt = new TH2D("pyMvPt",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMvPt = new TH2D("geMvPt",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMgvPt = new TH2D("pyMgvPt", ";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMgvPt = new TH2D("geMgvPt", ";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);

  // 2D responses                                                                                                                                              
  RooUnfoldResponse *m_pt_response = new RooUnfoldResponse(geMvPt, pyMvPt, "m_pt_response");
  RooUnfoldResponse *mg_pt_response = new RooUnfoldResponse(geMgvPt, pyMgvPt, "mg_pt_response");

  // Helpers
  // -------
  vector<PseudoJet> particles, particles_pyth, particles_pp, particles_pAu;
  
  unsigned nJets = 0, nJetsSD = 0, nJets_pyth = 0, nJetsSD_pyth = 0;

  unsigned count_evts = 0, count_jetevts = 0;

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
  //gen-level jets:
  Selector select_gen_jet_pt_min = fastjet::SelectorPtMin( jet_ptmin );
  Selector select_gen_jet_m_min = fastjet::SelectorMassMin( 0.0 );
  Selector sjet_gen = select_jet_rap && select_gen_jet_pt_min && select_jet_pt_max && select_gen_jet_m_min;

  // Choose a jet and area definition
  // --------------------------------
  JetDefinition jet_def = fastjet::JetDefinition(fastjet::antikt_algorithm, radius/*R*/);
  
  // create an area definition for the clustering
  //----------------------------------------------------------
  // ghosts should go up to the acceptance of the detector or
  // (with infinite acceptance) at least 2R beyond the region
  // where you plan to investigate jets.
  //GhostedAreaSpec area_spec = fastjet::GhostedAreaSpec( ghost_maxrap, ghost_repeat, ghost_area );
  //AreaDefinition  area_def = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts, area_spec);

  //Creating SoftDrop grooming object
  contrib::SoftDrop sd(Beta,z_cut,R0);
  cout << "SoftDrop groomer is: " << sd.description() << endl;
  
  //for later use looking up PDG masses using particle PID            
  TDatabasePDG *pdg = new TDatabasePDG();

  cout << "Performing analysis." << endl;
  // Cycle through events
  // --------------------  
  int nEventsUsed = 0;
  
  cout << "DEBUG: nentries_pp: " << chain_pp->GetEntries() << " nentries_pA: " << chain_pAu->GetEntries() << endl;
  cout << "DEBUG: nentries_pyth: " << chain_pyth->GetEntries() << endl;

  cout << "SHOULD BE FIXED!" << endl;

  unsigned npp = 0; unsigned npA = 0; unsigned npp_acc = 0; unsigned npA_acc = 0;

  try{
    //we check the pp event first. If it's bad, we just go to the next one. Once we find a good one, we check for a good pA event as well.
    //while ( reader_pp->NextEvent()) {
    for (int event = 0; event < chain_pyth->GetEntries(); ++ event) {
      reader_pyth->ReadEvent(event);

      reader_pp->ReadEvent(event);
      if ( reader_pp->ReadEvent(reader_pyth->GetNOfCurrentEvent()) != 1 ) {
	//cout << "no corresponding geant event...skipping event " << reader_pyth->GetNOfCurrentEvent() << endl;
        continue;//goes to the next pythia event 
      }
      if (reader_pyth->GetNOfCurrentEvent() != reader_pp->GetNOfCurrentEvent() ) { cerr << "ERROR: READING DIFFERENT EVENTS. EXITING." << endl; exit(1);}

      //after the above, should have a valid Pythia and corresponding Pythia+Geant+ZB2012 event

      //cout << "DEBUG: pp event #" << npp++ << endl;
      pp_weight = -9999;
      
      reader_pp->PrintStatus(10);
      //get the event header
      event_pp = reader_pp->GetEvent();
      header_pp = event_pp->GetHeader();
     
      // Get the output container from the reader
      // ----------------------------------------
      container_pp = reader_pp->GetOutputContainer();
      filename_pp =  reader_pp->GetInputChain()->GetCurrentFile()->GetName();
      pp_weight = LookupRun12Xsec( filename_pp );
      
      container_pyth = reader_pyth->GetOutputContainer();
      filename_pyth =  reader_pyth->GetInputChain()->GetCurrentFile()->GetName(); // should be the same as filename_pp
      
      //if the event lacks the desired trigger, skip it                                                                                                        
      //see function "SetTriggers()" for assignment of tID1, tID2 (or above until I write it)                                                                  
      if ( ! header_pp->HasTriggerId(tID_pp) ) {
	//cout << "DEBUG: skipping this event because the pp lacks appropriate triggers." << endl;
	continue;
      }

      //cout << "DEBUG: pp event #" << npp << " becomes accepted evt #" << npp_acc++ << endl;
	while ( reader_pAu->NextEvent() ) {
	//cout << "DEBUG: pA event #" << npA++ << endl;
	reader_pAu->PrintStatus(10);
      
      //get the event header
      event_pAu = reader_pAu->GetEvent();
      header_pAu = event_pAu->GetHeader();

      // Get the output container from the reader
      // ----------------------------------------
      container_pAu = reader_pAu->GetOutputContainer();
      
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~Skipping undesired events!~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                                                    

      //if the event lacks the desired trigger, skip it                                                                                                        
      //see function "SetTriggers()" for assignment of tID1, tID2 (or above until I write it)                                                                  
      if ( ! (header_pAu->HasTriggerId(tID1_pAu) || header_pAu->HasTriggerId(tID2_pAu) ) ) {
	//	cout << "DEBUG: skipping this event because the pAu lacks appropriate triggers." << endl;
	continue;
      }
      //removing some runs by hand in pA until we have bad run/tower lists
      //TEMPORARILY SKIPPING THESE RUNS for pA [should define these runIDs somewhere later so they're not magic numbers]                                   
      if (header_pAu->GetRunId() >= 16142059 && header_pAu->GetRunId() <= 16149001) {/*cout << "skipping bad pAu runs." << endl;*/ continue;}
      //something weird happened to the towers in run 16135032 (and it looks like it started at the end of run 16135031), so excluding both                  
      if (header_pAu->GetRunId() == 16135031 || header_pAu->GetRunId() == 16135032) {/*cout << "skipping bad pAu runs." << endl;*/ continue;}
      //Above 64000 seems like detectors saturate (tower multiplicity explodes).                                                                             
      if (header_pAu->GetBbcAdcSumEast() >= pAu_BBCE_ADC_sum_max) {/*cout << "DEBUG: BBCE ADC sum is too high! Skipping this bad pA event!" << endl;*/ continue;}


      //This cut loosely matches luminosity between pp and pA events:
      if (fabs(header_pAu->GetZdcCoincidenceRate() - header_pp->GetZdcCoincidenceRate()) > 10000) {
	continue; //lumis were too different; try the next pA event instead -- should have plenty to spare
      }
      /*
      //these cuts make the pA look more like pp:
      if (header_pAu->GetZdcCoincidenceRate() > 30000) {
	continue;//lumi is too high
      }
      */
      if ((header_pAu->GetBbcAdcSumEast() + header_pp->GetBbcAdcSumEast()) > 20000) {
	continue;//total "centrality" is too low (high) if < (>)
      }
      
      //cout << "DEBUG: this pA event is ok (unless pT[0] > 2*pt-hat-bin-upper-edge)! That means we have a pp-pA match!" << endl;
      //so clear vectors before we fill them:                    
      Pt.clear(); Eta.clear(); Phi.clear(); M.clear(); E.clear();
      zg.clear(); rg.clear(); mg.clear(); ptg.clear();
      ch_e_frac.clear(); mcd.clear();
      Pt_pyth.clear(); Eta_pyth.clear(); Phi_pyth.clear(); M_pyth.clear(); E_pyth.clear();
      zg_pyth.clear(); rg_pyth.clear(); mg_pyth.clear();  ptg_pyth.clear();
      ch_e_frac_pyth.clear();
      
      tau0.clear(); tau05.clear(); tau_05.clear(); tau_1.clear(); //for now, angularity observables may be wrong. Don't trust until further vetting
      tau0_g.clear(); tau05_g.clear(); tau_05_g.clear(); tau_1_g.clear();
      //initializing variables to -9999
      n_jets = -9999;
      refMult_pp = -9999; refMult_pAu = -9999;
      bbc_east_rate_pp = -9999; bbc_east_sum_pp = -9999;
      bbc_east_rate_pAu = -9999; bbc_east_sum_pAu = -9999;
      zdc_coinc_rate_pp = -9999; zdc_coinc_rate_pAu = -9999;
      bkgrho = -9999; bkgrho_m = -9999; bkgsigma = -9999; bkgsigma_m = -9999;
  
      ncons_neut.clear(); ncons_ch.clear();
      conspt_neut.clear(); conspt_ch.clear();
      conseta_neut.clear(); conseta_ch.clear();

      particles.clear(); particles_pyth.clear(); particles_pAu.clear(); particles_pp.clear();

      
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~// 

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~STARTING ANALYSIS ON THE EVENT!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      
      // Transform TStarJetVectors into (FastJet) PseudoJets; assign mass to particles
      // ----------------------------------------------------------
      GatherParticles(container_pyth, sv_pyth, particles_pyth, full, 1, pdg);
      GatherParticles(container_pp, sv_pp, particles_pp, full, 0, pdg); //"pdg" here finds the rest mass for particles with a given PID; 0 means detector-level.
      GatherParticles(container_pAu, sv_pAu, particles_pAu, full, 0, pdg); //"pdg" here finds the rest mass for particles with a given PID; 0 means detector-level.
      //adding both pp and pA particles to the vector of pseudojets containing the events' particles
      //particles.push_back(particles_pp); particles.push_back(particles_pAu);
      particles.insert(particles.end(), particles_pp.begin(), particles_pp.end());
      particles.insert(particles.end(), particles_pAu.begin(), particles_pAu.end());
      
      //calculating refMult (N_ch in |eta| < 0.5)
      double refmultpp = 0; double refmultpA = 0;
      for (int i = 0; i < particles_pp.size(); ++ i) {
	//testing a modified refmult that includes neutral particles as well
	if (particles_pp[i].user_index() != 0 && particles_pp[i].user_index() != -9999 && fabs(particles_pp[i].eta()) < 0.5) {
	  refmultpp ++;
	}
      }
      for (int i = 0; i < particles_pAu.size(); ++ i) {
	if (particles_pAu[i].user_index() != 0 && particles_pAu[i].user_index() != -9999 && fabs(particles_pAu[i].eta()) < 0.5) {
	  refmultpA ++;
	}
      }

      // refmultpp = particles_pp.size(); refmultpA = particles_pAu.size();
      // Analysis
      // --------
      // Apply selector, spart (defined above, under "Constituent selectors"), to the full particle set
      vector<PseudoJet> good_parts = spart( particles );
      vector<PseudoJet> good_parts_pyth = spart( particles_pyth );
      vector<PseudoJet> good_parts_MB = spart( particles_pAu  );

      // find corresponding jets with soft constituents
      // ----------------------------------------------
      ClusterSequence/*Area*/ csa ( good_parts, jet_def/*, area_def */); // WITHOUT background subtraction
      ClusterSequence csa_pyth ( good_parts_pyth, jet_def);
      ClusterSequence csa_MB ( good_parts_MB, jet_def);
      /*    
      // Background initialization - for later use with AA
      // -------------------------
      // Background selector - Exclude two hardest jets for background extermination
      Selector selector_bkgd = fastjet::SelectorAbsEtaMax( max_eta ) * (!fastjet::SelectorNHardest(2));
      // Area - same as for jets
      AreaDefinition area_def_bkgd ( area_def );
      // Jet definition - use kT instead of anti-kT algorithm here
      JetDefinition jet_def_bkgd (fastjet::kt_algorithm, R );
      // Energy density estimate from median ( pt_i / area_i )
      JetMedianBackgroundEstimator bkgd_estimator (selector_bkgd, jet_def_bkgd, area_def_bkgd);
      bkgd_estimator.set_particles( good_parts );
      // Subtract A*rho from the original pT & m
      Subtractor bkgd_subtractor (&bkgd_estimator);
      bkgd_subtractor.set_use_rho_m();
      
      //for plotting rho later:
      bkgrho = bkgd_estimator.rho();
      bkgsigma = bkgd_estimator.sigma();
      bkgrho_m = bkgd_estimator.rho_m();
      bkgsigma_m = bkgd_estimator.sigma_m();
*/
      
      vector<PseudoJet> initial = fastjet::sorted_by_pt(sjet(/*bkgd_subtractor(*/csa.inclusive_jets())/*)*/);//applies jet selector to clustered jets
      vector<PseudoJet> good_jets; //jets in population "initial" then have NEF selection applied (since this is data)
        
      vector<PseudoJet> initial_pyth = fastjet::sorted_by_pt(sjet_gen(csa_pyth.inclusive_jets()));//applies jet selector to clustered jets
      vector<PseudoJet> good_jets_pyth;

      vector<PseudoJet> initial_MB = fastjet::sorted_by_pt(sjet(csa_MB.inclusive_jets()));//applies jet selector to clustered jets
      vector<PseudoJet> good_jets_MB;


      //Implementing a neutral energy fraction cut of 90% on inclusive jets
      ApplyNEFSelection(initial, good_jets);
      good_jets_pyth = initial_pyth;
      good_jets_MB = initial_MB;
      /*      
      if (good_jets.size() != 0) {
	cout << "DEBUG: lead jet pt = " << good_jets[0].pt() << endl; 
      }
      */
      vector<PseudoJet> dummy;
      //used to do this for pp+pA but now that I pulled out the Pythia alone, can just ask for jets from both PYTHIA & PYTHIA+Ge+MB to be below the threshold
      if (DiscardEvent(filename_pp, good_jets_pyth, good_jets)) { break; } //jet pT is too high, will be overly weighted on the steeply falling spectrum. Try the next pp event since the pT is mostly due to the hard event from pp not the MB from pA
       
      //count_evts ++;
      //if (good_jets_MB.size() > 0) { if (good_jets_MB[0].pt() > 3.5) {count_jetevts ++; continue;}} //MB event has a semi-high-Q^2 process, so move on to next event, as a test

      //cout << "DEBUG: pA event #" << npA << " becomes accepted evt #" << npA_acc++ << endl;

      vector<PseudoJet> GroomedJets;
      //loop over the jets which passed cuts, groom them, and add to a vector (sorted by pt of the original jet)
      for (int i = 0; i < good_jets.size(); ++ i) {
	GroomedJets.push_back(sd(good_jets[i])); //in this step, we Soft Drop groom the jet and add it to the GroomedJet population [whether or not z > 0.1]
	if (sd(good_jets[i]).structure_of<SD>().symmetry() >= 0.1) { //incrementing total N of groomed jets with z > 0.1 for debug later
	  nJetsSD ++;
	}
      }

      vector<PseudoJet> GroomedJets_pyth;
      //loop over the jets which passed cuts, groom them, and add to a vector (sorted by pt of the original jet)                                               
      for (int i = 0; i < good_jets_pyth.size(); ++ i) {
        GroomedJets_pyth.push_back(sd(good_jets_pyth[i])); //in this step, we Soft Drop groom the jet and add it to the GroomedJet population [whether or not z > 0.1]   
        if (sd(good_jets_pyth[i]).structure_of<SD>().symmetry() >= 0.1) { //incrementing total N of groomed jets with z > 0.1 for debug later
          nJetsSD_pyth ++;
        }
      }
      
      //filling trees with GEANT+pAMB info for those events that had an embedding
      if (good_jets.size() != 0) {
	//cout << "DEBUG: ...and we have jets! Pt of 0th jet = " << good_jets[0].pt() << endl;
	//SoftDrop is a groomer not a tagger, so if we have at least one ungroomed jet, we should also have a SoftDrop'd jet.
	nJets += good_jets.size();
	n_jets = good_jets.size();
	refMult_pp = refmultpp;
	refMult_pAu = refmultpA;
	bbc_east_rate_pp = header_pp->GetBbcEastRate();
	bbc_east_sum_pp = header_pp->GetBbcAdcSumEast();
	bbc_east_rate_pAu = header_pAu->GetBbcEastRate();
	bbc_east_sum_pAu = header_pAu->GetBbcAdcSumEast(); 
	zdc_coinc_rate_pp = header_pp->GetZdcCoincidenceRate();
	zdc_coinc_rate_pAu = header_pAu->GetZdcCoincidenceRate();
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

      //filling trees with PYTHIA info for those events that had an embedding
      if (good_jets_pyth.size() != 0) {
	//cout << "DEBUG: ...and we have jets! Pt of 0th jet = " << good_jets[0].pt() << endl;
	//SoftDrop is a groomer not a tagger, so if we have at least one ungroomed jet, we should also have a SoftDrop'd jet.
	nJets_pyth += good_jets.size();
	//n_jets = good_jets_pyth.size();
	for (int i = 0; i < good_jets_pyth.size(); ++ i) {
	  Pt_pyth.push_back(good_jets_pyth[i].pt()); Eta_pyth.push_back(good_jets_pyth[i].eta()); Phi_pyth.push_back(good_jets_pyth[i].phi());
	  M_pyth.push_back(good_jets_pyth[i].m()); E_pyth.push_back(good_jets_pyth[i].e());
	  zg_pyth.push_back(GroomedJets_pyth[i].structure_of<SD>().symmetry()); rg_pyth.push_back(GroomedJets_pyth[i].structure_of<SD>().delta_R());
	  mg_pyth.push_back(GroomedJets_pyth[i].m()); ptg_pyth.push_back(GroomedJets_pyth[i].pt());
	  
	  //for calculating angularities - still haven't vetted these, so use with 15 grains of salt
	  double ch_e = 0; double tot_e = 0;
	  vector<PseudoJet> cons = good_jets_pyth[i].constituents(); //ungroomed jet's constituents
	  vector<PseudoJet> cons_g = GroomedJets_pyth[i].constituents(); //groomed jet's constituents

	  //using pT, not energy, so some names are misnomers
	  for (int j = 0; j < cons.size(); ++ j) { //looping over ungroomed jet's constituents
	    if (cons[j].user_index() != 0) {ch_e += cons[j].pt();}
	    tot_e += cons[j].pt();
	  }//ungroomed jet
	  ch_e_frac_pyth.push_back(ch_e/(double)tot_e); //CEF [charged energy fraction] in the jet
	}//for loop over jets
      }//if we had jets
      


      //embedded event has been accepted so let's geometrically match PYTHIA to PYTHIA+GEANT+ppZB2012+pAMB2015
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MATCHING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      //leaving the nomenclature from before this code segment was ported (laziness). So "g_" means not Geant anymore but Geant+pA_MB. "p_" still means just pythia.
      std::vector<fastjet::PseudoJet> g_matches; std::vector<fastjet::PseudoJet> p_matches;
      std::vector<fastjet::PseudoJet> g_matches_for_fakes; std::vector<fastjet::PseudoJet> p_matches_for_fakes; //only used to determine fakes               
      std::vector<fastjet::PseudoJet> fakes; std::vector<fastjet::PseudoJet> misses;
      std::vector<fastjet::PseudoJet> g_sd_matches; std::vector<fastjet::PseudoJet> p_sd_matches;
      std::vector<fastjet::PseudoJet> sd_fakes; std::vector<fastjet::PseudoJet> sd_misses;
      std::vector<int> match_indices;
      std::vector<int> miss_indices; std::vector<int> fake_indices;

      //matches & misses                                                                                                                                     
      if (good_jets_pyth.size() != 0) {
	g_matches.clear(); p_matches.clear(); misses.clear();
	match_indices.clear(); miss_indices.clear();
	match_indices = MatchJets(good_jets,good_jets_pyth, g_matches, p_matches); //find matches
	if (g_matches.size() != p_matches.size()) {std::cerr << "Somehow we have different-sized match vectors. Exiting!" <<std::endl; exit(1);}
	if (g_matches.size() < good_jets_pyth.size()) { //then we have misses
	  miss_indices = FakesandMisses(p_matches, good_jets_pyth, misses); //find misses                                                           
	}
      }

      //fakes                                                                                                                                                
      if (good_jets.size() != 0) {
	//clear the vectors to be used for determination of fakes (jets we find in Geant that don't have a match in Pythia)                                  
	fakes.clear(); fake_indices.clear();
	g_matches_for_fakes.clear(); p_matches_for_fakes.clear();
	MatchJets(good_jets_pyth, good_jets, p_matches_for_fakes, g_matches_for_fakes);
	if (g_matches_for_fakes.size() != p_matches_for_fakes.size()) {std::cerr << "Somehow we have different-sized match vectors. Exiting!" <<std::endl; exit(1);}
	if (p_matches_for_fakes.size() < good_jets.size()) { //then we have fakes
	  fake_indices = FakesandMisses(g_matches_for_fakes, good_jets, fakes);
	}
      }

      //need this to weight central embedding to look more like data
      double data_weight = 1;
      //if (bbc_east_sum_pAu < 20000) {data_weight_for_central = 1;}//don't weight periph because the pA is already pp-like
      //if (bbc_east_sum_pAu > 20000) {
      //now weighting peripheral too, so we don't need to check the bbc_east_sum first
      double refmult_weight = hrefmultrat->GetBinContent(hrefmultrat->GetXaxis()->FindBin(refMult_pp+refMult_pAu));
      double bbc_weight = hbbcrat->GetBinContent(hbbcrat->GetXaxis()->FindBin(bbc_east_sum_pAu));
      data_weight *= refmult_weight;
      data_weight *= bbc_weight;
      //data_weight = hrefmultbbcweight->GetBinContent(hrefmultbbcweight->GetXaxis()->FindBin(refMult_pp+refMult_pAu),hrefmultbbcweight->GetYaxis()->FindBin(bbc_east_sum_pAu));
	//}//weight because p+A has a different refMult and BBCE than PYTHIA+GEANT+pp2012ZB+pA2015MB
	
	//data_weight = 1; //for unweighted
	if (data_weight < 0) {cerr << "Weighting response by a negative number. Something's wrong! Exiting!" << endl; exit(1);}

      //we have matches and misses and fakes. Let's fill responses
     // if (bbc_east_sum_pAu < 20000) {
	for (int i = 0; i < misses.size(); ++ i) {
	  truth_unweight->Fill(misses[i].pt(),pp_weight);
	  truth_bbcweight->Fill(misses[i].pt(),pp_weight*bbc_weight);
	  truth_refweight->Fill(misses[i].pt(),pp_weight*refmult_weight);
	  truth_weight->Fill(misses[i].pt(),pp_weight*data_weight);

	  weight_v_pt->Fill(bbc_weight,misses[i].pt(),pp_weight);
	  
	  m_pt_response->Miss(misses[i].m(), misses[i].pt(), pp_weight*data_weight);
	  mg_pt_response->Miss(GroomedJets_pyth[miss_indices[i]].m(), good_jets_pyth[miss_indices[i]].pt(), pp_weight*data_weight);
	}
	//IMPORTANT NOTE:
	//match_indices contains the indices of pairs of geant and pythia matched jets. So if we loop over a list of
	//e.g. geant matches, the index i of a given geant match will be 2*i+1 in the match_indices list.
	//I.e. if we want to access the 3rd geant matched jet (at index 2), we skip over the first 2 pairs of match indices,
	//and the corresponding pythia match, to index at position 5. When accessing pythia jets in the match_indices list,
	//since they come first in each pair, it is similar but with 2*i instead.
	for (int i = 0; i < g_matches.size(); ++ i) { //g_matches.size == p_matches.size == 1/2 (match_indices.size())
	  truth_unweight->Fill(p_matches[i].pt(),pp_weight);
          truth_bbcweight->Fill(p_matches[i].pt(),pp_weight*bbc_weight);
          truth_refweight->Fill(p_matches[i].pt(),pp_weight*refmult_weight);
          truth_weight->Fill(p_matches[i].pt(),pp_weight*data_weight);

	  truth_matches_unweight->Fill(p_matches[i].pt(),pp_weight);
          truth_matches_bbcweight->Fill(p_matches[i].pt(),pp_weight*bbc_weight);
          truth_matches_refweight->Fill(p_matches[i].pt(),pp_weight*refmult_weight);
          truth_matches_weight->Fill(p_matches[i].pt(),pp_weight*data_weight);

	  weight_v_pt->Fill(bbc_weight,p_matches[i].pt(),pp_weight);

	  respT_v_partpT->Fill(g_matches[i].pt() / (double) p_matches[i].pt(), p_matches[i].pt(), pp_weight);
	  resM_v_partpT->Fill(g_matches[i].m() / (double) p_matches[i].m(), p_matches[i].pt(), pp_weight);
	  m_pt_response->Fill(g_matches[i].m(), g_matches[i].pt(), p_matches[i].m(), p_matches[i].pt(), pp_weight*data_weight);
	  mg_pt_response->Fill(GroomedJets[match_indices[(2*i)+1]].m(), good_jets[match_indices[(2*i)+1]].pt(),
			       GroomedJets_pyth[match_indices[2*i]].m(), good_jets_pyth[match_indices[2*i]].pt(), pp_weight*data_weight);
	}
	for (int i = 0; i < fakes.size(); ++ i) {
	  m_pt_response->Fake(fakes[i].m(), fakes[i].pt(), pp_weight*data_weight);
	  mg_pt_response->Fake(GroomedJets[fake_indices[i]].m(), good_jets[fake_indices[i]].pt(), pp_weight*data_weight);
	}
      //}//bbc_east_sum selection for the response
      // And we're done!
      // -----------------------------
      if (good_jets.size() != 0) { // there will be cases where we fill the tree with undefined pythia jets because we only require a non-zero geant jet. this is ok.
	nEventsUsed++; //this event was accepted and had at least one jet passing all criteria
	eventTree->Fill();
      }
      
      break; //THIS IS IMPORTANT: means we don't continue embedding a single pp event into multiple pA events. Once we've used it, we move on to the next one.

      } //Event loop inner
      //cout << "DEBUG: ENDED PA LOOP!!!" << endl;

    } // Event loop outer
    
    //cout << "DEBUG: ENDED PP LOOP!!!" << endl;
    /* // pp should always have no events left because that's how the loop ends
    if (reader_pp->NextEvent()) {
      cout << "DEBUG: pA ENDED FIRST!" << endl;
    }*/
    //but pA should have events left if pp was the cause for ending
    if (reader_pAu->NextEvent()) {
    cout << "DEBUG: pp ENDED FIRST! GOOD!" << endl;
    }
    
  }catch ( std::exception& e) {
    std::cerr << "Caught " << e.what() << std::endl;
    return -1;
  }
  
  // Output
  // ------                                                                                                                                                                                               
  TFile* fout = new TFile((outputDir + outFileName).c_str(), "RECREATE");

  // Close up shop
  // -------------
  //write trees
  eventTree->Write();
  //write responses
  m_pt_response->Write();
  mg_pt_response->Write();
  respT_v_partpT->Write();
  resM_v_partpT->Write();

  truth_unweight->Write();
  truth_bbcweight->Write();
  truth_refweight->Write();
  truth_weight->Write();
  
  truth_matches_unweight->Write();
  truth_matches_bbcweight->Write();
  truth_matches_refweight->Write();
  truth_matches_weight->Write();

  weight_v_pt->Write();
  
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
    
  cout << "DEBUG: " << count_evts << " accepted events, with " << count_jetevts << " MB events skipped due to having > 3.5 GeV jets" << endl;
  
  return 0;
}
  
