//!sim_test.cxx

#include "params.hh"
#include "funcs.hh"
#include "TStarJetPicoDefinitions.h"

using namespace fastjet;
using namespace std;
using namespace Analysis;

//! -------------------------
//! Command line arguments: ( Defaults
//! Defined for debugging in main )
//! [0]: output directory
//! [1]: output filename
//! [2]: flag: ch/ch+ne jets. Options: "ch" (full = 0), or "full" (full = 1).
//! [3]: flag: require match? Options: "nomatch" (match = 0), or "match" (match = 1).
//! [4]: input data: can be a single .root or a .txt or .list of root files - should always be last argument

int main (int argc, const char ** argv) {
  //  TStarJetPicoDefinitions::SetDebugLevel(0);
  TH1::SetDefaultSumw2( );  // Histograms will calculate gaussian errors
  TH2::SetDefaultSumw2( );
  TH3::SetDefaultSumw2( );
  
  // Read in command line arguments
  // ------------------------------
  // Defaults
  std::string outputDir = "out/"; // directory where everything will be saved
  std::string outFileName = "test.root"; //output file
  std::string chainList = "simlist.txt"; // input file: can be .root, .txt, .list
  double radius = 0.4; //jet radius parameter; input value can range from 0.1 to 9.9.
  bool full = 1;  //full = 1 => ch+ne; full = 0 => ch only.
  bool match = 0; //match = 0 => no match between Pythia & Pythia+Geant events.
  
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
    if (arguments[3] == "ch") {full = 0;} else {full = 1;}
    if (arguments[4] == "matched") {match = 1;} else if (arguments[4] == "unmatched") {match = 0;} else {cerr << "Not a valid flag!" << endl; exit(1);}
    chainList         = arguments[5];
    
    //Printing settings:
    cout << "Outputting to: " << (outputDir+outFileName).c_str() << "\nSettings:\nR = " << radius << " " <<  arguments[3] << " jets;\n match pythia and geant? " << match << ";\n input file: " << chainList << "\n";
    break;
  }
  default: { // Error: invalid custom settings
    __ERR("Invalid number of command line arguments");
    return -1;
    break;
  }
  }
  
  //  Initialize readers and provide chains
  TStarJetPicoReader* P6Reader = new TStarJetPicoReader();
  TChain* P6Chain = new TChain( "JetTreeMc" ); // PURE PYTHIA (particle)
  TStarJetPicoReader* GEANTReader = new TStarJetPicoReader();
  TChain* GEANTChain = new TChain( "JetTree" ); // CORRESPONDING GEANT (detector)
  
  // Check to see if the input is a .root file or a .txt
  bool inputIsRoot = Analysis::HasEnding( chainList.c_str(), ".root" );
  bool inputIsTxt  = Analysis::HasEnding( chainList.c_str(), ".txt"  );
  bool inputIsList = Analysis::HasEnding( chainList.c_str(), ".list" );
  
  // If its a recognized file type, build the chain
  // If its not recognized, exit
  if ( inputIsRoot ) { P6Chain->Add( chainList.c_str()); GEANTChain->Add( chainList.c_str());}
  else if ( inputIsTxt )  { P6Chain = TStarJetPicoUtils::BuildChainFromFileList(chainList.c_str()); GEANTChain = TStarJetPicoUtils::BuildChainFromFileList(chainList.c_str());}
  else if ( inputIsList)  { P6Chain = TStarJetPicoUtils::BuildChainFromFileList(chainList.c_str()); GEANTChain = TStarJetPicoUtils::BuildChainFromFileList(chainList.c_str());}
  else { __ERR("data file is not recognized type: .root or .txt only.") return -1; }

  TFile *fout = new TFile( ( outputDir + outFileName ).c_str() ,"RECREATE");
  fout->cd();

  //define relevant data structures
  TString pythiaFilename;
  TStarJetPicoEventHeader* p_header;
  TStarJetPicoEvent* p_event;
  TStarJetVectorContainer<TStarJetVector>* p_container;
  TStarJetVector* p_sv;
  
  TString geantFilename;
  TStarJetPicoEventHeader* g_header;
  TStarJetPicoEvent* g_event;
  TStarJetVectorContainer<TStarJetVector> * g_container;
  TStarJetVector* g_sv;
  
  //defining local containers to be linked to tree branches
  int p_EventID;
  double p_wt;
  vector<double> p_Pt; vector<double> p_Eta; vector<double> p_Phi; vector<double> p_M; vector<double> p_E;
    vector<double> p_Q; vector<double> p_px;
  
  int g_EventID;
  double g_wt;
  vector<double> g_Pt; vector<double> g_Eta; vector<double> g_Phi; vector<double> g_M; vector<double> g_E;
    vector<double> g_Q; vector<double> g_px;
    
    TTree *p_PJTree = new TTree("p_PJTree","p_PJTree");
    p_PJTree->Branch("p_Eta",&p_Eta); p_PJTree->Branch("p_Phi",&p_Phi);
    p_PJTree->Branch("p_px",&p_px); p_PJTree->Branch("p_E",&p_E);
    p_PJTree->Branch("p_Q", &p_Q);
        
    TTree *g_PJTree = new TTree("g_PJTree","g_PJTree");
    g_PJTree->Branch("g_Eta",&g_Eta); g_PJTree->Branch("g_Phi",&g_Phi);
    g_PJTree->Branch("g_px",&g_px); g_PJTree->Branch("g_E",&g_E);
    g_PJTree->Branch("g_Q", &g_Q);
  
  //for later use looking up PDG masses using particle PID
  TDatabasePDG *pdg = new TDatabasePDG();
  
  //  gRandom->SetSeed(1);//for consistency check with v1. Not needed otherwise.
  
  // Particle containers & counters
  vector<PseudoJet> p_Particles, g_Particles;
  int p_n_accepted = 0, g_n_accepted = 0;
  int counter_debug = 0;
  double mc_weight = -1;
  
  double hc = 0.9999; //to be varied in the systematic uncertainty variation
    
  p_n_accepted = 0; g_n_accepted = 0; counter_debug = 0;
    
    //initialize both readers
    InitReader(P6Reader, P6Chain, nEvents, "All", truth_absMaxVz, truth_vZDiff, truth_evPtMax, truth_evEtMax, truth_evEtMin, truth_DCA, truth_NFitPts, truth_FitOverMaxPts, sim_maxEtTow, hc, false, sim_badTowers, sim_bad_run_list);
    //InitReader(GEANTReader, GEANTChain, nEvents, det_triggerString, det_absMaxVz, det_vZDiff, det_evPtMax, det_evEtMax, det_evEtMin, det_DCA, det_NFitPts, det_FitOverMaxPts, sim_maxEtTow, hc, false, /*sim_badTowers, sim_bad_run_list);*/det_badTowers, dat_bad_run_list);
        InitReader(GEANTReader, GEANTChain, nEvents, "All", truth_absMaxVz, truth_vZDiff, truth_evPtMax, truth_evEtMax, truth_evEtMin, truth_DCA, truth_NFitPts, truth_FitOverMaxPts, sim_maxEtTow, hc, false, sim_badTowers, sim_bad_run_list);
    
    // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
    
    for (int event = 0; event < P6Chain->GetEntries(); ++ event) {
      P6Reader->ReadEvent(event);
      GEANTReader->ReadEvent(event);
      
      g_EventID = GEANTReader->GetNOfCurrentEvent();
      p_EventID = P6Reader->GetNOfCurrentEvent();
      
      //clearing vectors; initializing variables to -9999
      mc_weight = -9999;
      p_wt = -9999;
      p_Pt.clear(); p_Eta.clear(); p_Phi.clear(); p_M.clear(); p_E.clear();
      p_Q.clear(); p_px.clear();
      
      g_wt = -9999;
      g_Pt.clear(); g_Eta.clear(); g_Phi.clear(); g_M.clear(); g_E.clear();
      g_Q.clear(); g_px.clear();
      
      p_Particles.clear(); g_Particles.clear();
      
      if (P6Reader->ReadEvent(p_EventID) == 1) {p_n_accepted++;} // == 1 => event loaded and passed selections
      if (GEANTReader->ReadEvent(g_EventID) == 1) {g_n_accepted++;}
            
      P6Reader->PrintStatus(10); GEANTReader->PrintStatus(10);     // Print out reader status every 10 seconds
      
      //filling the data structures that were defined before the event loop:
      p_event = P6Reader->GetEvent();
      p_header = p_event->GetHeader();
      g_event = GEANTReader->GetEvent();
      g_header = g_event->GetHeader();
      
      p_container = P6Reader->GetOutputContainer();
      g_container = GEANTReader->GetOutputContainer();
      
      pythiaFilename =  P6Reader->GetInputChain()->GetCurrentFile()->GetName();
      geantFilename =  GEANTReader->GetInputChain()->GetCurrentFile()->GetName();
      
      p_wt = LookupRun12Xsec( pythiaFilename );
      g_wt = LookupRun12Xsec( geantFilename );
      mc_weight = p_wt; //arbitrarily setting it to pythia's but can be either
      
      // converts TStarJetVectors to PseudoJets, carefully assigning the proper particle mass in either case
      GatherParticles ( p_container, p_sv, p_Particles, full, 1, pdg); //Pythia; full = 0 => charged-only, 1 => ch+ne
      GatherParticles ( g_container, g_sv, g_Particles, full, 1, pdg); //GEANT //CHANGED THIS TO BE PYTHIA-LIKE IN MASS ASSIGNMENT FOR TESTING PURPOSES
            
        for (int i = 0; i < p_Particles.size(); ++ i) {
            p_px.push_back(p_Particles[i].px());
            p_E.push_back(p_Particles[i].E());
            p_Q.push_back(p_Particles[i].user_index());
            p_Eta.push_back(p_Particles[i].eta());
            p_Phi.push_back(p_Particles[i].phi());
        }
        
        if (p_Particles.size() != 0) {
            p_PJTree->Fill(); //when !match, will fill sometimes with empty geant vectors in the geant branches
        }
        
        for (int i = 0; i < g_Particles.size(); ++ i) {
            g_px.push_back(g_Particles[i].px());
            g_E.push_back(g_Particles[i].E());
            g_Q.push_back(g_Particles[i].user_index());
            g_Eta.push_back(g_Particles[i].eta());
            g_Phi.push_back(g_Particles[i].phi());
        }
        
        if (g_Particles.size() != 0) {
            g_PJTree->Fill();
        }

    }//event loop
    
    //~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ END EVENT LOOP! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
    
    fout->cd();
    
    cout << endl << endl << p_n_accepted << " pythia events and " << g_n_accepted << " geant events" << endl;
    
    //trees
    p_PJTree->Write();
    g_PJTree->Write();
      
    cout << endl << "Writing to:  " << fout->GetName() << endl;
    
    fout->Close();
    
  return 0;
}//main
