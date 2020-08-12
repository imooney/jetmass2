//Isaac Mooney, WSU - November 2019
//This file embeds ppJP2 data into pAuMB data
#include <TFile.h>

#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TString.h>
#include <TChain.h>
#include <TBranch.h>
#include <TMath.h>
#include <TRandom.h>
#include <TSystem.h>

#include "fastjet/PseudoJet.hh"

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

  // Read in command line arguments
  // ------------------------------
  // Defaults
  std::string executable = "./bin/data"; // placeholder                                                                                    
  std::string chainList_pyth = "list.txt"; // input file: can be .root, .txt, .list                                                                           
 
  std::string chainNameMc = "JetTreeMc";
  
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
    chainList_pyth = arguments[4];
    
    break;
  }
  default: { // Error: invalid custom settings
    __ERR("Invalid number of command line arguments");
      return -1;
    break;
  }
  }
 
    TChain* chain_pyth = new TChain( "JetTreeMc" );//chainNameMc.c_str() );
  cout << "CHAINNAMEMC! " << chainNameMc << endl;

  //chain_pyth->Add("/tier2/home/groups/rhi/STAR/Data/AddedEmbedPythiaRun12pp200/Cleanpp12Pico_pt25_35_g1.root");
  chain_pyth = TStarJetPicoUtils::BuildChainFromFileList( chainList_pyth.c_str(), "JetTreeMc", -1, 0, chain_pyth );
  
  // Build the event structure w/ cuts
  // ---------------------------------
  TStarJetPicoReader * reader_pyth = new TStarJetPicoReader();
  InitReader(reader_pyth, chain_pyth, nEvents, "All", truth_absMaxVz, truth_vZDiff, truth_evPtMax, truth_evEtMax, truth_evEtMin, truth_DCA, truth_NFitPts, truth_FitOverMaxPts, sim_maxEtTow, 0.9999, false, sim_badTowers, sim_bad_run_list);

  // Data classes
  // ------------
  TStarJetVectorContainer<TStarJetVector>* container_pyth;
  TStarJetVector* sv_pyth; // TLorentzVector* would be sufficient
  TString filename_pyth;
  TStarJetPicoEventHeader* header_pyth;
  TStarJetPicoEvent* event_pyth;
  TStarJetPicoEventCuts* EventCuts_pyth = reader_pyth->GetEventCuts();//will use this for hardcoding checks if events passed selections
  
  double pp_weight;

  try{
 
    for (int event = 0; event < chain_pyth->GetEntries(); ++ event) {
      reader_pyth->ReadEvent(event);
      
      pp_weight = -9999;
      container_pyth = reader_pyth->GetOutputContainer();

      cout << "DEBUG: smoking gun: " << reader_pyth->GetEvent()->GetHeader()->GetBbcEastRate() << endl;
      
      cout << "DEBUG: EVENT# " << reader_pyth->GetNOfCurrentEvent() << endl;
      for (int i = 0; i < container_pyth->GetEntries(); ++ i) {
	cout << "DEBUG: " << container_pyth->Get(i)->perp() << endl;
      }
    } // Event loop
    
  }catch ( std::exception& e) {
    std::cerr << "Caught " << e.what() << std::endl;
    return -1;
  }
  return 0;
}
  
