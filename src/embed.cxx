//!embed.cxx
//!Isaac Mooney, WSU - March 2020
//!This file runs the (pA) analysis on the STAR embedding of PYTHIA-6 into p+Au events for the jet mass project.
//!It takes in the Picos, performs selections, clusters particles, performs selections on the resulting jets,
//!applies the Soft Drop grooming procedure to a copy of the jet population, fills QA trees, jet trees, and responses, and writes to files.

#include "params.hh"
#include "funcs.hh"
#include "TStarJetPicoDefinitions.h"

using namespace fastjet;
using namespace std;
using namespace Analysis;
typedef fastjet::contrib::SoftDrop SD;

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
  std::string trigger = "pAuJP2";
  double radius = 0.4; //jet radius parameter; input value can range from 0.1 to 9.9.
  bool full = 1;  //full = 1 => ch+ne; full = 0 => ch only.
  bool match = 0; //match = 0 => no match between Pythia & Pythia+Geant events.
  
  // Now check to see if we were given modifying arguments
  switch ( argc ) {
  case 1: // Default case
    __OUT("Using Default Settings");
    break;
  case 8: { // Custom case
    __OUT("Using Custom Settings");
    std::vector<std::string> arguments( argv+1, argv+argc );
    
    // Set non-default values
    // ----------------------
    // output and file names
    outputDir         = arguments[0];
    outFileName       = arguments[1];
    radius            = radius_str_to_double (arguments[2]);
    trigger           = arguments[3]; //ppJP2, ppHT2, ppVPDMB, pAuJP2, pAuHT2, pAuBBCMB, or AA [TBD]     
    if (arguments[4] == "ch") {full = 0;} else {full = 1;}
    if (arguments[5] == "matched") {match = 1;} else if (arguments[5] == "unmatched") {match = 0;} else {cerr << "Not a valid flag!" << endl; exit(1);}
    chainList         = arguments[6];
    
    //test:
    for (int i = 0; i < 7; ++ i) {
      cout << "no argument there! " << arguments[i] << endl;
    }

    //Printing settings:
    cout << "Outputting to: " << (outputDir+outFileName).c_str() << "\nSettings:\nR = " << radius << " " <<  arguments[4] << " jets;\n match pythia and geant? " << match << ";\n input file: " << chainList << "\n";
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


  bool lowEA = 0;//sets activity to low or high (currently 60-90% and 0-30% respectively).

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
  
  /*
  //~~~
  TFile *fweight = new TFile(("~/jetmass2/out/embed/weight_for_embedding_R"+(string) argv[3]+".root").c_str(),"READ");
  TH1D *refmult_rat_lowEA = (TH1D*) fweight->Get("hrat_ref_lowEA");
  TH1D *refmult_rat_highEA = (TH1D*) fweight->Get("hrat_ref_highEA");
  TH1D *pt_rat_lowEA = (TH1D*) fweight->Get("hrat_pt_lowEA");
  TH1D* pt_rat_highEA = (TH1D*) fweight->Get("hrat_pt_highEA");
  TF1* frefmult_rat_lowEA = (TF1*) fweight->Get("frat_ref_lowEA");
  TF1* frefmult_rat_highEA = (TF1*) fweight->Get("frat_ref_highEA");
  refmult_rat_lowEA->SetDirectory(0);
  refmult_rat_highEA->SetDirectory(0);
  pt_rat_lowEA->SetDirectory(0);
  pt_rat_highEA->SetDirectory(0);
  fweight->Close();
  //~~~
  */
  /*
  //TEST! UNCOMMENT ABOVE AFTER
  TFile *fweight = new TFile(("~/jetmass2/out/embed/weight_for_embedding_mpt_R"+(string) argv[3]+".root").c_str(),"READ");
  TH1D *pt_rat_lowEA = (TH1D*) fweight->Get("hrat_pt_lowEA");
  TH1D* pt_rat_highEA = (TH1D*) fweight->Get("hrat_pt_highEA");
  pt_rat_lowEA->SetDirectory(0);
  pt_rat_highEA->SetDirectory(0);
  TH1D *m1520_rat_lowEA = (TH1D*) fweight->Get("hrat_m1520_lowEA");
  TH1D* m1520_rat_highEA = (TH1D*) fweight->Get("hrat_m1520_highEA");
  TH1D *m2025_rat_lowEA = (TH1D*) fweight->Get("hrat_m2025_lowEA");
  TH1D* m2025_rat_highEA = (TH1D*) fweight->Get("hrat_m2025_highEA");
  TH1D *m2530_rat_lowEA = (TH1D*) fweight->Get("hrat_m2530_lowEA");
  TH1D* m2530_rat_highEA = (TH1D*) fweight->Get("hrat_m2530_highEA");
  TH1D *m3040_rat_lowEA = (TH1D*) fweight->Get("hrat_m3040_lowEA");
  TH1D* m3040_rat_highEA = (TH1D*) fweight->Get("hrat_m3040_highEA");
  TH1D *m40up_rat_lowEA = (TH1D*) fweight->Get("hrat_m40up_lowEA");
  TH1D* m40up_rat_highEA = (TH1D*) fweight->Get("hrat_m40up_highEA");
  m1520_rat_lowEA->SetDirectory(0);
  m1520_rat_highEA->SetDirectory(0);
  m2025_rat_lowEA->SetDirectory(0);
  m2025_rat_highEA->SetDirectory(0);
  m2530_rat_lowEA->SetDirectory(0);
  m2530_rat_highEA->SetDirectory(0);
  m3040_rat_lowEA->SetDirectory(0);
  m3040_rat_highEA->SetDirectory(0);
  m40up_rat_lowEA->SetDirectory(0);
  m40up_rat_highEA->SetDirectory(0);
  
  fweight->Close();
  */

  TFile *fweightrefptfit = new TFile(("~/jetmass2/out/embed/weight_for_embedding_allmethods_R"+(string) argv[3]+".root").c_str(),"READ");
  TF1* frefmult_rat_lowEA = (TF1*) fweightrefptfit->Get("frat_ref_lowEA");                                                        
  TF1* frefmult_rat_highEA = (TF1*) fweightrefptfit->Get("frat_ref_highEA");    
  TF2* frefmult_pt_rat_lowEA = (TF2*) fweightrefptfit->Get("frat_ref_pt_lowEA");
  TF2* frefmult_pt_rat_highEA = (TF2*) fweightrefptfit->Get("frat_ref_pt_highEA");
  TH2D* pt_m_rat_lowEA = (TH2D*) fweightrefptfit->Get("hrat_pt_m_lowEA");
  TH2D* pt_m_rat_highEA = (TH2D*) fweightrefptfit->Get("hrat_pt_m_highEA");
  pt_m_rat_lowEA->SetDirectory(0);
  pt_m_rat_highEA->SetDirectory(0);
  fweightrefptfit->Close();
  


  TFile *fout = new TFile( ( outputDir + outFileName ).c_str() ,"RECREATE");
  fout->cd();
  
  
  //define relevant data structures
  TString pythiaFilename;
  TStarJetPicoEventHeader* p_header;
  TStarJetPicoEvent* p_event;
  TStarJetVectorContainer<TStarJetVector> * p_container;
  TStarJetVector* p_sv;
  /*
  TClonesArray* p_prim_trks;
  TClonesArray* g_prim_trks;
  TClonesArray* p_tows;
  TClonesArray* g_tows;
  */
  TString geantFilename;
  TStarJetPicoEventHeader* g_header;
  TStarJetPicoEvent* g_event;
  TStarJetVectorContainer<TStarJetVector> * g_container;
  TStarJetVector* g_sv;
  
  //defining local containers to be linked to tree branches
  int p_EventID;
  double p_n_jets, p_wt;
  vector<vector<double> > p_conspT;
  vector<double> p_jetMult;
  vector<double> p_Pt; vector<double> p_Eta; vector<double> p_Phi; vector<double> p_M; vector<double> p_E;
  vector<double> p_ch_e_frac;
  vector<double> p_zg; vector<double> p_rg; vector<double> p_mg; vector<double> p_ptg;
  vector<double> p_mcd;
  
  int g_EventID;
  double g_n_jets, g_wt;
  double event_weight;
  vector<vector<double> > g_conspT;
  vector<double> g_jetMult;
  vector<double> g_Pt; vector<double> g_Eta; vector<double> g_Phi; vector<double> g_M; vector<double> g_E;
  vector<double> g_ch_e_frac;
  vector<double> g_zg; vector<double> g_rg; vector<double> g_mg; vector<double> g_ptg;
  vector<double> g_mcd;

    //QA containers
    //EventHeader
    int p_EventId; int p_RunId;
    int p_RefMult; double p_UE_avgpt; double p_UE_npart;
    int p_NOfGlobalTracks; int p_NOfTowers; int p_NOfPrimaryTracks;
    double p_PVx; double p_PVy; double p_PVz;
    double p_vpdVz;
    double p_BbcEastRate; double p_BbcCoincidenceRate; double p_BbcAdcSumEast;
    //Primaries
    vector<int> p_prim_Charge;
    vector<double> p_prim_Px; vector<double> p_prim_Py; vector<double> p_prim_Pz; vector<double> p_prim_Pt;
    vector<double> p_prim_DCA; vector<double> p_prim_dEdx;
    //Towers
    vector<int> p_tow_Id;
    vector<double> p_tow_Energy; vector<double> p_tow_Et;
    vector<double> p_tow_Eta; vector<double> p_tow_Phi; vector<double> p_tow_EtaCorrected; vector<double> p_tow_PhiCorrected;
    
    //EventHeader
    int g_EventId; int g_RunId;
    int g_RefMult; double g_UE_avgpt; double g_UE_npart;
    int g_NOfGlobalTracks; int g_NOfTowers; int g_NOfPrimaryTracks;
    double g_PVx; double g_PVy; double g_PVz;
    double g_vpdVz;
    double g_BbcEastRate; double g_BbcCoincidenceRate; double g_BbcAdcSumEast;
    double g_ZdcCoincidenceRate;
    //Primaries
    vector<int> g_prim_Charge;
    vector<double> g_prim_Px; vector<double> g_prim_Py; vector<double> g_prim_Pz; vector<double> g_prim_Pt;
    vector<double> g_prim_DCA; vector<double> g_prim_dEdx;
    //Towers
    vector<int> g_tow_Id;
    vector<double> g_tow_Energy; vector<double> g_tow_Et;
    vector<double> g_tow_Eta; vector<double> g_tow_Phi; vector<double> g_tow_EtaCorrected; vector<double> g_tow_PhiCorrected;

    vector<vector<double> > suspect_pt; vector<vector<double> > suspect_eta; vector<vector<double> > suspect_phi; vector<vector<double> > suspect_Q;
    
    vector<double> g_part_pt; vector<double> g_part_eta; vector<double> p_part_pt; vector<double> p_part_eta;//temp

    vector<double> py_part_eta;//temp
    
    TTree *temptree = new TTree("temptree","temptree");
    temptree->Branch("p_RefMult",&p_RefMult);
    temptree->Branch("g_RefMult",&g_RefMult);
    temptree->Branch("g_BbcAdcSumEast",&g_BbcAdcSumEast);
    temptree->Branch("weight",&g_wt);
    temptree->Branch("g_part_eta",&g_part_eta);
    temptree->Branch("g_part_pt",&g_part_pt);
    temptree->Branch("p_part_eta",&p_part_eta);
    temptree->Branch("p_part_pt",&p_part_pt);
    

    //DEBUG tree
    TTree *suspects = new TTree("suspects","suspects");
    suspects->Branch("weight",&g_wt);
    suspects->Branch("g_EventId",&g_EventId);
    suspects->Branch("g_RunId",&g_RunId);
    suspects->Branch("g_RefMult",&g_RefMult);
    suspects->Branch("g_BbcAdcSumEast",&g_BbcAdcSumEast);
    suspects->Branch("pt",&suspect_pt);
    suspects->Branch("eta",&suspect_eta);
    suspects->Branch("phi",&suspect_phi);
    suspects->Branch("Q",&suspect_Q);
    
    
    //tree for QA of primary tracks, towers, event quantities
    TTree *p_QATree = new TTree("p_QA","p_QA");
    //pure PYTHIA:
    p_QATree->Branch("p_weight", &p_wt);
    p_QATree->Branch("p_EventId",&p_EventId);
    p_QATree->Branch("p_RunId",&p_RunId);
    p_QATree->Branch("p_RefMult",&p_RefMult);
    p_QATree->Branch("p_NOfGlobalTracks",&p_NOfGlobalTracks);
    p_QATree->Branch("p_NOfTowers",&p_NOfTowers);
    p_QATree->Branch("p_NOfPrimaryTracks",&p_NOfPrimaryTracks);
    p_QATree->Branch("p_PVx",&p_PVx);
    p_QATree->Branch("p_PVy",&p_PVy);
    p_QATree->Branch("p_PVz",&p_PVz);
    p_QATree->Branch("p_vpdVz",&p_vpdVz);
    p_QATree->Branch("p_BbcEastRate",&p_BbcEastRate);
    p_QATree->Branch("p_BbcCoincidenceRate",&p_BbcCoincidenceRate);
    p_QATree->Branch("p_BbcAdcSumEast",&p_BbcAdcSumEast);
    
    p_QATree->Branch("p_prim_Charge",&p_prim_Charge);
    p_QATree->Branch("p_prim_Px",&p_prim_Px);
    p_QATree->Branch("p_prim_Py",&p_prim_Py);
    p_QATree->Branch("p_prim_Pz",&p_prim_Pz);
    p_QATree->Branch("p_prim_Pt",&p_prim_Pt);
    p_QATree->Branch("p_prim_DCA",&p_prim_DCA);
    p_QATree->Branch("p_prim_dEdx",&p_prim_dEdx);
    
    p_QATree->Branch("p_tow_Id",&p_tow_Id);
    p_QATree->Branch("p_tow_Energy",&p_tow_Energy);
    p_QATree->Branch("p_tow_Et",&p_tow_Et);
    p_QATree->Branch("p_tow_Eta",&p_tow_Eta);
    p_QATree->Branch("p_tow_Phi",&p_tow_Phi);
    p_QATree->Branch("p_tow_EtaCorrected",&p_tow_EtaCorrected);
    p_QATree->Branch("p_tow_PhiCorrected",&p_tow_PhiCorrected);
    
    p_QATree->Branch("py_part_eta",&py_part_eta);

    TTree *g_QATree = new TTree("g_QA","g_QA");
    //PYTHIA embedded into min-bias pAu data:
    g_QATree->Branch("g_weight", &g_wt);
    g_QATree->Branch("g_EventId",&g_EventId);
    g_QATree->Branch("g_RunId",&g_RunId);
    g_QATree->Branch("g_RefMult",&g_RefMult);
    g_QATree->Branch("g_NOfGlobalTracks",&g_NOfGlobalTracks);
    g_QATree->Branch("g_NOfTowers",&g_NOfTowers);
    g_QATree->Branch("g_NOfPrimaryTracks",&g_NOfPrimaryTracks);
    g_QATree->Branch("g_PVx",&g_PVx);
    g_QATree->Branch("g_PVy",&g_PVy);
    g_QATree->Branch("g_PVz",&g_PVz);
    g_QATree->Branch("g_vpdVz",&g_vpdVz);
    g_QATree->Branch("g_BbcEastRate",&g_BbcEastRate);
    g_QATree->Branch("g_ZdcCoincidenceRate",&g_ZdcCoincidenceRate);
    g_QATree->Branch("g_BbcCoincidenceRate",&g_BbcCoincidenceRate);
    g_QATree->Branch("g_BbcAdcSumEast",&g_BbcAdcSumEast);
    
    g_QATree->Branch("g_prim_Charge",&g_prim_Charge);
    g_QATree->Branch("g_prim_Px",&g_prim_Px);
    g_QATree->Branch("g_prim_Py",&g_prim_Py);
    g_QATree->Branch("g_prim_Pz",&g_prim_Pz);
    g_QATree->Branch("g_prim_Pt",&g_prim_Pt);
    g_QATree->Branch("g_prim_DCA",&g_prim_DCA);
    g_QATree->Branch("g_prim_dEdx",&g_prim_dEdx);
    
    g_QATree->Branch("g_tow_Id",&g_tow_Id);
    g_QATree->Branch("g_tow_Energy",&g_tow_Energy);
    g_QATree->Branch("g_tow_Et",&g_tow_Et);
    g_QATree->Branch("g_tow_Eta",&g_tow_Eta);
    g_QATree->Branch("g_tow_Phi",&g_tow_Phi);
    g_QATree->Branch("g_tow_EtaCorrected",&g_tow_EtaCorrected);
    g_QATree->Branch("g_tow_PhiCorrected",&g_tow_PhiCorrected);
    
    
    
    
  //tree to hold jet and constituent quantites
  TTree *eventTree = new TTree("event","event");
  eventTree->Branch("p_n_jets", &p_n_jets);
  eventTree->Branch("p_conspT",&p_conspT);
  eventTree->Branch("p_jetMult",&p_jetMult);
  eventTree->Branch("p_UE_avgpt",&p_UE_avgpt);
  eventTree->Branch("p_UE_npart",&p_UE_npart);
  eventTree->Branch("p_Pt", &p_Pt); eventTree->Branch("p_Eta",&p_Eta); eventTree->Branch("p_Phi",&p_Phi); eventTree->Branch("p_M",&p_M); eventTree->Branch("p_E",&p_E);
  eventTree->Branch("p_ch_e_frac", &p_ch_e_frac);
  eventTree->Branch("p_zg", &p_zg); eventTree->Branch("p_rg", &p_rg); eventTree->Branch("p_mg", &p_mg); eventTree->Branch("p_ptg",&p_ptg);
  eventTree->Branch("p_mcd",&p_mcd);
  eventTree->Branch("p_weight", &p_wt); eventTree->Branch("p_EventID", &p_EventID);
  
  eventTree->Branch("g_n_jets", &g_n_jets);
  eventTree->Branch("g_BbcAdcSumEast",&g_BbcAdcSumEast);
  eventTree->Branch("g_conspT",&g_conspT);
  eventTree->Branch("g_jetMult",&g_jetMult);
  eventTree->Branch("g_RefMult",&g_RefMult);
  eventTree->Branch("g_UE_avgpt",&g_UE_avgpt);
  eventTree->Branch("g_UE_npart",&g_UE_npart);
  eventTree->Branch("g_Pt", &g_Pt); eventTree->Branch("g_Eta",&g_Eta); eventTree->Branch("g_Phi",&g_Phi); eventTree->Branch("g_M",&g_M); eventTree->Branch("g_E",&g_E);
  eventTree->Branch("g_ch_e_frac", &g_ch_e_frac);
  eventTree->Branch("g_zg", &g_zg); eventTree->Branch("g_rg", &g_rg); eventTree->Branch("g_mg", &g_mg); eventTree->Branch("g_ptg",&g_ptg);
  eventTree->Branch("g_mcd",&g_mcd);
  eventTree->Branch("g_weight", &g_wt); eventTree->Branch("g_EventID", &g_EventID);
  eventTree->Branch("event_weight",&event_weight);//we try weighting the events based on the refmult so embedding looks more like data


  //hist for statistical error correction - for matched events, but unmatched jets                                                                             
  TH1D* pt_gen_match_plus_miss = new TH1D("pt_gen_match_plus_miss","",15,5,80);
  
  
  //hists for use in responses - note:
  //"hist_measured and hist_truth are used to specify the dimensions of the distributions (the histogram contents are not used here), eg. for 2D or 3D distributions or non-uniform binning." - http://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html
  //I.e. these histograms are just used as a template, they're not actually what is being filled when the responses are constructed.
  TH2D *pyMvPt = new TH2D("pyMvPt",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMvPt = new TH2D("geMvPt",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyZgvPt = new TH2D("pyZgvPt", ";z_{g};p_{T} [GeV/c]",20,0,1,15,5,80);
  TH2D *geZgvPt = new TH2D("geZgvPt", ";z_{g};p_{T} [GeV/c]",20,0,1,9,15,60);
  TH2D *pyRgvPt = new TH2D("pyRgvPt", ";R_{g};p_{T} [GeV/c]",20,0,1,15,5,80);
  TH2D *geRgvPt = new TH2D("geRgvPt", ";R_{g};p_{T} [GeV/c]",20,0,1,9,15,60);
  TH2D *pyPtgvPt = new TH2D("pyPtgvPt", ";p_{T,g} [GeV/c];p_{T} [GeV/c]",15,5,80,15,5,80);
  TH2D *gePtgvPt = new TH2D("gePtgvPt", ";p_{T,g} [GeV/c];p_{T} [GeV/c]",9,15,60,9,15,60);
  TH2D *pyMgvPt = new TH2D("pyMgvPt", ";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMgvPt = new TH2D("geMgvPt", ";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
    
  // 1D responses
  RooUnfoldResponse *pt_response = new RooUnfoldResponse(9,15,60,15,5,80,"pt_response","");
  RooUnfoldResponse *m_response = new RooUnfoldResponse(14,0,14,14,0,14,"m_response","");
  RooUnfoldResponse *zg_response = new RooUnfoldResponse(20,0,1,20,0,1, "zg_response","");
  RooUnfoldResponse *rg_response = new RooUnfoldResponse(20,0,1,20,0,1, "rg_response","");
  RooUnfoldResponse *ptg_response = new RooUnfoldResponse(9,15,60,15,5,80, "ptg_response","");
  RooUnfoldResponse *mg_response = new RooUnfoldResponse(14,0,14,14,0,14, "mg_response","");
    
  // 2D responses
  RooUnfoldResponse *m_pt_response_lax = new RooUnfoldResponse(geMvPt, pyMvPt, "m_pt_response_lax");
  RooUnfoldResponse *m_pt_response_med = new RooUnfoldResponse(geMvPt, pyMvPt, "m_pt_response_med");
  RooUnfoldResponse *m_pt_response_strict = new RooUnfoldResponse(geMvPt, pyMvPt, "m_pt_response_strict");
  RooUnfoldResponse *m_pt_response_counts = new RooUnfoldResponse(geMvPt, pyMvPt, "m_pt_response_counts");
  RooUnfoldResponse *zg_pt_response = new RooUnfoldResponse(geZgvPt, pyZgvPt, "zg_pt_response");
  RooUnfoldResponse *rg_pt_response = new RooUnfoldResponse(geRgvPt, pyRgvPt, "rg_pt_response");
  RooUnfoldResponse *ptg_pt_response = new RooUnfoldResponse(gePtgvPt, pyPtgvPt, "ptg_pt_response");
  RooUnfoldResponse *mg_pt_response = new RooUnfoldResponse(geMgvPt, pyMgvPt, "mg_pt_response");
  RooUnfoldResponse *mg_pt_response_counts = new RooUnfoldResponse(geMgvPt, pyMgvPt, "mg_pt_response_counts");
  
  RooUnfoldResponse *m_pt_res_nom = new RooUnfoldResponse(geMvPt, pyMvPt, "m_pt_res_nom"); //nominal
  RooUnfoldResponse *m_pt_res_DS = new RooUnfoldResponse(geMvPt, pyMvPt, "m_pt_res_DS"); //smear detector spectrum  

  //vectors of responses & hists for easy writing to file later
  std::vector<RooUnfoldResponse*> res = {pt_response,m_response,zg_response,rg_response,ptg_response,mg_response,m_pt_response_lax,m_pt_response_counts,zg_pt_response,rg_pt_response,ptg_pt_response,mg_pt_response,mg_pt_response_counts, m_pt_response_med, m_pt_response_strict};    
    
  std::vector<RooUnfoldResponse*> res_syst = {m_pt_res_nom,m_pt_res_DS};
    
  //defining the algorithm and radius parameter for clustering jets
  JetDefinition jet_def(antikt_algorithm, radius/*R*/);
  
  //Creating SoftDrop grooming object
  contrib::SoftDrop sd(Beta,z_cut,R0);
  cout << "SoftDrop groomer is: " << sd.description() << endl;
  
  //for later use looking up PDG masses using particle PID
  TDatabasePDG *pdg = new TDatabasePDG();
  
  //  gRandom->SetSeed(1);//for consistency check with v1. Not needed otherwise.
  
  //SELECTORS
  // Constituent selectors
  // ---------------------
  Selector select_track_eta = fastjet::SelectorAbsEtaMax(max_track_eta);
  Selector select_lopt      = fastjet::SelectorPtMin( partMinPt );
  Selector select_loptmax   = fastjet::SelectorPtMax( partMaxPt );
  Selector spart = select_track_eta * select_lopt * select_loptmax;
  
  //for tests
  //  Selector spart_nocuts = fastjet::SelectorPtMax(1000);
    Selector spart_test = select_lopt;

  // Jet candidate selectors
  // -----------------------
  Selector select_jet_eta     = fastjet::SelectorAbsEtaMax(max_eta);
  Selector select_det_jet_pt_min  = fastjet::SelectorPtMin( det_jet_ptmin );
  Selector select_gen_jet_pt_min = fastjet::SelectorPtMin( jet_ptmin );
  Selector select_jet_pt_max  = fastjet::SelectorPtMax( jet_ptmax );
  Selector select_det_jet_m_min = fastjet::SelectorMassMin( mass_min );
  Selector select_gen_jet_m_min = fastjet::SelectorMassMin( 0.0 );
 

  //TEMP FOR COMPARISON! CHANGE BACK LATER!
  Selector select_temp_jetmin = fastjet::SelectorPtMin( 15.0 );
  Selector select_temp_jetmax = fastjet::SelectorPtMax( 30.0 );
  
  Selector sjet_gen = select_jet_eta && select_gen_jet_pt_min && select_jet_pt_max /*select_temp_jetmin && select_temp_jetmax*/ && select_gen_jet_m_min;
 
  Selector sjet_det = select_jet_eta && select_det_jet_pt_min && select_jet_pt_max && select_det_jet_m_min;
  
  //for tests
  //Selector sjet_nocuts = fastjet::SelectorPtMax(1000);
  //Selector sjet_test = select_gen_jet_pt_min;

  // Particle containers & counters
  vector<PseudoJet> p_Particles, g_Particles, p_JetsInitial, g_JetsInitial;
  vector<PseudoJet> p_Particles_nocuts, g_Particles_nocuts;
  vector<PseudoJet> p_Jets_nocuts, g_Jets_nocuts;
  int p_n_accepted = 0, g_n_accepted = 0; int p_NJets = 0, g_NJets = 0;
  int counter_debug = 0;
  double mc_weight = -1;
  double nEvts_for_weight = 0; //use this to divide total cross sections to make them per-event
  
  //debug info:
  vector<int> debug_jets_Py = {0,0,0,0,0};
  vector<int> debug_jets_Ge = {0,0,0,0,0};
  vector<int> debug_events_Py = {0,0,0,0,0};
  vector<int> debug_events_Ge = {0,0,0,0,0};
  
  double hc = 0.9999; //to be varied in the systematic uncertainty variation
  const int nSources = 1; //includes the nominal settings as a "systematic".//TEMP!
  for (int iSyst = 0; iSyst < nSources; ++ iSyst) {
    if (iSyst == 0) {cout << endl << "RUNNING WITH NOMINAL SETTINGS!" << endl << endl;}
    if (iSyst == 1) {cout << endl << "RUNNING WITH ADJUSTED DET-LEVEL pT SPECTRUM!" << endl << endl;}

    p_NJets = 0; g_NJets = 0; p_n_accepted = 0; g_n_accepted = 0; counter_debug = 0;
    //set parameters back to their nominal values after the previous iteration changed them.
    //hc = 0.9999;
    
    debug_events_Py[0] = P6Chain->GetEntries();
    debug_events_Ge[0] = GEANTChain->GetEntries();
    
    //initialize both readers
    InitReader(P6Reader, P6Chain, nEvents, "All", truth_absMaxVz, truth_vZDiff, truth_evPtMax, truth_evEtMax, truth_evEtMin, truth_DCA, truth_NFitPts, truth_FitOverMaxPts, sim_maxEtTow, hc, false, sim_badTowers, sim_bad_run_list);
    InitReader(GEANTReader, GEANTChain, nEvents, "All", det_absMaxVz, /*det_vZDiff*/ vzdiff, det_evPtMax, det_evEtMax, det_evEtMin, det_DCA, det_NFitPts, det_FitOverMaxPts, sim_maxEtTow, hc, false, /*pAu_badTowers*/badtows, pAu_bad_run_list);
    
    // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
    
    debug_events_Py[1] = P6Chain->GetEntries();
    debug_events_Ge[1] = GEANTChain->GetEntries();

    cout << "DEBUG: " << P6Chain->GetEntries() << " total py events" << endl;
    cout << "DEBUG: " << GEANTChain->GetEntries() << " total ge events" << endl;
    for (int event = 0; event < P6Chain->GetEntries(); ++ event) {
      P6Reader->ReadEvent(event);
     
      p_EventID = P6Reader->GetNOfCurrentEvent();
      //clearing vectors; initializing variables to -9999
        p_EventId = -9999; p_RunId = -9999;
        p_RefMult = -9999;
        p_NOfGlobalTracks = -9999; p_NOfTowers = -9999; p_NOfPrimaryTracks = -9999;
        p_PVx = -9999; p_PVy = -9999; p_PVz = -9999;
        p_vpdVz = -9999;
        p_BbcEastRate = -9999; p_BbcCoincidenceRate = -9999; p_BbcAdcSumEast = -9999;
        //Primaries
        p_prim_Charge.clear();
        p_prim_Px.clear(); p_prim_Py.clear(); p_prim_Pz.clear(); p_prim_Pt.clear();
        p_prim_DCA.clear(); p_prim_dEdx.clear();
        //Towers
        p_tow_Id.clear();
        p_tow_Energy.clear(); p_tow_Et.clear();
        p_tow_Eta.clear(); p_tow_Phi.clear(); p_tow_EtaCorrected.clear(); p_tow_PhiCorrected.clear();
        
        g_EventId = -9999; g_RunId = -9999;
        g_RefMult = -9999;
        g_NOfGlobalTracks = -9999; g_NOfTowers = -9999; g_NOfPrimaryTracks = -9999;
        g_PVx = -9999; g_PVy = -9999; g_PVz = -9999;
        g_vpdVz = -9999;
        g_BbcEastRate = -9999; g_BbcCoincidenceRate = -9999; g_BbcAdcSumEast = -9999;
	g_ZdcCoincidenceRate = -9999;
        //Primaries
        g_prim_Charge.clear();
        g_prim_Px.clear(); g_prim_Py.clear(); g_prim_Pz.clear(); g_prim_Pt.clear();
        g_prim_DCA.clear(); g_prim_dEdx.clear();
        //Towers
        g_tow_Id.clear();
        g_tow_Energy.clear(); g_tow_Et.clear();
        g_tow_Eta.clear(); g_tow_Phi.clear(); g_tow_EtaCorrected.clear(); g_tow_PhiCorrected.clear();
    
	suspect_pt.clear(); suspect_eta.clear(); suspect_phi.clear(); suspect_Q.clear();
        
	p_part_pt.clear(); //temp
	g_part_pt.clear();
	p_part_eta.clear();
	g_part_eta.clear();

	py_part_eta.clear();
	
	//temp
	p_UE_avgpt = -9999;
	p_UE_npart = -9999;
	g_UE_avgpt = -9999;
	g_UE_npart = -9999;
	
      mc_weight = -9999;
      p_n_jets = -9999; p_wt = -9999;
      p_conspT.clear();
      p_jetMult.clear();
      p_Pt.clear(); p_Eta.clear(); p_Phi.clear(); p_M.clear(); p_E.clear();
      p_ch_e_frac.clear();
      p_zg.clear(); p_rg.clear(); p_mg.clear(); p_ptg.clear();
      p_mcd.clear();
      
      event_weight = -9999;
      g_n_jets = -9999; g_wt = -9999;
      g_conspT.clear();
      g_jetMult.clear();
      g_Pt.clear(); g_Eta.clear(); g_Phi.clear(); g_M.clear(); g_E.clear();
      g_ch_e_frac.clear();
      g_zg.clear(); g_rg.clear(); g_mg.clear(); g_ptg.clear();
      g_mcd.clear();
      
      p_Particles_nocuts.clear(); g_Particles_nocuts.clear();
      p_Particles.clear(); g_Particles.clear();
      p_JetsInitial.clear(); g_JetsInitial.clear();
      p_Jets_nocuts.clear(); g_Jets_nocuts.clear();

      if (P6Reader->ReadEvent(p_EventID) == 1) {p_n_accepted++;} // == 1 => event loaded and passed selections
      if (P6Reader->ReadEvent(p_EventID) != 1) {cout << "DEBUG: SHOULDN'T SEE THIS? BAD PYTHIA EVENTS AREN'T READ, CORRECT?" << endl;}
      P6Reader->PrintStatus(10);     // Print out reader status every 10 seconds
      
      //filling the data structures that were defined before the event loop:
      p_event = P6Reader->GetEvent();
      p_header = p_event->GetHeader();  
      p_container = P6Reader->GetOutputContainer();
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~Skipping undesired events!~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                                                         
      pythiaFilename =  P6Reader->GetInputChain()->GetCurrentFile()->GetName();
      p_wt = LookupRun15Xsec( pythiaFilename );
     
      mc_weight = p_wt; //arbitrarily setting it to pythia's but can be either

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~QA~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      debug_events_Py[2] ++;
      
        p_EventId = p_header->GetEventId();
        p_RunId = p_header->GetRunId();
        //p_RefMult = p_header->GetReferenceMultiplicity();
        p_NOfGlobalTracks = p_header->GetNGlobalTracks();
        p_NOfTowers = p_header->GetNOfTowers();
        p_NOfPrimaryTracks = p_header->GetNOfPrimaryTracks();
        p_PVx = p_header->GetPrimaryVertexX();
        p_PVy = p_header->GetPrimaryVertexY();
        p_PVz = p_header->GetPrimaryVertexZ();
        p_vpdVz = p_header->GetVpdVz();
        p_BbcEastRate = p_header->GetBbcEastRate();
        p_BbcCoincidenceRate = p_header->GetBbcCoincidenceRate();
        p_BbcAdcSumEast = p_header->GetBbcAdcSumEast();
        
        //if (p_NOfPrimaryTracks > 0) {
        //for (int i = 0; i < p_NOfPrimaryTracks; ++ i) {}

	unsigned n_trks_in_evt = 0;
	//now looping over the tracks which pass selections (including hardcoded eta and pT cuts)
	
	TList *p_prim_trks = P6Reader->GetListOfSelectedTracks();
	//p_prim_trks = p_event->GetPrimaryTracks();
	TIter nextptrk(p_prim_trks);
	
        while (TStarJetPicoPrimaryTrack* p_track = (TStarJetPicoPrimaryTrack *) nextptrk()) {
	  double track_pt = (double) sqrt(p_track->GetPx()*p_track->GetPx() + p_track->GetPy()*p_track->GetPy());
	  //the charge requirement is due to PYTHIA not having towers and therefore storing tracks and towers together
	  if ( fabs(p_track->GetEta()) < 1 && track_pt > 0.2 && track_pt < 30.0 && p_track->GetCharge() != 0 && p_track->GetCharge() != -9999) {
	    p_prim_Charge.push_back(p_track->GetCharge());
	    p_prim_Px.push_back(p_track->GetPx());
	    p_prim_Py.push_back(p_track->GetPy());
	    p_prim_Pz.push_back(p_track->GetPz());
	    p_prim_Pt.push_back(track_pt);
	    p_prim_DCA.push_back(p_track->GetDCA());
	    p_prim_dEdx.push_back(p_track->GetdEdx());
	  }
	  else if (fabs(p_track->GetEta()) < 1 && track_pt > 0.2 && track_pt < 30.0 && p_track->GetCharge() == 0) {
	    p_tow_Eta.push_back(p_track->GetEta());
	    p_tow_Phi.push_back(p_track->GetPhi());
	    p_tow_Et.push_back(track_pt);//obviously this is slightly off but no way to access a track's Et without PID, right?
	  }
        }
	/*
	if (p_NOfTowers > 0) {
	  p_tows = p_event->GetTowers();
	  TIter nextptow(p_tows);
	  while (TStarJetPicoTower* p_tow = (TStarJetPicoTower *) nextptow()) {
	    if (p_tow->GetEt() > 0.2 && p_tow->GetEt() < 30.0 ) {
	      p_tow_Id.push_back(p_tow->GetId());
	      p_tow_Et.push_back(p_tow->GetEt());
	      p_tow_Eta.push_back(p_tow->GetEta());
	      p_tow_Phi.push_back(p_tow->GetPhi());
	      p_tow_EtaCorrected.push_back(p_tow->GetEtaCorrected());
	      p_tow_PhiCorrected.push_back(p_tow->GetPhiCorrected());
	    }
	    }
	    }
	*/
	//p_QATree->Fill(); //once per event
	  
	  
	// if we are not in matched mode, this token tells us that although the Pythia event is fine, the Geant event is bad and we shouldn't fill its tree later
	bool token_Ge_bad = 0;
	
	
	GEANTReader->ReadEvent(event); 
	g_EventID = GEANTReader->GetNOfCurrentEvent();
	if (GEANTReader->ReadEvent(g_EventID) == 1) {g_n_accepted++;}
	if (GEANTReader->ReadEvent(g_EventID) != 1) {token_Ge_bad = 1; cout << "BAD EVENT!" << endl;}
	
	if (token_Ge_bad == 0) {debug_events_Ge[2] ++;} //good geant event
	
	GEANTReader->PrintStatus(10);
	g_event = GEANTReader->GetEvent();
	g_header = g_event->GetHeader();
	g_container = GEANTReader->GetOutputContainer();
	geantFilename =  GEANTReader->GetInputChain()->GetCurrentFile()->GetName();
	g_wt = LookupRun15Xsec( geantFilename );

	//if the event lacks the desired trigger, skip it
	//TEMP! NOT REQUIRING TRIGGER RIGHT NOW. Makes it MB-like for comparison to PYTHIA. CHANGE LATER!
	
	//see function "SetTriggers()" for assignment of tID1, tID2 (or above until I write it)
	if ( ! (g_header->HasTriggerId(tID1) || g_header->HasTriggerId(tID2) || g_header->HasTriggerId(tID3) ) ) {//cout << "DEBUG: skipping this event because it lacks appropriate triggers. Does it have trigger ID " << tID1 << "? " << header->HasTriggerId(tID1) << endl;
	  token_Ge_bad = 1;
	  if (match) {continue;} //biases the pythia during matched mode, since we drop pythia events where the matching geant isn't triggered
	}
	
	if (trigger.find("pA") != string::npos) {//removing some runs by hand in pA until we have bad run/tower lists                                            
	  //TEMPORARILY SKIPPING THESE RUNS for pA [should define these runIDs somewhere later so they're not magic numbers]
	  
	  if (g_header->GetRunId() >= 16142059 && g_header->GetRunId() <= 16149001) {
	    token_Ge_bad = 1;
	    if (match) {continue;}
	  }
	  //something weird happened to the towers in run 16135032 (and it looks like it started at the end of run 16135031), so excluding both                  
	  if (g_header->GetRunId() == 16135031 || g_header->GetRunId() == 16135032) {
	    token_Ge_bad = 1;
	    if (match) {continue;}
	  }
	  
	  //the event cuts don't check if the vzdiff is acceptable, so I have to hardcode it here. UPDATE: I believe Nick updated the eventstructuredAu to check this condition, so this line should be redundant now                                                                                                           
	  //      if (!EventCuts->IsVertexZDiffOK(event)) {cout << "DEBUG: shouldn't see this now!" << endl; if (match) { continue;}}   
	  //Above 64000 seems like detectors saturate (tower multiplicity explodes).             
	  if (g_header->GetBbcAdcSumEast() >= pAu_BBCE_ADC_sum_max) {
	    token_Ge_bad = 1;
	    if (match) {continue;}
	  }
	  //selecting an activity range for the response matrix and eventual unfolding (low)
	  if (lowEA && (g_header->GetBbcAdcSumEast() < lowEA_low || g_header->GetBbcAdcSumEast() > lowEA_high)) {
	    token_Ge_bad = 1;
	    if (match) {continue;}//total "centrality" is too low (high) if < (>)                                                                   
	  }
	  //selecting an activity range for the response matrix and eventual unfolding (high)
	  if ((!lowEA) && (g_header->GetBbcAdcSumEast() < highEA_low || g_header->GetBbcAdcSumEast() > highEA_high)) {
	    token_Ge_bad = 1;
	    if (match) {continue;}//total "centrality" is too low (high) if < (>)                                                                   
	  }
	  
	
	}
	
	debug_events_Py[3] ++;
	if (token_Ge_bad == 0) {
	  debug_events_Ge[3] ++;
	}
	
        g_EventId = g_header->GetEventId();
        g_RunId = g_header->GetRunId();
        //g_RefMult = g_header->GetReferenceMultiplicity();
        g_NOfGlobalTracks = g_header->GetNGlobalTracks();
        g_NOfTowers = g_header->GetNOfTowers();
        g_NOfPrimaryTracks = g_header->GetNOfPrimaryTracks();
        g_PVx = g_header->GetPrimaryVertexX();
        g_PVy = g_header->GetPrimaryVertexY();
        g_PVz = g_header->GetPrimaryVertexZ();
        g_vpdVz = g_header->GetVpdVz();
        g_BbcEastRate = g_header->GetBbcEastRate();
	g_ZdcCoincidenceRate = g_header->GetZdcCoincidenceRate();
        g_BbcCoincidenceRate = g_header->GetBbcCoincidenceRate();
        g_BbcAdcSumEast = g_header->GetBbcAdcSumEast();
        
        //if (g_NOfPrimaryTracks > 0) {
          //  cout << "so it's come to this..." << endl;
	  // cout << "nprims: " << g_NOfPrimaryTracks << endl;
        //for (int i = 0; i < g_NOfPrimaryTracks; ++ i) {}

	  
	TList *g_prim_trks = GEANTReader->GetListOfSelectedTracks();
	//g_prim_trks = g_event->GetPrimaryTracks();
        TIter nextgtrk(g_prim_trks);
	//cout << "nentries: " << g_prim_trks->GetEntriesFast() << endl;
	int countieboi = 0;
        while (TStarJetPicoPrimaryTrack* g_track = (TStarJetPicoPrimaryTrack *) nextgtrk()) {
	  double track_pt = (double) sqrt(g_track->GetPx()*g_track->GetPx() + g_track->GetPy()*g_track->GetPy());
	  //TEMP!! (testing effect of increased Nbin on observables)
	  /*
	  if (g_track->GetTofTime() != -999 && track_pt >= 0.5) {//yes it's -999 and not -9999
	    token_Ge_bad = 1;
	    break;
	  }
	  */
	  //TEMP!! ^^
	  if (fabs(g_track->GetEta()) < 1 && track_pt > 0.2 && track_pt < 30.0) {
	    //countieboi ++; cout << countieboi << "th track";
	    g_prim_Charge.push_back(g_track->GetCharge());
            g_prim_Px.push_back(g_track->GetPx());
            g_prim_Py.push_back(g_track->GetPy());
            g_prim_Pz.push_back(g_track->GetPz());
	    g_prim_Pt.push_back(track_pt);
            g_prim_DCA.push_back(g_track->GetDCA());
            g_prim_dEdx.push_back(g_track->GetdEdx());
            //cout << " has dEdx = " << g_track->GetdEdx() << endl;
	  }
        }
	//TEMP (catching the events with high pT tracks in the MB)
	if (match && token_Ge_bad == 1) { continue; }
	//}
    //if (g_NOfTowers > 0) {
    TList *g_tows = GEANTReader->GetListOfSelectedTowers();//GetTowers();
        TIter nextgtow(g_tows);
        while (TStarJetPicoTower* g_tow = (TStarJetPicoTower *) nextgtow()) {
	  if (fabs(g_tow->GetEta()) && g_tow->GetEt() > 0.2 && g_tow->GetEt() < 30.0 ) {
            g_tow_Id.push_back(g_tow->GetId());
            g_tow_Et.push_back(g_tow->GetEt());
            g_tow_Eta.push_back(g_tow->GetEta());
            g_tow_Phi.push_back(g_tow->GetPhi());
            g_tow_EtaCorrected.push_back(g_tow->GetEtaCorrected());
            g_tow_PhiCorrected.push_back(g_tow->GetPhiCorrected());
	  }
        }
        //}
	//if (token_Ge_bad == 0) {
	//g_QATree->Fill(); //once per event 
        //}
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
        
      
	//if we require matching, must have an event in both P+G and P
	if ( match && GEANTReader->ReadEvent(p_EventID) != 1 )  {
	  //cout << "no corresponding geant event...skipping event " << p_EventID <<endl;
	  continue;//goes to the next event
	}
	//sanity check:
	if (match && (p_EventID != g_EventID) ) { cerr << "ERROR: READING DIFFERENT EVENTS. EXITING." << endl; exit(1);}	
	if (match && (pythiaFilename != geantFilename)) {cout << "I" << endl; std::cerr << "FILES DON'T MATCH! EXITING." << std::endl; exit(1);}
	if (match && (p_wt != g_wt)) {std::cerr << "WEIGHTS DON'T MATCH! EXITING." << std::endl; exit(1);}
     
      // converts TStarJetVectors to PseudoJets, carefully assigning the proper particle mass in either case
      GatherParticles ( p_container, p_sv, p_Particles, full, 1, pdg); //Pythia; full = 0 => charged-only, 1 => ch+ne
      GatherParticles ( g_container, g_sv, g_Particles, full, 0, pdg); //GEANT

      /*      
      GatherParticles ( p_container, p_sv, p_Particles_nocuts, full, 1, pdg); //Pythia; full = 0 => charged-only, 1 => ch+ne
      GatherParticles ( g_container, g_sv, g_Particles_nocuts, full, 0, pdg); //GEANT
      */

      
      /*      
      if (event == 2820) {
	cout << "BEGIN DEBUG TIME: (event = 2820; run = " << p_RunId << ")" << endl;
	for (int i = 0; i < p_Particles_nocuts.size(); ++ i) {
	  if (p_Particles_nocuts[i].pt() < 0.2) {
	    cout << "4momentum: ";
	    for (int j = 0; j < 4; ++ j) {
	      cout << p_Particles_nocuts[i].four_mom()[j] << " ";
	    }
	    cout << endl;
	    cout << "pt, eta, phi, theta, m: ";
	    cout << p_Particles_nocuts[i].pt() << " " << p_Particles_nocuts[i].eta() << " " << p_Particles_nocuts[i].phi() << " " << p_Particles_nocuts[i].theta() << " " << p_Particles_nocuts[i].m();
	    cout << endl;
	  }
	}
	cout << "END DEBUG TIME" << endl;
      }
      */

      //pythia
      // applying particle-level cuts
      //p_Particles_nocuts = p_Particles;
      //cout << "part size: " << p_Particles_nocuts.size() << " " << p_Particles.size() << endl;
      // vector<PseudoJet> p_nocut_Particles = spart_test(p_Particles_nocuts);
      vector<PseudoJet> p_cut_Particles = spart(p_Particles);
      //cout << "part size dos: " << p_nocut_Particles.size() << " " << p_cut_Particles.size() << endl;
      //Clustering jets
      //cout << "seg fault here? Event " << event << endl;
      //cout << "curious: " << token_Ge_bad << endl;
      //ClusterSequence p_Cluster_nocut(p_nocut_Particles, jet_def);
      //cout << "not yet!" << endl;
      //p_Jets_nocuts = sorted_by_pt(p_Cluster_nocut.inclusive_jets());
      ClusterSequence p_Cluster(p_cut_Particles, jet_def);
      p_JetsInitial = sorted_by_pt(sjet_gen(p_Cluster.inclusive_jets()));
      vector<PseudoJet> p_Jets;
      //geant
      //      if (token_Ge_bad == 0) {
	// applying particle-level cuts
	//g_Particles_nocuts = g_Particles;
	//vector<PseudoJet> g_nocut_Particles = spart_test(g_Particles_nocuts);
	vector<PseudoJet> g_cut_Particles = spart(g_Particles);
	//Clustering jets
	//	ClusterSequence g_Cluster_nocut(g_nocut_Particles, jet_def);
	//g_Jets_nocuts = sorted_by_pt(g_Cluster_nocut.inclusive_jets());
	ClusterSequence g_Cluster(g_cut_Particles, jet_def);
	g_JetsInitial = sorted_by_pt(sjet_det(g_Cluster.inclusive_jets()));
	//}
      vector<PseudoJet> g_Jets;
      
      //calculate REFMULT!
      //calculating refMult (N_ch in |eta| < 0.5)                                                                                                               
      double p_refmult = 0;
      for (int i = 0; i < p_cut_Particles.size(); ++ i) {
	py_part_eta.push_back(p_cut_Particles[i].eta());
        if (p_cut_Particles[i].user_index() != 0 && p_cut_Particles[i].user_index() != -9999 && fabs(p_cut_Particles[i].eta()) < 0.5) {
          p_refmult ++;
        }
      }
      p_RefMult = p_refmult;
      double g_refmult = 0;
      for (int i = 0; i < g_cut_Particles.size(); ++ i) {
        if (g_cut_Particles[i].user_index() != 0 && g_cut_Particles[i].user_index() != -9999 && fabs(g_cut_Particles[i].eta()) < 0.5) {
          if (g_cut_Particles[i].pt() < 0.2) { cout << "uh oh, wrong ref mult?" << endl;}
	  g_refmult ++;
        }
      }
      g_RefMult = g_refmult;
      cout << "REFMULT TEST: " << g_header->GetReferenceMultiplicity() << " ?=? " << g_RefMult << endl;

      

      p_QATree->Fill(); //once per event
      if (token_Ge_bad == 0) {
	g_QATree->Fill(); //once per event 
      }
      
      



      //Implementing a neutral energy fraction cut (of 90% currently) on inclusive det-level jets
      p_Jets = p_JetsInitial; //just passing intermediate -> final vector (no NEF selection on Pythia)


      if (p_Jets.size() != 0) {
	nEvts_for_weight ++;
      }
      
      //debug_jets_Py[0] += p_Jets_nocuts.size();
      if (token_Ge_bad != 1) {
	//debug_jets_Ge[0] += g_Jets_nocuts.size();
      }
      debug_jets_Py[1] += p_JetsInitial.size();
      if (token_Ge_bad != 1) {
	debug_jets_Ge[1] += g_JetsInitial.size();
	ApplyNEFSelection(g_JetsInitial, g_Jets); // uncomment this later!
	//g_Jets = g_JetsInitial;
      }
      
      debug_jets_Py[2] += p_Jets.size();
      if (token_Ge_bad != 1) {
	debug_jets_Ge[2] += g_Jets.size();
      }
      if (token_Ge_bad == 1) { //zero out the event info & vector of jets so we don't have info from these bad events added to the trees
	g_Jets.clear();
	g_EventId = -9999;
	g_wt = -9999;
      }

      //get UE avg pT
      vector<fastjet::PseudoJet> g_UE, p_UE; 
      if (g_Jets.size() > 0) {
	GatherUE (g_Jets[0], g_container , g_UE );
      }
      if (p_Jets.size() > 0) {
	GatherUE (p_Jets[0], p_container, p_UE );
      }
      double g_avgUEpt = -9999;
      if (g_UE.size() != 0) {
	g_avgUEpt = 0; //prevents zero-pT entries for events with no particles in the away-region
	for (int i = 0; i < g_UE.size(); ++ i) {
	  g_avgUEpt += g_UE[i].pt();
	}
	g_avgUEpt /= (double) g_UE.size();
      }
      double p_avgUEpt = -9999;
      if (p_UE.size() != 0) {
	p_avgUEpt = 0; //prevents zero-pT entries for events with no particles in the away-region
	for (int i = 0; i < p_UE.size(); ++ i) {
	  p_avgUEpt += p_UE[i].pt();
	}
	p_avgUEpt /= (double) p_UE.size();
      }
      p_UE_avgpt = p_avgUEpt;
      g_UE_avgpt = g_avgUEpt;
      
      //get UE npart
      p_UE_npart = p_UE.size();
      g_UE_npart = g_UE.size();

      //at this point, all criteria have been met (except checking if the event is bad in DiscardEmbedEvent() below)
      //so we can start filling jet information e.g.:
      p_n_jets = p_Jets.size();
      g_n_jets = g_Jets.size();
      
      vector<PseudoJet> p_GroomedJets; vector<PseudoJet> g_GroomedJets;
      //loop over the jets which passed cuts, groom them, and add to a vector (sorted by pt of the original jet)
      
      //SD grooming doesn't drop jets failing the criteria, so size of ungroomed = size of groomed
      for (int i = 0; i < p_Jets.size(); ++ i) {
	p_GroomedJets.push_back(sd(p_Jets[i]));
      }
      for (int i = 0; i < g_Jets.size(); ++ i) {
	g_GroomedJets.push_back(sd(g_Jets[i]));
      }
      
      //if jet pT > 2*pt-hat bin upper edge, discard event
      if (DiscardEmbedEvent(pythiaFilename, p_Jets, g_Jets)) { counter_debug ++; continue; }
      
      vector<double> suspectPt; vector<double> suspectEta; vector<double> suspectPhi; vector<double> suspectQ;
      for (int i = 0; i < g_Jets.size(); ++ i) {
	if (g_Jets[i].pt() > 80) {
	  vector<PseudoJet> cons = g_Jets[i].constituents();
	  for (int j = 0; j < cons.size(); ++ j) {
	    suspectPt.push_back(cons[j].pt());
	    suspectEta.push_back(cons[j].eta());
	    suspectPhi.push_back(cons[j].phi());
	    suspectQ.push_back(cons[j].user_index());
	  }
	  break;
	}
      }
      suspect_pt.push_back(suspectPt);
      suspect_eta.push_back(suspectEta);
      suspect_phi.push_back(suspectPhi);
      suspect_Q.push_back(suspectQ);
      
      //for calculating charged energy fraction of the jets
      vector<double> pch_e_frac, gch_e_frac;
      cout << "DEBUG: jets size: " << p_Jets.size() << endl;
      for (int i = 0; i < p_Jets.size(); ++ i) {
	cout << "DEBUG: jet " << i << " size: " << p_Jets[i].constituents().size() << endl;
	//looping over constituents
	double ch_e = 0; double tot_e = 0;//names are misnomers here since we use pT, not E.
	vector<PseudoJet> cons = p_Jets[i].constituents();
	for (int j = 0; j < cons.size(); ++ j) {
	  if (cons[j].user_index() != 0 && cons[j].user_index() != -9999) {cout << "DEBUG: charged cons " << j << " has pT " << cons[j].pt() << endl; ch_e += cons[j].pt();}
	  if (cons[j].user_index() == 0) {cout << "DEBUG: neutral cons " << j << " has pT " << cons[j].pt() << endl;}
	  tot_e += cons[j].pt();
	}
	pch_e_frac.push_back(ch_e/(double)tot_e);
	cout << "DEBUG: ch e frac = " << /*ch_e/(double)tot_e*/ pch_e_frac[i] << endl;
      }
      for (int i = 0; i < g_Jets.size(); ++ i) {
	//looping over constituents
	double ch_e = 0; double tot_e = 0;//names are misnomers here since we use pT, not E.
	vector<PseudoJet> cons = g_Jets[i].constituents();
	for (int j = 0; j < cons.size(); ++ j) {
	  if (cons[j].user_index() != 0 && cons[j].user_index() != -9999) {ch_e += cons[j].pt();}
	  tot_e += cons[j].pt();
	}
	gch_e_frac.push_back(ch_e/(double)tot_e);
      }
      
      //for calculating collinear dropped mass - i.e. how much mass we lost to grooming
      vector<double> pmcd, gmcd;
      for (int i = 0; i < p_Jets.size(); ++ i) {
	double m2 = (p_Jets[i].m())*(p_Jets[i].m()); double gm2 = (p_GroomedJets[i].m())*(p_GroomedJets[i].m());
	double m_cd = (double) sqrt(m2 - gm2); if ((m2 - gm2) < 1e-11) {m_cd = 0;}
	pmcd.push_back(m_cd);
      }
      for (int i = 0; i < g_Jets.size(); ++ i) {
	double m2 = (g_Jets[i].m())*(g_Jets[i].m()); double gm2 = (g_GroomedJets[i].m())*(g_GroomedJets[i].m());
	double m_cd = (double) sqrt(m2 - gm2); if ((m2 - gm2) < 1e-11) {m_cd = 0;}
	gmcd.push_back(m_cd);
      }
      if (pmcd.size() != p_Jets.size()) {cerr << "DEBUG: shouldn't have that size of M_cd vector is different than n_jets!" << endl; exit(1);}
      if (gmcd.size() != g_Jets.size()) {cerr << "DEBUG: shouldn't have that size of M_cd vector is different than n_jets!" << endl; exit(1);}
      
      //first check that we have a jet at detector-level in this event which has passed all cuts by this point, and that the event is good, then assign the weight for this event to make the embedding more data-like:
      if (g_Jets.size() != 0 && token_Ge_bad == 0) {
	//int binnum = refmult_rat_lowEA->GetXaxis()->FindBin(g_RefMult);
	if (lowEA_low < g_BbcAdcSumEast && lowEA_high > g_BbcAdcSumEast) {
	  
	  double refmultweight = frefmult_rat_lowEA->Eval(g_RefMult);//refmult_rat_lowEA->GetBinContent(binnum);
	  if (refmultweight == 0) {refmultweight = 1;}
	  event_weight = (double) mc_weight * refmultweight;
	  
	}
	else if (highEA_low < g_BbcAdcSumEast && highEA_high > g_BbcAdcSumEast) {
	  
	  double refmultweight = frefmult_rat_highEA->Eval(g_RefMult);//refmult_rat_highEA->GetBinContent(binnum);
	  if (refmultweight == 0) {refmultweight = 1;}
	  event_weight = (double) mc_weight * refmultweight;
	  
	}  
	else {
	  // event_weight = mc_weight; //this won't affect the weighted response since it's outside the activity range we select
	  if (token_Ge_bad == 0) {
	    cout << "Should never see this: a good event that wasn't in the activity range we select. Exiting!" << endl;
	    exit(1);
	  }
	}
      }
      //TEMP! UNWEIGHTED RESPONSE!
      //event_weight = mc_weight;

      
      //We have two vectors to be filled with matched jets. If they aren't, when looping over pythia jets, we have misses. Same goes when looping over geant jets with fakes. And for matches, we just fill with however many entries there are in the matched vectors.
      //MatchJets returns a vector of pairs of indices (i,j). The first entry is the position of the jet to match, the second its match's position, the third the position of the next jet to match, the fourth its match's position, etc.
      //FakesandMisses returns a vector of indices (i) corresponding to the indices of misses or fakes from the original candidate vector.
      if (match) {
	if (iSyst == 0) {
	  for (int i = 0; i < p_Jets.size(); ++ i) {
	    pt_gen_match_plus_miss->Fill(p_Jets[i].pt(),mc_weight);
	  }
	}
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MATCHING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	std::vector<fastjet::PseudoJet> g_matches; std::vector<fastjet::PseudoJet> p_matches;
	std::vector<fastjet::PseudoJet> g_matches_for_fakes; std::vector<fastjet::PseudoJet> p_matches_for_fakes; //only used to determine fakes
	std::vector<fastjet::PseudoJet> fakes; std::vector<fastjet::PseudoJet> misses;
	std::vector<fastjet::PseudoJet> g_sd_matches; std::vector<fastjet::PseudoJet> p_sd_matches;
	std::vector<fastjet::PseudoJet> sd_fakes; std::vector<fastjet::PseudoJet> sd_misses;
	std::vector<int> match_indices;
	std::vector<int> miss_indices; std::vector<int> fake_indices;

	//matches & misses
	if (p_Jets.size() != 0) {
	  g_matches.clear(); p_matches.clear(); misses.clear();
	  match_indices.clear(); miss_indices.clear();
                
	  match_indices = MatchJets(g_Jets, p_Jets, g_matches, p_matches); //find matches
                
	  if (g_matches.size() != p_matches.size()) {std::cerr << "Somehow we have different-sized match vectors. Exiting!" <<std::endl; exit(1);}
                
	  if (g_matches.size() < p_Jets.size()) { //then we have misses
	    miss_indices = FakesandMisses(p_matches, p_Jets, misses); //find misses
	  }
	}
        
	//fakes
	if (g_Jets.size() != 0) {
	  //clear the vectors to be used for determination of fakes (jets we find in Geant that don't have a match in Pythia)
	  fakes.clear(); fake_indices.clear();
	  g_matches_for_fakes.clear(); p_matches_for_fakes.clear();
		
	  MatchJets(p_Jets, g_Jets, p_matches_for_fakes, g_matches_for_fakes);
                
	  if (g_matches_for_fakes.size() != p_matches_for_fakes.size()) {std::cerr << "Somehow we have different-sized match vectors. Exiting!" <<std::endl; exit(1);}
                
	  if (p_matches_for_fakes.size() < g_Jets.size()) { //then we have fakes
	    fake_indices = FakesandMisses(g_matches_for_fakes, g_Jets, fakes);
	  }
	}
        
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	//NOTE: we don't perform a separate matching for groomed jets, so just index them with the match/miss/fake indices
            
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~FILL RESPONSES/HISTS/TREES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//	
	double prior_adjust = 0, prior_adjust_g = 0;
	for (int i = 0; i < misses.size(); ++ i) {
	  pt_response->Miss(misses[i].pt(), mc_weight);//don't weight miss jets by the data/embedding weight
	  m_response->Miss(misses[i].m(), mc_weight);
	  mg_response->Miss(p_GroomedJets[miss_indices[i]].m(), mc_weight);
	  cout << "A" << endl;
	  m_pt_response_lax->Miss(misses[i].m(), misses[i].pt(), mc_weight);
	  m_pt_response_med->Miss(misses[i].m(), misses[i].pt(), mc_weight);
	  m_pt_response_strict->Miss(misses[i].m(), misses[i].pt(), mc_weight);
	  mg_pt_response->Miss(p_GroomedJets[miss_indices[i]].m(), p_Jets[miss_indices[i]].pt(), mc_weight);
	  
	  res_syst[iSyst]->Miss(misses[i].m(), misses[i].pt(), mc_weight);
	  cout << "B" << endl;
	    
	}//for loop over misses
	
	//temp!
	if (g_matches.size() > 0) {//equivalently, p_matches.size()
	  //temp
	  for (int i = 0; i < p_cut_Particles.size(); ++ i) {
	    p_part_pt.push_back(p_cut_Particles[i].pt());
	    p_part_eta.push_back(p_cut_Particles[i].eta());
	  }
	  for (int i = 0; i < g_cut_Particles.size(); ++ i) {
	    g_part_pt.push_back(g_cut_Particles[i].pt());
	    g_part_eta.push_back(g_cut_Particles[i].eta());
	  }
	}

	//IMPORTANT NOTE:
	//match_indices contains the indices of pairs of geant and pythia matched jets. So if we loop over a list of
	//e.g. geant matches, the index i of a given geant match will be 2*i+1 in the match_indices list.
	//I.e. if we want to access the 3rd geant matched jet (at index 2), we skip over the first 2 pairs of match indices,
	//and the corresponding pythia match, to index at position 5. When accessing pythia jets in the match_indices list,
	//since they come first in each pair, it is similar but with 2*i instead. 
	for (int i = 0; i < g_matches.size(); ++ i) { //g_matches.size == p_matches.size == 1/2 (match_indices.size())
	  cout << "g_matches size " << g_matches.size() << endl;
	  //matches should be at same index in respective vectors
	  //RESPONSES:
	  double jetweight_strict = 1;
	  double jetweight_med = 1;
	  cout << "a" << endl;
	  if (lowEA_low < g_BbcAdcSumEast && lowEA_high > g_BbcAdcSumEast) {
	    cout << "ai" << endl;
	    jetweight_med = frefmult_pt_rat_lowEA->Eval(g_matches[i].pt(), g_BbcAdcSumEast);
	    cout << "aii" << endl;
	  }
	  cout << "b" << endl;
	  if (highEA_low < g_BbcAdcSumEast && highEA_high > g_BbcAdcSumEast) {
	    cout << "bi" << endl;
	    jetweight_med = frefmult_pt_rat_highEA->Eval(g_matches[i].pt(), g_BbcAdcSumEast);
	    cout << "bii" << endl;
	  }
	  cout << "c" << endl;
	  
	  //if (iSyst == 1) {
	    int binnum = pt_m_rat_lowEA->FindBin(g_matches[i].pt(),g_matches[i].m());
	    cout << "d" << endl;
	    if (lowEA_low < g_BbcAdcSumEast && lowEA_high > g_BbcAdcSumEast) {
	      cout << "di" << endl;
	      jetweight_strict = pt_m_rat_lowEA->GetBinContent(binnum);
	      cout << "dii" << endl;
	    }
	    cout << "e" << endl;
	    if (highEA_low < g_BbcAdcSumEast && highEA_high > g_BbcAdcSumEast) {
	      cout << "ei" << endl;
	      jetweight_strict = pt_m_rat_highEA->GetBinContent(binnum);
	      cout << "eii" << endl;
	    }
	    cout << "f" << endl;
	    if (jetweight_med <= 0 ) {jetweight_med = 1;}//makes sure empty bins aren't messing with things
	    if (jetweight_strict <= 0) {jetweight_strict = 1;}//makes sure empty bins aren't messing with things
	    
	    //for unfolding stability:
	    jetweight_med += 0.00001;
	    jetweight_strict += 0.00001;
	    
	    cout << "LATER!" << endl;
	    //}
	  /*
	  //TEMP!vv
	  int binnumpt = pt_rat_lowEA->GetXaxis()->FindBin(g_matches[i].pt());
	  int binnumm = m1520_rat_lowEA->GetXaxis()->FindBin(g_matches[i].m());//binning should be the same in low and high, and diff. pT selections
	  if (lowEA_low < g_BbcAdcSumEast && lowEA_high > g_BbcAdcSumEast) {
	    if (g_matches[i].pt() >= 15 && g_matches[i].pt() < 20) { jetweight = pt_rat_lowEA->GetBinContent(binnumpt) * m1520_rat_lowEA->GetBinContent(binnumm);}
	    if (g_matches[i].pt() >= 20 && g_matches[i].pt() < 25) { jetweight = pt_rat_lowEA->GetBinContent(binnumpt) * m2025_rat_lowEA->GetBinContent(binnumm);}
	    if (g_matches[i].pt() >= 25 && g_matches[i].pt() < 30) { jetweight = pt_rat_lowEA->GetBinContent(binnumpt) * m2530_rat_lowEA->GetBinContent(binnumm);}
	    if (g_matches[i].pt() >= 30 && g_matches[i].pt() < 40) { jetweight = pt_rat_lowEA->GetBinContent(binnumpt) * m3040_rat_lowEA->GetBinContent(binnumm);}
	    if (g_matches[i].pt() >= 40) { jetweight = pt_rat_lowEA->GetBinContent(binnumpt) * m40up_rat_lowEA->GetBinContent(binnumm);}
	    
	  }
	  if (highEA_low < g_BbcAdcSumEast && highEA_high > g_BbcAdcSumEast) {   
	    if (g_matches[i].pt() >= 15 && g_matches[i].pt() < 20) { jetweight = pt_rat_highEA->GetBinContent(binnumpt) * m1520_rat_highEA->GetBinContent(binnumm);}
	    if (g_matches[i].pt() >= 20 && g_matches[i].pt() < 25) { jetweight = pt_rat_highEA->GetBinContent(binnumpt) * m2025_rat_highEA->GetBinContent(binnumm);}
	    if (g_matches[i].pt() >= 25 && g_matches[i].pt() < 30) { jetweight = pt_rat_highEA->GetBinContent(binnumpt) * m2530_rat_highEA->GetBinContent(binnumm);}
	    if (g_matches[i].pt() >= 30 && g_matches[i].pt() < 40) { jetweight = pt_rat_highEA->GetBinContent(binnumpt) * m3040_rat_highEA->GetBinContent(binnumm);}
	    if (g_matches[i].pt() >= 40) { jetweight = pt_rat_highEA->GetBinContent(binnumpt) * m40up_rat_highEA->GetBinContent(binnumm);}
	    
	  }
	  
	  if (jetweight == 0) { jetweight = 1; }
	  */
	      //TEMP!^
	  pt_response->Fill(g_matches[i].pt(), p_matches[i].pt(), event_weight);
	  m_response->Fill(g_matches[i].m(), p_matches[i].m(), event_weight);
	  mg_response->Fill(g_GroomedJets[match_indices[(2*i)+1]].m(), p_GroomedJets[match_indices[2*i]].m(), event_weight);
	  cout << "C" << endl;
	  m_pt_response_lax->Fill(g_matches[i].m(), g_matches[i].pt(), p_matches[i].m(), p_matches[i].pt(), event_weight);
	  m_pt_response_med->Fill(g_matches[i].m(), g_matches[i].pt(), p_matches[i].m(), p_matches[i].pt(), mc_weight*jetweight_med);
	  m_pt_response_strict->Fill(g_matches[i].m(), g_matches[i].pt(), p_matches[i].m(), p_matches[i].pt(), mc_weight*jetweight_strict);
	  mg_pt_response->Fill(g_GroomedJets[match_indices[(2*i)+1]].m(), g_Jets[match_indices[(2*i)+1]].pt(),
			       p_GroomedJets[match_indices[2*i]].m(), p_Jets[match_indices[2*i]].pt(), event_weight);
	 
	  res_syst[iSyst]->Fill(g_matches[i].m(), g_matches[i].pt(), p_matches[i].m(), p_matches[i].pt(), event_weight);
	  cout << "D" << endl;
	  //(matched) TREES:

	  //ungroomed
	  p_Pt.push_back(p_matches[i].pt());
	  p_M.push_back(p_matches[i].m());
	  p_Eta.push_back(p_matches[i].eta());
	  p_Phi.push_back(p_matches[i].phi());
	  p_E.push_back(p_matches[i].e());
	  p_jetMult.push_back(p_matches[i].constituents().size());
	  //groomed
	  p_mg.push_back(p_GroomedJets[match_indices[2*i]].m());
	  p_ptg.push_back(p_GroomedJets[match_indices[2*i]].pt());
	  p_zg.push_back(p_GroomedJets[match_indices[2*i]].structure_of<SD>().symmetry());
	  p_rg.push_back(p_GroomedJets[match_indices[2*i]].structure_of<SD>().delta_R());
	  //ungroomed
	  g_Pt.push_back(g_matches[i].pt());
	  g_M.push_back(g_matches[i].m());
	  g_Eta.push_back(g_matches[i].eta());
	  g_Phi.push_back(g_matches[i].phi());
	  g_E.push_back(g_matches[i].e());
          g_jetMult.push_back(g_matches[i].constituents().size());
	  //groomed
	  g_mg.push_back(g_GroomedJets[match_indices[(2*i)+1]].m());
	  g_ptg.push_back(g_GroomedJets[match_indices[(2*i)+1]].pt());
	  g_zg.push_back(g_GroomedJets[match_indices[(2*i)+1]].structure_of<SD>().symmetry());
	  g_rg.push_back(g_GroomedJets[match_indices[(2*i)+1]].structure_of<SD>().delta_R());
	  //assigning vectors of values we calculated earlier to the tree
	  p_mcd.push_back(pmcd[match_indices[(2*i)+1]]);
	  g_mcd.push_back(gmcd[match_indices[(2*i)+1]]);
	  p_ch_e_frac.push_back(pch_e_frac[match_indices[(2*i)+1]]);
	  g_ch_e_frac.push_back(gch_e_frac[match_indices[(2*i)+1]]);
	}//for loop over matches
            
	for (int i = 0; i < fakes.size(); ++ i) {

	  double jetweight_strict = 1;
	  double jetweight_med = 1;
	  
	  if (lowEA_low < g_BbcAdcSumEast && lowEA_high > g_BbcAdcSumEast) {
	    jetweight_med = frefmult_pt_rat_lowEA->Eval(fakes[i].pt(), g_BbcAdcSumEast);
	  }
	  if (highEA_low < g_BbcAdcSumEast && highEA_high > g_BbcAdcSumEast) {
	    jetweight_med = frefmult_pt_rat_highEA->Eval(fakes[i].pt(), g_BbcAdcSumEast);
	  }
	  
	  
	  //if (iSyst == 1) {
	    int binnum = pt_m_rat_lowEA->FindBin(fakes[i].pt(),fakes[i].m());
	    if (lowEA_low < g_BbcAdcSumEast && lowEA_high > g_BbcAdcSumEast) {
	      jetweight_strict = pt_m_rat_lowEA->GetBinContent(binnum);
	    }
	    if (highEA_low < g_BbcAdcSumEast && highEA_high > g_BbcAdcSumEast) {
	      jetweight_strict = pt_m_rat_highEA->GetBinContent(binnum);
	    }
	    if (jetweight_med <= 0 ) {jetweight_med = 1;}//makes sure empty bins aren't messing with things
	    if (jetweight_strict <= 0) {jetweight_strict = 1;}//makes sure empty bins aren't messing with things
	    //for unfolding stability:
	    jetweight_med += 0.00001;
	    jetweight_strict += 0.00001;


	    /*
	  //if (iSyst == 1) {
	  int binnum = pt_rat_lowEA->GetXaxis()->FindBin(fakes[i].pt(),fakes[i].m());
	    if (lowEA_low < g_BbcAdcSumEast && lowEA_high > g_BbcAdcSumEast) {
	      jetweight = pt_rat_lowEA->GetBinContent(binnum);
	    }
	    if (highEA_low < g_BbcAdcSumEast && highEA_high > g_BbcAdcSumEast) {
	      jetweight = pt_rat_highEA->GetBinContent(binnum);
	    }
	    if (jetweight == 0) {jetweight = 1;}//makes sure empty bins aren't messing with things
	    //}
	    */
	  //TEMP!
	  //TEMP!vv
	  /*
	  int binnumpt = pt_rat_lowEA->GetXaxis()->FindBin(fakes[i].pt());
	  int binnumm = m1520_rat_lowEA->GetXaxis()->FindBin(fakes[i].m());//binning should be the same in low and high, and diff. pT selections
	  if (lowEA_low < g_BbcAdcSumEast && lowEA_high > g_BbcAdcSumEast) {
	    if (fakes[i].pt() >= 15 && fakes[i].pt() < 20) { jetweight = pt_rat_lowEA->GetBinContent(binnumpt) * m1520_rat_lowEA->GetBinContent(binnumm);}
	    if (fakes[i].pt() >= 20 && fakes[i].pt() < 25) { jetweight = pt_rat_lowEA->GetBinContent(binnumpt) * m2025_rat_lowEA->GetBinContent(binnumm);}
	    if (fakes[i].pt() >= 25 && fakes[i].pt() < 30) { jetweight = pt_rat_lowEA->GetBinContent(binnumpt) * m2530_rat_lowEA->GetBinContent(binnumm);}
	    if (fakes[i].pt() >= 30 && fakes[i].pt() < 40) { jetweight = pt_rat_lowEA->GetBinContent(binnumpt) * m3040_rat_lowEA->GetBinContent(binnumm);}
	    if (fakes[i].pt() >= 40) { jetweight = pt_rat_lowEA->GetBinContent(binnumpt) * m40up_rat_lowEA->GetBinContent(binnumm);}
	    
	  }
	  if (highEA_low < g_BbcAdcSumEast && highEA_high > g_BbcAdcSumEast) {   
	    if (fakes[i].pt() >= 15 && fakes[i].pt() < 20) { jetweight = pt_rat_highEA->GetBinContent(binnumpt) * m1520_rat_highEA->GetBinContent(binnumm);}
	    if (fakes[i].pt() >= 20 && fakes[i].pt() < 25) { jetweight = pt_rat_highEA->GetBinContent(binnumpt) * m2025_rat_highEA->GetBinContent(binnumm);}
	    if (fakes[i].pt() >= 25 && fakes[i].pt() < 30) { jetweight = pt_rat_highEA->GetBinContent(binnumpt) * m2530_rat_highEA->GetBinContent(binnumm);}
	    if (fakes[i].pt() >= 30 && fakes[i].pt() < 40) { jetweight = pt_rat_highEA->GetBinContent(binnumpt) * m3040_rat_highEA->GetBinContent(binnumm);}
	    if (fakes[i].pt() >= 40) { jetweight = pt_rat_highEA->GetBinContent(binnumpt) * m40up_rat_highEA->GetBinContent(binnumm);}
	    
	  }
	  if (jetweight == 0) { jetweight = 1; }
	  */
	  //
          pt_response->Fake(fakes[i].pt(), event_weight);
	  m_response->Fake(fakes[i].m(), event_weight);
	  mg_response->Fake(g_GroomedJets[fake_indices[i]].m(), event_weight);
	  cout << "E" << endl;
	  m_pt_response_lax->Fake(fakes[i].m(), fakes[i].pt(), event_weight);
	  m_pt_response_med->Fake(fakes[i].m(), fakes[i].pt(), mc_weight*jetweight_med);
	  m_pt_response_strict->Fake(fakes[i].m(), fakes[i].pt(), mc_weight*jetweight_strict);
	  mg_pt_response->Fake(g_GroomedJets[fake_indices[i]].m(), g_Jets[fake_indices[i]].pt(), event_weight);
	  res_syst[iSyst]->Fake(fakes[i].m(), fakes[i].pt(), event_weight);
	  cout << "F" << endl;
	}//for loop over fakes
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      }//matching-required conditional
      
      else { //matching not required
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~FILL RESPONSES/HISTS/TREES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	//looping over pythia jets
	for (int i = 0; i < p_Jets.size(); ++ i) {
	  //ungroomed
	  p_Pt.push_back(p_Jets[i].pt());
	  p_M.push_back(p_Jets[i].m());
	  p_Eta.push_back(p_Jets[i].eta());
	  p_Phi.push_back(p_Jets[i].phi());
	  p_E.push_back(p_Jets[i].e());
	    
	  //looping over constituents
	  vector<double> cons_pt;
	  vector<PseudoJet> cons = p_Jets[i].constituents();
	  int p_jet_size = cons.size();
	  p_jetMult.push_back(p_jet_size);
	  for (int j = 0; j < p_jet_size; ++ j) {
	    cons_pt.push_back(cons[j].pt());
	  }
	  p_conspT.push_back(cons_pt);

	  //groomed
	  p_mg.push_back(p_GroomedJets[i].m());
	  p_ptg.push_back(p_GroomedJets[i].pt());
	  p_zg.push_back(p_GroomedJets[i].structure_of<SD>().symmetry());
	  p_rg.push_back(p_GroomedJets[i].structure_of<SD>().delta_R());  
	}//end pythia jet loop
	//assigning vectors of values calculated earlier to the tree
	p_mcd = pmcd;
	p_ch_e_frac = pch_e_frac;
	
	if (!token_Ge_bad) {//Geant event is okay -- redundant since we zero out the bad geant event's jet vectors
	for (int i = 0; i < g_Jets.size(); ++ i) {
	  //ungroomed
	  if (g_Jets[i].pt() > 100) {
	    cout << "DEBUG 6/29: " << g_Jets[i].pt() << " " << g_Jets[i].m() << " " << g_Jets[i].eta() << " " << g_Jets[i].phi() << " " << g_EventID << " " << g_RunId << " " << g_Jets[i].constituents().size() << endl;
	    cout << "and cons: ";
	    for (int j = 0; j < g_Jets[i].constituents().size(); ++ j) {
	      cout << g_Jets[i].constituents()[j].pt() << " " << g_Jets[i].constituents()[j].eta() << " " << g_Jets[i].constituents()[j].phi() << " " << g_Jets[i].constituents()[j].m() << " " << g_Jets[i].constituents()[j].user_index() << endl;
	    }
	  }
	  g_Pt.push_back(g_Jets[i].pt());
	  g_M.push_back(g_Jets[i].m());
	  g_Eta.push_back(g_Jets[i].eta());
	  g_Phi.push_back(g_Jets[i].phi());
	  g_E.push_back(g_Jets[i].e());
	    
	  //looping over constituents
	  vector<double> cons_pt;
	  vector<PseudoJet> cons = g_Jets[i].constituents();
	  int g_jet_size = cons.size();
	  g_jetMult.push_back(g_jet_size);
	  for (int j = 0; j < g_jet_size; ++ j) {
	    cons_pt.push_back(cons[j].pt());
	  }
	  g_conspT.push_back(cons_pt);

	  //groomed
	  g_mg.push_back(g_GroomedJets[i].m());
	  g_ptg.push_back(g_GroomedJets[i].pt());
	  g_zg.push_back(g_GroomedJets[i].structure_of<SD>().symmetry());
	  g_rg.push_back(g_GroomedJets[i].structure_of<SD>().delta_R());
	}//end geant jet loop
	//assigning vectors of values calculated earlier to the tree
	g_mcd = gmcd;
	g_ch_e_frac = gch_e_frac;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	}//geant event okay conditional
      }//matching-not-required conditional
      
      
      if (p_Jets.size() != 0 && iSyst == 0) {
	//if (token_Ge_bad == 1) { // earlier we zero out geant jets from bad events, so here we can fill as normal
	  
	//}
	if (!token_Ge_bad) { suspects->Fill(); }//TEMP. Debugging too-high-pT embedded jets
	eventTree->Fill(); //when !match, will fill sometimes with empty geant vectors in the geant branches
	temptree->Fill();//temp
      }
        
      p_NJets += p_Jets.size(); g_NJets += g_Jets.size(); // add jets to total
    }//event loop
    
    //~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ END EVENT LOOP! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
    
    fout->cd();
    
    cout << endl << endl << "For systematic variation " << iSyst << endl
	 << "Of " << p_n_accepted << " pythia events and " << g_n_accepted << " geant events" << endl
	 << p_NJets << " gen jets have been found" << endl
	 << g_NJets << " det jets have been found" << endl << endl
	 << "Discarded " << counter_debug << " events on grounds of the found jets being too much higher than the pT-hat range" << endl;
    
    cout << "DEBUG: " << endl;
    for (int i = 0; i < debug_events_Ge.size(); ++ i) {
      cout << "N_EVENTS: PY = " << debug_events_Py[i] << " GE = " << debug_events_Ge[i] << endl;
    }
    for (int i = 0; i < debug_jets_Ge.size(); ++ i) {
      cout << "N_JETS: PY = " << debug_jets_Py[i] << " GE = " << debug_jets_Ge[i] << endl;
    }

    if (iSyst == 0) {//nominal case, so write closure and regular responses
      //trees
      eventTree->Write();
      p_QATree->Write();
      g_QATree->Write();
      suspects->Write();
 
      temptree->Write();//temp

      if (match) {
	pt_gen_match_plus_miss->Write();
	//hists

	//responses
	for (int i = 0; i < res.size(); ++ i) {
	  res[i]->Write();
	}
	
      }//match conditional
      
      
      cout << endl << "nevts " << nEvts_for_weight << endl; //temp to get n-events for cross section weights
      
    }//iSyst==0 conditional
    res_syst[iSyst]->Write();
  }//systematic variation loop

  cout << endl << "Writing to:  " << fout->GetName() << endl;

  fout->Close();
  
  //free the tclonesarray memory -- I think this has already been done behind the scenes because it causes a seg fault
  //p_prim_trks->Clear("C");
  //g_prim_trks->Clear("C");
  //p_tows->Clear("C");
  //g_tows->Clear("C");
  
  
  return 0;
}//main
