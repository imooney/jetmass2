//!sim.cxx
//!Isaac Mooney, WSU - July 2019
//!This file will run the (pp) analysis on the STAR-tuned Pythia6 and Pythia6+Geant for the jet mass project.
//!It takes in the Picos, performs selections, clusters particles, performs selections on the resulting jets,
//!applies the Soft Drop grooming procedure to a copy of the jet population, fills jet trees, and writes to files.

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
    bool full = 1;  //full = 1 => ch+ne; full = 0 => ch only.
    bool match = 0; //match = 0 => no match between Pythia & Pythia+Geant events.
    
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
            if (arguments[2] == "ch") {full = 0;} else {full = 1;}
            if (arguments[3] == "matched") {match = 1;} else if (arguments[3] == "unmatched") {match = 0;} else {cerr << "Not a valid flag!" << endl; exit(1);}
            chainList         = arguments[4];
            
            //Printing settings:
            cout << "Outputting to: " << (outputDir+outFileName).c_str() << "\nSettings:\n" << arguments[2] << " jets;\n match pythia and geant? " << match << ";\n input file: " << chainList << "\n";
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
    
    //initialize both readers
    InitReader(P6Reader, P6Chain, nEvents, "All", truth_absMaxVz, truth_vZDiff, truth_evPtMax, truth_evEtMax, truth_evEtMin, truth_DCA, truth_NFitPts, truth_FitOverMaxPts, sim_maxEtTow, 0.9999, false, sim_badTowers, sim_bad_run_list);
    InitReader(GEANTReader, GEANTChain, nEvents, det_triggerString, det_absMaxVz, det_vZDiff, det_evPtMax, det_evEtMax, det_evEtMin, det_DCA, det_NFitPts, det_FitOverMaxPts, sim_maxEtTow, 0.9999, false, /*sim_badTowers, sim_bad_run_list);*/det_badTowers, dat_bad_run_list);

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
    double p_n_jets, p_wt;
    vector<vector<double> > p_conspT;
    vector<double> p_jetMult;
    vector<double> p_Pt; vector<double> p_Eta; vector<double> p_Phi; vector<double> p_M; vector<double> p_E;
    vector<double> p_ch_e_frac;
    vector<double> p_zg; vector<double> p_rg; vector<double> p_mg; vector<double> p_ptg;
    vector<double> p_mcd;
    
    int g_EventID;
    double g_n_jets, g_wt;
    vector<vector<double> > g_conspT;
    vector<double> g_jetMult;
    vector<double> g_Pt; vector<double> g_Eta; vector<double> g_Phi; vector<double> g_M; vector<double> g_E;
    vector<double> g_ch_e_frac;
    vector<double> g_zg; vector<double> g_rg; vector<double> g_mg; vector<double> g_ptg;
    vector<double> g_mcd;

    //tree to hold jet and constituent quantites
    TTree *eventTree = new TTree("event","event");
    eventTree->Branch("p_n_jets", &p_n_jets);
    eventTree->Branch("p_conspT",&p_conspT);
    eventTree->Branch("p_jetMult",&p_jetMult);
    eventTree->Branch("p_Pt", &p_Pt); eventTree->Branch("p_Eta",&p_Eta); eventTree->Branch("p_Phi",&p_Phi); eventTree->Branch("p_M",&p_M); eventTree->Branch("p_E",&p_E);
    eventTree->Branch("p_ch_e_frac", &p_ch_e_frac);
    eventTree->Branch("p_zg", &p_zg); eventTree->Branch("p_rg", &p_rg); eventTree->Branch("p_mg", &p_mg); eventTree->Branch("p_ptg",&p_ptg);
    eventTree->Branch("p_mcd",&p_mcd);
    eventTree->Branch("p_weight", &p_wt); eventTree->Branch("p_EventID", &p_EventID);
    
    eventTree->Branch("g_n_jets", &g_n_jets);
    eventTree->Branch("g_conspT",&g_conspT);
    eventTree->Branch("g_jetMult",&g_jetMult);
    eventTree->Branch("g_Pt", &g_Pt); eventTree->Branch("g_Eta",&g_Eta); eventTree->Branch("g_Phi",&g_Phi); eventTree->Branch("g_M",&g_M); eventTree->Branch("g_E",&g_E);
    eventTree->Branch("g_ch_e_frac", &g_ch_e_frac);
    eventTree->Branch("g_zg", &g_zg); eventTree->Branch("g_rg", &g_rg); eventTree->Branch("g_mg", &g_mg); eventTree->Branch("g_ptg",&g_ptg);
    eventTree->Branch("g_mcd",&g_mcd);
    eventTree->Branch("g_weight", &g_wt); eventTree->Branch("g_EventID", &g_EventID);

    //1D hists for closure test
    TH1D* sampleA_pt_gen = new TH1D("sampleA_pt_gen","",15,5,80);
    TH1D* sampleA_pt_det = new TH1D("sampleA_pt_det","",9,15,60);
    TH1D* sampleA_m_gen = new TH1D("sampleA_m_gen","",14,0,14);
    TH1D* sampleA_m_det = new TH1D("sampleA_m_det","",14,0,14);
    TH1D* sampleA_mg_gen = new TH1D("sampleA_mg_gen","",14,0,14);
    TH1D* sampleA_mg_det = new TH1D("sampleA_mg_det","",14,0,14);

    TH1D* sampleB_pt_gen = new TH1D("sampleB_pt_gen","",15,5,80);
    TH1D* sampleB_pt_det = new TH1D("sampleB_pt_det","",9,15,60);
    TH1D* sampleB_m_gen = new TH1D("sampleB_m_gen","",14,0,14);
    TH1D* sampleB_m_det = new TH1D("sampleB_m_det","",14,0,14);
    TH1D* sampleB_mg_gen = new TH1D("sampleB_mg_gen","",14,0,14);
    TH1D* sampleB_mg_det = new TH1D("sampleB_mg_det","",14,0,14);

    //2D hists for closure test
    TH2D* sampleA_m_pt_gen = new TH2D("sampleA_m_pt_gen","",14,0,14,15,5,80);
    TH2D* sampleA_m_pt_det = new TH2D("sampleA_m_pt_det","",14,0,14,9,15,60);
    TH2D* sampleA_mg_pt_gen = new TH2D("sampleA_mg_pt_gen","",14,0,14,15,5,80);
    TH2D* sampleA_mg_pt_det = new TH2D("sampleA_mg_pt_det","",14,0,14,9,15,60);
    
    TH2D* sampleB_m_pt_gen = new TH2D("sampleB_m_pt_gen","",14,0,14,15,5,80);
    TH2D* sampleB_m_pt_det = new TH2D("sampleB_m_pt_det","",14,0,14,9,15,60);
    TH2D* sampleB_mg_pt_gen = new TH2D("sampleB_mg_pt_gen","",14,0,14,15,5,80);
    TH2D* sampleB_mg_pt_det = new TH2D("sampleB_mg_pt_det","",14,0,14,9,15,60);

    //~~~counts histograms for bin_drop, later~~~//
    
    //1D hists for closure test
    TH1D* sampleA_pt_gen_counts = new TH1D("sampleA_pt_gen_counts","",15,5,80);
    TH1D* sampleA_pt_det_counts = new TH1D("sampleA_pt_det_counts","",9,15,60);
    TH1D* sampleA_m_gen_counts = new TH1D("sampleA_m_gen_counts","",14,0,14);
    TH1D* sampleA_m_det_counts = new TH1D("sampleA_m_det_counts","",14,0,14);
    TH1D* sampleA_mg_gen_counts = new TH1D("sampleA_mg_gen_counts","",14,0,14);
    TH1D* sampleA_mg_det_counts = new TH1D("sampleA_mg_det_counts","",14,0,14);

    TH1D* sampleB_pt_gen_counts = new TH1D("sampleB_pt_gen_counts","",15,5,80);
    TH1D* sampleB_pt_det_counts = new TH1D("sampleB_pt_det_counts","",9,15,60);
    TH1D* sampleB_m_gen_counts = new TH1D("sampleB_m_gen_counts","",14,0,14);
    TH1D* sampleB_m_det_counts = new TH1D("sampleB_m_det_counts","",14,0,14);
    TH1D* sampleB_mg_gen_counts = new TH1D("sampleB_mg_gen_counts","",14,0,14);
    TH1D* sampleB_mg_det_counts = new TH1D("sampleB_mg_det_counts","",14,0,14);

    //2D hists for closure test
    TH2D* sampleA_m_pt_gen_counts = new TH2D("sampleA_m_pt_gen_counts","",14,0,14,15,5,80);
    TH2D* sampleA_m_pt_det_counts = new TH2D("sampleA_m_pt_det_counts","",14,0,14,9,15,60);
    TH2D* sampleA_mg_pt_gen_counts = new TH2D("sampleA_mg_pt_gen_counts","",14,0,14,15,5,80);
    TH2D* sampleA_mg_pt_det_counts = new TH2D("sampleA_mg_pt_det_counts","",14,0,14,9,15,60);
    
    TH2D* sampleB_m_pt_gen_counts = new TH2D("sampleB_m_pt_gen_counts","",14,0,14,15,5,80);
    TH2D* sampleB_m_pt_det_counts = new TH2D("sampleB_m_pt_det_counts","",14,0,14,9,15,60);
    TH2D* sampleB_mg_pt_gen_counts = new TH2D("sampleB_mg_pt_gen_counts","",14,0,14,15,5,80);
    TH2D* sampleB_mg_pt_det_counts = new TH2D("sampleB_mg_pt_det_counts","",14,0,14,9,15,60);

    //~~~~~~//
    
    //hists for use in responses
    TH2D *pyMvPt = new TH2D("pyMvPt",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
    TH2D *geMvPt = new TH2D("geMvPt",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
    TH2D *pyMvPt_counts = new TH2D("pyMvPt_counts",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
    TH2D *geMvPt_counts = new TH2D("geMvPt_counts",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
    TH2D *pyZgvPt = new TH2D("pyZgvPt", ";z_{g};p_{T} [GeV/c]",20,0,1,15,5,80);
    TH2D *geZgvPt = new TH2D("geZgvPt", ";z_{g};p_{T} [GeV/c]",20,0,1,9,15,60);
    TH2D *pyRgvPt = new TH2D("pyRgvPt", ";R_{g};p_{T} [GeV/c]",20,0,1,15,5,80);
    TH2D *geRgvPt = new TH2D("geRgvPt", ";R_{g};p_{T} [GeV/c]",20,0,1,9,15,60);
    TH2D *pyPtgvPt = new TH2D("pyPtgvPt", ";p_{T,g} [GeV/c];p_{T} [GeV/c]",15,5,80,15,5,80);
    TH2D *gePtgvPt = new TH2D("gePtgvPt", ";p_{T,g} [GeV/c];p_{T} [GeV/c]",9,15,60,9,15,60);
    TH2D *pyMgvPt = new TH2D("pyMgvPt", ";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
    TH2D *geMgvPt = new TH2D("geMgvPt", ";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
    TH2D *pyMgvPt_counts = new TH2D("pyMgvPt_counts",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
    TH2D *geMgvPt_counts = new TH2D("geMgvPt_counts",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
    
    // 1D responses
    RooUnfoldResponse *pt_response = new RooUnfoldResponse(9,15,60,15,5,80,"pt_response","");
    RooUnfoldResponse *m_response = new RooUnfoldResponse(14,0,14,14,0,14,"m_response","");
    RooUnfoldResponse *zg_response = new RooUnfoldResponse(20,0,1,20,0,1, "zg_response","");
    RooUnfoldResponse *rg_response = new RooUnfoldResponse(20,0,1,20,0,1, "rg_response","");
    RooUnfoldResponse *ptg_response = new RooUnfoldResponse(9,15,60,15,5,80, "ptg_response","");
    RooUnfoldResponse *mg_response = new RooUnfoldResponse(14,0,14,14,0,14, "mg_response","");
    
    //responses for training in the closure test
    RooUnfoldResponse *sampleA_pt_response = new RooUnfoldResponse(9,15,60,15,5,80,"sampleA_pt_response","");
    RooUnfoldResponse *sampleA_m_response = new RooUnfoldResponse(14,0,14,14,0,14,"sampleA_m_response","");
    RooUnfoldResponse *sampleA_zg_response = new RooUnfoldResponse(20,0,1,20,0,1, "sampleA_zg_response","");
    RooUnfoldResponse *sampleA_rg_response = new RooUnfoldResponse(20,0,1,20,0,1, "sampleA_rg_response","");
    RooUnfoldResponse *sampleA_ptg_response = new RooUnfoldResponse(9,15,60,15,5,80, "sampleA_ptg_response","");
    RooUnfoldResponse *sampleA_mg_response = new RooUnfoldResponse(14,0,14,14,0,14, "sampleA_mg_response","");

    //responses for validation in the closure test (unfolding pseudo-data with response constructed from same sample)
    RooUnfoldResponse *sampleB_pt_response = new RooUnfoldResponse(9,15,60,15,5,80,"sampleB_pt_response","");
    RooUnfoldResponse *sampleB_m_response = new RooUnfoldResponse(14,0,14,14,0,14,"sampleB_m_response","");
    RooUnfoldResponse *sampleB_zg_response = new RooUnfoldResponse(20,0,1,20,0,1, "sampleB_zg_response","");
    RooUnfoldResponse *sampleB_rg_response = new RooUnfoldResponse(20,0,1,20,0,1, "sampleB_rg_response","");
    RooUnfoldResponse *sampleB_ptg_response = new RooUnfoldResponse(9,15,60,15,5,80, "sampleB_ptg_response","");
    RooUnfoldResponse *sampleB_mg_response = new RooUnfoldResponse(14,0,14,14,0,14, "sampleB_mg_response","");
    
    // 2D responses
    RooUnfoldResponse *m_pt_response = new RooUnfoldResponse(geMvPt, pyMvPt, "m_pt_response");
    RooUnfoldResponse *m_pt_response_counts = new RooUnfoldResponse(geMvPt_counts, pyMvPt_counts, "m_pt_response_counts");
    RooUnfoldResponse *zg_pt_response = new RooUnfoldResponse(geZgvPt, pyZgvPt, "zg_pt_response");
    RooUnfoldResponse *rg_pt_response = new RooUnfoldResponse(geRgvPt, pyRgvPt, "rg_pt_response");
    RooUnfoldResponse *ptg_pt_response = new RooUnfoldResponse(gePtgvPt, pyPtgvPt, "ptg_pt_response");
    RooUnfoldResponse *mg_pt_response = new RooUnfoldResponse(geMgvPt, pyMgvPt, "mg_pt_response");
    RooUnfoldResponse *mg_pt_response_counts = new RooUnfoldResponse(geMgvPt_counts, pyMgvPt_counts, "mg_pt_response_counts");
    
    //responses for training in the closure test
    RooUnfoldResponse *sampleA_m_pt_response = new RooUnfoldResponse(geMvPt, pyMvPt, "sampleA_m_pt_response");
    RooUnfoldResponse *sampleA_m_pt_response_counts = new RooUnfoldResponse(geMvPt_counts, pyMvPt_counts, "sampleA_m_pt_response_counts");
    RooUnfoldResponse *sampleA_zg_pt_response = new RooUnfoldResponse(geZgvPt, pyZgvPt, "sampleA_zg_pt_response");
    RooUnfoldResponse *sampleA_rg_pt_response = new RooUnfoldResponse(geRgvPt, pyRgvPt, "sampleA_rg_pt_response");
    RooUnfoldResponse *sampleA_ptg_pt_response = new RooUnfoldResponse(gePtgvPt, pyPtgvPt, "sampleA_ptg_pt_response");
    RooUnfoldResponse *sampleA_mg_pt_response = new RooUnfoldResponse(geMgvPt, pyMgvPt, "sampleA_mg_pt_response");
    RooUnfoldResponse *sampleA_mg_pt_response_counts = new RooUnfoldResponse(geMgvPt_counts, pyMgvPt_counts, "sampleA_mg_pt_response_counts");
    
    //responses for validation in the closure test (unfolding pseudo-data with response constructed from same sample)
    RooUnfoldResponse *sampleB_m_pt_response = new RooUnfoldResponse(geMvPt, pyMvPt, "sampleB_m_pt_response");
    RooUnfoldResponse *sampleB_m_pt_response_counts = new RooUnfoldResponse(geMvPt_counts, pyMvPt_counts, "sampleB_m_pt_response_counts");
    RooUnfoldResponse *sampleB_zg_pt_response = new RooUnfoldResponse(geZgvPt, pyZgvPt, "sampleB_zg_pt_response");
    RooUnfoldResponse *sampleB_rg_pt_response = new RooUnfoldResponse(geRgvPt, pyRgvPt, "sampleB_rg_pt_response");
    RooUnfoldResponse *sampleB_ptg_pt_response = new RooUnfoldResponse(gePtgvPt, pyPtgvPt, "sampleB_ptg_pt_response");
    RooUnfoldResponse *sampleB_mg_pt_response = new RooUnfoldResponse(geMgvPt, pyMgvPt, "sampleB_mg_pt_response");
    RooUnfoldResponse *sampleB_mg_pt_response_counts = new RooUnfoldResponse(geMgvPt_counts, pyMgvPt_counts, "sampleB_mg_pt_response_counts");
    
    
    RooUnfoldResponse *m_response1520 = new RooUnfoldResponse(14,0,14,14,0,14,"m_response1520","");
    RooUnfoldResponse *m_response2025 = new RooUnfoldResponse(14,0,14,14,0,14,"m_response2025","");
    RooUnfoldResponse *m_response2530 = new RooUnfoldResponse(14,0,14,14,0,14,"m_response2530","");
    RooUnfoldResponse *m_response3040 = new RooUnfoldResponse(14,0,14,14,0,14,"m_response3040","");
    RooUnfoldResponse *m_response4060 = new RooUnfoldResponse(14,0,14,14,0,14,"m_response4060","");
    
    //vectors of responses & hists for easy writing to file later
    std::vector<RooUnfoldResponse*> res = {pt_response,m_response,zg_response,rg_response,ptg_response,mg_response,m_pt_response,m_pt_response_counts,zg_pt_response,rg_pt_response,ptg_pt_response,mg_pt_response,mg_pt_response_counts,m_response1520,m_response2025,m_response2530,m_response3040,m_response4060};
    
    std::vector<RooUnfoldResponse*> sampleA_res = {sampleA_pt_response,sampleA_m_response,sampleA_zg_response,sampleA_rg_response,sampleA_ptg_response,sampleA_mg_response,sampleA_m_pt_response,sampleA_m_pt_response_counts,sampleA_zg_pt_response,sampleA_rg_pt_response,sampleA_ptg_pt_response,sampleA_mg_pt_response,sampleA_mg_pt_response_counts};
    
    std::vector<RooUnfoldResponse*> sampleB_res = {sampleB_pt_response,sampleB_m_response,sampleB_zg_response,sampleB_rg_response,sampleB_ptg_response,sampleB_mg_response,sampleB_m_pt_response,sampleB_m_pt_response_counts,sampleB_zg_pt_response,sampleB_rg_pt_response,sampleB_ptg_pt_response,sampleB_mg_pt_response,sampleB_mg_pt_response_counts};
    
    std::vector<TH1D*> sampleA_h1Ds = {sampleA_pt_gen,sampleA_pt_det,sampleA_m_gen,sampleA_m_det,sampleA_mg_gen,sampleA_mg_det,sampleA_pt_gen_counts,sampleA_pt_det_counts,sampleA_m_gen_counts,sampleA_m_det_counts,sampleA_mg_gen_counts,sampleA_mg_det_counts};
    std::vector<TH2D*> sampleA_h2Ds = {sampleA_m_pt_gen,sampleA_m_pt_det,sampleA_mg_pt_gen,sampleA_mg_pt_det,sampleA_m_pt_gen_counts,sampleA_m_pt_det_counts,sampleA_mg_pt_gen_counts,sampleA_mg_pt_det_counts};
    std::vector<TH1D*> sampleB_h1Ds = {sampleB_pt_gen,sampleB_pt_det,sampleB_m_gen,sampleB_m_det,sampleB_mg_gen,sampleB_mg_det,sampleB_pt_gen_counts,sampleB_pt_det_counts,sampleB_m_gen_counts,sampleB_m_det_counts,sampleB_mg_gen_counts,sampleB_mg_det_counts};
    std::vector<TH2D*> sampleB_h2Ds = {sampleB_m_pt_gen,sampleB_m_pt_det,sampleB_mg_pt_gen,sampleB_mg_pt_det,sampleB_m_pt_gen_counts,sampleB_m_pt_det_counts,sampleB_mg_pt_gen_counts,sampleB_mg_pt_det_counts};
    
    
    //defining the algorithm and radius parameter for clustering jets
    JetDefinition jet_def(antikt_algorithm, R);

    //Creating SoftDrop grooming object
    contrib::SoftDrop sd(Beta,z_cut,R0);
    cout << "SoftDrop groomer is: " << sd.description() << endl;
    
    //for later use looking up PDG masses using particle PID
    TDatabasePDG *pdg = new TDatabasePDG();

    //SELECTORS
    // Constituent selectors
    // ---------------------
    Selector select_track_rap = fastjet::SelectorAbsRapMax(max_track_rap);
    Selector select_lopt      = fastjet::SelectorPtMin( partMinPt );
    Selector select_loptmax   = fastjet::SelectorPtMax( partMaxPt );
    Selector spart = select_track_rap * select_lopt * select_loptmax;
    
    // Jet candidate selectors
    // -----------------------
    Selector select_jet_rap     = fastjet::SelectorAbsRapMax(max_rap);
    Selector select_det_jet_pt_min  = fastjet::SelectorPtMin( det_jet_ptmin );
    Selector select_gen_jet_pt_min = fastjet::SelectorPtMin( jet_ptmin );
    Selector select_jet_pt_max  = fastjet::SelectorPtMax( jet_ptmax );
    Selector select_det_jet_m_min = fastjet::SelectorMassMin( mass_min );
    Selector select_gen_jet_m_min = fastjet::SelectorMassMin( 0.0 );
    
    Selector sjet_gen = select_jet_rap && select_gen_jet_pt_min && select_jet_pt_max && select_gen_jet_m_min;
    Selector sjet_det = select_jet_rap && select_det_jet_pt_min && select_jet_pt_max && select_det_jet_m_min;

    // Particle containers & counters
    vector<PseudoJet> p_Particles, g_Particles, p_JetsInitial, g_JetsInitial;
    int p_n_accepted = 0, g_n_accepted = 0; int p_NJets = 0, g_NJets = 0;
    int counter_debug = 0;
    double mc_weight = -1;
    // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
    
    for (int event = 0; event < P6Chain->GetEntries(); ++ event) {
        P6Reader->ReadEvent(event);
        GEANTReader->ReadEvent(event);
        
        g_EventID = GEANTReader->GetNOfCurrentEvent();
        p_EventID = P6Reader->GetNOfCurrentEvent();
	//NOTE: even # events will be used for training (response), odd will be used for validation (pseudo-data)
        
        //clearing vectors; initializing variables to -9999
        mc_weight = -9999;
        p_n_jets = -9999; p_wt = -9999;
        p_conspT.clear();
        p_jetMult.clear();
        p_Pt.clear(); p_Eta.clear(); p_Phi.clear(); p_M.clear(); p_E.clear();
        p_ch_e_frac.clear();
        p_zg.clear(); p_rg.clear(); p_mg.clear(); p_ptg.clear();
        p_mcd.clear();
        
	g_n_jets = -9999; g_wt = -9999;
        g_conspT.clear();
        g_jetMult.clear();
        g_Pt.clear(); g_Eta.clear(); g_Phi.clear(); g_M.clear(); g_E.clear();
        g_ch_e_frac.clear();
        g_zg.clear(); g_rg.clear(); g_mg.clear(); g_ptg.clear();
        g_mcd.clear();
        
        p_Particles.clear(); g_Particles.clear();
        p_JetsInitial.clear(); g_JetsInitial.clear();

        //if we require matching, must have an event in both P+G and P
        if ( match && GEANTReader->ReadEvent(p_EventID) != 1 ) {
            /*cout << "no corresponding geant event...skipping event " << p_eventID <<endl;*/
            continue;//goes to the next event
        }
        //sanity check:
        if ( match && (p_EventID != g_EventID) ) { cerr << "ERROR: READING DIFFERENT EVENTS. EXITING." << endl; exit(1);}
        
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
        
        if (match && (pythiaFilename != geantFilename)) {std::cerr << "FILES DON'T MATCH! EXITING." << std::endl; exit(1);}
        
        p_wt = LookupRun12Xsec( pythiaFilename );
        g_wt = LookupRun12Xsec( geantFilename );
        if (match && (p_wt != g_wt)) {std::cerr << "WEIGHTS DON'T MATCH! EXITING." << std::endl; exit(1);}
        mc_weight = p_wt; //arbitrarily setting it to pythia's but can be either
        
        // converts TStarJetVectors to PseudoJets, carefully assigning the proper particle mass in either case
        GatherParticles ( p_container, p_sv, p_Particles, full, 1, pdg); //Pythia; full = 0 => charged-only, 1 => ch+ne
        GatherParticles ( g_container, g_sv, g_Particles, full, 0, pdg); //GEANT
        
        // applying particle-level cuts
        vector<PseudoJet> p_cut_Particles = spart(p_Particles);
        vector<PseudoJet> g_cut_Particles = spart(g_Particles);
        
        //Clustering jets
        ClusterSequence p_Cluster(p_cut_Particles, jet_def);
        ClusterSequence g_Cluster(g_cut_Particles, jet_def);
        p_JetsInitial = sorted_by_pt(sjet_gen(p_Cluster.inclusive_jets()));
        g_JetsInitial = sorted_by_pt(sjet_det(g_Cluster.inclusive_jets()));
        
        vector<PseudoJet> p_Jets;
        vector<PseudoJet> g_Jets;
        
        //Implementing a neutral energy fraction cut (of 90% currently) on inclusive det-level jets
        p_Jets = p_JetsInitial; //just passing intermediate -> final vector (no NEF selection on Pythia)
        ApplyNEFSelection(g_JetsInitial, g_Jets);
        
        //at this point, all criteria have been met (except checking if the event is bad in DiscardEvent() below)
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
        
        if (DiscardEvent(pythiaFilename, p_Jets, g_Jets)) { counter_debug ++; continue; }
	
	//for calculating charged energy fraction of the jets
	vector<double> pch_e_frac, gch_e_frac;
	for (int i = 0; i < p_Jets.size(); ++ i) {
	  //looping over constituents
	  double ch_e = 0; double tot_e = 0;//names are misnomers here since we use pT, not E.
	  vector<PseudoJet> cons = p_Jets[i].constituents();
	  for (int j = 0; j < cons.size(); ++ j) {
	    if (cons[j].user_index() != 0) {ch_e += cons[j].pt();}
	    tot_e += cons[j].pt();
	  }
	  pch_e_frac.push_back(ch_e/(double)tot_e);
	}
	for (int i = 0; i < g_Jets.size(); ++ i) {
	  //looping over constituents
	  double ch_e = 0; double tot_e = 0;//names are misnomers here since we use pT, not E.
	  vector<PseudoJet> cons = g_Jets[i].constituents();
	  for (int j = 0; j < cons.size(); ++ j) {
	    if (cons[j].user_index() != 0) {ch_e += cons[j].pt();}
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

        //We have two vectors to be filled with matched jets. If they aren't, when looping over pythia jets, we have misses. Same goes when looping over geant jets with fakes. And for matches, we just fill with however many entries there are in the matched vectors.
        //MatchJets returns a vector of pairs of indices (i,j). The first entry is the position of the jet to match, the second its match's position, the third the position of the next jet to match, the fourth its match's position, etc.
        //FakesandMisses returns a vector of indices (i) corresponding to the indices of misses or fakes from the original candidate vector.
        
        if (match) {
	  //before the actual matching, fill the spectra for validation in the closure test
	  //(we only require there to be a matching Geant & Pythia event so our event reconstruction efficiency isn't folded unnecessarily into the closure)
	  if (p_EventID % 2 == 0) { //even events are sampleA
	    for (int i = 0; i < p_Jets.size(); ++ i) {
	      sampleA_pt_gen->Fill(p_Jets[i].pt(),mc_weight);
	      sampleA_m_gen->Fill(p_Jets[i].m(),mc_weight);
	      sampleA_mg_gen->Fill(p_GroomedJets[i].m(),mc_weight);
	      sampleA_m_pt_gen->Fill(p_Jets[i].m(),p_Jets[i].pt(),mc_weight);
	      sampleA_mg_pt_gen->Fill(p_GroomedJets[i].m(),p_Jets[i].pt(),mc_weight);
	      
	      sampleA_pt_gen_counts->Fill(p_Jets[i].pt());
	      sampleA_m_gen_counts->Fill(p_Jets[i].m());
	      sampleA_mg_gen_counts->Fill(p_GroomedJets[i].m());
	      sampleA_m_pt_gen_counts->Fill(p_Jets[i].m(),p_Jets[i].pt());
	      sampleA_mg_pt_gen_counts->Fill(p_GroomedJets[i].m(),p_Jets[i].pt());
	    }
	    for (int i = 0; i < g_Jets.size(); ++ i) {
	      sampleA_pt_det->Fill(g_Jets[i].pt(),mc_weight);
	      sampleA_m_det->Fill(g_Jets[i].m(),mc_weight);
	      sampleA_mg_det->Fill(g_GroomedJets[i].m(),mc_weight);
	      sampleA_m_pt_det->Fill(g_Jets[i].m(),g_Jets[i].pt(),mc_weight);
	      sampleA_mg_pt_det->Fill(g_GroomedJets[i].m(),g_Jets[i].pt(),mc_weight);

	      sampleA_pt_det_counts->Fill(g_Jets[i].pt());
	      sampleA_m_det_counts->Fill(g_Jets[i].m());
	      sampleA_mg_det_counts->Fill(g_GroomedJets[i].m());
	      sampleA_m_pt_det_counts->Fill(g_Jets[i].m(),g_Jets[i].pt());
	      sampleA_mg_pt_det_counts->Fill(g_GroomedJets[i].m(),g_Jets[i].pt());
	    }
	  }//end sampleA pseudo-data filling
	  if (p_EventID % 2 != 0) { //odd events are sampleB
	    for (int i = 0; i < p_Jets.size(); ++ i) {
	      sampleB_pt_gen->Fill(p_Jets[i].pt(),mc_weight);
	      sampleB_m_gen->Fill(p_Jets[i].m(),mc_weight);
	      sampleB_mg_gen->Fill(p_GroomedJets[i].m(),mc_weight);
	      sampleB_m_pt_gen->Fill(p_Jets[i].m(),p_Jets[i].pt(),mc_weight);
	      sampleB_mg_pt_gen->Fill(p_GroomedJets[i].m(),p_Jets[i].pt(),mc_weight);

	      sampleB_pt_gen_counts->Fill(p_Jets[i].pt());
	      sampleB_m_gen_counts->Fill(p_Jets[i].m());
	      sampleB_mg_gen_counts->Fill(p_GroomedJets[i].m());
	      sampleB_m_pt_gen_counts->Fill(p_Jets[i].m(),p_Jets[i].pt());
	      sampleB_mg_pt_gen_counts->Fill(p_GroomedJets[i].m(),p_Jets[i].pt());
	    }
	    for (int i = 0; i < g_Jets.size(); ++ i) {
	      sampleB_pt_det->Fill(g_Jets[i].pt(),mc_weight);
	      sampleB_m_det->Fill(g_Jets[i].m(),mc_weight);
	      sampleB_mg_det->Fill(g_GroomedJets[i].m(),mc_weight);
	      sampleB_m_pt_det->Fill(g_Jets[i].m(),g_Jets[i].pt(),mc_weight);
	      sampleB_mg_pt_det->Fill(g_GroomedJets[i].m(),g_Jets[i].pt(),mc_weight);

	      sampleB_pt_det_counts->Fill(g_Jets[i].pt());
	      sampleB_m_det_counts->Fill(g_Jets[i].m());
	      sampleB_mg_det_counts->Fill(g_GroomedJets[i].m());
	      sampleB_m_pt_det_counts->Fill(g_Jets[i].m(),g_Jets[i].pt());
	      sampleB_mg_pt_det_counts->Fill(g_GroomedJets[i].m(),g_Jets[i].pt());
	    }
	  }//end sampleB pseudo-data filling
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
                //clear the vectors to be used for determination of fakes (jets we find in Geant that don't have a match in Pythia
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
            
            for (int i = 0; i < misses.size(); ++ i) {
	      m_pt_response->Miss(misses[i].m(), misses[i].pt(), mc_weight);
	      mg_pt_response->Miss(p_GroomedJets[miss_indices[i]].m(), p_Jets[miss_indices[i]].pt(), mc_weight);
	      //closure sampleA responses
	      if (match && p_EventID % 2 == 0) { //throughout this section, checking if (match) is unnecessary because of the overall conditional, but I'll leave it.
		sampleA_pt_response->Miss(misses[i].pt(), mc_weight);
		sampleA_m_response->Miss(misses[i].m(), mc_weight);
		sampleA_m_pt_response->Miss(misses[i].m(), misses[i].pt(), mc_weight);
		sampleA_mg_pt_response->Miss(p_GroomedJets[miss_indices[i]].m(), p_Jets[miss_indices[i]].pt(), mc_weight);
	      }
	      //closure sampleB responses
	      if (match && p_EventID % 2 != 0) {
		sampleB_pt_response->Miss(misses[i].pt(), mc_weight);
		sampleB_m_response->Miss(misses[i].m(), mc_weight);
		sampleB_m_pt_response->Miss(misses[i].m(), misses[i].pt(), mc_weight);
		sampleB_mg_pt_response->Miss(p_GroomedJets[miss_indices[i]].m(), p_Jets[miss_indices[i]].pt(), mc_weight);
	      }
            }//for loop over misses
            
            for (int i = 0; i < g_matches.size(); ++ i) { //g_matches.size == p_matches.size == match_indices.size()
	      //matches should be at same index in respective vectors
	      //RESPONSES:
	      m_pt_response->Fill(g_matches[i].m(), g_matches[i].pt(), p_matches[i].m(), p_matches[i].pt(), mc_weight);
	      mg_pt_response->Fill(g_GroomedJets[match_indices[i]].m(), g_Jets[match_indices[i]].pt(),
				   p_GroomedJets[match_indices[i]].m(), p_Jets[match_indices[i]].pt());
	      //closure sampleA responses
	      if (match && p_EventID % 2 == 0) {
		sampleA_pt_response->Fill(g_matches[i].pt(), p_matches[i].pt(), mc_weight);
		sampleA_m_response->Fill(g_matches[i].m(), p_matches[i].m(), mc_weight);
		sampleA_m_pt_response->Fill(g_matches[i].m(), g_matches[i].pt(), p_matches[i].m(), p_matches[i].pt(), mc_weight);
		sampleA_mg_pt_response->Fill(g_GroomedJets[match_indices[i]].m(), g_Jets[match_indices[i]].pt(),
					     p_GroomedJets[match_indices[i]].m(), p_Jets[match_indices[i]].pt());
	      }
	      //closure sampleB responses
	      if (match && p_EventID % 2 != 0) {
		sampleB_pt_response->Fill(g_matches[i].pt(), p_matches[i].pt(), mc_weight);
		sampleB_m_response->Fill(g_matches[i].m(), p_matches[i].m(), mc_weight);
		sampleB_m_pt_response->Fill(g_matches[i].m(), g_matches[i].pt(), p_matches[i].m(), p_matches[i].pt(), mc_weight);
		sampleB_mg_pt_response->Fill(g_GroomedJets[match_indices[i]].m(), g_Jets[match_indices[i]].pt(),
					     p_GroomedJets[match_indices[i]].m(), p_Jets[match_indices[i]].pt());
	      }
	      
              
	      //(matched) TREES:
	      //ungroomed
	      p_Pt.push_back(p_matches[i].pt());
	      p_M.push_back(p_matches[i].m());
	      p_Eta.push_back(p_matches[i].eta());
	      p_Phi.push_back(p_matches[i].phi());
	      p_E.push_back(p_matches[i].e());
	      //groomed
	      p_mg.push_back(p_GroomedJets[match_indices[i]].m());
	      p_ptg.push_back(p_GroomedJets[match_indices[i]].pt());
	      p_zg.push_back(p_GroomedJets[match_indices[i]].structure_of<SD>().symmetry());
	      p_rg.push_back(p_GroomedJets[match_indices[i]].structure_of<SD>().delta_R());
	      //ungroomed
	      g_Pt.push_back(g_matches[i].pt());
	      g_M.push_back(g_matches[i].m());
	      g_Eta.push_back(g_matches[i].eta());
	      g_Phi.push_back(g_matches[i].phi());
	      g_E.push_back(g_matches[i].e());
	      //groomed
	      g_mg.push_back(g_GroomedJets[match_indices[i]].m());
	      g_ptg.push_back(g_GroomedJets[match_indices[i]].pt());
	      g_zg.push_back(g_GroomedJets[match_indices[i]].structure_of<SD>().symmetry());
	      g_rg.push_back(g_GroomedJets[match_indices[i]].structure_of<SD>().delta_R());
	      
	      //assigning vectors of values we calculated earlier to the tree
	      p_mcd.push_back(pmcd[match_indices[i]]);
	      g_mcd.push_back(gmcd[match_indices[i]]);
	      
	      p_ch_e_frac.push_back(pch_e_frac[match_indices[i]]);
	      g_ch_e_frac.push_back(gch_e_frac[match_indices[i]]);
	      
            }//for loop over matches
            
            for (int i = 0; i < fakes.size(); ++ i) {
	      m_pt_response->Fake(fakes[i].m(), fakes[i].pt(), mc_weight);
	      mg_pt_response->Fake(g_GroomedJets[fake_indices[i]].m(), g_Jets[fake_indices[i]].pt());
	      //closure sampleA responses
	      if (match && p_EventID % 2 == 0) {
		sampleA_pt_response->Fake(fakes[i].pt(), mc_weight);
		sampleA_m_response->Fake(fakes[i].m(), mc_weight);
		sampleA_m_pt_response->Fake(fakes[i].m(), fakes[i].pt(), mc_weight);
		sampleA_mg_pt_response->Fake(g_GroomedJets[fake_indices[i]].m(), g_Jets[fake_indices[i]].pt());
	      }
	      //closure sampleB responses
	      if (match && p_EventID % 2 != 0) {
		sampleB_pt_response->Fake(fakes[i].pt(), mc_weight);
		sampleB_m_response->Fake(fakes[i].m(), mc_weight);
		sampleB_m_pt_response->Fake(fakes[i].m(), fakes[i].pt(), mc_weight);
		sampleB_mg_pt_response->Fake(g_GroomedJets[fake_indices[i]].m(), g_Jets[fake_indices[i]].pt());
	      }
	      
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
	  
	  for (int i = 0; i < g_Jets.size(); ++ i) {
	    //ungroomed
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
        }//matching-not-required conditional
	
        if (p_Jets.size() != 0) {
	  eventTree->Fill(); //when !match, will fill sometimes with empty geant vectors in the geant branches
        }
        
        p_NJets += p_Jets.size(); g_NJets += g_Jets.size(); // add jets to total
    }//event loop
    
    //~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ END EVENT LOOP! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ //
    
    TFile *fout = new TFile( ( outputDir + outFileName ).c_str() ,"RECREATE");
    
    cout << endl << endl << "Of " << p_n_accepted << " pythia events and " << g_n_accepted << " geant events" << endl;
    cout << p_NJets << " gen jets have been found" << endl;
    cout << g_NJets << " det jets have been found" << endl << endl;
    cout << endl << "Writing to:  " << fout->GetName() << endl << endl;
    cout << "Discarded " << counter_debug << " events on grounds of the found jets being too much higher than the pT-hat range" << endl;
    
    //trees
    eventTree->Write();
    
    if (match) {
      //hists
      for (int i = 0; i < sampleA_h1Ds.size(); ++ i) {
	sampleA_h1Ds[i]->Write();
      }
      for (int i = 0; i < sampleA_h2Ds.size(); ++ i) {
	sampleA_h2Ds[i]->Write();
      }
      for (int i = 0; i < sampleB_h1Ds.size(); ++ i) {
	sampleB_h1Ds[i]->Write();
      }
      for (int i = 0; i < sampleB_h2Ds.size(); ++ i) {
	sampleB_h2Ds[i]->Write();
      } 
      
      //responses
      for (int i = 0; i < res.size(); ++ i) {
	res[i]->Write();
      }
      for (int i = 0; i < sampleA_res.size(); ++ i) {
	sampleA_res[i]->Write();
      }
      for (int i = 0; i < sampleB_res.size(); ++ i) {
	sampleB_res[i]->Write();
      }
      
    }//match conditional
    
    fout->Close();

    return 0;
}//main
