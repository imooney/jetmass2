//! Isaac Mooney, WSU - June, 2019OAOAOA
//! Jet mass analysis at RHIC kinematics

#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FinalPartons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "HepMC/GenParticle.h"
#include "HepMC/GenEvent.h"
//!ROOT headers
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"

#include "TDatabasePDG.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>

double m_chpion = 0.13957018; //!GeV/c^2

namespace Rivet {

  class STAR_MASS : public Analysis {
  public:
    
    /// Constructor
    STAR_MASS(string name = "STAR_MASS")
      : Analysis(name)
    {
      setNeedsCrossSection(true);
      _mode=0;
    }
    
    //! http://journals.aps.org/prd/pdf/10.1103/PhysRevD.91.111501
    //! Anti kT jet is reclustered using CA and min(pT1, pT2)/(pT1+pT2) > z_cut * (R12/Rj)^\beta
    //! for sub-jets that satisfy this condition, z_g = min(pT1,pT2)/(pT1+pT2)
    std::pair<double,double> SoftDrop(const fastjet::PseudoJet &j, double jetR, double z_cut, double beta)
    {
      //! give the soft drop groomer a short name
      //! Use a symmetry cut z > z_cut R^beta
      //! By default, there is no mass-drop requirement
      fastjet::contrib::RecursiveSymmetryCutBase::SymmetryMeasure symmetry_measure = fastjet::contrib::RecursiveSymmetryCutBase::scalar_z;
      fastjet::contrib::SoftDrop sd(beta, z_cut, symmetry_measure, jetR);
      PseudoJet sd_jet = sd(j);    
      if(sd_jet != 0){
	double z = sd_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry();
	double r = sd_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R();
	return std::pair<double,double>(z,r);
      }else 
	return std::pair<double,double>(-1,-1);
    }

    std::vector<int> MatchJets(const std::vector<fastjet::PseudoJet> candidates_safe, const std::vector<fastjet::PseudoJet> toMatch, std::vector<fastjet::PseudoJet> & c_matches, std::vector<fastjet::PseudoJet> & t_matches) {
      std::vector<int> match_indices;
      if (candidates_safe.size() == 0 || toMatch.size() == 0) {
	return match_indices; //later, match_indices being empty will tell us there were no matches                                                                      
      }
      
      //define candidates outside the loop so list continually dwindles as we remove matched candidates                                                                   
      std::vector<fastjet::PseudoJet> candidates = candidates_safe;
      for (unsigned i = 0; i < toMatch.size(); ++ i) { //for each jet in toMatch, we try to find a match from candidates_copy                                           
	//defined inside the loop so that for each toMatch jet there's a new set of candidates                                                                            
	std::vector<fastjet::PseudoJet> candidates_copy = candidates;
	fastjet::Selector selectMatchedJets = fastjet::SelectorCircle( _jetR );
	selectMatchedJets.set_reference( toMatch[i] );
	//note: matchedToJet and candidates_copy are equivalent, assuming candidates_safe was already sorted by pT                                                        
	std::vector<fastjet::PseudoJet> matchedToJet = sorted_by_pt( selectMatchedJets( candidates_copy ));
	if (matchedToJet.size() == 0) { continue; } //means no match to this jet. Remove none from candidates. Continuing on to the next one.                             
	else { //found at least one match. Need to remove the highest pT one from candidates and add the respective jets to the match vectors.                            
	  match_indices.push_back(i); //push back the toMatch match position                                                                                             
	  t_matches.push_back(toMatch[i]);
	  c_matches.push_back(matchedToJet[0]); //highest pT match                                                                                                        
	  for (unsigned j = 0; j < candidates.size(); ++ j) { //finding which one to delete from candidates before next toMatch iteration.
	    if (matchedToJet[0].delta_R(candidates[j]) < 0.0001) { //is probably the same jet                                                                            
	      candidates.erase(candidates.begin() + j); //removing the jet from the overall list of candidates so it can't be considered next time                        
	      match_indices.push_back(j); //push back the candidate match position                                                                                        
	      break; //should exit only the c_matches loop.                                                                                                              
	    }
	  }
	}
      }
      return match_indices;
    }
       
    //! Book histograms and initialise projections before the run
    void init() {

      //! ~~~~~~~~~~~~~~~Read in flags/parameters from a file~~~~~~~~~~~~~~~~~~~ //
      int nParams = 3; //change whenever more parameters are added
      py_input = true; //to determine whether we're analyzing pythia or herwig events
      std::string sim_name = "pythia8"; //default
      double jetRad = 0.4; //default
      std::string radiusText = "04";
      std::string decays = "undecayed"; //default
      std::string line;
      ifstream flag_file ("rivet_flags.txt");
      if (flag_file.is_open()) {
	for (int i = 0; i < 2*nParams; ++i) {//Each parameter has a comment above, so reading n parameters means iterating over 2n lines.
	  getline (flag_file,line);
	  if (line.find("inputIsPythia") != std::string::npos) {
	    if (line.find("TRUE") != std::string::npos) {
	      std::cout << "ANALYZING PYTHIA EVENTS!" << std::endl;
	      py_input = true;
	      sim_name = "pythia8";
	    }
	    else if (line.find("FALSE") != std::string::npos) {
	      std::cout << "ANALYZING HERWIG EVENTS!" << std::endl;
	      py_input = false;
	      sim_name = "herwig7";
	    }
	    else {
	      std::cerr << "Error. Cannot determine MC type to analyze. Exiting." << std::endl;
	      exit(1);
	    }
	  }
	  else if (line.find("jetRadius") != std::string::npos) {
	    //will now have to do some string magic to turn e.g. 04 into 0.4. 
	    //by the way, this should treat R > 1 correctly until R >= 10 [however the jet eta cut of 1 - R renders R > 1 meaningless]
	    std::size_t pos = line.find("=");
	    if (line.find("=") == std::string::npos) {
	      std::cerr << "Error. Cannot determine the jet radius for clustering. Exiting." << std::endl;
	    }
	    //finds the location of e.g. 04 in the line, then converts it to a double: 0.4
	    radiusText = line.substr(pos+2);//pos is location of =. Then there's a space, then the number. So pos+2 is position of e.g. 0 in 04.
	    std::cout << "DEBUG: " << radiusText << std::endl;
	    std::string dot = ".";
	    std::cout << "DEBUG: " << radiusText.substr(0,1) << " " << radiusText.substr(1,1) << std::endl;
	    std::string radiusNum = radiusText.substr(0,1)+dot+radiusText.substr(1,1); //e.g. 0.4
	    std::cout << "DEBUG: " << radiusNum << std::endl;
	    jetRad = (double) std::stod(radiusNum); //converting the string 0.4 to a double
	    std::cout << "DEBUG: radius text = " << radiusText << " and double = " << jetRad << std::endl;
	  }
	  else if (line.find("decayFlag") != std::string::npos) {
	    int pos = line.find("=");
	    if (line.find("=") == std::string::npos) {
	      std::cerr << "Error. Cannot determine final state information. Exiting." << std::endl;
            }
	    decays = line.substr(pos+2);
	    std::cout << "DEBUG: decays = " << decays << std::endl; 
	  }
	}
	flag_file.close();
      }
      else {
	std::cerr << "Unable to open file. Exiting." << std::endl;
	exit(1);
      }

      //! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
      
      //! for later use looking up PDG masses using particle PID
      pdg = new TDatabasePDG();

      //! defaults
      //_ptedges = {5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 60.0};
      //_nptbins = _ptedges.size()-1;
      _pTCut = 0.2;
      _pTMax = 1000.0;
      _pTconsMax = 30.0;
      _max_track_rap = 1.0;
      
      //! jet eta max is calculated from the radius given.
      //RADIUS = 4;
      _jetR = jetRad;//(float)RADIUS/10;

      //! debug flag 
      printDebug = false;
      
      //! parameters for splitting function
      z_cut=0.1;
      beta = 0;      
      //! final state and jet projection
      //! generic final state
      FinalState fs(-1.0, 1.0, 0.*GeV);
      addProjection(fs, "FS");
      FastJets fj(fs, FastJets::ANTIKT, _jetR);
      addProjection(fj, "Jets");
    
      //! skips hadronization
      FinalPartons fp(Cuts::abseta < 1.0);
      addProjection(fp, "FP");
      FastJets pj(fp, FastJets::ANTIKT, _jetR);
      addProjection(pj, "PartonJets");

      //handles the output for each pt-hat bin
      if(_mode == 0){
	//! For P6: "../out/py6_decayed_jewel_pthatbin580_R%s.root"
	fout = new TFile(("../out/star_mass_"+sim_name+"_IsaacsHepmcs_FSR_only_5pthatbin10_R"+radiusText+"_"+decays+".root").c_str(),"RECREATE");
	fout->cd();
      }    
      else if(_mode == 1){
	fout = new TFile(("../out/star_mass_"+sim_name+"_IsaacsHepmcs_FSR_only_10pthatbin15_R"+radiusText+"_"+decays+".root").c_str(),"RECREATE");
	fout->cd();
      }else if(_mode == 2){
	fout = new TFile(("../out/star_mass_"+sim_name+"_IsaacsHepmcs_FSR_only_15pthatbin20_R"+radiusText+"_"+decays+".root").c_str(),"RECREATE");
	fout->cd();
      }else if(_mode == 3){
	fout = new TFile(("../out/star_mass_"+sim_name+"_IsaacsHepmcs_FSR_only_20pthatbin25_R"+radiusText+"_"+decays+".root").c_str(),"RECREATE");
	fout->cd();
      }else if(_mode == 4){
	fout = new TFile(("../out/star_mass_"+sim_name+"_IsaacsHepmcs_FSR_only_25pthatbin30_R"+radiusText+"_"+decays+".root").c_str(),"RECREATE");
	fout->cd();
      }else if(_mode == 5){
	fout = new TFile(("../out/star_mass_"+sim_name+"_IsaacsHepmcs_FSR_only_30pthatbin35_R"+radiusText+"_"+decays+".root").c_str(),"RECREATE");
	fout->cd();
      }else if(_mode == 6){
	fout = new TFile(("../out/star_mass_"+sim_name+"_IsaacsHepmcs_FSR_only_35pthatbin40_R"+radiusText+"_"+decays+".root").c_str(),"RECREATE");
	fout->cd();
      }else if(_mode == 7){
	fout = new TFile(("../out/star_mass_"+sim_name+"_IsaacsHepmcs_FSR_only_40pthatbin45_R"+radiusText+"_"+decays+".root").c_str(),"RECREATE");
	fout->cd();
      }else if(_mode == 8){
	fout = new TFile(("../out/star_mass_"+sim_name+"_IsaacsHepmcs_FSR_only_45pthatbin50_R"+radiusText+"_"+decays+".root").c_str(),"RECREATE");
	fout->cd();
      }else if(_mode == 9){
	fout = new TFile(("../out/star_mass_"+sim_name+"_IsaacsHepmcs_FSR_only_50pthatbin60_R"+radiusText+"_"+decays+".root").c_str(),"RECREATE");
	fout->cd();
      }else if(_mode == 10){
	fout = new TFile(("../out/star_mass_"+sim_name+"_IsaacsHepmcs_FSR_only_60pthatbin80_R"+radiusText+"_"+decays+".root").c_str(),"RECREATE");
	fout->cd();
      }
        
      //! book desired histograms here
      //...

      //! partonic final state jets
      PartonTree=new TTree("PartonTree","Parton Jets");
      //! clear the vectors for each event!
      PLpt.clear();
      PLm.clear();
      PLptg.clear();
      PLmg.clear();
      PLeta.clear();
      PLphi.clear();
      PLzg.clear();
      PLrg.clear();
      PLncons.clear();
      PLnconsg.clear();
      PLconsM.clear();
      PLconspT.clear();
      PLconsEta.clear();
      
      pthat = -9999;
      mcweight = -9999;
      
      PartonTree->Branch("PLpt",&PLpt);
      PartonTree->Branch("PLm",&PLm);
      PartonTree->Branch("PLptg",&PLptg);
      PartonTree->Branch("PLmg",&PLmg);
      PartonTree->Branch("PLeta",&PLeta);
      PartonTree->Branch("PLphi",&PLphi);
      PartonTree->Branch("PLzg",&PLzg);
      PartonTree->Branch("PLrg",&PLrg);
      PartonTree->Branch("PLncons",&PLncons);
      PartonTree->Branch("PLnconsg",&PLnconsg);
      PartonTree->Branch("PLconsM",&PLconsM);
      PartonTree->Branch("PLconspT",&PLconspT);
      PartonTree->Branch("PLconsEta",&PLconsEta);
      PartonTree->Branch("pthat",&pthat,"pthat/d");
      PartonTree->Branch("mcweight",&mcweight,"mcweight/d");

      //! partonic final state jets matched geometrically to hadronic final state jets
      MatchTree = new TTree("MatchTree","Matched partonic/hadronic jets");
      PL_pt_match.clear();
      PL_m_match.clear();
      PL_mg_match.clear();
      PL_zg_match.clear();
      HL_pt_match.clear();
      HL_m_match.clear();
      HL_mg_match.clear();
      HL_zg_match.clear();

      MatchTree->Branch("PL_pt_match",&PL_pt_match);
      MatchTree->Branch("PL_m_match",&PL_m_match);
      MatchTree->Branch("PL_mg_match",&PL_mg_match);
      MatchTree->Branch("PL_zg_match",&PL_zg_match);
      MatchTree->Branch("HL_pt_match",&HL_pt_match);
      MatchTree->Branch("HL_m_match",&HL_m_match);
      MatchTree->Branch("HL_mg_match",&HL_mg_match);
      MatchTree->Branch("HL_zg_match",&HL_zg_match);
      MatchTree->Branch("pthat",&pthat,"pthat/d");
      MatchTree->Branch("mcweight",&mcweight,"mcweight/d");


      //! hadronic final state jets
      ResultTree=new TTree("ResultTree","Result Jets");
      //! clear the vectors for each event! 
      conspT.clear();
      consDist.clear();
      consGirth.clear();
      consM.clear();
      consEta.clear();
      jetTau1.clear();
      jetTau05.clear();
      jetTau0.clear();
      jetTau_05.clear();
      jetTau_1.clear();
      jetGirth.clear();
      qvg.clear();
      jetpT.clear();
      jetM.clear();
      jetMult.clear();
      sdjetpT.clear();
      sdjetM.clear();
      jeteta.clear();
      jetphi.clear();
      zg.clear();
      rg.clear();
      nCons.clear();
      
      ResultTree->Branch("mcweight", &mcweight,"mcweight/d");
      ResultTree->Branch("pthat", &pthat,"pthat/d");
      ResultTree->Branch("eventID",&eventID,"eventID/d");
      ResultTree->Branch("nCons",&nCons);
      ResultTree->Branch("conspT",&conspT);
      ResultTree->Branch("consDist",&consDist);
      ResultTree->Branch("consGirth",&consGirth);
      ResultTree->Branch("consM",&consM);
      ResultTree->Branch("consEta",&consEta);
      ResultTree->Branch("jetTau1",&jetTau1);
      ResultTree->Branch("jetTau05",&jetTau05);
      ResultTree->Branch("jetTau0",&jetTau0);
      ResultTree->Branch("jetTau_05",&jetTau_05);
      ResultTree->Branch("jetTau_1",&jetTau_1);
      ResultTree->Branch("jetGirth",&jetGirth);
      ResultTree->Branch("jetMult",&jetMult);
      ResultTree->Branch("qvg",&qvg);
      ResultTree->Branch("jetpT", &jetpT);
      ResultTree->Branch("jeteta", &jeteta);
      ResultTree->Branch("jetphi", &jetphi);
      ResultTree->Branch("jetM", &jetM);
      ResultTree->Branch("sdjetpT", &sdjetpT);
      ResultTree->Branch("sdjetM", &sdjetM);
      ResultTree->Branch("zg", &zg);
      ResultTree->Branch("rg", &rg);
      
    }//! init

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = handler().crossSection();

      //! clear the vectors for each event!
      
      //parton jets
      PLpt.clear();
      PLm.clear();
      PLptg.clear();
      PLmg.clear();
      PLeta.clear();
      PLphi.clear();
      PLzg.clear();
      PLrg.clear();
      PLncons.clear();
      PLnconsg.clear();
      PLconsM.clear();
      PLconsEta.clear();
      PLconspT.clear();
      
      //matched jets
      PL_pt_match.clear();
      PL_m_match.clear();
      PL_mg_match.clear();
      PL_zg_match.clear();
      HL_pt_match.clear();
      HL_m_match.clear();
      HL_mg_match.clear();
      HL_zg_match.clear();
      
      //hadron jets
      nCons.clear();
      conspT.clear();
      consDist.clear();
      consGirth.clear();
      consM.clear();
      consEta.clear();
      jetTau1.clear();
      jetTau05.clear();
      jetTau0.clear();
      jetTau_05.clear();
      jetTau_1.clear();
      jetGirth.clear();
      jetMult.clear();
      qvg.clear();
      jetpT.clear();
      jetM.clear();
      sdjetpT.clear();
      sdjetM.clear();
      jeteta.clear();
      jetphi.clear();
      zg.clear();
      rg.clear();
      
      pthat = -9999;
      mcweight = -9999;
      
      const double qscale = event.genEvent()->event_scale();
      //SHOULD BE CHANGED IF A DIFFERENT NUMBER OF EVENTS ARE RUN!
      int nEvents = -1;
      const char * fname = fout->GetName();
      std::string fstr = fname;
      if (fstr.find("py6") != std::string::npos) {
	nEvents = 10000000;
      }
      else {
	nEvents = 1000000;//100; for testing
      }
      eventID = event.genEvent()->event_number() + _mode*nEvents;
      // std::cout<<"pthat = "<<qscale<<" ; xsec-weight = "<<weight<<std::endl;

      // Constituent selectors
      // ---------------------
      fastjet::Selector select_track_rap = fastjet::SelectorAbsRapMax(_max_track_rap);
      fastjet::Selector select_lopt      = fastjet::SelectorPtMin( _pTCut );
      fastjet::Selector select_loptmax   = fastjet::SelectorPtMax( _pTconsMax );
      fastjet::Selector spart = select_track_rap * select_lopt * select_loptmax;

      fastjet::Selector select_eta  = fastjet::SelectorAbsEtaMax(1.0 - _jetR);
      fastjet::Selector select_pt_hi   = fastjet::SelectorPtMin(5.0);
      //fastjet::Selector select_pt_hi_max = fastjet::SelectorPtMax(_pTMax);
      fastjet::Selector select_pt_lo   = fastjet::SelectorPtMin(0.001);
      fastjet::Selector select_both_hi = select_pt_hi && select_eta;
      fastjet::Selector select_both_lo = select_pt_lo && select_eta;
      const ParticleVector& FS = applyProjection<FinalState>(event, "FS").particlesByPt();     
      const ParticleVector& FP = applyProjection<FinalPartons>(event, "FP").particlesByPt();
      
      Cut cuts = Cuts::etaIn(-1, 1) & (Cuts::pT > _pTCut*GeV); 

      fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, _jetR);
      
      //grooming                                                                                                                                                         
      fastjet::contrib::RecursiveSymmetryCutBase::SymmetryMeasure  symmetry_measure = fastjet::contrib::RecursiveSymmetryCutBase::scalar_z;
      fastjet::contrib::SoftDrop sd(beta, z_cut, symmetry_measure, _jetR);

      //QUARK v. GLUON PREP: figuring out what/where the hard scattered partons are
      vector<int> identity; vector<PseudoJet> hardparton;
      //for (int i = 0; i < event.genEvent()->particles_size(); ++ i) {
      foreach (const Particle &p, event.allParticles()) {
	if (!py_input) {
	  std::cout << "WARNING! Quark v. gluon identification is not a feature of the Herwig analysis. Do not use output pertaining to q v. g discrimination!" << std::endl;
	}
	if (/*event.genEvent()->statusHepMC(i)*/ p.genParticle()->status() == 23) {//PYTHIA only! NEED TO CHANGE FOR HERWIG?
	  if (p.pid() <=6 && p.pid() >= 1) {
	    identity.push_back(0);//0 MEANS QUARK!
	  }
	  else if (p.pid() == 21) {
	    identity.push_back(1);//1 MEANS GLUON!
	  }
	  else {
	    identity.push_back(-1);//MEANS NEITHER!
	  }
	  hardparton.push_back(PseudoJet(p.px(), p.py(), p.pz(), p.E()));	  
	}
      }
      
      //debugging!
      /*
      cout << "identity: " << identity.size() << " hardparton: " << hardparton.size() << std::endl;
      if ((int) identity.size() != (int) hardparton.size()) { std::cerr << "Number of hard partons should match size of the list of the identities of said hard partons, but that does not seem to be the case. Exiting." << std::endl; exit(1);}
      std::cout << "NUMBER OF INITIAL HARD PARTONS: " << hardparton.size() << std::endl;
      std::cout << "THEIR IDENTITIES ARE: " << std::endl;
      for (int i = 0; i < (int) identity.size(); ++ i) {
	if (identity[i] == 0) {std::cout << "     quark!" << std::endl;}
	if (identity[i] == 1) {std::cout << "     gluon!" << std::endl;}
	if (identity[i] == -1) {std::cout << "     neither quark nor gluon!" << std::endl;}
      }
      */

      //! Get all the hadrons
      PseudoJets pJet_sub; PseudoJets pJet_ch; PseudoJet intermediate; 	      
      foreach ( const Particle& p, FS) {
	intermediate = PseudoJet(p.px(), p.py(), p.pz(), p.E());
	
	//cout << "DEBUG: final state hadron = " << p.pid() << ": " << p.mass() << " " << pdg->GetParticle(p.pid())->Mass() << endl;

	intermediate.reset_PtYPhiM(p.pt(),p.rap(),p.phi(),(double) pdg->GetParticle(p.pid())->Mass()); //PDG MASSES!!!
       
	pJet_sub.push_back(intermediate);
	//! add only charged particles to charged jets
	if (p.charge() != 0) {
	  pJet_ch.push_back(PseudoJet(p.px(), p.py(), p.pz(), p.E()));
	}
	
      }     
      //! get all the partons
      PseudoJets PLjet;
      foreach ( const Particle & p, FP) {
	intermediate = PseudoJet(p.px(), p.py(), p.pz(), p.E());

	//cout << setprecision(5) << "DEBUG: final state parton = " << p.pid() << ": " << p.mass() << " " << pdg->GetParticle(p.pid())->Mass() << endl;
	
	//! the following two lines are possibilities for adjusting the mass assignment per particle. By default, parton's off-shell mass is used.
	//intermediate.reset_PtYPhiM(p.pt(),p.rap(),p.phi(),pdg->GetParticle(p.pid())->Mass()); //PDG MASSES!!!
	//intermediate.reset_PtYPhiM(p.pt(),p.rap(),p.phi(),0); //MASSLESS!!!
	
        PLjet.push_back(intermediate);
      }
      
      //discard bad particles -> cluster remaining particles into jets -> select on the jets
      
      vector<PseudoJet> pJet_cut = spart(pJet_sub);
      vector<PseudoJet> pJet_ch_cut = spart(pJet_ch);
      vector<PseudoJet> PL_cut = spart(PLjet);
      
      fastjet::ClusterSequence cs_sub(pJet_cut, jet_def);
      fastjet::ClusterSequence cs_ch(pJet_ch_cut, jet_def);
      fastjet::ClusterSequence cs_PL(PL_cut, jet_def);

      PseudoJets recoJets_bs = sorted_by_pt(cs_sub.inclusive_jets());
      PseudoJets recoJets_bs_ch = sorted_by_pt(cs_ch.inclusive_jets());
      PseudoJets recoJets_bs_pl = sorted_by_pt(cs_PL.inclusive_jets());

      PseudoJets recoJets_bkgsub = select_both_hi(recoJets_bs);
      PseudoJets recoJets_bkgsub_ch = select_both_hi(recoJets_bs_ch);
      PseudoJets recoJets_PL = select_both_hi(recoJets_bs_pl);
      
      if(printDebug){
	std::cout<<"bkg subtracted Jets: "<<std::endl;
	foreach (const PseudoJet jet, recoJets_bkgsub){
	  std::cout<<"    "<<jet.pt()<<std::endl;
	}
      }

      mcweight = weight;
      pthat = qscale;

      //discards events on the grounds of them having jets of pT > double the high end of the pT-hat bin from which they came.
      bool bad_event = 0;
      int upper_edge;
      const char * str = fout->GetName();
      std::string s = str;
      if (_mode == 0 && s.find("py6") == std::string::npos) {upper_edge = 10;} if (_mode == 1) {upper_edge = 15;} if (_mode == 2) {upper_edge = 20;} if (_mode == 3) {upper_edge = 25;}
      if (_mode == 4) {upper_edge = 30;} if (_mode == 5) {upper_edge = 35;} if (_mode == 6) {upper_edge = 40;} if (_mode == 7) {upper_edge = 45;}
      if (_mode == 8) {upper_edge = 50;} if (_mode == 9) {upper_edge = 60;} if (_mode == 10){upper_edge = 80;}
      if (s.find("py6") != std::string::npos) {upper_edge = 80;}
      if (recoJets_bkgsub.size() != 0) {
	if (recoJets_bkgsub[0].pt() > 2*upper_edge) {
	  std::cout << "removing event from bin " << _mode << " due to bad jet with pt, eta, phi, and m: " << recoJets_bkgsub[0].pt() << " " << recoJets_bkgsub[0].eta() << " " << recoJets_bkgsub[0].phi() << " " << recoJets_bkgsub[0].m() << std::endl;
	  bad_event = 1;
	}
      }
      if (bad_event) {return;}
      
      
      //loop over good parton-level jets, filling vectors with observables
      foreach(const PseudoJet jet, recoJets_PL) {
	PLpt.push_back(jet.pt());
	PLm.push_back(jet.m());
	PLeta.push_back(jet.eta());
	PLphi.push_back(jet.phi());
	PLncons.push_back(jet.constituents().size());
	
	vector<double> cons_M; vector<double> cons_Pt; vector<double> cons_Eta;
        foreach (const PseudoJet cons, jet.constituents()) {
	  cons_M.push_back(cons.m());
	  cons_Eta.push_back(cons.eta());
	  cons_Pt.push_back(cons.pt());
	}
	PLconsM.push_back(cons_M);
	PLconspT.push_back(cons_Pt);
	PLconsEta.push_back(cons_Eta);

	PseudoJet sd_jet = sd(jet);
	
	if(sd_jet != 0) { //! unless I'm mistaken, this condition should never evaluate to false, since grooming by default doesn't remove jets 

          PLptg.push_back(sd_jet.pt());
          PLmg.push_back(sd_jet.m());
          double z = sd_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry();
          double r = sd_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R();

          PLzg.push_back(z);
          PLrg.push_back(r);
	  PLnconsg.push_back(sd_jet.constituents().size());

	}	
	else {
	  std::cout << "Somehow don't have a groomed jet for the given ungroomed jet..." << std::endl;
	}	
      }

      
      //! DEBUG: std::cout << "number of good jets in the event: " << recoJets_bkgsub.size() << std::endl;
      
      //looping over good hadron-level jets, filling vectors with observables
      foreach(const PseudoJet jet, recoJets_bkgsub){// CHANGE WHEN GOING BETWEEN FULL JETS AND CHARGED JETS
	nCons.push_back(jet.constituents().size());
	
	//!~~~~~~~~~~~~~~~~~~~~~~quark vs. gluon jet studies~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	//CODE: qvg = 0 -> quark jet, qvg = 1 -> gluon jet, qvg = 2 -> neither, qvg = -1 -> no match
	//looping over each jet in the event and matching each to zero or one hard partons
     	
	//GEOMETRICALLY MATCH JET AND PARTON
	bool match = 0;
	//looping over hard partons in the event to find zero or one to match to this jet
	for (int i = 0; i < (int) hardparton.size(); ++ i) {
	  if (jet.delta_R(hardparton[i]) < _jetR) { //found a match
	    match = 1;
	    if (identity[i] == 0) { //means it's a quark!
	      /*std::cout << "matched this jet to this initial hard quark!" << std::endl;*/
	      qvg.push_back(0);
	    }
	    else if (identity[i] == 1) { //means it's a gluon!
	      /*std::cout << "matched this jet to this initial hard gluon!" << std::endl;*/
	      qvg.push_back(1);
	    }
	    else if (identity[i] == -1) { //means it's neither!
	      /*std::cout << "matched this jet to this initial hard parton that is neither quark nor gluon!" << std::endl;*/
	      qvg.push_back(2);
	    }
	    else {
	      std::cout << "WARNING: somehow the hard initiating parton does not have an identity..." << std::endl;
	    }
	    identity.erase(identity.begin() + i); //removes parton from future matching consideration
	    hardparton.erase(hardparton.begin() + i); //removes parton from future matching consideration
	  }
	  else { /*std::cout << "didn't match this jet with this hard parton" << std::endl;*/ }
	  if (match == 1) {break;} //ensures that the same jet can't be matched to multiple hard partons
	}
	if (match == 0) {/*std::cout << "didn't find a matching hard parton for this jet :(" << std::endl;*/ qvg.push_back(-1);} //no match!

	double jet_girth = 0;
	double a0 = 0;
	double a05 = 0;
	double a1 = 0;
	double a_05 = 0;
	double a_1 = 0;
	vector<double> cons_girth;
	vector<double> cons_dist;
	vector<double> cons_pt;
	vector<double> cons_mass;
	vector<double> cons_eta;
	//vector<PseudoJet> jet_cons = jet.constituents();
	foreach (const PseudoJet cons, jet.constituents()) {
	  cons_pt.push_back(cons.pt());
	  cons_mass.push_back(cons.m());
	  cons_eta.push_back(cons.eta());
	  double consR = fabs(cons.delta_R(jet));
	  cons_dist.push_back(consR);
	  double part_girth = consR * cons.pt() / (double) jet.pt();
	  cons_girth.push_back(part_girth);
	  jet_girth += part_girth;
	  //~~~angularity observables~~~//
	  a_1 += pow(consR,2+1)*cons.pt();
	  a_05 += pow(consR,2+0.5)*cons.pt();
	  a0 += pow(consR,2-0)*cons.pt();
	  a05 += pow(consR,2-0.5)*cons.pt();
	  a1 += pow(consR,2-1)*cons.pt();
	}
	
	conspT.push_back(cons_pt);
	consDist.push_back(cons_dist);
	consGirth.push_back(cons_girth);
	consM.push_back(cons_mass);
	consEta.push_back(cons_eta);
	//~~~
	//formula: tau_a = 1/pT sum (pT,i * deltaR^a) (as opposed to ^(2 - a))
	jetTau_1.push_back(a_1 / (double) jet.pt());
	jetTau_05.push_back(a_05 / (double) jet.pt());
	jetTau0.push_back(a0 / (double) jet.pt());
	jetTau05.push_back(a05 / (double) jet.pt());
	jetTau1.push_back(a1 / (double) jet.pt());
	jetGirth.push_back(jet_girth);
	jetMult.push_back(jet.constituents().size());
	jetpT.push_back(jet.pt());
	jetM.push_back(jet.m());
	jeteta.push_back(jet.eta());
	jetphi.push_back(jet.phi());

	//grooming:
	PseudoJet sd_jet = sd(jet);

	
	if(sd_jet != 0) { //! unless I'm mistaken, this condition should never evaluate to false, since grooming by default doesn't remove jets
	  
	  sdjetpT.push_back(sd_jet.pt());
	  sdjetM.push_back(sd_jet.m());
	  double z = sd_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry();
	  double r = sd_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R();
	  
	  zg.push_back(z);
	  rg.push_back(r);

	}
	else {
	  std::cout << "Somehow don't have a groomed jet for the given ungroomed jet..." << std::endl;
	  sdjetpT.push_back(-999);
	  sdjetM.push_back(-999);
	  zg.push_back(-999);
	  rg.push_back(-999);
	}
      }
      
      std::vector<fastjet::PseudoJet> HL_matches; std::vector<fastjet::PseudoJet> PL_matches;
      std::vector<int> match_indices;//alternate PL match 0, HL match 0, PL match 1, HL match 1, ...
      
      //PERFORM PARTON-TO-HADRON JET MATCHING!
      if (recoJets_PL.size() != 0) {
	HL_matches.clear(); PL_matches.clear();
	match_indices.clear();
	//match indices are a holdover from my jetmass2/src/sim code where I had to index a list of already groomed jets.
	//Here it is redundant because I groom each matched jet individually.
	match_indices = MatchJets(recoJets_bkgsub, recoJets_PL, HL_matches, PL_matches); //find matches
	if (HL_matches.size() != PL_matches.size()) {std::cerr << "Somehow we have different-sized match vectors. Exiting!" <<std::endl; exit(1);}
      }
      //FILL VECTORS WITH MATCHES!
      for (unsigned i = 0; i < HL_matches.size(); ++ i) {//HL_matches.size() == PL_matches.size() == 1/2 match_indices.size()
	//matches should be at the same index in respective vectors
	HL_pt_match.push_back(HL_matches[i].pt());
	HL_m_match.push_back(HL_matches[i].m());
	
	PL_pt_match.push_back(PL_matches[i].pt());
	PL_m_match.push_back(PL_matches[i].m());
	
	//grooming                                                                                                                                                        
	PseudoJet sd_HL = sd(HL_matches[i]);
	PseudoJet sd_PL = sd(PL_matches[i]);
	
	if (sd_HL != 0) { //should always evaluate to true
	  double z = sd_HL.structure_of<fastjet::contrib::SoftDrop>().symmetry();
	  HL_mg_match.push_back(sd_HL.m());
	  HL_zg_match.push_back(z);
	}
	if (sd_PL != 0) { //should always evaluate to true
	  double z = sd_PL.structure_of<fastjet::contrib::SoftDrop>().symmetry();
	  PL_mg_match.push_back(sd_PL.m());
	  PL_zg_match.push_back(z);
	}
      }
      
      if (PL_matches.size() != 0) {//should also be equivalent to asking if HL_matches.size() != 0
	MatchTree->Fill();
      }

      if (recoJets_PL.size() != 0) {
	PartonTree->Fill();
      }
      if (recoJets_bkgsub.size() != 0) {
	ResultTree->Fill();
      }
    }//! analyze function
    
    /// Normalise histograms etc., after the run
    void finalize() {
      
      //!Scale and write histograms here
      
      // ResultTree->Write();
      fout->Write();
      fout->Close();      
      
    }//! finalize
    
  protected:
    
    size_t _mode;
    
    //used to determine if analysis is over Pythia or Herwig HepMCs
    bool py_input;
    
    //params
    double _jetR;
    int RADIUS;
    double _pTCut;
    double _pTMax;
    double _pTconsMax;
    double _max_track_rap;
    bool printDebug;
    double z_cut;
    double beta;
    
    //event-level
    double mcweight;
    double pthat;
    double eventID;
    
    //parton-level jets
    vector<double> PLpt;
    vector<double> PLm;
    vector<double> PLptg;
    vector<double> PLmg;
    vector<double> PLeta;
    vector<double> PLphi;
    vector<double> PLzg;
    vector<double> PLrg;
    vector<double> PLncons;
    vector<double> PLnconsg;
    
    //parton-level jet constituents
    vector<vector<double> > PLconspT;
    vector<vector<double> > PLconsM;
    vector<vector<double> > PLconsEta;
    
    //matched hadron- & parton-level jets
    vector<double> PL_pt_match;
    vector<double> PL_m_match;
    vector<double> PL_mg_match;
    vector<double> PL_zg_match;
    vector<double> HL_pt_match;
    vector<double> HL_m_match;
    vector<double> HL_mg_match;
    vector<double> HL_zg_match;

    //hadron-level jet constituents
    vector<vector<double> > conspT;
    vector<vector<double> > consDist;
    vector<vector<double> > consGirth;
    vector<vector<double> > consM;
    vector<vector<double> > consEta;
    
    //hadron-level jets
    vector<int> qvg;
    vector<double> jetTau1;
    vector<double> jetTau05;
    vector<double> jetTau0;
    vector<double> jetTau_05;
    vector<double> jetTau_1;
    vector<double> jetGirth;
    vector<double> jetMult;
    vector<double> jetpT;
    vector<double> jetM;
    vector<double> sdjetpT;
    vector<double> sdjetM;
    vector<double> jeteta;
    vector<double> jetphi;
    vector<double> zg;
    vector<double> rg;
    vector<double> nCons;
    
  private:

    TFile * fout;    
    // Histo1DPtr _h_zg[2][4];
    // Scatter2DPtr _h_ratio[4];
    //! Declare histograms here, e.g.: TH1D * hJetpT[2];
    
    TTree * PartonTree;
    TTree * MatchTree;
    TTree * ResultTree;
  
    TDatabasePDG * pdg;
    
  };
  
  class STAR_MASS_10PTHAT15 : public STAR_MASS {
  public:
    STAR_MASS_10PTHAT15()
      : STAR_MASS("STAR_MASS_10PTHAT15")
    {
      _mode = 1;
    }
  };
  class STAR_MASS_15PTHAT20 : public STAR_MASS {
  public:
    STAR_MASS_15PTHAT20()
      : STAR_MASS("STAR_MASS_15PTHAT20")
    {
      _mode = 2;
    }
  };
  class STAR_MASS_20PTHAT25 : public STAR_MASS {
  public:
    STAR_MASS_20PTHAT25()
      : STAR_MASS("STAR_MASS_20PTHAT25")
    {
      _mode = 3;
    }
  };
  class STAR_MASS_25PTHAT30 : public STAR_MASS {
  public:
    STAR_MASS_25PTHAT30()
      : STAR_MASS("STAR_MASS_25PTHAT30")
    {
      _mode = 4;
    }
  };
  class STAR_MASS_30PTHAT35 : public STAR_MASS {
  public:
    STAR_MASS_30PTHAT35()
      : STAR_MASS("STAR_MASS_30PTHAT35")
    {
      _mode = 5;
    }
  };
  class STAR_MASS_35PTHAT40 : public STAR_MASS {
  public:
    STAR_MASS_35PTHAT40()
      : STAR_MASS("STAR_MASS_35PTHAT40")
    {
      _mode = 6;
    }
  };
  class STAR_MASS_40PTHAT45 : public STAR_MASS {
  public:
    STAR_MASS_40PTHAT45()
      : STAR_MASS("STAR_MASS_40PTHAT45")
    {
      _mode = 7;
    }
  };
  class STAR_MASS_45PTHAT50 : public STAR_MASS {
  public:
    STAR_MASS_45PTHAT50()
      : STAR_MASS("STAR_MASS_45PTHAT50")
    {
      _mode = 8;
    }
  };
  class STAR_MASS_50PTHAT60 : public STAR_MASS {
  public:
    STAR_MASS_50PTHAT60()
      : STAR_MASS("STAR_MASS_50PTHAT60")
    {
      _mode = 9;
    }
  };
  class STAR_MASS_60PTHAT80 : public STAR_MASS {
  public:
    STAR_MASS_60PTHAT80()
      : STAR_MASS("STAR_MASS_60PTHAT80")
    {
      _mode = 10;
    }
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(STAR_MASS);
  DECLARE_RIVET_PLUGIN(STAR_MASS_10PTHAT15);
  DECLARE_RIVET_PLUGIN(STAR_MASS_15PTHAT20);
  DECLARE_RIVET_PLUGIN(STAR_MASS_20PTHAT25);
  DECLARE_RIVET_PLUGIN(STAR_MASS_25PTHAT30);
  DECLARE_RIVET_PLUGIN(STAR_MASS_30PTHAT35);
  DECLARE_RIVET_PLUGIN(STAR_MASS_35PTHAT40);
  DECLARE_RIVET_PLUGIN(STAR_MASS_40PTHAT45);
  DECLARE_RIVET_PLUGIN(STAR_MASS_45PTHAT50);
  DECLARE_RIVET_PLUGIN(STAR_MASS_50PTHAT60);
  DECLARE_RIVET_PLUGIN(STAR_MASS_60PTHAT80);
  				       
  }
