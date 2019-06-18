//!  functions.cxx
//!  Isaac Mooney, WSU - June 2019

#include "params.hh"
#include "funcs.hh"
//!#include "ktTrackEff.hh"

typedef fastjet::contrib::SoftDrop SD;

namespace Analysis {

  //! -------------------------                                                                                                                                                                            
  //! IO/OS Manip functionality                                                                                                                                                                            
  //! -------------------------                                                                                                                                                                             
  //! Used to understand which format of input file is being used                                                                                                                                          
  //! ( .root file, .txt, .list, etc )                                                                                                                                                                     
  //! ---------------------------------------------------------------------                                                                                                                                 
  bool HasEnding (std::string const &full_string, std::string const &ending) {
    if (full_string.length() >= ending.length()) {
      return (0 == full_string.compare (full_string.length() - ending.length(), ending.length(), ending) );
    } else {
      return false;
    }
  }

  //! parse a CSV file to a set of unique entries.
  //! All comments must start on their own line, and be proceeded   
  //! by a pound sign (#)                                                                                                                                                              
  template <typename T>
  std::set<T> ParseCSV(std::string csv) {
    //! return set                                                                                                        
    std::set<T> ret;
    std::ifstream fs(csv);
    std::string line;
    //! first, split by line
    while (std::getline(fs, line)) {
      if (line.size() == 0) //!reject empty lines   
	continue;
      if (line[0] == '#') //!reject comments                                                                                                                                        
	continue;
      //! split the string by commas                                                                                                                                                   
      std::istringstream ss(line);
      while (ss) {
	std::string str_value;
	std::getline(ss, str_value, ',');
	if (CanCast<T>(str_value)) {
	  ret.insert(CastTo<T>(str_value));
	}
      }
    }
    return ret;
  }

  template<typename T>
  bool CanCast(std::string s) {
    std::istringstream iss(s);
    T dummy;
    iss >> std::skipws >> dummy;
    return iss && iss.eof();
  }

  template<typename T>
  T CastTo(std::string s) {
    std::istringstream iss(s);
    T dummy;
    iss >> std::skipws >> dummy;
    return dummy;
  }

  //  INITIATE READER with some selections
  void InitReader( TStarJetPicoReader * reader, TChain* chain, int nEvents, const std::string trig, const double vZ, const double vZDiff, const double Pt, const double Et, const double Etmin, const double DCA, const double NFit, const double NFitRatio, const double maxEtTow, const double hc, const bool mip_correction, const std::string badTows, const std::string bad_run_list) {
    // set the chain
    if (chain != nullptr) {reader->SetInputChain( chain );}
    // apply hadronic correction - subtract 100% of charged track energy from towers
    reader->SetApplyFractionHadronicCorrection( true );
    reader->SetFractionHadronicCorrection( hc ); //0.9999
    reader->SetRejectTowerElectrons( kFALSE );
    reader->SetApplyMIPCorrection( mip_correction );
    // if bad run list is specified, add to reader
    if (!bad_run_list.empty()) {
      std::set<int> bad_runs = ParseCSV<int>(bad_run_list);
      for (auto run : bad_runs) {
	reader->AddMaskedRun(run);
      }
    }
    // Event and track selection - all explained in params.hh
    // -------------------------
    TStarJetPicoEventCuts* evCuts = reader->GetEventCuts();
    evCuts->SetTriggerSelection( trig.c_str() ); //All, MB, HT, pp, ppHT, ppJP2. For pAu I think trigger strings have now been hardcoded, but haven't tested
    evCuts->SetVertexZCut ( vZ );
    evCuts->SetVertexZDiffCut( vZDiff );
    evCuts->SetRefMultCut( refMultCut );
    evCuts->SetMaxEventPtCut( Pt );
    evCuts->SetMaxEventEtCut( Et );
    evCuts->SetMinEventEtCut( Etmin );
    
    std::cout << "Using these event cuts:" << std::endl;
    std::cout << " trigger: " << evCuts->GetTriggerSelection() << std::endl;
    std::cout << " vz: " << evCuts->GetVertexZCut() << std::endl;
    std::cout << " |vpdvz - tpcvz|: " << evCuts->GetVertexZDiffCut() << std::endl;
    std::cout << " event max pT: " << evCuts->GetMaxEventPtCut() << std::endl;
    std::cout << " event max Et: " << evCuts->GetMaxEventEtCut() << std::endl;
    std::cout << " event min Et: " << evCuts->GetMinEventEtCut() << std::endl;

    // Tracks cuts - all explained in params.hh
    TStarJetPicoTrackCuts* trackCuts = reader->GetTrackCuts();
    trackCuts->SetDCACut( DCA );
    trackCuts->SetMinNFitPointsCut( NFit );
    trackCuts->SetFitOverMaxPointsCut( NFitRatio );    
    
    std::cout << "Using these track cuts:" << std::endl;
    std::cout << " dca : " << trackCuts->GetDCACut(  ) << std::endl;
    std::cout << " nfit : " <<   trackCuts->GetMinNFitPointsCut( ) << std::endl;
    std::cout << " nfitratio : " <<   trackCuts->GetFitOverMaxPointsCut( ) << std::endl;
    
    // Towers - all explained in params.hh
    TStarJetPicoTowerCuts* towerCuts = reader->GetTowerCuts();
    towerCuts->SetMaxEtCut( maxEtTow );
    towerCuts->AddBadTowers( badTows );
    std::cout << badTows << std::endl;
    
    std::cout << "Using these tower cuts:" << std::endl;
    std::cout << "  GetMaxEtCut = " << towerCuts->GetMaxEtCut() << std::endl;
    std::cout << "  Gety8PythiaCut = " << towerCuts->Gety8PythiaCut() << std::endl;
    
    // V0s: Turn off
    reader->SetProcessV0s(false);
    
    // Initialize the reader
    reader->Init( nEvents ); //runs through all events with -1
  }

  //This function simply converts TStarJetVectors from the event into PseudoJets for later clustering into jets with FastJet.
  //It also very importantly assigns a mass to each particle after the conversion.
  //The assigned mass is dependent on whether the function call is for particle-level (e.g. Pythia) [where we know the rest masses] or detector-level (e.g. Geant, data) [where we don't know them].
  void GatherParticles ( TStarJetVectorContainer<TStarJetVector> * container, TStarJetVector *sv, std::vector<fastjet::PseudoJet> & Particles, const bool full, const bool py, TDatabasePDG *pdg){
    for ( int i = 0; i < container->GetEntries() ; ++i ) {
      sv = container->Get(i);
      fastjet::PseudoJet current = fastjet::PseudoJet( *sv );

      if (sv->GetCharge() != 0 && !py) { // charged at detector-level -> charged pion mass
        current.reset_PtYPhiM(sqrt(current.perp2()),current.rap(),current.phi(), chPionMass); //assigning pion mass to charged particles
      }
      else if (sv->GetCharge() == 0 && !py) { // neutral at detector-level -> zero mass
        current.reset_PtYPhiM(sqrt(current.perp2()),current.rap(),current.phi(), 0); //neutral particles massless!
      }
      else if (py) { // at particle-level -> PDG mass
        current.reset_PtYPhiM(sqrt(current.perp2()), current.rap(), current.phi(), pdg->GetParticle(sv->mc_pdg_pid())->Mass());
      }
      //DEBUG:
      //std::cout << "Particle-level? " << py << std::endl << "ASSIGNING MASS " << current.m() << " TO CONSTITUENT " << i << " WITH CHARGE " << sv->GetCharge() << " AND MASS " << sv->GetPicoMass() << std::endl;
      //std::cout << "THE PARTICLE, WITH PID " << sv->mc_pdg_pid() << ", HAS MASS " << pdg->GetParticle(sv->mc_pdg_pid())->Mass() << "[should be the same as the last number above]" << std::endl;               

      if ((sv->GetCharge() == 0) && (full == 0)) { continue; } //!if we don't want full jets, skip neutrals                                                     
      current.set_user_index( sv->GetCharge() );
      
      //DEBUG:
      //std::cout << "ending with charge: " << current.user_index() << " for particle " << i << std::endl;
      
      Particles.push_back(current);
    }
    return;
  }

  //this function takes in a jet sample, and keeps whichever ones pass a neutral energy fraction selection (i.e. have >= x% of their energy in tracks) 
  //!the name is a bit of a misnomer, we're actually using pT, not energy
  void ApplyNEFSelection(const std::vector<fastjet::PseudoJet> init, std::vector<fastjet::PseudoJet> & result) {
    //!Implementing a neutral energy fraction cut on inclusive jets 
    for (int i = 0; i < init.size(); ++ i) { //looping over the jet
      double towersum = 0; double ptsum = 0; //start with zero energy, neutral or otherwise
      for (int j = 0; j < init[i].constituents().size(); ++ j) { //loop over the jet constituents
        if (init[i].constituents()[j].user_index() == 0) { //if it's neutral, add its pT to the towers
          towersum += init[i].constituents()[j].pt();
        }
        ptsum += init[i].constituents()[j].pt(); //whether or not it's neutral, add its pT to the total
      }//constituent loop
      if (towersum / (double) ptsum < NEF_max) { //if the fraction in the towers is less than the threshold, accept it
        result.push_back(init[i]);
      }//if
    }//jet loop
    return;
  }
  
  //~~~~~~~FILL HISTS [not used currently, but might be useful in the future so saving it here]~~~~~~~//
  
  void FillHistsHelper(Collection<std::string, TH1D> & c1D, Collection<std::string, TH2D> &c2D, Collection<std::string, TH3D> & c3D, const std::string flag, const fastjet::PseudoJet jet, const double weight) {
    c1D.fill(("m_" + flag + "_" + "jet").c_str(), jet.m(), weight);
    c2D.fill(("m_v_pt_" + flag + "_" + "jet").c_str(), jet.m(), jet.pt(), weight);
    c2D.fill(("m_v_pt_rebin_" + flag + "_" + "jet").c_str(), jet.m(), jet.pt(), weight);
    c3D.fill(("PtEtaPhi_" + flag + "_" + "jet").c_str(), jet.pt(), jet.eta(), jet.phi(), weight);
    for (int cons = 0; cons < jet.constituents().size(); ++ cons) {
      if (jet.constituents()[cons].pt() < partMinPt) {continue;} //ignores contributions from ghosts                       
      c3D.fill(("PtEtaPhi_" + flag + "_" + "cons").c_str(), jet.constituents()[cons].pt(), jet.constituents()[cons].eta(), jet.constituents()[cons].phi(), weight);
    }
    return;
  }
  
  void FillHists(Collection<std::string, TH1D> & c1D, Collection<std::string, TH2D> & c2D, Collection<std::string, TH3D> & c3D, const std::vector<fastjet::PseudoJet> jets, const double weight) {
    //leading
    FillHistsHelper(c1D, c2D, c3D, "lead", jets[0], weight);
    //subleading
    if (jets.size() > 1) {
      FillHistsHelper(c1D, c2D, c3D, "sublead", jets[1], weight);
    }
    //inclusive
    for (int i = 0; i < jets.size(); ++ i) {
      FillHistsHelper(c1D, c2D, c3D, "incl", jets[i], weight);
    }
    /*
    //trigger & recoil
    std::vector<fastjet::PseudoJet> candidates;
    bool which_one = GetTriggerJet(candidates, jets);
    if (candidates.size() == 0 || jets.size() < 2) { // means there isn't a trigger or there isn't a recoil
      return;
    }
    if (candidates.size() == 1 && jets.size() > 1) { //potential trigger
      if (fabs(fabs(jets[which_one].delta_phi_to(jets[(which_one + 1) % 2])) - Pi) < R) { //found a recoil
	FillHistsHelper(c1D, c2D, c3D, flag1, flag2, "trig", jets[which_one], weight); //filling hists for trigger
	FillHistsHelper(c1D, c2D, c3D, flag1, flag2, "rec", jets[(which_one + 1) % 2], weight); //filling hists for recoil
	return;
      }
    }
    if (candidates.size() == 2) {
      if (fabs(fabs(candidates[0].delta_phi_to(candidates[1])) - Pi) < R) { //trigger & recoil found!
	FillHistsHelper(c1D, c2D, c3D, flag1, flag2, "trig", candidates[0], weight); //filling hists for trigger
	FillHistsHelper(c1D, c2D, c3D, flag1, flag2, "rec", candidates[1], weight); //filling hists for recoil
	return;
      }
    }
    */
    return;
  }

  void FillSDHistsHelper(Collection<std::string, TH1D> & c1D, Collection<std::string, TH2D> &c2D, Collection<std::string, TH3D> & c3D, const std::string flag, const fastjet::PseudoJet jet, const double weight) {
    c1D.fill(("m_" + flag + "_" + "sd").c_str(), jet.m(), weight);
    c1D.fill(("zg_" + flag + "_" + "sd").c_str(), jet.structure_of<fastjet::contrib::SoftDrop>().symmetry(), weight);
    c1D.fill(("thetag_" + flag + "_" + "sd").c_str(), jet.structure_of<fastjet::contrib::SoftDrop>().delta_R(), weight);
    c2D.fill(("m_v_pt_" + flag + "_" + "sd").c_str(), jet.m(), jet.pt(), weight);
    c3D.fill(("PtEtaPhi_" + flag + "_" + "sd").c_str(), jet.pt(), jet.eta(), jet.phi(), weight);
    return;
  }
   

  void FillSDHists(Collection<std::string, TH1D> & c1D, Collection<std::string, TH2D> & c2D, Collection<std::string, TH3D> & c3D, const std::vector<fastjet::PseudoJet> jets, const double weight) {
    //leading
    FillSDHistsHelper(c1D, c2D, c3D, "lead", jets[0], weight);
    //subleading
    if (jets.size() > 1) {
      FillSDHistsHelper(c1D, c2D, c3D, "sublead", jets[1], weight);
    }
    //inclusive
    for (int i = 0; i < jets.size(); ++ i) {
      FillSDHistsHelper(c1D, c2D, c3D, "incl", jets[i], weight);
    }

  }

}
