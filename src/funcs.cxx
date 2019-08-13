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

  //gets the cross section for each pT bin for run 6   
  double LookupRun6Xsec(TString currentfile ) {

    static const Double_t Xsec[12] = {
      1.0,      //!Placeholder for 2-3                                                                                                                      
      1.30E+09, //!3-4                                                                                                                                      
      3.15E+08, //!4-5                                                                                                                                      
      1.37E+08, //!5-7                                                                                                                                      
      2.30E+07, //!7-9                                                                                                                                      
      5.53E+06, //!9-11                                                                                                                                     
      2.22E+06, //!11-15                                                                                                                                    
      3.90E+05, //!15-25                                                                                                                                    
      1.02E+04, //!25-35                                                                                                                                    
      5.01E+02, //!35-45                                                                                                                                    
      2.86E+01, //!45-55                                                                                                                                    
      1.46E+00  //!55-65                                                                                                                                    
    };

    static const Double_t Nmc[12] = {
      1,      //!2-3                                                                                                                                        
      672518, //!3-4                                                                                                                                        
      672447, //!4-5                                                                                                                                        
      393498, //!5-7                                                                                                                                        
      417659, //!7-9                                                                                                                                        
      412652, //!9-11                                                                                                                                       
      419030, //!11-15                                                                                                                                      
      396744, //!15-25                                                                                                                                      
      399919, //!25-35                                                                                                                                      
      119995, //!35-45                                                                                                                                      
      117999, //!45-55                                                                                                                                      
      119999  //!55-65                                                                                                                                      
    };

    Double_t w[12];
    for ( int i=0; i<12 ; ++i ){
      w[i] = Xsec[i] / Nmc[i];
    }

    if ( currentfile.Contains("3_4") ) return w[1];
    if ( currentfile.Contains("4_5") ) return w[2];
    if ( currentfile.Contains("5_7") ) return w[3];
    if ( currentfile.Contains("7_9") ) return w[4];
    if ( currentfile.Contains("9_11") ) return w[5];
    if ( currentfile.Contains("11_15") ) return w[6];
    if ( currentfile.Contains("15_25") ) return w[7];
    if ( currentfile.Contains("25_35") ) return w[8];
    if ( currentfile.Contains("35_45") ) return w[9];
    if ( currentfile.Contains("45_55") ) return w[10];
    if ( currentfile.Contains("55_65") ) return w[11];
    return 1;
  }

  //gets the cross section for each pT bin for run 12
  double LookupRun12Xsec( TString filename ){
    const int NUMBEROFPT = 11;
    //! const char *PTBINS[NUMBEROFPT]={"2_3","3_4","4_5","5_7","7_9","9_11","11_15","15_20","20_25","25_35","35_-1"};                                      
    const static float XSEC[NUMBEROFPT] = {9.00581646, 1.461908221, 0.3544350863, 0.1513760388, 0.02488645725, 0.005845846143, 0.002304880181, 0.000342661835, 4.562988397e-05, 9.738041626e-06, 5.019978175e-07};
    const static float NUMBEROFEVENT[NUMBEROFPT] = {2100295, 600300, 600300, 300289, 300289, 300289, 160295, 100302, 80293, 76303, 23307};
    const static std::vector<std::string> vptbins={"pp12Pico_pt2_3","pp12Pico_pt3_4","pp12Pico_pt4_5","pp12Pico_pt5_7","pp12Pico_pt7_9","pp12Pico_pt9_11","pp12Pico_pt11_15","pp12Pico_pt15_20","pp12Pico_pt20_25","pp12Pico_pt25_35","35_-1"};
    
    for ( int i=0; i<vptbins.size(); ++i ){
      if ( filename.Contains(vptbins.at(i).data())) return XSEC[i] / NUMBEROFEVENT[i];
    }

    throw std::runtime_error("Not a valid filename");
    return -1;
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
  
  //discards events on the grounds of them having jets of pT > double the high end of the pT-hat bin from which they came. Both the Py & Py+Ge event will be thrown out.
  //NOTE: I don't think the Picos save the pT-hat bin edges from the production, so I have to use the filenames instead
  bool DiscardEvent(const TString Filename, const std::vector<fastjet::PseudoJet> p_Jets, const std::vector<fastjet::PseudoJet> g_Jets) {
    bool bad_event = 0;
    //have to do some manipulations to get the upper edge of the pT-hard bin from the file name, e.g. isolating "25" from 20_25
    //unfortunately this only works for the naming convention of the Y12 Pythia6 and Pythia6+Geant files on disk. Will need to be rewritten slightly if
    //used on different files.
    std::string tail = ((std::string) Filename).substr(((std::string) Filename).size() - 10);
    std::string upstring = tail.substr(0,2);
    std::string upstring_copy = upstring;
    //since we search from the end of the string backward, there is a difference between files with two 2-digit pt-hat edges, with one, and with zero.
    if (upstring.find("_") != std::string::npos || upstring.find("-") != std::string::npos) {
      if (upstring.substr(1,1) != "_") {
	upstring = upstring.substr(1,1);
      }
      else {
	upstring = upstring.substr(0,1);
      }
    }
    int upbin = std::stoi(upstring);
    if (p_Jets.size() != 0) {
      if ((p_Jets[0].pt() > 2*upbin) && upstring_copy != "-1") {
	std::cout << "Removing event from file " << Filename << /*<< " with weight " << mc_weight*/ " due to bad jet with pt, eta, phi, and m: " << p_Jets[0].pt() << " " << p_Jets[0].eta() << " " << p_Jets[0].phi() << " " << p_Jets[0].m() << std::endl;
	bad_event = 1;
      }
    }
    if (g_Jets.size() != 0) {
      if ((g_Jets[0].pt() > 2*upbin) && upstring_copy != "-1") {
	std::cout << "Removing event from file " << Filename << /*<< " with weight " << mc_weight*/ " due to bad jet with pt, eta, phi, and m: " << g_Jets[0].pt() << " " << g_Jets[0].eta() << " " << g_Jets[0].phi() << " " << g_Jets[0].m() << std::endl;
	bad_event = 1;
      }
    }
    
    return bad_event;
  }
  
  //This function takes a vector of jets to be geometrically matched to another vector of "candidate" matches. Once matching occurs, a vector of indices is returned allowing one to index the original two vectors to fill responses etc. with the matched jet pairs. Redundantly (for debugging, etc.) the vectors of matches themselves are also updated for later use.
  //Note: We need to be able to remove jets from the "candidates" vector after they've been matched, so we make a copy in the function. Also make a copy candidates vector for each iteration on toMatch since this vector has selections applied to it
  //Note: In finding which jets were the matches, we know the toMatch jet match will be the 'i'th jet since we are iterating. The candidate_copy jet should be the highest pT match, so the first one in the candidate_copy list. Geometrically match the candidate_copy jet to the nearest candidate jet, since jets have been removed so they don't index to the same jet anymore
  std::vector<int> MatchJets(const std::vector<fastjet::PseudoJet> candidates_safe, const std::vector<fastjet::PseudoJet> toMatch, std::vector<fastjet::PseudoJet> & c_matches, std::vector<fastjet::PseudoJet> & t_matches) {
    std::vector<int> match_indices;
    if (candidates_safe.size() == 0 || toMatch.size() == 0) {
      return match_indices; //later, match_indices being empty will tell us there were no matches
    }
    //define candidates outside the loop so list continually dwindles as we remove matched candidates
    std::vector<fastjet::PseudoJet> candidates = candidates_safe;
    for (int i = 0; i < toMatch.size(); ++ i) { //for each jet in toMatch, we try to find a match from candidates_copy
      //defined inside the loop so that for each toMatch jet there's a new set of candidates
      std::vector<fastjet::PseudoJet> candidates_copy = candidates;
      fastjet::Selector selectMatchedJets = fastjet::SelectorCircle( R );
      selectMatchedJets.set_reference( toMatch[i] );
      //note: matchedToJet and candidates_copy are equivalent, assuming candidates_safe was already sorted by pT
      std::vector<fastjet::PseudoJet> matchedToJet = sorted_by_pt( selectMatchedJets( candidates_copy ));
      if (matchedToJet.size() == 0) { continue; } //means no match to this jet. Remove none from candidates. Continuing on to the next one.
      else { //found at least one match. Need to remove the highest pT one from candidates and add the respective jets to the match vectors.
	match_indices.push_back(i); //push back the toMatch match position
	t_matches.push_back(toMatch[i]);
	c_matches.push_back(matchedToJet[0]); //highest pT match
	for (int j = 0; j < candidates.size(); ++ j) { //finding which one to delete from candidates before next toMatch iteration.
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
  
  //!this function is similar to the "MatchJets" function, but it instead takes a list of jets which have already been matched ("candidates_safe") and finds the jets to which they correspond in the list of unmatched jets ("toMatch") by (exact) geometrical matching. The remaining jets with no corresponding already-matched jets are either misses or fakes depending on the call.
  std::vector<int> FakesandMisses(const std::vector<fastjet::PseudoJet> candidates_safe, const std::vector<fastjet::PseudoJet> toMatch, std::vector<fastjet::PseudoJet> & unmatched) {
    std::vector<int> miss_fake_index;
    std::vector<fastjet::PseudoJet> candidates = candidates_safe;
    for (int i = 0; i < toMatch.size(); ++ i) {
      std::vector<fastjet::PseudoJet> candidates_copy = candidates;
      fastjet::Selector selectMatchedJets = fastjet::SelectorCircle( 0.0001 ); //a "match" is now if we found the same exact jet.
      selectMatchedJets.set_reference( toMatch[i] );
      std::vector<fastjet::PseudoJet> matchedToJet = sorted_by_pt( selectMatchedJets( candidates_copy ));
      if (matchedToJet.size() == 0) { //means no match to this jet. Remove none from candidates. Add it to unmatched & continue to next one.
	miss_fake_index.push_back(i); //the ith jet in toMatch is a miss or a fake
	unmatched.push_back(toMatch[i]);
	continue;
      }
      else { //found at least one match. Need to remove the highest pT one from candidates. [Should just be 1 match]
	for (int j = 0; j < candidates.size(); ++ j) { //finding which one to delete from candidates before next toMatch iteration.
	  if (matchedToJet[0].delta_R(candidates[j]) < 0.0001) { //is probably the same jet
	    candidates.erase(candidates.begin() + j); //removing the jet from the overall list of candidates so it can't be considered next time
	    break; //should exit only the candidates loop.
	  }
	}
      }
    }
    return miss_fake_index;
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
