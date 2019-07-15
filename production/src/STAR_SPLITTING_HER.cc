//! -*- C++ -*-OA

//! Raghav Kunnawalkam Elayavalli
//! Wayne State University
//! Friday Jan 19th, 2018 (modified from earlier code to do zg at star)
//!
//! Analysis to calculate the splitting function at STAR for leading/subleading jets. 
//! The splitting function estimation is very similar to the CMS analysis except this will be at much lower pT 
//! 
//! List of plots:  
//! 1) zg for one centrality bins and leading/sub-leading jet pT bins
//! 2) rg for one centrlaity bins and leading/sub-leading jet pT bins 
//! 2) Aj as in asymmetry of leading/subleading jets 
//! 3) ratios for PbPb and pp
//!
//! Need to add plots that compare zg with the z from 2 subjets and 3 subjets 
//!

//! herwig pthat bin - 3-4, 4-5, 5-7, 7-9, 9-11, 11-15, 15-25, 25-35, 35-45, 45-55, 55-65  

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

double m_chpion = 0.13957018; //GeV/c^2

namespace Rivet {

  class STAR_SPLITTING_HER : public Analysis {
  public:
    
    /// Constructor
    STAR_SPLITTING_HER(string name = "STAR_SPLITTING_HER")
      : Analysis(name)
    {
      setNeedsCrossSection(true);
      _mode=0;
    }
    
    /*
    //! first setup only for pp jets without any centrlaity bins - so no BKG SUB

    //! candidates inside jets. need to draw a grid inside each jet to fill with candidates 
    struct candidate{
      //! these are the boundaries 
      double etaMin;
      double phiMin;
      double etaMax;
      double phiMax;
      //! Pseudojets will have the objects inside 
      PseudoJets objects;
      vector <int> jetID;
      //! Four momeentum to get useful informations like eta, phi, mass, pT etc... 
      FourMomentum candMom;
      FourMomentum bkgMom;
      double sumnegpT;
      double sumpT;
    };

    std::pair<int,int> FindCandidateBin(double eta, double phi){
      int etaloc = -9;
      int philoc = -9;      
      for(unsigned x = 0; x<_etabins.size()-1; ++x){
	for(unsigned y = 0; y<_phibins.size()-1; ++y){
	  if(eta >= _etabins[x] && eta < _etabins[x+1] &&
	     phi >= _phibins[y] && phi < _phibins[y+1]){
	    etaloc = x;
	    philoc = y;
	  }	  
	}
      }
      if(etaloc == -9 || philoc == -9)
	return std::pair<int,int>(-1,-1);
      else
	return std::pair<int,int>(etaloc,philoc);
    }
    */

    //! http://journals.aps.org/prd/pdf/10.1103/PhysRevD.91.111501
    //! Anti kT jet is reclustered using CA and min(pT1, pT2)/(pT1+pT2) > z_cut * (R12/Rj)^\beta
    //! for sub-jets that satisfy this condition, z_g = min(pT1,pT2)/(pT1+pT2)
    std::pair<double,double>  SoftDrop(const fastjet::PseudoJet &j, double jetR, double z_cut, double beta)
    {
      //! give the soft drop groomer a short name
      //! Use a symmetry cut z > z_cut R^beta
      //! By default, there is no mass-drop requirement
      fastjet::contrib::RecursiveSymmetryCutBase::SymmetryMeasure  symmetry_measure = fastjet::contrib::RecursiveSymmetryCutBase::scalar_z;
      fastjet::contrib::SoftDrop sd(beta, z_cut, symmetry_measure, jetR);
      PseudoJet sd_jet = sd(j);    
      if(sd_jet != 0){
	double z = sd_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry();
	double r = sd_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R();
	return std::pair<double,double>(z,r);
      }else 
	return std::pair<double,double>(-1,-1);
    }
       
    /// Book histograms and initialise projections before the run
    void init() {

      
      //for later use looking up PDG masses using particle PID
      pdg = new TDatabasePDG();
    

      //! jet pT bins, lead and sublead 
      // _ptleadedges = {10.0, 20.0, 30.0, 40.0, 60.0};
      // _nptleadbins = _ptleadedges.size()-1;
      // _ptsubleadedges = {5.0, 10.0, 20.0, 30.0, 40.0};
      // _nptsubleadbins = _ptsubleadedges.size()-1;

      //! Total ptbins 
      _ptedges = {5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 60.0};
      _nptbins = _ptedges.size()-1;
      _pTCut = 0.2;
      _pTMax = 1000.0;
      _pTconsMax = 30.0;
      _max_track_rap = 1.0;
      
      //! jet eta max = 0.6
      RADIUS = 4;
      _jetR = (float)RADIUS/10;

      // _jetparams += 0.2, 0.3, 0.4, 0.5, 0.6;
      // _njetbins= _jetparams.size();

      //! debug flag 
      printDebug = false;

      // zg bins
      // double _zgbins[] = {0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55};
      // _zgcenters += 0.0625, 0.0875, 0.1125, 0.1375, 0.1625, 0.1875, 0.2125, 0.2375, 0.2625, 0.2875, 0.3125, 0.3375, 0.3625, 0.3875, 0.4125, 0.4375, 0.4625, 0.4875, 0.5125, 0.5375,;
      
      // fine centrality bins 
      _centedges = {-10., 0.0, 0.2};
      _ncentbins = _centedges.size()-1;
      _Nj.resize(_ncentbins);
      _Nj_HC.resize(_ncentbins);
      _Nj_MC.resize(_ncentbins);
      // _Nlj.resize(_ncentbins);//! leading jet # Jets 
      // _Nslj.resize(_ncentbins);//! subleading jets 
      // _NevHardcore.resize(_ncentbins);
      // _NevMatched.resize(_ncentbins);
      for (size_t i = 0; i < _ncentbins; ++i) {
        _Nj[i].resize(_nptbins,0.0);
	_Nj_HC[i].resize(_nptbins, 0.0);
        _Nj_MC[i].resize(_nptbins, 0.0);
        // _Nlj[i].resize(_nptleadbins,0.0);
        // _Nslj[i].resize(_nptsubleadbins,0.0);
      }

      /*
      _delRMin=0.05;
      _etaMax=1.0;
      _phiMax    = M_PI;
      _Nbounds_eta = (int)2*_etaMax/(_delRMin);
      _Nbounds_phi = (int)2*_phiMax/(_delRMin);
      //! fill the grid: 
      for(int i = 0; i<=_Nbounds_eta; ++i){
	double etax = -1*_etaMax + i*(double)2*_etaMax/_Nbounds_eta;
	_etabins.push_back(etax);
      }
      for(int i = 0; i<=_Nbounds_phi; ++i){
	double phix = -1*_phiMax + i*(double)2*_phiMax/_Nbounds_phi;
	_phibins.push_back(phix);	
      }
      */
      
      //! parameters for splitting function
      z_cut=0.1;
      beta = 0;      
      //! final state and jet projection
      //! generic final state
      FinalState fs(-1.0, 1.0, 0.*GeV);
      addProjection(fs, "FS");
      FastJets fj(fs, FastJets::ANTIKT, _jetR);
      addProjection(fj, "Jets");

      FinalPartons fp(Cuts::abseta < 1.0);
      addProjection(fp, "FP");
      FastJets pj(fp, FastJets::ANTIKT, _jetR);
      addProjection(pj, "PartonJets");


      
      if(_mode == 0){
	fout = new TFile(Form(/*"Results/py6_decayed_jewel_pthatbin580_R0%d.root"*/"Results/star_subjetvars_herwig7_5pthatbin10_R0%d_decays_off_noHad_test.root", RADIUS),"RECREATE");
	fout->cd();
      }    
      else if(_mode == 1){
	fout = new TFile(Form("Results/star_subjetvars_herwig7_10pthatbin15_R0%d_decays_off_noHad_test.root", RADIUS),"RECREATE");
	fout->cd();
      }else if(_mode == 2){
	fout = new TFile(Form("Results/star_subjetvars_herwig7_15pthatbin20_R0%d_decays_off_noHad_test.root", RADIUS),"RECREATE");
	fout->cd();
      }else if(_mode == 3){
	fout = new TFile(Form("Results/star_subjetvars_herwig7_20pthatbin25_R0%d_decays_off_noHad_test.root", RADIUS),"RECREATE");
	fout->cd();
      }else if(_mode == 4){
	fout = new TFile(Form("Results/star_subjetvars_herwig7_25pthatbin30_R0%d_decays_off_noHad_test.root", RADIUS),"RECREATE");
	fout->cd();
      }else if(_mode == 5){
	fout = new TFile(Form("Results/star_subjetvars_herwig7_30pthatbin35_R0%d_decays_off_noHad_test.root", RADIUS),"RECREATE");
	fout->cd();
      }else if(_mode == 6){
	fout = new TFile(Form("Results/star_subjetvars_herwig7_35pthatbin40_R0%d_decays_off_noHad_test.root", RADIUS),"RECREATE");
	fout->cd();
      }else if(_mode == 7){
	fout = new TFile(Form("Results/star_subjetvars_herwig7_40pthatbin45_R0%d_decays_off_noHad_test.root", RADIUS),"RECREATE");
	fout->cd();
      }else if(_mode == 8){
	fout = new TFile(Form("Results/star_subjetvars_herwig7_45pthatbin50_R0%d_decays_off_noHad_test.root", RADIUS),"RECREATE");
	fout->cd();
      }else if(_mode == 9){
	fout = new TFile(Form("Results/star_subjetvars_herwig7_50pthatbin60_R0%d_decays_off_noHad_test.root", RADIUS),"RECREATE");
	fout->cd();
      }else if(_mode == 10){
	fout = new TFile(Form("Results/star_subjetvars_herwig7_60pthatbin80_R0%d_decays_off_noHad_test.root", RADIUS),"RECREATE");
	fout->cd();
      }
     
      
      //! histograms 
      for(size_t i = 0; i<_ncentbins; ++i){

	hJetpT[i] = new TH1F(Form("hJetpT_centbin%d", i), "", 100, 0, 100);
	hJetpT[i]->Sumw2();
	hJetEta[i] = new TH1F(Form("hJetEta_centbin%d", i), "", 50, -1, 1);
	hJetEta[i]->Sumw2();
	hJetPhi[i] = new TH1F(Form("hJetPhi_centbin%d", i), "", 50, M_PI, PI);
	hJetPhi[i]->Sumw2();

	hHCJetpT[i] = new TH1F(Form("hHCJetpT_centbin%d", i), "", 100, 0, 100);
	hHCJetpT[i]->Sumw2();
	hHCJetEta[i] = new TH1F(Form("hHCJetEta_centbin%d", i), "", 50, -1, 1);
	hHCJetEta[i]->Sumw2();
	hHCJetPhi[i] = new TH1F(Form("hHCJetPhi_centbin%d", i), "", 50, M_PI, PI);
	hHCJetPhi[i]->Sumw2();
	
	hMCJetpT[i] = new TH1F(Form("hMCJetpT_centbin%d", i), "", 100, 0, 100);
	hMCJetpT[i]->Sumw2();
	hMCJetEta[i] = new TH1F(Form("hMCJetEta_centbin%d", i), "", 50, -1, 1);
	hMCJetEta[i]->Sumw2();
	hMCJetPhi[i] = new TH1F(Form("hMCJetPhi_centbin%d", i), "", 50, M_PI, PI);
	hMCJetPhi[i]->Sumw2();

	//! centrality bin 
	for (size_t j = 0; j < _nptbins; ++j){
	  //! ptbins
	  // _h_zg[i][j] = bookHisto1D(1, i+1, j+1);
	  hzg[i][j] = new TH1F(Form("zg_centbin%d_ptbin%d", i, j), "", 20, 0.05, 0.55);
	  hzg[i][j]->Sumw2();
	  hRg[i][j] = new TH1F(Form("Rg_centbin%d_ptbin%d", i, j), "", 30, 0.00, 0.6);
	  hRg[i][j]->Sumw2();
	  hMj[i][j] = new TH1F(Form("Mj_centbin%d_ptbin%d", i, j), "", 20, 0, 10);
	  hMj[i][j]->Sumw2();
	  if(i==1){
	    // _h_ratio[j] = bookScatter2D(2, i, j+1);
	    hRatio[j] = new TH1F(Form("zg_ratio_ptbin%d", j), "", 20, 0.05, 0.55);
	    hRatio[j]->Sumw2();
	    hRatio_Rg[j] = new TH1F(Form("Rg_ratio_ptbin%d", j), "", 30, 0.0, 0.6);
	    hRatio_Rg[j]->Sumw2();
	    hRatio_Mj[j] = new TH1F(Form("Mj_ratio_ptbin%d", j), "", 20, 0, 10);
	    hRatio_Mj[j]->Sumw2();
	  }
	  hHCzg[i][j] = new TH1F(Form("HCzg_centbin%d_ptbin%d", i, j), "", 20, 0.05, 0.55);
	  hHCzg[i][j]->Sumw2();
	  hHCRg[i][j] = new TH1F(Form("HCRg_centbin%d_ptbin%d", i, j), "", 30, 0.00, 0.6);
	  hHCRg[i][j]->Sumw2();
	  hHCMj[i][j] = new TH1F(Form("HCMj_centbin%d_ptbin%d", i, j), "", 20, 0, 10);
	  hHCMj[i][j]->Sumw2();
	  if(i==1){
	    // _h_ratio[j] = bookScatter2D(2, i, j+1);
	    hHCRatio[j] = new TH1F(Form("HCzg_ratio_ptbin%d", j), "", 20, 0.05, 0.55);
	    hHCRatio[j]->Sumw2();
	    hHCRatio_Rg[j] = new TH1F(Form("HCRg_ratio_ptbin%d", j), "", 30, 0.0, 0.6);
	    hHCRatio_Rg[j]->Sumw2();
	    hHCRatio_Mj[j] = new TH1F(Form("HCMj_ratio_ptbin%d", j), "", 20, 0, 10);
	    hHCRatio_Mj[j]->Sumw2();
	  }
	  hMCzg[i][j] = new TH1F(Form("MCzg_centbin%d_ptbin%d", i, j), "", 20, 0.05, 0.55);
	  hMCzg[i][j]->Sumw2();
	  hMCRg[i][j] = new TH1F(Form("MCRg_centbin%d_ptbin%d", i, j), "", 30, 0.00, 0.6);
	  hMCRg[i][j]->Sumw2();
	  hMCMj[i][j] = new TH1F(Form("MCMj_centbin%d_ptbin%d", i, j), "", 20, 0, 10);
	  hMCMj[i][j]->Sumw2();
	  if(i==1){
	    // _h_ratio[j] = bookScatter2D(2, i, j+1);
	    hMCRatio[j] = new TH1F(Form("MCzg_ratio_ptbin%d", j), "", 20, 0.05, 0.55);
	    hMCRatio[j]->Sumw2();
	    hMCRatio_Rg[j] = new TH1F(Form("MCRg_ratio_ptbin%d", j), "", 30, 0.0, 0.6);
	    hMCRatio_Rg[j]->Sumw2();
	    hMCRatio_Mj[j] = new TH1F(Form("MCMj_ratio_ptbin%d", j), "", 20, 0, 10);
	    hMCRatio_Mj[j]->Sumw2();
	  }
	}
	
	// for (size_t j = 0; j < _nptleadbins; ++j){
	//   // _h_zg_lead[i][j] = bookHisto1D(3, i+1, j+1);	  
	//   hzg_lead[i][j] = new TH1F(Form("zg_lead_centbin%d_ptbin%d", i, j), "", 20, 0.05, 0.55);
	//   hzg_lead[i][j]->Sumw2();
	//   hRg_lead[i][j] = new TH1F(Form("Rg_lead_centbin%d_ptbin%d", i, j), "", 30, 0., 0.6);
	//   hRg_lead[i][j]->Sumw2();
	//   hMj_lead[i][j] = new TH1F(Form("Mj_lead_centbin%d_ptbin%d", i, j), "", 20, 0, 10);
	//   hMj_lead[i][j]->Sumw2();
	//   if(i==1){
	//     // _h_ratio_lead[j] = bookScatter2D(4, i, j+1);
	//     hRatio_lead[j] = new TH1F(Form("zg_ratio_lead_ptbin%d", j), "", 20, 0.05, 0.55);
	//     hRatio_lead_Rg[j] = new TH1F(Form("Rg_ratio_lead_ptbin%d", j), "", 30, 0., 0.6);
	//     hRatio_lead_Mj[j] = new TH1F(Form("Mj_ratio_lead_ptbin%d", j), "", 20, 0, 10);
	//   }
	// }
	// for (size_t j = 0; j < _nptsubleadbins; ++j){
	//   // _h_zg_sublead[i][j] = bookHisto1D(5, i+1, j+1);	  
	//   hzg_sublead[i][j] = new TH1F(Form("zg_sublead_centbin%d_ptbin%d", i, j), "", 20, 0.05, 0.55);
	//   hzg_sublead[i][j]->Sumw2();
	//   hRg_sublead[i][j] = new TH1F(Form("Rg_sublead_centbin%d_ptbin%d", i, j), "", 30, 0., 0.6);
	//   hRg_sublead[i][j]->Sumw2();
	//   hMj_sublead[i][j] = new TH1F(Form("Mj_sublead_centbin%d_ptbin%d", i, j), "", 20, 0, 10);
	//   hMj_sublead[i][j]->Sumw2();
	//   if(i==1){
	//     // _h_ratio_sublead[j] = bookScatter2D(6, i, j+1);
	//     hRatio_sublead[j] = new TH1F(Form("zg_ratio_sublead_ptbin%d", j), "", 20, 0.05, 0.55);
	//     hRatio_sublead_Rg[j] = new TH1F(Form("Rg_ratio_sublead_ptbin%d", j), "", 30, 0., 0.6);
	//     hRatio_sublead_Mj[j] = new TH1F(Form("Mj_ratio_sublead_ptbin%d", j), "", 20, 0, 10);
	//   }
	// }
	// //!AJ plots 
	// // _h_aj_hardcore[i] = bookHisto1D(8, i+1, 1);
	// hAJ_hardcore[i] = new TH1F(Form("hAJ_hardcore_centbin%d", i), "", 30, 0, 0.9);
	// // _h_aj_matched[i] = bookHisto1D(9, i+1, 1);
	// hAJ_matched[i] = new TH1F(Form("hAJ_matched_centbin%d", i), "", 30, 0, 0.9);
      }
      
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
      PartonTree->Branch("pthat",&pthat,"pthat/d");
      PartonTree->Branch("mcweight",&mcweight,"mcweight/d");

      m0PartonTree=new TTree("m0PartonTree","m0Parton Jets");
      //! clear the vectors for each event!
      m0PLpt.clear();
      m0PLm.clear();
      m0PLptg.clear();
      m0PLmg.clear();
      m0PLeta.clear();
      m0PLphi.clear();
      m0PLzg.clear();
      m0PLrg.clear();
      m0PLncons.clear();
      m0PLnconsg.clear();
      
      m0PartonTree->Branch("m0PLpt",&m0PLpt);
      m0PartonTree->Branch("m0PLm",&m0PLm);
      m0PartonTree->Branch("m0PLptg",&m0PLptg);
      m0PartonTree->Branch("m0PLmg",&m0PLmg);
      m0PartonTree->Branch("m0PLeta",&m0PLeta);
      m0PartonTree->Branch("m0PLphi",&m0PLphi);
      m0PartonTree->Branch("m0PLzg",&m0PLzg);
      m0PartonTree->Branch("m0PLrg",&m0PLrg);
      m0PartonTree->Branch("m0PLncons",&m0PLncons);
      m0PartonTree->Branch("m0PLnconsg",&m0PLnconsg);
      m0PartonTree->Branch("pthat",&pthat,"pthat/d");
      m0PartonTree->Branch("mcweight",&mcweight,"mcweight/d");


      TrigRecTree = new TTree("TrigRecTree","Recoil spectra for various triggers");
      trigpT.clear();
      recpT.clear();
      
      TrigRecTree->Branch("trigpT",&trigpT);
      TrigRecTree->Branch("recpT",&recpT);
      TrigRecTree->Branch("pthat", &pthat,"pthat/d");
      TrigRecTree->Branch("mcweight",&mcweight,"mcweight/d");

      ResultTree=new TTree("ResultTree","Result Jets");
      //! clear the vectors for each event! 
      conspT.clear();
      consDist.clear();
      consGirth.clear();
      consM.clear();
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
      rec_zg.clear();
      rec_rg.clear();
      rec_pt.clear();
      rec_mj.clear();
      rec_sdpt.clear();
      rec_sdmj.clear();
      rec_hard_mass.clear();
      rec_soft_mass.clear();
      
      ResultTree->Branch("cent", &cent,"cent/d");
      ResultTree->Branch("mcweight", &mcweight,"mcweight/d");
      ResultTree->Branch("pthat", &pthat,"pthat/d");
      ResultTree->Branch("eventID",&eventID,"eventID/d");
      ResultTree->Branch("conspT",&conspT);
      ResultTree->Branch("consDist",&consDist);
      ResultTree->Branch("consGirth",&consGirth);
      ResultTree->Branch("consM",&consM);
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
      ResultTree->Branch("rec_zg", &rec_zg);
      ResultTree->Branch("rec_rg", &rec_rg);      
      ResultTree->Branch("rec_pt", &rec_pt);
      ResultTree->Branch("rec_mj", &rec_mj);      
      ResultTree->Branch("rec_sdpt", &rec_sdpt);
      ResultTree->Branch("rec_sdmj", &rec_sdmj);
      ResultTree->Branch("rec_hard_mass", &rec_hard_mass);
      ResultTree->Branch("rec_soft_mass", &rec_soft_mass);      
      twosubjet_z.clear();
      twosubjet_theta.clear();
      threesubjet_z.clear();
      threesubjet_theta.clear();
      ResultTree->Branch("twosubjet_z", &twosubjet_z);
      ResultTree->Branch("twosubjet_theta", &twosubjet_theta);
      ResultTree->Branch("threesubjet_z", &threesubjet_z);
      ResultTree->Branch("threesubjet_theta", &threesubjet_theta);


      HCResultTree=new TTree("HCResultTree","Result Jets");
      //! clear the vectors for each event! 
      HCjetpT.clear();
      HCjetM.clear();
      HCsdjetpT.clear();
      HCsdjetM.clear();
      HCjeteta.clear();
      HCjetphi.clear();
      HCzg.clear();
      HCrg.clear();
      
      HCResultTree->Branch("cent", &cent,"cent/d");
      HCResultTree->Branch("mcweight", &mcweight,"mcweight/d");
      HCResultTree->Branch("pthat", &pthat,"pthat/d");
      HCResultTree->Branch("jetpT", &HCjetpT);
      HCResultTree->Branch("jeteta", &HCjeteta);
      HCResultTree->Branch("jetphi", &HCjetphi);
      HCResultTree->Branch("jetM", &HCjetM);
      HCResultTree->Branch("sdjetpT", &HCsdjetpT);
      HCResultTree->Branch("sdjetM", &HCsdjetM);
      HCResultTree->Branch("zg", &HCzg);
      HCResultTree->Branch("rg", &HCrg);
      HCtwosubjet_z.clear();
      HCtwosubjet_theta.clear();
      HCthreesubjet_z.clear();
      HCthreesubjet_theta.clear();
      HCResultTree->Branch("twosubjet_z", &HCtwosubjet_z);
      HCResultTree->Branch("twosubjet_theta", &HCtwosubjet_theta);
      HCResultTree->Branch("threesubjet_z", &HCthreesubjet_z);
      HCResultTree->Branch("threesubjet_theta", &HCthreesubjet_theta);


      MCResultTree=new TTree("MCResultTree","Result Jets");
      //! clear the vectors for each event! 
      MCjetpT.clear();
      MCjetM.clear();
      MCsdjetpT.clear();
      MCsdjetM.clear();
      MCjeteta.clear();
      MCjetphi.clear();
      MCzg.clear();
      MCrg.clear();
      
      MCResultTree->Branch("cent", &cent,"cent/d");
      MCResultTree->Branch("mcweight", &mcweight,"mcweight/d");
      MCResultTree->Branch("pthat", &pthat,"pthat/d");
      MCResultTree->Branch("jetpT", &MCjetpT);
      MCResultTree->Branch("jeteta", &MCjeteta);
      MCResultTree->Branch("jetphi", &MCjetphi);
      MCResultTree->Branch("jetM", &MCjetM);
      MCResultTree->Branch("sdjetpT", &MCsdjetpT);
      MCResultTree->Branch("sdjetM", &MCsdjetM);
      MCResultTree->Branch("zg", &MCzg);
      MCResultTree->Branch("rg", &MCrg);
      MCtwosubjet_z.clear();
      MCtwosubjet_theta.clear();
      MCthreesubjet_z.clear();
      MCthreesubjet_theta.clear();
      MCResultTree->Branch("twosubjet_z", &MCtwosubjet_z);
      MCResultTree->Branch("twosubjet_theta", &MCtwosubjet_theta);
      MCResultTree->Branch("threesubjet_z", &MCthreesubjet_z);
      MCResultTree->Branch("threesubjet_theta", &MCthreesubjet_theta);

      
    }//! init

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      
      const double weight = handler().crossSection();

      //! clear the vectors for each event! 
      recpT.clear();
      trigpT.clear();

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
      
      m0PLpt.clear();
      m0PLm.clear();
      m0PLptg.clear();
      m0PLmg.clear();
      m0PLeta.clear();
      m0PLphi.clear();
      m0PLzg.clear();
      m0PLrg.clear();
      m0PLncons.clear();
      m0PLnconsg.clear();
      

      conspT.clear();
      consDist.clear();
      consGirth.clear();
      consM.clear();
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
      twosubjet_z.clear();
      twosubjet_theta.clear();
      threesubjet_z.clear();
      threesubjet_theta.clear();
      rec_zg.clear();
      rec_rg.clear();
      rec_pt.clear();
      rec_mj.clear();
      rec_sdpt.clear();
      rec_sdmj.clear();
      rec_hard_mass.clear();
      rec_soft_mass.clear();

      HCjetpT.clear();
      HCjetM.clear();
      HCsdjetpT.clear();
      HCsdjetM.clear();
      HCjeteta.clear();
      HCjetphi.clear();
      HCzg.clear();
      HCrg.clear();
      HCtwosubjet_z.clear();
      HCtwosubjet_theta.clear();
      HCthreesubjet_z.clear();
      HCthreesubjet_theta.clear();

      MCjetpT.clear();
      MCjetM.clear();
      MCsdjetpT.clear();
      MCsdjetM.clear();
      MCjeteta.clear();
      MCjetphi.clear();
      MCzg.clear();
      MCrg.clear();
      MCtwosubjet_z.clear();
      MCtwosubjet_theta.clear();
      MCthreesubjet_z.clear();
      MCthreesubjet_theta.clear();
      
      const double qscale = event.genEvent()->event_scale();
      //MUST BE CHANGED IF MORE PYTHIA EVENTS ARE RUN!!!!!!!!!!!!!!!!!!
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
      
      // const double ecent = (event.genEvent()->heavy_ion()?event.genEvent()->heavy_ion()->impact_parameter():-1.);

      // int centbin(-1);
      // for (size_t i = 0; i < _ncentbins; ++i) {
      // 	if (ecent >= _centedges[i] && ecent < _centedges[i+1]) centbin=i;
      // }
      int centbin=0;
      
      // Cut cuts = Cuts::etaIn(-0.6, 0.6) & (Cuts::pT > _pTCut*GeV);
      // const FastJets& Ajets = applyProjection<FastJets>(event, "Jets");
      // const Jets ajets = Ajets.jetsByPt(cuts);
      // if(ajets.size()==0) vetoEvent;
      // if(printDebug)
      // 	std::cout<<"There are "<<ajets.size()<<" inclusive jets in this event"<<std::endl;

      /*
      //! get the scattered particles from the gen event:
      vector <HepMC::GenParticle*> pscat; 
      foreach (const HepMC::GenParticle* p, particles(event.genEvent())) {
	if(p->status() == 3)
	  pscat.push_back((HepMC::GenParticle*)p);
      }
      
      bool doSubtraction = true;
      if(pscat.size() == 0)
       	doSubtraction = false;

      if(printDebug)
	std::cout<<"doSubtraction = "<<doSubtraction<<std::endl;
      
      //! initialize the grid
      vector <vector<candidate> > grid_hardcore;
      for(int x = 0; x<_Nbounds_eta; ++x){
	vector <candidate> xgrid;
	for(int y = 0; y<_Nbounds_phi; ++y){
	  candidate test;
	  test.etaMin = _etabins[x];
	  test.phiMin = _phibins[y];
	  test.etaMax = _etabins[x+1];
	  test.phiMax = _phibins[y+1];
	  test.sumnegpT = 0.0;
	  xgrid.push_back(test);
	}
	grid_hardcore.push_back(xgrid);
	xgrid.clear();
      }

      if(printDebug)
	std::cout<<"after grid initialization"<<std::endl;
      
      
      const ParticleVector& FS = applyProjection<FinalState>(event, "FS").particlesByPt();     

      //! first add the hard core objects in the event. 	      
      foreach ( const Particle& p, FS) {
	if(p.pt() >= 2.*GeV){
	  double ceta = p.eta();
	  double cphi = p.phi(MINUSPI_PLUSPI);
	  std::pair<int,int> candpos = FindCandidateBin(ceta, cphi);
	  if(candpos.first!=-1 && candpos.second!=-1){
	    grid_hardcore[candpos.first][candpos.second].candMom+=p.momentum();
	    grid_hardcore[candpos.first][candpos.second].sumpT+=p.pt();
	  }
	}
      }

      if(printDebug)
	std::cout<<"filling the grid"<<std::endl;
      
      //! build up the pseudojets to do the clustering
      PseudoJets pJet_hardcore;
      for(int x = 0; x<_Nbounds_eta; ++x){
	for(int y = 0; y<_Nbounds_phi; ++y){
	  FourMomentum jet_hardcore = grid_hardcore[x][y].candMom;
	  double px, py, pz, E;
	  px = jet_hardcore.px();
	  py = jet_hardcore.py();
	  pz = jet_hardcore.pz();
	  E  = jet_hardcore.E();	  
	  PseudoJet part_hardcore(px,py,pz,E);
	  if(part_hardcore.pt() != 0 && part_hardcore.E() != 0)
	    pJet_hardcore.push_back(part_hardcore);
	}
      }//!grid loop

      grid_hardcore.clear();
      
      if(printDebug)
	std::cout<<""<<std::endl;
      */

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
      
      //shouldn't have to check for triggers in the event since we do that in running over the pythia!
      /*
      bool trig_in_event = 0;
      foreach ( const Particle & p, FS) {
	if (trig_in_event) { break; } //if we have a trigger, no need to keep looking
	if(p.pt() >= 5.4*GeV && p.pid() == 111){ //found a pi0 trigger
	  cout << "HAVE A TRIGGER!" << endl;
	  trig_in_event = 1;
	  break;
	}
      }
      cout << "This should always say 1: " << trig_in_event << endl;
      if (!trig_in_event) {return;}
      */
      Cut cuts = Cuts::etaIn(-1, 1) & (Cuts::pT > _pTCut*GeV); 
      
      double highest_pt_pi0 = -1; double trig_phi = -9999;
      //      cout << "starting event loop" << endl;
      foreach (const Particle &p, event.allParticles(cuts)/*FS*/) {//now that we know there was a trigger, time to find the highest pT one.
	//cout << "PID: " << p.pid() << endl;
        if (p.pid() == 111/*p.charge() == 0*/ && p.pt() > highest_pt_pi0) {
	  highest_pt_pi0 = p.pt(); trig_phi = p.phi();
	} 
      }
      trigpT.push_back(highest_pt_pi0); //we've guaranteed the event has a neutral trigger, so should never be pushing_back with anything less than 5.4.

      fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, _jetR);
      
      
      //QUARK v. GLUON PREP: figuring out what/where the hard scattered partons are
      vector<int> identity; vector<PseudoJet> hardparton;
      //for (int i = 0; i < event.genEvent()->particles_size(); ++ i) {
      foreach (const Particle &p, event.allParticles()) {
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
      // if(recoJets_hardcore.size()>=2){
      // 	//! fill Aj_hardcore
      // 	if(recoJets_hardcore.at(0).pt() >= 20 &&
      // 	   recoJets_hardcore.at(1).pt() >= 10 &&
      // 	   recoJets_hardcore.at(0).delta_phi_to(recoJets_hardcore.at(1)) > (M_PI - _jetR)){
      // 	  // _h_aj_hardcore[centbin]->fill((recoJets_hardcore.at(0).pt() - recoJets_hardcore.at(1).pt())/
      // 	  // 			      (recoJets_hardcore.at(0).pt() + recoJets_hardcore.at(1).pt()), weight);
      // 	  hAJ_hardcore[centbin]->Fill((recoJets_hardcore.at(0).pt() - recoJets_hardcore.at(1).pt())/
      // 				      (recoJets_hardcore.at(0).pt() + recoJets_hardcore.at(1).pt()), weight);
      // 	  _NevHardcore[centbin]+=weight;
      // 	}
      // }
      
      /*
      //! initialize the low core grid
      vector <vector<candidate> > grid_lowcore;
      for(int x = 0; x<_Nbounds_eta; ++x){
	vector <candidate> xgrid;
	for(int y = 0; y<_Nbounds_phi; ++y){
	  candidate test;
	  test.etaMin = _etabins[x];
	  test.phiMin = _phibins[y];
	  test.etaMax = _etabins[x+1];
	  test.phiMax = _phibins[y+1];
	  test.sumnegpT = 0.0;
	  xgrid.push_back(test);
	}
	grid_lowcore.push_back(xgrid);
	xgrid.clear();
      }
      if(printDebug)
	std::cout<<"Going to fill the grid with low core constituents"<<std::endl;
      
      foreach ( const Particle& p, FS) {
	if(p.pt() >= 0.2*GeV){
	  double ceta = p.eta();
	  double cphi = p.phi(MINUSPI_PLUSPI);
	  std::pair<int,int> candpos = FindCandidateBin(ceta, cphi);
	  if(candpos.first!=-1 && candpos.second!=-1){
	    grid_lowcore[candpos.first][candpos.second].candMom+=p.momentum();
	    grid_lowcore[candpos.first][candpos.second].sumpT+=p.pt();
	  }
	}
      }

      if(printDebug)
	std::cout<<"doing subtraction in the grid"<<std::endl;
      
      //! Now loop over the scattering centers that are near jets and add that as BKG
      if(doSubtraction){
	for(unsigned ip = 0; ip<pscat.size(); ++ip){
	  double px,py,pz,E;
	  px = pscat[ip]->momentum().px();
	  py = pscat[ip]->momentum().py();
	  pz = pscat[ip]->momentum().pz();
	  E  = pscat[ip]->momentum().e();
	  PseudoJet part(px,py,pz,E);
	  bool iJ = false;
	  foreach(const PseudoJet j, ajets){
	    double delR = deltaR(j.eta(), j.phi_std(), part.eta(), part.phi_std());
	    if(delR < _jetR){
	      iJ = true;
	      break;
	    }
	  }
	  if(iJ){
	    std::pair<int,int> candpos = FindCandidateBin(part.eta(), part.phi_std());
	    if(candpos.first!=-1 && candpos.second!=-1){
	      grid_lowcore[candpos.first][candpos.second].bkgMom+=FourMomentum(E, px, py, pz);
	      grid_lowcore[candpos.first][candpos.second].sumpT-=part.pt();
	    }
	  }
	}// scattering center loop
      }// dosub

      //! build up the pseudojets to do the clustering
      PseudoJets pJet_sub;
      for(int x = 0; x<_Nbounds_eta; ++x){
	for(int y = 0; y<_Nbounds_phi; ++y){
	  FourMomentum jet_sub = grid_lowcore[x][y].candMom;
	  if(doSubtraction)
	    jet_sub -= grid_lowcore[x][y].bkgMom;
	  double px, py, pz, E;
	  if(grid_lowcore[x][y].candMom.pT() - grid_lowcore[x][y].bkgMom.pT() < 0.0){
	    grid_lowcore[x][y].sumnegpT += grid_lowcore[x][y].bkgMom.pT() - grid_lowcore[x][y].candMom.pT();
	    px = 0.0; py = 0.0; pz = 0.0; E = 0.0;
	  }else{
	    px = jet_sub.px();
	    py = jet_sub.py();
	    pz = jet_sub.pz();
	    E  = jet_sub.E();
	  }
	  PseudoJet part_sub(px,py,pz,E);
	  if(part_sub.pt() != 0 && part_sub.E() != 0)
	    pJet_sub.push_back(part_sub);
	}
      }//!grid loop

      grid_lowcore.clear();

      if(printDebug)
	std::cout<<"just found all the subtracted jets "<<std::endl;
      */

      //! Get all the objects > 0.2 GeV for the matched jets
      PseudoJets pJet_sub; PseudoJets pJet_ch; PseudoJet intermediate; PseudoJet m0intermediate;
      //! first add the hard core objects in the event. 	      
      foreach ( const Particle& p, FS) {
	//if(p.pt() >= 0.2*GeV && fabs(p.eta()) < 1){
	intermediate = PseudoJet(p.px(), p.py(), p.pz(), p.E());
	//if (p.charge() != 0) {
	//	  std::cout << p.mass() << std::endl;
	//intermediate.reset_PtYPhiM(p.pt(),p.rap(),p.phi(), p.mass());//,m_chpion);
	
	//if (p.mass() > 1) {cout << "MASSIVE PARTICLE BELOW!" << endl;}
	//if (p.mass() > 9) {cout << "CHONKYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY v" << endl;}
	//cout << p.pid() << ": " << p.mass() << " " << pdg->GetParticle(p.pid())->Mass() << endl;
	//if (p.mass() > 1) {cout << "MASSIVE PARTICLE ABOVE!" << endl;}
	//if (p.mass() > 9) {cout << "CHONKYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY ^" << endl;}

	intermediate.reset_PtYPhiM(p.pt(),p.rap(),p.phi(),(double) pdg->GetParticle(p.pid())->Mass()); //PDG MASSES!!!
       
	//}
	pJet_sub.push_back(intermediate);
	if (p.charge() != 0) {
	  pJet_ch.push_back(PseudoJet(p.px(), p.py(), p.pz(), p.E()));
	}
	//}
      }     

      PseudoJets PLjet; PseudoJets m0PLjet;
      foreach ( const Particle & p, FP) {
	intermediate = PseudoJet(p.px(), p.py(), p.pz(), p.E());
	m0intermediate = PseudoJet(p.px(), p.py(), p.pz(), p.E());
	if (p.mass() > 1) {cout << "MASSIVE PARTICLE BELOW!" << endl;}
	if (p.mass() > 9) {cout << "CHONKYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY v" << endl;}
	cout << p.pid() << ": " << p.mass() << " " << pdg->GetParticle(p.pid())->Mass() << endl;
	if (p.mass() > 1) {cout << "MASSIVE PARTICLE ABOVE!" << endl;}
	if (p.mass() > 9) {cout << "CHONKYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY ^" << endl;}
	
	intermediate.reset_PtYPhiM(p.pt(),p.rap(),p.phi(),pdg->GetParticle(p.pid())->Mass()); //PDG MASSES!!!
	m0intermediate.reset_PtYPhiM(p.pt(),p.rap(),p.phi(),0); //MASSLESS!!!
	
        //intermediate.reset_PtYPhiM(p.pt(),p.rap(),p.phi(), 0);//p.mass());//,m_chpion);//!MASSLESS PARTONS!!!
        PLjet.push_back(intermediate);
	m0PLjet.push_back(m0intermediate);
      }
      
      vector<PseudoJet> pJet_cut = spart(pJet_sub);
      vector<PseudoJet> pJet_ch_cut = spart(pJet_ch);
      vector<PseudoJet> PL_cut = spart(PLjet);
      vector<PseudoJet> m0PL_cut = spart(m0PLjet);
      
      
      fastjet::ClusterSequence cs_sub(pJet_cut, jet_def);
      fastjet::ClusterSequence cs_ch(pJet_ch_cut, jet_def);
      fastjet::ClusterSequence cs_PL(PL_cut, jet_def);
      fastjet::ClusterSequence cs_m0PL(m0PL_cut, jet_def);

      PseudoJets recoJets_bs = sorted_by_pt(cs_sub.inclusive_jets());
      PseudoJets recoJets_bs_ch = sorted_by_pt(cs_ch.inclusive_jets());
      PseudoJets recoJets_bs_pl = sorted_by_pt(cs_PL.inclusive_jets());
      PseudoJets recoJets_bs_m0pl = sorted_by_pt(cs_m0PL.inclusive_jets());

      PseudoJets recoJets_bkgsub = select_both_hi(recoJets_bs); //want to use these for response construction
      PseudoJets recoJets_bkgsub_ch = select_both_hi(recoJets_bs_ch);
      PseudoJets recoJets_PL = select_both_hi(recoJets_bs_pl);
      PseudoJets recoJets_m0PL = select_both_hi(recoJets_bs_m0pl);
      
      /*
      if(recoJets_bkgsub.size()==0){
	// std::cout<<"NO low core RECONSTRUCTED JETS or only 1 jet"<<std::endl;
	vetoEvent;
      }
      if(recoJets_bkgsub_ch.size()==0) {
	//if there are no charged jets, veto the event.
	vetoEvent;
      }
      */

      if(printDebug){
	std::cout<<"bkg subtracted Jets: "<<std::endl;
	foreach (const PseudoJet jet, recoJets_bkgsub){
	  std::cout<<"    "<<jet.pt()<<std::endl;
	}
      }
      
      // vector<int> jetMatch;
      // //!Now we have recoJets_hardcore and recoJets_bkgsub. We need to match the appropriate jets to get AJ_matched.
      // if(recoJets_bkgsub.size() >= recoJets_hardcore.size()){
      // 	jetMatch.resize(recoJets_bkgsub.size());
      // 	int leadjetCounter = 0;
      // 	foreach (const PseudoJet jet1, recoJets_bkgsub){
      // 	  int counter = 0;
      // 	  foreach (const PseudoJet jet2, recoJets_hardcore){	    
      // 	    double delR = deltaR(jet1.eta(), jet1.phi_std(), jet2.eta(), jet2.phi_std());
      // 	    if(delR < _jetR){
      // 	      jetMatch.push_back(counter);
      // 	      if(printDebug){
      // 		std::cout<<"bkg subtracted jet (pt, eta, phi): # "<<leadjetCounter<<" ("<<jet1.pt()<<", "<<jet1.eta()<<", "<<jet1.phi()
      // 			 <<") is Matched with hard core jet: # "<<counter<<" ("<<jet2.pt()<<", "<<jet2.eta()<<", "<<jet2.phi()<<") "<<std::endl;
      // 	      }
      // 	    }
      // 	    counter++;
      // 	  }
      // 	  leadjetCounter++;
      // 	}
      // }//! making sure that subtracted jet collection is same or larger than hard core jet collection
      // // if(jetMatch.size()>1){
      // // 	//! give possibility for the jet matching to change
      // // 	if((jetMatch[0] == 0 || jetMatch[0] == 1) && (jetMatch[1] == 0 || jetMatch[1] == 1)){
      // // _h_aj_matched[centbin]->fill((recoJets_bkgsub.at(0).pt() - recoJets_bkgsub.at(1).pt())/
      // // 				   (recoJets_bkgsub.at(0).pt() - recoJets_bkgsub.at(1).pt()), weight);
      // if(recoJets_bkgsub.size()>=2){
      // 	hAJ_matched[centbin]->Fill((recoJets_bkgsub.at(0).pt() - recoJets_bkgsub.at(1).pt())/
      // 				   (recoJets_bkgsub.at(0).pt() + recoJets_bkgsub.at(1).pt()), weight);
      // 	_NevMatched[centbin]+=weight;
      // }
      // // 	}
      // // }
      
      //! loop over the jets now:
      // int jetcounter = 0;

      // if(recoJets_bkgsub.size() >= 2) {

      // 	pt1 = recoJets_bkgsub.at(0).pt();
      // 	const std::pair<double,double> sd1 = SoftDrop(recoJets_bkgsub.at(0), _jetR, z_cut, beta);
      // 	zg1 = sd1.first;
      // 	delR1 = sd1.second;
      // 	m1 = recoJets_bkgsub.at(0).m();
      // 	pt2 = recoJets_bkgsub.at(1).pt();
      // 	const std::pair<double,double> sd2 = SoftDrop(recoJets_bkgsub.at(1), _jetR, z_cut, beta);
      // 	zg2 = sd2.first;
      // 	delR2 = sd2.second;
      // 	m2 = recoJets_bkgsub.at(1).m();
      cent = centbin;
      mcweight = weight;
      pthat = qscale;
      // 	if(pt1 > 20.0 && pt2 > 10.0){
      // 	  Nsplitting->Fill(cent, pt1, zg1, delR1, m1, pt2, zg2, delR2, m2);
      // 	}
	
      // }
      //cout << "TEST!" << endl;
      foreach(const PseudoJet jet, recoJets_bkgsub_ch) {
	//get highest pT pi0 (or neutral particle) from the list of particles above. Find its phi. Here, only take jets +- 1/4 from back-to-back from it.
	double recoil_side = trig_phi + pi;
	double recoil_side_modded = fmod(recoil_side, (2.0*pi)); //trig phi is from 0 to 2pi
	double recoil_side_modded_sub = recoil_side_modded - 0.25;
	double recoil_side_modded_add = recoil_side_modded + 0.25;
	double recoil_side_modded_sub_modded = fmod(recoil_side_modded_sub,(2.0*pi));
	double recoil_side_modded_add_modded = fmod(recoil_side_modded_add,(2.0*pi));
	double smaller = -9999; double bigger = -9999;
	/*
	cout << "trig phi: " << trig_phi << endl;
	cout << "trig phi + pi: " << recoil_side << endl;
	cout << "(trig phi + pi) % 2pi: " << recoil_side_modded << endl;
	cout << "     - 0.25: " << recoil_side_modded_sub << endl;
	cout << "     + 0.25: " << recoil_side_modded_add << endl;
	cout << "          % 2pi: " << recoil_side_modded_sub_modded << endl;
	cout << "          % 2pi: " << recoil_side_modded_add_modded << endl;
	cout << "jet phi: " << jet.phi() << endl;
	*/
	if (recoil_side_modded_sub_modded != recoil_side_modded_sub || recoil_side_modded_add_modded != recoil_side_modded_add) {
	  if (recoil_side_modded_sub_modded > recoil_side_modded_add_modded) {
	    smaller = recoil_side_modded_add_modded; bigger = recoil_side_modded_sub_modded;
	  }
	  else {
	    smaller = recoil_side_modded_sub_modded; bigger = recoil_side_modded_add_modded;
	  }
	  if ((0 < jet.phi() && jet.phi() < smaller) || (bigger < jet.phi() && jet.phi() < 2.0*pi)) {
	    //  cout << "pushing back" << endl;
	    recpT.push_back(jet.pt());
	  }
	}
	else {
	  if ( (jet.phi() > recoil_side_modded - 0.25) && (jet.phi() < recoil_side_modded + 0.25)) {
	    //cout << "pushing back" << endl;
	    recpT.push_back(jet.pt()); //can take more than one recoil if requirements are satisfied
	  }
	}
      }
      if (recpT.size() == 0) {
	recpT.push_back(-9999);
      }
      
      TrigRecTree->Fill();


      //HEY LOOK HERE!!!!!!!!!!!!!!!!!!!!!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
      
      foreach(const PseudoJet jet, recoJets_m0PL) {
	m0PLpt.push_back(jet.pt());
	m0PLm.push_back(jet.m());
	m0PLeta.push_back(jet.eta());
	m0PLphi.push_back(jet.phi());
	m0PLncons.push_back(jet.constituents().size());
	
	fastjet::contrib::RecursiveSymmetryCutBase::SymmetryMeasure  symmetry_measure = fastjet::contrib::RecursiveSymmetryCutBase::scalar_z;
	fastjet::contrib::SoftDrop sd(beta, z_cut, symmetry_measure, _jetR);
        PseudoJet sd_jet = sd(jet);
	
	if(sd_jet != 0){

          m0PLptg.push_back(sd_jet.pt());
          m0PLmg.push_back(sd_jet.m());
          double z = sd_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry();
          double r = sd_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R();

          m0PLzg.push_back(z);
          m0PLrg.push_back(r);
	  m0PLnconsg.push_back(sd_jet.constituents().size());

	}
	
      }

      foreach(const PseudoJet jet, recoJets_PL) {
	PLpt.push_back(jet.pt());
	PLm.push_back(jet.m());
	PLeta.push_back(jet.eta());
	PLphi.push_back(jet.phi());
	PLncons.push_back(jet.constituents().size());
	
	fastjet::contrib::RecursiveSymmetryCutBase::SymmetryMeasure  symmetry_measure = fastjet::contrib::RecursiveSymmetryCutBase::scalar_z;
	fastjet::contrib::SoftDrop sd(beta, z_cut, symmetry_measure, _jetR);
        PseudoJet sd_jet = sd(jet);
	
	if(sd_jet != 0){

          PLptg.push_back(sd_jet.pt());
          PLmg.push_back(sd_jet.m());
          double z = sd_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry();
          double r = sd_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R();

          PLzg.push_back(z);
          PLrg.push_back(r);
	  PLnconsg.push_back(sd_jet.constituents().size());

	}
	
      }
      
      //debugging:
      //std::cout << "number of jets above 0.2 GeV in the event: " << recoJets_bkgsub.size() << std::endl;
      foreach(const PseudoJet jet, recoJets_bkgsub){// CHANGE WHEN GOING BETWEEN FULL JETS AND CHARGED JETS
	//! Inclusive
	hJetpT[centbin]->Fill(jet.pt(), weight);
	hJetEta[centbin]->Fill(jet.eta(), weight);
	hJetPhi[centbin]->Fill(jet.phi(), weight);
	
	//GEOMETRICALLY MATCH JET AND PARTON
	bool match = 0;
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
	//vector<PseudoJet> jet_cons = jet.constituents();
	foreach (const PseudoJet cons, jet.constituents()) {
	  cons_pt.push_back(cons.pt());
	  double consR = fabs(cons.delta_R(jet));
	  cons_dist.push_back(consR);
	  double part_girth = consR * cons.pt() / (double) jet.pt();
	  cons_girth.push_back(part_girth);
	  cons_mass.push_back(cons.m());
	  jet_girth += part_girth;
	  //~~~
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

	fastjet::contrib::RecursiveSymmetryCutBase::SymmetryMeasure  symmetry_measure = fastjet::contrib::RecursiveSymmetryCutBase::scalar_z;
	fastjet::contrib::SoftDrop sd(beta, z_cut, symmetry_measure, _jetR);
	PseudoJet sd_jet = sd(jet);

	//! Loop over the jet constituent and make three smaller jets
	//! apply a zcut on the smaller of the subjets
	
	vector<PseudoJet> jetconsts = jet.constituents();
	
	fastjet::JetDefinition jetd(fastjet::antikt_algorithm, _jetR*0.5);
	fastjet::ClusterSequence clust_seq_full(jetconsts, jetd);

	// cout<<"at the const sub jets section"<<endl;

	vector<PseudoJet> subjets= sorted_by_pt(clust_seq_full.inclusive_jets());
	if(subjets.size() >= 3){
	  
	  //! three subjets - perform some kinda grooming
	  // vector<PseudoJet> threesubjets = sorted_by_pt(clust_seq_full.exclusive_jets(3));
	  PseudoJet comblead = subjets.at(0) + subjets.at(1);
	  double zsj = subjets.at(2).pt()/(comblead.pt()+subjets.at(2).pt());
	  double rsj = comblead.delta_R(subjets.at(2));
	  if(zsj > 0.1){
	    threesubjet_z.push_back(zsj);
	    threesubjet_theta.push_back(rsj);
	  }else {
	    threesubjet_z.push_back(-999);
	    threesubjet_theta.push_back(-999);
	  }
	}else{
	  threesubjet_z.push_back(-999);
	  threesubjet_theta.push_back(-999);
	}

	if(subjets.size() >= 2){
	  //! three subjets - no grooming. just straight up two subjets 
	  double zsj = subjets.at(1).pt()/(subjets.at(0).pt()+subjets.at(1).pt());
	  double rsj = subjets.at(0).delta_R(subjets.at(1));
	  if(zsj > 0.1){
	    twosubjet_z.push_back(zsj);
	    twosubjet_theta.push_back(rsj);
	  }else {
	    twosubjet_z.push_back(-999);
	    twosubjet_theta.push_back(-999);
	  }
	}else{
	  twosubjet_z.push_back(-999);
	  twosubjet_theta.push_back(-999);
	}
       	

	if(sd_jet != 0){
	  
	  sdjetpT.push_back(sd_jet.pt());
	  sdjetM.push_back(sd_jet.m());
	  double z = sd_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry();
	  double r = sd_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R();
	  
	  zg.push_back(z);
	  rg.push_back(r);
	  
	  vector<double> reczg;
	  vector<double> recrg;
	  vector<double> recpt;
	  vector<double> recmj;
	  vector<double> recsdpt;
	  vector<double> recsdmj;
	  vector<double> rechardmass;
	  vector<double> recsoftmass;
          
	  PseudoJet jj, j1, j2;
	  jj = jet;
	  
	  // while(jj.fastjet::PseudoJet::has_parents(j1, j2)){
	  bool continueRec = true;

	  while(continueRec){

	    PseudoJet sd_rec = sd(jj);

	    if(sd_rec == 0){
	      if(printDebug)std::cout<<" *******BREAK******* Recursive branch grooming fails, z< 0.1  *******BREAK*******"<<std::endl;
	      break;
	    }

	    if(sd_rec.fastjet::PseudoJet::has_parents(j1, j2))
	      continueRec = true;
	    else{
	      continueRec = false;
	      break;
	    }
	    //! sort the subjets to look at the leading prong 
	    if(j1.pt() < j2.pt()) {
	      PseudoJet test;
	      test = j1;
	      j1 = j2;
	      j2 = test;
	    }

	    double dR = sd_rec.structure_of<fastjet::contrib::SoftDrop>().delta_R();
	    double zrel = sd_rec.structure_of<fastjet::contrib::SoftDrop>().symmetry();
	    double sd_pt = sd_rec.pt();
	    double sd_mj = sd_rec.m();

	    // if ( jj.pt() < j1.pt() ) std::cout << "WTF! " << jj.pt() << "  " << j1.pt() << endl;

	    
	    if(zrel > 0.1) {
	    
	      recpt.push_back(jj.pt());
	      recmj.push_back(jj.m());
	      rechardmass.push_back(j1.m());
	      recsoftmass.push_back(j2.m());
	      
	      // collect info and fill in the vectors
	      // double dR = j1.fastjet::PseudoJet::delta_R(j2);
	      // double zrel = j2.pt()/(j1.pt() + j2.pt());
	      
	      reczg.push_back(zrel);
	      recrg.push_back(dR);
	      recsdpt.push_back(sd_pt);
	      recsdmj.push_back(sd_mj);
	      
	      // if(zrel > 0.1){
	      
	      //   reczg.push_back(zrel);
	      //   recrg.push_back(dR);
	      
	      // }else {
	      
	      //   reczg.push_back(-999);
	      //   recrg.push_back(-999);
	      
	      //   break;
	      
	      // }
	    }else{
	      recpt.push_back(-999);
	      recmj.push_back(-999);
	      rechardmass.push_back(-999);
	      recsoftmass.push_back(-999);
	      
	      // collect info and fill in the vectors
	      // double dR = j1.fastjet::PseudoJet::delta_R(j2);
	      // double zrel = j2.pt()/(j1.pt() + j2.pt());
	      
	      reczg.push_back(-999);
	      recrg.push_back(-999);
	      recsdpt.push_back(-999);
	      recsdmj.push_back(-999);	      

	    }

	    jj = j1;

	  }
	  
	  rec_zg.push_back(reczg);
	  rec_rg.push_back(recrg);
	  rec_pt.push_back(recpt);
	  rec_mj.push_back(recmj);
	  rec_sdpt.push_back(recsdpt);
	  rec_sdmj.push_back(recsdmj);
	  rec_hard_mass.push_back(rechardmass);
	  rec_soft_mass.push_back(recsoftmass);

	}else {
	  sdjetpT.push_back(-999);
	  sdjetM.push_back(-999);
	  zg.push_back(-999);
	  rg.push_back(-999);
	  vector<double> reczg;
	  vector<double> recrg;
	  reczg.push_back(-999);
	  recrg.push_back(-999);
	  rec_zg.push_back(reczg);
	  rec_rg.push_back(recrg);

	  vector<double> recpt;
	  vector<double> recmj;
	  recpt.push_back(-999);
	  recmj.push_back(-999);
	  rec_pt.push_back(recmj);
	  rec_mj.push_back(recpt);	  

	  vector<double> recsdpt;
	  vector<double> recsdmj;
	  recsdpt.push_back(-999);
	  recsdmj.push_back(-999);
	  rec_sdpt.push_back(recsdmj);
	  rec_sdmj.push_back(recsdpt);	  
	  
	  vector<double> rechardmass;
	  vector<double> recsoftmass;
	  rechardmass.push_back(-999);
	  recsoftmass.push_back(-999);
	  rec_hard_mass.push_back(rechardmass);
	  rec_soft_mass.push_back(recsoftmass);	  
	}
	
	int ptbin = -1;
	for (size_t j = 0; j < _nptbins; ++j) {
	  if (jet.pt() >= _ptedges[j] && jet.pt() < _ptedges[j+1]) {
	    ptbin = j;
	  }
	}
	if(ptbin == -1) continue;
	const std::pair<double,double> zg = SoftDrop(jet, _jetR, z_cut, beta);
	if(zg.second > 0.0){
	  // _h_zg[centbin][ptbin]->fill(zg.first, weight);
	  hzg[centbin][ptbin]->Fill(zg.first, weight);
	  hRg[centbin][ptbin]->Fill(zg.second, weight);
	  hMj[centbin][ptbin]->Fill(jet.m(), weight);
	  _Nj[centbin][ptbin]+=weight;
	}
	//! Leading Jet 
	// if(jetcounter == 0){
	//   int ptleadbin = -1;
	//   for (size_t j = 0; j < _nptleadbins; ++j) {
	//     if (jet.pt() >= _ptleadedges[j] && jet.pt() < _ptleadedges[j+1]) {
	//       ptleadbin = j;
	//     }	
	//   }
	//   if(ptleadbin!=-1){
	//     if(zg.second > 0.0){
	//       // _h_zg_lead[centbin][ptleadbin]->fill(zg.first, weight);
	//       hzg_lead[centbin][ptleadbin]->Fill(zg.first, weight);
	//       hRg_lead[centbin][ptleadbin]->Fill(zg.second, weight);
	//       hMj_lead[centbin][ptleadbin]->Fill(jet.m(), weight);
	//       _Nlj[centbin][ptleadbin]+=weight;
	//     }
	//   }
	// }
	// //! SubLeading Jet 
	// if(jetcounter == 1){
	//   int ptsubleadbin = -1;
	//   for (size_t j = 0; j < _nptsubleadbins; ++j) {
	//     if (jet.pt() >= _ptsubleadedges[j] && jet.pt() < _ptsubleadedges[j+1]) {
	//       ptsubleadbin = j;
	//     }	
	//   }
	//   if(ptsubleadbin!=-1){
	//     if(zg.second > 0.0){
	//       // _h_zg_sublead[centbin][ptsubleadbin]->fill(zg.first, weight);
	//       hzg_sublead[centbin][ptsubleadbin]->Fill(zg.first, weight);
	//       hRg_sublead[centbin][ptsubleadbin]->Fill(zg.second, weight);
	//       hMj_sublead[centbin][ptsubleadbin]->Fill(jet.m(), weight);
	//       _Nslj[centbin][ptsubleadbin]+=weight;
	//     }
	//   }
	// }
	// jetcounter++;


	//! print the rec_sdpt vectors to make sure that pT of successive legs are smaller than previous legs
	// if(printDebug){
	//   std::cout<<std::endl<<"size of rec_sdpt = "<<rec_sdpt.size() << "which is basically number of jets "<<std::endl;
	//   if(rec_sdpt.size() >0){
	    
	//     for(unsigned it = 0; it<rec_sdpt.size(); ++it){
	      
	//       std::cout<<"     splits per jet = rec_sdpt.at(it).size() = " << rec_sdpt.at(it).size()<<std::endl;
	      
	//       if(rec_sdpt.at(it).size() > 0){

	// 	for(unsigned is = 0; is < rec_sdpt.at(it).size(); ++is){

	// 	  std::cout<<"             split #"<<is+1<<"  rec pT = "<<rec_sdpt.at(it).at(is)<<std::endl<<std::endl;
		  
	// 	}
		
	//       }

	//     }

	//   }
	// }
      }
      PartonTree->Fill();
      m0PartonTree->Fill();
      ResultTree->Fill();
      
      // jetcounter = 0;
      
      PseudoJets pJet_hardcore;
      //! first add the hard core objects in the event. 	      
      foreach ( const Particle& p, FS) {
      	if(p.pt() >= 2.*GeV){
      	  pJet_hardcore.push_back(PseudoJet(p.px(), p.py(), p.pz(), p.E()));
      	}
      }
      if ( printDebug )
      	std::cout<<"list of hard core objects in the event = "<<pJet_hardcore.size()<<std::endl;
      // //! Hard core Jets 	
      fastjet::ClusterSequence cs_hardcore(pJet_hardcore, jet_def);
      PseudoJets recoJets_HT = sorted_by_pt(cs_hardcore.inclusive_jets(_pTCut));
      PseudoJets recoJets_hardcore = select_both_hi(recoJets_HT);
      
      
      if(recoJets_hardcore.size()==0){
      	// std::cout<<"NO hard core RECONSTRUCTED JETS "<<std::endl;
      	vetoEvent;
      }

      if(printDebug){
      	std::cout<<"Hard core Jets: "<<std::endl;
      	foreach (const PseudoJet jet, recoJets_hardcore){
      	  std::cout<<"    "<<jet.pt()<<std::endl;
      	}
      }
      
      foreach(const PseudoJet jet, recoJets_hardcore){
	//! Inclusive
	hHCJetpT[centbin]->Fill(jet.pt(), weight);
	hHCJetEta[centbin]->Fill(jet.eta(), weight);
	hHCJetPhi[centbin]->Fill(jet.phi(), weight);
	HCjetpT.push_back(jet.pt());
	HCjetM.push_back(jet.m());
	HCjeteta.push_back(jet.eta());
	HCjetphi.push_back(jet.phi());

	fastjet::contrib::RecursiveSymmetryCutBase::SymmetryMeasure  symmetry_measure = fastjet::contrib::RecursiveSymmetryCutBase::scalar_z;
	fastjet::contrib::SoftDrop sd(beta, z_cut, symmetry_measure, _jetR);
	PseudoJet sd_jet = sd(jet);    



	
	vector<PseudoJet> jetconsts = jet.constituents();
	fastjet::JetDefinition jetd(fastjet::antikt_algorithm, _jetR*0.5);
	fastjet::ClusterSequence clust_seq_full(jetconsts, jetd);

	// cout<<"at the const sub jets section"<<endl;

	vector<PseudoJet> subjets= sorted_by_pt(clust_seq_full.inclusive_jets());
	if(subjets.size() >= 3){
	  
	  //! three subjets - perform some kinda grooming
	  // vector<PseudoJet> threesubjets = sorted_by_pt(clust_seq_full.exclusive_jets(3));
	  PseudoJet comblead = subjets.at(0) + subjets.at(1);
	  double zsj = subjets.at(2).pt()/(comblead.pt()+subjets.at(2).pt());
	  double rsj = comblead.delta_R(subjets.at(2));
	  if(zsj > 0.1){
	  HCthreesubjet_z.push_back(zsj);
	  HCthreesubjet_theta.push_back(rsj);
	  }else {
	  HCthreesubjet_z.push_back(-999);
	  HCthreesubjet_theta.push_back(-999);
	  }
	}else{
	HCthreesubjet_z.push_back(-999);
	HCthreesubjet_theta.push_back(-999);
	}

	if(subjets.size() >= 2){
	  //! three subjets - no grooming. just straight up two subjets 
	  double zsj = subjets.at(1).pt()/(subjets.at(0).pt()+subjets.at(1).pt());
	  double rsj = subjets.at(0).delta_R(subjets.at(1));
	  if(zsj > 0.1){
	  HCtwosubjet_z.push_back(zsj);
	  HCtwosubjet_theta.push_back(rsj);
	  }else {
	  HCtwosubjet_z.push_back(-999);
	  HCtwosubjet_theta.push_back(-999);
	  }
	}else{
	HCtwosubjet_z.push_back(-999);
	HCtwosubjet_theta.push_back(-999);
	}


	
	if(sd_jet != 0){

	  HCsdjetpT.push_back(sd_jet.pt());
	  HCsdjetM.push_back(sd_jet.m());
	  double z = sd_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry();
	  double r = sd_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R();
	  
	  HCzg.push_back(z);
	  HCrg.push_back(r);

	  int ptbin = -1;
	  for (size_t j = 0; j < _nptbins; ++j) {
	    if (jet.pt() >= _ptedges[j] && jet.pt() < _ptedges[j+1]) {
	      ptbin = j;
	    }
	  }
	  if(ptbin == -1) continue;
	  const std::pair<double,double> zg = SoftDrop(jet, _jetR, z_cut, beta);
	  if(zg.second > 0.0){
	    // _h_zg[centbin][ptbin]->fill(zg.first, weight);
	    hHCzg[centbin][ptbin]->Fill(zg.first, weight);
	    hHCRg[centbin][ptbin]->Fill(zg.second, weight);
	    hHCMj[centbin][ptbin]->Fill(jet.m(), weight);
	    _Nj_HC[centbin][ptbin]+=weight;
	  }
	}
	else{
	  HCsdjetpT.push_back(-999);
	  HCsdjetM.push_back(-999);
	  HCzg.push_back(-999);
	  HCrg.push_back(-999);	  	  
	}
	// jetcounter++;
      }

      HCResultTree->Fill();



      //! Compare the HC jets - recoJets_hardcore and Inclusive Jets - recoJets_bkgsub
      //! first check the sizes to make sure they exist and then loop over the smaller jets 
      if(recoJets_hardcore.size() > 0 && recoJets_bkgsub.size() > 0) {

	if(recoJets_hardcore.size() <= recoJets_bkgsub.size()) {

	  PseudoJets recoJets_matched;
	  
	  foreach(const PseudoJet jet_hc, recoJets_hardcore){

	    foreach(const PseudoJet jet_inc, recoJets_bkgsub){

	      if(deltaR(jet_inc.eta(), jet_inc.phi_std(), jet_hc.eta(), jet_hc.phi_std()) < _jetR) {

		PseudoJet jet_mc = jet_inc;
		recoJets_matched.push_back(jet_mc);
		
	      }

	    }

	  }

	  if(recoJets_matched.size() >= 1){

	    if(printDebug){
	      std::cout<<"matched Jets: "<<std::endl;
	      foreach (const PseudoJet jet, recoJets_matched){
		std::cout<<"    "<<jet.pt()<<std::endl;
	      }
	    }
	    
	    foreach(const PseudoJet jet, recoJets_matched){
	      //! Inclusive
	      hMCJetpT[centbin]->Fill(jet.pt(), weight);
	      hMCJetEta[centbin]->Fill(jet.eta(), weight);
	      hMCJetPhi[centbin]->Fill(jet.phi(), weight);
	      MCjetpT.push_back(jet.pt());
	      MCjetM.push_back(jet.m());
	      MCjeteta.push_back(jet.eta());
	      MCjetphi.push_back(jet.phi());
	      
	      fastjet::contrib::RecursiveSymmetryCutBase::SymmetryMeasure  symmetry_measure = fastjet::contrib::RecursiveSymmetryCutBase::scalar_z;
	      fastjet::contrib::SoftDrop sd(beta, z_cut, symmetry_measure, _jetR);
	      PseudoJet sd_jet = sd(jet);    


	
	      vector<PseudoJet> jetconsts = jet.constituents();
	      fastjet::JetDefinition jetd(fastjet::antikt_algorithm, _jetR*0.5);
	      fastjet::ClusterSequence clust_seq_full(jetconsts, jetd);

	      // cout<<"at the const sub jets section"<<endl;

	      vector<PseudoJet> subjets= sorted_by_pt(clust_seq_full.inclusive_jets());
	      if(subjets.size() >= 3){
	  
		//! three subjets - perform some kinda grooming
		// vector<PseudoJet> threesubjets = sorted_by_pt(clust_seq_full.exclusive_jets(3));
		PseudoJet comblead = subjets.at(0) + subjets.at(1);
		double zsj = subjets.at(2).pt()/(comblead.pt()+subjets.at(2).pt());
		double rsj = comblead.delta_R(subjets.at(2));
		if(zsj > 0.1){
		  MCthreesubjet_z.push_back(zsj);
		  MCthreesubjet_theta.push_back(rsj);
		}else {
		  MCthreesubjet_z.push_back(-999);
		  MCthreesubjet_theta.push_back(-999);
		}
	      }else{
		MCthreesubjet_z.push_back(-999);
		MCthreesubjet_theta.push_back(-999);
	      }

	      if(subjets.size() >= 2){
		//! three subjets - no grooming. just straight up two subjets 
		double zsj = subjets.at(1).pt()/(subjets.at(0).pt()+subjets.at(1).pt());
		double rsj = subjets.at(0).delta_R(subjets.at(1));
		if(zsj > 0.1){
		  MCtwosubjet_z.push_back(zsj);
		  MCtwosubjet_theta.push_back(rsj);
		}else {
		  MCtwosubjet_z.push_back(-999);
		  MCtwosubjet_theta.push_back(-999);
		}
	      }else{
		MCtwosubjet_z.push_back(-999);
		MCtwosubjet_theta.push_back(-999);
	      }



	      
	      if(sd_jet != 0){
		
		MCsdjetpT.push_back(sd_jet.pt());
		MCsdjetM.push_back(sd_jet.m());
		double z = sd_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry();
		double r = sd_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R();
		
		MCzg.push_back(z);
		MCrg.push_back(r);
		
		int ptbin = -1;
		for (size_t j = 0; j < _nptbins; ++j) {
		  if (jet.pt() >= _ptedges[j] && jet.pt() < _ptedges[j+1]) {
		    ptbin = j;
		  }
		}
		if(ptbin == -1) continue;
		const std::pair<double,double> zg = SoftDrop(jet, _jetR, z_cut, beta);
		if(zg.second > 0.0){
		  // _h_zg[centbin][ptbin]->fill(zg.first, weight);
		  hMCzg[centbin][ptbin]->Fill(zg.first, weight);
		  hMCRg[centbin][ptbin]->Fill(zg.second, weight);
		  hMCMj[centbin][ptbin]->Fill(jet.m(), weight);
		  _Nj_MC[centbin][ptbin]+=weight;
		}
	      }	else{
		MCsdjetpT.push_back(-999);
		MCsdjetM.push_back(-999);
		MCzg.push_back(-999);
		MCrg.push_back(-999);	  	  
	      }
	      
	    }
	  
	  }
	    
	} else if(recoJets_bkgsub.size() < recoJets_hardcore.size()) {

	  std::cout<<"You have way more hardcore jets than inclusive jets! probably shoudnt happen!"<<std::endl;
	  /*
	  foreach (const PseudoJet jet, recoJets_hardcore){
	    std::cout<<"    HC : "<<jet.pt()<<", "<<jet.eta()<<", "<<jet.phi()<<std::endl;
	  }
	  std::cout<<"---------------"<<std::endl;
	  foreach (const PseudoJet jet, recoJets_bkgsub){
	    std::cout<<"    IC : "<<jet.pt()<<", "<<jet.eta()<<", "<<jet.phi()<<std::endl;
	  }
	  
	  std::cout<<"**************"<<std::endl;
	  */
	}

	MCResultTree->Fill();

      }
      
      
      
    }//! analyze function

    /// Normalise histograms etc., after the run
    void finalize() {

      // TH1F * hRatio[4];
      // TH1F * hRatio_lead[3];
      // TH1F * hRatio_sublead[3];

      
      for (size_t i = 0; i < _ncentbins; ++i) {

	// // scale(_h_aj_hardcore[i], 1./(_NevHardcore[i]>0.?_NevHardcore[i]:1.));
	// hAJ_hardcore[i]->Scale(1./(_NevHardcore[i]>0.?_NevHardcore[i]:1.));
	// hAJ_hardcore[i]->Write();
	// // scale(_h_aj_matched[i], 1./(_NevMatched[i]>0.?_NevMatched[i]:1.));
	// hAJ_matched[i]->Scale(1./(_NevMatched[i]>0.?_NevMatched[i]:1.));
	// hAJ_matched[i]->Write();
	
	for (size_t j = 0; j < _nptbins; ++j) {
	  // scale(_h_zg[i][j], 1./(_Nj[i][j]>0.?_Nj[i][j]:1.));
	  hzg[i][j]->Scale(1./(_Nj[i][j]>0.?_Nj[i][j]:1.));
	  hzg[i][j]->Write();
	  hRg[i][j]->Scale(1./(_Nj[i][j]>0.?_Nj[i][j]:1.));
	  hRg[i][j]->Write();
	  hMj[i][j]->Scale(1./(_Nj[i][j]>0.?_Nj[i][j]:1.));
	  hMj[i][j]->Write();
	  if(i >= 1){
	    hRatio[j] = (TH1F*)hzg[i][j]->Clone(Form("Ratio_centbin%d_ptbin%d", i, j));
	    hRatio[j]->Divide(hzg[0][j]);
	    hRatio[j]->Write();
	    hRatio_Rg[j] = (TH1F*)hRg[i][j]->Clone(Form("Ratio_Rg_centbin%d_ptbin%d", i, j));
	    hRatio_Rg[j]->Divide(hRg[0][j]);
	    hRatio_Rg[j]->Write();
	    hRatio_Mj[j] = (TH1F*)hMj[i][j]->Clone(Form("Ratio_Mj_centbin%d_ptbin%d", i, j));
	    hRatio_Mj[j]->Divide(hMj[0][j]);
	    hRatio_Mj[j]->Write();
	    // for(size_t k = 0; k<_zgcenters.size(); ++k){
	    //   _h_ratio[j]->addPoint(_zgcenters[k],
	    // 			    hRatio[j]->GetBinContent(hRatio[j]->FindBin(_zgcenters[k])),
	    // 			    0.025,
	    // 			    hRatio[j]->GetBinError(hRatio[j]->FindBin(_zgcenters[k])));
	    // }	    
	  }	  


	  hHCzg[i][j]->Scale(1./(_Nj_HC[i][j]>0.?_Nj_HC[i][j]:1.));
	  hHCzg[i][j]->Write();
	  hHCRg[i][j]->Scale(1./(_Nj_HC[i][j]>0.?_Nj_HC[i][j]:1.));
	  hHCRg[i][j]->Write();
	  hHCMj[i][j]->Scale(1./(_Nj_HC[i][j]>0.?_Nj_HC[i][j]:1.));
	  hHCMj[i][j]->Write();
	  if(i >= 1){
	    hHCRatio[j] = (TH1F*)hHCzg[i][j]->Clone(Form("HCRatio_centbin%d_ptbin%d", i, j));
	    hHCRatio[j]->Divide(hHCzg[0][j]);
	    hHCRatio[j]->Write();
	    hHCRatio_Rg[j] = (TH1F*)hHCRg[i][j]->Clone(Form("HCRatio_Rg_centbin%d_ptbin%d", i, j));
	    hHCRatio_Rg[j]->Divide(hHCRg[0][j]);
	    hHCRatio_Rg[j]->Write();
	    hHCRatio_Mj[j] = (TH1F*)hHCMj[i][j]->Clone(Form("HCRatio_Mj_centbin%d_ptbin%d", i, j));
	    hHCRatio_Mj[j]->Divide(hHCMj[0][j]);
	    hHCRatio_Mj[j]->Write();
	    // for(size_t k = 0; k<_zgcenters.size(); ++k){
	    //   _h_ratio[j]->addPoint(_zgcenters[k],
	    // 			    hRatio[j]->GetBinContent(hRatio[j]->FindBin(_zgcenters[k])),
	    // 			    0.025,
	    // 			    hRatio[j]->GetBinError(hRatio[j]->FindBin(_zgcenters[k])));
	    // }	    
	  }	  


	  hMCzg[i][j]->Scale(1./(_Nj_MC[i][j]>0.?_Nj_MC[i][j]:1.));
	  hMCzg[i][j]->Write();
	  hMCRg[i][j]->Scale(1./(_Nj_MC[i][j]>0.?_Nj_MC[i][j]:1.));
	  hMCRg[i][j]->Write();
	  hMCMj[i][j]->Scale(1./(_Nj_MC[i][j]>0.?_Nj_MC[i][j]:1.));
	  hMCMj[i][j]->Write();
	  if(i >= 1){
	    hMCRatio[j] = (TH1F*)hMCzg[i][j]->Clone(Form("MCRatio_centbin%d_ptbin%d", i, j));
	    hMCRatio[j]->Divide(hMCzg[0][j]);
	    hMCRatio[j]->Write();
	    hMCRatio_Rg[j] = (TH1F*)hMCRg[i][j]->Clone(Form("MCRatio_Rg_centbin%d_ptbin%d", i, j));
	    hMCRatio_Rg[j]->Divide(hMCRg[0][j]);
	    hMCRatio_Rg[j]->Write();
	    hMCRatio_Mj[j] = (TH1F*)hMCMj[i][j]->Clone(Form("MCRatio_Mj_centbin%d_ptbin%d", i, j));
	    hMCRatio_Mj[j]->Divide(hMCMj[0][j]);
	    hMCRatio_Mj[j]->Write();
	    // for(size_t k = 0; k<_zgcenters.size(); ++k){
	    //   _h_ratio[j]->addPoint(_zgcenters[k],
	    // 			    hRatio[j]->GetBinContent(hRatio[j]->FindBin(_zgcenters[k])),
	    // 			    0.025,
	    // 			    hRatio[j]->GetBinError(hRatio[j]->FindBin(_zgcenters[k])));
	    // }	    
	  }	  

	  
	}


	// for (size_t j = 0; j < _nptleadbins; ++j) {
	//   // scale(_h_zg_lead[i][j], 1./(_Nlj[i][j]>0.?_Nlj[i][j]:1.));
	//   hzg_lead[i][j]->Scale(1./(_Nlj[i][j]>0.?_Nlj[i][j]:1.));
	//   hzg_lead[i][j]->Write();
	//   hRg_lead[i][j]->Scale(1./(_Nlj[i][j]>0.?_Nlj[i][j]:1.));
	//   hRg_lead[i][j]->Write();
	//   hMj_lead[i][j]->Scale(1./(_Nlj[i][j]>0.?_Nlj[i][j]:1.));
	//   hMj_lead[i][j]->Write();
	//   if(i >= 1){
	//     hRatio_lead[j] = (TH1F*)hzg_lead[i][j]->Clone(Form("Ratio_lead_centbin%d_ptbin%d", i, j));
	//     hRatio_lead[j]->Divide(hzg_lead[0][j]);
	//     hRatio_lead[j]->Write();
	//     hRatio_lead_Rg[j] = (TH1F*)hRg_lead[i][j]->Clone(Form("Ratio_Rg_lead_centbin%d_ptbin%d", i, j));
	//     hRatio_lead_Rg[j]->Divide(hRg_lead[0][j]);
	//     hRatio_lead_Rg[j]->Write();
	//     hRatio_lead_Mj[j] = (TH1F*)hMj_lead[i][j]->Clone(Form("Ratio_Mj_lead_centbin%d_ptbin%d", i, j));
	//     hRatio_lead_Mj[j]->Divide(hMj_lead[0][j]);
	//     hRatio_lead_Mj[j]->Write();
	//     // for(size_t k = 0; k<_zgcenters.size(); ++k){
	//     //   _h_ratio_lead[j]->addPoint(_zgcenters[k],
	//     // 				 hRatio_lead[j]->GetBinContent(hRatio_lead[j]->FindBin(_zgcenters[k])),
	//     // 				 0.025,
	//     // 				 hRatio_lead[j]->GetBinError(hRatio_lead[j]->FindBin(_zgcenters[k])));
	//     // }	    
	//   }	  
	// }
	// for (size_t j = 0; j < _nptsubleadbins; ++j) {
	//   // scale(_h_zg_sublead[i][j], 1./(_Nslj[i][j]>0.?_Nslj[i][j]:1.));
	//   hzg_sublead[i][j]->Scale(1./(_Nslj[i][j]>0.?_Nslj[i][j]:1.));
	//   hzg_sublead[i][j]->Write();
	//   hRg_sublead[i][j]->Scale(1./(_Nslj[i][j]>0.?_Nslj[i][j]:1.));
	//   hRg_sublead[i][j]->Write();
	//   hMj_sublead[i][j]->Scale(1./(_Nslj[i][j]>0.?_Nslj[i][j]:1.));
	//   hMj_sublead[i][j]->Write();
	//   if(i >= 1){
	//     hRatio_sublead[j] = (TH1F*)hzg_sublead[i][j]->Clone(Form("Ratio_sublead_centbin%d_ptbin%d", i, j));
	//     hRatio_sublead[j]->Divide(hzg_sublead[0][j]);
	//     hRatio_sublead[j]->Write();
	//     hRatio_sublead_Rg[j] = (TH1F*)hRg_sublead[i][j]->Clone(Form("Ratio_Rg_sublead_centbin%d_ptbin%d", i, j));
	//     hRatio_sublead_Rg[j]->Divide(hRg_sublead[0][j]);
	//     hRatio_sublead_Rg[j]->Write();
	//     hRatio_sublead_Mj[j] = (TH1F*)hMj_sublead[i][j]->Clone(Form("Ratio_Mj_sublead_centbin%d_ptbin%d", i, j));
	//     hRatio_sublead_Mj[j]->Divide(hMj_sublead[0][j]);
	//     hRatio_sublead_Mj[j]->Write();
	//     // for(size_t k = 0; k<_zgcenters.size(); ++k){
	//     //   _h_ratio_sublead[j]->addPoint(_zgcenters[k],
	//     // 				    hRatio_sublead[j]->GetBinContent(hRatio_sublead[j]->FindBin(_zgcenters[k])),
	//     // 				    0.025,
	//     // 				    hRatio_sublead[j]->GetBinError(hRatio_sublead[j]->FindBin(_zgcenters[k])));
	//     // }	    
	//   }	  
	// }
      }//! centrality loop
      // ResultTree->Write();
      fout->Write();
      fout->Close();      
      
    }//! finalize

  protected:

    size_t _mode;
    
    vector<double> _centedges, _ptedges, _ptleadedges, _ptsubleadedges;
    size_t _ncentbins, _nptbins, _nptleadbins, _nptsubleadbins;
    double _jetR;
    int RADIUS;
    double _pTCut;
    double _pTMax;
    double _pTconsMax;
    double _max_track_rap;
    bool printDebug;
    vector<vector<double> > _Nj;
    vector<vector<double> > _Nj_HC;
    vector<vector<double> > _Nj_MC;
    // vector<vector<double> > _Nlj;
    // vector<vector<double> > _Nslj;
    // vector<double> _NevHardcore;
    // vector<double> _NevMatched;
    // vector<double> _zgcenters;
    vector<double> _etabins;
    vector<double> _phibins;
    int _Nbounds_eta;
    int _Nbounds_phi;
    double _delRMin;
    double _etaMax;
    double _phiMax;    
    double z_cut;
    double beta;

    double cent;
    double mcweight;
    double pthat;
    double eventID;
    
    vector<double> trigpT;
    vector<double> recpT;

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
    
    vector<double> m0PLpt;
    vector<double> m0PLm;
    vector<double> m0PLptg;
    vector<double> m0PLmg;
    vector<double> m0PLeta;
    vector<double> m0PLphi;
    vector<double> m0PLzg;
    vector<double> m0PLrg;
    vector<double> m0PLncons;
    vector<double> m0PLnconsg;
    
    vector<vector<double> > conspT;
    vector<vector<double> > consDist;
    vector<vector<double> > consGirth;
    vector<vector<double> > consM;

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
    vector<double> twosubjet_z;
    vector<double> twosubjet_theta;
    vector<double> threesubjet_z;
    vector<double> threesubjet_theta;
    vector<vector<double> > rec_zg;
    vector<vector<double> > rec_rg;
    vector<vector<double> > rec_pt;
    vector<vector<double> > rec_mj;
    vector<vector<double> > rec_sdpt;
    vector<vector<double> > rec_sdmj;
    vector<vector<double> > rec_hard_mass;
    vector<vector<double> > rec_soft_mass;

    vector<double> HCjetpT;
    vector<double> HCjetM;
    vector<double> HCsdjetpT;
    vector<double> HCsdjetM;
    vector<double> HCjeteta;
    vector<double> HCjetphi;
    vector<double> HCzg;
    vector<double> HCrg;
    vector<double> HCtwosubjet_z;
    vector<double> HCtwosubjet_theta;
    vector<double> HCthreesubjet_z;
    vector<double> HCthreesubjet_theta;

    vector<double> MCjetpT;
    vector<double> MCjetM;
    vector<double> MCsdjetpT;
    vector<double> MCsdjetM;
    vector<double> MCjeteta;
    vector<double> MCjetphi;
    vector<double> MCzg;
    vector<double> MCrg;
    vector<double> MCtwosubjet_z;
    vector<double> MCtwosubjet_theta;
    vector<double> MCthreesubjet_z;
    vector<double> MCthreesubjet_theta;
    
    
  private:

    TFile * fout;    
    // Histo1DPtr _h_zg[2][4];
    // Scatter2DPtr _h_ratio[4];
    TH1F * hJetpT[2];
    TH1F * hJetEta[2];
    TH1F * hJetPhi[2];
    TH1F * hzg[2][7];    
    TH1F * hRatio[7];
    TH1F * hRg[2][7];
    TH1F * hRatio_Rg[7];
    TH1F * hMj[2][7];
    TH1F * hRatio_Mj[7];

    TH1F * hHCJetpT[2];
    TH1F * hHCJetEta[2];
    TH1F * hHCJetPhi[2];
    TH1F * hHCzg[2][7];    
    TH1F * hHCRatio[7];
    TH1F * hHCRg[2][7];
    TH1F * hHCRatio_Rg[7];
    TH1F * hHCMj[2][7];
    TH1F * hHCRatio_Mj[7];

    TH1F * hMCJetpT[2];
    TH1F * hMCJetEta[2];
    TH1F * hMCJetPhi[2];
    TH1F * hMCzg[2][7];    
    TH1F * hMCRatio[7];
    TH1F * hMCRg[2][7];
    TH1F * hMCRatio_Rg[7];
    TH1F * hMCMj[2][7];
    TH1F * hMCRatio_Mj[7];

    
    // Histo1DPtr _h_zg_lead[2][3];
    // Scatter2DPtr _h_ratio_lead[3];
    // TH1F * hzg_lead[2][4];
    // TH1F * hRatio_lead[4];
    // TH1F * hRg_lead[2][4];
    // TH1F * hRatio_lead_Rg[4];
    // TH1F * hMj_lead[2][4];
    // TH1F * hRatio_lead_Mj[4];
    // // Histo1DPtr _h_zg_sublead[2][3];
    // // Scatter2DPtr _h_ratio_sublead[3];
    // TH1F * hzg_sublead[2][4];
    // TH1F * hRatio_sublead[4];
    // TH1F * hRg_sublead[2][4];
    // TH1F * hRatio_sublead_Rg[4];
    // TH1F * hMj_sublead[2][4];
    // TH1F * hRatio_sublead_Mj[4];
    // // Histo1DPtr _h_aj_hardcore[2];
    // // Histo1DPtr _h_aj_matched[2];
    // TH1F * hAJ_hardcore[2];
    // TH1F * hAJ_matched[2];

    TTree * TrigRecTree;
    TTree * PartonTree;
    TTree * m0PartonTree;
    TTree * ResultTree;
    TTree * HCResultTree;
    TTree * MCResultTree;
  
    TDatabasePDG * pdg;
    
  };
  
  class STAR_SPLITTING_HER_10PTHAT15 : public STAR_SPLITTING_HER {
  public:
    STAR_SPLITTING_HER_10PTHAT15()
      : STAR_SPLITTING_HER("STAR_SPLITTING_HER_10PTHAT15")
    {
      _mode = 1;
    }
  };
  class STAR_SPLITTING_HER_15PTHAT20 : public STAR_SPLITTING_HER {
  public:
    STAR_SPLITTING_HER_15PTHAT20()
      : STAR_SPLITTING_HER("STAR_SPLITTING_HER_15PTHAT20")
    {
      _mode = 2;
    }
  };
  class STAR_SPLITTING_HER_20PTHAT25 : public STAR_SPLITTING_HER {
  public:
    STAR_SPLITTING_HER_20PTHAT25()
      : STAR_SPLITTING_HER("STAR_SPLITTING_HER_20PTHAT25")
    {
      _mode = 3;
    }
  };
  class STAR_SPLITTING_HER_25PTHAT30 : public STAR_SPLITTING_HER {
  public:
    STAR_SPLITTING_HER_25PTHAT30()
      : STAR_SPLITTING_HER("STAR_SPLITTING_HER_25PTHAT30")
    {
      _mode = 4;
    }
  };
  class STAR_SPLITTING_HER_30PTHAT35 : public STAR_SPLITTING_HER {
  public:
    STAR_SPLITTING_HER_30PTHAT35()
      : STAR_SPLITTING_HER("STAR_SPLITTING_HER_30PTHAT35")
    {
      _mode = 5;
    }
  };
  class STAR_SPLITTING_HER_35PTHAT40 : public STAR_SPLITTING_HER {
  public:
    STAR_SPLITTING_HER_35PTHAT40()
      : STAR_SPLITTING_HER("STAR_SPLITTING_HER_35PTHAT40")
    {
      _mode = 6;
    }
  };
  class STAR_SPLITTING_HER_40PTHAT45 : public STAR_SPLITTING_HER {
  public:
    STAR_SPLITTING_HER_40PTHAT45()
      : STAR_SPLITTING_HER("STAR_SPLITTING_HER_40PTHAT45")
    {
      _mode = 7;
    }
  };
  class STAR_SPLITTING_HER_45PTHAT50 : public STAR_SPLITTING_HER {
  public:
    STAR_SPLITTING_HER_45PTHAT50()
      : STAR_SPLITTING_HER("STAR_SPLITTING_HER_45PTHAT50")
    {
      _mode = 8;
    }
  };
  class STAR_SPLITTING_HER_50PTHAT60 : public STAR_SPLITTING_HER {
  public:
    STAR_SPLITTING_HER_50PTHAT60()
      : STAR_SPLITTING_HER("STAR_SPLITTING_HER_50PTHAT60")
    {
      _mode = 9;
    }
  };
  class STAR_SPLITTING_HER_60PTHAT80 : public STAR_SPLITTING_HER {
  public:
    STAR_SPLITTING_HER_60PTHAT80()
      : STAR_SPLITTING_HER("STAR_SPLITTING_HER_60PTHAT80")
    {
      _mode = 10;
    }
  };


  
  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(STAR_SPLITTING_HER);
  DECLARE_RIVET_PLUGIN(STAR_SPLITTING_HER_10PTHAT15);
  DECLARE_RIVET_PLUGIN(STAR_SPLITTING_HER_15PTHAT20);
  DECLARE_RIVET_PLUGIN(STAR_SPLITTING_HER_20PTHAT25);
  DECLARE_RIVET_PLUGIN(STAR_SPLITTING_HER_25PTHAT30);
  DECLARE_RIVET_PLUGIN(STAR_SPLITTING_HER_30PTHAT35);
  DECLARE_RIVET_PLUGIN(STAR_SPLITTING_HER_35PTHAT40);
  DECLARE_RIVET_PLUGIN(STAR_SPLITTING_HER_40PTHAT45);
  DECLARE_RIVET_PLUGIN(STAR_SPLITTING_HER_45PTHAT50);
  DECLARE_RIVET_PLUGIN(STAR_SPLITTING_HER_50PTHAT60);
  DECLARE_RIVET_PLUGIN(STAR_SPLITTING_HER_60PTHAT80);
  				       
  }
