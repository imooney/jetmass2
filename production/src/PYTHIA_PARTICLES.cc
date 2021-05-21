//! -*- C++ -*-

//! Raghav Kunnawalkam Elayavalli
//! Wayne State University
//! Wednesday October 8th, 2018 
//!
//! Analysis to read in pythia particles which will be used to generate a ttree of particles 
//!
//!
//! herwig pthat bin - 3-4, 4-5, 5-7, 7-9, 9-11, 11-15, 15-25, 25-35, 35-45, 45-55, 55-65  

#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "HepMC/GenParticle.h"
#include "HepMC/GenEvent.h"
//!ROOT headers
#include "TTree.h"
#include "TFile.h"

namespace Rivet {

  class PYTHIA_PARTICLES : public Analysis {
  public:

    /// Constructor
    PYTHIA_PARTICLES(string name = "PYTHIA_PARTICLES")
      : Analysis(name)
    {
      setNeedsCrossSection(true);
    }

       
    /// Book histograms and initialise projections before the run
    void init() {

      //! final state and jet projection
      //! generic final state
      FinalState fs(-1.0, 1.0, 0.*GeV);
      addProjection(fs, "FS");

      fout = new TFile("pthia8_particleTree.root", "RECREATE");

      ParticleTree = new TTree("ParticleTree", "final state particle tree");
      partpT.clear();
      parteta.clear();
      partphi.clear();
      partm.clear();
      ParticleTree->Branch("mcweight", &mcweight,"mcweight/d");
      ParticleTree->Branch("pthat", &pthat,"pthat/d");
      ParticleTree->Branch("partpT", &partpT);
      ParticleTree->Branch("parteta", &parteta);
      ParticleTree->Branch("partphi", &partphi);
      ParticleTree->Branch("partm", &partm);
      ParticleTree->Branch("partc", &partc);
      
    }//! init

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      mcweight = handler().crossSection();
      // mcweight = event.weight();
      
      //! clear the vectors for each event! 
      partpT.clear();
      parteta.clear();
      partphi.clear();
      partm.clear();
      partc.clear();
      
      pthat = event.genEvent()->event_scale();
      
      const ParticleVector& FS = applyProjection<FinalState>(event, "FS").particlesByPt();     
      foreach ( const Particle& p, FS) {
	PseudoJet part = PseudoJet(p.px(), p.py(), p.pz(), p.E());
	partpT.push_back(part.pt());
	parteta.push_back(part.eta());
	partphi.push_back(part.phi());
	partm.push_back(part.m());
	partc.push_back(p.charge());
      }

      ParticleTree->Fill();
      
    }//! analyze function

    /// Normalise histograms etc., after the run
    void finalize() {

      fout->Write();
      fout->Close();      
      
    }//! finalize

  protected:

    double mcweight;
    double pthat;
    vector<double> partpT;    
    vector<double> parteta;    
    vector<double> partphi;    
    vector<double> partm;    
    vector<double> partc;    
    
  private:

    TFile * fout;    

    TTree * ParticleTree;
    
  };

  
  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(PYTHIA_PARTICLES);


}
