//Isaac Mooney, WSU, September 2019
//This file implements the MC closure test
//lately I've been using a more up-to-date version of this macro on my desktop, and this one is now deprecated (not irreparably)
//if I have time later, I'll update it with the new version

#include <ctime>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <algorithm>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TRandom.h>

#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

using namespace std;

int main (int argc, const char** argv) {
  //intro                                                                                                                    
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                                                        
  //basic argument checking.                                                                                                                                   
  if (argc < 2) {
    cerr << "Should be at least one argument: jet radius (can also pass a grooming flag, after). Received "
         << argc-1 << ". Exiting." << endl;
    exit(1);
  }

  //argv[1] should be the jet radius e.g. "04".                                                                                                                
  string radius = (string) argv[1];
  radius = "_R"+radius;//appending the "_R" used in file naming.                                                                                               
  bool groomflag = 0;
  string hname = "m"; //for internal histogram name calls: default                                                                                             
  string htitle = "M"; //for axis titles of histograms: default                                                                                                
  if (argc > 2) {
    string groom_str = (string) argv[2];
    if (groom_str == "g") {
      groomflag = 1; //a logic switch                                                                                                                          
      hname = "mg"; //for internal histogram name calls                                                                                                        
      htitle = "M_{g}"; //for axis titles of histograms                                                                                                        
      cout << "Running Mg unfolding on groomed jets!" << endl;
    }
  }

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();

  const string match_path = "~/jetmass2/out/sim/";
  const string match_file = "sim_matched";

  //input files                                                                                                                                                
  TFile *f = new TFile((match_path+match_file+radius+"_paper_bindropped_new.root").c_str(),"READ");
  //responses
  RooUnfoldResponse *r1D_m_A = (RooUnfoldResponse*) f->Get(("sampleA_"+hname+"_response").c_str());
  RooUnfoldResponse *r1D_pt_A = (RooUnfoldResponse*) f->Get("sampleA_pt_response");
  RooUnfoldResponse *r2D_A = (RooUnfoldResponse*) f->Get(("sampleA_"+hname+"_pt_response").c_str());
  RooUnfoldResponse *r1D_m_B = (RooUnfoldResponse*) f->Get(("sampleB_"+hname+"_response").c_str());
  RooUnfoldResponse *r1D_pt_B = (RooUnfoldResponse*) f->Get("sampleB_pt_response");
  RooUnfoldResponse *r2D_B = (RooUnfoldResponse*) f->Get(("sampleB_"+hname+"_pt_response").c_str());
  
  //pseudo-data                           
  //pop A
  TH2D* m_pt_det_A = (TH2D*) f->Get(("sampleA_"+hname+"_pt_det").c_str());
  TH2D* m_pt_gen_A = (TH2D*) f->Get(("sampleA_"+hname+"_pt_gen").c_str());
  TH1D* m_det_A = (TH1D*) f->Get(("sampleA_"+hname+"_det").c_str());
  TH1D* m_gen_A = (TH1D*) f->Get(("sampleA_"+hname+"_gen").c_str());
  TH1D* pt_det_A = (TH1D*) f->Get("sampleA_pt_det");
  TH1D* pt_gen_A = (TH1D*) f->Get("sampleA_pt_gen");
  //pop B
  TH2D* m_pt_det_B = (TH2D*) f->Get(("sampleB_"+hname+"_pt_det").c_str());
  TH2D* m_pt_gen_B = (TH2D*) f->Get(("sampleB_"+hname+"_pt_gen").c_str());
  TH1D* m_det_B = (TH1D*) f->Get(("sampleB_"+hname+"_det").c_str());
  TH1D* m_gen_B = (TH1D*) f->Get(("sampleB_"+hname+"_gen").c_str());
  TH1D* pt_det_B = (TH1D*) f->Get("sampleB_pt_det");
  TH1D* pt_gen_B = (TH1D*) f->Get("sampleB_pt_gen");
  
  cout << "B" << endl;
  RooUnfoldBayes *unfold_opp_m1D = new RooUnfoldBayes(r1D_m_A, m_det_B, 4, false, "unfold_opp_m1D","");
  RooUnfoldBayes *unfold_opp_pt1D = new RooUnfoldBayes(r1D_pt_A, pt_det_B, 4, false, "unfold_opp_pt1D","");
  RooUnfoldBayes *unfold_opp_2D = new RooUnfoldBayes(r2D_A, m_pt_det_B, 4, false, "unfold_opp_2D","");
  RooUnfoldBayes *unfold_same_m1D = new RooUnfoldBayes(r1D_m_A, m_det_A, 4, false, "unfold_same_m1D","");
  RooUnfoldBayes *unfold_same_pt1D = new RooUnfoldBayes(r1D_pt_A, pt_det_A, 4, false, "unfold_same_pt1D","");
  RooUnfoldBayes *unfold_same_2D = new RooUnfoldBayes(r2D_A, m_pt_det_A, 4, false, "unfold_same_2D","");
  
  cout << "C" << endl;
  
  TH1D *reco_opp_m1D = (TH1D*) unfold_opp_m1D->Hreco((RooUnfold::ErrorTreatment) 3);
  TH1D *reco_opp_pt1D = (TH1D*) unfold_opp_pt1D->Hreco((RooUnfold::ErrorTreatment) 3);
  TH2D *reco_opp_2D = (TH2D*) unfold_opp_2D->Hreco((RooUnfold::ErrorTreatment) 3);
  TH1D *reco_same_m1D = (TH1D*) unfold_same_m1D->Hreco((RooUnfold::ErrorTreatment) 3);
  TH1D *reco_same_pt1D = (TH1D*) unfold_same_pt1D->Hreco((RooUnfold::ErrorTreatment) 3);
  TH2D *reco_same_2D = (TH2D*) unfold_same_2D->Hreco((RooUnfold::ErrorTreatment) 3);

  reco_opp_m1D->Divide(m_gen_A);
  reco_opp_pt1D->Divide(pt_gen_A);
  TH1D* opp_pt2D = (TH1D*) reco_opp_2D->ProjectionY("opp_pt2D");
  vector<TH1D*> opp_m_projs = {(TH1D*) reco_opp_2D->ProjectionX("opp_m2030",reco_opp_2D->GetYaxis()->FindBin(20),reco_opp_2D->GetYaxis()->FindBin(30)-1),
  (TH1D*) reco_opp_2D->ProjectionX("opp_m3045",reco_opp_2D->GetYaxis()->FindBin(30),reco_opp_2D->GetYaxis()->FindBin(45)-1)};
  reco_same_m1D->Divide(m_gen_A);
  reco_same_pt1D->Divide(pt_gen_A);
  TH1D* same_pt2D = (TH1D*) reco_same_2D->ProjectionY("same_pt2D");
  vector<TH1D*> same_m_projs = {(TH1D*) reco_same_2D->ProjectionX("same_m2030",reco_same_2D->GetYaxis()->FindBin(20),reco_same_2D->GetYaxis()->FindBin(30)-1),
  (TH1D*) reco_same_2D->ProjectionX("same_m3045",reco_same_2D->GetYaxis()->FindBin(30),reco_same_2D->GetYaxis()->FindBin(45)-1)};
  vector<TH1D*> gen_m_A_projs = {(TH1D*) m_pt_gen_A->ProjectionX("gen_m2030",m_pt_gen_A->GetYaxis()->FindBin(20),m_pt_gen_A->GetYaxis()->FindBin(30)-1),
  (TH1D*) m_pt_gen_A->ProjectionX("gen_m3045",m_pt_gen_A->GetYaxis()->FindBin(30),m_pt_gen_A->GetYaxis()->FindBin(45)-1)};
  
  const int nBins = 3; //20-25, 25-30,30-40
  for (int i = 0; i < nBins; ++ i) {
    opp_m_projs[i]->Divide(gen_m_A_projs[i]); //det B unfolded w/ A, divided by gen A
    same_m_projs[i]->Divide(gen_m_A_projs[i]);//det A unfolded w/ A, divided by gen A
  }
  


  string ftitle = "closure";
  if (groomflag) {
    ftitle = "groomed_"+ftitle;
  }

  TFile *fout = new TFile(("~/jetmass2/out/closure/"+ftitle+radius+"_new.root").c_str(),"RECREATE");
  fout->cd();
  
  reco_opp_m1D->Write();
  reco_opp_pt1D->Write();
  reco_same_m1D->Write();
  reco_same_pt1D->Write();
  for (int i = 0; i < nBins; ++ i) {
    opp_m_projs[i]->Write();
    same_m_projs[i]->Write();
  }
  
  fout->Close();
  
  return 0;
  
}
