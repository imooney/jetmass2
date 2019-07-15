//Isaac Mooney, WSU, June 2019
//This file takes in a root file, pulls the tree(s) out, fills histograms with the entries, and writes to an output root file

#include <ctime>
#include <iostream>
#include <iomanip>
#include <math.h>

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

using namespace std;

//self explanatory: takes branches in a tree, "treestr", in file f, and fills histograms with them
void TreetoHist (TFile *f, string treestr, vector<TH1D*> h1Ds, vector<TH2D*> h2Ds, vector<TH3D*> h3Ds, bool dataflag) {  
  //initializing the variables that will be filled by the values in the branches later
  vector<double> *Pt = 0; vector<double> *Eta = 0; vector<double> *Phi = 0;
  vector<double> *M = 0; vector<double> *zg = 0; vector<double> *rg = 0; vector<double> *mg = 0;
  vector<double> *ptg = 0; vector<double> *E = 0; vector<double> *ch_e_frac = 0; vector<double> *mcd = 0;
  vector<double> *tau0 = 0; vector<double> *tau05 = 0; vector<double> *tau_05 = 0; vector<double> *tau_1 = 0;
  vector<double> *tau0_g = 0; vector<double> *tau05_g = 0; vector<double> *tau_05_g = 0; vector<double> *tau_1_g = 0;

  double weight = 1, n_jets = -1, bbc_east_rate = -1, bbc_east_sum = -1;
    
  //getting the tree and linking the branches with the variables
  TTree *t = (TTree*) f->Get(treestr.c_str());
  t->SetBranchAddress("n_jets",&n_jets);
  t->SetBranchAddress("bbc_east_rate",&bbc_east_rate);
  t->SetBranchAddress("bbc_east_sum",&bbc_east_sum);
  t->SetBranchAddress("Pt",&Pt);
  t->SetBranchAddress("Eta",&Eta);
  t->SetBranchAddress("Phi",&Phi);
  t->SetBranchAddress("E",&E);
  t->SetBranchAddress("M",&M);
  t->SetBranchAddress("zg",&zg);
  t->SetBranchAddress("rg",&rg);
  t->SetBranchAddress("mg",&mg);
  t->SetBranchAddress("ptg",&ptg);
  t->SetBranchAddress("ch_e_frac",&ch_e_frac);
  t->SetBranchAddress("mcd",&mcd);
  t->SetBranchAddress("tau0",&tau0);
  t->SetBranchAddress("tau05",&tau05);
  t->SetBranchAddress("tau_05",&tau_05);
  t->SetBranchAddress("tau_1",&tau_1);
  t->SetBranchAddress("tau0_g",&tau0_g); 
  t->SetBranchAddress("tau05_g",&tau05_g); 
  t->SetBranchAddress("tau_05_g",&tau_05_g);
  t->SetBranchAddress("tau_1_g",&tau_1_g);
  
  if (dataflag == 0) {t->SetBranchAddress("weight", &weight);}

  cout << ("RUNNING OVER TREE "+treestr+"! Entries: ").c_str() << t->GetEntries() << endl;
  const clock_t begin_time = clock(); //timing - for debugging and for fun
  for (unsigned i = 0; i < t->GetEntries(); ++ i) { //"event" loop
    if (i % 1000 == 0 && i != 0) { //can change this to a frequency of your preference (for real data I use 1e5 or 1e6)
      cout << "Still chuggin. On event " << i << endl;
      cout << "Total time passed: " << fixed << setprecision(5) << double(clock() - begin_time) /(double) CLOCKS_PER_SEC << " secs" << endl;
    }
    t->GetEntry(i);
    
    //filling "event" observables
    h1Ds[13]->Fill(n_jets, weight);
    //looping over "tracks" and filling histograms
    for (unsigned j = 0; j < Pt->size(); ++ j) { //all vectors of doubles in the branches should have the same size 
      //inclusive plots (1D)                                                                                                                                   
      h1Ds[0]->Fill(Pt->at(j), weight);
      h1Ds[1]->Fill(Pt->at(j), weight); //different binnings between these two ^
      h1Ds[2]->Fill(Eta->at(j), weight);
      h1Ds[3]->Fill(Phi->at(j), weight);
      h1Ds[4]->Fill(E->at(j), weight);
      h1Ds[5]->Fill(M->at(j), weight);
      h1Ds[6]->Fill(ch_e_frac->at(j), weight);
      if (zg->at(j) >= 0.1) { //otherwise, we tag and drop the groomed jet
        h1Ds[7]->Fill(zg->at(j), weight);
        h1Ds[8]->Fill(rg->at(j), weight);
        h1Ds[9]->Fill(mg->at(j), weight);
        h1Ds[10]->Fill(ptg->at(j), weight);
        h1Ds[11]->Fill(ptg->at(j) / (double) Pt->at(j));
      }
      h1Ds[12]->Fill(mcd->at(j), weight); //collinear dropped mass is filled whether or not the jet passes the SD criterion. Can be changed if necessary.
      
      //2Ds!                                                                                                                                           
      h2Ds[0]->Fill(M->at(j), Pt->at(j), weight);
      h2Ds[1]->Fill(ch_e_frac->at(j), Pt->at(j), weight);
      if (zg->at(j) >= 0.1) {                                                                                                                           
	h2Ds[2]->Fill(zg->at(j), Pt->at(j), weight);
	h2Ds[3]->Fill(rg->at(j), Pt->at(j), weight);
	h2Ds[4]->Fill(mg->at(j), Pt->at(j), weight);
	h2Ds[5]->Fill(ptg->at(j), Pt->at(j), weight);
	h2Ds[6]->Fill(ptg->at(j) / (double) Pt->at(j), Pt->at(j), weight);
	h2Ds[11]->Fill(mg->at(j), Pt->at(j), 1); //to use later to drop stats-limited bins                                                                   
	h2Ds[16]->Fill(log10(tau0_g->at(j)), Pt->at(j), weight);
	h2Ds[17]->Fill(log10(tau05_g->at(j)), Pt->at(j), weight);
	h2Ds[18]->Fill(log10(tau_05_g->at(j)), Pt->at(j), weight);
	h2Ds[19]->Fill(log10(tau_1_g->at(j)), Pt->at(j), weight);
      }                                                                                                                                                     
      h2Ds[7]->Fill(mcd->at(j), Pt->at(j), weight);
      h2Ds[8]->Fill(Phi->at(j), Pt->at(j), weight);
      h2Ds[9]->Fill(Eta->at(j), Pt->at(j), weight);
      h2Ds[10]->Fill(M->at(j), Pt->at(j), 1); //to use later to drop stats-limited bins                                                                      
      h2Ds[12]->Fill(log10(tau0->at(j)), Pt->at(j), weight);
      h2Ds[13]->Fill(log10(tau05->at(j)), Pt->at(j), weight);
      h2Ds[14]->Fill(log10(tau_05->at(j)), Pt->at(j), weight);
      h2Ds[15]->Fill(log10(tau_1->at(j)), Pt->at(j), weight);
      
      h3Ds[0]->Fill(bbc_east_rate, Pt->at(j), M->at(j), weight);
      h3Ds[1]->Fill(bbc_east_sum, Pt->at(j), M->at(j), weight);
    }//!jet loop
  }//!event loop
  
  //!needs to be outside the event loop; not sure exactly what it does
  t->ResetBranchAddresses();
  return;
}

int main (int argc, const char ** argv) {
  //intro
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  //basic argument checking.
  if (argc != 4) {
    cerr << "Should be three arguments: output location, output name, input name. Received "
	 << argc << ". Exiting." << endl;
    exit(1);
  }

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();

  //opening file containing some example trees
  //argv[3] should be the name of the input file
  string fin_name = (string) argv[3];
  TFile *fin = new TFile(fin_name.c_str(),"READ");
  cout << "DEBUG: input file name is " << fin->GetName() << endl;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~hists~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

  //event observables:
  TH1D* n_jets = new TH1D("n_jets","",20,-0.5,19.5);
  
  //pT range for particle-level is different from detector-level so define it here for ease of change later if necessary:
  const int ptbins = 9;
  const double ptlow = 15; const double pthigh = 60;
  
  //making variable bin size histogram                                                                                                                        
  const int nBins_pt = 8;
  double edges[nBins_pt + 1] = {5,10,15,20,25,30,40,60,100};
  
  //1D inclusives (including below 15 GeV pT jets)                                                                                                             
  TH1D* pt_var_bin = new TH1D("pt_var_bin","",nBins_pt,edges);
  TH1D* pt = new TH1D("pt","",15,5,80);
  TH1D* eta = new TH1D("eta","",50,-1,1);
  TH1D* phi = new TH1D("phi","",50,0,2*M_PI);
  TH1D* e = new TH1D("e","",100,0,100);
  TH1D* m = new TH1D("m","",14,0,14);
  TH1D* ch_e_frac = new TH1D("ch_e_frac","",10,0,1);
  TH1D* zg = new TH1D("zg","",20,0,1);
  TH1D* rg = new TH1D("rg","",20,0,1);
  TH1D* mg = new TH1D("mg","",14,0,14);
  TH1D* ptg = new TH1D("ptg","",15,5,80);
  TH1D* ratio_ptg_pt = new TH1D("ratio_ptg_pt","",21,0,2);
  TH1D* mcd = new TH1D("mcd","",20,0,5);

  //1D substructure
  TH1D* conspT = new TH1D("conspT","",30,0.2,30.2);
  TH1D* consDist = new TH1D("consDist","",20,0,0.4);
  TH1D* consGirth = new TH1D("consGirth","",20,0,0.12);
  TH1D* jetMult = new TH1D("jetMult","",25,0,25);
  TH1D* jetGirth = new TH1D("jetGirth","",20,0,0.3);

  //2D hists mostly for taking slices in pT                                 
  TH2D* conspT_v_pt = new TH2D("conspT_v_pt","",30,0.2,30.2,ptbins,ptlow,pthigh);
  TH2D* consDist_v_pt = new TH2D("consDist_v_pt","",20,0,0.4,ptbins,ptlow,pthigh);
  TH2D* consGirth_v_pt = new TH2D("consGirth_v_pt","",20,0,0.12,ptbins,ptlow,pthigh);
  TH2D* jetMult_v_pt = new TH2D("jetMult_v_pt","",25,0,25,ptbins,ptlow,pthigh);
  TH2D* jetGirth_v_pt = new TH2D("jetGirth_v_pt","",20,0,0.3,ptbins,ptlow,pthigh);

  TH2D* phi_v_pt = new TH2D("phi_v_pt","",24,0,2*M_PI,ptbins,ptlow,pthigh);
  TH2D* eta_v_pt = new TH2D("eta_v_pt","",25,-1,1,ptbins,ptlow,pthigh);
  TH2D* m_v_pt = new TH2D("m_v_pt","",14,0,14,ptbins,ptlow,pthigh);
  TH2D* ch_frac_v_pt = new TH2D("ch_frac_v_pt","",10,0,1,ptbins,ptlow,pthigh);
  TH2D* zg_v_pt = new TH2D("zg_v_pt","",20,0,1,ptbins,ptlow,pthigh);
  TH2D* rg_v_pt = new TH2D("rg_v_pt","",20,0,1,ptbins,ptlow,pthigh);
  TH2D* mg_v_pt = new TH2D("mg_v_pt","",14,0,14,ptbins,ptlow,pthigh);
  TH2D* ptg_v_pt = new TH2D("ptg_v_pt","",ptbins,ptlow,pthigh,ptbins,ptlow,pthigh);
  TH2D* ratio_ptgpt_v_pt = new TH2D("ratio_ptgpt_v_pt","",21,0,2,ptbins,ptlow,pthigh);
  TH2D* mcd_v_pt = new TH2D("mcd_v_pt","",20,0,5,ptbins,ptlow,pthigh);

  //jet mass dependence on event activity
  TH3D* bbc_east_rate_v_pt_v_m = new TH3D("bbc_east_rate_v_pt_v_m","",70,0,7e6,ptbins,ptlow,pthigh,14,0,14); //luminosity dependence
  TH3D* bbc_east_sum_v_pt_v_m = new TH3D("bbc_east_sum_v_pt_v_m","",100,0,6.5e4,ptbins,ptlow,pthigh,14,0,14); //centrality dependence
  
  //angularities
  TH2D* tau0_v_pt = new TH2D("tau0_v_pt","",30,-8,0,ptbins,ptlow,pthigh);
  TH2D* tau05_v_pt = new TH2D("tau05_v_pt","",30,-8,0,ptbins,ptlow,pthigh);
  TH2D* tau_05_v_pt = new TH2D("tau_05_v_pt","",30,-8,0,ptbins,ptlow,pthigh);
  TH2D* tau_1_v_pt = new TH2D("tau_1_v_pt","",30,-8,0,ptbins,ptlow,pthigh);
  TH2D* tau0_g_v_pt = new TH2D("tau0_g_v_pt","",30,-8,0,ptbins,ptlow,pthigh);
  TH2D* tau05_g_v_pt = new TH2D("tau05_g_v_pt","",30,-8,0,ptbins,ptlow,pthigh);
  TH2D* tau_05_g_v_pt = new TH2D("tau_05_g_v_pt","",30,-8,0,ptbins,ptlow,pthigh);
  TH2D* tau_1_g_v_pt = new TH2D("tau_1_g_v_pt","",30,-8,0,ptbins,ptlow,pthigh);

  TH2D* m_v_pt_counts = new TH2D("m_v_pt_counts","",14,0,14,9,15,60);
  TH2D* mg_v_pt_counts = new TH2D("mg_v_pt_counts","",14,0,14,9,15,60);

  //putting them in a vector to more easily shuttle them back and forth in the function. Drawback: have to know their order.
  vector<TH1D*> h1Ds = {pt_var_bin,pt,eta,phi,e,m,ch_e_frac,zg,rg,mg,ptg,ratio_ptg_pt,mcd,n_jets};
  vector<TH2D*> h2Ds = {m_v_pt,ch_frac_v_pt,zg_v_pt,rg_v_pt,mg_v_pt,ptg_v_pt,ratio_ptgpt_v_pt,mcd_v_pt,phi_v_pt,eta_v_pt,m_v_pt_counts,mg_v_pt_counts,tau0_v_pt,tau05_v_pt,tau_05_v_pt,tau_1_v_pt,tau0_g_v_pt,tau05_g_v_pt,tau_05_g_v_pt,tau_1_g_v_pt};
  vector<TH3D*> h3Ds = {bbc_east_rate_v_pt_v_m, bbc_east_sum_v_pt_v_m};
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

  //calling analysis function(s)! "event" here is the internal name of the tree in "fin"  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  TreetoHist (fin, "event", h1Ds, h2Ds, h3Ds, 1); //last argument (bool) tells the function this is data, so no weighting necessary
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

  //outro
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  //creating output file in which to deposit histograms
  //argv[1] should be the desired location of the output file, argv[2] should be the desired name
  TFile *fout = new TFile(((string) argv[1]+(string) argv[2]).c_str(),"RECREATE");
  cout << "DEBUG: output file name is " << fout->GetName() << endl;
  
  //writing hists to file
  for (unsigned i = 0; i < h1Ds.size(); ++ i) {
    h1Ds[i]->Write();
  }
  for (unsigned i = 0; i < h2Ds.size(); ++ i) {
    h2Ds[i]->Write();
  }
  for (unsigned i = 0; i < h3Ds.size(); ++ i) {
    h3Ds[i]->Write();
  }
  cout << "Wrote to " << fout->GetName() << endl;
  
  //closing file
  fout->Close();
  cout << "Closed " << fout->GetName() << endl;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  return 0;
}
