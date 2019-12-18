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
void TreetoHist (TFile *f/*, string treestr, vector<string> branchnames, vector<TH1D*> h1Ds*/, vector<TH2D*> h2Ds, vector<TH3D*> h3Ds, bool hadronic) {  
  cout << "STARTING TREETOHIST" << endl;
  
  //initializing the variables that will be filled by the values in the branches later
  vector<double> *Pt = 0; vector<double> *Eta = 0; vector<double> *Phi = 0;
  vector<double> *M = 0; vector<double> *zg = 0; vector<double> *rg = 0; vector<double> *mg = 0;
  vector<int> *qvg = 0;
  
  double weight = -1;
  
  cout << "INITIALIZED BRANCH VARIABLES!" << endl;
  
  //getting the tree and linking the branches with the variables
  TTree *t = (TTree*) f->Get("ResultTree"/*treestr.c_str()*/);
  t->SetBranchAddress("jetpT"/*branchnames[0].c_str()*/,&Pt);
  t->SetBranchAddress("jeteta"/*branchnames[1].c_str()*/,&Eta);
  t->SetBranchAddress("jetphi"/*branchnames[2].c_str()*/,&Phi);
  t->SetBranchAddress("jetM"/*branchnames[3].c_str()*/,&M);
  t->SetBranchAddress("sdjetM"/*branchnames[4].c_str()*/,&mg);
  t->SetBranchAddress("zg"/*branchnames[5].c_str()*/,&zg);
  t->SetBranchAddress("rg"/*branchnames[6].c_str()*/,&rg);

  if (hadronic) {
    t->SetBranchAddress("qvg",&qvg);
  }
  
  t->SetBranchAddress("mcweight", &weight);

  cout << "ASSIGNED VARIABLES TO BRANCHES!" << endl;
  
  cout << /*(*/"RUNNING OVER TREE ResultTree"/*+treestr+*/ << "! Entries: "/*).c_str()*/ << t->GetEntries() << endl;
  const clock_t begin_time = clock(); //timing - for debugging and for fun
  for (unsigned i = 0; i < t->GetEntries(); ++ i) { //"event" loop
    //cout << "IN EVENT LOOP" << endl;
    if (i % 1000 == 0 && i != 0) { //can change this to a frequency of your preference (for real data I use 1e5 or 1e6)
      cout << "Still chuggin. On event " << i << endl;
      cout << "Total time passed: " << fixed << setprecision(5) << double(clock() - begin_time) /(double) CLOCKS_PER_SEC << " secs" << endl;
    }
    //    cout << "GETTING ENTRY " << i << endl;
    t->GetEntry(i);
    cout << "GOT ENTRY " << i << endl;
    
    //filling "event" observables
    //
    //looping over jets and filling histograms
    for (unsigned j = 0; j < Pt->size(); ++ j) { //all vectors of doubles in the branches should have the same size
      //cout << "IN JET LOOP" << endl;
      //2Ds!                                                                                                                                           
      h2Ds[0]->Fill(M->at(j), Pt->at(j), weight);
      //quark v. gluon jets
      if (hadronic && qvg->at(j) == 0) {//q jet!                                                                                                                           
        h2Ds[1]->Fill(M->at(j),Pt->at(j),weight);
      }
      else if (hadronic && qvg->at(j) == 1) {//g jet!                                                                                                                      
        h2Ds[2]->Fill(M->at(j),Pt->at(j),weight);
      }
      else {//neither (or no match to hard parton)!
        h2Ds[3]->Fill(M->at(j),Pt->at(j),weight);
      }

      h2Ds[4]->Fill(Eta->at(j), Pt->at(j), weight);
      h2Ds[5]->Fill(Phi->at(j), Pt->at(j), weight);
      if (zg->at(j) >= 0.1) { //otherwise, we tag and drop the groomed jet
	h2Ds[6]->Fill(mg->at(j), Pt->at(j), weight);
	h2Ds[7]->Fill(zg->at(j), Pt->at(j), weight);
	h2Ds[8]->Fill(rg->at(j), Pt->at(j), weight);
      
	//3Ds!
	h3Ds[0]->Fill(M->at(j),mg->at(j),Pt->at(j),weight);
	h3Ds[1]->Fill(M->at(j),zg->at(j),Pt->at(j),weight);
	h3Ds[2]->Fill(M->at(j),rg->at(j),Pt->at(j),weight);
	h3Ds[3]->Fill(mg->at(j),zg->at(j),Pt->at(j),weight);
	h3Ds[4]->Fill(mg->at(j),rg->at(j),Pt->at(j),weight);
	h3Ds[5]->Fill(zg->at(j),rg->at(j),Pt->at(j),weight);
      } 
    }//!jet loop
    //    cout << "DONE WITH JET LOOP" << endl;
  }//!event loop
  cout << "DONE WITH EVENT LOOP" << endl;
  cout << "HeLP!" << endl;
  //!needs to be outside the event loop; not sure exactly what it does
  t->ResetBranchAddresses();
  cout << "HeLP2!" << endl;
  //f->Close();
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
  //string fin_name = (string) argv[3];
  cout << "WHAT ARE YOUUUUUUUUU? " << (string) argv[3] << endl;
  string fin_name; fin_name.assign((string) argv[3]);//strcpy(fin_name,(string) argv[3]);
  TFile *fin = new TFile(fin_name.c_str(),"READ");
  cout << "DEBUG: input file name is " << fin->GetName() << endl;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~hists~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

  //event observables:
  //
  
  //pT range for particle-level is different from detector-level so define it here for ease of change later if necessary:
  const int ptbins = 15;
  const double ptlow = 5; const double pthigh = 80;
  
  //~~~HADRON LEVEL~~~
  //2D correlations with pT
  TH2D* mvpt = new TH2D("mvpt",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,ptbins,ptlow,pthigh);
  
  TH2D* mvpt_q = new TH2D("mvpt_q",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,ptbins,ptlow,pthigh);//quark jets
  TH2D* mvpt_g = new TH2D("mvpt_g",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,ptbins,ptlow,pthigh);//gluon jets
  TH2D* mvpt_n = new TH2D("mvpt_n",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,ptbins,ptlow,pthigh);//neither
  
  TH2D* mgvpt = new TH2D("mgvpt",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,ptbins,ptlow,pthigh);
  TH2D* etavpt = new TH2D("etavpt",";#eta;p_{T} [GeV/c]",50,-1,1,ptbins,ptlow,pthigh);
  TH2D* phivpt = new TH2D("phivpt",";#phi;p_{T} [GeV/c]",50,0,2*M_PI,ptbins,ptlow,pthigh);
  TH2D* zgvpt = new TH2D("zgvpt",";z_{g} [GeV/c^{2}];p_{T} [GeV/c]",10,0,0.5,ptbins,ptlow,pthigh);
  TH2D* rgvpt = new TH2D("rgvpt",";R_{g} [GeV/c^{2}];p_{T} [GeV/c]",10,0,0.5,ptbins,ptlow,pthigh);

  //3D correlations
  TH3D* mvmgvpt = new TH3D("mvmgvpt",";M [GeV/c^{2}];M_{g} [GeV/c^{2}]; p_{T} [GeV/c]",14,0,14,14,0,14,ptbins,ptlow,pthigh);
  TH3D* mvzgvpt = new TH3D("mvzgvpt",";M [GeV/c^{2}];z_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,ptbins,ptlow,pthigh);
  TH3D* mvrgvpt = new TH3D("mvrgvpt",";M [GeV/c^{2}];R_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,ptbins,ptlow,pthigh);
  TH3D* mgvzgvpt = new TH3D("mgvzgvpt",";M_{g} [GeV/c^{2}];z_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,ptbins,ptlow,pthigh);
  TH3D* mgvrgvpt = new TH3D("mgvrgvpt",";M_{g} [GeV/c^{2}];R_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,ptbins,ptlow,pthigh);
  TH3D* zgvrgvpt = new TH3D("zgvrgvpt",";z_{g};R_{g}; p_{T} [GeV/c]",10,0,0.5,10,0,0.5,ptbins,ptlow,pthigh);

  //putting them in a vector to more easily shuttle them back and forth in the function. Drawback: have to know their order.
  vector<TH2D*> h2Ds = {mvpt,mvpt_q,mvpt_g,mvpt_n,etavpt,phivpt,mgvpt,zgvpt,rgvpt};
  vector<TH3D*> h3Ds = {mvmgvpt,mvzgvpt,mvrgvpt,mgvzgvpt,mgvrgvpt,zgvrgvpt}; 
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  
  cout << "HISTS DONE!" << endl;
  
  //calling analysis function(s)! "event" here is the internal name of the tree in "fin"  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  TreetoHist (fin, h2Ds, h3Ds, 1); //1 = hadronic
  cout << "HeLP3!" << endl;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  
  //outro
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  //creating output file in which to deposit histograms
  //argv[1] should be the desired location of the output file, argv[2] should be the desired name
  TFile *fout = new TFile(((string) argv[1]+(string) argv[2]).c_str(),"RECREATE");
  cout << "DEBUG: output file name is " << fout->GetName() << endl;
  
  //writing hists to file
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
  fin->Close();
  cout << "Closed " << fin->GetName() <<endl;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  return 0;
}
