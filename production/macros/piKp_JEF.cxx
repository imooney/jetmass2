//Isaac Mooney, WSU, September 2019
//this file takes in a jet tree and calculates the fraction of jet energy carried by pions, Kaons, and protons as a function of pT.

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
void TreetoHist (/*TFile *fin,*/ vector<TH1D*> h1Ds, const char* argv1, const char* argv2, const char* argv3) {  
   //opening file containing some example trees
  //argv[3] should be the name of the input file
  string fin_name = (string) argv3;
  TFile *fin = new TFile(fin_name.c_str(),"READ");
  cout << "DEBUG: input file name is " << fin->GetName() << endl;

  TTree *t = (TTree*) fin->Get("ResultTree"); 
  
  //initializing the variables that will be filled by the values in the branches later
  vector<double> *Pt = 0;
  vector<vector<double> > *PID = 0;
  vector<vector<double> > *consPt = 0;
  
  double weight = -1;
  
  //linking the branches with the variables
  t->SetBranchAddress("consPID",&PID);
  t->SetBranchAddress("conspT",&consPt);  
  t->SetBranchAddress("mcweight", &weight);
  t->SetBranchAddress("jetpT",&Pt);

  
  vector<double> JEF_pi(15),JEF_K(15),JEF_p(15),JEF_other(15),JEF_tot(15);
  
  cout << "RUNNING OVER TREE ResultTree! Entries: " << t->GetEntries() << endl;
  const clock_t begin_time = clock(); //timing - for debugging and for fun
  for (unsigned i = 0; i < t->GetEntries(); ++ i) { //"event" loop
    if (i % 1000000 == 0 && i != 0) { //can change this to a frequency of your preference (for real data I use 1e5 or 1e6)
      cout << "Still chuggin. On event " << i << endl;
      cout << "Total time passed: " << fixed << setprecision(5) << double(clock() - begin_time) /(double) CLOCKS_PER_SEC << " secs" << endl;
    }
    t->GetEntry(i);
    
    //looping over jets and filling histograms
    for (unsigned j = 0; j < Pt->size(); ++ j) { //all vectors of doubles in the branches should have the same size
      
      for (unsigned k = 0; k < PID->at(j).size(); ++ k) { //constituent loop
	//~~~this block will handle the pi/K/p jet energy fraction calculation~~~//
	if (fabs(PID->at(j).at(k)) == (double) 211 || PID->at(j).at(k) == (double) 111) {
	  //cout << "DEBUG: PIONS: PID = " << PID->at(j).at(k) << " should be |211| or 111" << endl;
	  //for each jet pT we have a running count of the JEF of pions. We therefore index the vector of pions depending on the jet's pT. E.g. for 5 GeV bins starting at 5 GeV, a 20 GeV jet is in the 3rd bin (indexing from 0) or: (20 - 5)/5. 
	  JEF_pi[(int)((Pt->at(j)-5)/(double)5)] += consPt->at(j).at(k)*weight;
	  //cout << "DEBUG: PIONS: " << weight << " " << consPt->at(j).at(k)*weight << endl;
	}
	else if (fabs(PID->at(j).at(k)) == (double) 321 || PID->at(j).at(k) == (double) 311 || PID->at(j).at(k) == (double) 310 || PID->at(j).at(k) == (double) 130) {
	  //cout << "DEBUG: KAONS: PID = " << PID->at(j).at(k) << " should be |321|, 311, 310, or 130" << endl;
	  JEF_K[(int)((Pt->at(j)-5)/(double)5)] += consPt->at(j).at(k)*weight;
	  //cout << "DEBUG: KAONS: " << weight << " " << consPt->at(j).at(k)*weight << endl;
	}
	else if (fabs(PID->at(j).at(k)) == (double) 2212) {
	  //cout << "DEBUG: PROTONS: PID = " << PID->at(j).at(k) << " should be |2212|" << endl;
	  JEF_p[(int)((Pt->at(j)-5)/(double)5)] += consPt->at(j).at(k)*weight;
	  //cout << "DEBUG: PROTONS: " << weight << " " << consPt->at(j).at(k)*weight << endl;
	}
	else {
	  //cout << "DEBUG: OTHERS: PID = " << PID->at(j).at(k) << " should be anything but |211|, 111, |321|, 311, 310, 130, or |2212|" << endl;
	  JEF_other[(int)((Pt->at(j)-5)/(double)5)] += consPt->at(j).at(k)*weight;
	  //cout << "DEBUG: OTHERS: " << weight << " " << consPt->at(j).at(k)*weight << endl;
	}
	JEF_tot[(int)((Pt->at(j)-5)/(double)5)] += consPt->at(j).at(k)*weight;
	//cout << "DEBUG: TOTAL: " << weight << " " << consPt->at(j).at(k)*weight << endl;
	//~~~          ~~~//
      }//!constituent loop
    }//!jet loop
  }//!event loop

  //normalizing the JEFs
  cout << "jef_tot.size(): " << JEF_tot.size() << endl;
  for (unsigned i = 0; i < JEF_tot.size(); ++ i) {
    cout << "DEBUG: " << JEF_pi[i] << " " << JEF_K[i] << " " << JEF_p[i] << " " << JEF_other[i] << " " << JEF_tot[i] << endl;
    if (JEF_tot[i] != 0) {
      cout << "A1" << endl;
      JEF_pi[i] /= (double) JEF_tot[i];
      JEF_K[i] /= (double) JEF_tot[i];
      JEF_p[i] /= (double) JEF_tot[i];
      JEF_other[i] /= (double) JEF_tot[i];
      cout << "A2" << endl;
    }
    else { //to avoid dividing by zero when there are no jets in the bin
      cout << "A3" << endl;
      JEF_pi[i] = 0;
      JEF_K[i] = 0;
      JEF_p[i] = 0;
      JEF_other[i] = 0;
      cout << "A4" << endl;
    }
    cout << "B" << endl;
    h1Ds[0]->Fill(5*i+5,JEF_pi[i]);//re-converts the bin number to a jet pT value for filling the hist
    h1Ds[1]->Fill(5*i+5,JEF_K[i]);
    h1Ds[2]->Fill(5*i+5,JEF_p[i]);
    h1Ds[3]->Fill(5*i+5,JEF_other[i]);
    cout << "C" << endl;
  }
  
  cout << "D" << endl;
  
  //!needs to be outside the event loop; not sure exactly what it does
  //t->ResetBranchAddresses();
  
  
  //outro
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  //creating output file in which to deposit histograms
  //argv[1] should be the desired location of the output file, argv[2] should be the desired name
  cout << "DEBUG: " << endl;
  cout << argv1 << " " << argv2 << endl;
  TFile *fout = new TFile(((string) argv1+(string) argv2).c_str(),"RECREATE");
  cout << "DEBUG: output file name is " << fout->GetName() << endl;
  
  //writing hists to file
  for (unsigned i = 0; i < h1Ds.size(); ++ i) {
    h1Ds[i]->Write();
  }
 
  cout << "Wrote to " << fout->GetName() << endl;
  
  //closing file
  fout->Close();
  cout << "Closed " << fout->GetName() << endl;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

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
  /*
  //opening file containing some example trees
  //argv[3] should be the name of the input file
  string fin_name = (string) argv[3];
  TFile *fin = new TFile(fin_name.c_str(),"READ");
  cout << "DEBUG: input file name is " << fin->GetName() << endl;
*/
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~hists~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

  //jet energy fraction by species, and delta_R & z as function of species
  TH1D* ptvJEF_pi = new TH1D("ptvJEF_pi","",15,5,80);
  TH1D* ptvJEF_K = new TH1D("ptvJEF_K","",15,5,80);
  TH1D* ptvJEF_p = new TH1D("ptvJEF_p","",15,5,80);
  TH1D* ptvJEF_other = new TH1D("ptvJEF_other","",15,5,80);
  TH2D* deltaRvpt_pi = new TH2D("deltaRvpt_pi","",10,0,0.5,15,5,80);
  TH2D* deltaRvpt_K = new TH2D("deltaRvpt_K","",10,0,0.5,15,5,80);
  TH2D* deltaRvpt_p = new TH2D("deltaRvpt_p","",10,0,0.5,15,5,80);
  TH2D* deltaRvpt_other = new TH2D("deltaRvpt_other","",10,0,0.5,15,5,80);
  TH2D* zvpt_pi = new TH2D("zvpt_pi","",10,0,1,15,5,80);
  TH2D* zvpt_K = new TH2D("zvpt_K","",10,0,1,15,5,80);
  TH2D* zvpt_p = new TH2D("zvpt_p","",10,0,1,15,5,80);
  TH2D* zvpt_other = new TH2D("zvpt_other","",10,0,1,15,5,80);
  
  //need 4 1D hists of jet pT weighted by the JEF. How is JEF calculated? Add up all the pTs of the constituents and divide by the jet pT. Then weight that number with the MC weight. Keep a running tally of these numbers. Add the next jet's weighted number to this first number, and also to the running total. At the end, divide each of the four weighted numbers by the weighted total to get the JEF for each of the species. Stack them here. Or create the stack in another plotting macro if ownership is weird w/r/t prettification later.
  //do the same thing at the constituent level with z and delta_R?
 
  //putting them in a vector to more easily shuttle them back and forth in the function. Drawback: have to know their order.
  vector<TH1D*> h1Ds = {ptvJEF_pi, ptvJEF_K, ptvJEF_p, ptvJEF_other};
  //  vector<TH2D*> h2Ds = {deltaRvpt_pi,deltaRvpt_K,deltaRvpt_p,deltaRvpt_other,zvpt_pi,zvpt_K,zvpt_p,zvpt_other};
  
  //calling analysis function(s)! "event" here is the internal name of the tree in "fin"  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  //  TTree *t = (TTree*) fin->Get("ResultTree");
  TreetoHist (/*fin,*/ h1Ds, argv[1], argv[2], argv[3]);
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
/*
  //outro
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  //creating output file in which to deposit histograms
  //argv[1] should be the desired location of the output file, argv[2] should be the desired name
  cout << "DEBUG: " << endl;
  cout << argv[1] << " " << argv[2] << endl;
  TFile *fout = new TFile(((string) argv[1]+(string) argv[2]).c_str(),"RECREATE");
  cout << "DEBUG: output file name is " << fout->GetName() << endl;
  
  //writing hists to file
  for (unsigned i = 0; i < h1Ds.size(); ++ i) {
    h1Ds[i]->Write();
  }
 
  cout << "Wrote to " << fout->GetName() << endl;
  
  //closing file
  fout->Close();
  cout << "Closed " << fout->GetName() << endl;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  */
  return 0;
}
