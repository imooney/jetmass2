//Isaac Mooney, WSU, August 2019
//This file addresses a problem in RooUnfold pertaining to the statistical errors.
//RooUnfold seems to assume that the "truth" jets which were not matched to "measured" jets still count toward the number of unfolded jets.
//So the statistical errors are calculated basically by the square root of the number of unfolded jets, which has a contribution from these "misses".
//This is bad statistics. We should really be scaling the error from the data jets by some factor to account for inefficiency in the truth -> measured process.
//We should then take the errors given by RooUnfold and bin-by-bin scale them up by x = the square root of (N_matches + N_misses) / N_matches.
//This comes from a simple algebraic rearrangement of sqrt(N_matches + N_misses) * x = sqrt(N_matches) * (N_matches + N_misses) / N_matches;
//where the LHS is asking the question "By what do we need to scale the errors given by RooUnfold to match the RHS in which we properly scale the data counts by the inefficiency factor.
//So, in this macro, we take in the appropriate histograms (from a file produced from MC requiring matched Pythia and Pythia+Geant events, to eliminate the unnecessary contribution of event matching efficiency to the overall efficiency), and calculate this number, x, for each bin, outputting it as a histogram for use after the unfolding.

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

#include "RooUnfold.h"

using namespace std;

int main (int argc, const char** argv) {
  //intro                                                                                                                                                      
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                                                       
  //basic argument checking.                                                                                                                                   
  if (argc != 2) {
    cerr << "Should be one argument: jet radius. Received "
         << argc-1 << ". Exiting." << endl;
    exit(1);
  }  
   
  //argv[1] should be the jet radius e.g. "04".                                                                                                               
  string radius = (string) argv[1];
  radius = "_R"+radius;//appending the "_R" used in file naming.
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();

  const string path = "~/jetmass2/out/sim/";
  const string file_in = "sim_matched";

  TFile *fin = new TFile((path+file_in+radius+".root").c_str(),"READ");
  cout << "DEBUG: input file name is " << fin->GetName() << endl;
  
  RooUnfoldResponse* res = (RooUnfoldResponse*) fin->Get("pt_response");
  TH1D* match_plus_miss = (TH1D*) fin->Get("pt_gen_match_plus_miss");
  TH1D* match = (TH1D*) res->Hresponse()->ProjectionY("match");
  
  TH1D* hratio = (TH1D*) match_plus_miss->Clone("hratio");
  TH1D* efficiency = (TH1D*) match->Clone("efficiency");
  
  hratio->Divide(match);
  efficiency->Divide(match_plus_miss);
  
  const int nBins = hratio->GetNbinsX();
  
  for (int i = 1; i <= nBins; ++ i) {
    hratio->SetBinContent(i,sqrt(hratio->GetBinContent(i)));
  }
  
  TFile *fout = new TFile((path+"stat_err_scaling"+radius+".root").c_str(),"RECREATE");

  fout->cd();
  hratio->Write();
  efficiency->Write();

  cout << endl
       << "Wrote statistical error scaling to " << fout->GetName() << endl;

  //closing files                                                                                                                                              
  fout->Close();                                                              
  cout << endl << "Closed " << fout->GetName() << endl;
  
  return 0;
}
