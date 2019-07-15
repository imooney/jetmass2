//Isaac Mooney, WSU - June 2019
//This file determines a preliminary bad tower list for any dataset [just change the file/hist name(s) to adapt]
//It takes in histograms of frequency and energy per tower from the output of the QA,
//finds the mean* and the RMS in order to drop towers outside of mean+-3*sigma.
//First towers are dropped by checking the per-event frequency, then the
//average energy, and finally the average energy above 2 GeV.
//A tower in the normal range for all three is a good tower.
//NOTE: this is currently used for pAu 2015, which has bad calibration so half the towers are higher average frequency/energy than the others and must be averaged separately. To use the file in general, just remove all of this doubling in the remove_towers() function.

#include <iostream>
#include <algorithm>
#include <set>
#include <fstream>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

using namespace std;

void remove_towers (set<int> & bad_tows, TH1D* h) {
  set<int>::iterator it;
  
  //oddity of pAu2015 dataset specifically: calibration differs for tower ID < 2400 and > 2400. So need two copies of everything.
  unsigned midpoint = (unsigned) h->GetNbinsX()/(double) 2;
  unsigned endpoint = (unsigned) h->GetNbinsX();
  double lowavg = 0, highavg = 0, lowsigma = 0, highsigma = 0;
  double nonzero_tows_low = 0, nonzero_tows_high = 0; //to count the nonzero towers

  cout << "DEBUG: midpoint = " << midpoint << "; endpoint = " << endpoint << endl;
  for (unsigned i = 1; i <= midpoint; ++ i) {//NOTE: towers index from 1, so do histograms. So it aligns perfectly!
    if (h->GetBinContent(i) && (bad_tows.find(i) == bad_tows.end())) {//empty or bad towers are ignored in calculating the average
      lowavg += h->GetBinContent(i);
      lowsigma += (h->GetBinContent(i)*h->GetBinContent(i));
      nonzero_tows_low ++;
    }
  }//tower loop low
  
  for (unsigned i = midpoint+1; i <= endpoint; ++ i) {
    if (h->GetBinContent(i) && (bad_tows.find(i) == bad_tows.end())) {//empty or bad towers are ignored in calculating the average
      highavg += h->GetBinContent(i);
      highsigma += (h->GetBinContent(i)*h->GetBinContent(i));
      nonzero_tows_high ++;
    }
  }//tower loop high
  
  //cout << "DEBUG: before dividing by total counts, lowavg = " << lowavg << " and highavg = " << highavg << endl;
  //cout << "DEBUG: before dividing by total counts and adjusting by mean^2, lowsigma = " << lowsigma << " and highsigma = " << highsigma << endl;
  //cout << "DEBUG: total counts low = " << nonzero_tows_low << " and high = " << nonzero_tows_high << endl;
  lowavg /= (double) nonzero_tows_low;
  highavg /= (double) nonzero_tows_high;
  lowsigma /= (double) nonzero_tows_low;
  highsigma /= (double) nonzero_tows_high;
  //cout << "DEBUG: before adjusting by mean^2, lowsigma = " << lowsigma << " and highsigma = " << highsigma << endl;
  lowsigma -= (lowavg*lowavg); //must be done after the mean is done being calculated
  highsigma -= (highavg*highavg); //ditto
  lowsigma = TMath::Abs(lowsigma); //shouldn't be negative
  highsigma = TMath::Abs(highsigma); //ditto
  lowsigma = sqrt(lowsigma); //standard deviation is square root of variance
  highsigma = sqrt(highsigma); //ditto
  cout << "DEBUG: final means: lowavg = " << lowavg << " and highavg = " << highavg << endl;
  cout << "DEBUG: final std. devs.: lowsigma = " << lowsigma << " and highsigma = " << highsigma << endl;

  cout << "DEBUG: final limits, low: " << lowavg - 3*lowsigma << " - " << lowavg + 3*lowsigma << endl
       << "                and high: " << highavg - 3*highsigma << " - " << highavg + 3*highsigma << endl;

  //fill the set of bad towers with the bad towers below 2400
  for (unsigned i = 1; i <= midpoint; ++ i) {
    if (h->GetBinContent(i) > lowavg + 3*lowsigma || h->GetBinContent(i) < lowavg - 3*lowsigma) {
      bad_tows.insert(i); //inserts the tower which was outside the window for a good tower
    }
  }//low tower ID bad tower fill loop

  //fill the set of bad towers with the bad towers above 2400
  for (unsigned i = midpoint+1; i <= endpoint; ++ i) {
    if (h->GetBinContent(i) > highavg + 3*highsigma || h->GetBinContent(i) < highavg - 3*highsigma) {
      bad_tows.insert(i); //inserts the tower which was outside the window for a good tower
    }
  }//high tower ID bad tower fill loop
  
  cout << "DEBUG: In the function, this is the current list: ";
  for (it=bad_tows.begin(); it!=bad_tows.end(); ++it) {cout << ' ' << *it << "[" << h->GetBinContent(*it) << "]";}
  cout << '\n';
  
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

  //opening file containing some example trees                                                                                                                 
  //argv[3] should be the name of the input file                                                                                                               
  string fin_name = (string) argv[3];
  TFile *fin = new TFile(fin_name.c_str(),"READ");
  cout << "DEBUG: input file name is " << fin->GetName() << endl;
  
  //pulling necessary histograms from file
  TH1D* hglob = (TH1D*) fin->Get("hn_trks"); //want a global observable so that the hist's integral is N_events, for scaling the frequency hist
  TH1D* hfreq_above2 = (TH1D*) fin->Get("htowerId_Eabove2"); //for scaling the weighted energy 1D by the counts to get the average E per tower
  TH1D* hfreq = (TH1D*) fin->Get("htowerId");
  TH1D* henergy = (TH1D*) fin->Get("htowerId_weighted");
  TH1D* henergy_above2 = (TH1D*) fin->Get("htowerId_weighted_above2");
  //to be filled without the bad towers:
  TH1D* hfreq_bad_removed = new TH1D("hfreq_bad_removed","",4800,0.5,4800.5);
  TH1D* henergy_bad_removed = new TH1D("henergy_bad_removed","",4800,0.5,4800.5);
  TH1D* henergy_above2_bad_removed = new TH1D("henergy_above2_bad_removed","",4800,0.5,4800.5);

  const int nTows = 4800;
   
  //turns the histogram y-axis into a per-tower average Et
  //  gStyle->SetErrorY(0);
  double errors[nTows] = {0};
  henergy->Divide(hfreq);
  henergy_above2->Divide(hfreq_above2);
  henergy->SetError(errors);
  henergy_above2->SetError(errors);
  
  //turns the histogram y-axis into per-event average - this must come after the henergy->Divide(hfreq) line
  hfreq->Scale(1/(double)hglob->Integral());

  //this set contains towers from the status tables marked as bad status at any point in the run period we are using
  //[any tower that only has a bad status in the runs that we omit does not appear in this list]
  set<int> bad_status_tows = {106, 139, 483, 562, 563, 564, 585, 603, 627, 680, 721, 722, 725, 916, 1023, 1154, 1158, 1204, 1221, 1237, 1283, 1301, 1388, 1434, 1440, 1563, 1564, 1567, 1602, 1612, 1654, 1753, 1762, 1776, 1952, 2092, 2415, 2569, 2697, 2929, 2969, 3005, 3186, 3360, 3436, 3498, 3504, 3682, 3702, 3739, 3834, 4006, 4047, 4104, 4177, 4331, 4496, 4500, 4659, 4712};
  
  set<int> bad_tows; //running list of the bad towers
  set<int>::iterator it;
  
  cout << "DEBUG: SIZE of bad tower list: " << bad_tows.size() << endl;
  //first, for towers with bad frequencies
  remove_towers(bad_tows, hfreq);
  cout << "DEBUG: SIZE of bad tower list: " << bad_tows.size() << endl;
  cout << "DEBUG: After the first function, this is the current list: ";
  for (it=bad_tows.begin(); it!=bad_tows.end(); ++it) {cout << ' ' << *it << "[" << hfreq->GetBinContent(*it) << "]";}
  cout << '\n';
  //next, for towers with bad energies
  remove_towers(bad_tows, henergy);
  cout << "DEBUG: SIZE of bad tower list: " << bad_tows.size() << endl;
  cout << "DEBUG: After the second function, this is the current list: ";
  for (it=bad_tows.begin(); it!=bad_tows.end(); ++it) {cout << ' ' << *it << "[" << henergy->GetBinContent(*it) << "]";}
  cout << '\n';
  //next, for towers with bad energies above 2 GeV
  remove_towers(bad_tows, henergy_above2);
  cout << "DEBUG: SIZE of bad tower list: " << bad_tows.size() << endl;
  cout << "DEBUG: After the third function, this is the current list: ";
  for (it=bad_tows.begin(); it!=bad_tows.end(); ++it) {cout << ' ' << *it << "[" << henergy_above2->GetBinContent(*it) << "]";}
  cout << '\n';

  //now we merge the bad towers with the bad status towers:
  set<int> combined_bad_tows;
  merge(bad_tows.begin(), bad_tows.end(), bad_status_tows.begin(), bad_status_tows.end(),
	inserter(combined_bad_tows, combined_bad_tows.begin()));

  //removes the bad towers from the histograms for visual aid
  for (unsigned i = 1; i <= nTows; ++ i) {
    if (combined_bad_tows.find(i) == combined_bad_tows.end()) {//means this is a good tower
      hfreq_bad_removed->Fill(i,hfreq->GetBinContent(i));//skips bad towers and fills what is otherwise a copy of the first histogram
      hfreq_bad_removed->SetBinError(i,hfreq->GetBinError(i));//keeps the same errors, too
      henergy_bad_removed->Fill(i,henergy->GetBinContent(i));
      henergy_bad_removed->SetBinError(i,henergy->GetBinError(i));
      henergy_above2_bad_removed->Fill(i,henergy_above2->GetBinContent(i));
      henergy_above2_bad_removed->SetBinError(i,henergy_above2->GetBinError(i));
    }
  }
  
  //decouples them from the current (input) file since they should be associated with the final output file
  hfreq->SetDirectory(0);
  henergy->SetDirectory(0);
  henergy_above2->SetDirectory(0);
  hfreq_bad_removed->SetDirectory(0);
  henergy_bad_removed->SetDirectory(0);
  henergy_above2_bad_removed->SetDirectory(0);
  fin->Close();
    
  //outro: writing the bad tower list and histograms without the bad towers
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~bad tower list~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  //opening file "bad_towers_pAu2015.list"
  fstream file("../lists/bad_towers_pAu2015.list", fstream::in | fstream::out | fstream::trunc);
  if (!file) {
    cerr << "Error in creating file!!!" << endl; exit(1);
  }
  else {
    cout << "bad_towers_pAu2015.list created successfully." << endl;
  }
  
  file << "#Bad (hot+dead+bad_status) towers from 2015 pAu BBCMB data. " << combined_bad_tows.size() << " total bad towers. (iam; july 2, 2019)" << endl;
  file << "#including " << bad_status_tows.size() << " bad towers from the status tables for the 2015 pAu data which were not caught by the 3-sigma removal procedure:" << endl
       << "#563, 564, 585, 603, 627, 680, 721, 722, 725, 1023, 1154, 1221, 1301, 1563, 1564, 1602, 1776, 1952, 2697, 3504, 3682, 3702, 3834, 4104, 4496, 4500, 4659, 4712" << endl;
  //writing each bad tower ID to the file
  for (it=combined_bad_tows.begin(); it!= --combined_bad_tows.end(); ++it) {file << *it << ", ";}
  it = --combined_bad_tows.end(); //grabs the last element in the list separately, so we don't have a comma after the last element
  file << *it; //writes that element
  
  file.close();
  cout << "Closed bad_towers_pAu2015.list" << endl;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                                                        
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~histograms~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  //creating output file in which to deposit histograms                                                                                                        
  //argv[1] should be the desired location of the output file, argv[2] should be the desired name                                                              
  TFile *fout = new TFile(((string) argv[1]+(string) argv[2]).c_str(),"RECREATE");
  cout << "DEBUG: output file name is " << fout->GetName() << endl;

  //writing hists to file         
  hfreq->Write(); henergy->Write(); henergy_above2->Write();
  hfreq_bad_removed->Write(); henergy_bad_removed->Write(); henergy_above2_bad_removed->Write();
  
  cout << "Wrote to " << fout->GetName() << endl;

  //closing file                                                                                                                                               
  fout->Close();
  cout << "Closed " << fout->GetName() << endl;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//  
  
  return 0;
  
}
