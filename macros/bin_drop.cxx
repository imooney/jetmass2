//Isaac Mooney, WSU, August 2019                                                                                                         
//This file takes in all output of the main analysis to this point, drops low statistics bins, and outputs the results to similarly named files

//First, we drop data spectra and closure pseudo-data spectra with fewer than 20 jets.
//Next, we drop the bins in all responses corresponding to the bins dropped in the closure pseudo-data spectra

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

//loops over the 2D spectra to find bins with bin content < 20. These bins are passed to the response, to be dropped there as well.
//PS this assumes mass (or whatever observable) is on the x-axis!
vector<int> DropLowStatsBins(RooUnfoldResponse* weighted, TH2D* det_unweighted) {//, TH2D* gen) {
  //TH2D* unweighted_res = (TH2D*) unweighted->Hresponse();

  vector<int> bins_to_drop;

  //this method requires the binning for the TH2 to be the same as for the RooUnfoldResponse
  for (int i = 1; i <= det_unweighted->GetXaxis()->GetNbins(); ++ i) {//mass bins in det-level spectrum
    for (int j = 1; j <= det_unweighted->GetYaxis()->GetNbins(); ++ j) {//pT bins in det-level spectrum
      if (det_unweighted->GetBinContent(i,j) < 20) {
	cout << "For mass bin " << i << " and pT bin " << j << ": " << endl;
	cout << "Removing column: " << (det_unweighted->GetXaxis()->GetNbins()*(j-1))+i << endl;
	//for a given det-level (M,pT) with lacking statistics, we have to remove the whole gen-level ~column~ -> no PL jets can wind up here at DL
	bins_to_drop.push_back((det_unweighted->GetXaxis()->GetNbins()*(j-1))+i);
	for (int k = 1; k <= weighted->Hresponse()->GetYaxis()->GetNbins(); ++ k) {//gen bins in response
	  weighted->Hresponse()->SetBinContent((det_unweighted->GetXaxis()->GetNbins()*(j-1))+i,k,0);
	  weighted->Hresponse()->SetBinError((det_unweighted->GetXaxis()->GetNbins()*(j-1))+i,k,0);                                      
	}                                                                                                                      
      }                                                                                                                       
    }
  }
  return bins_to_drop;
}


//loops over the unweighted 1D histogram and finds any bins where there are fewer than 20 counts. These statistics-limited bins are dropped from the corresponding weighted histogram.
vector<int> DropLowStatsBins(TH1D* weighted, TH1D* unweighted) {     
  
  vector<int> bins_to_drop;
  
  for (int i = 1; i <= unweighted->GetXaxis()->GetNbins(); ++ i) {                                                         
    if (unweighted->GetBinContent(i) < 20) {                                                                                  
      bins_to_drop.push_back(i);
      weighted->SetBinContent(i,0);                                                                                          
      weighted->SetBinError(i,0);                                                                                             
    }                                                                                                                            
  }                                                                                                                               

  return bins_to_drop;
}                                                                                                                                                         
                                                                                                                                                            
//loops over the unweighted 2D histogram and finds any bins where there are fewer than 20 counts. These statistics-limited bins are dropped from the corresponding weighted histogram.
void DropLowStatsBins(TH2D* weighted, TH2D* unweighted) {                                                                
  for (int i = 1; i <= unweighted->GetXaxis()->GetNbins(); ++ i) {                                                        
    for (int j = 1; j <= unweighted->GetYaxis()->GetNbins(); ++ j) {                                                            
      if (unweighted->GetBinContent(i,j) < 20) {                                                                             
	weighted->SetBinContent(i,j,0);                                                                                      
	weighted->SetBinError(i,j,0);                                                                                        
      }                                                                                                                  
    }                                                                                                                            
  }                                                                                                                                                       
  return;                                                                                                                                                 
} 

void DropBins(RooUnfoldResponse *response, vector<int> to_drop) {
  for (unsigned i = 0; i < to_drop.size(); ++ i) { //make sure you index from 0 to size - 1 here since you're using the vector.
    for (int j = 1; j <= response->Hresponse()->GetYaxis()->GetNbins(); ++ j) { //indexing response here so 1 to size.
      response->Hresponse()->SetBinContent(to_drop[i], j, 0);
      response->Hresponse()->SetBinError(to_drop[i], j, 0);
    }
  }
  
  return;
}

int main () {
  //intro                                                                                                                                                      
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                                                       
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();
  
  const string sim_path = "~/jetmass2/out/sim/hists/";
  const string match_path = "~/jetmass2/out/sim/";
  const string data_path = "~/jetmass2/out/data/hists/";
  const string match_file = "sim_matched_closureconsistency";
  const string data_file = "data_hists_ppJP2";
  const string sim_file = "unmatched_hists";
  
  TFile *match_in = new TFile((match_path+match_file+".root").c_str(),"READ");
  TFile *data_in = new TFile((data_path+data_file+".root").c_str(),"READ");
  TFile *sim_in = new TFile((sim_path+sim_file+".root").c_str(),"READ");
  
  //~~~closure - must be done first
  RooUnfoldResponse *sampleA_m_pt_response = (RooUnfoldResponse*) match_in->Get("sampleA_m_pt_response");
  RooUnfoldResponse *sampleA_mg_pt_response = (RooUnfoldResponse*) match_in->Get("sampleA_mg_pt_response");
  RooUnfoldResponse *sampleB_m_pt_response = (RooUnfoldResponse*) match_in->Get("sampleB_m_pt_response");
  RooUnfoldResponse *sampleB_mg_pt_response = (RooUnfoldResponse*) match_in->Get("sampleB_mg_pt_response");
  RooUnfoldResponse *sampleA_pt_response = (RooUnfoldResponse*) match_in->Get("sampleA_pt_response");
  RooUnfoldResponse *sampleB_pt_response = (RooUnfoldResponse*) match_in->Get("sampleB_pt_response");
  RooUnfoldResponse *sampleA_m_response = (RooUnfoldResponse*) match_in->Get("sampleA_m_response");
  RooUnfoldResponse *sampleB_m_response = (RooUnfoldResponse*) match_in->Get("sampleB_m_response");
  RooUnfoldResponse *sampleA_mg_response = (RooUnfoldResponse*) match_in->Get("sampleA_mg_response");
  RooUnfoldResponse *sampleB_mg_response = (RooUnfoldResponse*) match_in->Get("sampleB_mg_response");
 
  //1D closure spectra to be pruned and to be used to prune responses
  TH1D *sampleA_pt_gen = (TH1D*) match_in->Get("sampleA_pt_gen");
  TH1D *sampleA_pt_det = (TH1D*) match_in->Get("sampleA_pt_det");
  TH1D *sampleB_pt_gen = (TH1D*) match_in->Get("sampleB_pt_gen");
  TH1D *sampleB_pt_det = (TH1D*) match_in->Get("sampleB_pt_det");
  TH1D *sampleA_pt_gen_counts = (TH1D*) match_in->Get("sampleA_pt_gen_counts");
  TH1D *sampleA_pt_det_counts = (TH1D*) match_in->Get("sampleA_pt_det_counts");
  TH1D *sampleB_pt_gen_counts = (TH1D*) match_in->Get("sampleB_pt_gen_counts");
  TH1D *sampleB_pt_det_counts = (TH1D*) match_in->Get("sampleB_pt_det_counts");
  
  TH1D *sampleA_m_gen = (TH1D*) match_in->Get("sampleA_m_gen");
  TH1D *sampleA_m_det = (TH1D*) match_in->Get("sampleA_m_det");
  TH1D *sampleB_m_gen = (TH1D*) match_in->Get("sampleB_m_gen");
  TH1D *sampleB_m_det = (TH1D*) match_in->Get("sampleB_m_det");
  TH1D *sampleA_m_gen_counts = (TH1D*) match_in->Get("sampleA_m_gen_counts");
  TH1D *sampleA_m_det_counts = (TH1D*) match_in->Get("sampleA_m_det_counts");
  TH1D *sampleB_m_gen_counts = (TH1D*) match_in->Get("sampleB_m_gen_counts");
  TH1D *sampleB_m_det_counts = (TH1D*) match_in->Get("sampleB_m_det_counts");

  TH1D *sampleA_mg_gen = (TH1D*) match_in->Get("sampleA_mg_gen");
  TH1D *sampleA_mg_det = (TH1D*) match_in->Get("sampleA_mg_det");
  TH1D *sampleB_mg_gen = (TH1D*) match_in->Get("sampleB_mg_gen");
  TH1D *sampleB_mg_det = (TH1D*) match_in->Get("sampleB_mg_det");
  TH1D *sampleA_mg_gen_counts = (TH1D*) match_in->Get("sampleA_mg_gen_counts");
  TH1D *sampleA_mg_det_counts = (TH1D*) match_in->Get("sampleA_mg_det_counts");
  TH1D *sampleB_mg_gen_counts = (TH1D*) match_in->Get("sampleB_mg_gen_counts");
  TH1D *sampleB_mg_det_counts = (TH1D*) match_in->Get("sampleB_mg_det_counts");

  //2D closure spectra to be pruned and to be used to prune responses
  TH2D *sampleA_m_pt_gen = (TH2D*) match_in->Get("sampleA_m_pt_gen");
  TH2D *sampleA_m_pt_det = (TH2D*) match_in->Get("sampleA_m_pt_det");
  TH2D *sampleB_m_pt_gen = (TH2D*) match_in->Get("sampleB_m_pt_gen");
  TH2D *sampleB_m_pt_det = (TH2D*) match_in->Get("sampleB_m_pt_det");
  TH2D *sampleA_m_pt_gen_counts = (TH2D*) match_in->Get("sampleA_m_pt_gen_counts");
  TH2D *sampleA_m_pt_det_counts = (TH2D*) match_in->Get("sampleA_m_pt_det_counts");
  TH2D *sampleB_m_pt_gen_counts = (TH2D*) match_in->Get("sampleB_m_pt_gen_counts");
  TH2D *sampleB_m_pt_det_counts = (TH2D*) match_in->Get("sampleB_m_pt_det_counts");
  TH2D *sampleA_mg_pt_gen = (TH2D*) match_in->Get("sampleA_mg_pt_gen");
  TH2D *sampleA_mg_pt_det = (TH2D*) match_in->Get("sampleA_mg_pt_det");
  TH2D *sampleB_mg_pt_gen = (TH2D*) match_in->Get("sampleB_mg_pt_gen");
  TH2D *sampleB_mg_pt_det = (TH2D*) match_in->Get("sampleB_mg_pt_det");
  TH2D *sampleA_mg_pt_gen_counts = (TH2D*) match_in->Get("sampleA_mg_pt_gen_counts");
  TH2D *sampleA_mg_pt_det_counts = (TH2D*) match_in->Get("sampleA_mg_pt_det_counts");
  TH2D *sampleB_mg_pt_gen_counts = (TH2D*) match_in->Get("sampleB_mg_pt_gen_counts");
  TH2D *sampleB_mg_pt_det_counts = (TH2D*) match_in->Get("sampleB_mg_pt_det_counts");

  //~~~matched MC
  RooUnfoldResponse *m_pt_response = (RooUnfoldResponse*) match_in->Get("m_pt_response");
  RooUnfoldResponse *mg_pt_response = (RooUnfoldResponse*) match_in->Get("mg_pt_response");
  
  //~~~systematics
  RooUnfoldResponse *m_pt_res_nom = (RooUnfoldResponse*) match_in->Get("m_pt_res_nom");
  RooUnfoldResponse *m_pt_res_TS = (RooUnfoldResponse*) match_in->Get("m_pt_res_TS");
  RooUnfoldResponse *m_pt_res_TU = (RooUnfoldResponse*) match_in->Get("m_pt_res_TU");
  RooUnfoldResponse *m_pt_res_HC50 = (RooUnfoldResponse*) match_in->Get("m_pt_res_HC50");
  RooUnfoldResponse *m_pt_res_DS = (RooUnfoldResponse*) match_in->Get("m_pt_res_DS");
  RooUnfoldResponse *m_pt_res_GS = (RooUnfoldResponse*) match_in->Get("m_pt_res_GS");
  RooUnfoldResponse *m_pt_res_MS = (RooUnfoldResponse*) match_in->Get("m_pt_res_MS");
  RooUnfoldResponse *mg_pt_res_nom = (RooUnfoldResponse*) match_in->Get("mg_pt_res_nom");
  RooUnfoldResponse *mg_pt_res_TS = (RooUnfoldResponse*) match_in->Get("mg_pt_res_TS");
  RooUnfoldResponse *mg_pt_res_TU = (RooUnfoldResponse*) match_in->Get("mg_pt_res_TU");
  RooUnfoldResponse *mg_pt_res_HC50 = (RooUnfoldResponse*) match_in->Get("mg_pt_res_HC50");
  RooUnfoldResponse *mg_pt_res_DS = (RooUnfoldResponse*) match_in->Get("mg_pt_res_DS");
  RooUnfoldResponse *mg_pt_res_GS = (RooUnfoldResponse*) match_in->Get("mg_pt_res_GS");
  RooUnfoldResponse *mg_pt_res_MS = (RooUnfoldResponse*) match_in->Get("mg_pt_res_MS");

  //~~~data
  TH2D *m_v_pt_d = (TH2D*) data_in->Get("m_v_pt");
  TH2D *m_v_pt_g = (TH2D*) sim_in->Get("m_v_pt");
  TH2D *m_v_pt_p = (TH2D*) sim_in->Get("PL_m_v_pt");
  TH2D *m_v_pt_d_counts = (TH2D*) data_in->Get("m_v_pt_counts");
  TH2D *m_v_pt_g_counts = (TH2D*) sim_in->Get("m_v_pt_counts");
  TH2D *m_v_pt_p_counts = (TH2D*) sim_in->Get("PL_m_v_pt_counts");
  TH2D *mg_v_pt_d = (TH2D*) data_in->Get("mg_v_pt");
  TH2D *mg_v_pt_g = (TH2D*) sim_in->Get("mg_v_pt");
  TH2D *mg_v_pt_p = (TH2D*) sim_in->Get("PL_mg_v_pt");
  TH2D *mg_v_pt_d_counts = (TH2D*) data_in->Get("mg_v_pt_counts");
  TH2D *mg_v_pt_g_counts = (TH2D*) sim_in->Get("mg_v_pt_counts");
  TH2D *mg_v_pt_p_counts = (TH2D*) sim_in->Get("PL_mg_v_pt_counts");

  
  vector<int> to_drop_m; vector<int> to_drop_mg;
  vector<int> to_drop_m1D; vector<int> to_drop_mg1D; vector<int> to_drop_pt1D;
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//  
  //dropping low-stats closure spectra bins
  //pt                                                        
  DropLowStatsBins(sampleA_pt_gen, sampleA_pt_gen_counts);                                                                                
  to_drop_pt1D = DropLowStatsBins(sampleA_pt_det, sampleA_pt_det_counts);                                                                                
  DropLowStatsBins(sampleB_pt_gen, sampleB_pt_gen_counts);                                                                              
  DropLowStatsBins(sampleB_pt_det, sampleB_pt_det_counts);        
  //m
  DropLowStatsBins(sampleA_m_gen, sampleA_m_gen_counts);                                                                                
  to_drop_m1D = DropLowStatsBins(sampleA_m_det, sampleA_m_det_counts);                                                                                
  DropLowStatsBins(sampleB_m_gen, sampleB_m_gen_counts);                                                                              
  DropLowStatsBins(sampleB_m_det, sampleB_m_det_counts);        
  //m
  DropLowStatsBins(sampleA_mg_gen, sampleA_mg_gen_counts);                                                                                
  to_drop_mg1D = DropLowStatsBins(sampleA_mg_det, sampleA_mg_det_counts);                                                                                
  DropLowStatsBins(sampleB_mg_gen, sampleB_mg_gen_counts);                                                                              
  DropLowStatsBins(sampleB_mg_det, sampleB_mg_det_counts);        
  //m-pt
  DropLowStatsBins(sampleA_m_pt_gen, sampleA_m_pt_gen_counts);                                                                            
  DropLowStatsBins(sampleA_m_pt_det, sampleA_m_pt_det_counts);                                                                            
  DropLowStatsBins(sampleB_m_pt_gen, sampleB_m_pt_gen_counts);                                                                         
  DropLowStatsBins(sampleB_m_pt_det, sampleB_m_pt_det_counts); 
  //mg-pt
  DropLowStatsBins(sampleA_mg_pt_gen, sampleA_mg_pt_gen_counts);                                                                            
  DropLowStatsBins(sampleA_mg_pt_det, sampleA_mg_pt_det_counts);                                                                            
  DropLowStatsBins(sampleB_mg_pt_gen, sampleB_mg_pt_gen_counts);                                                                          
  DropLowStatsBins(sampleB_mg_pt_det, sampleB_mg_pt_det_counts); 
  //~~~
  
  //dropping low-stats spectra bins
  DropLowStatsBins(m_v_pt_d, m_v_pt_d_counts);                                                                             
  DropLowStatsBins(m_v_pt_g, m_v_pt_g_counts);                                                                             
  DropLowStatsBins(m_v_pt_p, m_v_pt_p_counts);
  DropLowStatsBins(mg_v_pt_d, mg_v_pt_d_counts);                                                                             
  DropLowStatsBins(mg_v_pt_g, mg_v_pt_g_counts);                                                                             
  DropLowStatsBins(mg_v_pt_p, mg_v_pt_p_counts);
  //~~~


  //dropping bins from the closure responses using the closure spectra
  //remember to call this block of after spectra have been pruned already
  to_drop_m = DropLowStatsBins(sampleA_m_pt_response, sampleA_m_pt_det_counts);
  to_drop_mg = DropLowStatsBins(sampleA_mg_pt_response, sampleA_mg_pt_det_counts);
  //we determined the bins to drop using sampleA, now we drop the same bins from sampleB (should be roughly the same statistics)
  DropBins(sampleB_m_pt_response, to_drop_m);
  DropBins(sampleB_mg_pt_response, to_drop_mg);  
  //1Ds:
  DropBins(sampleA_pt_response, to_drop_pt1D);
  DropBins(sampleA_m_response, to_drop_m1D);
  DropBins(sampleA_mg_response, to_drop_mg1D);
  DropBins(sampleB_pt_response, to_drop_pt1D);
  DropBins(sampleB_m_response, to_drop_m1D);
  DropBins(sampleB_mg_response, to_drop_mg1D);
  //~~~

  //dropping bins from the responses using the closure spectra
  DropBins(m_pt_response, to_drop_m);   
  DropBins(mg_pt_response, to_drop_mg);
  //~~~

  //dropping bins from the systematics responses using the closure spectra
  DropBins(m_pt_res_nom,to_drop_m);
  DropBins(mg_pt_res_nom,to_drop_mg);
  DropBins(m_pt_res_TS,to_drop_m);
  DropBins(mg_pt_res_TS,to_drop_mg);
  DropBins(m_pt_res_TU,to_drop_m);
  DropBins(mg_pt_res_TU,to_drop_mg); 
  DropBins(m_pt_res_HC50,to_drop_m);
  DropBins(mg_pt_res_HC50,to_drop_mg);
  DropBins(m_pt_res_DS,to_drop_m);
  DropBins(mg_pt_res_DS,to_drop_mg);
  DropBins(m_pt_res_GS,to_drop_m);
  DropBins(mg_pt_res_GS,to_drop_mg);
  DropBins(m_pt_res_MS,to_drop_m);
  DropBins(mg_pt_res_MS,to_drop_mg);   
  
  //~~~

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

  //output files - same as inputs but with "_bindropped" appended.
  TFile *match_out = new TFile((match_path+match_file+"_bindropped.root").c_str(),"RECREATE");
  TFile *data_out = new TFile((data_path+data_file+"_bindropped.root").c_str(),"RECREATE");
  TFile *sim_out = new TFile((sim_path+sim_file+"_bindropped.root").c_str(),"RECREATE");

  match_out->cd();
  //nominal responses
  m_pt_response->Write(); mg_pt_response->Write();
  //responses for closure test
  sampleA_pt_response->Write(); sampleB_pt_response->Write();
  sampleA_m_response->Write(); sampleB_m_response->Write();
  sampleA_mg_response->Write(); sampleB_mg_response->Write();
  sampleA_m_pt_response->Write(); sampleA_mg_pt_response->Write();
  sampleB_m_pt_response->Write(); sampleB_mg_pt_response->Write();
  sampleA_pt_gen->Write(); sampleB_pt_gen->Write(); sampleA_pt_det->Write(); sampleB_pt_det->Write();
  sampleA_m_gen->Write(); sampleB_m_gen->Write(); sampleA_m_det->Write(); sampleB_m_det->Write();
  sampleA_mg_gen->Write(); sampleB_mg_gen->Write(); sampleA_mg_det->Write(); sampleB_mg_det->Write();
  sampleA_m_pt_gen->Write(); sampleB_m_pt_gen->Write(); sampleA_m_pt_det->Write(); sampleB_m_pt_det->Write();
  sampleA_mg_pt_gen->Write(); sampleB_mg_pt_gen->Write(); sampleA_mg_pt_det->Write(); sampleB_mg_pt_det->Write();
  //responses for systematic uncertainty variation
  m_pt_res_nom->Write(); mg_pt_res_nom->Write();
  m_pt_res_TS->Write(); mg_pt_res_TS->Write();
  m_pt_res_TU->Write(); mg_pt_res_TU->Write();
  m_pt_res_HC50->Write(); mg_pt_res_HC50->Write();
  m_pt_res_DS->Write(); mg_pt_res_DS->Write();
  m_pt_res_GS->Write(); mg_pt_res_GS->Write();
  m_pt_res_MS->Write(); mg_pt_res_MS->Write();

  data_out->cd();
  //data spectra
  m_v_pt_d->Write();
  mg_v_pt_d->Write();

  sim_out->cd();
  //(unmatched) simulation spectra
  m_v_pt_g->Write(); m_v_pt_p->Write();
  mg_v_pt_g->Write(); mg_v_pt_p->Write();
  
  cout << endl
       << "Wrote data to " << data_out->GetName() << endl
       << "Wrote sim spectra to " << sim_out->GetName() << endl
       << "Wrote sim responses to " << match_out->GetName() << endl
       << "Wrote closure to " << match_out->GetName() << endl
       << "Wrote systematics responses to " << match_out->GetName() << endl;

  
  cout << endl << "Dropped mass bins: " << endl;
  for (unsigned i = 0; i < to_drop_m.size() - 1; ++ i) {
    cout << to_drop_m[i] << ", ";
  }
  cout << to_drop_m[to_drop_m.size() - 1] << endl;
  
  cout << "Dropped Mg bins: " << endl;
  for (unsigned i = 0; i < to_drop_mg.size() - 1; ++ i) {
    cout << to_drop_mg[i] << ", ";
  }
  cout << to_drop_mg[to_drop_mg.size() - 1] << endl;


  //closing files
  data_out->Close(); sim_out->Close(); match_out->Close();
  cout << endl << "Closed \t" << data_out->GetName() << ",\n\t" << sim_out->GetName() 
       << ",\n\t" << match_out->GetName() << endl;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//      
  
  return 0;
}
