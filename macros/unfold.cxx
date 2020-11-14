//Isaac Mooney, WSU, August 2019
//This file takes in the data and systematics responses, RooUnfolds using the data and the nominal and varied responses, appropriately adjusts the statistical errors using an input file (produced/described by macros/stat_err_scaling.cxx), calculates the systematic uncertainties on the mass result, and saves the result with these uncertainties as output.

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

//THESE FUNCTIONS ARE HERE ONLY TEMPORARILY UNTIL I HAVE TIME TO LINK THE PLOTS LIBRARY/MAKE A NEW ONE

//projects a 2D histogram in desired ranges and returns an array of the (1D) projections on desired axis                                                        
std::vector<TH1D*> Projection2D (TH2D * hist2D, const int nBins, double * ranges, const std::string axis) {
  std::vector<TH1D*> proj1Ds;
  for (int i = 0; i < nBins; ++ i) {
    std::string low = std::to_string(ranges[i]);
    std::string high = std::to_string(ranges[i+1]);
    std::string low_rough = low.substr(0,2);
    std::string high_rough = high.substr(0,2);
    if (low_rough.substr(1,2) == ".") {low_rough = low_rough.substr(0,1);}
    if (high_rough.substr(1,2) == ".") {high_rough = high_rough.substr(0,1);}
    if (axis == "x" || axis == "X" || axis == "1") {
      proj1Ds.push_back(hist2D->ProjectionX((hist2D->GetName() + axis + low_rough + high_rough).c_str(),ranges[i],ranges[i+1]- 1));
    }
    else if (axis == "y" || axis == "Y" || axis == "2") {
      proj1Ds.push_back(hist2D->ProjectionY((hist2D->GetName() + axis + low_rough + high_rough).c_str(),ranges[i],ranges[i+1] - 1));
    }
    else {
      std::cerr << "Improper axis given for projections. Exiting." << std::endl; exit(1);
    }
    proj1Ds[i]->SetTitle("");
  }
  return proj1Ds;
}

//prettify 1D histogram                                                                                                                                        
void Prettify1D (TH1D * hist, const Color_t markColor, const Style_t markStyle, const double markSize, const Color_t lineColor,
                 const std::string xTitle, const std::string yTitle, const double lowx, const double highx, const double lowy, const double highy) {
  if (yTitle.find("1/N") != std::string::npos && yTitle.find("N_{cons}") == std::string::npos) {
    hist->Scale(1/(double)hist->Integral());
    double binwidth = (hist->GetXaxis()->GetXmax() - hist->GetXaxis()->GetXmin()) / (double) hist->GetXaxis()->GetNbins();
    hist->Scale(1/(double)binwidth);
  }
  else if (yTitle.find("1/N") != std::string::npos && yTitle.find("N_{cons}") != std::string::npos) {
    double binwidth = (hist->GetXaxis()->GetXmax() - hist->GetXaxis()->GetXmin()) / (double) hist->GetXaxis()->GetNbins();
    hist->Scale(1/(double)binwidth);
  }
  else {
    cout << "Not scaling by histogram " << hist->GetName() << "'s integral & bin width! If this is a cross section measurement, this will be a problem! Fix by passing histogram y-axis title containing anything but 'arb'" << endl;
  }

  hist->SetMarkerColor(markColor); hist->SetMarkerStyle(markStyle); hist->SetMarkerSize(markSize); hist->SetLineColor(lineColor);
  hist->GetXaxis()->SetTitle((xTitle).c_str()); hist->GetYaxis()->SetTitle((yTitle).c_str());
  if (highx != -1) {
    hist->GetXaxis()->SetRangeUser(lowx, highx);
  }
  if (highy != -1) {
    hist->GetYaxis()->SetRangeUser(lowy, highy);
  }
  return;
}

//prettify 1D histogram to be drawn as line with "C" option                                                                                                     
void Prettify1DwLineStyle(TH1D * hist, const Color_t lineColor, const Style_t lineStyle, const double lineWidth,
                          const std::string xTitle, const std::string yTitle, const double lowx, const double highx, const double lowy, const double highy) {
  if (yTitle.find("1/N") != std::string::npos && yTitle.find("N_{cons}") == std::string::npos) { 
    hist->Scale(1/(double)hist->Integral());
    double binwidth = (hist->GetXaxis()->GetXmax() - hist->GetXaxis()->GetXmin()) / (double) hist->GetXaxis()->GetNbins();
    hist->Scale(1/(double)binwidth);
  }
  else if (yTitle.find("1/N") != std::string::npos && yTitle.find("N_{cons}") != std::string::npos) {
    double binwidth = (hist->GetXaxis()->GetXmax() - hist->GetXaxis()->GetXmin()) / (double) hist->GetXaxis()->GetNbins();
    hist->Scale(1/(double)binwidth);
  }
  else {
    cout << "Not scaling by histogram " << hist->GetName() << "'s integral & bin width! If this is a cross section measurement, this will be a problem! Fix by passing histogram y-axis title containing anything but 'arb'" << endl;
  }

  hist->SetLineColor(lineColor); hist->SetLineStyle(lineStyle); hist->SetLineWidth(lineWidth);
  hist->GetXaxis()->SetTitle((xTitle).c_str()); hist->GetYaxis()->SetTitle((yTitle).c_str());
  if (highx != -1) {
    hist->GetXaxis()->SetRangeUser(lowx, highx);
  }
  if (highy != -1) {
    hist->GetYaxis()->SetRangeUser(lowy, highy);
  }
  hist->Sumw2(0);
  return;
}

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
  const string data_path = "~/jetmass2/out/data/hists/";
  const string match_file = "FINAL_sim_matched";
  const string data_file = "FINAL_data_hists_ppJP2";
  
  //input files
  TFile *fres = new TFile((match_path+match_file+radius+"_paper_bindropped.root").c_str(),"READ");
  TFile *fdat = new TFile((data_path+data_file+radius+"_bindropped.root").c_str(),"READ");
  cout << "A" << endl;
  //systematics responses (rnom is nominal response for unfolding)
  RooUnfoldResponse *rnom = (RooUnfoldResponse*) fres->Get((hname+"_pt_res_nom").c_str());
  RooUnfoldResponse *rTS = (RooUnfoldResponse*) fres->Get((hname+"_pt_res_TS").c_str()); 
  RooUnfoldResponse *rTU = (RooUnfoldResponse*) fres->Get((hname+"_pt_res_TU").c_str()); 
  RooUnfoldResponse *rHC50 = (RooUnfoldResponse*) fres->Get((hname+"_pt_res_HC50").c_str()); 
  RooUnfoldResponse *rDS = (RooUnfoldResponse*) fres->Get((hname+"_pt_res_DS").c_str()); 
  RooUnfoldResponse *rGS = (RooUnfoldResponse*) fres->Get((hname+"_pt_res_GS").c_str()); 
  //  RooUnfoldResponse *rMS = (RooUnfoldResponse*) fres->Get((hname+"_pt_res_MS").c_str());
  RooUnfoldResponse *rMS_2025_nom = (RooUnfoldResponse*) fres->Get((hname+"_res2025_nom").c_str());
  RooUnfoldResponse *rMS_2025_h7smear = (RooUnfoldResponse*) fres->Get((hname+"_res2025_h7smear").c_str());
  RooUnfoldResponse *rMS_2025_p8smear = (RooUnfoldResponse*) fres->Get((hname+"_res2025_p8smear").c_str());
  RooUnfoldResponse *rMS_2530_nom = (RooUnfoldResponse*) fres->Get((hname+"_res2530_nom").c_str());
  RooUnfoldResponse *rMS_2530_h7smear = (RooUnfoldResponse*) fres->Get((hname+"_res2530_h7smear").c_str());
  RooUnfoldResponse *rMS_2530_p8smear = (RooUnfoldResponse*) fres->Get((hname+"_res2530_p8smear").c_str());
  RooUnfoldResponse *rMS_3040_nom = (RooUnfoldResponse*) fres->Get((hname+"_res3040_nom").c_str());
  RooUnfoldResponse *rMS_3040_h7smear = (RooUnfoldResponse*) fres->Get((hname+"_res3040_h7smear").c_str());
  RooUnfoldResponse *rMS_3040_p8smear = (RooUnfoldResponse*) fres->Get((hname+"_res3040_p8smear").c_str());
  
  //data spectrum
  TH2D* m_pt_dat = (TH2D*) fdat->Get((hname+"_v_pt").c_str());
  TH1D* m_dat_2025 = (TH1D*) fdat->Get((hname+"_2025_d").c_str());
  TH1D* m_dat_2530 = (TH1D*) fdat->Get((hname+"_2530_d").c_str());
  TH1D* m_dat_3040 = (TH1D*) fdat->Get((hname+"_3040_d").c_str());
  
  cout << "B" << endl;
  
  RooUnfoldBayes *unfold_nom = new RooUnfoldBayes(rnom, m_pt_dat, 4, false, "unfold_nom","");
  RooUnfoldBayes *unfold_IP2 = new RooUnfoldBayes(rnom, m_pt_dat, 2, false, "unfold_IP2","");
  RooUnfoldBayes *unfold_IP6 = new RooUnfoldBayes(rnom, m_pt_dat, 6, false, "unfold_IP6","");
  RooUnfoldBayes *unfold_TS = new RooUnfoldBayes(rTS, m_pt_dat, 4, false, "unfold_TS","");
  RooUnfoldBayes *unfold_TU = new RooUnfoldBayes(rTU, m_pt_dat, 4, false, "unfold_TU","");
  RooUnfoldBayes *unfold_HC50 = new RooUnfoldBayes(rHC50, m_pt_dat, 4, false, "unfold_HC50","");
  RooUnfoldBayes *unfold_DS = new RooUnfoldBayes(rDS, m_pt_dat, 4, false, "unfold_DS","");
  RooUnfoldBayes *unfold_GS = new RooUnfoldBayes(rGS, m_pt_dat, 4, false, "unfold_GS","");
  //RooUnfoldBayes *unfold_MS = new RooUnfoldBayes(rMS, m_pt_dat, 4, false, "unfold_MS","");
  RooUnfoldBayes *unfold_nom1D_2025 = new RooUnfoldBayes(rMS_2025_nom,m_dat_2025,4,false,"unfold_nom1D_2025","");
  RooUnfoldBayes *unfold_nom1D_2530 = new RooUnfoldBayes(rMS_2530_nom,m_dat_2530,4,false,"unfold_nom1D_2530",""); 
  RooUnfoldBayes *unfold_nom1D_3040 = new RooUnfoldBayes(rMS_3040_nom,m_dat_3040,4,false,"unfold_nom1D_3040","");
  RooUnfoldBayes *unfold_h7smear1D_2025 = new RooUnfoldBayes(rMS_2025_h7smear,m_dat_2025,4,false,"unfold_h7smear1D_2025","");
  RooUnfoldBayes *unfold_h7smear1D_2530 = new RooUnfoldBayes(rMS_2530_h7smear,m_dat_2530,4,false,"unfold_h7smear1D_2530","");
  RooUnfoldBayes *unfold_h7smear1D_3040 = new RooUnfoldBayes(rMS_3040_h7smear,m_dat_3040,4,false,"unfold_h7smear1D_3040","");
  RooUnfoldBayes *unfold_p8smear1D_2025 = new RooUnfoldBayes(rMS_2025_p8smear,m_dat_2025,4,false,"unfold_p8smear1D_2025","");
  RooUnfoldBayes *unfold_p8smear1D_2530 = new RooUnfoldBayes(rMS_2530_p8smear,m_dat_2530,4,false,"unfold_p8smear1D_2530","");
  RooUnfoldBayes *unfold_p8smear1D_3040 = new RooUnfoldBayes(rMS_3040_p8smear,m_dat_3040,4,false,"unfold_p8smear1D_3040","");


  cout << "C" << endl;
  TH2D *reco_nom = (TH2D*) unfold_nom->Hreco((RooUnfold::ErrorTreatment) 3);
  TH2D *reco_IP2 = (TH2D*) unfold_IP2->Hreco((RooUnfold::ErrorTreatment) 3);
  TH2D *reco_IP6 = (TH2D*) unfold_IP6->Hreco((RooUnfold::ErrorTreatment) 3);
  TH2D *reco_TS = (TH2D*) unfold_TS->Hreco((RooUnfold::ErrorTreatment) 3);
  TH2D *reco_TU = (TH2D*) unfold_TU->Hreco((RooUnfold::ErrorTreatment) 3);
  TH2D *reco_HC50 = (TH2D*) unfold_HC50->Hreco((RooUnfold::ErrorTreatment) 3);
  TH2D *reco_DS = (TH2D*) unfold_DS->Hreco((RooUnfold::ErrorTreatment) 3);
  TH2D *reco_GS = (TH2D*) unfold_GS->Hreco((RooUnfold::ErrorTreatment) 3);
  //TH2D *reco_MS = (TH2D*) unfold_MS->Hreco((RooUnfold::ErrorTreatment) 3);
  TH1D *reco_nom1D_2025 = (TH1D*) unfold_nom1D_2025->Hreco((RooUnfold::ErrorTreatment) 3);
  TH1D *reco_nom1D_2530 = (TH1D*) unfold_nom1D_2530->Hreco((RooUnfold::ErrorTreatment) 3);
  TH1D *reco_nom1D_3040 = (TH1D*) unfold_nom1D_3040->Hreco((RooUnfold::ErrorTreatment) 3);
  TH1D *reco_h7smear1D_2025 = (TH1D*) unfold_h7smear1D_2025->Hreco((RooUnfold::ErrorTreatment) 3);
  TH1D *reco_h7smear1D_2530 = (TH1D*) unfold_h7smear1D_2530->Hreco((RooUnfold::ErrorTreatment) 3);
  TH1D *reco_h7smear1D_3040 = (TH1D*) unfold_h7smear1D_3040->Hreco((RooUnfold::ErrorTreatment) 3);
  TH1D *reco_p8smear1D_2025 = (TH1D*) unfold_p8smear1D_2025->Hreco((RooUnfold::ErrorTreatment) 3);
  TH1D *reco_p8smear1D_2530 = (TH1D*) unfold_p8smear1D_2530->Hreco((RooUnfold::ErrorTreatment) 3);
  TH1D *reco_p8smear1D_3040 = (TH1D*) unfold_p8smear1D_3040->Hreco((RooUnfold::ErrorTreatment) 3);
  
  //These ranges are okay now that Projection2D has the plotting bug removed.
  const int nBins = 3; //change back to 3 if you show a different range
  const int nBins_m = 14;
  TAxis* reco_axis = reco_nom->GetYaxis(); TAxis* det_axis = m_pt_dat->GetYaxis();
  
  double ranges[nBins + 1] = {(double) reco_axis->FindBin(20), (double) reco_axis->FindBin(25), (double) reco_axis->FindBin(30), (double) reco_axis->FindBin(40)};
  double ranges_d[nBins + 1] = {(double) det_axis->FindBin(20), (double) det_axis->FindBin(25), (double) det_axis->FindBin(30), (double) det_axis->FindBin(40)};
  string pts[nBins + 1] = {"20","25","30","40"};
  
  vector<TH1D*> reco_noms = Projection2D(reco_nom,nBins,ranges,"x");
  vector<TH1D*> reco_IP2s = Projection2D(reco_IP2,nBins,ranges,"x");
  vector<TH1D*> reco_IP6s = Projection2D(reco_IP6,nBins,ranges,"x");
  vector<TH1D*> reco_TSs = Projection2D(reco_TS,nBins,ranges,"x");
  vector<TH1D*> reco_TUs = Projection2D(reco_TU,nBins,ranges,"x");
  vector<TH1D*> reco_HC50s = Projection2D(reco_HC50,nBins,ranges,"x");
  vector<TH1D*> reco_DSs = Projection2D(reco_DS,nBins,ranges,"x");
  vector<TH1D*> reco_GSs = Projection2D(reco_GS,nBins,ranges,"x");
  //vector<TH1D*> reco_MSs = Projection2D(reco_MS,nBins,ranges,"x");

  vector<TH1D*> reco_noms_copy;
  vector<TH1D*> reco_IP2s_copy;
  vector<TH1D*> reco_IP6s_copy;
  vector<TH1D*> reco_TSs_copy;
  vector<TH1D*> reco_TUs_copy;
  vector<TH1D*> reco_HC50s_copy;
  vector<TH1D*> reco_DSs_copy;
  vector<TH1D*> reco_GSs_copy;
  //vector<TH1D*> reco_MSs_copy;

  
  for (int i = 0; i < nBins; ++ i) {
    reco_noms_copy.push_back((TH1D*) reco_noms[i]->Clone(("nom_"+to_string(i)).c_str()));
    reco_IP2s_copy.push_back((TH1D*) reco_IP2s[i]->Clone(("IP2_"+to_string(i)).c_str()));
    reco_IP6s_copy.push_back((TH1D*) reco_IP6s[i]->Clone(("IP6_"+to_string(i)).c_str()));
    reco_TSs_copy.push_back((TH1D*) reco_TSs[i]->Clone(("TS_"+to_string(i)).c_str()));
    reco_TUs_copy.push_back((TH1D*) reco_TUs[i]->Clone(("TU_"+to_string(i)).c_str()));
    reco_HC50s_copy.push_back((TH1D*) reco_HC50s[i]->Clone(("HC50_"+to_string(i)).c_str()));
    reco_DSs_copy.push_back((TH1D*) reco_DSs[i]->Clone(("DS_"+to_string(i)).c_str()));
    reco_GSs_copy.push_back((TH1D*) reco_GSs[i]->Clone(("GS_"+to_string(i)).c_str()));
    //reco_MSs_copy.push_back((TH1D*) reco_MSs[i]->Clone(("MS_"+to_string(i)).c_str()));

  }

  for (int i = 0; i < nBins; ++ i) {
    reco_noms[i]->Scale(1/(double)reco_noms[i]->Integral());
    reco_IP2s[i]->Scale(1/(double)reco_IP2s[i]->Integral());
    reco_IP6s[i]->Scale(1/(double)reco_IP6s[i]->Integral());
    reco_TSs[i]->Scale(1/(double)reco_TSs[i]->Integral());
    reco_TUs[i]->Scale(1/(double)reco_TUs[i]->Integral());
    reco_HC50s[i]->Scale(1/(double)reco_HC50s[i]->Integral());
    reco_DSs[i]->Scale(1/(double)reco_DSs[i]->Integral());
    reco_GSs[i]->Scale(1/(double)reco_GSs[i]->Integral());
    //reco_MSs[i]->Scale(1/(double)reco_MSs[i]->Integral());
    
    reco_IP2s[i]->Divide(reco_noms[i]);
    reco_IP6s[i]->Divide(reco_noms[i]);
    reco_TSs[i]->Divide(reco_noms[i]);
    reco_TUs[i]->Divide(reco_noms[i]);
    reco_HC50s[i]->Divide(reco_noms[i]);
    reco_DSs[i]->Divide(reco_noms[i]);
    reco_GSs[i]->Divide(reco_noms[i]);
    //reco_MSs[i]->Divide(reco_noms[i]);
  }
  
  //1D mass smearing systematic ratio and bookkeeping
  
  reco_nom1D_2025->Scale(1/(double)reco_nom1D_2025->Integral());
  reco_nom1D_2530->Scale(1/(double)reco_nom1D_2530->Integral());
  reco_nom1D_3040->Scale(1/(double)reco_nom1D_3040->Integral());
  reco_h7smear1D_2025->Scale(1/(double)reco_h7smear1D_2025->Integral());
  reco_h7smear1D_2530->Scale(1/(double)reco_h7smear1D_2530->Integral());
  reco_h7smear1D_3040->Scale(1/(double)reco_h7smear1D_3040->Integral());
  reco_p8smear1D_2025->Scale(1/(double)reco_p8smear1D_2025->Integral());
  reco_p8smear1D_2530->Scale(1/(double)reco_p8smear1D_2530->Integral());
  reco_p8smear1D_3040->Scale(1/(double)reco_p8smear1D_3040->Integral());

  reco_h7smear1D_2025->Divide(reco_nom1D_2025);
  reco_p8smear1D_2025->Divide(reco_nom1D_2025);
  reco_h7smear1D_2530->Divide(reco_nom1D_2530);
  reco_p8smear1D_2530->Divide(reco_nom1D_2530);
  reco_h7smear1D_3040->Divide(reco_nom1D_3040);
  reco_p8smear1D_3040->Divide(reco_nom1D_3040);

  Double_t stats[5] ={0,0,0,0,0};
  
  reco_h7smear1D_2025->PutStats(stats);
  reco_p8smear1D_2025->PutStats(stats);
  reco_h7smear1D_2530->PutStats(stats);
  reco_p8smear1D_2530->PutStats(stats);
  reco_h7smear1D_3040->PutStats(stats);
  reco_p8smear1D_3040->PutStats(stats);

  reco_h7smear1D_2025->Sumw2(0);
  reco_p8smear1D_2025->Sumw2(0);
  reco_h7smear1D_2530->Sumw2(0);
  reco_p8smear1D_2530->Sumw2(0);
  reco_h7smear1D_3040->Sumw2(0);
  reco_p8smear1D_3040->Sumw2(0);
  
  for (int j = 1; j <= nBins_m; ++ j) {
    reco_h7smear1D_2025->SetBinContent(j,fabs(reco_h7smear1D_2025->GetBinContent(j) - 1));
    reco_p8smear1D_2025->SetBinContent(j,fabs(reco_p8smear1D_2025->GetBinContent(j) - 1));
    reco_h7smear1D_2530->SetBinContent(j,fabs(reco_h7smear1D_2530->GetBinContent(j) - 1));
    reco_p8smear1D_2530->SetBinContent(j,fabs(reco_p8smear1D_2530->GetBinContent(j) - 1));
    reco_h7smear1D_3040->SetBinContent(j,fabs(reco_h7smear1D_3040->GetBinContent(j) - 1));
    reco_p8smear1D_3040->SetBinContent(j,fabs(reco_p8smear1D_3040->GetBinContent(j) - 1));
  }
  
  for (int i = 0; i < nBins; ++ i) {
    reco_IP2s[i]->PutStats(stats);
    reco_IP6s[i]->PutStats(stats);
    reco_TSs[i]->PutStats(stats);
    reco_TUs[i]->PutStats(stats);
    reco_HC50s[i]->PutStats(stats);
    reco_DSs[i]->PutStats(stats);
    reco_GSs[i]->PutStats(stats);
    //reco_MSs[i]->PutStats(stats);
    
    reco_IP2s[i]->Sumw2(0);
    reco_IP6s[i]->Sumw2(0);
    reco_TSs[i]->Sumw2(0);
    reco_TUs[i]->Sumw2(0);
    reco_HC50s[i]->Sumw2(0);
    reco_DSs[i]->Sumw2(0);
    reco_GSs[i]->Sumw2(0);
    //reco_MSs[i]->Sumw2(0);
    
    
    for (int j = 1; j <= nBins_m; ++ j) {
      //turning the ratio into a percentage.
      reco_IP2s[i]->SetBinContent(j,fabs(reco_IP2s[i]->GetBinContent(j) - 1));
      reco_IP6s[i]->SetBinContent(j,fabs(reco_IP6s[i]->GetBinContent(j) - 1));
      reco_TSs[i]->SetBinContent(j,fabs(reco_TSs[i]->GetBinContent(j) - 1));
      reco_TUs[i]->SetBinContent(j,fabs(reco_TUs[i]->GetBinContent(j) - 1));
      reco_HC50s[i]->SetBinContent(j,fabs(reco_HC50s[i]->GetBinContent(j) - 1));
      reco_DSs[i]->SetBinContent(j,fabs(reco_DSs[i]->GetBinContent(j) - 1));
      reco_GSs[i]->SetBinContent(j,fabs(reco_GSs[i]->GetBinContent(j) - 1));
      //reco_MSs[i]->SetBinContent(j,fabs(reco_MSs[i]->GetBinContent(j) - 1));
    }
    
  }
  
  cout << "D" << endl;
  
  for (int i = 0; i < nBins; ++ i) {
    Prettify1DwLineStyle(reco_IP2s[i],2, kSolid, 2,(htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1); 
    reco_IP2s[i]->SetFillColor(2); reco_IP2s[i]->SetFillStyle(3305);
    Prettify1DwLineStyle(reco_IP6s[i],3, kSolid, 2,(htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1);
    reco_IP6s[i]->SetFillColor(3); reco_IP6s[i]->SetFillStyle(3395);
    Prettify1DwLineStyle(reco_TSs[i],4, kSolid, 2,(htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1);
    reco_TSs[i]->SetFillColor(4); reco_TSs[i]->SetFillStyle(3490);
    Prettify1DwLineStyle(reco_TUs[i],5, kSolid, 2,(htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1);
    reco_TUs[i]->SetFillColor(5); reco_TUs[i]->SetFillStyle(3436);
    Prettify1DwLineStyle(reco_HC50s[i],6, kSolid, 2,(htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1);
    reco_HC50s[i]->SetFillColor(6); reco_HC50s[i]->SetFillStyle(3335);
    Prettify1DwLineStyle(reco_DSs[i],8, kSolid, 2,(htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1); 
    reco_DSs[i]->SetFillColor(8); reco_DSs[i]->SetFillStyle(3944);
    Prettify1DwLineStyle(reco_GSs[i],9, kSolid, 2,(htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1);
    reco_GSs[i]->SetFillColor(9); reco_GSs[i]->SetFillStyle(3544);
    //  Prettify1DwLineStyle(reco_MSs[i],11, kSolid, 2,(htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1);
    //reco_MSs[i]->SetFillColor(11); reco_MSs[i]->SetFillStyle(3690);
  }
  Prettify1DwLineStyle(reco_h7smear1D_2025, 11, kSolid, 2,(htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1);
  reco_h7smear1D_2025->SetFillColor(11); reco_h7smear1D_2025->SetFillStyle(3690);
  Prettify1DwLineStyle(reco_h7smear1D_2530, 11, kSolid, 2,(htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1);
  reco_h7smear1D_2530->SetFillColor(11); reco_h7smear1D_2530->SetFillStyle(3690);
  Prettify1DwLineStyle(reco_h7smear1D_3040, 11, kSolid, 2,(htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1);
  reco_h7smear1D_3040->SetFillColor(11); reco_h7smear1D_3040->SetFillStyle(3690);
 
  Prettify1DwLineStyle(reco_p8smear1D_2025, 11, kSolid, 2,(htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1);
  reco_p8smear1D_2025->SetFillColor(11); reco_p8smear1D_2025->SetFillStyle(3444);
  Prettify1DwLineStyle(reco_p8smear1D_2530, 11, kSolid, 2,(htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1);
  reco_p8smear1D_2530->SetFillColor(11); reco_p8smear1D_2530->SetFillStyle(3444);
  Prettify1DwLineStyle(reco_p8smear1D_3040, 11, kSolid, 2,(htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1);
  reco_p8smear1D_3040->SetFillColor(11); reco_p8smear1D_3040->SetFillStyle(3444);
  

  //taking maximum envelopes!
  TH2D* env_HC = new TH2D("env_HC","",14,0,14,15,5,80);
  vector<TH1D*> env_HCs = Projection2D(env_HC,nBins,ranges,"x");
  TH2D* env_un = new TH2D("env_un","",14,0,14,15,5,80);
  vector<TH1D*> env_uns = Projection2D(env_un,nBins,ranges,"x");
  TH2D* net = new TH2D("net","",14,0,14,15,5,80);
  vector<TH1D*> nets = Projection2D(net,nBins,ranges,"x");

  vector<vector< double> > syst_errs2D;
  
  for (int i = 0; i < nBins; ++ i) {
    vector<double> syst_errs1D; 
    reco_noms_copy[i]->Scale(1/(double)reco_noms_copy[i]->Integral());
    //for each mass bin, determine the largest effect in each category
    for (int j = 1; j < nBins_m + 1; ++ j) {
      //hadronic correction envelope
      double hcs [1] = {reco_HC50s[i]->GetBinContent(j)};
      set<double> hc_sort (hcs, hcs+1);
      set<double>::iterator hc = hc_sort.end(); hc --;
      double hc_envelope = *hc;
      env_HCs[i]->SetBinContent(j, hc_envelope);
      //a clunky way of doing this for both current ranges of 1D mass responses: 20-30, 30-45:
      TH1D* reco_h7_1D = new TH1D("reco_h7_1D","",14,0,14);
      TH1D* reco_p8_1D = new TH1D("reco_p8_1D","",14,0,14);
      if (i == 0) {reco_h7_1D = reco_h7smear1D_2025; reco_p8_1D = reco_p8smear1D_2025;}
      if (i == 1) {reco_h7_1D = reco_h7smear1D_2530; reco_p8_1D = reco_p8smear1D_2530;}
      if (i == 2) {reco_h7_1D = reco_h7smear1D_3040; reco_p8_1D = reco_p8smear1D_3040;}
      //unfolding envelope - using an ordered set here to automatically get the largest value
      double uns [6] = {reco_DSs[i]->GetBinContent(j), reco_GSs[i]->GetBinContent(j), /*reco_MSs[i]->GetBinContent(j),*/ reco_h7_1D->GetBinContent(j), reco_p8_1D->GetBinContent(j), reco_IP2s[i]->GetBinContent(j), reco_IP6s[i]->GetBinContent(j)};
      set<double> un_sort (uns, uns+6);
      set<double>::iterator un = un_sort.end(); un --;
      double un_envelope = *un;
      env_uns[i]->SetBinContent(j, un_envelope);
      //total uncertainty = TU + TS + un envelope + hc envelope
      double square = pow(hc_envelope,2) + pow(un_envelope,2) + pow(reco_TUs[i]->GetBinContent(j),2) + pow(reco_TSs[i]->GetBinContent(j),2);
      nets[i]->SetBinContent(j,sqrt(square));
      syst_errs1D.push_back(nets[i]->GetBinContent(j));
    }
    
    syst_errs2D.push_back(syst_errs1D);
  }

  //prettification
  for (int i = 0; i < nBins; ++ i) {
    Prettify1DwLineStyle(env_HCs[i],7, kSolid, 2,(htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1);
    env_HCs[i]->SetFillColor(7); env_HCs[i]->SetFillStyle(3353);
    Prettify1DwLineStyle(env_uns[i],2, kSolid, 2,(htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1);
    env_uns[i]->SetFillColor(2); env_uns[i]->SetFillStyle(3305);
    Prettify1DwLineStyle(nets[i], kBlack, kSolid, 2, (htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1);
  }

  cout << "E" << endl;
  
  //unfolded result with systematic errors!
  vector<TH1D*> w_systs;
  for (int i = 0; i < nBins; ++ i) { 
    w_systs.push_back((TH1D*) reco_noms_copy[i]->Clone(("w_systs_"+to_string(i)).c_str()));
    for (int j = 1; j < nBins_m + 1; ++ j) {
      w_systs[i]->SetBinError(j, syst_errs2D[i][j-1]*w_systs[i]->GetBinContent(j));
    }
  }

  vector<TH1D*> dats = Projection2D(m_pt_dat,nBins,ranges_d,"x");

  for (int i = 0; i < nBins; ++ i) {
    Prettify1D(reco_noms_copy[i],kRed,kFullStar,4,kRed,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,10,0,0.5);
    Prettify1D(dats[i], kBlack, kOpenStar, 4, kBlack, (htitle+" [GeV/c^{2}]").c_str(), ("1/N dN/d"+htitle).c_str(),0,10,0,0.5);
    Prettify1D(w_systs[i],kRed,kFullStar,0,kRed,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,10,0,0.5);
    w_systs[i]->SetFillColor(kRed - 10); w_systs[i]->SetFillStyle(1001);
  }
   
  //scaling errors
  TFile *fstats = new TFile(("~/jetmass2/out/sim/FINAL_stat_err_scaling"+radius+".root").c_str(),"READ");
  TH1D* scalefactors = (TH1D*) fstats->Get("hratio");

  for (int i = 0; i < nBins; ++ i) {
    cout << "i = " << i << endl << endl;
    for (int j = 1; j <= reco_noms[i]->GetNbinsX(); ++ j) {
      double scaling = -1;
      TAxis *scalex = scalefactors->GetXaxis();
      //will access the scalefactors histogram for pt values 20, 25, 30, 35, so if you change the ranges you show, change the FindBin calls, too.   
      if (i == 0) {scaling = scalefactors->GetBinContent(scalex->FindBin(20));/*old: 1.122;*/} //these numbers are calculated using the bin content of the ratio of gen. matched spectrum to gen. inclusive (unmatched). See macros/stat_err_scaling.cxx.
      if (i == 1) {scaling = scalefactors->GetBinContent(scalex->FindBin(25));/*old: 1.082;*/}
      //function from algorithm.h is max(a,b). Can't take 3 arguments, so have to do max(max(a,b),c).
      if (i == 2) {scaling = max(scalefactors->GetBinContent(scalex->FindBin(30)),scalefactors->GetBinContent(scalex->FindBin(35)));/*old: 1.062;*/}
      double binerror = reco_noms_copy[i]->GetBinError(j);
      reco_noms_copy[i]->SetBinError(j,(double) binerror*scaling);
      cout << "bin " << j << " error: " << binerror << endl;
      cout << "corresponding raw data error: " << dats[i]->GetBinError(j) << endl;
      cout << "scaling bin " << j << " error by: " << scaling << endl;
      cout << "scaled error: " << binerror*scaling << endl;
    }
  }
  
  cout << "F" << endl;
  
  string ftitle = "unfolded";
  if (groomflag) {
    ftitle = "groomed_"+ftitle;
  }

  TFile *fout = new TFile(("~/jetmass2/out/unfold/FINAL_"+ftitle+radius+"_paper.root").c_str(),"RECREATE");
  fout->cd();
    
  for (int i = 0; i < nBins; ++ i) {
    reco_IP2s[i]->Write(); //systematics
    reco_IP6s[i]->Write();
    reco_TSs[i]->Write();
    reco_TUs[i]->Write();
    reco_HC50s[i]->Write();
    reco_DSs[i]->Write();
    reco_GSs[i]->Write();
    //reco_MSs[i]->Write();
    reco_h7smear1D_2025->Write();
    reco_h7smear1D_2530->Write();
    reco_h7smear1D_3040->Write();
    reco_p8smear1D_2025->Write();
    reco_p8smear1D_2530->Write();
    reco_p8smear1D_3040->Write();
    
    env_HCs[i]->Write(); //systematic envelopes
    env_uns[i]->Write();
    nets[i]->Write();
        
      
    reco_noms_copy[i]->Write(); //unfolded datapoints
    dats[i]->Write(); //raw datapoints
    w_systs[i]->Write(); //unfolded data with errors equal to net systematic uncertainty
  }
  
  cout << "G" << endl;
  
  fout->Close();
  

  return 0;
}
