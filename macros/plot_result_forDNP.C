//Isaac Mooney, WSU, August 2019
//This file plots the unfolded (groomed) mass results for the DNP talk

#include "Plots_old.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

using namespace std;

//when running functions with arguments in an interactive root session, do e.g. root "plot_result(4, 0, 0)" - the quotation marks are necessary

/*
//this is a temp function to produce a plot for Raghav -- delete it!
void plot_result_forDNP(int radius, bool groom, int anim) {
  gROOT->ForceStyle(); //forces use of ~/rootlogon.C's style settings 
  
  string radstring = to_string(radius);
  string animstring = to_string(anim);
  
  radstring = "_R0"+radstring;//to obtain e.g. _R04 as in filenames
  string fstart = "";
  string hname = "m";
  string htitle = "M";
  if (groom) {
    fstart = "groomed_";
    hname = "mg";
    htitle = "M_{g}";
  }
  
  const int nRows = 1; //0.4
  const int nCols = 1;//2; //20-30,30-45 //!@!
  
  TCanvas *cws = new TCanvas("cws","cws",600,600);
  //  DivideCanvas(cws,"0",nCols,nRows);
  
  //open files from which to pull lots of curves                                                        
  //will use this file for a side-by-side comparison of groomed to ungroomed
  TFile *funfold_ungroom = new TFile(("~/jetmass2/out/unfold/unfolded"+radstring+".root").c_str(),"READ");
  TFile *funfold;
  if (radstring == "_R04") {funfold = new TFile(("~/jetmass2/out/unfold/"+fstart+"unfolded"+radstring+".root").c_str(),"READ");}//keeping the systematics for R = 0.4 consistent with the old method, since it already had preliminary status
  else {funfold = new TFile(("~/jetmass2/out/unfold/"+fstart+"unfolded"+radstring+"_test.root").c_str(),"READ");}
  TFile *fp6 = new TFile(("~/jetmass2/out/sim/hists/unmatched_hists"+radstring+"_bindropped.root").c_str(),"READ");
  TFile *fp8 = new TFile(("~/jetmass2/production/out/pythia/hists/pythia8"+radstring+"_undecayed_hists.root").c_str(),"READ");
  TFile *fh7 = new TFile(("~/jetmass/production/macros/hists/hists_allsim_lowzgremoved"+radstring+".root").c_str(),"READ");//temporarily using the old files (they're the same, but will point to new ones later)

  //get histograms from the files
  vector<TH1D*> recos = {(TH1D*) funfold->Get("nom_0"),(TH1D*) funfold->Get("nom_1")};
  vector<TH1D*> w_systs = {(TH1D*) funfold->Get("w_systs_0"),(TH1D*) funfold->Get("w_systs_1")};
  vector<TH1D*> recos_ungroom = {(TH1D*) funfold_ungroom->Get("nom_0"),(TH1D*) funfold_ungroom->Get("nom_1")};
  vector<TH1D*> w_systs_ungroom = {(TH1D*) funfold_ungroom->Get("w_systs_0"),(TH1D*) funfold_ungroom->Get("w_systs_1")};
  TH2D* p6_2D = (TH2D*) fp6->Get(("PL_"+hname+"_v_pt").c_str());
  TH2D* p8_2D = (TH2D*) fp8->Get((hname+"vpt").c_str());
  TH2D* h7_2D = (TH2D*) fh7->Get((hname+"vpt_h7off").c_str());//internal name will be changed later
  TH2D* PL_2D = (TH2D*) fp8->Get("PLmvpt");//fPL->Get("PLmvHLpt");
  TH2D* PL_g_2D = (TH2D*) fp8->Get("PLmgvpt");//fPL->Get("PLmgvHLpt");
  TH2D* p8_q_2D = (TH2D*) fp8->Get("mvpt_q");
  TH2D* p8_g_2D = (TH2D*) fp8->Get("mvpt_g");
  TH2D* p8_n_2D = (TH2D*) fp8->Get("mvpt_n");

  TH1D* dummy_left = new TH1D("dummy_left","",10,0,9.999);

  TAxis* p8_axis = p8_2D->GetYaxis();
  double ranges[nCols + 1] = {(double) p8_axis->FindBin(20), (double) p8_axis->FindBin(30)}; //!@!
  string pts[nCols + 1] = {"20","30"}; //!@!
  
  vector<TH1D*> p6s = Projection2D(p6_2D,nCols,ranges,"x");
  vector<TH1D*> p8s = Projection2D(p8_2D,nCols,ranges,"x");
  vector<TH1D*> h7s = Projection2D(h7_2D,nCols,ranges,"x");
  vector<TH1D*> PLs = Projection2D(PL_2D,nCols,ranges,"x");
  vector<TH1D*> PL_gs = Projection2D(PL_g_2D,nCols,ranges,"x");
  vector<TH1D*> p8_qs = Projection2D(p8_q_2D,nCols,ranges,"x");
  vector<TH1D*> p8_gs = Projection2D(p8_g_2D,nCols,ranges,"x");
  vector<TH1D*> p8_ns = Projection2D(p8_n_2D,nCols,ranges,"x");
  
  //get qvg content fractions
  //vector<double> qs = {(double) p8_qs[0]->Integral()/(double)(p8_qs[0]->Integral()+p8_gs[0]->Integral()+p8_ns[0]->Integral()),
  //(double) p8_qs[1]->Integral()/(double)(p8_qs[1]->Integral()+p8_gs[1]->Integral()+p8_ns[1]->Integral())}; //!@!
  //vector<double> gs = {(double) p8_gs[0]->Integral()/(double)(p8_qs[0]->Integral()+p8_gs[0]->Integral()+p8_ns[0]->Integral()),
  //(double) p8_gs[1]->Integral()/(double)(p8_qs[1]->Integral()+p8_gs[1]->Integral()+p8_ns[1]->Integral())}; //!@!
  //cout << "Q AND G FRACS: " << qs[0] << " " << qs[1] << " " << gs[0] << " " << gs[1] << endl; //!@!

  //clones for proper error bars for the mean mass later
  vector<TH1D*> p6_clones;
  vector<TH1D*> p8_clones;
  vector<TH1D*> h7_clones;
  
  vector<TH1D*> p6_ratios;
  vector<TH1D*> p8_ratios;
  vector<TH1D*> h7_ratios;
  vector<TH1D*> PL_ratios;
  vector<TH1D*> PL_g_ratios;
  vector<TH1D*> p8_q_ratios;
  vector<TH1D*> p8_g_ratios;

  vector<TH1D*> p6_ratios_bar;
  vector<TH1D*> p8_ratios_bar;
  vector<TH1D*> h7_ratios_bar;
  vector<TH1D*> PL_ratios_bar;
  vector<TH1D*> PL_g_ratios_bar;
  vector<TH1D*> p8_q_ratios_bar;
  vector<TH1D*> p8_g_ratios_bar;

  vector<TH1D*> ratio_systematics;

  //this will be the first "hist" drawn so needs to have correct axis properties
  dummy_left->GetYaxis()->SetRangeUser(0,1.999);
  
  dummy_left->GetYaxis()->SetTitle("MC / data");
  dummy_left->GetYaxis()->SetTitleSize(0.15);
  dummy_left->GetYaxis()->SetTitleOffset(0.6);
  dummy_left->GetYaxis()->SetNdivisions(505);
  dummy_left->GetYaxis()->SetLabelSize(0.17);

  dummy_left->GetXaxis()->SetTitle((htitle+" [GeV/c^{2}]").c_str());
  dummy_left->GetXaxis()->SetTitleSize(0.2);
  dummy_left->GetXaxis()->SetTitleOffset(0.8);
  dummy_left->GetXaxis()->SetLabelSize(0.17);
  dummy_left->GetXaxis()->SetNdivisions(505);
  
  for (unsigned i = 0; i < nCols; ++ i) {    
    cout << "UNFOLDED DATA MEAN " << i << ": " << w_systs[i]->GetMean() << endl;
    Prettify1D(recos[i],kRed,kFullStar,2,kRed,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,10,0,0.5);
    Prettify1D(w_systs[i],kRed,kFullStar,0,kRed,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,10,0.001,0.4999);
    Prettify1D(recos_ungroom[i],kBlue,kFullStar,2,kBlue,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,10,0,0.5);
    Prettify1D(w_systs_ungroom[i],kBlue,kFullStar,0,kBlue,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,10,0.001,0.4999);

    w_systs[i]->SetFillColor(kRed - 10); w_systs[i]->SetFillStyle(1001);
    
    Prettify1DwLineStyle(p6s[i], kBlue, kSolid,2, (htitle+" [GeV/c^{2}]").c_str(), ("1/N dN/d"+htitle).c_str(),0,10,0,0.5);
    Prettify1DwLineStyle(p8s[i],kBlack,kSolid,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,10,0,0.5);
    Prettify1DwLineStyle(p8_qs[i],kOrange-1,kSolid,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,10,0,0.5);
    Prettify1DwLineStyle(p8_gs[i],kGreen+2,kSolid,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,10,0,0.5);
    Prettify1DwLineStyle(h7s[i],kMagenta,kSolid,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,10,0,0.5);
    Prettify1DwLineStyle(PLs[i],kBlack,kDotted,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,10,0,0.5);
    Prettify1DwLineStyle(PL_gs[i],kBlack,kDashed,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,10,0,0.5);
    
    p6_ratios.push_back((TH1D*) p6s[i]->Clone(("p6ratio"+to_string(i)).c_str()));
    p8_ratios.push_back((TH1D*) p8s[i]->Clone(("p8ratio"+to_string(i)).c_str()));
    h7_ratios.push_back((TH1D*) h7s[i]->Clone(("h7ratio"+to_string(i)).c_str()));
    PL_ratios.push_back((TH1D*) PLs[i]->Clone(("PLratio"+to_string(i)).c_str()));
    PL_g_ratios.push_back((TH1D*) PL_gs[i]->Clone(("PL_gratio"+to_string(i)).c_str()));
    p8_q_ratios.push_back((TH1D*) p8_qs[i]->Clone(("p8_qratio"+to_string(i)).c_str()));
    p8_g_ratios.push_back((TH1D*) p8_gs[i]->Clone(("p8_gratio"+to_string(i)).c_str()));
    
    p6_clones.push_back((TH1D*) p6s[i]->Clone(("p6clone"+to_string(i)).c_str()));
    p8_clones.push_back((TH1D*) p8s[i]->Clone(("p8clone"+to_string(i)).c_str()));
    h7_clones.push_back((TH1D*) h7s[i]->Clone(("h7clone"+to_string(i)).c_str()));
    
    ratio_systematics.push_back((TH1D*) w_systs[i]->Clone(("ratio_systs"+to_string(i)).c_str()));
    
    p6s[i]->Sumw2(0);
    p8s[i]->Sumw2(0);
    h7s[i]->Sumw2(0);
    PLs[i]->Sumw2(0);
    PL_gs[i]->Sumw2(0);
    p8_qs[i]->Sumw2(0);
    p8_gs[i]->Sumw2(0);
    recos[i]->Sumw2();

    p6_ratios[i]->Divide(recos[i]);
    p8_ratios[i]->Divide(recos[i]);
    h7_ratios[i]->Divide(recos[i]);
    PL_ratios[i]->Divide(recos[i]);
    PL_g_ratios[i]->Divide(recos[i]);    
    p8_q_ratios[i]->Divide(recos[i]);
    p8_g_ratios[i]->Divide(recos[i]);
    
    ratio_systematics[i]->Divide(recos[i]);
    
    p6_ratios[i]->GetYaxis()->SetRangeUser(0,2);
    p8_ratios[i]->GetYaxis()->SetRangeUser(0,2);
    h7_ratios[i]->GetYaxis()->SetRangeUser(0,2);
    PL_ratios[i]->GetYaxis()->SetRangeUser(0,2);
    PL_g_ratios[i]->GetYaxis()->SetRangeUser(0,2);
    p8_q_ratios[i]->GetYaxis()->SetRangeUser(0,2);
    p8_g_ratios[i]->GetYaxis()->SetRangeUser(0,2);
    ratio_systematics[i]->GetYaxis()->SetRangeUser(0,2);

    w_systs[i]->GetYaxis()->SetLabelSize(0.065);
    
    p6_ratios_bar.push_back((TH1D*) p6_ratios[i]->Clone(("p6ratio_bar"+to_string(i)).c_str()));
    p8_ratios_bar.push_back((TH1D*) p8_ratios[i]->Clone(("p8ratio_bar"+to_string(i)).c_str()));
    h7_ratios_bar.push_back((TH1D*) h7_ratios[i]->Clone(("h7ratio_bar"+to_string(i)).c_str()));
    PL_ratios_bar.push_back((TH1D*) PL_ratios[i]->Clone(("PLratio_bar"+to_string(i)).c_str()));
    PL_g_ratios_bar.push_back((TH1D*) PL_g_ratios[i]->Clone(("PL_gratio_bar"+to_string(i)).c_str()));
    p8_q_ratios_bar.push_back((TH1D*) p8_q_ratios[i]->Clone(("p8_qratio_bar"+to_string(i)).c_str()));
    p8_g_ratios_bar.push_back((TH1D*) p8_g_ratios[i]->Clone(("p8_gratio_bar"+to_string(i)).c_str()));
    
    p6_ratios_bar[i]->Sumw2(0);
    p8_ratios_bar[i]->Sumw2(0);
    h7_ratios_bar[i]->Sumw2(0);
    PL_ratios_bar[i]->Sumw2(0);
    PL_g_ratios_bar[i]->Sumw2(0);
    p8_q_ratios_bar[i]->Sumw2(0);
    p8_g_ratios_bar[i]->Sumw2(0);

    //    p8_qs[i]->Scale((double)qs[i]); p8_gs[i]->Scale((double)gs[i]); //!@!
    p8_qs[i]->Sumw2(0); p8_gs[i]->Sumw2(0);
    
  }

  TLegend *tleg = new TLegend(0.55,0.3,0.75,0.4); tleg->SetBorderSize(0);
  TLegend *tleg2 = new TLegend(0.45,0.45,0.7,0.6); tleg2->SetBorderSize(0);
  TLegend *tleg3 = new TLegend(0.45,0.4,0.7,0.45); tleg3->SetBorderSize(0);
  TH1D* for_legend = (TH1D*) w_systs[0]->Clone("for_legend"); for_legend->SetMarkerSize(2);
  tleg->AddEntry(for_legend,"STAR","pf");
  
  tleg->SetTextSize(0.04); tleg2->SetTextSize(0.04); tleg3->SetTextSize(0.04);
  
  if (anim == 3 || anim == 4) {tleg2->AddEntry(p6s[0],"PYTHIA-6 Perugia","l");} //!@!
  if (anim == 4) {
    tleg2->AddEntry(h7s[0],"HERWIG-7 EE4C","l");
    tleg2->AddEntry(p8s[0],"PYTHIA-8 Monash","l");
    tleg3->AddEntry(PLs[0],"PYTHIA-8 parton jets","l");
  }
  
  TLatex *ttitle = new TLatex(); TLatex *slice = new TLatex();
  ttitle->SetTextSize(0.05); slice->SetTextSize(0.05);
  TLatex *tprelim = new TLatex();
  tprelim->SetTextSize(0.05); tprelim->SetTextColor(kRed);
  
  TLine *one = new TLine(0,1,10,1); one->SetLineStyle(kDashed);
  
  for (int i = 0; i < nCols; ++ i) {
    //cws->cd(i+1); //!@!
    TPad *padup = new TPad(("padup"+to_string(i)).c_str(),("padup"+to_string(i)).c_str(),0,0.3,1,1.0);
    padup->SetBottomMargin(0);
    padup->SetTopMargin(0.1);
    padup->SetRightMargin(0.2);
    padup->SetLeftMargin(0.2);
    padup->Draw();
    padup->cd();
    w_systs[i]->Draw("E3same");
    if (anim == 3 || anim == 4) {p6s[i]->Draw("Csame");} //!@!                                                 
    if (anim == 4) {
      p8s[i]->Draw("Csame");
      h7s[i]->Draw("Csame");
      PLs[i]->Draw("Csame");
    }                   
    
    recos[i]->Draw("same");
    tleg->Draw("same");
    tleg2->Draw("same");
    tleg3->Draw("same");
    
    slice->DrawLatex(0.5,0.33,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
    
    ttitle->DrawLatex(0.5,0.41, "pp #sqrt{s} = 200 GeV");
    ttitle->DrawLatex(0.5,0.37, "anti-k_{T}, R = 0.4, |#eta_{jet}| < 1-R");
    tprelim->DrawLatex(0.5,0.45, "STAR Preliminary");
    
    
    cws->cd(); //!@!
    TPad *paddown = new TPad(("paddown"+to_string(i)).c_str(),("paddown"+to_string(i)).c_str(),0,0.05,1,0.3);
    paddown->SetTopMargin(0);
    paddown->SetBottomMargin(0.4);
    paddown->SetRightMargin(0.2);
    paddown->SetLeftMargin(0.2);
    
    paddown->Draw();
    paddown->cd();
    if (i==0) {dummy_left->Draw();}
    ratio_systematics[i]->Draw("E3same");
    
    if (anim == 3 || anim == 4) { //!@!
      p6_ratios[i]->Draw("E1same");
      p6_ratios_bar[i]->Draw("same");
    }
    if (anim == 4) {
      p8_ratios[i]->Draw("E1same");
      p8_ratios_bar[i]->Draw("same");
      h7_ratios[i]->Draw("E1same");
      h7_ratios_bar[i]->Draw("same");
      PL_ratios[i]->Draw("E1same");
      PL_ratios_bar[i]->Draw("same");
    }
    
    one->Draw("same");
    
  }
  
  cws->SaveAs("~/jetmass2/plots/single_pt_groomed_mass_for_raghav.pdf");
  
  return;
}
*/

//this function plots a 1x2 plot of mass as a function of pT, but with ratios of MC to data
TH1D* plot_result_forDNP(int radius, bool groom, int anim) {
  gROOT->ForceStyle(); //forces use of ~/rootlogon.C's style settings 
  
  string radstring = to_string(radius);
  string animstring = to_string(anim);
  
  radstring = "_R0"+radstring;//to obtain e.g. _R04 as in filenames
  string fstart = "";
  string hname = "m";
  string htitle = "M";
  if (groom) {
    fstart = "groomed_";
    hname = "mg";
    htitle = "M_{g}";
  }
  
  const int nRows = 1; //0.4
  const int nCols = 2; //20-30,30-45
  
  TCanvas *cws = new TCanvas("cws","cws",800,500);
  DivideCanvas(cws,"0",nCols,nRows);
  
  //open files from which to pull lots of curves                                                        
  //will use this file for a side-by-side comparison of groomed to ungroomed
  TFile *funfold_ungroom = new TFile(("~/jetmass2/out/unfold/unfolded"+radstring+".root").c_str(),"READ");
  TFile *funfold;
  if (radstring == "_R04") {funfold = new TFile(("~/jetmass2/out/unfold/"+fstart+"unfolded"+radstring+".root").c_str(),"READ");}//keeping the systematics for R = 0.4 consistent with the old method, since it already had preliminary status
  else {funfold = new TFile(("~/jetmass2/out/unfold/"+fstart+"unfolded"+radstring+"_test.root").c_str(),"READ");}
  TFile *fp6 = new TFile(("~/jetmass2/out/sim/hists/unmatched_hists"+radstring+"_bindropped.root").c_str(),"READ");
  TFile *fp8 = new TFile(("~/jetmass2/production/out/pythia/hists/pythia8"+radstring+"_undecayed_hists.root").c_str(),"READ");
  TFile *fh7 = new TFile(("~/jetmass/production/macros/hists/hists_allsim_lowzgremoved"+radstring+".root").c_str(),"READ");//temporarily using the old files (they're the same, but will point to new ones later)
  //TFile *fPL = new TFile(("~/jetmass2/production/out/pythia/hists/pythia8"+radstring+"_FSR_only_undecayed_hists.root").c_str(),"READ");

  TFile *fang = new TFile("~/jetmass2/out/sim/hists/Angantyr_for_HP.root","READ");
  
  //get histograms from the files
  vector<TH1D*> recos = {(TH1D*) funfold->Get("nom_0"),(TH1D*) funfold->Get("nom_1")};
  vector<TH1D*> w_systs = {(TH1D*) funfold->Get("w_systs_0"),(TH1D*) funfold->Get("w_systs_1")};
  vector<TH1D*> recos_ungroom = {(TH1D*) funfold_ungroom->Get("nom_0"),(TH1D*) funfold_ungroom->Get("nom_1")};
  vector<TH1D*> w_systs_ungroom = {(TH1D*) funfold_ungroom->Get("w_systs_0"),(TH1D*) funfold_ungroom->Get("w_systs_1")};
  TH2D* p6_2D = (TH2D*) fp6->Get(("PL_"+hname+"_v_pt").c_str());
  TH2D* p8_2D = (TH2D*) fp8->Get((hname+"vpt").c_str());
  TH2D* h7_2D = (TH2D*) fh7->Get((hname+"vpt_h7off").c_str());//internal name will be changed later
  TH2D* PL_2D = (TH2D*) fp8->Get("PLmvpt");//fPL->Get("PLmvHLpt");
  TH2D* PL_g_2D = (TH2D*) fp8->Get("PLmgvpt");//fPL->Get("PLmgvHLpt");
  TH2D* p8_q_2D = (TH2D*) fp8->Get("mvpt_q");
  TH2D* p8_g_2D = (TH2D*) fp8->Get("mvpt_g");
  TH2D* p8_n_2D = (TH2D*) fp8->Get("mvpt_n");

  TH1D* ang0_1D = (TH1D*) fang->Get("Angantyr0");
  TH1D* ang1_1D = (TH1D*) fang->Get("Angantyr1");
  vector<TH1D*> angs = {ang0_1D,ang1_1D};
  
  TH1D* dummy_left = new TH1D("dummy_left","",10,0,9.999);
  TH1D* dummy_right = new TH1D("dummy_right","",10,0.001,9.999);

  TAxis* p8_axis = p8_2D->GetYaxis();
  double ranges[nCols + 1] = {(double) p8_axis->FindBin(20), (double) p8_axis->FindBin(30), (double) p8_axis->FindBin(45)};
  string pts[nCols + 1] = {"20","30","45"};
  
  vector<TH1D*> p6s = Projection2D(p6_2D,nCols,ranges,"x");
  vector<TH1D*> p8s = Projection2D(p8_2D,nCols,ranges,"x");
  vector<TH1D*> h7s = Projection2D(h7_2D,nCols,ranges,"x");
  vector<TH1D*> PLs = Projection2D(PL_2D,nCols,ranges,"x");
  vector<TH1D*> PL_gs = Projection2D(PL_g_2D,nCols,ranges,"x");
  vector<TH1D*> p8_qs = Projection2D(p8_q_2D,nCols,ranges,"x");
  vector<TH1D*> p8_gs = Projection2D(p8_g_2D,nCols,ranges,"x");
  vector<TH1D*> p8_ns = Projection2D(p8_n_2D,nCols,ranges,"x");
  
  //get qvg content fractions
  vector<double> qs = {(double) p8_qs[0]->Integral()/(double)(p8_qs[0]->Integral()+p8_gs[0]->Integral()+p8_ns[0]->Integral()),
		       (double) p8_qs[1]->Integral()/(double)(p8_qs[1]->Integral()+p8_gs[1]->Integral()+p8_ns[1]->Integral())};
  vector<double> gs = {(double) p8_gs[0]->Integral()/(double)(p8_qs[0]->Integral()+p8_gs[0]->Integral()+p8_ns[0]->Integral()),
		       (double) p8_gs[1]->Integral()/(double)(p8_qs[1]->Integral()+p8_gs[1]->Integral()+p8_ns[1]->Integral())};
  cout << "Q AND G FRACS: " << qs[0] << " " << qs[1] << " " << gs[0] << " " << gs[1] << endl;
  
  //clones for proper error bars for the mean mass later
  vector<TH1D*> p6_clones;
  vector<TH1D*> p8_clones;
  vector<TH1D*> h7_clones;
  
  vector<TH1D*> p6_ratios;
  vector<TH1D*> p8_ratios;
  vector<TH1D*> h7_ratios;
  vector<TH1D*> PL_ratios;
  vector<TH1D*> PL_g_ratios;
  vector<TH1D*> p8_q_ratios;
  vector<TH1D*> p8_g_ratios;

  vector<TH1D*> p6_ratios_bar;
  vector<TH1D*> p8_ratios_bar;
  vector<TH1D*> h7_ratios_bar;
  vector<TH1D*> PL_ratios_bar;
  vector<TH1D*> PL_g_ratios_bar;
  vector<TH1D*> p8_q_ratios_bar;
  vector<TH1D*> p8_g_ratios_bar;

  vector<TH1D*> ratio_systematics;

  //this will be the first "hist" drawn so needs to have correct axis properties
  dummy_left->GetYaxis()->SetRangeUser(0,1.999);
  dummy_right->GetYaxis()->SetRangeUser(0,1.999);
  
  dummy_left->GetYaxis()->SetTitle("MC / data");
  dummy_left->GetYaxis()->SetTitleSize(0.15);
  dummy_left->GetYaxis()->SetTitleOffset(0.6);
  dummy_left->GetYaxis()->SetNdivisions(505);
  dummy_left->GetYaxis()->SetLabelSize(0.13);
  dummy_right->GetYaxis()->SetTitleSize(0.15);
  dummy_right->GetYaxis()->SetTitleOffset(0.6);
  dummy_right->GetYaxis()->SetNdivisions(505);
  dummy_right->GetYaxis()->SetLabelSize(0.13);

  dummy_left->GetXaxis()->SetTitle((htitle+" [GeV/c^{2}]").c_str());
  //  dummy_left->GetXaxis()->SetTitle("M_{(g)} [GeV/c^{2}]");
  dummy_left->GetXaxis()->SetTitleSize(0.2);
  dummy_left->GetXaxis()->SetTitleOffset(0.8);
  dummy_left->GetXaxis()->SetLabelSize(0.13);
  dummy_left->GetXaxis()->SetNdivisions(505);
  dummy_right->GetXaxis()->SetTitle((htitle+" [GeV/c^{2}]").c_str());
  //dummy_right->GetXaxis()->SetTitle("M_{(g)} [GeV/c^{2}]");
  dummy_right->GetXaxis()->SetTitleSize(0.2);
  dummy_right->GetXaxis()->SetTitleOffset(0.8);
  dummy_right->GetXaxis()->SetLabelSize(0.13);
  dummy_right->GetXaxis()->SetNdivisions(505);
  
  
  for (unsigned i = 0; i < nCols; ++ i) {    
    cout << "UNFOLDED DATA MEAN " << i << ": " << w_systs[i]->GetMean() << endl;
    Prettify1D(recos[i],kRed,kFullStar,2,kRed,(htitle+" [GeV/c^{2}]").c_str(),("1/N_{jet} dN_{jet}/d"+htitle+" [c^{2}/GeV]").c_str(),0,10,0,0.5);
    Prettify1D(w_systs[i],kRed,kFullStar,0,kRed,(htitle+" [GeV/c^{2}]").c_str(),("1/N_{jet} dN_{jet}/d"+htitle+" [c^{2}/GeV]").c_str(),0,10,0.001,0.4999);
    Prettify1D(recos_ungroom[i],kBlue,kFullStar,2,kBlue,(htitle+" [GeV/c^{2}]").c_str(),("1/N_{jet} dN_{jet}/d"+htitle+" [c^{2}/GeV]").c_str(),0,10,0,0.5);
    Prettify1D(w_systs_ungroom[i],kBlue,kFullStar,0,kBlue,(htitle+" [GeV/c^{2}]").c_str(),("1/N_{jet} dN_{jet}/d"+htitle+" [c^{2}/GeV]").c_str(),0,10,0.001,0.4999);

    w_systs[i]->SetFillColor(kRed - 10); w_systs[i]->SetFillStyle(1001);
    w_systs_ungroom[i]->SetFillColor(kBlue - 10); w_systs_ungroom[i]->SetFillStyle(1001);
    
    Prettify1DwLineStyle(p6s[i], kBlue, kSolid,2, (htitle+" [GeV/c^{2}]").c_str(), ("1/N_{jet} dN_{jet}/d"+htitle+" [c^{2}/GeV]").c_str(),0,10,0,0.5);
    Prettify1DwLineStyle(p8s[i],kBlack,kSolid,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N_{jet} dN_{jet}/d"+htitle+" [c^{2}/GeV]").c_str(),0,10,0,0.5);
    Prettify1DwLineStyle(angs[i],kViolet,kSolid,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N_{jet} dN_{jet}/d"+htitle+" [c^{2}/GeV]").c_str(),0,10,0,0.5);
    Prettify1DwLineStyle(p8_qs[i],kOrange-1,kSolid,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N_{jet} dN_{jet}/d"+htitle+" [c^{2}/GeV]").c_str(),0,10,0,0.5);
    Prettify1DwLineStyle(p8_gs[i],kGreen+2,kSolid,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N_{jet} dN_{jet}/d"+htitle+" [c^{2}/GeV]").c_str(),0,10,0,0.5);
    Prettify1DwLineStyle(h7s[i],kMagenta,kSolid,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N_{jet} dN_{jet}/d"+htitle+" [c^{2}/GeV]").c_str(),0,10,0,0.5);
    Prettify1DwLineStyle(PLs[i],kBlack,kDotted,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N_{jet} dN_{jet}/d"+htitle+" [c^{2}/GeV]").c_str(),0,10,0,0.5);
    Prettify1DwLineStyle(PL_gs[i],kBlack,kDashed,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N_{jet} dN_{jet}/d"+htitle+" [c^{2}/GeV]").c_str(),0,10,0,0.5);
    
    p6_ratios.push_back((TH1D*) p6s[i]->Clone(("p6ratio"+to_string(i)).c_str()));
    p8_ratios.push_back((TH1D*) p8s[i]->Clone(("p8ratio"+to_string(i)).c_str()));
    h7_ratios.push_back((TH1D*) h7s[i]->Clone(("h7ratio"+to_string(i)).c_str()));
    PL_ratios.push_back((TH1D*) PLs[i]->Clone(("PLratio"+to_string(i)).c_str()));
    PL_g_ratios.push_back((TH1D*) PL_gs[i]->Clone(("PL_gratio"+to_string(i)).c_str()));
    p8_q_ratios.push_back((TH1D*) p8_qs[i]->Clone(("p8_qratio"+to_string(i)).c_str()));
    p8_g_ratios.push_back((TH1D*) p8_gs[i]->Clone(("p8_gratio"+to_string(i)).c_str()));
    
    p6_clones.push_back((TH1D*) p6s[i]->Clone(("p6clone"+to_string(i)).c_str()));
    p8_clones.push_back((TH1D*) p8s[i]->Clone(("p8clone"+to_string(i)).c_str()));
    h7_clones.push_back((TH1D*) h7s[i]->Clone(("h7clone"+to_string(i)).c_str()));
    
    ratio_systematics.push_back((TH1D*) w_systs[i]->Clone(("ratio_systs"+to_string(i)).c_str()));
    
    p6s[i]->Sumw2(0);
    p8s[i]->Sumw2(0);
    h7s[i]->Sumw2(0);
    PLs[i]->Sumw2(0);
    PL_gs[i]->Sumw2(0);
    p8_qs[i]->Sumw2(0);
    p8_gs[i]->Sumw2(0);
    angs[i]->Sumw2(0);
    recos[i]->Sumw2();

    p6_ratios[i]->Divide(recos[i]);
    p8_ratios[i]->Divide(recos[i]);
    h7_ratios[i]->Divide(recos[i]);
    PL_ratios[i]->Divide(recos[i]);
    PL_g_ratios[i]->Divide(recos[i]);    
    p8_q_ratios[i]->Divide(recos[i]);
    p8_g_ratios[i]->Divide(recos[i]);
    
    ratio_systematics[i]->Divide(recos[i]);
    
    p6_ratios[i]->GetYaxis()->SetRangeUser(0,2);
    p8_ratios[i]->GetYaxis()->SetRangeUser(0,2);
    h7_ratios[i]->GetYaxis()->SetRangeUser(0,2);
    PL_ratios[i]->GetYaxis()->SetRangeUser(0,2);
    PL_g_ratios[i]->GetYaxis()->SetRangeUser(0,2);
    p8_q_ratios[i]->GetYaxis()->SetRangeUser(0,2);
    p8_g_ratios[i]->GetYaxis()->SetRangeUser(0,2);
    ratio_systematics[i]->GetYaxis()->SetRangeUser(0,2);

    w_systs[i]->GetYaxis()->SetLabelSize(0.075);
    
    p6_ratios_bar.push_back((TH1D*) p6_ratios[i]->Clone(("p6ratio_bar"+to_string(i)).c_str()));
    p8_ratios_bar.push_back((TH1D*) p8_ratios[i]->Clone(("p8ratio_bar"+to_string(i)).c_str()));
    h7_ratios_bar.push_back((TH1D*) h7_ratios[i]->Clone(("h7ratio_bar"+to_string(i)).c_str()));
    PL_ratios_bar.push_back((TH1D*) PL_ratios[i]->Clone(("PLratio_bar"+to_string(i)).c_str()));
    PL_g_ratios_bar.push_back((TH1D*) PL_g_ratios[i]->Clone(("PL_gratio_bar"+to_string(i)).c_str()));
    p8_q_ratios_bar.push_back((TH1D*) p8_q_ratios[i]->Clone(("p8_qratio_bar"+to_string(i)).c_str()));
    p8_g_ratios_bar.push_back((TH1D*) p8_g_ratios[i]->Clone(("p8_gratio_bar"+to_string(i)).c_str()));
    
    p6_ratios_bar[i]->Sumw2(0);
    p8_ratios_bar[i]->Sumw2(0);
    h7_ratios_bar[i]->Sumw2(0);
    PL_ratios_bar[i]->Sumw2(0);
    PL_g_ratios_bar[i]->Sumw2(0);
    p8_q_ratios_bar[i]->Sumw2(0);
    p8_g_ratios_bar[i]->Sumw2(0);

    p8_qs[i]->Scale((double)qs[i]); p8_gs[i]->Scale((double)gs[i]);
    p8_qs[i]->Sumw2(0); p8_gs[i]->Sumw2(0);
    
  }

  TLegend *tleg = new TLegend(0.6,0.45,0.9,0.6); tleg->SetBorderSize(0);
  //TLegend *tleg2 = new TLegend(0.1,0.5,0.5,0.65); tleg2->SetBorderSize(0);
  TLegend *tleg2; if (!groom) {tleg2 = new TLegend(0.2,0.52,0.55,0.75);} 
  if (groom) {tleg2 = new TLegend(0.2,0.46,0.55,0.67);}
  tleg2->SetBorderSize(0);
  
  //TLegend *tleg2 = new TLegend(0.248965, 0.454372, 0.598997, 0.652933); tleg2->SetBorderSize(0);
  TLegend *tleg3 = new TLegend(0.25,0.53,0.65,0.65); tleg3->SetBorderSize(0);
  TH1D* for_legend = (TH1D*) w_systs[0]->Clone("for_legend"); for_legend->SetMarkerSize(2);
  TH1D* for_legend_ungroom = (TH1D*) w_systs_ungroom[0]->Clone("for_legend_ungroom"); for_legend_ungroom->SetMarkerSize(2);
  tleg->AddEntry(for_legend,"STAR","pf");
  //  tleg->AddEntry(for_legend_ungroom,"STAR M","pf");
  if (anim == 1) {tleg3->AddEntry(PLs[0],"PYTHIA-8 parton jets","l");}
  
  // if ((!groom) && (anim == 2)) {
  // tleg2->AddEntry(p8_qs[0],"PYTHIA-8 q jets","l");
  // tleg2->AddEntry(p8_gs[0],"PYTHIA-8 g jets","l");
  //}
  
  if (anim == 3) {
    tleg2->AddEntry(p6s[0],"PYTHIA-6 Perugia 2012","l");
  }
  if (anim == 3 || anim == 4) {
    tleg2->AddEntry(h7s[0],"HERWIG-7 EE4C","l");
    tleg2->AddEntry(p8s[0],"PYTHIA-8 Monash","l");
    //tleg2->AddEntry(angs[0],"PYTHIA-8 Monash Angantyr","l");
  }
  
  if (groom) {
    //    tleg3->AddEntry(PL_gs[0],"PYTHIA-8 SD parton jets","l");
  }
  
  TLatex *ttitle = new TLatex(); TLatex *slice = new TLatex();
  ttitle->SetTextSize(0.07); slice->SetTextSize(0.07);
  TLatex *tprelim = new TLatex();
  tprelim->SetTextSize(0.07); tprelim->SetTextColor(kRed);
  
  TLine *one = new TLine(0,1,14,1); one->SetLineStyle(kDashed);
  
  for (int i = 0; i < nCols; ++ i) {
    cws->cd(i+1);
    TPad *padup = new TPad(("padup"+to_string(i)).c_str(),("padup"+to_string(i)).c_str(),0,0.4,1,1.0);
    padup->SetBottomMargin(0);
    padup->SetTopMargin(0);
    if (i!=0) {padup->SetLeftMargin(0);}
    padup->SetRightMargin(0);
    padup->Draw();
    padup->cd();
    w_systs[i]->Draw("E3same");
    //w_systs_ungroom[i]->Draw("E3same");
    if (anim == 1) {PLs[i]->Draw("Csame");}
    //if ((!groom) && (anim == 2)) {
    //p8_qs[i]->Draw("Csame");
    // p8_gs[i]->Draw("Csame");
    //}
    if (anim == 3) {
      p6s[i]->Draw("Csame");
    }                                                 
    if (anim == 3 || anim == 4) {
      p8s[i]->Draw("Csame");
      //angs[i]->Draw("Csame");
      h7s[i]->Draw("Csame");
    }                   
    if (groom) {
      //      PL_gs[i]->Draw("Csame");
    }
    recos[i]->Draw("same");
    //recos_ungroom[i]->Draw("same");
    if (i == 0) {tleg->Draw("same");}
    if (i == 1) {tleg2->Draw("same");}
    if (i == 1 && anim == 1) {tleg3->Draw("same");}

    if (i == 0) {
      slice->DrawLatexNDC(0.4,0.7,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
    }
    if (i == 1 && !groom) {
      slice->DrawLatexNDC(0.4,0.8,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
    }
    if (i == 1 && groom) {
      slice->DrawLatex(2.2,0.36,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
    }
    if (i == 0) {
      ttitle->DrawLatex(2.2,0.45, "p+p #sqrt{s} = 200 GeV");
      ttitle->DrawLatex(2.2,0.41, "anti-k_{T}, R = 0.4, |#eta_{jet}| < 1-R");
    }
    if (i == 1) {
      tprelim->DrawLatex(2.2,0.45, "STAR Preliminary");
    }
    if (i == 1 && groom) {
      ttitle->DrawLatex(2.2,0.41,"SoftDrop z_{cut} = 0.1, #beta = 0");
    }
    
    
    cws->cd(i+1);
    TPad *paddown = new TPad(("paddown"+to_string(i)).c_str(),("paddown"+to_string(i)).c_str(),0,0.05,1,0.4);
    paddown->SetTopMargin(0);
    paddown->SetBottomMargin(0.4);
    paddown->SetRightMargin(0);
    paddown->SetLeftMargin(0);
    if (i==0) {paddown->SetLeftMargin(0.2);}
    paddown->Draw();
    paddown->cd();
    if (i==0) {dummy_left->Draw();}
    else if (i!=0) {dummy_right->Draw();}
    ratio_systematics[i]->Draw("E3same");
    if (anim == 1) {
      PL_ratios[i]->Draw("E1same");
      PL_ratios_bar[i]->Draw("same");
    }
    // if ((!groom) && (anim == 2)) {
    // p8_q_ratios[i]->Draw("E1same");
    // p8_q_ratios_bar[i]->Draw("same");
    // p8_g_ratios[i]->Draw("E1same");
    // p8_g_ratios_bar[i]->Draw("same");
    // }
    if (anim == 3) {
      p6_ratios[i]->Draw("E1same");
      p6_ratios_bar[i]->Draw("same");
    }
    if (anim == 3 || anim == 4) {
      p8_ratios[i]->Draw("E1same");
      p8_ratios_bar[i]->Draw("same");
      h7_ratios[i]->Draw("E1same");
      h7_ratios_bar[i]->Draw("same");
    }
     
    //if(!groom) {//}
    //if(groom) {PL_g_ratios[i]->Draw("E1same");PL_g_ratios_bar[i]->Draw("same");}
    
    one->Draw("same");
    
  }
  //cws->SaveAs(("~/jetmass2/plots/DNP_talk/"+fstart+"jet_mass_2panel_w_ratios_"+animstring+"_prelim.pdf").c_str());
  cws->SaveAs(("~/jetmass2/plots/HP_proceedings/"+fstart+"pp_jetmass.pdf").c_str());
  
  //MEAN v. pT:

  TCanvas *cm = new TCanvas("cm","cm",800,600);
  cm->cd();
  
  const int nBins = nCols;
  double bin_edges[nBins+1] = {20,30,45};
  
  TH1D *reco_m = new TH1D("reco_m","",nBins,bin_edges);
  //BLOCK COMMENT STARTS HERE
  /*TH1D *reco_m_copy = new TH1D("reco_m_copy","",nBins,bin_edges);
  TH1D *p6_m = new TH1D("p6_m","",nBins,bin_edges);
  TH1D *p8_m = new TH1D("p8_m","",nBins,bin_edges);
  TH1D *h7_m = new TH1D("h7_m","",nBins,bin_edges);

  for (int i = 0; i < nCols; ++ i) {
    reco_m->SetBinContent(i+1, w_systs[i]->GetMean());
    reco_m->SetBinError(i+1, w_systs[i]->GetMeanError());
    
    p6_m->SetBinContent(i+1, p6_clones[i]->GetMean());
    p6_m->SetBinError(i+1, p6_clones[i]->GetMeanError());
    
    p8_m->SetBinContent(i+1, p8_clones[i]->GetMean());
    p8_m->SetBinError(i+1, p8_clones[i]->GetMeanError());
    
    h7_m->SetBinContent(i+1, h7_clones[i]->GetMean());
    h7_m->SetBinError(i+1, h7_clones[i]->GetMeanError());

  }
  
  reco_m_copy = (TH1D*) reco_m->Clone("reco_m_copy");
  for (int i = 0; i < nCols; ++ i) {
    reco_m_copy->SetBinError(i,0);
  }

  reco_m_copy->Sumw2();
  p6_m->Sumw2();
  p8_m->Sumw2();
  h7_m->Sumw2();
  
  Prettify1D(reco_m_copy,kRed,kFullStar,2.5,kRed,"p_{T} [GeV/c]",("<"+htitle+">"+" [GeV/c^{2}]").c_str(),20,40,0,10);
  Prettify1D(reco_m,kRed,kFullStar,0,kRed,"p_{T} [GeV/c]",("<"+htitle+">"+" [GeV/c^{2}]").c_str(),20,40,0,10);
  reco_m->SetFillColor(kRed - 10); reco_m->SetFillStyle(1001);
  Prettify1D(p6_m, kBlue, kOpenSquare,2,kBlue,"p_{T} [GeV/c]",("<"+htitle+">"+" [GeV/c^{2}]").c_str(),20,40,0,10);
  Prettify1D(p8_m,kBlack,kFullCircle,2,kBlack,"p_{T} [GeV/c]",("<"+htitle+">"+" [GeV/c^{2}]").c_str(),20,40,0,10);
  Prettify1D(h7_m,kMagenta,kOpenCircle,2,kMagenta,"p_{T} [GeV/c]",("<"+htitle+">"+" [GeV/c^{2}]").c_str(),20,40,0,10);
  
  //  gStyle->SetErrorX();

  reco_m->GetXaxis()->SetTitleSize(0.05);
  reco_m->GetYaxis()->SetTitleSize(0.05);
  
  TLegend *tm = new TLegend(0.25,0.25,0.45,0.35); tm->SetBorderSize(0);
  TLegend *tm2 = new TLegend(0.55,0.25,0.75,0.35); tm2->SetBorderSize(0);
  TH1D* for_legend_m = (TH1D*) reco_m->Clone("for_legend_m"); for_legend_m->SetMarkerSize(2);
  tm->AddEntry(for_legend_m,"Unfolded data","pfe");
  tm->AddEntry(p6_m,"PYTHIA-6","pe");
  tm2->AddEntry(h7_m,"HERWIG-7","pe");
  tm2->AddEntry(p8_m,"PYTHIA-8","pe");
  tm->SetTextSize(0.05); tm2->SetTextSize(0.05);
  ttitle->SetTextSize(0.048);
  
  reco_m->Draw("E2");
  reco_m_copy->Draw("same");
  p6_m->Draw("same");
  p8_m->Draw("same");
  h7_m->Draw("same");
  tm->Draw("same");
  tm2->Draw("same");
  
  ttitle->DrawLatex(22,9, "pp 200 GeV run12 JP2");
  ttitle->DrawLatex(22,8, "anti-k_{T} jets, |#eta| < 1 - R");
  if (groom) {
    ttitle->DrawLatex(22,7,"SoftDrop z_{cut} = 0.1, #beta = 0");
  }

  //  cm->SaveAs(("~/jetmass2/plots/DNP_talk/"+fstart+"mean_jet_mass.pdf").c_str());
*/  //BLOCK COMMENT ENDS HERE
  return reco_m;
}

//this function plots the mean mass and groomed mass as a function of jet R for two different pT ranges
void plot_result_forDNP() {
  gROOT->ForceStyle();

  const int nCols = 3; //(number of radius bins)
 
  //open files from which to pull lots of curves                                                                                              
  TFile *funfold_02 = new TFile("~/jetmass2/out/unfold/unfolded_R02.root","READ");
  TFile *funfold_04 = new TFile("~/jetmass2/out/unfold/unfolded_R04.root","READ");
  TFile *funfold_06 = new TFile("~/jetmass2/out/unfold/unfolded_R06.root","READ");
  TFile *fgunfold_02 = new TFile("~/jetmass2/out/unfold/groomed_unfolded_R02.root","READ");
  TFile *fgunfold_04 = new TFile("~/jetmass2/out/unfold/groomed_unfolded_R04.root","READ");
  TFile *fgunfold_06 = new TFile("~/jetmass2/out/unfold/groomed_unfolded_R06.root","READ");

  //get histograms from the files
  vector<TH1D*> w_systs_02 = {(TH1D*) funfold_02->Get("w_systs_0"),(TH1D*) funfold_02->Get("w_systs_1")};
  vector<TH1D*> w_systs_04 = {(TH1D*) funfold_04->Get("w_systs_0"),(TH1D*) funfold_04->Get("w_systs_1")};
  vector<TH1D*> w_systs_06 = {(TH1D*) funfold_06->Get("w_systs_0"),(TH1D*) funfold_06->Get("w_systs_1")};
  vector<TH1D*> gw_systs_02 = {(TH1D*) fgunfold_02->Get("w_systs_0"),(TH1D*) fgunfold_02->Get("w_systs_1")};
  vector<TH1D*> gw_systs_04 = {(TH1D*) fgunfold_04->Get("w_systs_0"),(TH1D*) fgunfold_04->Get("w_systs_1")};
  vector<TH1D*> gw_systs_06 = {(TH1D*) fgunfold_06->Get("w_systs_0"),(TH1D*) fgunfold_06->Get("w_systs_1")};
  
  vector<TH1D*> w_systs_low = {(TH1D*) funfold_02->Get("w_systs_0"),(TH1D*) funfold_04->Get("w_systs_0"),(TH1D*) funfold_06->Get("w_systs_0")};
  vector<TH1D*> gw_systs_low = {(TH1D*) fgunfold_02->Get("w_systs_0"),(TH1D*) fgunfold_04->Get("w_systs_0"),(TH1D*) fgunfold_06->Get("w_systs_0")};
  vector<TH1D*> w_systs_high = {(TH1D*) funfold_02->Get("w_systs_1"),(TH1D*) funfold_04->Get("w_systs_1"),(TH1D*) funfold_06->Get("w_systs_1")};
  vector<TH1D*> gw_systs_high = {(TH1D*) fgunfold_02->Get("w_systs_1"),(TH1D*) fgunfold_04->Get("w_systs_1"),(TH1D*) fgunfold_06->Get("w_systs_1")};

  
  TH1D *reco_m_low = new TH1D("reco_m_low","",3,0.1,0.7);
  TH1D *reco_m_low_copy = new TH1D("reco_m_low_copy","",3,0.1,0.7);
  TH1D *reco_m_high = new TH1D("reco_m_high","",3,0.1,0.7);
  TH1D *reco_m_high_copy = new TH1D("reco_m_high_copy","",3,0.1,0.7);
  TH1D *greco_m_low = new TH1D("greco_m_low","",3,0.1,0.7);
  TH1D *greco_m_low_copy = new TH1D("greco_m_low_copy","",3,0.1,0.7);
  TH1D *greco_m_high = new TH1D("greco_m_high","",3,0.1,0.7);
  TH1D *greco_m_high_copy = new TH1D("greco_m_high_copy","",3,0.1,0.7);

  for (int i = 0; i < nCols; ++ i) {
    reco_m_low->SetBinContent(i+1, w_systs_low[i]->GetMean());
    reco_m_low->SetBinError(i+1, w_systs_low[i]->GetMeanError());
    reco_m_high->SetBinContent(i+1, w_systs_high[i]->GetMean());
    reco_m_high->SetBinError(i+1, w_systs_high[i]->GetMeanError());
    greco_m_low->SetBinContent(i+1, gw_systs_low[i]->GetMean());
    greco_m_low->SetBinError(i+1, gw_systs_low[i]->GetMeanError());
    greco_m_high->SetBinContent(i+1, gw_systs_high[i]->GetMean());
    greco_m_high->SetBinError(i+1, gw_systs_high[i]->GetMeanError());
  }
 
  reco_m_low_copy = (TH1D*) reco_m_low->Clone("reco_m_low_copy");
  reco_m_high_copy = (TH1D*) reco_m_high->Clone("reco_m_high_copy");
  greco_m_low_copy = (TH1D*) greco_m_low->Clone("greco_m_low_copy");
  greco_m_high_copy = (TH1D*) greco_m_high->Clone("greco_m_high_copy");
  for (int i = 0; i < nCols; ++ i) {
    reco_m_low_copy->SetBinError(i,0);
    reco_m_high_copy->SetBinError(i,0);
    greco_m_low_copy->SetBinError(i,0);
    greco_m_high_copy->SetBinError(i,0);
  }
  
  reco_m_low_copy->Sumw2();
  reco_m_high_copy->Sumw2();
  greco_m_low_copy->Sumw2();
  greco_m_high_copy->Sumw2();
  
  Prettify1D(reco_m_low_copy,kRed,kFullStar,2.5,kRed,"jet radius, R","<M> [GeV/c^{2}]",0.1,0.7,2,7);
  Prettify1D(reco_m_high_copy,kBlue,kFullStar,2.5,kBlue,"jet radius, R","<M> [GeV/c^{2}]",0.1,0.7,2,7);
  Prettify1D(greco_m_low_copy,kRed,kFullStar,2.5,kRed,"jet radius, R","<M_{g}> [GeV/c^{2}]",0.1,0.7,2,7);
  Prettify1D(greco_m_high_copy,kBlue,kFullStar,2.5,kBlue,"jet radius, R","<M_{g}> [GeV/c^{2}]",0.1,0.7,2,7);

  Prettify1D(reco_m_low,kRed,kFullStar,0,kRed,"jet radius, R","<M> [GeV/c^{2}]",0.1,0.7,2,7);
  reco_m_low->SetFillColor(kRed - 10); reco_m_low->SetFillStyle(1001);
  Prettify1D(greco_m_low,kRed,kFullStar,0,kRed,"jet radius, R","<M_{g}> [GeV/c^{2}]",0.1,0.7,2,7);
  greco_m_low->SetFillColor(kRed-10); greco_m_low->SetFillStyle(1001);
  Prettify1D(reco_m_high,kBlue,kFullStar,0,kBlue,"jet radius, R","<M> [GeV/c^{2}]",0.1,0.7,2,7);
  reco_m_high->SetFillColor(kBlue - 10); reco_m_high->SetFillStyle(1001);
  Prettify1D(greco_m_high,kBlue,kFullStar,0,kBlue,"jet radius, R","<M_{g}> [GeV/c^{2}]",0.1,0.7,2,7);
  greco_m_high->SetFillColor(kBlue-10); greco_m_high->SetFillStyle(1001);

  //define dummies and prettify them                                                                                                                            
  TH1D* hdummy = new TH1D("hdummy","",10,0.1,0.6999);
  TH1D* hdummyg = new TH1D("hdummyg","",10,0.1001,0.6999);
  hdummy->GetYaxis()->SetRangeUser(2,7);
  hdummyg->GetYaxis()->SetRangeUser(2,7);
  hdummy->GetXaxis()->SetTitle("jet radius, R");
  hdummyg->GetXaxis()->SetTitle("jet radius, R");
  hdummy->GetYaxis()->SetTitle("<M> [GeV/c^{2}]");
  hdummyg->GetYaxis()->SetTitle("<M_{g}> [GeV/c^{2}]");
  
  hdummyg->GetYaxis()->SetTitleOffset(1.1);
  hdummyg->GetXaxis()->SetTitleSize(0.07);
  hdummy->GetXaxis()->SetTitleSize(0.07);
  
  gStyle->SetPadTickY(0);
  gStyle->SetErrorX();
  
  TGaxis *raxis = new TGaxis(0.7,2,0.7,7,2,7,505,"+L");  
  raxis->SetLabelFont(42);
  raxis->SetLabelSize(0.05);
  raxis->SetTitleFont(42);
  raxis->SetTitleSize(0.05);
  raxis->SetTitle("<M_{g}> [GeV/c^{2}]");

  hdummy->GetXaxis()->SetNdivisions(504);
  hdummy->GetYaxis()->SetNdivisions(505);
  hdummyg->GetXaxis()->SetNdivisions(504);
  hdummyg->GetYaxis()->SetNdivisions(505);

  hdummy->GetXaxis()->SetTitleSize(0.05);
  hdummyg->GetXaxis()->SetTitleSize(0.05);
  hdummy->GetYaxis()->SetTitleSize(0.05);
  hdummy->GetXaxis()->SetTitleOffset(1.2);
  hdummyg->GetXaxis()->SetTitleOffset(1.2);
  
  //legends
  TLatex *title = new TLatex();
  TLegend *tm = new TLegend(0.2,0.6,0.5,0.8); tm->SetBorderSize(0);
  TH1D* for_legend_m_low = (TH1D*) reco_m_low->Clone("for_legend_m_low"); for_legend_m_low->SetMarkerSize(2);
  TH1D* for_legend_mg_low = (TH1D*) greco_m_low->Clone("for_legend_mg_low"); for_legend_mg_low->SetMarkerSize(2);
  TH1D* for_legend_m_high = (TH1D*) reco_m_high->Clone("for_legend_m_high"); for_legend_m_high->SetMarkerSize(2);
  TH1D* for_legend_mg_high = (TH1D*) greco_m_high->Clone("for_legend_mg_high"); for_legend_mg_high->SetMarkerSize(2);
  tm->AddEntry(for_legend_m_low,"p_{T} #in (20,30) GeV/c","pfe");
  tm->AddEntry(for_legend_m_high,"p_{T} #in (30,45) GeV/c","pfe");

  tm->SetTextSize(0.04);
  title->SetTextSize(0.048);
  
  TCanvas *cm = new TCanvas("cm","cm",800,500);
  cm->cd();
  TPad *pleft = new TPad("pleft","",0,0,0.5,1);
  TPad *pright = new TPad("pright","",0.5,0,1,1);
  pleft->SetRightMargin(0);
  pleft->SetLeftMargin(0.2);
  pleft->SetBottomMargin(0.17);
  pright->SetLeftMargin(0);
  pright->SetRightMargin(0.2);
  pright->SetBottomMargin(0.17);
  pleft->Draw();
  pright->Draw();
  
  pleft->cd();  
  hdummy->Draw();
  reco_m_low->Draw("E2same");
  reco_m_low_copy->Draw("same");
  reco_m_high->Draw("E2same");
  reco_m_high_copy->Draw("same");
  title->DrawLatex(0.12,6.6, "pp 200 GeV run12 JP2");
  title->DrawLatex(0.12,6.1, "anti-k_{T} jets, |#eta| < 1 - R");
  title->DrawLatex(0.12,5.6,"SoftDrop z_{cut} = 0.1, #beta = 0");
  
  pright->cd();
  
  hdummyg->Draw();
  greco_m_low->Draw("E2sameY+");
  greco_m_low_copy->Draw("same");
  greco_m_high->Draw("E2sameY+");
  greco_m_high_copy->Draw("same");
  raxis->Draw("same");
  tm->Draw("same");

  //cm->SaveAs("~/jetmass2/plots/DNP_talk/mass_mean_radial_scan_2panel_w_groomed.pdf");

}

/*
//this function calls the previous function in order to plot the mass and groomed mass individually, then pull the mean mass as a function of pT from each, and plot them together. Note: when using this call, the latter helper function call will overwrite the plots of the former - so only use this to get the mean histograms plotted together, and nothing more.
void plot_result_forDNP() {
  const int nCols = 2; //(number of pT bins)
  
  TH1D* reco_m = (TH1D*) plot_result_forDNP(4,0,0);
  TH1D* reco_mg = (TH1D*) plot_result_forDNP(4,1,0);
  
  TCanvas *ccomp = new TCanvas("ccomp","ccomp",800,500);
  ccomp->cd();

  TH1D* reco_m_copy = (TH1D*) reco_m->Clone("reco_m_copy");
  TH1D* reco_mg_copy = (TH1D*) reco_mg->Clone("reco_mg_copy");
  for (int i = 0; i < nCols; ++ i) {
    reco_m_copy->SetBinError(i,0);
    reco_mg_copy->SetBinError(i,0);
  }

  reco_m_copy->Sumw2();
  reco_mg_copy->Sumw2();
  Prettify1D(reco_m_copy,kRed,kFullStar,2.5,kRed,"p_{T} [GeV/c]","<M> [GeV/c^{2}]",20,45,3,6);
  Prettify1D(reco_mg_copy,kBlue,kFullStar,2.5,kBlue,"p_{T} [GeV/c]","<M_{g}> [GeV/c^{2}]",20,45,3,6);
  reco_m->GetYaxis()->SetRangeUser(3,6);
  reco_mg->GetYaxis()->SetRangeUser(3,6);
  reco_mg->SetMarkerColor(kBlue);
  reco_mg->SetLineColor(kBlue);
  reco_mg->SetFillColor(kBlue-10);
  
  gStyle->SetPadTickY(0);
  
  TGaxis *raxis = new TGaxis(40,3,40,6,3,6,510,"+L");  
  raxis->SetLabelFont(42);
  raxis->SetLabelSize(0.05);
  raxis->SetTitleFont(42);
  raxis->SetTitleSize(0.05);
  raxis->SetTitle("<M_{g}> [GeV/c^{2}]");

  reco_m->GetXaxis()->SetNdivisions(504);
  
  //legends
  TLatex *title = new TLatex();
  TLegend *tm = new TLegend(0.5,0.2,0.8,0.35); tm->SetBorderSize(0);
  TH1D* for_legend_m = (TH1D*) reco_m->Clone("for_legend_m"); for_legend_m->SetMarkerSize(2);
  TH1D* for_legend_mg = (TH1D*) reco_mg->Clone("for_legend_mg"); for_legend_mg->SetMarkerSize(2);
  tm->AddEntry(for_legend_m,"Mass","pfe");
  tm->AddEntry(for_legend_mg,"Groomed mass","pfe");
  tm->SetTextSize(0.05);
  title->SetTextSize(0.048);
  
  reco_m->Draw("E2");
  reco_m_copy->Draw("same");
  reco_mg->Draw("E2sameY+");
  reco_mg_copy->Draw("same");
  tm->Draw("same");
  raxis->Draw("same");
  
  title->DrawLatex(21,5.75, "pp 200 GeV run12 JP2");
  title->DrawLatex(21,5.55, "anti-k_{T} jets, |#eta| < 1 - R");
  title->DrawLatex(21,5.35,"SoftDrop z_{cut} = 0.1, #beta = 0");

  //  ccomp->SaveAs("~/jetmass2/plots/DNP_talk/mean_jetmass_groom_v_ungroom.pdf");
  
  return;
}
*/
//this function plots a 3x2 plot of mass as a function of radius and pT
void plot_result_forDNP(bool groom, int anim) {
  gROOT->ForceStyle();

  string fstart = "";
  string hname = "m";
  string htitle = "M";
  if (groom) {
    fstart = "groomed_";
    hname = "mg";
    htitle = "M_{g}";
  }
  
  string radstring;
  string animstring = to_string(anim);
  
  const int nRows = 3; //0.2, 0.4, 0.6
  const int nCols = 2/*3*/;//now: 20-30,30-45 //other times: 20-25,25-30,30-40
  
  TCanvas *cws = new TCanvas("cws","cws",800,700);
  DivideCanvas(cws,"0",nCols,nRows);
  
  for (int iRow = 0; iRow < 3; ++ iRow) {
    if (iRow == 0) {radstring = "_R02";}
    if (iRow == 1) {radstring = "_R04";}
    if (iRow == 2) {radstring = "_R06";}
    //open files from which to pull lots of curves                                   
    TFile *funfold;
    if (iRow == 1) {funfold = new TFile(("~/jetmass2/out/unfold/"+fstart+"unfolded"+radstring+".root").c_str(),"READ"); }//this keeps the R = 0.4 net systematic consistent with the mass preliminary for now
    else {funfold = new TFile(("~/jetmass2/out/unfold/"+fstart+"unfolded"+radstring+"_test.root").c_str(),"READ");}//change back if you change pt ranges!
    TFile *fp6 = new TFile(("~/jetmass2/out/sim/hists/unmatched_hists"+radstring+"_bindropped.root").c_str(),"READ");
    TFile *fp8 = new TFile(("~/jetmass2/production/out/pythia/hists/pythia8"+radstring+"_undecayed_hists.root").c_str(),"READ");
    TFile *fh7 = new TFile(("~/jetmass/production/macros/hists/hists_allsim_lowzgremoved"+radstring+".root").c_str(),"READ");//temporarily using the old files (they're the same, but will point to new ones later)
    //TFile *fPL = new TFile(("~/jetmass2/production/out/pythia/hists/pythia8"+radstring+"_FSR_only_undecayed_hists.root").c_str(),"READ");

    //get histograms from the files
    vector<TH1D*> recos = {(TH1D*) funfold->Get("nom_0"),(TH1D*) funfold->Get("nom_1")/*,(TH1D*) funfold->Get("nom_2")*/};
    vector<TH1D*> w_systs = {(TH1D*) funfold->Get("w_systs_0"),(TH1D*) funfold->Get("w_systs_1")/*,(TH1D*) funfold->Get("w_systs_2")*/};
    TH2D* p6_2D = (TH2D*) fp6->Get(("PL_"+hname+"_v_pt").c_str());
    TH2D* p8_2D = (TH2D*) fp8->Get((hname+"vpt").c_str());
    TH2D* h7_2D = (TH2D*) fh7->Get((hname+"vpt_h7off").c_str());//internal name will be changed later
    TH2D* PL_2D = (TH2D*) fp8->Get("PLmvpt");//fPL->Get("PLmvHLpt");
    TH2D* PL_g_2D = (TH2D*) fp8->Get("PLmgvpt");//fPL->Get("PLmgvHLpt");
    TH2D* p8_q_2D = (TH2D*) fp8->Get("mvpt_q");
    TH2D* p8_g_2D = (TH2D*) fp8->Get("mvpt_g");
    TH2D* p8_n_2D = (TH2D*) fp8->Get("mvpt_n");
    
    TAxis* p8_axis = p8_2D->GetYaxis();
    double ranges[nCols + 1] = {(double) p8_axis->FindBin(20)/*, (double) p8_axis->FindBin(25)*/, (double) p8_axis->FindBin(30), (double) p8_axis->FindBin(45)/*(40)*/};
    string pts[nCols + 1] = {"20"/*,"25"*/,"30","45"/*"40"*/};
  
    vector<TH1D*> p6s = {(TH1D*)p6_2D->ProjectionX(("p6_0"+to_string(iRow)).c_str(),ranges[0],ranges[1]-1),(TH1D*)p6_2D->ProjectionX(("p6_1"+to_string(iRow)).c_str(),ranges[1],ranges[2]-1)/*,(TH1D*)p6_2D->ProjectionX(("p6_2"+to_string(iRow)).c_str(),ranges[2],ranges[3]-1)*/};
    vector<TH1D*> p8s = {(TH1D*)p8_2D->ProjectionX(("p8_0"+to_string(iRow)).c_str(),ranges[0],ranges[1]-1),(TH1D*)p8_2D->ProjectionX(("p8_1"+to_string(iRow)).c_str(),ranges[1],ranges[2]-1)/*,(TH1D*)p8_2D->ProjectionX(("p8_2"+to_string(iRow)).c_str(),ranges[2],ranges[3]-1)*/};
    vector<TH1D*> h7s = {(TH1D*)h7_2D->ProjectionX(("h7_0"+to_string(iRow)).c_str(),ranges[0],ranges[1]-1),(TH1D*)h7_2D->ProjectionX(("h7_1"+to_string(iRow)).c_str(),ranges[1],ranges[2]-1)/*,(TH1D*)h7_2D->ProjectionX(("h7_2"+to_string(iRow)).c_str(),ranges[2],ranges[3]-1)*/};
    vector<TH1D*> PLs = {(TH1D*)PL_2D->ProjectionX(("PL_0"+to_string(iRow)).c_str(),ranges[0],ranges[1]-1),(TH1D*)PL_2D->ProjectionX(("PL_1"+to_string(iRow)).c_str(),ranges[1],ranges[2]-1)/*,(TH1D*)PL_2D->ProjectionX(("PL_2"+to_string(iRow)).c_str(),ranges[2],ranges[3]-1)*/};
    vector<TH1D*> PL_gs = {(TH1D*)PL_g_2D->ProjectionX(("PL_g_0"+to_string(iRow)).c_str(),ranges[0],ranges[1]-1),(TH1D*)PL_g_2D->ProjectionX(("PL_g_1"+to_string(iRow)).c_str(),ranges[1],ranges[2]-1)/*,(TH1D*)PL_g_2D->ProjectionX(("PL_g_2"+to_string(iRow)).c_str(),ranges[2],ranges[3]-1)*/};
    vector<TH1D*> p8_qs = {(TH1D*)p8_q_2D->ProjectionX(("p8_q_0"+to_string(iRow)).c_str(),ranges[0],ranges[1]-1),(TH1D*)p8_q_2D->ProjectionX(("p8_q_1"+to_string(iRow)).c_str(),ranges[1],ranges[2]-1)/*,(TH1D*)p8_q_2D->ProjectionX(("p8_q_2"+to_string(iRow)).c_str(),ranges[2],ranges[3]-1)*/};
    vector<TH1D*> p8_gs = {(TH1D*)p8_g_2D->ProjectionX(("p8_g_0"+to_string(iRow)).c_str(),ranges[0],ranges[1]-1),(TH1D*)p8_g_2D->ProjectionX(("p8_g_1"+to_string(iRow)).c_str(),ranges[1],ranges[2]-1)/*,(TH1D*)p8_g_2D->ProjectionX(("p8_g_2"+to_string(iRow)).c_str(),ranges[2],ranges[3]-1)*/};
    vector<TH1D*> p8_ns = {(TH1D*)p8_n_2D->ProjectionX(("p8_n_0"+to_string(iRow)).c_str(),ranges[0],ranges[1]-1),(TH1D*)p8_n_2D->ProjectionX(("p8_n_1"+to_string(iRow)).c_str(),ranges[1],ranges[2]-1)/*,(TH1D*)p8_n_2D->ProjectionX(("p8_n_2"+to_string(iRow)).c_str(),ranges[2],ranges[3]-1)*/};  
    //vector<TH1D*> p6s = Projection2D(p6_2D,nCols,ranges,"x");
    //vector<TH1D*> p8s = Projection2D(p8_2D,nCols,ranges,"x");
    //vector<TH1D*> h7s = Projection2D(h7_2D,nCols,ranges,"x");
    //vector<TH1D*> PLs = Projection2D(PL_2D,nCols,ranges,"x");
    //vector<TH1D*> PL_gs = Projection2D(PL_g_2D,nCols,ranges,"x");
    //vector<TH1D*> p8_qs = Projection2D(p8_q_2D,nCols,ranges,"x");
    //vector<TH1D*> p8_gs = Projection2D(p8_g_2D,nCols,ranges,"x");
  
    //get qvg content fractions
    vector<double> qs = {(double) p8_qs[0]->Integral()/(double)(p8_qs[0]->Integral()+p8_gs[0]->Integral()+p8_ns[0]->Integral()),
			 (double) p8_qs[1]->Integral()/(double)(p8_qs[1]->Integral()+p8_gs[1]->Integral()+p8_ns[1]->Integral())};
    vector<double> gs = {(double) p8_gs[0]->Integral()/(double)(p8_qs[0]->Integral()+p8_gs[0]->Integral()+p8_ns[0]->Integral()),
			 (double) p8_gs[1]->Integral()/(double)(p8_qs[1]->Integral()+p8_gs[1]->Integral()+p8_ns[1]->Integral())};
    cout << "Q AND G FRACS: q0, q1, g0, g1: " << qs[0] << " " << qs[1] << " " << gs[0] << " " << gs[1] << ", " << radstring <<  endl;
    
    
    TH1D* ldummy = new TH1D("ldummy","",14,0,13.999);
    TH1D* rdummy = new TH1D("rdummy","",14,0.0001,13.999);
    //ldummy->GetXaxis()->SetTitle((htitle+" [GeV/c^{2}]").c_str());
    //rdummy->GetXaxis()->SetTitle((htitle+" [GeV/c^{2}]").c_str());
    ldummy->GetYaxis()->SetTitle(("1/N dN/d"+htitle).c_str());
    ldummy->GetYaxis()->SetRangeUser(0.0001,0.4999);
    rdummy->GetYaxis()->SetRangeUser(0.0001,0.4999);
    ldummy->GetYaxis()->SetLabelSize(0.075);
    ldummy->GetXaxis()->SetLabelSize(0.075);
    rdummy->GetXaxis()->SetLabelSize(0.075);
    ldummy->GetYaxis()->SetNdivisions(505);
    rdummy->GetYaxis()->SetNdivisions(505);
    /* ldummy->GetXaxis()->SetTitleOffset(1.2);
    rdummy->GetXaxis()->SetTitleOffset(1.2);
    ldummy->GetXaxis()->SetTitleSize(0.1);
    rdummy->GetXaxis()->SetTitleSize(0.1);
    */
    for (unsigned i = 0; i < nCols; ++ i) {
      cout << "UNFOLDED DATA MEAN " << i << ": " << w_systs[i]->GetMean() << endl;
      Prettify1D(recos[i],kRed,kFullStar,2,kRed,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.75);
      Prettify1D(w_systs[i],kRed,kFullStar,0,kRed,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.75);
      w_systs[i]->SetFillColor(kRed - 10); w_systs[i]->SetFillStyle(1001);
      Prettify1DwLineStyle(p6s[i], kBlue, kSolid,2, (htitle+" [GeV/c^{2}]").c_str(), ("1/N dN/d"+htitle).c_str(),0,14,0,0.75);
      Prettify1DwLineStyle(p8s[i],kBlack,kSolid,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.75);
      Prettify1DwLineStyle(p8_qs[i],kOrange-1,kSolid,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.75);
      Prettify1DwLineStyle(p8_gs[i],kGreen+2,kSolid,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.75);
      Prettify1DwLineStyle(h7s[i],kMagenta,kSolid,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.75);
      Prettify1DwLineStyle(PLs[i],kBlack,kDotted,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.75);
      Prettify1DwLineStyle(PL_gs[i],kBlack,kDashed,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.75);
      
      p6s[i]->Sumw2(0);
      p8s[i]->Sumw2(0);
      p8_qs[i]->Sumw2(0);
      p8_gs[i]->Sumw2(0);
      h7s[i]->Sumw2(0);
      PLs[i]->Sumw2(0);
      PL_gs[i]->Sumw2(0);
      
      p8_qs[i]->Scale((double)qs[i]);
      p8_gs[i]->Scale((double)gs[i]);
      p8_qs[i]->Sumw2(0); p8_gs[i]->Sumw2(0);
    }
    
    TLegend *tleg = new TLegend(0.4,0.45,0.8,0.6); tleg->SetBorderSize(0);
    TLegend *tleg2 = new TLegend(0.4,0.45,0.8,0.6); tleg2->SetBorderSize(0);
    TLegend *tleg3 = new TLegend(0.1,0.5,0.55,0.65); tleg3->SetBorderSize(0);
    TH1D* for_legend = (TH1D*) w_systs[0]->Clone("for_legend"); for_legend->SetMarkerSize(2);
    tleg->AddEntry(for_legend,"STAR","pf");
    if (anim == 1) {
      tleg3->AddEntry(PLs[0],"PYTHIA-8 parton jets","l");
    }
    if ((!groom) && (anim == 2)) {
      tleg2->AddEntry(p8_qs[0],"PYTHIA-8 q jets","l");
      tleg2->AddEntry(p8_gs[0],"PYTHIA-8 g jets","l");
    }
    if (anim == 3/* || anim == 5*/) {
      tleg->AddEntry(p6s[0],"PYTHIA-6 Perugia","l");
    }
    if (anim == 4 /*|| anim == 5*/) {
      tleg2->AddEntry(h7s[0],"HERWIG-7 EE4C","l");
    } 
    if (anim == 5) {
      tleg2->AddEntry(p8s[0],"PYTHIA-8 Monash","l");
    }
    
    if (groom) {
      //tleg3->AddEntry(PL_gs[0],"PYTHIA-8 SD parton jets","l");
    }
    

    TLatex *ttitle = new TLatex(); TLatex *slice = new TLatex();
    ttitle->SetTextSize(0.08); slice->SetTextSize(0.07);
    TLatex *tprelim = new TLatex();
    tprelim->SetTextSize(0.08); tprelim->SetTextColor(kRed);

    for (int i = 0; i < nCols; ++ i) {
      cws->cd((i+1)+nCols*iRow);
      if (i == 0) {ldummy->Draw();} if (i == 1) {rdummy->Draw();}
      w_systs[i]->Draw("E3same");
      if (anim == 1) {PLs[i]->Draw("Csame");}
      if ((!groom) && (anim == 2)) {
	p8_qs[i]->Draw("Csame");
	p8_gs[i]->Draw("Csame");
      }
      if (anim == 3 /*|| anim == 5*/) {p6s[i]->Draw("Csame");}
      if (anim == 4 /*|| anim == 5*/) {h7s[i]->Draw("Csame");}                                                        
      if (anim == 5) {p8s[i]->Draw("Csame");}
      
      if (groom) {
	//PL_gs[i]->Draw("Csame");
      }
      
      recos[i]->Draw("same");

      int rad_int = (iRow+1)*2;//iRow = 0,1,2 representating radius = 0.2,0.4,0.6 -> 2,4,6
      slice->DrawLatexNDC(0.4,0.77,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
      slice->DrawLatexNDC(0.4,0.67,("R = 0."+to_string(rad_int)).c_str());

      if (i == 0 && iRow == 0) {
	ttitle->DrawLatex(4,0.45, "pp #sqrt{s} = 200 GeV");
      }
      if (i == 0 && iRow == 1) {
	ttitle->DrawLatex(4,0.45, "anti-k_{T}, |#eta_{jet}| < 1-R");
      }
      if (i == 1 && iRow == 1 && groom) {
	ttitle->DrawLatex(4,0.45,"SoftDrop z_{cut} = 0.1, #beta = 0");
      }
      if (i == 1 && iRow == 0) {
	tprelim->DrawLatex(4,0.45,"STAR Preliminary");
      }
      if (i == 1 && iRow == 0) {
	tleg->Draw("same");
      }
      if (i == 1 && iRow == 1) {
	tleg2->Draw("same");
      }
      if (i == 1 && iRow == 2 && anim == 1) {
	tleg3->Draw("same");
      }
    }
  }

  //cws->SaveAs(("~/jetmass2/plots/DNP_talk/"+fstart+"jet_mass_6panel_"+animstring+"_prelim.pdf").c_str());

  return;
}
/*
//this file plots a 1x3 of the (groomed) mass result with additional MC curves (e.g. quark & gluon jet discrimination)
//radius format: 4 = 0.4, ...
void plot_result_forDNP (int radius, bool groom) {
  //  gStyle->SetPalette(kPastel);
  gROOT->ForceStyle();
  
  string radstring = to_string(radius);
  
  radstring = "_R0"+radstring;//to obtain e.g. _R04 as in filenames
  string fstart = "";
  string hname = "m";
  string htitle = "M";
  if (groom) {
    fstart = "groomed_";
    hname = "mg";
    htitle = "M_{g}";
  }

  //open files from which to pull lots of curves
  TFile *funfold = new TFile(("~/jetmass2/out/unfold/"+fstart+"unfolded"+radstring+".root").c_str(),"READ");//change this if you change nBins!
  TFile *fp6 = new TFile(("~/jetmass2/out/sim/hists/unmatched_hists"+radstring+"_bindropped.root").c_str(),"READ");
  TFile *fp8 = new TFile(("~/jetmass2/production/out/pythia/hists/pythia8"+radstring+"_undecayed_hists.root").c_str(),"READ");
  TFile *fh7 = new TFile(("~/jetmass/production/macros/hists/hists_allsim_lowzgremoved"+radstring+".root").c_str(),"READ");//temporarily using the old files (they're the same, but will point to new ones later)
  //TFile *fPL = new TFile(("~/jetmass2/production/out/pythia/hists/pythia8"+radstring+"_FSR_only_undecayed_hists.root").c_str(),"READ");
  
  //get histograms from the files
  vector<TH1D*> recos = {(TH1D*) funfold->Get("nom_0"),(TH1D*) funfold->Get("nom_1")};
  vector<TH1D*> w_systs = {(TH1D*) funfold->Get("w_systs_0"),(TH1D*) funfold->Get("w_systs_1")};
  TH2D* p6_2D = (TH2D*) fp6->Get(("PL_"+hname+"_v_pt").c_str());
  TH2D* p8_2D = (TH2D*) fp8->Get((hname+"vpt").c_str());
  TH2D* h7_2D = (TH2D*) fh7->Get((hname+"vpt_h7off").c_str());//internal name will be changed later
  TH2D* PL_2D = (TH2D*) fp8->Get("PLmvpt");//fPL->Get("PLmvHLpt");
  TH2D* PL_g_2D = (TH2D*) fp8->Get("PLmgvpt");//fPL->Get("PLmgvHLpt");
  TH2D* p8_q_2D = (TH2D*) fp8->Get("mvpt_q");
  TH2D* p8_g_2D = (TH2D*) fp8->Get("mvpt_g");

  const int nBins = 2;
  TAxis* p8_axis = p8_2D->GetYaxis();
  double ranges[nBins + 1] = {(double) p8_axis->FindBin(20), (double) p8_axis->FindBin(30), (double) p8_axis->FindBin(45)};
  string pts[nBins + 1] = {"20","30","45"};
  
  vector<TH1D*> p6s = Projection2D(p6_2D,nBins,ranges,"x");
  vector<TH1D*> p8s = Projection2D(p8_2D,nBins,ranges,"x");
  vector<TH1D*> h7s = Projection2D(h7_2D,nBins,ranges,"x");
  vector<TH1D*> PLs = Projection2D(PL_2D,nBins,ranges,"x");
  vector<TH1D*> PL_gs = Projection2D(PL_g_2D,nBins,ranges,"x");
  vector<TH1D*> p8_qs = Projection2D(p8_q_2D,nBins,ranges,"x");
  vector<TH1D*> p8_gs = Projection2D(p8_g_2D,nBins,ranges,"x");
  
  TH1D* ldummy = new TH1D("ldummy","",14,0,13.999);
  TH1D* rdummy = new TH1D("rdummy","",14,0.0001,13.999);
  ldummy->GetXaxis()->SetTitle((htitle+" [GeV/c^{2}]").c_str());
  rdummy->GetXaxis()->SetTitle((htitle+" [GeV/c^{2}]").c_str());
  ldummy->GetYaxis()->SetTitle(("1/N dN/d"+htitle).c_str());
  ldummy->GetYaxis()->SetRangeUser(0,0.4999);
  rdummy->GetYaxis()->SetRangeUser(0,0.4999);
  
  
  for (unsigned i = 0; i < nBins; ++ i) {
    Prettify1D(recos[i],kRed,kFullStar,4,kRed,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.5);
    Prettify1D(w_systs[i],kRed,kFullStar,0,kRed,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.5);
    w_systs[i]->SetFillColor(kRed - 10); w_systs[i]->SetFillStyle(1001);
    Prettify1DwLineStyle(p6s[i], kBlue, kSolid,3, (htitle+" [GeV/c^{2}]").c_str(), ("1/N dN/d"+htitle).c_str(),0,14,0,0.5);
    Prettify1DwLineStyle(p8s[i],kBlack,kSolid,3,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.5);
    Prettify1DwLineStyle(p8_qs[i],kOrange-1,kSolid,3,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.5);
    Prettify1DwLineStyle(p8_gs[i],kGreen+2,kSolid,3,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.5);                        
    Prettify1DwLineStyle(h7s[i],kMagenta,kSolid,3,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.5);
    Prettify1DwLineStyle(PLs[i],kBlack,kDotted,3,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.5);
    Prettify1DwLineStyle(PL_gs[i],kBlack,kDashed,3,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.5); 
    
    p6s[i]->Sumw2(0);
    p8s[i]->Sumw2(0);
    p8_qs[i]->Sumw2(0);
    p8_gs[i]->Sumw2(0);
    h7s[i]->Sumw2(0);
    PLs[i]->Sumw2(0);
    PL_gs[i]->Sumw2(0);
  }

  TCanvas *cws = new TCanvas("cws","cws",800,500);
  DivideCanvas(cws,"0",nBins,1);

  TLegend *tleg = new TLegend(0.55,0.5,0.75,0.7); tleg->SetBorderSize(0);
  TLegend *tleg2 = new TLegend(0.1,0.58,0.3,0.75); tleg2->SetBorderSize(0);
  TH1D* for_legend = (TH1D*) w_systs[0]->Clone("for_legend"); for_legend->SetMarkerSize(2);
  tleg->AddEntry(for_legend,"STAR","pf");
  tleg->AddEntry(p6s[0],"PYTHIA-6","l");
  tleg2->AddEntry(h7s[0],"HERWIG-7","l");
  tleg->AddEntry(p8s[0],"PYTHIA-8","l");
  
  tleg2->AddEntry(PLs[0],"PYTHIA-8 parton jets","l");
  if (groom) {
    tleg2->AddEntry(PL_gs[0],"PYTHIA-8 SD parton jets","l");
  }
  
  if (!groom) {
    tleg2->AddEntry(p8_qs[0],"PYTHIA-8 q jets","l");
    tleg2->AddEntry(p8_gs[0],"PYTHIA-8 g jets","l");
  }
  
  
  TLatex *ttitle = new TLatex(); ttitle->SetTextAlign(11); ttitle->SetTextSize(0.05);
  TLatex *slice = new TLatex();

  for (int i = 0; i < nBins; ++ i) {
    cws->cd(i+1);
    if (i == 0) {ldummy->Draw();} if (i == 1) {rdummy->Draw();}
    w_systs[i]->Draw("E3same");
    p6s[i]->Draw("Csame");
    p8s[i]->Draw("Csame");
    h7s[i]->Draw("Csame");
    PLs[i]->Draw("Csame");
    if (groom) {
      PL_gs[i]->Draw("Csame");
    }
    
    if (!groom) {
      p8_qs[i]->Draw("Csame");
      p8_gs[i]->Draw("Csame");
    }
    
    recos[i]->Draw("same");
    slice->DrawLatexNDC(0.5,0.77,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
    if (i == 0) {
      ttitle->DrawLatex(1,0.46, "pp 200 GeV run12 JP2");
      ttitle->DrawLatex(1,0.43, ("anti-k_{T}, R = 0."+to_string(radius)).c_str());
      ttitle->DrawLatex(1,0.4, ("Ch+Ne jets, |#eta| < 0."+to_string(10-radius)).c_str());
      tleg->Draw("same");
    }
    if (i == 1) { tleg2->Draw("same");}
  }

  cws->SaveAs(("~/jetmass2/plots/DNP_talk/"+fstart+"mass_result"+radstring+".pdf").c_str());

  return;
}
*/
//this file plots a 1x3 of the (groomed) mass result for various jet radii (and fixed pT)
void plot_result_forDNP (bool groom,int dummy1,int dummy2,int dummy3) {
  gROOT->ForceStyle();
  //gStyle->SetPalette(kPastel);
 
  string fstart = "";
  string hname = "m";
  string htitle = "M";
  if (groom) {
    fstart = "groomed_";
    hname = "mg";
    htitle = "M_{g}";
  }

  //open files from which to pull lots of curves
  TFile *funfold_R02 = new TFile(("~/jetmass2/out/unfold/"+fstart+"unfolded_R02.root").c_str(),"READ");
  TFile *fp6_R02 = new TFile("~/jetmass2/out/sim/hists/unmatched_hists_R02_bindropped.root","READ");
  TFile *fp8_R02 = new TFile("~/jetmass2/production/out/pythia/hists/pythia8_R02_undecayed_hists.root","READ");
  TFile *fh7_R02 = new TFile("~/jetmass/production/macros/hists/hists_allsim_lowzgremoved_R02.root","READ");//temporarily using the old files (they're the same, but will point to new ones later)
  TFile *funfold_R04 = new TFile(("~/jetmass2/out/unfold/"+fstart+"unfolded_R04.root").c_str(),"READ");
  TFile *fp6_R04 = new TFile("~/jetmass2/out/sim/hists/unmatched_hists_R04_bindropped.root","READ");
  TFile *fp8_R04 = new TFile("~/jetmass2/production/out/pythia/hists/pythia8_R04_undecayed_hists.root","READ");
  TFile *fh7_R04 = new TFile("~/jetmass/production/macros/hists/hists_allsim_lowzgremoved_R04.root","READ");//temporarily using the old files (they're the same, but will point to new ones later)
  TFile *funfold_R06 = new TFile(("~/jetmass2/out/unfold/"+fstart+"unfolded_R06.root").c_str(),"READ");
  TFile *fp6_R06 = new TFile("~/jetmass2/out/sim/hists/unmatched_hists_R06_bindropped.root","READ");
  TFile *fp8_R06 = new TFile("~/jetmass2/production/out/pythia/hists/pythia8_R06_undecayed_hists.root","READ");
  TFile *fh7_R06 = new TFile("~/jetmass/production/macros/hists/hists_allsim_lowzgremoved_R06.root","READ");//temporarily using the old files (they're the same, but will point to new ones later)
  

  //get histograms from the files
  vector<TH1D*> recos = {(TH1D*) ((TH1D*) funfold_R02->Get("nom_2"))->Clone("recos0"), (TH1D*) ((TH1D*) funfold_R04->Get("nom_2"))->Clone("recos1"), (TH1D*) ((TH1D*) funfold_R06->Get("nom_2"))->Clone("recos2")};
  vector<TH1D*> w_systs = {(TH1D*) ((TH1D*) funfold_R02->Get("w_systs_2"))->Clone("w_systs0"), (TH1D*) ((TH1D*) funfold_R04->Get("w_systs_2"))->Clone("w_systs1"), (TH1D*) ((TH1D*) funfold_R06->Get("w_systs_2"))->Clone("w_systs2")}; 
  vector<TH2D*> p6_2D = {(TH2D*) fp6_R02->Get(("PL_"+hname+"_v_pt").c_str()), (TH2D*) fp6_R04->Get(("PL_"+hname+"_v_pt").c_str()), (TH2D*) fp6_R06->Get(("PL_"+hname+"_v_pt").c_str())};
  vector<TH2D*> p8_2D = {(TH2D*) fp8_R02->Get((hname+"vpt").c_str()), (TH2D*) fp8_R04->Get((hname+"vpt").c_str()), (TH2D*) fp8_R06->Get((hname+"vpt").c_str())};
  vector<TH2D*> h7_2D = {(TH2D*) fh7_R02->Get((hname+"vpt_h7off").c_str()),(TH2D*) fh7_R04->Get((hname+"vpt_h7off").c_str()),(TH2D*) fh7_R06->Get((hname+"vpt_h7off").c_str())};
  vector<TH2D*> PL_2D = {(TH2D*) fp8_R02->Get("PLmvpt"),(TH2D*) fp8_R04->Get("PLmvpt"),(TH2D*) fp8_R06->Get("PLmvpt")};
  vector<TH2D*> PL_g_2D = {(TH2D*) fp8_R02->Get("PLmgvpt"),(TH2D*) fp8_R04->Get("PLmgvpt"),(TH2D*) fp8_R06->Get("PLmgvpt")};

  const int nPts = 1;
  const int nRads = 3;
  TAxis* p6_axis = p6_2D[0]->GetYaxis();
  double ranges[nPts + 1] = {(double) p6_axis->FindBin(30), (double) p6_axis->FindBin(40)};
  string pts[nPts + 1] = {"30","40"};
  
  vector<TH1D*> p6s, p8s, h7s, PLs, PL_gs; 
  for (int i = 0; i < p6_2D.size(); ++ i) {
    p6s.push_back(p6_2D[i]->ProjectionX(("p6"+to_string(i)).c_str(),p6_2D[i]->GetYaxis()->FindBin(30),p6_2D[i]->GetYaxis()->FindBin(40)));
    p8s.push_back(p8_2D[i]->ProjectionX(("p8"+to_string(i)).c_str(),p8_2D[i]->GetYaxis()->FindBin(30),p8_2D[i]->GetYaxis()->FindBin(40)));
    h7s.push_back(h7_2D[i]->ProjectionX(("h7"+to_string(i)).c_str(),h7_2D[i]->GetYaxis()->FindBin(30),h7_2D[i]->GetYaxis()->FindBin(40)));
    PLs.push_back(PL_2D[i]->ProjectionX(("PL"+to_string(i)).c_str(),PL_2D[i]->GetYaxis()->FindBin(30),PL_2D[i]->GetYaxis()->FindBin(40)));
    PL_gs.push_back(PL_g_2D[i]->ProjectionX(("PL_g"+to_string(i)).c_str(),PL_g_2D[i]->GetYaxis()->FindBin(30),PL_g_2D[i]->GetYaxis()->FindBin(40)));
  }

  TH1D* hdummyl = new TH1D("hdummyl","",14,0,13.999);
  TH1D* hdummym = new TH1D("hdummym","",14,0.001,13.999);
  TH1D* hdummyr = new TH1D("hdummyr","",14,0.001,13.999);
  hdummyl->GetYaxis()->SetRangeUser(0,0.4999);
  hdummym->GetYaxis()->SetRangeUser(0,0.4999);
  hdummyr->GetYaxis()->SetRangeUser(0,0.4999);
  //  gStyle->SetTitleOffset(1.02, "x");
  hdummyl->GetXaxis()->SetTitle((htitle+" [GeV/c^{2}]").c_str());
  hdummym->GetXaxis()->SetTitle((htitle+" [GeV/c^{2}]").c_str());
  hdummyr->GetXaxis()->SetTitle((htitle+" [GeV/c^{2}]").c_str());
  hdummyl->GetYaxis()->SetTitle(("1/N dN/d"+htitle).c_str());
  vector<TH1D*> dumvec = {hdummyl,hdummym,hdummyr};
  
  for (unsigned i = 0; i < nRads; ++ i) {
    Prettify1D(recos[i],kRed,kFullStar,2.5,kRed,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.5);
    Prettify1D(w_systs[i],kRed,kFullStar,0,kRed,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.5);
    w_systs[i]->SetFillColor(kRed - 10); w_systs[i]->SetFillStyle(1001);
    Prettify1DwLineStyle(p6s[i], kBlue, kSolid,2, (htitle+" [GeV/c^{2}]").c_str(), ("1/N dN/d"+htitle).c_str(),0,14,0,0.5);
    Prettify1DwLineStyle(p8s[i],kBlack,kSolid,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.5);
    Prettify1DwLineStyle(h7s[i],kMagenta,kSolid,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.5);
    Prettify1DwLineStyle(PLs[i],kBlack,kDotted,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.5);
    Prettify1DwLineStyle(PL_gs[i],kBlack,kDashed,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.5);
    
    p6s[i]->Sumw2(0);
    p8s[i]->Sumw2(0);
    h7s[i]->Sumw2(0);
    PLs[i]->Sumw2(0);
    PL_gs[i]->Sumw2(0);
  }
  
  TCanvas *cws = new TCanvas("cws","cws",800,400);
  DivideCanvas(cws,"0",nRads,1);

  TLegend *tleg = new TLegend(0.1,0.8,0.65,0.95); tleg->SetBorderSize(0);
  TLegend *tleg2 = new TLegend(0.025,0.8,0.6,0.95); tleg2->SetBorderSize(0);
  TH1D* for_legend = (TH1D*) w_systs[0]->Clone("for_legend"); for_legend->SetMarkerSize(2);
  tleg->AddEntry(for_legend,"STAR","pf");
  tleg->AddEntry(p6s[0],"PYTHIA-6","l");
  tleg2->AddEntry(h7s[0],"HERWIG-7","l");
  tleg->AddEntry(p8s[0],"PYTHIA-8","l");
  
  tleg2->AddEntry(PLs[0],"PYTHIA-8 parton jets","l");
  if (groom) {
    tleg2->AddEntry(PL_gs[0],"PYTHIA-8 SD parton jets","l");
  }
  
  TLatex *ttitle = new TLatex(); ttitle->SetTextAlign(11); ttitle->SetTextSize(0.05);
  TLatex *slice = new TLatex();
  
  for (int i = 0; i < nRads; ++ i) {
    cws->cd(i+1);
    dumvec[i]->Draw();
    w_systs[i]->Draw("E3same");
    p6s[i]->Draw("Csame");
    p8s[i]->Draw("Csame");
    h7s[i]->Draw("Csame");
    PLs[i]->Draw("Csame");
    if (groom) {
      PL_gs[i]->Draw("Csame");
    }
    recos[i]->Draw("same");
    if (i == 0) {
      slice->DrawLatexNDC(0.5,0.77,"30 < p_{T} < 40 GeV/c");
      slice->DrawLatexNDC(0.5,0.72,"R = 0.2");
      ttitle->DrawLatex(4.2,0.46, "pp 200 GeV run12 JP2");
      ttitle->DrawLatex(4.2,0.43, "anti-k_{T} jets, |#eta| < 1 - R");
      if (groom) {
	ttitle->DrawLatex(4.2,0.4,"SoftDrop z_{cut} = 0.1, #beta = 0");
      }
    }
    if (i == 1) {
      tleg2->Draw("same");
      slice->DrawLatexNDC(0.5,0.72,"R = 0.4");
    }
    if (i == 2) {
      slice->DrawLatexNDC(0.5,0.72,"R = 0.6");
      tleg->Draw("same");
    }
  }
  
  //cws->SaveAs(("~/jetmass2/plots/pp_paper/"+fstart+"mass_result_radscan.pdf").c_str());

  return;
}
