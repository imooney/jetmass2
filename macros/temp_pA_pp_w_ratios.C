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

//this function plots a 1x2 plot of mass as a function of pT, but with ratios of MC to data
void temp_pA_pp_w_ratios() {
  gROOT->ForceStyle(); //forces use of ~/rootlogon.C's style settings 
  
  const int nRows = 1; //0.4
  const int nCols = 2; //20-30,30-45
  
  //open files from which to pull curves                                                        
  TFile *fpAu = new TFile("~/jetmass2/out/data/data_pAuJP2.root","READ");
  //get histograms from the files
  TTree *tpAu = (TTree*) fpAu->Get("event");
  //  TTree *tpAu = (TTree*) tpAu_pre->CloneTree();//("tpAu");
  
  TFile *fpp = new TFile("~/jetmass2/out/data/data_ppJP2_R04.root","READ");
  TTree *tpp = (TTree*) fpp->Get("event");
  //TTree *tpp = (TTree*) tpp_pre->CloneTree();//("tpp");
  
  TH1D* ppEast_low = new TH1D("ppEast_low","",10,0,10);
  TH1D* ppEast_high = new TH1D("ppEast_high","",10,0,10);
  TH1D* ppWest_low = new TH1D("ppWest_low","",10,0,10);
  TH1D* ppWest_high = new TH1D("ppWest_high","",10,0,10);

  TH1D* pAulowEA_East_low = new TH1D("pAulowEA_East_low","",10,0,10);
  TH1D* pAulowEA_East_high = new TH1D("pAulowEA_East_high","",10,0,10);
  TH1D* pAulowEA_West_low = new TH1D("pAulowEA_West_low","",10,0,10);
  TH1D* pAulowEA_West_high = new TH1D("pAulowEA_West_high","",10,0,10);
  
  TH1D* pAuhighEA_East_low = new TH1D("pAuhighEA_East_low","",10,0,10);
  TH1D* pAuhighEA_East_high = new TH1D("pAuhighEA_East_high","",10,0,10);
  TH1D* pAuhighEA_West_low = new TH1D("pAuhighEA_West_low","",10,0,10);
  TH1D* pAuhighEA_West_high = new TH1D("pAuhighEA_West_high","",10,0,10);
  
  tpp->Draw("M>>ppEast_low","Eta < 0 && Pt > 20 && Pt < 25");
  tpp->Draw("M>>ppWest_low","Eta > 0 && Pt > 20 && Pt < 25");
  tpp->Draw("M>>ppEast_high","Eta < 0 && Pt > 30");
  tpp->Draw("M>>ppWest_high","Eta > 0 && Pt > 30");
 
  tpAu->Draw("M");
  tpAu->Draw("M>>pAulowEA_East_low","Eta < 0 && Pt > 20 && Pt < 25 && bbc_east_sum < 20000");
  tpAu->Draw("M>>pAulowEA_West_low","Eta > 0 && Pt > 20 && Pt < 25 && bbc_east_sum < 20000");
  tpAu->Draw("M>>pAulowEA_East_high","Eta < 0 && Pt > 30 && bbc_east_sum < 20000");
  tpAu->Draw("M>>pAulowEA_West_high","Eta > 0 && Pt > 30 && bbc_east_sum < 20000");
  
  tpAu->Draw("M>>pAuhighEA_East_low","Eta < 0 && Pt > 20 && Pt < 25 && bbc_east_sum > 30000");
  tpAu->Draw("M>>pAuhighEA_West_low","Eta > 0 && Pt > 20 && Pt < 25 && bbc_east_sum > 30000");
  tpAu->Draw("M>>pAuhighEA_East_high","Eta < 0 && Pt > 30 && bbc_east_sum > 30000");
  tpAu->Draw("M>>pAuhighEA_West_high","Eta > 0 && Pt > 30 && bbc_east_sum > 30000");
  
  TCanvas *clow = new TCanvas("clow","clow",800,500);
  DivideCanvas(clow,"0",nCols,nRows);

  TCanvas *chigh = new TCanvas("chigh","chigh",800,500);
  DivideCanvas(chigh,"0",nCols,nRows);
  TH1D* dummy_left = new TH1D("dummy_left","",10,0,9.999);
  TH1D* dummy_right = new TH1D("dummy_right","",10,0.001,9.999);

  TH1D* ppEW_low;
  TH1D* pAulowEA_EW_low;
  TH1D* pAulowEA_pp_E_low;
  TH1D* pAulowEA_pp_W_low;

  TH1D* ppEW_high;
  TH1D* pAulowEA_EW_high;
  TH1D* pAulowEA_pp_E_high;
  TH1D* pAulowEA_pp_W_high;

  // TH1D* pAulowEA_EW_low;
  TH1D* pAuhighEA_EW_low;
  TH1D* pAuhighEA_pAulowEA_E_low;
  TH1D* pAuhighEA_pAulowEA_W_low;

  //TH1D* pAulowEA_EW_high;
  TH1D* pAuhighEA_EW_high;
  TH1D* pAuhighEA_pAulowEA_E_high;
  TH1D* pAuhighEA_pAulowEA_W_high;

  //this will be the first "hist" drawn so needs to have correct axis properties
  dummy_left->GetYaxis()->SetRangeUser(0.5,1.499);
  dummy_right->GetYaxis()->SetRangeUser(0.5,1.499);
  
  dummy_left->GetYaxis()->SetTitle("ratios");//"E/W || pA/pp || high-EA/low-EA");
  dummy_left->GetYaxis()->SetTitleSize(0.15);
  dummy_left->GetYaxis()->SetTitleOffset(0.6);
  dummy_left->GetYaxis()->SetNdivisions(505);
  dummy_left->GetYaxis()->SetLabelSize(0.13);
  dummy_right->GetYaxis()->SetTitleSize(0.15);
  dummy_right->GetYaxis()->SetTitleOffset(0.6);
  dummy_right->GetYaxis()->SetNdivisions(505);
  dummy_right->GetYaxis()->SetLabelSize(0.13);

  dummy_left->GetXaxis()->SetTitle("M [GeV/c^{2}]");
  dummy_left->GetXaxis()->SetTitleSize(0.2);
  dummy_left->GetXaxis()->SetTitleOffset(0.8);
  dummy_left->GetXaxis()->SetLabelSize(0.13);
  dummy_left->GetXaxis()->SetNdivisions(505);
  dummy_right->GetXaxis()->SetTitle("M [GeV/c^{2}]");
  dummy_right->GetXaxis()->SetTitleSize(0.2);
  dummy_right->GetXaxis()->SetTitleOffset(0.8);
  dummy_right->GetXaxis()->SetLabelSize(0.13);
  dummy_right->GetXaxis()->SetNdivisions(505);
    
  
  Prettify1D(ppEast_low,kGreen+2,kOpenCircle,2,kGreen+2,"M [GeV/c^{2}]","1/N_{j}dN/dM",0,10,0,0.45);
  Prettify1D(ppEast_high,kGreen+2,kOpenCircle,2,kGreen+2,"M [GeV/c^{2}]","1/N_{j}dN/dM",0,10,0,0.45);
  Prettify1D(ppWest_low,kRed+2,kOpenSquare,2,kRed+2,"M [GeV/c^{2}]","1/N_{j}dN/dM",0,10,0,0.45);
  Prettify1D(ppWest_high,kRed+2,kOpenSquare,2,kRed+2,"M [GeV/c^{2}]","1/N_{j}dN/dM",0,10,0,0.45);
     
  Prettify1D(pAulowEA_East_low,kGreen+2,kFullCircle,2,kGreen+2,"M [GeV/c^{2}]","1/N_{j}dN/dM",0,10,0,0.45);
  Prettify1D(pAulowEA_East_high,kGreen+2,kFullCircle,2,kGreen+2,"M [GeV/c^{2}]","1/N_{j}dN/dM",0,10,0,0.45);
  Prettify1D(pAulowEA_West_low,kRed+2,kFullSquare,2,kRed+2,"M [GeV/c^{2}]","1/N_{j}dN/dM",0,10,0,0.45);
  Prettify1D(pAulowEA_West_high,kRed+2,kFullSquare,2,kRed+2,"M [GeV/c^{2}]","1/N_{j}dN/dM",0,10,0,0.45);
  
  Prettify1D(pAuhighEA_East_low,kGreen+2,kFullCross,2,kGreen+2,"M [GeV/c^{2}]","1/N_{j}dN/dM",0,10,0,0.45);
  Prettify1D(pAuhighEA_East_high,kGreen+2,kFullCross,2,kGreen+2,"M [GeV/c^{2}]","1/N_{j}dN/dM",0,10,0,0.45);
  Prettify1D(pAuhighEA_West_low,kRed+2,kFullDiamond,2,kRed+2,"M [GeV/c^{2}]","1/N_{j}dN/dM",0,10,0,0.45);
  Prettify1D(pAuhighEA_West_high,kRed+2,kFullDiamond,2,kRed+2,"M [GeV/c^{2}]","1/N_{j}dN/dM",0,10,0,0.45);
  
  
  ppEW_low = (TH1D*) ppEast_low->Clone("ppEW_low");
  ppEW_high = (TH1D*) ppEast_high->Clone("ppEW_high");
  
  pAulowEA_EW_low = (TH1D*) pAulowEA_East_low->Clone("pAulowEA_EW_low");
  pAulowEA_EW_high = (TH1D*) pAulowEA_East_high->Clone("pAulowEA_EW_high");

  pAulowEA_pp_E_low = (TH1D*) pAulowEA_East_low->Clone("pAulowEA_pp_E_low");
  pAulowEA_pp_E_high = (TH1D*) pAulowEA_East_high->Clone("pAulowEA_pp_E_high");

  pAulowEA_pp_W_low = (TH1D*) pAulowEA_West_low->Clone("pAulowEA_pp_W_low");
  pAulowEA_pp_W_high = (TH1D*) pAulowEA_West_high->Clone("pAulowEA_pp_W_high");

  pAuhighEA_EW_low = (TH1D*) pAuhighEA_East_low->Clone("pAuhighEA_EW_low");
  pAuhighEA_EW_high = (TH1D*) pAuhighEA_East_high->Clone("pAuhighEA_EW_high");

  pAuhighEA_pAulowEA_E_low = (TH1D*) pAuhighEA_East_low->Clone("pAuhighEA_pAulowEA_E_low");
  pAuhighEA_pAulowEA_E_high = (TH1D*) pAuhighEA_East_high->Clone("pAuhighEA_pAulowEA_E_high");
  
  pAuhighEA_pAulowEA_W_low = (TH1D*) pAuhighEA_West_low->Clone("pAuhighEA_pAulowEA_W_low");
  pAuhighEA_pAulowEA_W_high = (TH1D*) pAuhighEA_West_high->Clone("pAuhighEA_pAulowEA_W_high");

  
  ppEW_low->Divide(ppWest_low);
  ppEW_high->Divide(ppWest_high);
  
  pAulowEA_EW_low->Divide(pAulowEA_West_low);
  pAulowEA_EW_high->Divide(pAulowEA_West_high);
  
  pAulowEA_pp_E_low->Divide(ppEast_low);
  pAulowEA_pp_E_high->Divide(ppEast_high);
  
  pAulowEA_pp_W_low->Divide(ppWest_low);
  pAulowEA_pp_W_high->Divide(ppWest_high);
  
  pAuhighEA_EW_low->Divide(pAuhighEA_West_low);
  pAuhighEA_EW_high->Divide(pAuhighEA_West_high);
  
  pAuhighEA_pAulowEA_E_low->Divide(pAulowEA_East_low);
  pAuhighEA_pAulowEA_E_high->Divide(pAulowEA_East_high);
  
  pAuhighEA_pAulowEA_W_low->Divide(pAulowEA_West_low);
  pAuhighEA_pAulowEA_W_high->Divide(pAulowEA_West_high);
  
  pAulowEA_EW_low->GetYaxis()->SetRangeUser(0.5,1.5);
  pAulowEA_EW_high->GetYaxis()->SetRangeUser(0.5,1.5);
  
  
  TLegend *tleg = new TLegend(0.525,0.45,0.85,0.6); tleg->SetBorderSize(0);
  tleg->AddEntry(pAulowEA_East_low,"pAu low-EA east","p");
  tleg->AddEntry(pAulowEA_West_low,"pAu low-EA west","p");
  tleg->SetTextSize(0.06);
  TLegend *tleg2 = new TLegend(0.6,0.45,0.9,0.6); tleg2->SetBorderSize(0);
  tleg2->AddEntry(ppEast_low,"pp east","p");
  tleg2->AddEntry(ppWest_low,"pp west","p");
  tleg2->SetTextSize(0.06);
  TLegend *tleg3 = new TLegend(0.525,0.45,0.85,0.6); tleg3->SetBorderSize(0);
  tleg3->AddEntry(pAuhighEA_East_low,"pAu high-EA east","p");
  tleg3->AddEntry(pAuhighEA_West_low,"pAu high-EA west","p");
  tleg3->SetTextSize(0.06);
  
  TLatex *ttitle = new TLatex(); TLatex *slice = new TLatex();
  ttitle->SetTextSize(0.07); slice->SetTextSize(0.07);
  
  TLine *one = new TLine(0,1,10,1); one->SetLineStyle(kDashed);
  
  //PANEL 1
  clow->cd(1);
  TPad *padup0 = new TPad("padup0","padup0",0,0.4,1,1.0);
  padup0->SetBottomMargin(0);
  padup0->SetTopMargin(0);
  padup0->SetRightMargin(0);
  padup0->Draw();
  padup0->cd();
  pAulowEA_West_low->Draw("same");
  pAulowEA_East_low->Draw("same");
  ppWest_low->Draw("same");
  ppEast_low->Draw("same");
  tleg->Draw("same");
  
  slice->DrawLatexNDC(0.55,0.75,"20 < p_{T} < 25 GeV/c");
  
  ttitle->DrawLatex(2.2,0.41, "anti-k_{T}, R = 0.4, |#eta_{jet}| < 1-R");  
  
  clow->cd(1);
  TPad *paddown0 = new TPad("paddown0","paddown0",0,0.05,1,0.4);
  paddown0->SetTopMargin(0);
  paddown0->SetBottomMargin(0.4);
  paddown0->SetRightMargin(0);
  paddown0->SetLeftMargin(0);
  paddown0->SetLeftMargin(0.2);
  paddown0->Draw();
  paddown0->cd();
  dummy_left->Draw();
  ppEW_low->SetMarkerColor(kOrange);
  pAulowEA_EW_low->SetMarkerColor(kOrange);
  pAulowEA_EW_low->Draw("E1same");
  pAulowEA_pp_E_low->Draw("E1same");
  pAulowEA_pp_W_low->Draw("E1same");
  ppEW_low->Draw("E1same");
  one->Draw("same");
  
  //PANEL 2
  clow->cd(2);
  TPad *padup1 = new TPad("padup1","padup1",0,0.4,1,1.0);
  padup1->SetBottomMargin(0);
  padup1->SetTopMargin(0);
  padup1->SetLeftMargin(0);
  padup1->SetRightMargin(0);
  padup1->Draw();
  padup1->cd();
  pAulowEA_West_high->Draw("same");
  pAulowEA_East_high->Draw("same");
  ppWest_high->Draw("same");
  ppEast_high->Draw("same");
  tleg2->Draw("same");
  
  slice->DrawLatexNDC(0.55,0.75,"p_{T} > 30 GeV/c");
  
  clow->cd(2);
  TPad *paddown1 = new TPad("paddown1","paddown1",0,0.05,1,0.4);
  paddown1->SetTopMargin(0);
  paddown1->SetBottomMargin(0.4);
  paddown1->SetRightMargin(0);
  paddown1->SetLeftMargin(0);
  paddown1->Draw();
  paddown1->cd();
  dummy_right->Draw();
  ppEW_high->SetMarkerColor(kOrange);
  pAulowEA_EW_high->SetMarkerColor(kOrange);
  pAulowEA_EW_high->Draw("E1same");
  ppEW_high->Draw("E1same");
  pAulowEA_pp_E_high->Draw("E1same");
  pAulowEA_pp_W_high->Draw("E1same");  
  one->Draw("same");
  

  ///!!!!!!!!!!!!!!!!!!!///


  //PANEL 1
  chigh->cd(1);
  TPad *padup2 = new TPad("padup2","padup2",0,0.4,1,1.0);
  padup2->SetBottomMargin(0);
  padup2->SetTopMargin(0);
  padup2->SetRightMargin(0);
  padup2->Draw();
  padup2->cd();
  pAuhighEA_West_low->Draw("same");
  pAuhighEA_East_low->Draw("same");
  pAulowEA_West_low->Draw("same");
  pAulowEA_East_low->Draw("same");
  tleg->Draw("same");
  
  slice->DrawLatexNDC(0.55,0.75,"20 < p_{T} < 25 GeV/c");
  
  ttitle->DrawLatex(2.2,0.41, "anti-k_{T}, R = 0.4, |#eta_{jet}| < 1-R");  
  
  chigh->cd(1);
  TPad *paddown2 = new TPad("paddown2","paddown2",0,0.05,1,0.4);
  paddown2->SetTopMargin(0);
  paddown2->SetBottomMargin(0.4);
  paddown2->SetRightMargin(0);
  paddown2->SetLeftMargin(0);
  paddown2->SetLeftMargin(0.2);
  paddown2->Draw();
  paddown2->cd();
  dummy_left->Draw();
  pAuhighEA_EW_low->SetMarkerColor(kOrange);
  pAuhighEA_EW_low->Draw("E1same");
  pAuhighEA_pAulowEA_E_low->Draw("E1same");
  pAuhighEA_pAulowEA_W_low->Draw("E1same");
  pAulowEA_EW_low->Draw("E1same");
  one->Draw("same");
  
  //PANEL 2
  chigh->cd(2);
  TPad *padup3 = new TPad("padup3","padup3",0,0.4,1,1.0);
  padup3->SetBottomMargin(0);
  padup3->SetTopMargin(0);
  padup3->SetLeftMargin(0);
  padup3->SetRightMargin(0);
  padup3->Draw();
  padup3->cd();
  pAuhighEA_West_high->Draw("same");
  pAuhighEA_East_high->Draw("same");
  pAulowEA_West_high->Draw("same");
  pAulowEA_East_high->Draw("same");
  tleg3->Draw("same");
  
  slice->DrawLatexNDC(0.55,0.75,"p_{T} > 30 GeV/c");
  
  chigh->cd(2);
  TPad *paddown3 = new TPad("paddown3","paddown3",0,0.05,1,0.4);
  paddown3->SetTopMargin(0);
  paddown3->SetBottomMargin(0.4);
  paddown3->SetRightMargin(0);
  paddown3->SetLeftMargin(0);
  paddown3->Draw();
  paddown3->cd();
  dummy_right->Draw();
  pAuhighEA_EW_high->SetMarkerColor(kOrange);
  pAuhighEA_EW_high->Draw("E1same");
  pAulowEA_EW_high->Draw("E1same");
  pAuhighEA_pAulowEA_E_high->Draw("E1same");
  pAuhighEA_pAulowEA_W_high->Draw("E1same");  
  one->Draw("same");

  


  
  return;
}
