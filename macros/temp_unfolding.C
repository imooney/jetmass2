#include "RooUnfoldResponse.h"
#include "Plots_old.h"
#include <algorithm>

using namespace std;

void temp_unfolding () {
  gStyle->SetPalette(kPastel);
  
  TFile *fdat = new TFile("~/jetmass2/out/data/hists/data_hists_ppJP2_bindropped.root","READ");
  TFile *fres = new TFile("~/jetmass2/out/sim/hists/matched_hists_bindropped.root","READ");
  TFile *fold = new TFile("~/jetmass2/out/unfold/unfolded_v1.root","READ");

  vector<TH1D*> m_pt_unf_old = {(TH1D*) fold->Get("nom_0"),(TH1D*) fold->Get("nom_1"),(TH1D*) fold->Get("nom_2")};
  vector<TH1D*> dats_old = {(TH1D*) fold->Get("m_v_pt_dx23"),(TH1D*) fold->Get("m_v_pt_dx34"),(TH1D*) fold->Get("m_v_pt_dx46")};
  TH2D* m_pt_dat = (TH2D*) fdat->Get("m_v_pt");
  RooUnfoldResponse *rnom = (RooUnfoldResponse*) fres->Get("m_pt_response");
  
  RooUnfoldBayes *unfold_nom = new RooUnfoldBayes(rnom, m_pt_dat, 4, false, "unfold_nom","");  
  TH2D *reco_nom = (TH2D*) unfold_nom->Hreco((RooUnfold::ErrorTreatment) 3);

  
  const int nBins = 3;
  TAxis* reco_axis = reco_nom->GetYaxis(); TAxis* det_axis = m_pt_dat->GetYaxis();  
  double ranges[nBins + 1] = {(double) reco_axis->FindBin(20), (double) reco_axis->FindBin(25), (double) reco_axis->FindBin(30), (double) reco_axis->FindBin(40)};//3,4,5,6,8,12};
  double ranges_d[nBins + 1] = {(double) det_axis->FindBin(20), (double) det_axis->FindBin(25), (double) det_axis->FindBin(30), (double) det_axis->FindBin(40)};//1,2,3,4,6,10};
  string pts[nBins + 1] = {"20","25","30","40"};
  
  vector<TH1D*> reco_noms = Projection2D(reco_nom,nBins,ranges,"x");
  vector<TH1D*> dats = Projection2D(m_pt_dat,nBins,ranges_d,"x");
  
  for (int i = 0; i < nBins; ++ i) {
    Prettify1D(reco_noms[i],kRed,kFullStar,4,kRed,"M [GeV/c^{2}]","1/N dN/dM",0,14,/*0.9,1.1*/0,0.5);
    Prettify1D(m_pt_unf_old[i],kRed,kOpenStar,4,kRed,"M [GeV/c^{2}]","1/N dN/dM",0,14,/*0.9,1.1*/0,0.5);
    Prettify1D(dats[i], kBlack, kFullStar, 4, kBlack, "M [GeV/c^{2}]", "1/N dN/dM",0,14,/*0.9,1.1*/0,0.5);
    Prettify1D(dats_old[i], kBlack, kOpenStar, 4, kBlack, "M [GeV/c^{2}]", "1/N dN/dM",0,14,/*0.9,1.1*/0,0.5);
  }
  /*
  for (int i = 0; i < nBins; ++ i) {
    reco_noms[i]->Divide(m_pt_unf_old[i]);
  }
  */
  //scaling errors
  TFile *fstats = new TFile("~/jetmass2/out/sim/stat_err_scaling.root","READ");
  TH1D* scalefactors = (TH1D*) fstats->Get("hratio");
  
  for (int i = 0; i < nBins; ++ i) {
    cout << "i = " << i << endl << endl;
    for (int j = 1; j <= reco_noms[i]->GetNbinsX(); ++ j) {
      double scaling = -1;
      TAxis *scalex = scalefactors->GetXaxis();
      //will access the scalefactors histogram for pt values 20, 25, 30, 35, so if you change the ranges you show, change the FindBin calls, too.
      if (i == 0) {scaling = scalefactors->GetBinContent(scalex->FindBin(20));/*old: 1.122;*/} //these numbers are calculated using the bin content of the ratio of gen. matched spectrum to gen. inclusive (unmatched). See macros/stat_err_scaling.cxx.
      if (i == 1) {scaling = scalefactors->GetBinContent(scalex->FindBin(25));/*old: 1.082;*/}
      if (i == 2) {scaling = max(scalefactors->GetBinContent(scalex->FindBin(30)),scalefactors->GetBinContent(scalex->FindBin(35)));/*old: 1.062;*/}
      double binerror = reco_noms[i]->GetBinError(j);
      reco_noms[i]->SetBinError(j,(double) binerror*scaling);
      cout << "bin " << j << " error: " << binerror << endl;
      cout << "corresponding raw data error: " << dats[i]->GetBinError(j) << endl;
      cout << "scaling bin " << j << " error by: " << scaling << endl;
      cout << "scaled error: " << binerror*scaling << endl;
    }
  }
  
  TCanvas *cws = new TCanvas("cws","cws",1200,500);
  DivideCanvas(cws,"0",3,1);

  TLegend *twsysts2 = new TLegend(0.5,0.5,0.75,0.7); twsysts2->SetBorderSize(0);
  twsysts2->AddEntry(dats[0],"Raw data -new","p");
  twsysts2->AddEntry(dats_old[0],"Raw data -old","p");
  twsysts2->AddEntry(reco_noms[0],"Unfolded data -new","p");
  twsysts2->AddEntry(m_pt_unf_old[0],"Unfolded data -old","p");
  twsysts2->SetTextSize(0.04);

  TLatex *ttitle = new TLatex(); ttitle->SetTextAlign(11); ttitle->SetTextSize(0.05);
  TLatex *slice = new TLatex();
  
  TH1D* hdummycws = new TH1D("hdummycws",";;1/N dN/dM",1,0,14);
  hdummycws->GetYaxis()->SetRangeUser(0,0.4);
  hdummycws->GetYaxis()->SetTitleSize(0.06);

  for (int i = 0; i < nBins; ++ i) {

    cws->cd(i+1); dats[i]->Draw("9"); dats_old[i]->Draw("9same"); m_pt_unf_old[i]->Draw("9same"); reco_noms[i]->Draw("same9"); slice->DrawLatexNDC(0.5,0.77,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
    if (i == 0) {ttitle->DrawLatex(4.2,0.46, "pp 200 GeV run12 JP2");ttitle->DrawLatex(4.2,0.43, "anti-k_{T}, R = 0.4");ttitle->DrawLatex(4.2,0.4, "Ch+Ne jets, |#eta| < 0.6");}
    if (i == 2) { twsysts2->Draw("same9");}
  }

  //  cws->SaveAs("~/jetmass2/plots/compare_code_v1_v2_w_bugs.pdf");

  return;
}
