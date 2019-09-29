//Isaac Mooney, WSU, September 2019
//This file plots the analysis note plots.

#include "Plots_old.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

using namespace std;

void rawdata(TFile *fdat, TFile *fsim ) {
  //get 2Ds
  TH2D* hdat = (TH2D*) fdat->Get("m_v_pt");
  TH2D* hge = (TH2D*) fsim->Get("m_v_pt");
  TH2D* hpy = (TH2D*) fsim->Get("PL_m_v_pt");
  TH2D* hdatg = (TH2D*) fdat->Get("mg_v_pt");
  TH2D* hgeg = (TH2D*) fsim->Get("mg_v_pt");
  TH2D* hpyg = (TH2D*) fsim->Get("PL_mg_v_pt");

  //project onto 1Ds
  TH1D* hdat_1D = (TH1D*) hdat->ProjectionX("hdat_1D",hdat->GetYaxis()->FindBin(20),hdat->GetYaxis()->FindBin(25)-1);
  TH1D* hge_1D = (TH1D*) hge->ProjectionX("hge_1D",hge->GetYaxis()->FindBin(20),hge->GetYaxis()->FindBin(25)-1);
  TH1D* hpy_1D = (TH1D*) hpy->ProjectionX("hpy_1D",hpy->GetYaxis()->FindBin(20),hpy->GetYaxis()->FindBin(25)-1);
  TH1D* hdatg_1D = (TH1D*) hdatg->ProjectionX("hdatg_1D",hdatg->GetYaxis()->FindBin(20),hdatg->GetYaxis()->FindBin(25)-1);
  TH1D* hgeg_1D = (TH1D*) hgeg->ProjectionX("hgeg_1D",hgeg->GetYaxis()->FindBin(20),hgeg->GetYaxis()->FindBin(25)-1);
  TH1D* hpyg_1D = (TH1D*) hpyg->ProjectionX("hpyg_1D",hpyg->GetYaxis()->FindBin(20),hpyg->GetYaxis()->FindBin(25)-1);
  
  //define dummies and prettify them
  TH1D* hdummy = new TH1D("hdummy","",10,0,9.999);
  TH1D* hdummyg = new TH1D("hdummyg","",10,0.001,10);
  hdummy->GetYaxis()->SetRangeUser(0,0.5);
  hdummyg->GetYaxis()->SetRangeUser(0,0.5);
  hdummy->GetXaxis()->SetTitle("M [GeV/c^{2}]");
  hdummyg->GetXaxis()->SetTitle("M_{g} [GeV/c^{2}]");
  hdummy->GetYaxis()->SetTitle("1/N dN/dM");
  hdummyg->GetYaxis()->SetTitle("1/N dN/dM_{g}");
  
  
  //prettify
  Prettify1D(hdat_1D,kBlack,kFullStar,2.5,kBlack,"M [GeV/c^{2}]","1/N dN/dM",0,10,0,0.5);
  Prettify1D(hge_1D,kBlue,kOpenCircle,2,kBlue,"M [GeV/c^{2}]","1/N dN/dM",0,10,0,0.5);
  Prettify1DwLineStyle(hpy_1D, kBlue, kSolid,2, "M [GeV/c^{2}]", "1/N dN/dM",0,10,0,0.5);
  Prettify1D(hdatg_1D,kBlack,kFullStar,2.5,kBlack,"M_{g} [GeV/c^{2}]","1/N dN/dM_{g}",0,10,0,0.5);
  Prettify1D(hgeg_1D,kBlue,kOpenCircle,2,kBlue,"M_{g} [GeV/c^{2}]","1/N dN/dM_{g}",0,10,0,0.5);
  Prettify1DwLineStyle(hpyg_1D, kBlue, kSolid,2, "M_{g} [GeV/c^{2}]", "1/N dN/dM_{g}",0,10,0,0.5); 

  hpy_1D->Sumw2(0);
  hpyg_1D->Sumw2(0);
  
  hdummyg->GetYaxis()->SetTitleOffset(1.1);
  hdummyg->GetXaxis()->SetTitleSize(0.07);
  hdummy->GetXaxis()->SetTitleSize(0.07);
  
  
  //legends
  TLatex *lleft = new TLatex();
  TLatex *lright = new TLatex();
  TLegend *t = new TLegend(0.3,0.5,0.7,0.7);
  t->SetBorderSize(0);
  t->AddEntry(hpy_1D,"PYTHIA-6","l");
  t->AddEntry(hge_1D,"PYTHIA-6+GEANT","p");
  t->AddEntry(hdat_1D,"Raw STAR data","p");
  t->SetTextSize(0.045);
  
  
  //draw
  TCanvas *craw = new TCanvas("craw","craw",800,500);
  craw->cd();
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
  hdummy->Draw();//so the axes will be slightly shifted so 0/10 don't overlap
  hge_1D->Draw("same");
  //axis->Draw("same");
  hdat_1D->Draw("same");
  hpy_1D->Draw("Csame");
  lleft->DrawLatex(0.5,0.45,"pp #sqrt{s_{NN}} = 200 GeV, JP2");
  lleft->DrawLatex(0.5,0.41,"anti-k_{t} R = 0.4 jets, |#eta|<1-R_{jet}");
  pright->cd();
  hdummyg->Draw("Y+");
  hgeg_1D->Draw("same");
  hdatg_1D->Draw("same");
  hpyg_1D->Draw("Csame");
  lright->DrawLatex(0.5,0.45,"20<p_{T}^{jet}<25 GeV/c");
  lright->DrawLatex(0.5,0.41,"SoftDrop z_{cut} = 0.1, #beta = 0");
  t->Draw("same");

  cout << "Finished drawing raw data plots" << endl;
  
  craw->SaveAs("~/jetmass2/plots/pp_AN/raw_mass.pdf");

  return;
}

//plots the mass resolution for different pT ranges on the same panel
void resolution(TFile* fp6match) {
    //get 2Ds
  TH2D* hmres = (TH2D*) fp6match->Get("ratioMvGePt");
  TH2D* hmresg = (TH2D*) fp6match->Get("ratioMgvGePt");

  //project onto 1Ds
  const int nCols = 3;
  TAxis *ptaxis = (TAxis*) hmres->GetYaxis();
  double ranges[nCols + 1] = {(double) ptaxis->FindBin(20), (double) ptaxis->FindBin(25), (double) ptaxis->FindBin(30), (double) ptaxis->FindBin(40)};
  vector<TH1D*> hmres_1D = Projection2D(hmres,nCols,ranges,"x");
  vector<TH1D*> hmresg_1D = Projection2D(hmresg,nCols,ranges,"x");

  //define dummies and prettify them
  TH1D* hdummy = new TH1D("hdummy","",20,0,1.999);
  TH1D* hdummyg = new TH1D("hdummyg","",20,0.001,2);
  hdummy->GetYaxis()->SetRangeUser(0.0001,0.25);
  hdummyg->GetYaxis()->SetRangeUser(0.0001,0.25);
  hdummy->GetXaxis()->SetTitle("r = M^{det} / M^{part}");
  hdummyg->GetXaxis()->SetTitle("r = M_{g}^{det} / M_{g}^{part}");
  hdummy->GetYaxis()->SetTitle("1/N_{pair} dN/dr");
  //hdummyg->GetYaxis()->SetTitle("a.u.");
  
  for (int i = 0; i < nCols; ++ i) {
    hmres_1D[i]->Scale(1/(double)hmres_1D[i]->Integral("w"));
    hmresg_1D[i]->Scale(1/(double)hmresg_1D[i]->Integral("w"));
  }
  
  //prettify
  Prettify1D(hmres_1D[0],kBlue,kOpenCircle,2,kBlue,"M^{det} / M^{part}","a.u.",0,2,0.0001,0.25);
  Prettify1D(hmres_1D[1],kMagenta+1,kOpenSquare,2,kMagenta+1,"M^{det} / M^{part}","a.u.",0,2,0.0001,0.25);
  Prettify1D(hmres_1D[2],kGreen+2, kOpenCross,2,kGreen+2, "M^{det} / M^{part}", "a.u.",0,2,0.0001,0.25);
  Prettify1D(hmresg_1D[0],kBlue,kOpenCircle,2,kBlue,"M_{g}^{det} / M_{g}^{part}","a.u.",0,2,0.0001,0.25);
  Prettify1D(hmresg_1D[1],kMagenta+1,kOpenSquare,2,kBlue,"M_{g}^{det} / M_{g}^{part}","a.u.",0,2,0.0001,0.25);
  Prettify1D(hmresg_1D[2],kGreen+2, kOpenCross,2,kGreen+2, "M_{g}^{det} / M_{g}^{part}", "a.u.",0,2,0.0001,0.25); 

  hdummy->GetXaxis()->SetNdivisions(505);
  hdummyg->GetXaxis()->SetNdivisions(505);
  hdummy->GetYaxis()->SetTitleOffset(0.7);
  hdummyg->GetYaxis()->SetTitleOffset(1.1);
  hdummyg->GetXaxis()->SetTitleSize(0.07);
  hdummy->GetXaxis()->SetTitleSize(0.07);
  
  
  //legends
  TLatex *lleft = new TLatex();
  TLatex *lright = new TLatex();
  TLegend *t = new TLegend(0.55,0.45,0.75,0.65);
  t->SetBorderSize(0);
  t->AddEntry(hmres_1D[0],"p_{T}^{det.} #in (20,25) GeV/c","p");
  t->AddEntry(hmres_1D[1],"p_{T}^{det.} #in (25,30) GeV/c","p");
  t->AddEntry(hmres_1D[2],"p_{T}^{det.} #in (30,40) GeV/c","p");
  t->SetTextSize(0.035);
  
  
  //draw
  TCanvas *cres = new TCanvas("cres","cres",800,500);
  cres->cd();
  /*  TPad *pleft = new TPad("pleft","",0,0,0.5,1);
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
  */hdummy->Draw();//so the axes will be slightly shifted so 0/10 don't overlap
  hmres_1D[0]->Draw("same");
  hmres_1D[1]->Draw("same");
  hmres_1D[2]->Draw("Csame");
  lleft->DrawLatex(0.1,0.23,"pp #sqrt{s} = 200 GeV, 2012");
  lleft->DrawLatex(0.1,0.21,"PYTHIA-6+GEANT");
  lright->DrawLatex(0.1,0.19,"anti-k_{t} R = 0.4 jets, |#eta|<1-R_{jet}");
  t->Draw("same");
  /* pright->cd();
  hdummyg->Draw();
  hmresg_1D[0]->Draw("same");
  hmresg_1D[1]->Draw("same");
  hmresg_1D[2]->Draw("Csame");
  lright->DrawLatex(0.1,0.23,"SoftDrop z_{cut} = 0.1, #beta = 0");
  */
  cout << "Finished drawing resolution plots" << endl;
  
  cres->SaveAs("~/jetmass2/plots/DNP_talk/mass_resolution_v_detpt.pdf");
  
  return;
}

void closure(TFile *fclos) {
  TH1D* opp_m1D = (TH1D*) fclos->Get("unfold_opp_m1D");
  TH1D* opp_pt1D = (TH1D*) fclos->Get("unfold_opp_pt1D");
  TH1D* opp_m2030 = (TH1D*) fclos->Get("opp_m2030");
  TH1D* opp_m3045 = (TH1D*) fclos->Get("opp_m3045");
  TH1D* same_m1D = (TH1D*) fclos->Get("unfold_same_m1D");
  TH1D* same_pt1D = (TH1D*) fclos->Get("unfold_same_pt1D");
  TH1D* same_m2030 = (TH1D*) fclos->Get("same_m2030");
  TH1D* same_m3045 = (TH1D*) fclos->Get("same_m3045");
  TH1D* ldummy_m = new TH1D("ldummy_m","",10,0,9.999);
  TH1D* rdummy_m = new TH1D("rdummy_m","",10,0.0001,9.999);

  ldummy_m->GetXaxis()->SetTitle("M [GeV/c^{2}]");
  rdummy_m->GetXaxis()->SetTitle("M [GeV/c^{2}]");
  ldummy_m->GetYaxis()->SetTitle("reco / gen");
  rdummy_m->GetYaxis()->SetTitle("reco / gen");
  //  ldummy_m->SetTitle("STAR Simulation");
  ldummy_m->GetYaxis()->SetRangeUser(0.5,1.4999);
  rdummy_m->GetYaxis()->SetRangeUser(0.5,1.4999);
  
  opp_m1D->SetTitle("");
  same_m1D->SetTitle("");
  opp_pt1D->SetTitle("");
  same_pt1D->SetTitle("");
  opp_m2030->SetTitle("");
  same_m2030->SetTitle("");
  opp_m3045->SetTitle("");
  same_m3045->SetTitle("");
  same_m1D->GetYaxis()->SetTitleOffset(0.5);
  same_pt1D->GetYaxis()->SetTitleOffset(0.5);

  Prettify1D(opp_m1D,kBlue,kOpenSquare,2,kBlue,"M [GeV/c^{2}]","reco / gen",0,10,0.5,1.5);
  Prettify1D(opp_pt1D,kBlue,kOpenSquare,2,kBlue,"p_{T} [GeV/c]","reco / gen",5,60,0.5,1.5);
  Prettify1D(opp_m2030,kBlue,kOpenSquare,2,kBlue,"M [GeV/c^{2}]","reco / gen",0,10,0.5,1.5);
  Prettify1D(opp_m3045,kBlue,kOpenSquare,2,kBlue,"M [GeV/c^{2}]","reco / gen",0,10,0.5,1.5);
  Prettify1D(same_m1D,kMagenta,kFullCircle,2,kMagenta,"M [GeV/c^{2}]","reco / gen",0,10,0.5,1.5);
  Prettify1D(same_pt1D,kMagenta,kFullCircle,2,kMagenta,"p_{T} [GeV/c]","reco / gen",5,60,0.5,1.5);
  Prettify1D(same_m2030,kMagenta,kFullCircle,2,kMagenta,"M [GeV/c^{2}]","reco / gen",0,10,0.5,1.5);
  Prettify1D(same_m3045,kMagenta,kFullCircle,2,kMagenta,"M [GeV/c^{2}]","reco / gen",0,10,0.5,1.5);
  
  TLegend *t = new TLegend(0.3,0.65,0.6,0.8);
  t->SetBorderSize(0);
  t->AddEntry(same_m1D,"Training (4 iter.)","p");
  t->AddEntry(opp_m1D,"Validation (4 iter.)","p");
  TLegend *t2 = new TLegend(0.3,0.2,0.6,0.35);
  t2->SetBorderSize(0);
  t2->AddEntry(same_m1D,"Training (4 iter.)","p");
  t2->AddEntry(opp_m1D,"Validation (4 iter.)","p");

  
  TLatex *ptrange = new TLatex();
  TLatex *ltitle = new TLatex();
  
  TLine* one_m = new TLine(0,1,10,1);
  TLine* one_pt = new TLine(5,1,60,1);
  TLine* plus_m = new TLine(0,1.05,10,1.05);
  TLine* plus_pt = new TLine(5,1.05,60,1.05);
  TLine* minus_m = new TLine(0,0.95,10,0.95);
  TLine* minus_pt = new TLine(5,0.95,60,0.95);
  one_m->SetLineStyle(kSolid);
  one_pt->SetLineStyle(kSolid);
  plus_m->SetLineStyle(kDashed);
  plus_pt->SetLineStyle(kDashed);
  minus_m->SetLineStyle(kDashed);
  minus_pt->SetLineStyle(kDashed);
  
  TCanvas *cm1D = new TCanvas("cm1D","cm1D",800,500);
  TCanvas *cpt1D = new TCanvas("cpt1D","cpt1D",800,500);
  TCanvas *cm2D = new TCanvas("cm2D","cm2D",800,500);
  cm2D->Divide(2,1,0,0);
  /*
  cm1D->cd();
  same_m1D->Draw();
  opp_m1D->Draw("same");
  t->Draw("same");
  one_m->Draw("same");
  plus_m->Draw("same");
  minus_m->Draw("same");
  ltitle->DrawLatex(0.5,0.78,"pp 200 GeV run12");
  ltitle->DrawLatex(0.5,0.7,"anti-k_{T}, R = 0.4");
  ltitle->DrawLatex(0.5,0.62,"Ch+Ne jets, |#eta| < 0.6");
  ltitle->DrawLatex(0.5,0.54,"1D Bayesian unfolding closure");
  cm1D->SaveAs("~/jetmass2/plots/pp_AN/mass_closure_1D.pdf");
  
  cpt1D->cd();
  same_pt1D->Draw();
  opp_pt1D->Draw("same");
  t->Draw("same");
  one_pt->Draw("same");
  plus_pt->Draw("same");
  minus_pt->Draw("same");
  ltitle->DrawLatex(7.5,0.78,"pp 200 GeV run12");
  ltitle->DrawLatex(7.5,0.7,"anti-k_{T}, R = 0.4");
  ltitle->DrawLatex(7.5,0.62,"Ch+Ne jets, |#eta| < 0.6");
  ltitle->DrawLatex(7.5,0.54,"1D Bayesian unfolding closure");
  cpt1D->SaveAs("~/jetmass2/plots/pp_AN/pt_closure_1D.pdf");
*/  
  cm2D->cd();
  cm2D->cd(1);
  ldummy_m->Draw();
  same_m2030->Draw("same");
  opp_m2030->Draw("same");
  one_m->Draw("same");
  plus_m->Draw("same");
  minus_m->Draw("same");
  ptrange->DrawLatex(3,1.3,"20 < p_{T} < 30 GeV/c");
  
  cm2D->cd(2);
  rdummy_m->Draw();
  same_m3045->Draw("same");
  opp_m3045->Draw("same");
  t2->Draw("same");
  one_m->Draw("same");
  plus_m->Draw("same");
  minus_m->Draw("same");
  ptrange->DrawLatex(3,1.3,"30 < p_{T} < 45 GeV/c");
  
  cm2D->SaveAs("~/jetmass2/plots/pp_AN/m_pt_closure_2D.pdf");
  
  return;
}

//There is basically a duplicate of this function in unfold.cxx. This is just for plotting things for the analysis note. The other one is where I get my actual numbers for the systematics - not that they should be different in any way.
void systematics (TFile *fres, TFile *fdat, string hname, string htitle, string fstart) {
  gROOT->ForceStyle();
  gStyle->SetPalette(kPastel);
  
  RooUnfoldResponse *rnom = (RooUnfoldResponse*) fres->Get((hname+"_pt_res_nom").c_str());                                                              
  RooUnfoldResponse *rTS = (RooUnfoldResponse*) fres->Get((hname+"_pt_res_TS").c_str());                                                           
  RooUnfoldResponse *rTU = (RooUnfoldResponse*) fres->Get((hname+"_pt_res_TU").c_str());                                                          
  RooUnfoldResponse *rHC50 = (RooUnfoldResponse*) fres->Get((hname+"_pt_res_HC50").c_str());                                                          
  RooUnfoldResponse *rDS = (RooUnfoldResponse*) fres->Get((hname+"_pt_res_DS").c_str());                                                           
  RooUnfoldResponse *rGS = (RooUnfoldResponse*) fres->Get((hname+"_pt_res_GS").c_str());                                                           
  RooUnfoldResponse *rMS = (RooUnfoldResponse*) fres->Get((hname+"_pt_res_MS").c_str());                                                                       

  TH2D* m_pt_dat = (TH2D*) fdat->Get((hname+"_v_pt").c_str());

  RooUnfoldBayes *unfold_nom = new RooUnfoldBayes(rnom, m_pt_dat, 4, false, "unfold_nom","");
  RooUnfoldBayes *unfold_IP2 = new RooUnfoldBayes(rnom, m_pt_dat, 2, false, "unfold_IP2","");
  RooUnfoldBayes *unfold_IP6 = new RooUnfoldBayes(rnom, m_pt_dat, 6, false, "unfold_IP6","");
  RooUnfoldBayes *unfold_TS = new RooUnfoldBayes(rTS, m_pt_dat, 4, false, "unfold_TS","");
  RooUnfoldBayes *unfold_TU = new RooUnfoldBayes(rTU, m_pt_dat, 4, false, "unfold_TU","");
  RooUnfoldBayes *unfold_HC50 = new RooUnfoldBayes(rHC50, m_pt_dat, 4, false, "unfold_HC50","");
  RooUnfoldBayes *unfold_DS = new RooUnfoldBayes(rDS, m_pt_dat, 4, false, "unfold_DS","");
  RooUnfoldBayes *unfold_GS = new RooUnfoldBayes(rGS, m_pt_dat, 4, false, "unfold_GS","");
  RooUnfoldBayes *unfold_MS = new RooUnfoldBayes(rMS, m_pt_dat, 4, false, "unfold_MS","");

  TH2D *reco_nom = (TH2D*) unfold_nom->Hreco((RooUnfold::ErrorTreatment) 3);
  TH2D *reco_IP2 = (TH2D*) unfold_IP2->Hreco((RooUnfold::ErrorTreatment) 3);
  TH2D *reco_IP6 = (TH2D*) unfold_IP6->Hreco((RooUnfold::ErrorTreatment) 3);
  TH2D *reco_TS = (TH2D*) unfold_TS->Hreco((RooUnfold::ErrorTreatment) 3);
  TH2D *reco_TU = (TH2D*) unfold_TU->Hreco((RooUnfold::ErrorTreatment) 3);
  TH2D *reco_HC50 = (TH2D*) unfold_HC50->Hreco((RooUnfold::ErrorTreatment) 3);
  TH2D *reco_DS = (TH2D*) unfold_DS->Hreco((RooUnfold::ErrorTreatment) 3);
  TH2D *reco_GS = (TH2D*) unfold_GS->Hreco((RooUnfold::ErrorTreatment) 3);
  TH2D *reco_MS = (TH2D*) unfold_MS->Hreco((RooUnfold::ErrorTreatment) 3);
  
  //These ranges are okay now that Projection2D has the plotting bug removed.                                                                                   
  const int nBins = 2/*3*/;
  TAxis* reco_axis = reco_nom->GetYaxis(); TAxis* det_axis = m_pt_dat->GetYaxis();
  double ranges[nBins + 1] = {(double) reco_axis->FindBin(20), /*(double) reco_axis->FindBin(25),*/ (double) reco_axis->FindBin(30), (double) reco_axis->FindBin(45)/*(40)*/};                                                                                                                                           
  double ranges_d[nBins + 1] = {(double) det_axis->FindBin(20),/* (double) det_axis->FindBin(25),*/ (double) det_axis->FindBin(30), (double) det_axis->FindBin(45)/*(40)*/};                                                                                                                                              
  string pts[nBins + 1] = {"20",/*"25",*/"30","45"/*"40"*/};

  vector<TH1D*> reco_noms = Projection2D(reco_nom,nBins,ranges,"x");
  vector<TH1D*> reco_IP2s = Projection2D(reco_IP2,nBins,ranges,"x");
  vector<TH1D*> reco_IP6s = Projection2D(reco_IP6,nBins,ranges,"x");
  vector<TH1D*> reco_TSs = Projection2D(reco_TS,nBins,ranges,"x");
  vector<TH1D*> reco_TUs = Projection2D(reco_TU,nBins,ranges,"x");
  vector<TH1D*> reco_HC50s = Projection2D(reco_HC50,nBins,ranges,"x");
  vector<TH1D*> reco_DSs = Projection2D(reco_DS,nBins,ranges,"x");
  vector<TH1D*> reco_GSs = Projection2D(reco_GS,nBins,ranges,"x");
  vector<TH1D*> reco_MSs = Projection2D(reco_MS,nBins,ranges,"x");

  vector<TH1D*> reco_noms_copy;
  
  for (int i = 0; i < nBins; ++ i) {
    reco_noms_copy.push_back((TH1D*) reco_noms[i]->Clone(("nom_"+to_string(i)).c_str()));
    
    reco_noms[i]->Scale(1/(double)reco_noms[i]->Integral());
    reco_IP2s[i]->Scale(1/(double)reco_IP2s[i]->Integral());
    reco_IP6s[i]->Scale(1/(double)reco_IP6s[i]->Integral());
    reco_TSs[i]->Scale(1/(double)reco_TSs[i]->Integral());
    reco_TUs[i]->Scale(1/(double)reco_TUs[i]->Integral());
    reco_HC50s[i]->Scale(1/(double)reco_HC50s[i]->Integral());
    reco_DSs[i]->Scale(1/(double)reco_DSs[i]->Integral());
    reco_GSs[i]->Scale(1/(double)reco_GSs[i]->Integral());
    reco_MSs[i]->Scale(1/(double)reco_MSs[i]->Integral());

    reco_IP2s[i]->Divide(reco_noms[i]);
    reco_IP6s[i]->Divide(reco_noms[i]);
    reco_TSs[i]->Divide(reco_noms[i]);
    reco_TUs[i]->Divide(reco_noms[i]);
    reco_HC50s[i]->Divide(reco_noms[i]);
    reco_DSs[i]->Divide(reco_noms[i]);
    reco_GSs[i]->Divide(reco_noms[i]);
    reco_MSs[i]->Divide(reco_noms[i]);
    
    Double_t stats[5] = {0,0,0,0,0};
    reco_IP2s[i]->PutStats(stats);
    reco_IP6s[i]->PutStats(stats);
    reco_TSs[i]->PutStats(stats);
    reco_TUs[i]->PutStats(stats);
    reco_HC50s[i]->PutStats(stats);
    reco_DSs[i]->PutStats(stats);
    reco_GSs[i]->PutStats(stats);
    reco_MSs[i]->PutStats(stats);

    reco_IP2s[i]->Sumw2(0);
    reco_IP6s[i]->Sumw2(0);
    reco_TSs[i]->Sumw2(0);
    reco_TUs[i]->Sumw2(0);
    reco_HC50s[i]->Sumw2(0);
    reco_DSs[i]->Sumw2(0);
    reco_GSs[i]->Sumw2(0);
    reco_MSs[i]->Sumw2(0);
    
    for (int j = 1; j <= 14; ++ j) {
      //turning the ratio into a percentage.                                                                                                                   
      reco_IP2s[i]->SetBinContent(j,fabs(reco_IP2s[i]->GetBinContent(j) - 1));
      reco_IP6s[i]->SetBinContent(j,fabs(reco_IP6s[i]->GetBinContent(j) - 1));
      reco_TSs[i]->SetBinContent(j,fabs(reco_TSs[i]->GetBinContent(j) - 1));
      reco_TUs[i]->SetBinContent(j,fabs(reco_TUs[i]->GetBinContent(j) - 1));
      reco_HC50s[i]->SetBinContent(j,fabs(reco_HC50s[i]->GetBinContent(j) - 1));
      reco_DSs[i]->SetBinContent(j,fabs(reco_DSs[i]->GetBinContent(j) - 1));
      reco_GSs[i]->SetBinContent(j,fabs(reco_GSs[i]->GetBinContent(j) - 1));
      reco_MSs[i]->SetBinContent(j,fabs(reco_MSs[i]->GetBinContent(j) - 1));
    }
    
    Prettify1DwLineStyle(reco_IP2s[i],2, kSolid, 2,(htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1);
    reco_IP2s[i]->SetFillColor(2); reco_IP2s[i]->SetFillStyle(3305);
    Prettify1DwLineStyle(reco_IP6s[i],3, kSolid, 2,(htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1);
    reco_IP6s[i]->SetFillColor(3); reco_IP6s[i]->SetFillStyle(3395);
    Prettify1DwLineStyle(reco_TSs[i],kGreen+2, kSolid, 2,(htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1);
    reco_TSs[i]->SetFillColor(kGreen+2); reco_TSs[i]->SetFillStyle(3490);
    Prettify1DwLineStyle(reco_TUs[i],kMagenta+1, kSolid, 2,(htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1);
    reco_TUs[i]->SetFillColor(kMagenta+1); reco_TUs[i]->SetFillStyle(3436);
    Prettify1DwLineStyle(reco_HC50s[i],6, kSolid, 2,(htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1);
    reco_HC50s[i]->SetFillColor(6); reco_HC50s[i]->SetFillStyle(3335);
    Prettify1DwLineStyle(reco_DSs[i],8, kSolid, 2,(htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1);
    reco_DSs[i]->SetFillColor(8); reco_DSs[i]->SetFillStyle(3944);
    Prettify1DwLineStyle(reco_GSs[i],9, kSolid, 2,(htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1);
    reco_GSs[i]->SetFillColor(9); reco_GSs[i]->SetFillStyle(3544);
    Prettify1DwLineStyle(reco_MSs[i],11, kSolid, 2,(htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1);
    reco_MSs[i]->SetFillColor(11); reco_MSs[i]->SetFillStyle(3690);
  }
  
  TLegend *l1 = new TLegend(0.5,0.3,0.7,0.6); l1->SetBorderSize(0);
  l1->AddEntry(reco_IP2s[0],"IP2","f");
  l1->AddEntry(reco_IP6s[0],"IP6","f");
  l1->AddEntry(reco_TSs[0],"TS","f");
  l1->AddEntry(reco_TUs[0],"TU","f");

  TLegend *l2 = new TLegend(0.5,0.3,0.7,0.6); l2->SetBorderSize(0);
  l2->AddEntry(reco_HC50s[0],"HC50","f");
  l2->AddEntry(reco_DSs[0],"DS","f");
  l2->AddEntry(reco_GSs[0],"GS","f");
  l2->AddEntry(reco_MSs[0],"MS","f");


  TLatex *slice = new TLatex();

  TCanvas *csys = MakeCanvas ("csys","0",800,500);
  DivideCanvas(csys,"0", nBins,1);

  TH1D* ldummy = new TH1D("ldummy",(";"+htitle+" [GeV/c^{2}];relative uncertainty").c_str(),10,0,9.999);
  TH1D* rdummy = new TH1D("rdummy",(";"+htitle+" [GeV/c^{2}];relative uncertainty").c_str(),10,0.001,9.999);

  ldummy->GetYaxis()->SetRangeUser(0,0.999);
  rdummy->GetYaxis()->SetRangeUser(0,0.999);

  TLatex *p = new TLatex();
  //p->SetTextAlign(11);
  //p->SetTextSize(0.07);

  TLatex *t = new TLatex();
  for (int i = 0; i < nBins; ++ i) {
    csys->cd(i+1);
    if (i == 0) {ldummy->Draw();} if (i == 1) {rdummy->Draw();}
    reco_HC50s[i]->Draw("LF2same"); 
    reco_TSs[i]->Draw("LF2same");
    reco_TUs[i]->Draw("LF2same");
    reco_DSs[i]->Draw("LF2same");
    reco_GSs[i]->Draw("LF2same");
    reco_IP2s[i]->Draw("LF2same");
    reco_IP6s[i]->Draw("LF2same");
    reco_MSs[i]->Draw("LF2same");
    if (i == 0) {
      p->DrawLatexNDC(0.3,0.9, "pp 200 GeV run12 JP2");
      p->DrawLatexNDC(0.3,0.85, "anti-k_{T}, R = 0.4");
      p->DrawLatexNDC(0.3,0.8, "Ch+Ne jets, |#eta| < 0.4");
    }
    if (i==0) {l1->Draw("same");} if (i==1) {l2->Draw("same");}
    slice->DrawLatexNDC(0.3,0.7,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
  }

  //csys->SaveAs(("~/jetmass2/plots/DNP_talk/"+fstart+"systematics.pdf").c_str());
  
  
  //taking maximum envelopes!                                                                                                                                  
  TH2D* env_HC = new TH2D("env_HC","",14,0,14,15,5,80);
  vector<TH1D*> env_HCs = Projection2D(env_HC,nBins,ranges,"x");
  TH2D* env_un = new TH2D("env_un","",14,0,14,15,5,80);
  vector<TH1D*> env_uns = Projection2D(env_un,nBins,ranges,"x");
  TH2D* net = new TH2D("net","",14,0,14,15,5,80);
  vector<TH1D*> nets = Projection2D(net,nBins,ranges,"x");
  TH2D* stat = new TH2D("stat","",14,0,14,15,5,80);
  vector<TH1D*> stats = Projection2D(stat,nBins,ranges,"x");
  
  const int nBins_m = 14;

  vector<vector< double> > syst_errs2D;

  for (int i = 0; i < nBins; ++ i) {
    vector<double> syst_errs1D;
    reco_noms_copy[i]->Scale(1/(double)reco_noms_copy[i]->Integral());

    for (int j = 1; j < nBins_m + 1; ++ j) {
      //hadronic correction envelope - using an ordered set here to automatically get the largest value.                                                        
      double hcs [1] = {reco_HC50s[i]->GetBinContent(j)};
      set<double> hc_sort (hcs, hcs+1);
      set<double>::iterator hc = hc_sort.end(); hc --;
      double hc_envelope = *hc;
      env_HCs[i]->SetBinContent(j, hc_envelope);
      //unfolding envelope                                                                                                                                     
      double uns [5] = {reco_DSs[i]->GetBinContent(j), reco_GSs[i]->GetBinContent(j), reco_MSs[i]->GetBinContent(j), reco_IP2s[i]->GetBinContent(j), reco_IP6s[i]->GetBinContent(j)};
      set<double> un_sort (uns, uns+5);
      set<double>::iterator un = un_sort.end(); un --;
      double un_envelope = *un;
      env_uns[i]->SetBinContent(j, un_envelope);
      //total uncertainty = TU + TS + un envelope + hc envelope   
      double square = pow(hc_envelope,2) + pow(un_envelope,2) + pow(reco_TUs[i]->GetBinContent(j),2) + pow(reco_TSs[i]->GetBinContent(j),2);
      nets[i]->SetBinContent(j,sqrt(square));
      stats[i]->SetBinContent(j,reco_noms_copy[i]->GetBinError(j));
      syst_errs1D.push_back(nets[i]->GetBinContent(j));
    }
    
    syst_errs2D.push_back(syst_errs1D);
  }
  
  for (int i = 0; i < nBins; ++ i) {
    Prettify1DwLineStyle(env_HCs[i],kRed, kSolid, 2,(htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1);
    env_HCs[i]->SetFillColor(kRed); env_HCs[i]->SetFillStyle(3353);
    Prettify1DwLineStyle(env_uns[i],kBlue, kSolid, 2,(htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1);
    env_uns[i]->SetFillColor(kBlue); env_uns[i]->SetFillStyle(3305);
    Prettify1DwLineStyle(nets[i], kBlack, kSolid, 2, (htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1);
    Prettify1DwLineStyle(stats[i], kBlack, kDashed, 2, (htitle+" [GeV/c^{2}]").c_str(),"relative uncertainty",0,14,0,1);
  }
  
  TLegend *tenvs = new TLegend(0.15,0.4,0.45,0.75); tenvs->SetBorderSize(0);
  tenvs->AddEntry(env_HCs[0],"Hadronic correction","f");
  tenvs->AddEntry(reco_TSs[0],"Tower scale","f");
  tenvs->AddEntry(reco_TUs[0],"Tracking","f");
  tenvs->AddEntry(env_uns[0],"Unfolding","f");
  tenvs->AddEntry(nets[0],"Total systematic","l");

  TCanvas *cenvs = new TCanvas("cenvs","cenvs",800,500);
  DivideCanvas(cenvs,"0",nBins,1);

  for (int i = 0; i < nBins; ++ i) {
    cenvs->cd(i+1);
    if (i == 0) {ldummy->Draw();} if (i == 1) {rdummy->Draw();}
    env_uns[i]->Draw("lf2same PLC PFC PMC");
    reco_TUs[i]->Draw("lf2same PLC PFC PMC");
    env_HCs[i]->Draw("lf2same PLC PFC PMC");
    reco_TSs[i]->Draw("lf2same PLC PFC PMC");
    nets[i]->Draw("lf2same");
    if (i == 0) {
      slice->DrawLatexNDC(0.3,0.7,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
      p->DrawLatexNDC(0.3,0.9, "pp 200 GeV 2012");
      p->DrawLatexNDC(0.3,0.85, "anti-k_{T}, R = 0.4");
      p->DrawLatexNDC(0.3,0.8, "Ch+Ne jets, |#eta| < 0.6");
    }
    if (i == 1) {
      slice->DrawLatexNDC(0.3,0.8,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
      tenvs->Draw("same");
    }
  }
  /*
  for (int i = 1; i <= nets[1]->GetNbinsX(); ++ i) {
    if (i == 2 || i == 5 || i == 8) {
      cout << 100*env_HCs[1]->GetBinContent(i) << "% " << 100*reco_TSs[1]->GetBinContent(i) << "% " << 100*reco_TUs[1]->GetBinContent(i) << "% " << 100*env_uns[1]->GetBinContent(i) << "% " << 100*nets[1]->GetBinContent(i) << endl;
    }
  }
  */
  
  cenvs->SaveAs(("~/jetmass2/plots/DNP_talk/"+fstart+"systematic_envelopes.pdf").c_str());
  
}

//when running functions with arguments in an interactive root session, do e.g. root "AN_plots.C(4, 0)" - the quotation marks are necessary 

void AN_plots(int radius, bool groom) {
  gROOT->ForceStyle(); //forces use of ~/rootlogon.C's style settings
  
  //~~~SELECT A JET RADIUS AND GROOMED/UNGROOMED JETS~~~//
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
  //~~~~~~//

  //~~~OPEN NECESSARY FILES~~~//
  TFile *fdat = new TFile(("~/jetmass2/out/data/hists/data_hists_ppJP2"+radstring+".root").c_str(),"READ");
  TFile *fp6unmatch = new TFile(("~/jetmass2/out/sim/hists/unmatched_hists"+radstring+".root").c_str(),"READ");
  TFile *fp6match = new TFile(("~/jetmass2/out/sim/hists/matched_hists"+radstring+".root").c_str(),"READ");
  TFile *fp8 = new TFile(("~/jetmass2/production/out/pythia/hists/pythia8"+radstring+"_undecayed_hists.root").c_str(),"READ");
  TFile *fh7 = new TFile(("~/jetmass/production/macros/hists/hists_allsim_lowzgremoved"+radstring+".root").c_str(),"READ");//temporarily using the old files (they're the same, but will point to new ones later)                                                                                                             
  //  TFile *fPL = new TFile(("~/jetmass2/production/out/pythia/hists/pythia8"+radstring+"_FSR_only_undecayed_hists.root").c_str(),"READ");
  
  TFile *fclos = new TFile(("~/jetmass2/out/closure/closure"+radstring+".root").c_str(),"READ");
  TFile *fsysts = new TFile(("~/jetmass2/out/sim/sim_matched"+radstring+"_bindropped.root").c_str(),"READ");

  //plots the uncorrected jet mass for a bin of pT and compares to MC
  //  rawdata(fdat, fp6unmatch);
  
  //plots the mass resolution for different pT ranges on the same panel
  //  resolution(fp6match);

  //plots the closure test results for 1D and 2D unfolding
  //  closure(fclos);
  
  //plots the systematic uncertainties
    systematics(fsysts, fdat, hname, htitle, fstart);
  
  //~~~~~~//
  


  return;
}
