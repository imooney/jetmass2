#include "RooUnfoldResponse.h"
#include "Plots_old.h"
#include <algorithm>

using namespace std;

void temp_systematics () {
  gStyle->SetPalette(kPastel);

  TFile *fold = new TFile("~/jetmass2_11-10-2020_11_10-2020/out/unfold/groomed_unfolded_v1.root","READ");
  vector<TH1D*> systs_old = {(TH1D*) fold->Get("w_systs_0"),(TH1D*) fold->Get("w_systs_1"),(TH1D*) fold->Get("w_systs_2")};
  vector<TH1D*> reco_old = {(TH1D*) fold->Get("nom_0"),(TH1D*) fold->Get("nom_1"),(TH1D*) fold->Get("nom_2")};
  
  TFile *fres = new TFile("~/jetmass2_11-10-2020_11_10-2020/out/sim/sim_matched_allbugsfixed_bindropped.root","READ");
  RooUnfoldResponse *rnom = (RooUnfoldResponse*) fres->Get("mg_pt_res_nom"); //!
  RooUnfoldResponse *rTS = (RooUnfoldResponse*) fres->Get("mg_pt_res_TS"); //!
  RooUnfoldResponse *rTU = (RooUnfoldResponse*) fres->Get("mg_pt_res_TU"); //!
  RooUnfoldResponse *rHC50 = (RooUnfoldResponse*) fres->Get("mg_pt_res_HC50"); //!
  RooUnfoldResponse *rDS = (RooUnfoldResponse*) fres->Get("mg_pt_res_DS"); //!
  RooUnfoldResponse *rGS = (RooUnfoldResponse*) fres->Get("mg_pt_res_GS"); //!
  RooUnfoldResponse *rMS = (RooUnfoldResponse*) fres->Get("mg_pt_res_MS"); //!
  
  TFile *fdat = new TFile("~/jetmass2_11-10-2020_11_10-2020/out/data/hists/data_hists_ppJP2_bindropped.root","READ");
  
  TH2D* mg_pt_dat = (TH2D*) fdat->Get("mg_v_pt");
  
  RooUnfoldBayes *unfold_nom = new RooUnfoldBayes(rnom, mg_pt_dat, 4, false, "unfold_nom","");
  RooUnfoldBayes *unfold_IP2 = new RooUnfoldBayes(rnom, mg_pt_dat, 2, false, "unfold_IP2","");
  RooUnfoldBayes *unfold_IP6 = new RooUnfoldBayes(rnom, mg_pt_dat, 6, false, "unfold_IP6","");
  RooUnfoldBayes *unfold_TS = new RooUnfoldBayes(rTS, mg_pt_dat, 4, false, "unfold_TS","");
  RooUnfoldBayes *unfold_TU = new RooUnfoldBayes(rTU, mg_pt_dat, 4, false, "unfold_TU","");
  RooUnfoldBayes *unfold_HC50 = new RooUnfoldBayes(rHC50, mg_pt_dat, 4, false, "unfold_HC50","");
  RooUnfoldBayes *unfold_DS = new RooUnfoldBayes(rDS, mg_pt_dat, 4, false, "unfold_DS","");
  RooUnfoldBayes *unfold_GS = new RooUnfoldBayes(rGS, mg_pt_dat, 4, false, "unfold_GS","");
  RooUnfoldBayes *unfold_MS = new RooUnfoldBayes(rMS, mg_pt_dat, 4, false, "unfold_MS","");
  
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
  const int nBins = 3;
  TAxis* reco_axis = reco_nom->GetYaxis(); TAxis* det_axis = mg_pt_dat->GetYaxis();  
  double ranges[nBins + 1] = {(double) reco_axis->FindBin(20), (double) reco_axis->FindBin(25), (double) reco_axis->FindBin(30), (double) reco_axis->FindBin(40)};//3,4,5,6,8,12};
  double ranges_d[nBins + 1] = {(double) det_axis->FindBin(20), (double) det_axis->FindBin(25), (double) det_axis->FindBin(30), (double) det_axis->FindBin(40)};//1,2,3,4,6,10};
  string pts[nBins + 1] = {"20","25","30","40"};
  
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
  vector<TH1D*> reco_IP2s_copy;
  vector<TH1D*> reco_IP6s_copy;
  vector<TH1D*> reco_TSs_copy;
  vector<TH1D*> reco_TUs_copy;
  vector<TH1D*> reco_HC50s_copy;
  vector<TH1D*> reco_DSs_copy;
  vector<TH1D*> reco_GSs_copy;
  vector<TH1D*> reco_MSs_copy;

  
  for (int i = 0; i < nBins; ++ i) {
    reco_noms_copy.push_back((TH1D*) reco_noms[i]->Clone(("nom_"+to_string(i)).c_str()));
    reco_IP2s_copy.push_back((TH1D*) reco_IP2s[i]->Clone(("IP2_"+to_string(i)).c_str()));
    reco_IP6s_copy.push_back((TH1D*) reco_IP6s[i]->Clone(("IP6_"+to_string(i)).c_str()));
    reco_TSs_copy.push_back((TH1D*) reco_TSs[i]->Clone(("TS_"+to_string(i)).c_str()));
    reco_TUs_copy.push_back((TH1D*) reco_TUs[i]->Clone(("TU_"+to_string(i)).c_str()));
    reco_HC50s_copy.push_back((TH1D*) reco_HC50s[i]->Clone(("HC50_"+to_string(i)).c_str()));
    reco_DSs_copy.push_back((TH1D*) reco_DSs[i]->Clone(("DS_"+to_string(i)).c_str()));
    reco_GSs_copy.push_back((TH1D*) reco_GSs[i]->Clone(("GS_"+to_string(i)).c_str()));
    reco_MSs_copy.push_back((TH1D*) reco_MSs[i]->Clone(("MS_"+to_string(i)).c_str()));

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
    reco_MSs[i]->Scale(1/(double)reco_MSs[i]->Integral());
    
    reco_IP2s[i]->Divide(reco_noms[i]);
    reco_IP6s[i]->Divide(reco_noms[i]);
    reco_TSs[i]->Divide(reco_noms[i]);
    reco_TUs[i]->Divide(reco_noms[i]);
    reco_HC50s[i]->Divide(reco_noms[i]);
    reco_DSs[i]->Divide(reco_noms[i]);
    reco_GSs[i]->Divide(reco_noms[i]);
    reco_MSs[i]->Divide(reco_noms[i]);
  }
  
  for (int i = 0; i < nBins; ++ i) {
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
    
  }
  
  for (int i = 0; i < nBins; ++ i) {
    Prettify1DwLineStyle(reco_IP2s[i],2, kSolid, 2,"M_{g} [GeV/c^{2}]","relative uncertainty",0,14,0,1); 
    reco_IP2s[i]->SetFillColor(2); reco_IP2s[i]->SetFillStyle(3305);
    Prettify1DwLineStyle(reco_IP6s[i],3, kSolid, 2,"M_{g} [GeV/c^{2}]","relative uncertainty",0,14,0,1);
    reco_IP6s[i]->SetFillColor(3); reco_IP6s[i]->SetFillStyle(3395);
    Prettify1DwLineStyle(reco_TSs[i],4, kSolid, 2,"M_{g} [GeV/c^{2}]","relative uncertainty",0,14,0,1);
    reco_TSs[i]->SetFillColor(4); reco_TSs[i]->SetFillStyle(3490);
    Prettify1DwLineStyle(reco_TUs[i],5, kSolid, 2,"M_{g} [GeV/c^{2}]","relative uncertainty",0,14,0,1);
    reco_TUs[i]->SetFillColor(5); reco_TUs[i]->SetFillStyle(3436);
    Prettify1DwLineStyle(reco_HC50s[i],6, kSolid, 2,"M_{g} [GeV/c^{2}]","relative uncertainty",0,14,0,1);
    reco_HC50s[i]->SetFillColor(6); reco_HC50s[i]->SetFillStyle(3335);
    Prettify1DwLineStyle(reco_DSs[i],8, kSolid, 2,"M_{g} [GeV/c^{2}]","relative uncertainty",0,14,0,1); 
    reco_DSs[i]->SetFillColor(8); reco_DSs[i]->SetFillStyle(3944);
    Prettify1DwLineStyle(reco_GSs[i],9, kSolid, 2,"M_{g} [GeV/c^{2}]","relative uncertainty",0,14,0,1);
    reco_GSs[i]->SetFillColor(9); reco_GSs[i]->SetFillStyle(3544);
    Prettify1DwLineStyle(reco_MSs[i],11, kSolid, 2,"M_{g} [GeV/c^{2}]","relative uncertainty",0,14,0,1);
    reco_MSs[i]->SetFillColor(11); reco_MSs[i]->SetFillStyle(3690);

  }
  
  TLegend *l1 = new TLegend(0.1,0.5,0.28,0.8); l1->SetBorderSize(0);
  l1->AddEntry(reco_IP2s[0],"IP2","f");
  l1->AddEntry(reco_IP6s[0],"IP6","f");
  l1->AddEntry(reco_TSs[0],"TS","f");
  l1->AddEntry(reco_TUs[0],"TU","f");

  TLegend *l2 = new TLegend(0.1,0.5,0.28,0.8); l2->SetBorderSize(0);
  l2->AddEntry(reco_HC50s[0],"HC50","f");
  l2->AddEntry(reco_DSs[0],"DS","f");
  l2->AddEntry(reco_GSs[0],"GS","f");
  l2->AddEntry(reco_MSs[0],"MS","f");
 

  TLatex *slice = new TLatex();
  
  TCanvas *csys = MakeCanvas ("csys","0",1200,500);
  DivideCanvas(csys,"0", 3,1);
  
  TH1D* hdummy = new TH1D("hdummy",";M_{g} [GeV/c^{2}];relative uncertainty",14,0,14);
  hdummy->GetYaxis()->SetRangeUser(0,1);
    
  TLatex *p = new TLatex();
  p->SetTextAlign(11);
  p->SetTextSize(0.07);
  
  /*csys->cd(1); hdummy->Draw();*/ TLatex *t = new TLatex();
  for (int i = 0; i < nBins; ++ i) {
    csys->cd(i+1); 
    reco_HC50s[i]->Draw("lf2same   "); reco_TSs[i]->Draw("lf2same   "); reco_TUs[i]->Draw("lf2same   ");
    reco_DSs[i]->Draw("lf2same   "); reco_GSs[i]->Draw("lf2same   "); reco_IP2s[i]->Draw("lf2same   "); reco_IP6s[i]->Draw("lf2same   "); reco_MSs[i]->Draw("lf2same   ");
    if (i == 0) {
      p->DrawLatexNDC(0.2,0.65, "pp 200 GeV run12 JP2");
      p->DrawLatexNDC(0.2,0.55, "anti-k_{T}, R = 0.4");
      p->DrawLatexNDC(0.2,0.45, "Ch+Ne jets, |#eta| < 0.4");
    }
    if (i==1) {l1->Draw("same");} if (i==2) {l2->Draw("same");} slice->DrawLatexNDC(0.3,0.8,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
  }
    
  //  csys->SaveAs("~/jetmass/plots/systematics/systematics_R04.pdf");


  

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
    Prettify1DwLineStyle(env_HCs[i],7, kSolid, 2,"M_{g} [GeV/c^{2}]","relative uncertainty",0,14,0,1);
    env_HCs[i]->SetFillColor(7); env_HCs[i]->SetFillStyle(3353);
    Prettify1DwLineStyle(env_uns[i],2, kSolid, 2,"M_{g} [GeV/c^{2}]","relative uncertainty",0,14,0,1);
    env_uns[i]->SetFillColor(2); env_uns[i]->SetFillStyle(3305);
    Prettify1DwLineStyle(nets[i], kBlack, kSolid, 2, "M_{g} [GeV/c^{2}]","relative uncertainty",0,14,0,1);
    Prettify1DwLineStyle(stats[i], kBlack, kDashed, 2, "M_{g} [GeV/c^{2}]","relative uncertainty",0,14,0,1);
  }
  
  TLegend *tenvs = new TLegend(0.1,0.4,0.4,0.75); tenvs->SetBorderSize(0);
  tenvs->AddEntry(env_HCs[0],"Hadronic correction","f");
  tenvs->AddEntry(reco_TSs[0],"Tower scale","f");
  tenvs->AddEntry(reco_TUs[0],"Tracking","f");
  tenvs->AddEntry(env_uns[0],"Unfolding","f");
  tenvs->AddEntry(nets[0],"Total systematic uncertainty","l");
  
  TCanvas *cenvs = new TCanvas("cenvs","cenvs",1200,500);
  DivideCanvas(cenvs,"0",3,1);
  
  TH1D* hdummyenvs = new TH1D("hdummyenvs",";;relative uncertainty",1,0,14);
  hdummyenvs->GetYaxis()->SetRangeUser(0,1);

  for (int i = 0; i < nBins; ++ i) {
    cenvs->cd(i+1); 
    env_HCs[i]->Draw("lf2same   "); reco_TSs[i]->Draw("lf2same   "); reco_TUs[i]->Draw("lf2same   "); env_uns[i]->Draw("lf2same   "); nets[i]->Draw("lf2same"); slice->DrawLatexNDC(0.3,0.8,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
    if (i == 0) {
      p->DrawLatexNDC(0.2,0.65, "pp 200 GeV run12 JP2");
      p->DrawLatexNDC(0.2,0.55, "anti-k_{T}, R = 0.4");
      p->DrawLatexNDC(0.2,0.45, "Ch+Ne jets, |#eta| < 0.6");
    }
    if (i == 1) {tenvs->Draw("same");}
  }
  
  for (int i = 1; i <= nets[1]->GetNbinsX(); ++ i) {
    if (i == 2 || i == 5 || i == 8) {
      cout << 100*env_HCs[1]->GetBinContent(i) << "% " << 100*reco_TSs[1]->GetBinContent(i) << "% " << 100*reco_TUs[1]->GetBinContent(i) << "% " << 100*env_uns[1]->GetBinContent(i) << "% " << 100*nets[1]->GetBinContent(i) << endl;
    }
  }
   
  // cenvs->SaveAs("~/jetmass/plots/systematics/systematic_envelopes_R04.pdf");
  
  //unfolded result with systematic errors!
  vector<TH1D*> w_systs;
  for (int i = 0; i < nBins; ++ i) { 
    w_systs.push_back((TH1D*) reco_noms_copy[i]->Clone(("w_systs_"+to_string(i)).c_str()));
    for (int j = 1; j < nBins_m + 1; ++ j) {
      w_systs[i]->SetBinError(j, syst_errs2D[i][j-1]*w_systs[i]->GetBinContent(j));
    }
  }

  vector<TH1D*> dats = Projection2D(mg_pt_dat,nBins,ranges_d,"x");

  for (int i = 0; i < nBins; ++ i) {
    Prettify1D(reco_noms_copy[i],kRed,kFullStar,4,kRed,"M_{g} [GeV/c^{2}]","1/N dN/dM_{g}",0,10,0,0.5);
    Prettify1D(dats[i], kBlack, kOpenStar, 4/*4*/, kBlack, "M_{g} [GeV/c^{2}]", "1/N dN/dM_{g}",0,10,0,0.5);
    Prettify1D(w_systs[i],kRed,kFullStar,0,kRed,"M_{g} [GeV/c^{2}]","1/N dN/dM_{g}",0,10,0,0.5);
    Prettify1D(systs_old[i],kBlue,kFullStar,0,kBlue,"M_{g} [GeV/c^{2}]","1/N dN/dM_{g}",0,10,0,0.5);
    Prettify1D(reco_old[i],kBlue,kFullStar,4,kBlue,"M_{g} [GeV/c^{2}]","1/N dN/dM_{g}",0,10,0,0.5);
    w_systs[i]->SetFillColor(kRed - 10); w_systs[i]->SetFillStyle(/*3352*/1001);
    systs_old[i]->SetFillColor(kBlue - 10); systs_old[i]->SetFillStyle(/*3352*/1001);   
  }
   
  //scaling errors
  TFile *fstats = new TFile("~/jetmass2_11-10-2020_11_10-2020/out/sim/stat_err_scaling.root","READ");
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
      double binerror = reco_noms_copy[i]->GetBinError(j);
      reco_noms_copy[i]->SetBinError(j,(double) binerror*scaling);
      cout << "bin " << j << " error: " << binerror << endl;
      cout << "corresponding raw data error: " << dats[i]->GetBinError(j) << endl;
      cout << "scaling bin " << j << " error by: " << scaling << endl;
      cout << "scaled error: " << binerror*scaling << endl;
    }
  }
  
  TCanvas *cws = new TCanvas("cws","cws",1200,500);
  DivideCanvas(cws,"0",3,1);

  TLegend *twsysts2 = new TLegend(0.5,0.55,0.75,0.7); twsysts2->SetBorderSize(0);
  //twsysts2->AddEntry(dats[0],"Raw data","p");
  TH1D* for_legend = (TH1D*) w_systs[0]->Clone("for_legend"); for_legend->SetMarkerSize(2);
  TH1D* for_legend_old = (TH1D*) systs_old[0]->Clone("for_legend_old"); for_legend_old->SetMarkerSize(2);
  twsysts2->AddEntry(for_legend,"Unfolded data (v2)","pf");
  twsysts2->AddEntry(for_legend_old,"Unfolded data (v1)","pf");
  //twsysts2->AddEntry(reco_old[0],"Unfolded data (v1)","p");
  
  TLatex *tpost = new TLatex(); tpost->SetTextColor(kRed);
  
  TLatex *ttitle = new TLatex(); ttitle->SetTextAlign(11); ttitle->SetTextSize(0.05);
  t->SetTextAlign(11);
  t->SetTextSize(0.07);
  
  TH1D* hdummycws = new TH1D("hdummycws",";;1/N dN/dM_{g}",1,0,10);
  hdummycws->GetYaxis()->SetRangeUser(0,0.4);
  hdummycws->GetYaxis()->SetTitleSize(0.06);

  //cws->cd(1); hdummycws->Draw(); t = PanelTitle();
  for (int i = 0; i < nBins; ++ i) {

    cws->cd(i+1); w_systs[i]->Draw("E3same9"); systs_old[i]->Draw("E3same9"); /*dats[i]->Draw("same");*/reco_old[i]->Draw("same9"); reco_noms_copy[i]->Draw("same9"); slice->DrawLatexNDC(0.5,0.77,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
    //w_systs[i]->GetXaxis()->SetTitleSize(0.08); w_systs[i]->GetYaxis()->SetTitleSize(0.08);
    if (i == 0) {ttitle->DrawLatex(4.2,0.46, "pp 200 GeV run12 JP2");ttitle->DrawLatex(4.2,0.43, "anti-k_{T}, R = 0.4");ttitle->DrawLatex(4.2,0.4, "Ch+Ne jets, |#eta| < 0.4");}
    if (i == 1) { /*tpost->DrawLatexNDC(0.15,0.87,"PREVIEW - Work in Progress");*/}
    if (i == 2) { twsysts2->Draw("same9");}
    
  }

  //  cws->SaveAs("~/jetmass/plots/systematics/unfolded_w_systs_etc_and_PL_R04.pdf");
  
  
  TCanvas *ccheck = new TCanvas("ccheck","ccheck",1200,500);
  DivideCanvas(ccheck,"0",3,1);
  
  vector<TH1D*> ratio_old_new = {(TH1D*) reco_old[0]->Clone("ratio0"), (TH1D*) reco_old[1]->Clone("ratio1"), (TH1D*) reco_old[2]->Clone("ratio2")};
  vector<TH1D*> systs_ratio = {(TH1D*) systs_old[0]->Clone("sratio0"), (TH1D*) systs_old[1]->Clone("sratio1"), (TH1D*) systs_old[2]->Clone("sratio2")};

  TLine *unity = new TLine (0,1,10,1); unity->SetLineStyle(kDashed);
  
  for (int i = 0; i < w_systs.size(); ++ i) {
    
    for (int j = 1; j < w_systs[i]->GetNbinsX(); ++ j) {
      //cout << "DEBUG: " << ratio_old_new[i]->GetBinError(j) << " " << w_systs[i]->GetBinError(j) << " " << ratio_old_new[i]->GetBinError(j)/(double) w_systs[i]->GetBinError(j) << endl;
      systs_ratio[i]->SetBinContent(j, systs_ratio[i]->GetBinError(j) / (double) w_systs[i]->GetBinError(j)); 
    }
    
    ratio_old_new[i]->Divide(reco_noms_copy[i]/*w_systs[i]*/);
    ratio_old_new[i]->GetYaxis()->SetRangeUser(0,2);
    ratio_old_new[i]->GetXaxis()->SetRangeUser(0,10);
    ratio_old_new[i]->SetMarkerStyle(kOpenCircle);
    ratio_old_new[i]->SetMarkerColor(kMagenta);
    ratio_old_new[i]->SetLineColor(kMagenta);
    ratio_old_new[i]->SetMarkerSize(2);
    systs_ratio[i]->SetMarkerStyle(kOpenCircle);
    systs_ratio[i]->SetMarkerColor(kGreen+2);
    systs_ratio[i]->SetLineColor(kGreen+2);
    systs_ratio[i]->SetMarkerSize(2);
    
    ccheck->cd(i+1); ratio_old_new[i]->Draw(); systs_ratio[i]->Draw("same");
    unity->Draw("same");
  }
  /*  
  TFile *fout = new TFile("~/jetmass2_11-10-2020_11_10-2020/out/unfold/unfolded_v2_bugsfixed.root","RECREATE");
  fout->cd();
  for (int i = 0; i < nBins; ++ i) {
    reco_noms_copy[i]->Write();
    dats[i]->Write();
    w_systs[i]->Write();
  }
  fout->Close();
  */

  return;
}
