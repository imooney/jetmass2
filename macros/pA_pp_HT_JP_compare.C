//Isaac Mooney, WSU, October 2019
//This file plots the comparison between pA HT2, pA JP2, pp HT2, and pp JP2 for the mass and the NEF 

#include "Plots_old.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

using namespace std;


void pA_pp_HT_JP_compare () {
  gROOT->ForceStyle();
  
  const int nPts = 5;
  const int nCols = 3;
  const int nRows = 2;
    
  TCanvas *cnef = new TCanvas("cnef","cnef",800,600);  
  DivideCanvas(cnef,"0",nCols,nRows);

  TCanvas *cm = new TCanvas("cm","cm",800,600);  
  DivideCanvas(cm,"0",nCols,nRows);


  TFile *fpAJP = new TFile("~/jetmass2/out/data/hists/data_hists_pAuJP2.root","READ");
  TFile *fpAHT = new TFile("~/jetmass2/out/data/hists/data_hists_pAuHT2.root","READ");
  TFile *fppJP = new TFile("~/jetmass2/out/data/hists/data_hists_ppJP2_R04.root","READ");
  TFile *fppHT = new TFile("~/jetmass2/out/data/hists/data_hists_ppHT.root","READ");

  TH2D* nef_pAJP_1 = (TH2D*) fpAJP->Get("ch_frac_v_pt");
  TH2D* nef_pAJP = (TH2D*) nef_pAJP_1->Clone("nef_pAJP");
  TH2D* nef_pAHT_1 = (TH2D*) fpAHT->Get("ch_frac_v_pt");
  TH2D* nef_pAHT = (TH2D*) nef_pAHT_1->Clone("nef_pAHT");
  TH2D* nef_ppJP_1 = (TH2D*) fppJP->Get("ch_frac_v_pt");
  TH2D* nef_ppJP = (TH2D*) nef_ppJP_1->Clone("nef_ppJP");
  TH2D* nef_ppHT_1 = (TH2D*) fppHT->Get("ch_frac_v_pt");
  TH2D* nef_ppHT = (TH2D*) nef_ppHT_1->Clone("nef_ppHT");

  TH2D* m_pAJP_1 = (TH2D*) fpAJP->Get("m_v_pt");
  TH2D* m_pAJP = (TH2D*) m_pAJP_1->Clone("m_pAJP");
  TH2D* m_pAHT_1 = (TH2D*) fpAHT->Get("m_v_pt");
  TH2D* m_pAHT = (TH2D*) m_pAHT_1->Clone("m_pAHT");
  TH2D* m_ppJP_1 = (TH2D*) fppJP->Get("m_v_pt");
  TH2D* m_ppJP = (TH2D*) m_ppJP_1->Clone("m_ppJP");
  TH2D* m_ppHT_1 = (TH2D*) fppHT->Get("m_v_pt");
  TH2D* m_ppHT = (TH2D*) m_ppHT_1->Clone("m_ppHT");
  
  TH3D* m_pAJP_1_BBC = (TH3D*) fpAJP->Get("bbc_east_sum_v_pt_v_m");
  TH3D* m_pAJP_BBC_lowEA = (TH3D*) m_pAJP_1_BBC->Clone("m_pAJP_BBC_lowEA");
  TH3D* m_pAJP_BBC_highEA = (TH3D*) m_pAJP_1_BBC->Clone("m_pAJP_BBC_highEA");
  m_pAJP_BBC_lowEA->GetXaxis()->SetRangeUser(8000,16000);
  TH2D* m_pAJP_lowEA = (TH2D*) m_pAJP_BBC_lowEA->Project3D("yz");
  m_pAJP_BBC_highEA->GetXaxis()->SetRangeUser(38000,100000);
  TH2D* m_pAJP_highEA = (TH2D*) m_pAJP_BBC_highEA->Project3D("yz");
  
  TAxis* pt_axis = nef_pAJP->GetYaxis();
  double ranges[nPts + 1] = {(double) pt_axis->FindBin(15), (double) pt_axis->FindBin(20), (double) pt_axis->FindBin(25), (double) pt_axis->FindBin(30), (double) pt_axis->FindBin(40), (double) pt_axis->FindBin(60)};                                                            
  string pts[nPts + 1] = {"15","20","25","30","40","60"};

  vector<TH1D*> nef_pAJPs = Projection2D(nef_pAJP,nPts,ranges,"x");
  vector<TH1D*> nef_pAHTs = Projection2D(nef_pAHT,nPts,ranges,"x");
  vector<TH1D*> nef_ppJPs = Projection2D(nef_ppJP,nPts,ranges,"x");
  vector<TH1D*> nef_ppHTs = Projection2D(nef_ppHT,nPts,ranges,"x");
  
  vector<TH1D*> m_pAJPs = Projection2D(m_pAJP,nPts,ranges,"x");
  vector<TH1D*> m_pAHTs = Projection2D(m_pAHT,nPts,ranges,"x");
  vector<TH1D*> m_ppJPs = Projection2D(m_ppJP,nPts,ranges,"x");
  vector<TH1D*> m_ppHTs = Projection2D(m_ppHT,nPts,ranges,"x");

  vector<TH1D*> m_pAJPs_lowEA = Projection2D(m_pAJP_lowEA,nPts,ranges,"x");
  vector<TH1D*> m_pAJPs_highEA = Projection2D(m_pAJP_highEA,nPts,ranges,"x");
  
  cout << "# of jets total: " << nef_pAJP->GetEntries() << " <--- pA JP" << endl;
  cout << "# of jets total: " << nef_pAHT->GetEntries() << " <--- pA HT" << endl;
  cout << "# of jets total: " << nef_ppJP->GetEntries() << " <--- pp JP" << endl;
  cout << "# of jets total: " << nef_ppHT->GetEntries() << " <--- pp HT" << endl;
  
  cout << "# of jets total: " << m_pAJP_lowEA->GetEntries() << " <---pA JP - low EA" << endl;
  cout << "# of jets total: " << m_pAJP_highEA->GetEntries() << " <---pA JP - high EA" << endl;
  cout << endl;
  
  for (int i = 0; i < nPts; ++ i) {
    cout << "# of jets in bin " << i << ": " << nef_pAJPs[i]->Integral() << " <--- pA JP" << endl;
    cout << "# of jets in bin " << i << ": " << nef_pAHTs[i]->Integral() << " <--- pA HT" << endl;
    cout << "# of jets in bin " << i << ": " << nef_ppJPs[i]->Integral() << " <--- pp JP" << endl;
    cout << "# of jets in bin " << i << ": " << nef_ppHTs[i]->Integral() << " <--- pp HT" << endl;

    cout << "# of jets in bin " << i << ": " << m_pAJPs_lowEA[i]->Integral() << " <--- pA JP - low EA" << endl;
    cout << "# of jets in bin " << i << ": " << m_pAJPs_highEA[i]->Integral() << " <--- pA JP - high EA" << endl;

    cout << endl;
    
    Prettify1D(nef_pAJPs[i],kViolet,kFullStar,2,kViolet,"F = E_{T,ch}/E_{T,tot}","1/N dN/dF",0,1,0,2.5);
    Prettify1D(nef_pAHTs[i],kViolet,kOpenStar,2,kViolet,"F = E_{T,ch}/E_{T,tot}","1/N dN/dF",0,1,0,2.5);
    Prettify1D(nef_ppJPs[i],kOrange,kFullStar,2,kOrange,"F = E_{T,ch}/E_{T,tot}","1/N dN/dF",0,1,0,2.5);
    Prettify1D(nef_ppHTs[i],kOrange,kOpenStar,2,kOrange,"F = E_{T,ch}/E_{T,tot}","1/N dN/dF",0,1,0,2.5);

    Prettify1D(m_pAJPs[i],kViolet,kFullStar,2,kViolet,"M [GeV/c^{2}]","1/N dN/dM",0,10,0,0.55);
    Prettify1D(m_pAHTs[i],kViolet,kOpenStar,2,kViolet,"M [GeV/c^{2}]","1/N dN/dM",0,10,0,0.55);
    Prettify1D(m_ppJPs[i],kOrange,kFullStar,2,kOrange,"M [GeV/c^{2}]","1/N dN/dM",0,10,0,0.55);
    Prettify1D(m_ppHTs[i],kOrange,kOpenStar,2,kOrange,"M [GeV/c^{2}]","1/N dN/dM",0,10,0,0.55);
    
    Prettify1D(m_pAJPs_lowEA[i],kViolet,kFullStar,2,kViolet,"M [GeV/c^{2}]","1/N dN/dM",0,10,0,0.55);
    Prettify1D(m_pAJPs_highEA[i],kViolet,kOpenStar,2,kViolet,"M [GeV/c^{2}]","1/N dN/dM",0,10,0,0.55);
 
    nef_pAJPs[i]->Sumw2();
    nef_pAHTs[i]->Sumw2();
    nef_ppJPs[i]->Sumw2();
    nef_ppHTs[i]->Sumw2();

    m_pAJPs[i]->Sumw2();
    m_pAHTs[i]->Sumw2();
    m_ppJPs[i]->Sumw2();
    m_ppHTs[i]->Sumw2();
    
    m_pAJPs_lowEA[i]->Sumw2();
    m_pAJPs_highEA[i]->Sumw2();
  }
  
  //statistics:
  for (int i = 0; i < nPts; ++ i) {
    cout << endl;
    cout << "mean of bin " << i << ": " << nef_pAJPs[i]->GetMean() << " <--- pA JP NEF" << endl;
    cout << "mean of bin " << i << ": " << nef_pAHTs[i]->GetMean() << " <--- pA HT NEF" << endl;
    cout << "mean of bin " << i << ": " << nef_ppJPs[i]->GetMean() << " <--- pp JP NEF" << endl;
    cout << "mean of bin " << i << ": " << nef_ppHTs[i]->GetMean() << " <--- pp HT NEF" << endl;
    cout << endl;
    cout << "RMS of bin " << i << ": " << nef_pAJPs[i]->GetRMS() << " <--- pA JP NEF" << endl;
    cout << "RMS of bin " << i << ": " << nef_pAHTs[i]->GetRMS() << " <--- pA HT NEF" << endl;
    cout << "RMS of bin " << i << ": " << nef_ppJPs[i]->GetRMS() << " <--- pp JP NEF" << endl;
    cout << "RMS of bin " << i << ": " << nef_ppHTs[i]->GetRMS() << " <--- pp HT NEF" << endl;
    cout << endl;
    cout << "mean of bin " << i << ": " << m_pAJPs[i]->GetMean() << " <--- pA JP M" << endl;
    cout << "mean of bin " << i << ": " << m_pAHTs[i]->GetMean() << " <--- pA HT M" << endl;
    cout << "mean of bin " << i << ": " << m_ppJPs[i]->GetMean() << " <--- pp JP M" << endl;
    cout << "mean of bin " << i << ": " << m_ppHTs[i]->GetMean() << " <--- pp HT M" << endl;

    cout << "mean of bin " << i << ": " << m_pAJPs_lowEA[i]->GetMean() << " <--- pA JP low EA M" << endl;
    cout << "mean of bin " << i << ": " << m_pAJPs_highEA[i]->GetMean() << " <--- pA JP high EA M" << endl;
    cout << endl;
    cout << "RMS of bin " << i << ": " << m_pAJPs[i]->GetRMS() << " <--- pA JP M" << endl;
    cout << "RMS of bin " << i << ": " << m_pAHTs[i]->GetRMS() << " <--- pA HT M" << endl;
    cout << "RMS of bin " << i << ": " << m_ppJPs[i]->GetRMS() << " <--- pp JP M" << endl;
    cout << "RMS of bin " << i << ": " << m_ppHTs[i]->GetRMS() << " <--- pp HT M" << endl;
    
    cout << "RMS of bin " << i << ": " << m_pAJPs_lowEA[i]->GetRMS() << " <--- pA JP low EA M" << endl;
    cout << "RMS of bin " << i << ": " << m_pAJPs_highEA[i]->GetRMS() << " <--- pA JP high EA M" << endl;
    cout << endl;
  }
  
   
  TLegend *tleg = new TLegend(0.55,0.15,0.75,0.35); tleg->SetBorderSize(0);
  tleg->AddEntry(m_pAJPs_lowEA[0],"pAu JP2 low EA","p");
  tleg->AddEntry(m_pAJPs_highEA[0],"pAu JP2 high EA","p");
  //tleg->AddEntry(nef_pAHTs[0],"pAu HT2","p");
  TLegend *tleg2 = new TLegend(0.55,0.15,0.75,0.35); tleg2->SetBorderSize(0);
  tleg2->AddEntry(m_ppJPs[0],"pp JP2","p");
  //tleg2->AddEntry(nef_ppHTs[0],"pp HT2","p");
  
  TLatex *ttitle = new TLatex(); TLatex *slice = new TLatex();
  /*
  cnef->cd(1);
  ttitle->DrawLatexNDC(0.3,0.5, "anti-k_{T}, R = 0.4, |#eta_{jet}| < 1-R");
  
  for (int i = 0; i < nPts; ++ i) {
    cnef->cd(i+2);
    nef_pAJPs[i]->Draw("same");
    nef_pAHTs[i]->Draw("same");
    nef_ppJPs[i]->Draw("same");
    nef_ppHTs[i]->Draw("same");
    if (i == 0) {tleg->Draw("same");}
    if (i == 1) {tleg2->Draw("same");}
    slice->DrawLatexNDC(0.3,0.85,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
  }
  */
  cm->cd(1);
  ttitle->DrawLatexNDC(0.3,0.5, "anti-k_{T}, R = 0.4, |#eta_{jet}| < 1-R");
  
  for (int i = 0; i < nPts; ++ i) {
    cm->cd(i+2);
    m_pAJPs_lowEA[i]->Draw("same");
    m_pAJPs_highEA[i]->Draw("same");
    //m_pAHTs[i]->Draw("same");
    m_ppJPs[i]->Draw("same");
    //m_ppHTs[i]->Draw("same");
    if (i == 0) {tleg->Draw("same");}
    if (i == 1) {tleg2->Draw("same");}
    slice->DrawLatexNDC(0.3,0.85,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
  }
  
  //cnef->SaveAs("~/jetmass2/plots/pA_pp_HT_JP_compare_nef.pdf");
  //  cm->SaveAs("~/jetmass2/plots/pA_pp_HT_JP_compare_m_by_bbcsum.pdf");
  
  
  
  return;
}
