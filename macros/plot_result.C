#include "Plots_old.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

using namespace std;

//radius format: 4 = 0.4, ...
void plot_result (int radius, bool groom) {
  gStyle->SetPalette(kPastel);
  
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
  TFile *funfold = new TFile(("~/jetmass2/out/unfold/"+fstart+"unfolded"+radstring+".root").c_str(),"READ");
  TFile *fp6 = new TFile(("~/jetmass2/out/sim/hists/unmatched_hists"+radstring+"_bindropped.root").c_str(),"READ");
  TFile *fp8 = new TFile(("~/jetmass2/production/out/pythia/hists/pythia8"+radstring+"_undecayed_hists.root").c_str(),"READ");
  //TFile *fp8 = new TFile(("~/jetmass/production/macros/hists/hists_allsim_lowzgremoved"+radstring+".root").c_str(),"READ");//temporarily using the old files (they're the same, but will point to new ones later)
  TFile *fh7 = new TFile(("~/jetmass/production/macros/hists/hists_allsim_lowzgremoved"+radstring+".root").c_str(),"READ");//temporarily using the old files (they're the same, but will point to new ones later)
  //TFile *fPL = new TFile(("~/jetmass2/production/out/pythia/hists/pythia8"+radstring+"_FSR_only_undecayed_hists.root").c_str(),"READ");
  TFile *fPL = new TFile("~/jetmass2/production/out/pythia/hists/test_nozcut_R04.root","READ");
  
  //get histograms from the files
  vector<TH1D*> recos = {(TH1D*) funfold->Get("nom_0"),(TH1D*) funfold->Get("nom_1")};
  vector<TH1D*> w_systs = {(TH1D*) funfold->Get("w_systs_0"),(TH1D*) funfold->Get("w_systs_1")};
  TH2D* p6_2D = (TH2D*) fp6->Get(("PL_"+hname+"_v_pt").c_str());
  TH2D* p8_2D = (TH2D*) fp8->Get((hname+"vpt").c_str());
  TH2D* h7_2D = (TH2D*) fh7->Get((hname+"vpt_h7off").c_str());//internal name will be changed later
  TH2D* PL_2D = (TH2D*) fPL->Get("PLmvHLpt");
  TH2D* PL_g_2D = (TH2D*) fPL->Get("PLmgvHLpt");
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
  
  for (unsigned i = 0; i < nBins; ++ i) {
    Prettify1D(recos[i],kRed,kFullStar,4,kRed,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.5);
    Prettify1D(w_systs[i],kRed,kFullStar,0,kRed,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.5);
    w_systs[i]->SetFillColor(kRed - 10); w_systs[i]->SetFillStyle(1001);
    Prettify1DwLineStyle(p6s[i], kBlue, kSolid,2, (htitle+" [GeV/c^{2}]").c_str(), ("1/N dN/d"+htitle).c_str(),0,14,0,0.5);
    Prettify1DwLineStyle(p8s[i],kBlack,kSolid,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.5);
    Prettify1DwLineStyle(p8_qs[i],kOrange-1,kSolid,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.5);
    Prettify1DwLineStyle(p8_gs[i],kGreen+2,kSolid,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.5);                        
    Prettify1DwLineStyle(h7s[i],kMagenta,kSolid,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.5);
    Prettify1DwLineStyle(PLs[i],kBlack,kDotted,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.5);
    Prettify1DwLineStyle(PL_gs[i],kBlack,kDashed,2,(htitle+" [GeV/c^{2}]").c_str(),("1/N dN/d"+htitle).c_str(),0,14,0,0.5); 
  }

  TCanvas *cws = new TCanvas("cws","cws",800,500);
  DivideCanvas(cws,"0",nBins,1);

  TLegend *tleg = new TLegend(0.5,0.55,0.75,0.7); tleg->SetBorderSize(0);
  TLegend *tleg2 = new TLegend(0.5,0.55,0.75,0.7); tleg2->SetBorderSize(0);
  TH1D* for_legend = (TH1D*) w_systs[0]->Clone("for_legend"); for_legend->SetMarkerSize(2);
  tleg->AddEntry(for_legend,"Unfolded data","pf");
  tleg->AddEntry(p6s[0],"Pythia6","l");
  tleg2->AddEntry(h7s[0],"Herwig7","l");
  tleg2->AddEntry(p8s[0],"Pythia8","l");
  
  //tleg2->AddEntry(PLs[0],"Pythia8 P jets","l");
  if (groom) {
    //tleg2->AddEntry(PL_gs[0],"Pythia8 groomed P jets","l");
  }
  if (!groom) {
    //tleg2->AddEntry(p8_qs[0],"Pythia8 q jets","l");
    //tleg2->AddEntry(p8_gs[0],"Pythia8 g jets","l");
  }
  
  TLatex *ttitle = new TLatex(); ttitle->SetTextAlign(11); ttitle->SetTextSize(0.05);
  TLatex *slice = new TLatex();

  for (int i = 0; i < nBins; ++ i) {
    cws->cd(i+1);
    w_systs[i]->Draw("E3same");
    recos[i]->Draw("same");
    p6s[i]->Draw("Csame");
    p8s[i]->Draw("Csame");
    h7s[i]->Draw("Csame");
    /*PLs[i]->Draw("Csame");
    if (groom) {
      PL_gs[i]->Draw("Csame");
    }
    if (!groom) {
      p8_qs[i]->Draw("Csame");
      p8_gs[i]->Draw("Csame");
      }*/
    slice->DrawLatexNDC(0.5,0.77,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
    if (i == 0) {ttitle->DrawLatex(4.2,0.46, "pp 200 GeV run12 JP2");ttitle->DrawLatex(4.2,0.43, "anti-k_{T}, R = 0.4");ttitle->DrawLatex(4.2,0.4, "Ch+Ne jets, |#eta| < 0.6"); tleg->Draw("same");}
    if (i == 1) { tleg2->Draw("same");}
  }

  cws->SaveAs(("~/jetmass2/plots/results/"+fstart+"mass_result"+radstring+".pdf").c_str());

  return;
}
