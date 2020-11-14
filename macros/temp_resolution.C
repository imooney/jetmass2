void temp_resolution () {
  
  TFile* matchFile = new TFile( "~/jetmass2_11-10-2020_11_10-2020/out/sim/sim_matched_R04_no_masscut.root", "READ");

  TTree *t = (TTree*) matchFile->Get("event");
  
  double weight = 1;
  vector<double> *gePt = 0; vector<double> *geM = 0;
  vector<double> *pyPt = 0; vector<double> *pyM = 0;
  
  t->SetBranchAddress("p_M",&pyM);
  t->SetBranchAddress("g_M",&geM);
  t->SetBranchAddress("p_Pt",&pyPt);
  t->SetBranchAddress("g_Pt",&gePt);

  TH3D* res_v_m_v_pt = new TH3D("res_v_m_v_pt","",50,0,2,14,0,14,9,15,60);
  
  for (int i = 0; i < t->GetEntries(); ++ i) {
    t->GetEntry(i);
    for(int j = 0; j < pyM->size(); ++ j) {
      res_v_m_v_pt->Fill(geM->at(j) / (double) pyM->at(j), geM->at(j), gePt->at(j), weight);
    }
  }
  
  t->ResetBranchAddresses();
  
  TFile *fout = new TFile("~/jetmass2_11-10-2020_11_10-2020/out/sim/hists/temp_resolution_no_mass_cut.root","RECREATE");
  
  fout->cd();
  
  res_v_m_v_pt->Write();
  
  fout->Close();
  
}
