#define load_cxx
#include "load.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "../DMRH/mt2.h"

bool GoodJet(Int_t size_e, Int_t size_m, vector<float>* Eta_e, vector<float>* Phi_e, vector<float>* Eta_m, vector<float>* Phi_m, Float_t jet_eta, Float_t jet_phi)
{
   bool flag = 1;
   Float_t dR_cut_lept = 0.2;
   for (size_t i = 0; i < size_e; i++)
   {
      double dR = std::sqrt((Eta_e->at(i) - jet_eta)*(Eta_e->at(i) - jet_eta)+(Phi_e->at(i) - jet_phi)*(Phi_e->at(i) - jet_phi));
      if (dR < dR_cut_lept) flag = 0;
   }
   for (size_t i = 0; i < size_m; i++)
   {
      double dR = std::sqrt((Eta_m->at(i) - jet_eta)*(Eta_m->at(i) - jet_eta)+(Phi_m->at(i) - jet_phi)*(Phi_m->at(i) - jet_phi));
      if (dR < dR_cut_lept) flag = 0;
   }
   return flag;
}

void load::Loop()
{
   asymm_mt2_lester_bisect::disableCopyrightMessage();
   std::string output_file_name = loaded_file_name.substr(0, loaded_file_name.find(".")) + "_output.root";
   TFile *outputFile = new TFile(output_file_name.c_str(), "RECREATE");

   auto h_mt2 = new TH1F("h_mt2", ";;", 50, 0., 1500.);
   // auto cutflow = new TH1F("cutflow", "Cutflow", 6, 0.5, 0.5 + 6);
   auto h_mt = new TH1F("h_mt", ";;", 50, 0., 1500.);
   auto h_met = new TH1F("h_met", ";;", 50, 0., 1500.);
   auto h_pt_l = new TH1F("h_pt_l", ";;", 50, 0., 1500.);
   auto h_n_bjet = new TH1F("h_bjet_n", ";;", 10, 0, 10);
   auto h_pt_bjet = new TH1F("h_pt_bjet", ";;", 50, 0., 1500.);
   auto h_pt_jet = new TH1F("h_pt_jet", ";;", 50, 0., 1500.);
   auto h_m_bl = new TH1F("h_m_bl", ";;", 50, 0., 1500.);
   auto h_dphi_min = new TH1F("h_dphi_min", ";;", 50, 0., 6.3);
   auto h_n_l = new TH1F("h_n_l", ";;", 10, 0, 10);
   auto h_n_jets = new TH1F("h_n_jets", ";;", 10, 0, 10);

   // cutflow->GetXaxis()->SetBinLabel(1, "No Cut");
   // cutflow->GetXaxis()->SetBinLabel(2, "Exactly 1 lepton detected");
   // cutflow->GetXaxis()->SetBinLabel(3, "At least 2 jets");
   // cutflow->GetXaxis()->SetBinLabel(4, "MET > 175GeV");
   // cutflow->GetXaxis()->SetBinLabel(5, "MT of l and MET > 175GeV");
   // cutflow->GetXaxis()->SetBinLabel(6, "At least one B-jet");
   double cs = 29720; //fb
   double u_cs = 87.57; //fb
   double N = 10000;
   double L = 140; //fb^-1
   double w = cs * L / N;




   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      for (auto& element : *ElectronE) element /= 1000.;
      for (auto& element : *ElectronPT) element /= 1000.;
      for (auto& element : *JetE) element /= 1000.;
      for (auto& element : *JetPT) element /= 1000.;
      MissingETMET /= 1000.;
      for (auto& element : *MuonE) element /= 1000.;
      for (auto& element : *MuonPT) element /= 1000.;
      float mt2 = -999.;
      int l_n = ElectronE->size() + MuonE->size();
      if ((l_n != 1) or (MissingETMET < 200.)) continue;
      TLorentzVector good_lept;
      if (ElectronE->size())
      {
         good_lept.SetPtEtaPhiM(ElectronPT->at(0), ElectronEta->at(0), ElectronPhi->at(0), 0.000511);
      }
      else 
      {
         good_lept.SetPtEtaPhiM(MuonPT->at(0), MuonEta->at(0), MuonPhi->at(0), 0.10566);
      }
      if (good_lept.Pt() < 30.) continue;
      int jet_count = 0;
      int bjet_count = 0;
      float dphi_min = 999.;
      TLorentzVector bjet;
      TLorentzVector lead_jet;
      bool var = 1;
      std::vector<float> j_pt;
      for (size_t i = 0; i < JetE->size(); i++)
      {
         float dphi = std::abs(JetPhi->at(i) - MissingETPHI);
         if (dphi < dphi_min)
         {
            dphi_min = dphi;
         }
         if (GoodJet(ElectronE->size(), MuonE->size(), ElectronEta, ElectronPhi, MuonEta, MuonPhi, JetEta->at(i), JetPhi->at(i)))
         {
            jet_count++;
            j_pt.push_back(JetPT->at(i));
            if (JetBTag->at(i))
            {
               bjet_count++;
               bjet.SetPtEtaPhiE(JetPT->at(i), JetEta->at(i), JetPhi->at(i), JetE->at(i));
            }
            else if (var)
            {
               lead_jet.SetPtEtaPhiE(JetPT->at(i), JetEta->at(i), JetPhi->at(i), JetE->at(i));
               var = 0;
            }
         }
      }
      if ((bjet_count != 1) or (jet_count < 2)) continue;
      for(float i: j_pt) h_pt_jet->Fill(i, weight_total_NOSYS);
      h_pt_bjet->Fill(bjet.Pt(), weight_total_NOSYS);
      h_m_bl->Fill((bjet+good_lept).M(), weight_total_NOSYS);
      h_n_bjet->Fill(bjet_count, weight_total_NOSYS);
      h_n_jets->Fill(jet_count, weight_total_NOSYS);
      h_dphi_min->Fill(dphi_min, weight_total_NOSYS);
      h_pt_l->Fill(good_lept.Pt(), weight_total_NOSYS);

      h_n_l->Fill(l_n, weight_total_NOSYS);
      h_met->Fill(MissingETMET, weight_total_NOSYS);
      float kos = std::cos(std::abs(good_lept.Phi() - MissingETPHI));
      float mt = std::sqrt(2 * MissingETMET * good_lept.Pt() * (1 - kos));
      h_mt->Fill(mt, weight_total_NOSYS);

      TLorentzVector tqark = bjet + good_lept;
      double mVisA = tqark.M();
      double mVisB = lead_jet.M();
      double pxA = tqark.Px();
      double pyA = tqark.Py();
      double pxB = lead_jet.Px();
      double pyB = lead_jet.Py();
      double pxMiss = MissingETMET * std::cos(MissingETPHI);
      double pyMiss = MissingETMET * std::sin(MissingETPHI);
      double chiA = 200.;
      double chiB = 200.;
      double desiredPrecisionOnMt2 = 0.;
      mt2 =  asymm_mt2_lester_bisect::get_mT2(
      mVisA, pxA, pyA,
      mVisB, pxB, pyB,
      pxMiss, pyMiss,
      chiA, chiB,
      desiredPrecisionOnMt2);
      h_mt2->Fill(mt2, weight_total_NOSYS);
   }
   h_mt2->Write();
   // cutflow->Write();
   h_mt->Write();
   h_met->Write();
   h_pt_l->Write();
   h_n_bjet->Write();
   h_pt_bjet->Write();
   h_pt_jet->Write();
   h_m_bl->Write();
   h_dphi_min->Write();
   h_n_l->Write();
   h_n_jets->Write();
   outputFile->Close();
}
