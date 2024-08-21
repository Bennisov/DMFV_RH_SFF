#define cuts_cxx
#include "cuts.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "mt2.h"
#include <iostream>
#include <fstream>

bool GoodJet(Int_t size_e, Int_t size_m, Int_t size_p, Float_t* Eta_e, Float_t* Phi_e, Float_t* Eta_m, Float_t* Phi_m, Float_t* Eta_p, Float_t* Phi_p, Float_t jet_eta, Float_t jet_phi)
{
   bool flag = 1;
   Float_t dR_cut_lept = 0.2;
   for (size_t i = 0; i < size_e; i++)
   {
      double dR = std::sqrt((Eta_e[i] - jet_eta)*(Eta_e[i] - jet_eta)+(Phi_e[i] - jet_phi)*(Phi_e[i] - jet_phi));
      if (dR < dR_cut_lept) flag = 0;
   }
   for (size_t i = 0; i < size_m; i++)
   {
      double dR = std::sqrt((Eta_m[i] - jet_eta)*(Eta_m[i] - jet_eta)+(Phi_m[i] - jet_phi)*(Phi_m[i] - jet_phi));
      if (dR < dR_cut_lept) flag = 0;
   }
   for (size_t i = 0; i < size_p; i++)
   {
      double dR = std::sqrt((Eta_p[i] - jet_eta)*(Eta_p[i] - jet_eta)+(Phi_p[i] - jet_phi)*(Phi_p[i] - jet_phi));
      if (dR < dR_cut_lept) flag = 0;
   }
   return flag;
}

void cuts::Loop()
{
   //std::ofstream out;
   //out.open("../NN.mt2/sig.csv");
   asymm_mt2_lester_bisect::disableCopyrightMessage();
   TFile *outputFile = new TFile("output_1.0_1000.root", "RECREATE");
   auto h_mt2 = new TH1F("h_mt2", ";;", 16, 200., 1000.);
   auto cutflow = new TH1F("cutflow", "Cutflow", 12, 0.5, 0.5 + 12);
   auto h_mt = new TH1F("h_mt", ";;", 20, 0., 1500.);
   auto h_met = new TH1F("h_met", ";;", 14, 200., 900.);
   auto h_pt_l = new TH1F("h_pt_l", ";;", 20, 0., 1500.);
   auto h_n_bjet = new TH1F("h_bjet_n", ";;", 10, 0, 10);
   auto h_pt_bjet = new TH1F("h_pt_bjet", ";;", 20, 0., 1500.);
   auto h_pt_jet = new TH1F("h_pt_jet", ";;", 20, 0., 1500.);
   auto h_m_bl = new TH1F("h_m_bl", ";;", 20, 0., 1500.);
   auto h_dphi_min = new TH1F("h_dphi_min", ";;", 50, 0., 6.3);
   auto h_n_l = new TH1F("h_n_l", ";;", 10, 0, 10);
   auto h_n_jets = new TH1F("h_n_jets", ";;", 10, 0, 10);

   cutflow->GetXaxis()->SetBinLabel(1, "No Cut");
   cutflow->GetXaxis()->SetBinLabel(2, "At least 1 lepton detected");
   cutflow->GetXaxis()->SetBinLabel(3, "Lepton PT>30GeV");
   cutflow->GetXaxis()->SetBinLabel(4, "At least 1 B-jet");
   cutflow->GetXaxis()->SetBinLabel(5, "Lead B-Jet PT>30GeV");
   cutflow->GetXaxis()->SetBinLabel(6, "At least 2 jets");
   cutflow->GetXaxis()->SetBinLabel(7, "Lead Jet PT>100GeV");
   cutflow->GetXaxis()->SetBinLabel(8, "Inv mass of B-jet and l < 160GeV");
   cutflow->GetXaxis()->SetBinLabel(9, "Trans mass of l and met > 160GeV");
   cutflow->GetXaxis()->SetBinLabel(10, "Dphi min > 0.6");
   cutflow->GetXaxis()->SetBinLabel(11, "Missing MET > 90GeV");
   cutflow->GetXaxis()->SetBinLabel(12, "mt2 of bjet,l,jet > X"); 

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   double cs = 0.001646; //pb
   double L = 300.; //fb^-1
   cs *= 1000; //fb
   double w = cs * L / nentries;

   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      int l_n = Electron_size + Muon_size;
      cutflow->Fill(1, w);

      if (l_n < 1) continue;
      cutflow->Fill(2, w);

      TLorentzVector good_lept;
      if (Electron_size)
      {
         if (Muon_size)
         {
            if (Electron_PT[0] > Muon_PT[0]) good_lept.SetPtEtaPhiM(Electron_PT[0], Electron_Eta[0], Electron_Phi[0], 0.000511);
            else good_lept.SetPtEtaPhiM(Muon_PT[0], Muon_Eta[0], Muon_Phi[0], 0.10566);
         }
         else good_lept.SetPtEtaPhiM(Electron_PT[0], Electron_Eta[0], Electron_Phi[0], 0.000511);
      }
      else 
      {
         good_lept.SetPtEtaPhiM(Muon_PT[0], Muon_Eta[0], Muon_Phi[0], 0.10566);
      }

      if (good_lept.Pt() < 30.) continue;
      cutflow->Fill(3, w);

      int jet_count = 0;
      int bjet_count = 0;
      float dphi_min = 999.;
      TLorentzVector bjet;
      TLorentzVector lead_jet;
      bool var = 1;
      std::vector<float> j_pt;
      for (size_t i = 0; i < Jet_size; i++)
      {
         float dphi = std::abs(Jet_Phi[i] - MissingET_Phi[0]);
         if (dphi < dphi_min)
         {
            dphi_min = dphi;
         }
         if (GoodJet(Electron_size, Muon_size, Photon_size, Electron_Eta, Electron_Phi, Muon_Eta, Muon_Phi, Photon_Eta, Photon_Phi, Jet_Eta[i], Jet_Phi[i]))
         {
            jet_count++;
            j_pt.push_back(Jet_PT[i]);
            if (Jet_BTag[i])
            {
               bjet_count++;
               bjet.SetPtEtaPhiM(Jet_PT[i], Jet_Eta[i], Jet_Phi[i], Jet_Mass[i]);
            }
            else if (var)
            {
               lead_jet.SetPtEtaPhiM(Jet_PT[i], Jet_Eta[i], Jet_Phi[i], Jet_Mass[i]);
               var = 0;
            }
         }
      }
      if (bjet_count < 1) continue;
      cutflow->Fill(4, w);

      if (bjet.Pt() < 30.) continue;
      cutflow->Fill(5, w);

      if (jet_count < 2) continue;
      cutflow->Fill(6, w);

      if (lead_jet.Pt() < 100.) continue;
      cutflow->Fill(7, w);

      float kos = std::cos(std::abs(good_lept.Phi() - MissingET_Phi[0]));
      float mt = std::sqrt(2 * MissingET_MET[0] * good_lept.Pt() * (1 - kos));
      float mbl = (bjet+good_lept).M();

      if (mbl < 160.) continue;
      cutflow->Fill(8, w);

      if (mt > 160.) continue;
      cutflow->Fill(9, w);

      if (dphi_min < 0.6) continue;
      cutflow->Fill(10, w);

      if (MissingET_MET[0] < 90.) continue;
      cutflow->Fill(11, w);


      TLorentzVector tqark = bjet + good_lept;
      double mVisA = tqark.M();
      double mVisB = lead_jet.M();
      double pxA = tqark.Px();
      double pyA = tqark.Py();
      double pxB = lead_jet.Px();
      double pyB = lead_jet.Py();
      double pxMiss = MissingET_MET[0] * std::cos(MissingET_Phi[0]);
      double pyMiss = MissingET_MET[0] * std::sin(MissingET_Phi[0]);
      double chiA = 200.;
      double chiB = 200.;
      double desiredPrecisionOnMt2 = 0.;
      float MT2 = -999.;
      MT2 =  asymm_mt2_lester_bisect::get_mT2(
      mVisA, pxA, pyA,
      mVisB, pxB, pyB,
      pxMiss, pyMiss,
      chiA, chiB,
      desiredPrecisionOnMt2);
      if (MT2 < 400.) continue;
      cutflow->Fill(12, w);

      h_mt2->Fill(MT2, w);
      h_pt_bjet->Fill(bjet.Pt(), w);
      h_m_bl->Fill(mbl, w);
      for (float i : j_pt) h_pt_jet->Fill(i, w);
      h_n_bjet->Fill(bjet_count, w);
      h_n_jets->Fill(jet_count, w);
      h_dphi_min->Fill(dphi_min, w);
      h_pt_l->Fill(good_lept.Pt(), w);
      h_n_l->Fill(l_n, w);
      h_met->Fill(MissingET_MET[0], w);
      h_mt->Fill(mt, w);
   }
   h_mt2->Write();
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
   cutflow->Write();
   outputFile->Close();
   //out.close();
}
