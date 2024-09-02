#define load_cxx
#include "load.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "../DMRH/mt2.h"

// float sig_calc(float s, float b, float b_err_abs)
// {
//    float tot = s+b;
//    float b2 = b*b;
//    float b_err2 = b_err_abs*b_err_abs;
//    float b_plus_err2 = b+b_err2;
//    float retval = std::sqrt(2*((tot) * std::log(tot*b_plus_err2/(b2+tot*b_err2))-b2/b_err2 * std::log(1+b_err2*s/(b*b_plus_err2))));
//    return retval;
// }


void load::Loop()
{
   asymm_mt2_lester_bisect::disableCopyrightMessage();
   //std::string output_file_name = loaded_file_name.substr(0, loaded_file_name.find(".")) + "_output.root";
   //TFile *outputFile = new TFile(output_file_name.c_str(), "RECREATE");
   TFile *outputFile = new TFile("300_2_histo.root", "RECREATE");

   auto h_mt2 = new TH1F("h_mt2", ";;", 16, 200., 1000.);
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
   auto h_pt_ljet = new TH1F("h_pt_ljet", ";;", 20, 0., 1500.);
   auto h_dphi_bl = new TH1F("h_dphi_bl", ";;", 50, 0., 6.3);
   auto h_dphi_bm = new TH1F("h_dphi_bm", ";;", 50, 0., 6.3);
   auto h_dr_bjetl = new TH1F("h_dr_bjetl", ";;", 50, 0., 6.3);


   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      float mt2 = -999.;

      int l_n = ElectronE->size() + MuonE->size();

      TLorentzVector good_lept;
      if (ElectronE->size())
      {
         if (MuonE->size())
         {
            if (ElectronPT->at(0)>MuonPT->at(0)) good_lept.SetPtEtaPhiM(ElectronPT->at(0), ElectronEta->at(0), ElectronPhi->at(0), 0.000511);
            else good_lept.SetPtEtaPhiM(MuonPT->at(0), MuonEta->at(0), MuonPhi->at(0), 0.10566);
         }
         else good_lept.SetPtEtaPhiM(ElectronPT->at(0), ElectronEta->at(0), ElectronPhi->at(0), 0.000511);
      }
      else 
      {
         good_lept.SetPtEtaPhiM(MuonPT->at(0), MuonEta->at(0), MuonPhi->at(0), 0.10566);
      }

      int jet_count = 0;
      int bjet_count = 0;
      float dphi_min = 999.;

      TLorentzVector bjet;
      TLorentzVector lead_jet;
      bool var1 = 1;
      bool var2 = 1;
      std::vector<float> j_pt;
      for (size_t i = 0; i < JetE->size(); i++)
      {
         float dphi = std::abs(JetPhi->at(i) - ETmissPhi_NOSYS);
         if (dphi < dphi_min)
         {
            dphi_min = dphi;
         }
         jet_count++;
         j_pt.push_back(JetPT->at(i));
         if (JetBTag->at(i))
         {
            bjet_count++;
            if (var2) 
            {
               bjet.SetPtEtaPhiE(JetPT->at(i), JetEta->at(i), JetPhi->at(i), JetE->at(i));
               var2 = 0;
            }
         }
         else if (var1)
         {
            lead_jet.SetPtEtaPhiE(JetPT->at(i), JetEta->at(i), JetPhi->at(i), JetE->at(i));
            var1 = 0;
         }
      }

      float kos = std::cos(std::abs(good_lept.Phi() - ETmissPhi_NOSYS));
      float mt = std::sqrt(2 * ETmiss_NOSYS * good_lept.Pt() * (1 - kos));
      float mbl = (bjet+good_lept).M();
      float dphi_bl = std::abs(bjet.Phi() - good_lept.Phi());
      float dphi_bm = std::abs(bjet.Phi() - ETmissPhi_NOSYS);
      float dr_bjetl = std::sqrt((bjet.Phi() - good_lept.Phi())*(bjet.Phi() - good_lept.Phi()) + (bjet.Eta() - good_lept.Eta())*(bjet.Eta() - good_lept.Eta()));

      


      TLorentzVector tqark = bjet + good_lept;
      double mVisA = tqark.M();
      double mVisB = lead_jet.M();
      double pxA = tqark.Px();
      double pyA = tqark.Py();
      double pxB = lead_jet.Px();
      double pyB = lead_jet.Py();
      double pxMiss = ETmiss_NOSYS * std::cos(ETmissPhi_NOSYS);
      double pyMiss = ETmiss_NOSYS * std::sin(ETmissPhi_NOSYS);
      double chiA = 0.;
      double chiB = 80.;
      double desiredPrecisionOnMt2 = 0.;
      mt2 =  asymm_mt2_lester_bisect::get_mT2(
      mVisA, pxA, pyA,
      mVisB, pxB, pyB,
      pxMiss, pyMiss,
      chiA, chiB,
      desiredPrecisionOnMt2);

      h_mt2->Fill(mt2, weight_total_NOSYS);
      h_pt_bjet->Fill(bjet.Pt(), weight_total_NOSYS);
      h_m_bl->Fill(mbl, weight_total_NOSYS);
      h_n_bjet->Fill(bjet_count, weight_total_NOSYS);
      h_n_jets->Fill(jet_count, weight_total_NOSYS);
      h_dphi_min->Fill(dphi_min, weight_total_NOSYS);
      h_pt_l->Fill(good_lept.Pt(), weight_total_NOSYS);
      h_n_l->Fill(l_n, weight_total_NOSYS);
      h_met->Fill(ETmiss_NOSYS, weight_total_NOSYS);
      h_mt->Fill(mt, weight_total_NOSYS);
      h_pt_ljet->Fill(lead_jet.Pt(), weight_total_NOSYS);
      h_dphi_bl->Fill(dphi_bl, weight_total_NOSYS);
      h_dphi_bm->Fill(dphi_bm, weight_total_NOSYS);
      h_dr_bjetl->Fill(dr_bjetl, weight_total_NOSYS);
   }


   // float cuts = cutflow->GetBinContent(11);
   // std::cout<<cuts;
   // float nocuts = cutflow->GetBinContent(7);

   

   // std::cout<<"Significance cuts: "<<sig_calc(cuts, 26044., 0.3 * 26044.)<<std::endl;
   // std::cout<<"Significance nocuts: "<<sig_calc(nocuts, 26044., 0.3 * 26044.)<<std::endl;

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
   h_pt_ljet->Write();
   h_dphi_bl->Write();
   h_dphi_bm->Write();
   h_dr_bjetl->Write();
   outputFile->Close();
}
