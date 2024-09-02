#define precuts_cxx
#include "precuts.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>

bool GoodJet(Int_t size_e, Float_t* Eta_e, Float_t* Phi_e, Float_t jet_eta, Float_t jet_phi)
{
   bool flag = 1;
   Float_t dR_cut_lept = 0.2;
   for (size_t i = 0; i < size_e; i++)
   {
      double dR = std::sqrt((Eta_e[i] - jet_eta)*(Eta_e[i] - jet_eta)+(Phi_e[i] - jet_phi)*(Phi_e[i] - jet_phi));
      if (dR < dR_cut_lept) flag = 0;
   }
   return flag;
}

bool GoodMuon(Int_t size_j, Float_t* Eta_j, Float_t* Phi_j, Float_t mu_eta, Float_t mu_phi)
{
   bool flag = 1;
   Float_t dR_cut_lept = 0.2;
   for (size_t i = 0; i < size_j; i++)
   {
      double dR = std::sqrt((Eta_j[i] - mu_eta)*(Eta_j[i] - mu_eta)+(Phi_j[i] - mu_phi)*(Phi_j[i] - mu_phi));
      if (dR < dR_cut_lept) flag = 0;
   }
   return flag;
}

void precuts::Loop()
{


   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   double cs = 0.7047;
   cs *= 1000;
   double L = 29.1;
   Long64_t nbytes = 0, nb = 0;
   double w = cs * L / nentries;

   TFile *file = new TFile("300_2_precut2.root", "RECREATE");
   TTree *tree = new TTree("reco", "reco");

   vector<float>   ElectronE;
   vector<float>   ElectronEta;
   vector<float>   ElectronPT;
   vector<float>   ElectronPhi;
   vector<char>    JetBTag;
   vector<float>   JetE;
   vector<float>   JetEta;
   vector<float>   JetPT;
   vector<float>   JetPhi;
   Float_t         ETmiss_NOSYS;
   Float_t         ETmissPhi_NOSYS;
   vector<float>   MuonE;
   vector<float>   MuonEta;
   vector<float>   MuonPT;
   vector<float>   MuonPhi;
   ULong64_t       eventNumber;
   Float_t         met_met_NOSYS;
   Float_t         met_phi_NOSYS;
   Float_t         met_significance_NOSYS;
   Float_t         met_sumet_NOSYS;
   Double_t        weight_total_NOSYS;

   tree->Branch("ElectronE", &ElectronE);
   tree->Branch("ElectronEta", &ElectronEta);
   tree->Branch("ElectronPT", &ElectronPT);
   tree->Branch("ElectronPhi", &ElectronPhi);
   tree->Branch("JetBTag", &JetBTag);
   tree->Branch("JetE", &JetE);
   tree->Branch("JetEta", &JetEta);
   tree->Branch("JetPT", &JetPT);
   tree->Branch("JetPhi", &JetPhi);
   tree->Branch("ETmiss_NOSYS", &ETmiss_NOSYS);
   tree->Branch("ETmissPhi_NOSYS", &ETmissPhi_NOSYS);
   tree->Branch("MuonE", &MuonE);
   tree->Branch("MuonEta", &MuonEta);
   tree->Branch("MuonPT", &MuonPT);
   tree->Branch("MuonPhi", &MuonPhi);
   tree->Branch("eventNumber", &eventNumber);
   tree->Branch("met_met_NOSYS", &met_met_NOSYS);
   tree->Branch("met_phi_NOSYS", &met_phi_NOSYS);
   tree->Branch("met_significance_NOSYS", &met_significance_NOSYS);
   tree->Branch("met_sumet_NOSYS", &met_sumet_NOSYS);
   tree->Branch("weight_total_NOSYS", &weight_total_NOSYS);





   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (MissingET_MET[0] < 100.) continue;
      if (Jet_size < 2) continue;
      int b_jets = 0;
      bool var = 1;
      bool var2 = 1;
      for (size_t i = 0; i < Jet_size; i++)
      {
         if (Jet_BTag[i]) 
         {
            b_jets++;
            if (var)
            {
               var = 0;
               if (Jet_PT[i] < 30.) continue;               
            }
         }
         else if (var2)
         {
            var2 = 0;
            if (Jet_PT[i] < 30.) continue;
         }
      }
      if (b_jets < 1) continue;
      if (!(Electron_size + Muon_size)) continue;
      if (Electron_size)
      {
         if (Muon_size)
         {
            if (Electron_PT[0]>Muon_PT[0])
            {
               if(Electron_PT[0] < 30.) continue;
            }
            else 
            {
               if(Muon_PT[0] < 30.) continue;
            }
         }
         else
         {
            if(Electron_PT[0] < 30.) continue;
         } 
      }
      else 
      {
         if(Muon_PT[0] < 30.) continue;
      }

      ElectronE.clear();
      ElectronEta.clear();
      ElectronPT.clear();
      ElectronPhi.clear();
      JetBTag.clear();
      JetE.clear();
      JetEta.clear();
      JetPT.clear();
      JetPhi.clear();
      ETmiss_NOSYS = 0.;
      ETmissPhi_NOSYS = 0.;
      MuonE.clear();
      MuonEta.clear();
      MuonPT.clear();
      MuonPhi.clear();
      eventNumber = 0;
      met_met_NOSYS = 0.;
      met_phi_NOSYS = 0.;
      met_significance_NOSYS = 0.;
      met_sumet_NOSYS = 0.;
      weight_total_NOSYS = 0.;
      for (size_t i = 0; i < Electron_size; i++)
      {
         TLorentzVector e;
         e.SetPtEtaPhiM(Electron_PT[i], Electron_Eta[i], Electron_Phi[i], 0.000511);
         ElectronE.push_back(e.E());
         ElectronPhi.push_back(e.Phi());
         ElectronPT.push_back(e.Pt());
         ElectronEta.push_back(e.Eta());
      }
      for (size_t i = 0; i < Muon_size; i++)
      {
         if (GoodMuon(Jet_size, Jet_Eta, Jet_Phi, Muon_Eta[i], Muon_Phi[i])) 
         {
            TLorentzVector mu;
            mu.SetPtEtaPhiM(Muon_PT[i], Muon_Eta[i], Muon_Phi[i], 0.10566);
            MuonE.push_back(mu.E());
            MuonPhi.push_back(mu.Phi());
            MuonPT.push_back(mu.Pt());
            MuonEta.push_back(mu.Eta());
         }
      }
      for (size_t i = 0; i < Jet_size; i++)
      {
         if (GoodJet(Electron_size, Electron_Eta, Electron_Phi, Jet_Eta[i], Jet_Phi[i]))
         {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(Jet_PT[i], Jet_Eta[i], Jet_Phi[i], Jet_Mass[i]);
            JetE.push_back(jet.E());
            JetPhi.push_back(jet.Phi());
            JetPT.push_back(jet.Pt());
            JetEta.push_back(jet.Eta());
            JetBTag.push_back(Jet_BTag[i]);
         }
      }
      ETmiss_NOSYS =  MissingET_MET[0];
      ETmissPhi_NOSYS = MissingET_Phi[0];

      eventNumber = jentry;
      met_met_NOSYS = MissingET_MET[0];
      met_phi_NOSYS = MissingET_Phi[0];
      met_significance_NOSYS = 0.;
      met_sumet_NOSYS = 0.;
      //std::cout<<w<<std::endl;
      weight_total_NOSYS = w;

      tree->Fill();
   }
   tree->Write();

   file->Close();
}
