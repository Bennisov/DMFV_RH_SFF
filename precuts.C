#define precuts_cxx
#include "precuts.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>

void precuts::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   std::cout<<nentries<<std::endl;
   double cs = 0.001646; //pb
   double L = 21.9; //fb^-1
   cs *= 1000; //fbs
   double w = cs * L / nentries;

   TFile *file = new TFile("../bck/signal.root", "RECREATE");
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
   Float_t         MissingETMET;
   Float_t         MissingETPHI;
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
   tree->Branch("MissingETMET", &MissingETMET);
   tree->Branch("MissingETPHI", &MissingETPHI);
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
      for (size_t i = 0; i < Jet_size; i++)
      {
         if (Jet_BTag[i]) b_jets++;
      }
      if (b_jets < 1) continue;
      ElectronE.clear();
      ElectronEta.clear();
      ElectronPT.clear();
      ElectronPhi.clear();
      JetBTag.clear();
      JetE.clear();
      JetEta.clear();
      JetPT.clear();
      JetPhi.clear();
      MissingETMET = 0.;
      MissingETPHI = 0.;
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
         e.SetPtEtaPhiM(1000*Electron_PT[i], Electron_Eta[i], Electron_Phi[i], 0.511);
         ElectronE.push_back(e.E());
         ElectronPhi.push_back(e.Phi());
         ElectronPT.push_back(e.Pt());
         ElectronEta.push_back(e.Eta());
      }
      for (size_t i = 0; i < Muon_size; i++)
      {
         TLorentzVector mu;
         mu.SetPtEtaPhiM(1000*Muon_PT[i], Muon_Eta[i], Muon_Phi[i], 105.66);
         MuonE.push_back(mu.E());
         MuonPhi.push_back(mu.Phi());
         MuonPT.push_back(mu.Pt());
         MuonEta.push_back(mu.Eta());
      }
      for (size_t i = 0; i < Jet_size; i++)
      {
         TLorentzVector jet;
         jet.SetPtEtaPhiM(1000*Jet_PT[i], Jet_Eta[i], Jet_Phi[i], 1000*Jet_Mass[i]);
         JetE.push_back(jet.E());
         JetPhi.push_back(jet.Phi());
         JetPT.push_back(jet.Pt());
         JetEta.push_back(jet.Eta());
         JetBTag.push_back(Jet_BTag[i]);
      }
      MissingETMET = 1000*MissingET_MET[0];
      MissingETPHI = 1000*MissingET_Phi[0];

      eventNumber = jentry;
      met_met_NOSYS = 1000*MissingET_MET[0];
      met_phi_NOSYS = 1000*MissingET_Phi[0];
      met_significance_NOSYS = 0.;
      met_sumet_NOSYS = 0.;
      //std::cout<<w<<std::endl;
      weight_total_NOSYS = w;

      tree->Fill();
   }
   tree->Write();

   file->Close();
}
