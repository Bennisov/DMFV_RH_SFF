//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Aug 13 12:23:22 2024 by ROOT version 6.30/04
// from TTree reco/reco
// found on file: Data_0_2022_data.root
//////////////////////////////////////////////////////////

#ifndef load_h
#define load_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class load {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   std::string loaded_file_name = "others.root";
      /* 
      "others.root"
      "singletop.root"
      "ttbar.root"
      "ttZ.root"
      "wjets.root"
      "zjets.root"  
      "data.root"
       */
   // Declaration of leaf types
   vector<float>   *ElectronE;
   vector<float>   *ElectronEta;
   vector<float>   *ElectronPT;
   vector<float>   *ElectronPhi;
   vector<char>    *JetBTag;
   vector<float>   *JetE;
   vector<float>   *JetEta;
   vector<float>   *JetPT;
   vector<float>   *JetPhi;
   Float_t         MissingETMET;
   Float_t         MissingETPHI;
   vector<float>   *MuonE;
   vector<float>   *MuonEta;
   vector<float>   *MuonPT;
   vector<float>   *MuonPhi;
   ULong64_t       eventNumber;
   Float_t         met_met_NOSYS;
   Float_t         met_phi_NOSYS;
   Float_t         met_significance_NOSYS;
   Float_t         met_sumet_NOSYS;
   Double_t        weight_total_NOSYS;

   // List of branches
   TBranch        *b_ElectronE;   //!
   TBranch        *b_ElectronEta;   //!
   TBranch        *b_ElectronPT;   //!
   TBranch        *b_ElectronPhi;   //!
   TBranch        *b_JetBTag;   //!
   TBranch        *b_JetE;   //!
   TBranch        *b_JetEta;   //!
   TBranch        *b_JetPT;   //!
   TBranch        *b_JetPhi;   //!
   TBranch        *b_MissingETMET;   //!
   TBranch        *b_MissingETPHI;   //!
   TBranch        *b_MuonE;   //!
   TBranch        *b_MuonEta;   //!
   TBranch        *b_MuonPT;   //!
   TBranch        *b_MuonPhi;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_met_met_NOSYS;   //!
   TBranch        *b_met_phi_NOSYS;   //!
   TBranch        *b_met_significance_NOSYS;   //!
   TBranch        *b_met_sumet_NOSYS;   //!
   TBranch        *b_weight_total_NOSYS;   //!

   load(TTree *tree=0);
   virtual ~load();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef load_cxx
load::load(TTree *tree) : fChain(0) 
{
   if (tree == 0) {
      #ifdef SINGLE_TREE
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("reco",tree);
#else 
      TChain * chain = new TChain("reco","");
      chain->Add(loaded_file_name.c_str());
      tree = chain;
#endif // SINGLE_TREE
   }
   Init(tree);
}

load::~load()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t load::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t load::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void load::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   ElectronE = 0;
   ElectronEta = 0;
   ElectronPT = 0;
   ElectronPhi = 0;
   JetBTag = 0;
   JetE = 0;
   JetEta = 0;
   JetPT = 0;
   JetPhi = 0;
   MuonE = 0;
   MuonEta = 0;
   MuonPT = 0;
   MuonPhi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ElectronE", &ElectronE, &b_ElectronE);
   fChain->SetBranchAddress("ElectronEta", &ElectronEta, &b_ElectronEta);
   fChain->SetBranchAddress("ElectronPT", &ElectronPT, &b_ElectronPT);
   fChain->SetBranchAddress("ElectronPhi", &ElectronPhi, &b_ElectronPhi);
   fChain->SetBranchAddress("JetBTag", &JetBTag, &b_JetBTag);
   fChain->SetBranchAddress("JetE", &JetE, &b_JetE);
   fChain->SetBranchAddress("JetEta", &JetEta, &b_JetEta);
   fChain->SetBranchAddress("JetPT", &JetPT, &b_JetPT);
   fChain->SetBranchAddress("JetPhi", &JetPhi, &b_JetPhi);
   fChain->SetBranchAddress("MissingETMET", &MissingETMET, &b_MissingETMET);
   fChain->SetBranchAddress("MissingETPHI", &MissingETPHI, &b_MissingETPHI);
   fChain->SetBranchAddress("MuonE", &MuonE, &b_MuonE);
   fChain->SetBranchAddress("MuonEta", &MuonEta, &b_MuonEta);
   fChain->SetBranchAddress("MuonPT", &MuonPT, &b_MuonPT);
   fChain->SetBranchAddress("MuonPhi", &MuonPhi, &b_MuonPhi);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("met_met_NOSYS", &met_met_NOSYS, &b_met_met_NOSYS);
   fChain->SetBranchAddress("met_phi_NOSYS", &met_phi_NOSYS, &b_met_phi_NOSYS);
   fChain->SetBranchAddress("met_significance_NOSYS", &met_significance_NOSYS, &b_met_significance_NOSYS);
   fChain->SetBranchAddress("met_sumet_NOSYS", &met_sumet_NOSYS, &b_met_sumet_NOSYS);
   fChain->SetBranchAddress("weight_total_NOSYS", &weight_total_NOSYS, &b_weight_total_NOSYS);
   Notify();
}

Bool_t load::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void load::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t load::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef load_cxx
