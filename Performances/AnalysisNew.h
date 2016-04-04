//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Sep 14 16:40:26 2015 by ROOT version 5.34/07
// from TTree myTree/myTree
// found on file: ../../../src/VFPixSimulation/TrkJetAnalyzer/python/histoOUTPUT_VBF.root
//////////////////////////////////////////////////////////

#ifndef AnalysisNew_h
#define AnalysisNew_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class AnalysisNew {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           ngenJet;
   Float_t         genJetPt[50];   //[ngenJet]
   Float_t         genJetMass[50];   //[ngenJet]
   Float_t         genJetEta[50];   //[ngenJet]
   Float_t         genJetPhi[50];   //[ngenJet]
   Float_t         genJetEmE[50];   //[ngenJet]
   Float_t         genJetHadrE[50];   //[ngenJet]
   Float_t         genJetInvE[50];   //[ngenJet]
   Float_t         genJetAuxE[50];   //[ngenJet]
   Int_t           genJetNconst[50];   //[ngenJet]
   Int_t           ntrkJet;
   Float_t         trkJetPt[40];   //[ntrkJet]
   Float_t         trkJetMass[40];   //[ntrkJet]
   Float_t         trkJetEta[40];   //[ntrkJet]
   Float_t         trkJetPhi[40];   //[ntrkJet]
   Int_t           ntrk;
   Float_t         trkPt[5799];   //[ntrk]
   Float_t         trkEta[5799];   //[ntrk]
   Float_t         trkPhi[5799];   //[ntrk]
   Int_t           ngenCand;
   Float_t         genCandPt[31];   //[ngenCand]
   Float_t         genCandEta[31];   //[ngenCand]
   Float_t         genCandPhi[31];   //[ngenCand]
   Float_t         genCandMass[31];   //[ngenCand]
   Float_t         genCandStatus[31];   //[ngenCand]
   Int_t           genCandPdgId[31];   //[ngenCand]
   Int_t           genCandMothPdgId[31];   //[ngenCand]
   Int_t           nvbfQuark;
   Float_t         vbfQuarkPt[6];   //[nvbfQuark]
   Float_t         vbfQuarkEta[6];   //[nvbfQuark]
   Float_t         vbfQuarkPhi[6];   //[nvbfQuark]
   Float_t         vbfQuarkMass[6];   //[nvbfQuark]
   Float_t         vbfQuarkStatus[6];   //[nvbfQuark]
   Int_t           vbfQuarkPdgId[6];   //[nvbfQuark]
   Int_t           vbfQuarkMothPdgId[6];   //[nvbfQuark]
   Int_t           nZ;
   Float_t         ZPt[2];   //[nZ]
   Float_t         ZEta[2];   //[nZ]
   Float_t         ZPhi[2];   //[nZ]
   Float_t         ZMass[2];   //[nZ]
   Float_t         ZStatus[2];   //[nZ]
   Int_t           ZPdgId[2];   //[nZ]
   Int_t           nMu;
   Float_t         MuPt[4];   //[nMu]
   Float_t         MuEta[4];   //[nMu]
   Float_t         MuPhi[4];   //[nMu]
   Float_t         MuMass[4];   //[nMu]
   Float_t         MuStatus[4];   //[nMu]
   Int_t           MuPdgId[4];   //[nMu]
   Int_t           npu;
   Float_t         pu_zpos[157];   //[npu]
   Float_t         pu_sumpt_lowpt[157];   //[npu]
   Float_t         pu_sumpt_highpt[157];   //[npu]
   Float_t         pu_ntrks_lowpt[157];   //[npu]
   Float_t         pu_ntrks_highpt[157];   //[npu]
   Float_t         vxMC;
   Float_t         vyMC;
   Float_t         vzMC;
   Int_t           nvertex;
   Float_t         vx[103];   //[nvertex]
   Float_t         vy[103];   //[nvertex]
   Float_t         vz[103];   //[nvertex]
   Float_t         vntracks[103];   //[nvertex]
   Float_t         vchi2[103];   //[nvertex]
   Float_t         vndof[103];   //[nvertex]
   Int_t           mH;

   // List of branches
   TBranch        *b_ngenJet;   //!
   TBranch        *b_genJetPt;   //!
   TBranch        *b_genJetMass;   //!
   TBranch        *b_genJetEta;   //!
   TBranch        *b_genJetPhi;   //!
   TBranch        *b_genJetEmE;   //!
   TBranch        *b_genJetHadrE;   //!
   TBranch        *b_genJetInvE;   //!
   TBranch        *b_genJetAuxE;   //!
   TBranch        *b_genJetNconst;   //!
   TBranch        *b_ntrkJet;   //!
   TBranch        *b_trkJetPt;   //!
   TBranch        *b_trkJetMass;   //!
   TBranch        *b_trkJetEta;   //!
   TBranch        *b_trkJetPhi;   //!
   TBranch        *b_ntrk;   //!
   TBranch        *b_trkPt;   //!
   TBranch        *b_trkEta;   //!
   TBranch        *b_trkPhi;   //!
   TBranch        *b_ngenCand;   //!
   TBranch        *b_genCandPt;   //!
   TBranch        *b_genCandEta;   //!
   TBranch        *b_genCandPhi;   //!
   TBranch        *b_genCandMass;   //!
   TBranch        *b_genCandStatus;   //!
   TBranch        *b_genCandPdgId;   //!
   TBranch        *b_genCandMothPdgId;   //!
   TBranch        *b_nvbfQuark;   //!
   TBranch        *b_vbfQuarkPt;   //!
   TBranch        *b_vbfQuarkEta;   //!
   TBranch        *b_vbfQuarkPhi;   //!
   TBranch        *b_vbfQuarkMass;   //!
   TBranch        *b_vbfQuarkStatus;   //!
   TBranch        *b_vbfQuarkPdgId;   //!
   TBranch        *b_vbfQuarkMothPdgId;   //!
   TBranch        *b_nZ;   //!
   TBranch        *b_ZPt;   //!
   TBranch        *b_ZEta;   //!
   TBranch        *b_ZPhi;   //!
   TBranch        *b_ZMass;   //!
   TBranch        *b_ZStatus;   //!
   TBranch        *b_ZPdgId;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_MuPt;   //!
   TBranch        *b_MuEta;   //!
   TBranch        *b_MuPhi;   //!
   TBranch        *b_MuMass;   //!
   TBranch        *b_MuStatus;   //!
   TBranch        *b_MuPdgId;   //!
   TBranch        *b_npu;   //!
   TBranch        *b_pu_zpos;   //!
   TBranch        *b_pu_sumpt_lowpt;   //!
   TBranch        *b_pu_sumpt_highpt;   //!
   TBranch        *b_pu_ntrks_lowpt;   //!
   TBranch        *b_pu_ntrks_highpt;   //!
   TBranch        *b_vxMC;   //!
   TBranch        *b_vyMC;   //!
   TBranch        *b_vzMC;   //!
   TBranch        *b_nvertex;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_vntracks;   //!
   TBranch        *b_vchi2;   //!
   TBranch        *b_vndof;   //!
   TBranch        *b_mH;   //!

   AnalysisNew(TTree *tree=0);
   virtual ~AnalysisNew();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(std::string isFit, std::string isDelphes);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef AnalysisNew_cxx
AnalysisNew::AnalysisNew(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/cmshome/lsoffi/CMSSW_6_2_0_SLHC23_patch1/src/rootfiles/output2.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/cmshome/lsoffi/CMSSW_6_2_0_SLHC23_patch1/src/rootfiles/output2.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/cmshome/lsoffi/CMSSW_6_2_0_SLHC23_patch1/src/rootfiles/output2.root:/TrkJetPhase2");
      dir->GetObject("myTree",tree);

   }
   Init(tree);
}

AnalysisNew::~AnalysisNew()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AnalysisNew::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AnalysisNew::LoadTree(Long64_t entry)
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

void AnalysisNew::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ngenJet", &ngenJet, &b_ngenJet);
   fChain->SetBranchAddress("genJetPt", genJetPt, &b_genJetPt);
   fChain->SetBranchAddress("genJetMass", genJetMass, &b_genJetMass);
   fChain->SetBranchAddress("genJetEta", genJetEta, &b_genJetEta);
   fChain->SetBranchAddress("genJetPhi", genJetPhi, &b_genJetPhi);
   fChain->SetBranchAddress("genJetEmE", genJetEmE, &b_genJetEmE);
   fChain->SetBranchAddress("genJetHadrE", genJetHadrE, &b_genJetHadrE);
   fChain->SetBranchAddress("genJetInvE", genJetInvE, &b_genJetInvE);
   fChain->SetBranchAddress("genJetAuxE", genJetAuxE, &b_genJetAuxE);
   fChain->SetBranchAddress("genJetNconst", genJetNconst, &b_genJetNconst);
   fChain->SetBranchAddress("ntrkJet", &ntrkJet, &b_ntrkJet);
   fChain->SetBranchAddress("trkJetPt", trkJetPt, &b_trkJetPt);
   fChain->SetBranchAddress("trkJetMass", trkJetMass, &b_trkJetMass);
   fChain->SetBranchAddress("trkJetEta", trkJetEta, &b_trkJetEta);
   fChain->SetBranchAddress("trkJetPhi", trkJetPhi, &b_trkJetPhi);
   fChain->SetBranchAddress("ntrk", &ntrk, &b_ntrk);
   fChain->SetBranchAddress("trkPt", trkPt, &b_trkPt);
   fChain->SetBranchAddress("trkEta", trkEta, &b_trkEta);
   fChain->SetBranchAddress("trkPhi", trkPhi, &b_trkPhi);
   fChain->SetBranchAddress("ngenCand", &ngenCand, &b_ngenCand);
   fChain->SetBranchAddress("genCandPt", genCandPt, &b_genCandPt);
   fChain->SetBranchAddress("genCandEta", genCandEta, &b_genCandEta);
   fChain->SetBranchAddress("genCandPhi", genCandPhi, &b_genCandPhi);
   fChain->SetBranchAddress("genCandMass", genCandMass, &b_genCandMass);
   fChain->SetBranchAddress("genCandStatus", genCandStatus, &b_genCandStatus);
   fChain->SetBranchAddress("genCandPdgId", genCandPdgId, &b_genCandPdgId);
   fChain->SetBranchAddress("genCandMothPdgId", genCandMothPdgId, &b_genCandMothPdgId);
   fChain->SetBranchAddress("nvbfQuark", &nvbfQuark, &b_nvbfQuark);
   fChain->SetBranchAddress("vbfQuarkPt", vbfQuarkPt, &b_vbfQuarkPt);
   fChain->SetBranchAddress("vbfQuarkEta", vbfQuarkEta, &b_vbfQuarkEta);
   fChain->SetBranchAddress("vbfQuarkPhi", vbfQuarkPhi, &b_vbfQuarkPhi);
   fChain->SetBranchAddress("vbfQuarkMass", vbfQuarkMass, &b_vbfQuarkMass);
   fChain->SetBranchAddress("vbfQuarkStatus", vbfQuarkStatus, &b_vbfQuarkStatus);
   fChain->SetBranchAddress("vbfQuarkPdgId", vbfQuarkPdgId, &b_vbfQuarkPdgId);
   fChain->SetBranchAddress("vbfQuarkMothPdgId", vbfQuarkMothPdgId, &b_vbfQuarkMothPdgId);
   fChain->SetBranchAddress("nZ", &nZ, &b_nZ);
   fChain->SetBranchAddress("ZPt", ZPt, &b_ZPt);
   fChain->SetBranchAddress("ZEta", ZEta, &b_ZEta);
   fChain->SetBranchAddress("ZPhi", ZPhi, &b_ZPhi);
   fChain->SetBranchAddress("ZMass", ZMass, &b_ZMass);
   fChain->SetBranchAddress("ZStatus", ZStatus, &b_ZStatus);
   fChain->SetBranchAddress("ZPdgId", ZPdgId, &b_ZPdgId);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("MuPt", MuPt, &b_MuPt);
   fChain->SetBranchAddress("MuEta", MuEta, &b_MuEta);
   fChain->SetBranchAddress("MuPhi", MuPhi, &b_MuPhi);
   fChain->SetBranchAddress("MuMass", MuMass, &b_MuMass);
   fChain->SetBranchAddress("MuStatus", MuStatus, &b_MuStatus);
   fChain->SetBranchAddress("MuPdgId", MuPdgId, &b_MuPdgId);
   fChain->SetBranchAddress("npu", &npu, &b_npu);
   fChain->SetBranchAddress("pu_zpos", pu_zpos, &b_pu_zpos);
   fChain->SetBranchAddress("pu_sumpt_lowpt", pu_sumpt_lowpt, &b_pu_sumpt_lowpt);
   fChain->SetBranchAddress("pu_sumpt_highpt", pu_sumpt_highpt, &b_pu_sumpt_highpt);
   fChain->SetBranchAddress("pu_ntrks_lowpt", pu_ntrks_lowpt, &b_pu_ntrks_lowpt);
   fChain->SetBranchAddress("pu_ntrks_highpt", pu_ntrks_highpt, &b_pu_ntrks_highpt);
   fChain->SetBranchAddress("vxMC", &vxMC, &b_vxMC);
   fChain->SetBranchAddress("vyMC", &vyMC, &b_vyMC);
   fChain->SetBranchAddress("vzMC", &vzMC, &b_vzMC);
   fChain->SetBranchAddress("nvertex", &nvertex, &b_nvertex);
   fChain->SetBranchAddress("vx", vx, &b_vx);
   fChain->SetBranchAddress("vy", vy, &b_vy);
   fChain->SetBranchAddress("vz", vz, &b_vz);
   fChain->SetBranchAddress("vntracks", vntracks, &b_vntracks);
   fChain->SetBranchAddress("vchi2", vchi2, &b_vchi2);
   fChain->SetBranchAddress("vndof", vndof, &b_vndof);
   fChain->SetBranchAddress("mH", &mH, &b_mH);
   Notify();
}

Bool_t AnalysisNew::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AnalysisNew::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AnalysisNew::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef AnalysisNew_cxx
