// -*- C++ -*-
//
// Package:    TrkJetAnalyzer
// Class:      TrkJetAnalyzer
// 
/**\class TrkJetAnalyzer TrkJetAnalyzer.cc TrkJetPhase2/TrkJetAnalyzer/plugins/TrkJetAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Livia Soffi
//         Created:  Thu, 23 Jul 2015 12:55:55 GMT
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/TrackJetCollection.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//---- for GenParticles
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
// #include "DataFormats/Candidate/interface/Candidate.h"

//---- for GenJets
#include "DataFormats/JetReco/interface/GenJet.h" 

#include "TTree.h"
#include "TH1.h"
#include "TLorentzVector.h"
//
// class declaration
//


class TrkJetAnalyzer : public edm::EDAnalyzer {
   public:
      explicit TrkJetAnalyzer(const edm::ParameterSet&);
      ~TrkJetAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
  edm::InputTag GenJetCollection_;
  edm::InputTag GenParticlesCollection_;
  edm::InputTag TrackJetCollection_;
  edm::InputTag TrackCollection_;
  edm::InputTag SimTrackCollection_;
  edm::InputTag VertexCollection_;
  edm::InputTag PUSummaryInfoCollection_;
  edm::InputTag PfCandidatesCollection_;

  TTree* myTree_;

  //---Vertex and Pileup
 
  Int_t npu_;
  Float_t pu_zpos_[300];
  Float_t pu_sumpt_lowpt_[300];
  Float_t pu_sumpt_highpt_[300];
  Float_t pu_ntrks_lowpt_[300];
  Float_t pu_ntrks_highpt_[300];

  Float_t vxMC_;
  Float_t vyMC_;
  Float_t vzMC_;
  Int_t nvertex_;
  Float_t vx_[300];
  Float_t vy_[300];
  Float_t vz_[300];
  Float_t vntracks_[300];
  Float_t vchi2_[300];
  Float_t vndof_[300];




  //---- gen jets
  int ngenJet_;
  float genJetPt_[250];
  float genJetEta_[250];
  float genJetPhi_[250];
  float genJetMass_[250]; 
  float genJetEmE_[250]   ;
  float genJetHadrE_[250] ;
  float genJetInvE_[250]  ;
  float genJetAuxE_[250]  ;
  int genJetNconst_[250];

  //------- trk jets
  int ntrkJet_;
  float trkJetPt_[250];
  float trkJetEta_[250];
  float trkJetPhi_[250];
  float trkJetMass_[250]; 
  float trkJetNtrk_[250]; 
 
  //------- tracks
  int ntrk_;
  float trkPt_[9000];
  float trkEta_[9000];
  float trkPhi_[9000];

  //----- gen particles
  int ngenCand_;
  float genCandPt_[30]   ;
  float genCandEta_[30]  ;
  float genCandPhi_[30]  ;
  float genCandMass_[30] ;
  float genCandStatus_[30];
  int genCandPdgId_[30];
  int genCandMothPdgId_[30];
 
  //----- VBF quarks
  int nvbfQuark_;
  float vbfQuarkPt_[30]   ;
  float vbfQuarkEta_[30]  ;
  float vbfQuarkPhi_[30]  ;
  float vbfQuarkMass_[30] ;
  float vbfQuarkStatus_[30];
  int vbfQuarkPdgId_[30];
  int vbfQuarkMothPdgId_[30];


 //----- GGH quarks
  int ngghQuark_;
  float gghQuarkPt_[30]   ;
  float gghQuarkEta_[30]  ;
  float gghQuarkPhi_[30]  ;
  float gghQuarkMass_[30] ;
  float gghQuarkStatus_[30];
  int gghQuarkPdgId_[30];
  int gghQuarkMothPdgId_[30];


  //----- Z
  int nZ_;
  float ZPt_[2]   ;
  float ZEta_[2]  ;
  float ZPhi_[2]  ;
  float ZMass_[2] ;
  float ZStatus_[2];
  int ZPdgId_[2];
 
 //----- Muons
  int nMu_;
  float MuPt_[4]   ;
  float MuEta_[4]  ;
  float MuPhi_[4]  ;
  float MuMass_[4] ;
  float MuStatus_[4];
  int MuPdgId_[4];
 
  int mH_; 
   
  unsigned int minTracks_;
  TH1D *ntrkhisto;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TrkJetAnalyzer::TrkJetAnalyzer(const edm::ParameterSet& iConfig):
  minTracks_(iConfig.getUntrackedParameter<unsigned int>("minTracks",0))
{
   //now do what ever initialization is needed
  // I want to make a histogram of number of tracks in a<dwn event     
  edm::Service<TFileService> fs;

  GenJetCollection_       = iConfig.getParameter<edm::InputTag>("GenJetCollection");
  GenParticlesCollection_ = iConfig.getParameter<edm::InputTag>("GenParticlesCollection");
  TrackJetCollection_       = iConfig.getParameter<edm::InputTag>("TrackJetCollection");
  TrackCollection_       = iConfig.getParameter<edm::InputTag>("TrackCollection");
  //  SimTrackCollection_       = iConfig.getParameter<edm::InputTag>("simTrackCollection");
  VertexCollection_       = iConfig.getParameter<edm::InputTag>("VertexCollection");
  PUSummaryInfoCollection_ = iConfig.getParameter<edm::InputTag>("PUSummaryInfoCollection");
  //PfCandidatesCollection_ = iConfig.getParameter<edm::InputTag>("PfCandidatesCollection");

  myTree_ = fs -> make <TTree>("myTree","myTree");
  myTree_ -> Branch("ngenJet", &ngenJet_, "ngenJet/I");
  myTree_ -> Branch("genJetPt", genJetPt_, "genJetPt[ngenJet]/F");
  myTree_ -> Branch("genJetMass", genJetMass_, "genJet[ngenJet]/F");
  myTree_ -> Branch("genJetEta", genJetEta_, "genJetEta[ngenJet]/F");
  myTree_ -> Branch("genJetPhi", genJetPhi_, "genJetPhi[ngenJet]/F");
  myTree_ -> Branch("genJetEmE", genJetEmE_, "genJeEmEt[ngenJet]/F");
  myTree_ -> Branch("genJetHadrE", genJetHadrE_, "genJetHadrE[ngenJet]/F");
  myTree_ -> Branch("genJetInvE", genJetInvE_, "genJetInvE[ngenJet]/F");
  myTree_ -> Branch("genJetAuxE", genJetAuxE_, "genJetAuxE[ngenJet]/F");
  myTree_ -> Branch("genJetNconst", genJetNconst_, "genJetNconst[ngenJet]/I");
  
  myTree_ -> Branch("ntrkJet", &ntrkJet_, "ntrkJet/I");
  myTree_ -> Branch("trkJetPt", trkJetPt_, "trkJetPt[ntrkJet]/F");
  myTree_ -> Branch("trkJetMass", trkJetMass_, "trkJet[ntrkJet]/F");
  myTree_ -> Branch("trkJetEta", trkJetEta_, "trkJetEta[ntrkJet]/F");
  myTree_ -> Branch("trkJetPhi", trkJetPhi_, "trkJetPhi[ntrkJet]/F");
  myTree_ -> Branch("trkJetNtrk", trkJetNtrk_, "trkJetNtrk[ntrkJet]/F");
 
  myTree_ -> Branch("ntrk", &ntrk_, "ntrk/I");
  myTree_ -> Branch("trkPt", trkPt_, "trkPt[ntrk]/F");
  myTree_ -> Branch("trkEta", trkEta_, "trkEta[ntrk]/F");
  myTree_ -> Branch("trkPhi", trkPhi_, "trkPhi[ntrk]/F");
 

  myTree_ -> Branch("ngenCand",   &ngenCand_,     "ngenCand/I")   ;
  myTree_ -> Branch("genCandPt",   genCandPt_,    "genCandPt[ngenCand]/F")   ;
  myTree_ -> Branch("genCandEta" , genCandEta_,    "genCandEta[ngenCand]/F")  ;
  myTree_ -> Branch("genCandPhi",  genCandPhi_,    "genCandPhi[ngenCand]/F")  ;
  myTree_ -> Branch("genCandMass", genCandMass_,   "genCandMass[ngenCand]/F") ;
  myTree_ -> Branch("genCandStatus", genCandStatus_,  "genCandStatus[ngenCand]/F");
  myTree_ -> Branch("genCandPdgId", genCandPdgId_,  "genCandPdgId[ngenCand]/I");
  myTree_ -> Branch("genCandMothPdgId", genCandMothPdgId_,  "genCandMothPdgId[ngenCand]/I");
 

  myTree_ -> Branch("nvbfQuark",   &nvbfQuark_,     "nvbfQuark/I")   ;
  myTree_ -> Branch("vbfQuarkPt",   vbfQuarkPt_,    "vbfQuarkPt[nvbfQuark]/F")   ;
  myTree_ -> Branch("vbfQuarkEta" , vbfQuarkEta_,    "vbfQuarkEta[nvbfQuark]/F")  ;
  myTree_ -> Branch("vbfQuarkPhi",  vbfQuarkPhi_,    "vbfQuarkPhi[nvbfQuark]/F")  ;
  myTree_ -> Branch("vbfQuarkMass", vbfQuarkMass_,   "vbfQuarkMass[nvbfQuark]/F") ;
  myTree_ -> Branch("vbfQuarkStatus", vbfQuarkStatus_,  "vbfQuarkStatus[nvbfQuark]/F");
  myTree_ -> Branch("vbfQuarkPdgId", vbfQuarkPdgId_,  "vbfQuarkPdgId[nvbfQuark]/I");
  myTree_ -> Branch("vbfQuarkMothPdgId", vbfQuarkMothPdgId_,  "vbfQuarkMothPdgId[nvbfQuark]/I");


  myTree_ -> Branch("ngghQuark",   &ngghQuark_,     "ngghQuark/I")   ;
  myTree_ -> Branch("gghQuarkPt",   gghQuarkPt_,    "gghQuarkPt[ngghQuark]/F")   ;
  myTree_ -> Branch("gghQuarkEta" , gghQuarkEta_,    "gghQuarkEta[ngghQuark]/F")  ;
  myTree_ -> Branch("gghQuarkPhi",  gghQuarkPhi_,    "gghQuarkPhi[ngghQuark]/F")  ;
  myTree_ -> Branch("gghQuarkMass", gghQuarkMass_,   "gghQuarkMass[ngghQuark]/F") ;
  myTree_ -> Branch("gghQuarkStatus", gghQuarkStatus_,  "gghQuarkStatus[ngghQuark]/F");
  myTree_ -> Branch("gghQuarkPdgId", gghQuarkPdgId_,  "gghQuarkPdgId[ngghQuark]/I");
  myTree_ -> Branch("gghQuarkMothPdgId", gghQuarkMothPdgId_,  "gghQuarkMothPdgId[ngghQuark]/I");
  
  myTree_ -> Branch("nZ",   &nZ_,     "nZ/I")   ;
  myTree_ -> Branch("ZPt",   ZPt_,    "ZPt[nZ]/F")   ;
  myTree_ -> Branch("ZEta" , ZEta_,    "ZEta[nZ]/F")  ;
  myTree_ -> Branch("ZPhi",  ZPhi_,    "ZPhi[nZ]/F")  ;
  myTree_ -> Branch("ZMass", ZMass_,   "ZMass[nZ]/F") ;
  myTree_ -> Branch("ZStatus", ZStatus_,  "ZStatus[nZ]/F");
  myTree_ -> Branch("ZPdgId", ZPdgId_,  "ZPdgId[nZ]/I");
  //myTree_ -> Branch("ZMothPdgId", ZMothPdgId_,  "ZMothPdgId[nZ]/I");
 
  myTree_ -> Branch("nMu",   &nMu_,     "nMu/I")   ;
  myTree_ -> Branch("MuPt",   MuPt_,    "MuPt[nMu]/F")   ;
  myTree_ -> Branch("MuEta" , MuEta_,    "MuEta[nMu]/F")  ;
  myTree_ -> Branch("MuPhi",  MuPhi_,    "MuPhi[nMu]/F")  ;
  myTree_ -> Branch("MuMass", MuMass_,   "MuMass[nMu]/F") ;
  myTree_ -> Branch("MuStatus", MuStatus_,  "MuStatus[nMu]/F");
  myTree_ -> Branch("MuPdgId", MuPdgId_,  "MuPdgId[nMu]/I");
  //myTree_ -> Branch("MuMothPdgId", MuMothPdgId_,  "MuMothPdgId[nMu]/I");
  
   
  myTree_->Branch("npu", &npu_, "npu/I");
  myTree_->Branch("pu_zpos", &pu_zpos_, "pu_zpos[npu]/F");
  myTree_->Branch("pu_sumpt_lowpt", &pu_sumpt_lowpt_, "pu_sumpt_lowpt[npu]/F");
  myTree_->Branch("pu_sumpt_highpt", &pu_sumpt_highpt_, "pu_sumpt_highpt[npu]/F");
  myTree_->Branch("pu_ntrks_lowpt", &pu_ntrks_lowpt_, "pu_ntrks_lowpt[npu]/F");
  myTree_->Branch("pu_ntrks_highpt", &pu_ntrks_highpt_, "pu_ntrks_highpt[npu]/F");


  myTree_->Branch("vxMC",&vxMC_,"vxMC/F");
  myTree_->Branch("vyMC",&vyMC_,"vyMC/F");
  myTree_->Branch("vzMC",&vzMC_,"vzMC/F");
  myTree_->Branch("nvertex",&nvertex_,"nvertex/I");
  myTree_->Branch("vx",&vx_,"vx[nvertex]/F");
  myTree_->Branch("vy",&vy_,"vy[nvertex]/F");
  myTree_->Branch("vz",&vz_,"vz[nvertex]/F");
  myTree_->Branch("vntracks",&vntracks_,"vntracks[nvertex]/F");
  myTree_->Branch("vchi2",&vchi2_,"vchi2[nvertex]/F");
  myTree_->Branch("vndof",&vndof_,"vndof[nvertex]/F");

  myTree_->Branch("mH", &mH_, "mH/I");

  ntrkhisto = fs->make<TH1D>("tracks" , "Tracks" , 320 , 2500 , 10000 );

}


TrkJetAnalyzer::~TrkJetAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void TrkJetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   
   
   edm::Handle<reco::GenParticleCollection> genParticles;
   iEvent.getByLabel(GenParticlesCollection_,genParticles);
   
   edm::Handle<reco::GenJetCollection> genJet;
   iEvent.getByLabel(GenJetCollection_,genJet);
   
   edm::Handle<reco::TrackJetCollection> trkJet;
   iEvent.getByLabel(TrackJetCollection_,trkJet);
  
   edm::Handle<reco::TrackCollection> trk;
   iEvent.getByLabel(TrackCollection_,trk);
  
  
   edm::Handle<reco::VertexCollection> vtx;
   iEvent.getByLabel(VertexCollection_, vtx);

   edm::Handle<std::vector<PileupSummaryInfo> > pileup;
   iEvent.getByLabel(PUSummaryInfoCollection_, pileup);

 
      
  //save gen particles info
   ngenCand_ = 0;
   for(int i = 0;i<30;i++){
     genCandPt_[i] =-999.;
     genCandEta_[i]=-999.;
     genCandPhi_[i]=-999.;
     genCandMass_[i]=-999.;
     genCandStatus_[i]=-999.;
     genCandPdgId_[i]=-999;
     genCandMothPdgId_[i]=-999;
   }
   // std::cout<<"----------------------------------------"<<std::endl;
   int k = 0;
   for (reco::GenParticleCollection::const_iterator genCandIter = genParticles->begin(); genCandIter != genParticles->end(); genCandIter++){
     if(k>30)continue;
     genCandPt_[k]= genCandIter->pt();
     genCandEta_[k]= genCandIter->eta();
     genCandPhi_[k]= genCandIter->phi();
     genCandMass_[k]= genCandIter->mass();
     genCandStatus_[k]= genCandIter->status();
     genCandPdgId_[k]= genCandIter->pdgId();
     //  genCandMothPdgId_[k]= genCandIter->mother()->pdgId();
     //     int   pid = genCandIter->pdgId ();
     //  if( (abs(pid) < 10 || abs(pid)==21)) for (unsigned int j=0;  j<genCandIter->numberOfMothers(); j++)  std::cout<< genCandStatus_[k]<<" "<<genCandPdgId_[k]<<"  "<<genCandPt_[k]<<std::endl;
     ngenCand_++;
     k++;
   }
   
   //save info from vbfQuarks
   nvbfQuark_ = 0;
   for(int i = 0;i<30;i++){
     vbfQuarkPt_[i] =-999.;
     vbfQuarkEta_[i]=-999.;
     vbfQuarkPhi_[i]=-999.;
     vbfQuarkMass_[i]=-999.;
     vbfQuarkStatus_[i]=-999.;
     vbfQuarkPdgId_[i]=-999;
     vbfQuarkMothPdgId_[i]=-999;
   }
   int r = 0;
   for (reco::GenParticleCollection::const_iterator genCandIter = genParticles->begin(); genCandIter != genParticles->end(); genCandIter++){
     int status = genCandIter->status ();
     int   pid = genCandIter->pdgId ();
     
     bool isVBFquark = false ; 
     if ( status == 3 && (abs(pid) < 10 || abs(pid)==21) && genCandIter->numberOfMothers() > 0 ) { 
        for (unsigned int j=0; !isVBFquark && j<genCandIter->numberOfMothers(); j++) {
          for (unsigned int k=0; k<genCandIter->mother(j)->numberOfDaughters(); k++) {
            if ( abs(genCandIter->mother(j)->daughter(k)->pdgId()) == 25 ) {
	      vbfQuarkPt_[r]= genCandIter->pt();
	      vbfQuarkEta_[r]= genCandIter->eta();
	      vbfQuarkPhi_[r]= genCandIter->phi();
	      vbfQuarkMass_[r]= genCandIter->mass();
	      vbfQuarkStatus_[r]= genCandIter->status();
	      vbfQuarkPdgId_[r]= genCandIter->pdgId();
	      //   vbfQuarkMothPdgId_[r]= genCandIter->mother()->pdgId();
	      nvbfQuark_++;
	      isVBFquark = true;
	      //	      std::cout<<vbfQuarkPdgId_[r]<<"  "<<vbfQuarkMass_[r]<<"   "<<vbfQuarkEta_[r]<<"   "<<genCandIter->mother(j)->pdgId()<<"  "<<genCandIter->mother(j)->mass()<<std::endl;
	      r++;
	    }
	  }
	}
     }
   }

//save info from gghQuarks
   ngghQuark_ = 0;
   for(int i = 0;i<30;i++){
     gghQuarkPt_[i] =-999.;
     gghQuarkEta_[i]=-999.;
     gghQuarkPhi_[i]=-999.;
     gghQuarkMass_[i]=-999.;
     gghQuarkStatus_[i]=-999.;
     gghQuarkPdgId_[i]=-999;
     gghQuarkMothPdgId_[i]=-999;
   }
   int rr = 0;
   for (reco::GenParticleCollection::const_iterator genCandIter = genParticles->begin(); genCandIter != genParticles->end(); genCandIter++){
     int status = genCandIter->status ();
     int   pid = genCandIter->pdgId ();
     
     //     bool isGGHquark = false ; 
     if ( status == 3 && (abs(pid) < 10 || abs(pid)==21) && genCandIter->numberOfMothers() > 0 ) { 
       //       for (unsigned int j=0; !isGGHquark && j<genCandIter->numberOfMothers(); j++) {
       // for (unsigned int k=0; k<genCandIter->mother(j)->numberOfDaughters(); k++) {
       //  if ( abs(genCandIter->mother(j)->daughter(k)->pdgId()) == 25 ) {

	      gghQuarkPt_[rr]= genCandIter->pt();
	      gghQuarkEta_[rr]= genCandIter->eta();
	      gghQuarkPhi_[rr]= genCandIter->phi();
	      gghQuarkMass_[rr]= genCandIter->mass();
	      gghQuarkStatus_[rr]= genCandIter->status();
	      gghQuarkPdgId_[rr]= genCandIter->pdgId();
	      //   gghQuarkMothPdgId_[r]= genCandIter->mother()->pdgId();
	      ngghQuark_++;
	      //   std::cout<<gghQuarkPdgId_[rr]<<"  "<<gghQuarkMass_[rr]<<"   "<<gghQuarkEta_[rr]<<"   "<<genCandIter->mother(j)->pdgId()<<"  "<<genCandIter->mother(j)->mass()<<std::endl;
	      //	      isGGHquark = true;
	      rr++;
	      //  }
	      //	 }
	      //}
     }
   }

   //reconstruct H->ZZ->4L
   nZ_ = 0;
   for(int i = 0;i<2;i++){
     ZPt_[i] =-999.;
     ZEta_[i]=-999.;
     ZPhi_[i]=-999.;
     ZMass_[i]=-999.;
     ZStatus_[i]=-999.;
     ZPdgId_[i]=-999;
    
   }

   nMu_ = 0;
   for(int i = 0;i<4;i++){
     MuPt_[i] =-999.;
     MuEta_[i]=-999.;
     MuPhi_[i]=-999.;
     MuMass_[i]=-999.;
     MuStatus_[i]=-999.;
     MuPdgId_[i]=-999;
    
   }
  int t = 0;
  int mu = 0;
   for (reco::GenParticleCollection::const_iterator genCandIter = genParticles->begin(); genCandIter != genParticles->end(); genCandIter++){
     int status = genCandIter->status ();
     int   pid = genCandIter->pdgId ();
     
     //save Zs infos
     bool isZ = false ; 
     if (status == 3 && pid == 23 && genCandIter->numberOfMothers() > 0 ) { 
        for (unsigned int j=0; !isZ && j<genCandIter->numberOfMothers(); j++) {
	  if ( abs(genCandIter->mother(j)->pdgId()) == 25 ) {
	      ZPt_[t]= genCandIter->pt();
	      ZEta_[t]= genCandIter->eta();
	      ZPhi_[t]= genCandIter->phi();
	      ZMass_[t]= genCandIter->mass();
	      ZStatus_[t]= genCandIter->status();
	      ZPdgId_[t]= genCandIter->pdgId();
	      //   ZMothPdgId_[t]= genCandIter->mother()->pdgId();
	      nZ_++;
	      t++;
	      isZ=true;
	    }
	  }
	}

     //save Leptons infos
     // bool isMu = false ; 
     if (status == 3 && (abs(pid) == 13||abs(pid)==11||abs(pid)==15)/*  && genCandIter->numberOfMothers() > 0*/ ) { 
       /*for (unsigned int j=0; !isMu && j<genCandIter->numberOfMothers(); j++) {
	  if ( abs(genCandIter->mother(j)->pdgId()) == 23 ) {
	    for (unsigned int k=0; !isMu && j<(genCandIter->mother(j))->numberOfMothers(); k++) {
	      if ( abs(genCandIter->mother(j)->mother(k)->pdgId()) == 25 ) {
       */  
	      MuPt_[mu]= genCandIter->pt();
	      MuEta_[mu]= genCandIter->eta();
	      MuPhi_[mu]= genCandIter->phi();
	      MuMass_[mu]= genCandIter->mass();
	      MuStatus_[mu]= genCandIter->status();
	      MuPdgId_[mu]= genCandIter->pdgId();
	      //   MuMothPdgId_[mu]= genCandIter->mother()->pdgId();
	      nMu_++;
	      mu++;
	      //   isMu =true;
	      /* }
	    }
	  }
	  }*/
     }
     
   }

   if(nMu_>3){
     TLorentzVector higgs;
    
       TLorentzVector* tmu0=new TLorentzVector();
       TLorentzVector* tmu1=new TLorentzVector();
       TLorentzVector* tmu2=new TLorentzVector();
       TLorentzVector* tmu3=new TLorentzVector();
       tmu0->SetPtEtaPhiM(MuPt_[0],MuEta_[0],MuPhi_[0], MuMass_[0]);
       tmu1->SetPtEtaPhiM(MuPt_[1],MuEta_[1],MuPhi_[1], MuMass_[1]);
       tmu2->SetPtEtaPhiM(MuPt_[2],MuEta_[2],MuPhi_[2], MuMass_[2]);
       tmu3->SetPtEtaPhiM(MuPt_[3],MuEta_[3],MuPhi_[3], MuMass_[3]);
       higgs= *tmu0+*tmu1+*tmu2+*tmu3;
       mH_ =higgs.M();
   } 


   //save pu infos
   npu_ = 0;
 
   if( pileup.isValid() ) 
     {
       std::vector<PileupSummaryInfo>::const_iterator PVI;       
       for(PVI = pileup->begin(); PVI != pileup->end(); ++PVI) 
	 {
	   if(PVI->getBunchCrossing() != 0) 
	     continue;
	   npu_ = PVI->getPU_NumInteractions();
	   int sv = PVI->getPU_zpositions().size() < 50 ? PVI->getPU_zpositions().size() : 300;
	   for (int iPU=0;iPU<sv;++iPU)
	     {
	       pu_zpos_[iPU]=PVI->getPU_zpositions()[iPU];
	       pu_sumpt_lowpt_[iPU]=PVI->getPU_sumpT_lowpT()[iPU];
	       pu_sumpt_highpt_[iPU]=PVI->getPU_sumpT_highpT()[iPU];
	       pu_ntrks_lowpt_[iPU]=PVI->getPU_ntrks_lowpT()[iPU];
	       pu_ntrks_highpt_[iPU]=PVI->getPU_ntrks_highpT()[iPU];
	     }
	 }
     }
   
   
  // Get the primary vertex coordinates
 
   nvertex_=0;
   for (reco::VertexCollection::const_iterator it = vtx->begin(); 
	it != vtx->end(); ++it) {
     
     vx_[nvertex_] = (it->isValid()) ? it->x() : 999.;
     vy_[nvertex_] = (it->isValid()) ? it->y() : 999.;
     vz_[nvertex_] = (it->isValid()) ? it->z() : 999.;
     
     vntracks_[nvertex_] = (it->isValid()) ? it->tracksSize() : 0;
     vchi2_[nvertex_] = (it->isValid()) ? it->normalizedChi2() : 100.;
     vndof_[nvertex_] = (it->isValid()) ? it->ndof() : 0.;
     
     nvertex_++;
   }
   
   //save gen jet informations
   ngenJet_ = 0;
   for(int i = 0;i<250;i++){
     genJetPt_[i] =-999.;
     genJetEta_[i]=-999.;
     genJetPhi_[i]=-999.;
     genJetMass_[i]=-999.;
     genJetEmE_[i]=-999.;
     genJetHadrE_[i]=-999.;
     genJetInvE_[i]=-999.;
     genJetAuxE_[i]=-999.;
     genJetNconst_[i]=-999;
 }
   int j=0;
   for (reco::GenJetCollection::const_iterator genJetIter=genJet->begin(); genJetIter!=genJet->end(); genJetIter++){
     if(genJetIter->pt()<5)continue; 
     genJetPt_[j]=genJetIter->pt();
     genJetEta_[j]=genJetIter->eta();
     genJetPhi_[j]=genJetIter->phi();
     genJetMass_[j]=genJetIter->mass();
     genJetEmE_[j]=genJetIter->emEnergy();
     genJetHadrE_[j]  =genJetIter->hadEnergy();
     genJetInvE_[j]=genJetIter->invisibleEnergy();
     genJetAuxE_[j]=genJetIter->auxiliaryEnergy();
     genJetNconst_[j]=genJetIter->nConstituents();
     ngenJet_++;
     
     //  std::cout<<genJetPt_[j]<<"  "<<genJetEta_[j]<<"   "<<genJetPhi_[j]<<std::endl;
     j++;
   }
   
   //save trk jet infos
   ntrkJet_=0;
   for(int l = 0;l<250;l++){
     trkJetPt_[l] =-999.;
     trkJetEta_[l]=-999.;
     trkJetPhi_[l]=-999.;
     trkJetMass_[l]=-999.;
     trkJetNtrk_[l]=-999.;
 }
   int z=0;
   for (reco::TrackJetCollection::const_iterator trkJetIter=trkJet->begin(); trkJetIter!=trkJet->end(); trkJetIter++){
     //     if(trkJetIter->pt()<10)continue; 
     trkJetPt_[z]=trkJetIter->pt();
     trkJetEta_[z]=trkJetIter->eta();
     trkJetPhi_[z]=trkJetIter->phi();
     trkJetMass_[z]=trkJetIter->mass();
     trkJetNtrk_[z]=trkJetIter->tracks().size();
     z++;
     ntrkJet_++;
   }
  
   //save trk  infos
   ntrk_=0;
   for(int ll = 0;ll<9000;ll++){
     trkPt_[ll] =-999.;
     trkEta_[ll]=-999.;
     trkPhi_[ll]=-999.;
   }
   int zz=0;
   for (reco::TrackCollection::const_iterator trkIter=trk->begin(); trkIter!=trk->end(); trkIter++){
     //if(trkIter->pt()<5)continue; 
     trkPt_[zz]=trkIter->pt();
     trkEta_[zz]=trkIter->eta();
     trkPhi_[zz]=trkIter->phi();
     zz++;
     ntrk_++;
   }
  
   //look only at H->ZZ->4Mu
   if(nMu_>=0) myTree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
TrkJetAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrkJetAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
TrkJetAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
TrkJetAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
TrkJetAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
TrkJetAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrkJetAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrkJetAnalyzer);
