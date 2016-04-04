#define AnalysisNew_cxx
#include "AnalysisNew.h"

#include "TCanvas.h"
#include "THStack.h"
#include "TStyle.h"
#include <stdio.h>
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "THStack.h"
#include "TProfile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"                                                                                                                                                         
#include "TLegend.h"                                                                                                                                                        
#include "TPaveText.h"                                                                                                                                                      
#include "TColor.h"                                                                                                                                                         
#include "TLatex.h"                                                                                                                                                         
#include "TLorentzVector.h"
#include "TEfficiency.h"
void AnalysisNew::Loop(std::string isFit, std::string isDelphes)
{ bool isMatched[300];
   Double_t trkJetDR[300];
   Double_t vbfQuarkAssocId[300];

   TH1F* h_nvbfQuark = new TH1F("h_nvbfQuark","", 20,-0.5,20.5);
   TH1F* h_nTrkJet = new TH1F("h_nTrkJet","", 20,-0.5,20.5);

   TH1F* h_nGenJet = new TH1F("h_nGenJet","", 20,-0.5,20.5);
   TH1F* h_genJetPt = new TH1F("h_genJetPt","", 30,0,400);
   TH1F* h_genJetEta = new TH1F("h_genJetEta","", 60,-6,6);
   TH1F* h_genJetPhi = new TH1F("h_genJetPhi","", 60,-4,4);

   TH1F* h_nMatchJet = new TH1F("h_nMatchJet","", 20,-0.5,20.5);

   TH1F* h_ntrk = new TH1F("h_ntrk","", 2000,-0.5,2000.5);
   

   //trk jet efficiency
   TH1F* h_den = new TH1F("h_den","", 20,0,4);
   TH1F* h_num = new TH1F("h_den","", 20,0,4);

   TH1F* h_vbfQuarkPt = new TH1F("h_vbfQuarkPt","", 30,0,400);
   TH1F* h_vbfQuarkEta = new TH1F("h_vbfQuarkEta","", 60,-6,6);
   TH1F* h_vbfQuarkPhi = new TH1F("h_vbfQuarkPhi","", 60,-4,4);
   TH1F* h_vbfQuarkMass = new TH1F("h_vbfQuarkMass","", 30,0,50);

   TH1F* h_trkJetPt = new TH1F("h_trkJetPt","", 30,0,400);
   TH1F* h_trkJetEta = new TH1F("h_trkJetEta","", 60,-6,6);
   TH1F* h_trkJetPhi = new TH1F("h_trkJetPhi","", 60,-4,4);
   TH1F* h_trkJetMass = new TH1F("h_trkJetMass","", 60,0,100);
   TH1F* h_trkJetDR = new TH1F("h_trkJetDR","", 30,0,5);

   TH2F* h2_pt_eta = new TH2F("h2_pt_eta","",60,0,6,30, 0., 400);
   TH2F* h2_phi_eta = new TH2F("h2_phi_eta","",60,0,6,60, -4., 4);
   TH2F* h2_resp_eta = new TH2F("h2_resp_eta","",10,0,4,20, 0., 2);
   

   TH1F* h_trkPt = new TH1F("h_trkPt","", 60,0,80);
   TH1F* h_trkEta = new TH1F("h_trkEta","", 60,-6,6);
   TH1F* h_trkPhi = new TH1F("h_trkPhi","", 60,-4,4);
   
   TH1F* h_trkInJetPt = new TH1F("h_trkInJetPt","", 60,0,80);
   TH2F* h_trkInJetPtEta = new TH2F("h_trkInJetPtEta","", 60,0,80, 60,-6,6);
   TH1F* h_trkInJetP = new TH1F("h_trkInJetP","", 60,0,80);
   TH1F* h_trkInJetEta = new TH1F("h_trkInJetEta","", 60,-6,6);
   TH1F* h_trkInJetPhi = new TH1F("h_trkInJetPhi","", 60,-4,4);
   TH1F* h_ntrkInJet = new TH1F("h_ntrkInJet","", 200,-0.5,200.5);


   //compute vbf/ggh separation variables
   TH1F* h_nJetsNoOverl = new TH1F("h_nJetsNoOverl","", 10,-0.5,10.5);
   TH1F* h_DeltaEtajj = new TH1F("h_deltaEtajj","", 20,0,6);
   TH1F* h_Mjj = new TH1F("h_Mjj","", 80,0,1000);
   TH1F* h_Jet1Pt = new TH1F("h_Jet1Pt","", 260,0,500);
   TH1F* h_Jet1Ntrk = new TH1F("h_Jet1Ntrk","", 60,0,50);
   TH1F* h_Jet2Ntrk = new TH1F("h_Jet2Ntrk","", 60,0,50);
   TH1F* h_Jet2Pt = new TH1F("h_Jet2Pt","", 260,0,500);
   TH1F* h_Jet1Eta = new TH1F("h_Jet1Eta","", 60,-5,5);
   TH1F* h_Jet2Eta = new TH1F("h_Jet2Eta","", 60,-5,5);
   TH1F* h_Jet1Phi = new TH1F("h_Jet1Phi","", 260,-4,4);
   TH1F* h_Jet2Phi = new TH1F("h_Jet2Phi","", 260,-4,4);
   TH1F* h_Djj = new TH1F("h_Djj","", 40,-1,4);


   int counter=0;
  
   TFile* f_out = new TFile(("output_"+isFit+"_"+isDelphes+".root").c_str(), "RECREATE");
   TTree* treeout = new TTree(("tree_"+isFit).c_str(),("tree_"+isFit).c_str());
   double deltaEta_;
   double M_;

   treeout->Branch("deltaEta", &deltaEta_, "deltaEta/D");
   treeout->Branch("M", &M_, "M/D");

   Long64_t nentries = fChain->GetEntriesFast();
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     
        if(jentry%1000==0)std::cout<<jentry<<std::endl;
	
      //match gen and trk jets

	double ngenJ=0;
	for(int pp = 0; pp< ngenJet;pp++){
	    if(genJetPt[pp]<20)continue;
	    ngenJ++;
	    h_genJetPt->Fill(genJetPt[pp]);
	    h_genJetEta->Fill(genJetEta[pp]);
	    h_genJetPhi->Fill(genJetPhi[pp]);
	  }
	h_nGenJet->Fill(ngenJ);

	int ntrkJ=0;

     	for(int pp = 0; pp< ntrkJet;pp++){
	    if(trkJetPt[pp]<10)continue;
	    h_trkJetPt->Fill(trkJetPt[pp]);
	    h_trkJetEta->Fill(trkJetEta[pp]);
	    h_trkJetPhi->Fill(trkJetPhi[pp]);
	    h_trkJetMass->Fill(trkJetMass[pp]);
	    h_trkJetDR->Fill(trkJetDR[pp]);
	    ntrkJ++;
	  }
	h_nTrkJet->Fill(ntrkJ);
	//	std::cout<<"----------------------------------------------"<<std::endl;
	bool matchWvbfQ;
	if(isFit=="GGH")matchWvbfQ = false;
	else  matchWvbfQ = false;
	for(int pp = 0; pp< ntrkJet;pp++){
	  if(trkJetPt[pp]<10)continue;
	  
	  double dRMin = 999.;
	  int idQuark=999;
	  double dRMinGen = 999.;
	  int idQuarkGen=999;
	  double dRMinGenJet = 999.;
	  int idQuarkGenJet=999;
	  isMatched[pp]=false;
	  trkJetDR[pp]=999.;
	  //	  std::cout<<"jet: "<<pp<<" Pt: "<<trkJetPt[pp]<<" Eta: "<<trkJetEta[pp]<<" Phi: "<<trkJetPhi[pp]<<" Mass:"<<trkJetMass[pp]<<std::endl;
	   
	  if(matchWvbfQ==true){	
	    for(int quark = 0; quark< nvbfQuark;quark++){
	      double dRJet = 999.;
	      TLorentzVector* trkP4 = new TLorentzVector();
	      trkP4->SetPtEtaPhiM(trkJetPt[pp],trkJetEta[pp],trkJetPhi[pp],trkJetMass[pp]);
	      TLorentzVector* quarkP4 = new TLorentzVector();
	      quarkP4->SetPtEtaPhiM(vbfQuarkPt[quark],vbfQuarkEta[quark],vbfQuarkPhi[quark],vbfQuarkMass[quark]);
	      dRJet = trkP4->DeltaR(*quarkP4);
	      
	      //compute alternative DR
	      double dRJetNew = sqrt(pow(vbfQuarkEta[quark] - trkJetEta[pp],2) + pow(vbfQuarkPhi[quark] - trkJetPhi[pp],2) );
	      
	      if (dRJet<dRMin){	    
		dRMin=dRJet;
		idQuark = quark; 
		//	if(dRJetNew<0.5)std::cout<<" jet: "<<pp<<" q: "<<idQuark<<" dr New : "<<dRJetNew<<" dr Old : "<<dRJet<<std::endl;
	      }
	      
	    }
	    
	  }else{
	    for(int q = 0; q< ngenCand;q++){
	      double dRJet = 999.;
	      if(genCandPt[q]<10||genCandStatus[q]!=3 || (abs(genCandPdgId[q])> 10 && abs(genCandPdgId[q])!=21) )continue;
	      //  std::cout<<"genCand: "<<q<<" Pt: "<<genCandPt[q]<<" Eta: "<<genCandEta[q]<<" Phi: "<<genCandPhi[q]<<" Mass: "<<genCandMass[q]<<std::endl;
	      TLorentzVector* trkP4 = new TLorentzVector();
	      trkP4->SetPtEtaPhiM(trkJetPt[pp],trkJetEta[pp],trkJetPhi[pp],trkJetMass[pp]);
	      TLorentzVector* quarkP4 = new TLorentzVector();
	      quarkP4->SetPtEtaPhiM(genCandPt[q],genCandEta[q],genCandPhi[q],genCandMass[q]);
	      dRJet = trkP4->DeltaR(*quarkP4);
	      if (dRJet<dRMinGen){
		dRMinGen=dRJet;
		idQuarkGen = q;
	      }
	    }
	    
	  }
	  //	    std::cout<<"------> "<<std::cout<<" jet: "<<pp<<" vbfq: "<<idQuarkGen<<" dr  : "<<dRMinGen<<std::endl;
	  //std::cout<<"------> "<<std::cout<<" jet: "<<pp<<" genq: "<<idQuark<<" dr  : "<<dRMin<<std::endl;
	  //std::cout<<" jet: "<<pp<<" q: "<<idQuark<<" dr: "<<dRMin<<std::endl;
	  bool matchToGenJet=true;
	  if(matchToGenJet==true){

	    for(int j = 0; j< ngenJet;j++){
	      double dRJet = 999.;
	      if(genJetPt[j]<10)continue;
	      //  std::cout<<"genCand: "<<q<<" Pt: "<<genCandPt[q]<<" Eta: "<<genCandEta[q]<<" Phi: "<<genCandPhi[q]<<" Mass: "<<genCandMass[q]<<std::endl;
	      TLorentzVector* trkP4 = new TLorentzVector();
	      trkP4->SetPtEtaPhiM(trkJetPt[pp],trkJetEta[pp],trkJetPhi[pp],trkJetMass[pp]);
	      TLorentzVector* genP4 = new TLorentzVector();
	      genP4->SetPtEtaPhiM(genJetPt[j],genJetEta[j],genJetPhi[j],genJetMass[j]);
	      dRJet = trkP4->DeltaR(*genP4);
	      if (dRJet<dRMinGenJet){
		dRMinGenJet=dRJet;
		idQuarkGenJet = j;
	      }
	    }


	  }


	  trkJetDR[pp]=dRMinGenJet;
	  vbfQuarkAssocId[pp]=idQuarkGenJet;
	  //	    if(idQuark>=0&&dRMin<0.3)std::cout<<" jet: "<<pp<<" q: "<<idQuark<<" dr: "<<dRMin<<std::endl;
	  h2_pt_eta->Fill(fabs(trkJetEta[pp]),trkJetPt[pp]);
	  h2_phi_eta->Fill(fabs(trkJetEta[pp]),trkJetPhi[pp]);
	  if(idQuark>=0&&dRMin<0.3)h2_resp_eta->Fill(fabs(trkJetEta[pp]),trkJetPt[pp]/vbfQuarkPt[idQuark]);
	  
	}
	
	//VBF ggh separation
	//	if(nMu<4)continue;
	
	//1. consider only jets that do not overlap with 4muons in the event
	std::vector<int> jetsOk;
	int nJetsOk = 0;
	for(int jet = 0; jet< ntrkJet;jet++){
	  bool isJetOk=true;
	  TLorentzVector* tjet = new TLorentzVector();
	  tjet->SetPtEtaPhiM(trkJetPt[jet], trkJetEta[jet], trkJetPhi[jet], trkJetMass[jet]);	  
	  for(int mu = 0;mu<nMu;mu++){
	    //  TLorentzVector* tmu = new TLorentzVector();
	    //tmu->SetPtEtaPhiM(MuPt[mu], MuEta[mu], MuPhi[mu], MuMass[mu]);
	    //double deltaR = tjet->DeltaR(*tmu);
	    double deltaRNew = sqrt(pow(MuEta[mu] - trkJetEta[jet],2) + pow(MuPhi[mu] - trkJetPhi[jet],2) );
	    if(deltaRNew<0.5) isJetOk = false;
	    }
	
	if(isJetOk==true&&trkJetPt[jet]>10.&&trkJetDR[jet]<0.1){
	  //  std::cout<<trkJetDR[jet]<<std::endl;
	  jetsOk.push_back(jet);
	  nJetsOk++;
	}
      }

      h_nJetsNoOverl->Fill(nJetsOk);
      

      //study the 2 jet category for the time being
      if(nJetsOk>1)   {  
	//	std::cout<<nJetsOk<<std::endl;     
	//look for 1st jet in the event
	double jet1Pt = 0;
	int jet1Id = 999;
	for(int jet = 0; jet< ntrkJet;jet++){
	  bool isGood = false;
	  for(int i = 0; i< jetsOk.size();i++) if(jet == jetsOk[i]) isGood = true;
	  if(!isGood)continue;
	  if(trkJetPt[jet]> jet1Pt){
	    jet1Pt = trkJetPt[jet];
	    jet1Id = jet;
	  }
	}
	//std::cout<<jet1Id<<std::endl;
	//look for 2nd jet in the event
	double jet2Pt = 0;
	int jet2Id = 999;
	for(int jet = 0; jet< ntrkJet;jet++){
	  if(jet == jet1Id)continue;
	  bool isGood = false;
	  for(int i = 0; i< jetsOk.size();i++) if(jet == jetsOk[i]) isGood = true;
	  if(!isGood)continue;
	  if(trkJetPt[jet]> jet2Pt && trkJetPt[jet]<jet1Pt){
	    jet2Pt = trkJetPt[jet];
	    jet2Id = jet;
	  }
	}
	//std::cout<<jet2Id<<std::endl;
	//save dijet system infos
	h_Jet1Pt->Fill(trkJetPt[jet1Id]);
	h_Jet1Eta->Fill(trkJetEta[jet1Id]);
	h_Jet1Phi->Fill(trkJetPhi[jet1Id]);
	
	h_Jet2Pt->Fill(trkJetPt[jet2Id]);
	h_Jet2Eta->Fill(trkJetEta[jet2Id]);
	h_Jet2Phi->Fill(trkJetPhi[jet2Id]);
	//	std::cout<<jet1Pt<<" "<<jet2Pt<<" "<<vbfQuarkAssocId[jet1Id]<< " "<<vbfQuarkAssocId[jet2Id]<<std::endl;
	if(jet1Pt>0 && jet2Pt>0&&vbfQuarkAssocId[jet1Id]!=vbfQuarkAssocId[jet2Id]){
	 

	  TLorentzVector* tjet1 = new TLorentzVector();
	  tjet1->SetPtEtaPhiM(trkJetPt[jet1Id], trkJetEta[jet1Id], trkJetPhi[jet1Id], trkJetMass[jet1Id]);	  
	  TLorentzVector* tjet2 = new TLorentzVector();
	  tjet2->SetPtEtaPhiM(trkJetPt[jet2Id], trkJetEta[jet2Id], trkJetPhi[jet2Id], trkJetMass[jet2Id]);	  

	  double deltaR12 = tjet1->DeltaR(*tjet2);
	  if(deltaR12<1)continue;
	  h_DeltaEtajj->Fill(fabs(trkJetEta[jet1Id]-trkJetEta[jet2Id]));
	  TLorentzVector sum;
	  sum = *tjet1+*tjet2;
	  double mjj = sum.M();
	  h_Mjj->Fill(mjj);
	  M_=mjj;
	  deltaEta_=fabs(trkJetEta[jet1Id]-trkJetEta[jet2Id]);
	  if(isDelphes=="true") h_Djj->Fill(0.206*deltaEta_+0.00011*M_-0.5799);
	  if(isDelphes=="false") h_Djj->Fill(0.1458*deltaEta_+0.00103*M_-0.497);
	  treeout->Fill();
	  //	  std::cout<<"dEta: "<<deltaEta_<<" trkJetEta[jet1Id]: "<<trkJetEta[jet1Id]<<" trkJetEta[jet2Id]: "<<trkJetEta[jet2Id]<<" trkJetPhi[jet1Id]: "<<trkJetPhi[jet1Id]<<" trkJetPhi[jet2Id]: "<<trkJetPhi[jet2Id]<<" q1: "<< vbfQuarkAssocId[jet1Id]<<" q2: "<< vbfQuarkAssocId[jet2Id]<<" dr: "<<deltaR12<<std::endl;
	}
	
      }


      
      
      //compute efficiency
      
      for(int quark = 0; quark< nvbfQuark;quark++){
	if(vbfQuarkPt[quark]<30) continue;
	double dRMinEff = 999.;
	h_den->Fill(vbfQuarkEta[quark]);
	for(int pp = 0; pp< ntrkJet;pp++){
	  double dRJet = 999.;	
	  if(trkJetPt[pp]<30)continue;
	  TLorentzVector* trkP4 = new TLorentzVector();
	  trkP4->SetPtEtaPhiM(trkJetPt[pp],trkJetEta[pp],trkJetPhi[pp],trkJetMass[pp]);
	  TLorentzVector* quarkP4 = new TLorentzVector();
	  quarkP4->SetPtEtaPhiM(vbfQuarkPt[quark],vbfQuarkEta[quark],vbfQuarkPhi[quark],vbfQuarkMass[quark]);
	  dRJet = trkP4->DeltaR(*quarkP4); 
	  if(dRJet<dRMinEff) dRMinEff = dRJet;	  
	}
	if(dRMinEff<0.3)h_num->Fill(vbfQuarkEta[quark]);	
      }
      
      
    
      if(abs(ntrk)<100000){
	int ntrkGood=0;
	for(int i = 0;i<ntrk;i++){
	  if(trkPt[i]<1)continue;
	  h_trkPt->Fill(trkPt[i]);
	  h_trkEta->Fill(trkEta[i]);
	  h_trkPhi->Fill(trkPhi[i]);
	  ntrkGood++;
	}
	h_ntrk->Fill(ntrkGood);
	for(int pp = 0; pp< ntrkJet;pp++){
	  int ntrkInJet=0;
	  if(trkJetPt[pp]<15)continue;
	  for(int i = 0;i<ntrk;i++){
	    if(trkPt[i]<1)continue;
	    double deltaR = sqrt(pow(trkJetEta[pp]-trkEta[i],2)+pow(trkJetPhi[pp]-trkPhi[i],2));       
	    if(deltaR<0.5){
	      h_trkInJetPt->Fill(trkPt[i]);
	      h_trkInJetPtEta->Fill(trkPt[i],trkEta[i]);
	      h_trkInJetP->Fill(trkPt[i]*TMath::CosH(trkEta[i]));
	      h_trkInJetEta->Fill(trkEta[i]);
	      h_trkInJetPhi->Fill(trkPhi[i]);
	      ntrkInJet++;
	    }	 
	  }
	h_ntrkInJet->Fill(ntrkInJet);
	}
      

      }
      if(nvbfQuark>0){
      h_nvbfQuark->Fill(nvbfQuark);
      for(int i = 0;i<nvbfQuark;i++){
	if(vbfQuarkPt[i]<30) continue;
	h_vbfQuarkPt->Fill(vbfQuarkPt[i]);
	h_vbfQuarkEta->Fill(vbfQuarkEta[i]);
	h_vbfQuarkPhi->Fill(vbfQuarkPhi[i]);
	h_vbfQuarkMass->Fill(vbfQuarkMass[i]);
      }
      }
	 
   }
   std::cout<<"out"<<std::endl;
   
   gStyle->SetOptStat(0);
   TCanvas* c1 = new TCanvas("c1", "c1", 1);
   c1->cd();
   
   TLatex* latex = new TLatex(1.01, 0.96,"CMS Phase II Simulation          #sqrt{s}=14 TeV, PU=140 ");
   latex->SetNDC();
   latex->SetTextAngle(0);
   latex->SetTextColor(kBlack);  
   latex->SetTextFont(42);
   latex->SetTextAlign(31); 
   latex->SetTextSize(0.05);    
 
  
   std::cout<<counter<<std::endl;
   //trkInJet pt
   h_trkInJetPt->SetLineColor(kBlue);
   h_trkInJetPt->GetYaxis()->SetTitle("Events");     
   h_trkInJetPt->GetXaxis()->SetTitle("Tracks p_{T} [GeV]");     
   
   h_trkInJetPt->DrawNormalized("hist");
   latex->Draw("same"); 
   c1->SetLogy(0); 
   c1->SaveAs("plots/Tracks_Pt.png");
   c1->SetLogy();
   c1->SaveAs("plots/Tracks_Pt_LOG.png");
  
   //trkInJet pt eta
   h_trkInJetPtEta->SetLineColor(kBlue);
   h_trkInJetPtEta->GetYaxis()->SetTitle("Tracks #eta");     
   h_trkInJetPtEta->GetXaxis()->SetTitle("Tracks p_{T} [GeV]");     
   
   h_trkInJetPtEta->Draw("COLZ");
   latex->Draw("same"); 
   c1->SetLogy(0); 
   c1->SaveAs("plots/Tracks_PtEta.png");
   
   //trkInJet p
   h_trkInJetP->SetLineColor(kBlue);
   h_trkInJetP->GetYaxis()->SetTitle("Events");     
   h_trkInJetP->GetXaxis()->SetTitle("Tracks |p| [GeV]");     
   
   h_trkInJetP->DrawNormalized("hist");
   latex->Draw("same"); 
   c1->SetLogy(0); 
   c1->SaveAs("plots/Tracks_P.png");
   c1->SetLogy();
   c1->SaveAs("plots/Tracks_P_LOG.png");
   
   //trkInJet eta
   h_trkInJetEta->SetLineColor(kBlue);
   h_trkInJetEta->GetYaxis()->SetTitle("Events");     
   h_trkInJetEta->GetXaxis()->SetTitle("Tracks #eta");     
  
   h_trkInJetEta->DrawNormalized("hist"); 
   latex->Draw("same"); 
   c1->SetLogy(0); 
   c1->SaveAs("plots/Tracks_Eta.png");
   c1->SetLogy();
   c1->SaveAs("plots/Tracks_Eta_LOG.png");
   
   //trkInJet phi
   h_trkInJetPhi->SetLineColor(kBlue);
   h_trkInJetPhi->GetYaxis()->SetTitle("Events");     
   h_trkInJetPhi->GetXaxis()->SetTitle("Tracks #phi");     
 
   h_trkInJetPhi->DrawNormalized("hist"); 
   latex->Draw("same"); 
   c1->SetLogy(0);   
   c1->SaveAs("plots/Tracks_Phi.png");
   c1->SetLogy();
   c1->SaveAs("plots/Tracks_Phi_LOG.png");
   


   /*
   //jet n
   h_ntrkJet->SetLineColor(kBlue);
   h_nvbfQuark->SetLineColor(kMagenta+7);
   h_nvbfQuark->SetLineWidth(2);
   h_ntrkJet->SetLineWidth(2);

   TLegend* legmc;
   legmc = new TLegend(0.65, 0.75, 0.95, 0.95, "", "brNDC"); 
   legmc->SetTextFont(42);
   legmc->SetBorderSize(0);
   legmc->SetFillStyle(0);
   legmc->AddEntry(h_ntrkJet,"Tracker Jets", "L");
   legmc->AddEntry(h_nvbfQuark,"VBF Quarks", "L");

   h_nvbfQuark->GetYaxis()->SetTitle("Events");     
   h_nvbfQuark->GetXaxis()->SetTitle("Number of Jets");     
 
   h_nvbfQuark->DrawNormalized("hist");
   h_ntrkJet->DrawNormalized("histsame"); 
   latex->Draw("same"); 
   legmc->Draw("same");
   c1->SetLogy(0); 
   c1->SaveAs("plots/JetN.png");
   c1->SetLogy();
   c1->SaveAs("plots/JetN_LOG.png");
  
   //jet pt
   h_trkJetPt->SetLineColor(kBlue);
   h_vbfQuarkPt->SetLineColor(kMagenta+7);
   h_vbfQuarkPt->SetLineWidth(2);
   h_trkJetPt->SetLineWidth(2);
   h_trkJetPt->GetYaxis()->SetTitle("Events");     
   h_trkJetPt->GetXaxis()->SetTitle("Jet p_{T} [GeV]");     
 
   h_trkJetPt->DrawNormalized("hist"); 
   h_vbfQuarkPt->DrawNormalized("histsame");
   latex->Draw("same"); 
   legmc->Draw("same");
   c1->SetLogy(0); 
   c1->SaveAs("plots/JetPt.png");
   c1->SetLogy();
   c1->SaveAs("plots/JetPt_LOG.png");
   
   //jet eta
   h_trkJetEta->SetLineColor(kBlue);
   h_vbfQuarkEta->SetLineColor(kMagenta+7);
   h_vbfQuarkEta->SetLineWidth(2);
			    h_trkJetEta->SetLineWidth(2);
   h_trkJetEta->GetYaxis()->SetTitle("Events");     
   h_trkJetEta->GetXaxis()->SetTitle("Jet #eta ");     
 
   h_trkJetEta->DrawNormalized("hist"); 
   h_vbfQuarkEta->DrawNormalized("histsame");
   latex->Draw("same"); 
   legmc->Draw("same");
   c1->SetLogy(0); 
   c1->SaveAs("plots/JetEta.png");
   c1->SetLogy();
   c1->SaveAs("plots/JetEta_LOG.png");
 
  //jet phi
   h_trkJetPhi->SetLineColor(kBlue);
   h_vbfQuarkPhi->SetLineColor(kMagenta+7);
   h_vbfQuarkPhi->SetLineWidth(2);
   h_trkJetPhi->SetLineWidth(2);
   h_trkJetPhi->GetYaxis()->SetTitle("Events");     
   h_trkJetPhi->GetXaxis()->SetTitle("Jet #phi ");     
 
   h_trkJetPhi->DrawNormalized("hist"); 
   h_vbfQuarkPhi->DrawNormalized("histsame");
   latex->Draw("same"); 
   legmc->Draw("same");
   c1->SetLogy(0); 
   c1->SaveAs("plots/JetPhi.png");
   c1->SetLogy();
   c1->SaveAs("plots/JetPhi_LOG.png");
 

   //jet mass
   h_trkJetMass->SetLineColor(kBlue);
   h_vbfQuarkMass->SetLineColor(kMagenta+7);
   h_vbfQuarkMass->SetLineWidth(2);
   h_trkJetMass->SetLineWidth(2);
   h_trkJetMass->GetYaxis()->SetTitle("Events");     
   h_trkJetMass->GetXaxis()->SetTitle("Jet Mass [GeV]");     
 
   h_trkJetMass->DrawNormalized("hist"); 
   h_vbfQuarkMass->DrawNormalized("histsame");
   latex->Draw("same"); 
   legmc->Draw("same");
   c1->SetLogy(0); 
   c1->SaveAs("plots/JetMass.png");
   c1->SetLogy();
   c1->SaveAs("plots/JetMass_LOG.png");
 

   //jet DR
   h_trkJetDR->SetLineColor(kBlue);
   h_trkJetDR->SetLineWidth(2);
   h_trkJetDR->GetYaxis()->SetTitle("Events");     
   h_trkJetDR->GetXaxis()->SetTitle("#Delta R Matching"); 
 
   h_trkJetDR->DrawNormalized("hist");
   latex->Draw("same"); 
   c1->SetLogy(0);     
   c1->SaveAs("plots/JetDR.png");
   c1->SetLogy();
   c1->SaveAs("plots/JetDR_LOG.png");
   
   
   //2D distributions
   c1->SetLogy(0);
  
   h2_pt_eta->GetYaxis()->SetTitle("Tracker Jet p_{T}");     
   h2_pt_eta->GetXaxis()->SetTitle("Tracker Jet |#eta|");    
   h2_pt_eta->Draw("COLZ");
   latex->Draw("same"); 
   c1->SaveAs("plots/Jet_pt_vs_eta.png");
   
  
   h2_phi_eta->GetYaxis()->SetTitle("Tracker Jet #phi");     
   h2_phi_eta->GetXaxis()->SetTitle("Tracker Jet |#eta|");     
   h2_phi_eta->Draw("COLZ");
   latex->Draw("same"); 
   c1->SaveAs("plots/Jet_phi_vs_eta.png");
   
   
   h2_resp_eta->GetYaxis()->SetTitle("Trk Jet p_{T}/Gen Jet p_{T}  ");     
   h2_resp_eta->GetXaxis()->SetTitle("Trk Jet |#eta|");    
   h2_resp_eta->Draw("COLZ");
   latex->Draw("same"); 
   c1->SaveAs("plots/Jet_response_vs_eta.png");
  
   //means and rms response vs eta
   int etabins=10;
   double mean[etabins];
   double rms[etabins];
   double mean_err[etabins];
   double rms_err[etabins];
 
   for(int k=1; k<=10;k++){
     std::cout<<k<<std::endl;
     TH1F* h1 = (TH1F*) h2_resp_eta->ProjectionY("",k,k,"e");
     h1->GetXaxis()->SetTitle("Trk Jet p_{T} /Gen Jet p_{T}");
     h1->GetYaxis()->SetTitle("Entries");
     bool doFit=false;
     if(doFit){
       if(h1->GetEntries()>0) h1->Fit("gaus", "", "",  0., 2.);
       h1->GetFunction("gaus")->SetLineColor(kBlue);
       h1->GetFunction("gaus")->SetLineWidth(2);       
       mean[k] = h1->GetFunction("gaus")->GetParameter(1);
       mean_err[k] = h1->GetFunction("gaus")->GetParError(1);
       rms[k]  = h1->GetFunction("gaus")->GetParameter(2);
       rms_err[k] = h1->GetFunction("gaus")->GetParError(2);
     }else{

       mean[k-1] = h1->GetMean();
       mean_err[k-1] = h1->GetMeanError();
       rms[k-1]  = h1->GetRMS();
       rms_err[k-1] = h1->GetRMSError();
       std::cout<<mean[k-1]<<std::endl;
     }
     h1->Draw("P");

     c1->SaveAs(TString::Format("plots/Jet_response_fit_bin_%d.png",k));
     }
   double eta[10] = {0.2, 0.6, 1,1.4, 1.8, 2.2, 2.6, 3.0, 3.4, 3.8};
   double eta_err[10] = {0.2, 0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2};
   TGraphErrors* means = new TGraphErrors(10, eta, mean,eta_err, mean_err);
   TGraphErrors* RMSs = new TGraphErrors(10, eta,  rms,eta_err, rms_err);
   means->SetTitle("");
   means->SetName("meanResponse");
   means->Draw("APE");
   means->GetYaxis()->SetRangeUser(0., 1.2);
   means->GetYaxis()->SetTitle("<Trk Jet p_{T} /VBF Quark p_{T}>");
   latex->Draw("same"); 
   c1->SaveAs("plots/RESP_mean_vs_eta.png");
   RMSs->Draw("APE");
   RMSs->SetName("rmsResponse");

   latex->Draw("same"); 
   c1->SaveAs("plots/RESP_RMS_vs_eta.png");


   //efficiency  
   TEfficiency* pEff =new TEfficiency(*h_num,*h_den);
   TLatex* latexY = new TLatex(0.06, 0.95,"Tracker Efficiency");
   latexY->SetNDC();
   latexY->SetTextAngle(90);
   latexY->SetTextColor(kBlack);  
   latexY->SetTextFont(42);
   latexY->SetTextAlign(31); 
   latexY->SetTextSize(0.05);    
   TLatex* latexX = new TLatex(0.95, 0.05,"VBF Quark |#eta|");
   latexX->SetNDC();
   latexX->SetTextAngle(0);
   latexX->SetTextColor(kBlack);  
   latexX->SetTextFont(42);
   latexX->SetTextAlign(31); 
   latexX->SetTextSize(0.05);    
  
 
   pEff->Draw("AP");
   pEff->SetName("efficiency");
   latex->Draw("same"); 
   latexY->Draw("same"); 
   latexX->Draw("same"); 
   c1->SaveAs("plots/eff_vs_eta.png");
   */
   f_out->cd();

   h_nJetsNoOverl->Write();

   h_DeltaEtajj->Write();

   h_Mjj->Write();
   h_Djj->Write();

   h_trkJetPt->Write();
   h_trkJetEta->Write();
   h_trkJetPhi->Write();
  
   h_genJetPt->Write();
   h_genJetEta->Write();
   h_genJetPhi->Write();
  
   h_nTrkJet->Write();
   h_nGenJet->Write();
   
   h_Jet1Pt->Write();
   h_Jet1Eta->Write();
   h_Jet1Phi->Write();

   h_Jet2Pt->Write();
   h_Jet2Eta->Write();
   h_Jet2Phi->Write();
  
   h_trkPt->Write();
   h_trkPhi->Write();
   h_trkEta->Write();
   h_ntrk->Write();
  
   h_vbfQuarkPt->Write();
   h_vbfQuarkPhi->Write();
   h_vbfQuarkEta->Write();
   h_nvbfQuark->Write();
  
   h_trkInJetPt->Write();
   h_trkInJetPhi->Write();
   h_trkInJetEta->Write();
   h_ntrkInJet->Write();
  
   // pEff->Write();
   // means->Write();
   //RMSs->Write();
   treeout->Fill(); 
   f_out->Write();
   f_out->Close();


}
