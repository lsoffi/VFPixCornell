/*#include <iostream> 
#include <algorithm>
#include <vector>
#include <TStyle.h>
#include <fstream>
#include <TROOT.h>
#include <TObjArray.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraphSmooth.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include "TGraph2D.h"
#include <TH1F.h>
#include <TString.h>
#include <TCut.h>
#include <TChain.h>
#include <TF1.h>
#include "TProfile.h"
#include <TPaveText.h>
#include <TLegend.h>
#include <TFile.h>
#include "TH2F.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooBreitWigner.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooHistFunc.h"
#include "TLatex.h"
#include <TPaletteAxis.h>
#include <RooWorkspace.h>
#include <RooAbsPdf.h>
#include <RooGenericPdf.h>
#include <RooCmdArg.h>
#include <RooGaussian.h>
using namespace RooFit;
*/

void compareFisher(){
 TCanvas* c1 = new TCanvas("c1", "c1", 1);
 c1->cd();
 TLatex* latex = new TLatex(0.8534483,0.965035,"CMS Phase II Simulation                             #sqrt{s}=14 TeV, PU=140 ");
 latex->SetNDC();
 latex->SetTextAngle(0);
 latex->SetTextColor(kBlack);  
 latex->SetTextFont(42);
 latex->SetTextAlign(31); 
 latex->SetTextSize(0.03);    
 TFile* f_d = TFile::Open("fisher_true.root");
 TFile* f_f = TFile::Open("fisher_false.root");
 TH1F* h_d = (TH1F*)f_d->Get("Method_Fisher/Fisher/MVA_Fisher_rejBvsS");
 TH1F* h_f = (TH1F*)f_f->Get("Method_Fisher/Fisher/MVA_Fisher_rejBvsS");
 h_d->SetLineColor(kBlue);
 h_f->SetLineColor(kBlack);
 h_d->Draw("hist");
 h_f->Draw("histsame");

 TLegend* legmc;
 legmc = new TLegend(0.302299,0.4377622,0.5405172,0.6388112, "", "brNDC"); 
 legmc->SetTextFont(42);
 legmc->SetBorderSize(0);
 legmc->SetFillStyle(0);
 legmc->AddEntry(h_d, "Delphes", "L");
 legmc->AddEntry(h_f, "Full-Sim", "L");
 latex->Draw("same");
 legmc->Draw("same");

 c1->SetLogy(0);
 c1->SaveAs("plots/Fisher_Comaprison.png");
 c1->SaveAs("plots/Fisher_Comaprison.pdf");


}

void PlotEfficiency() {


 TCanvas* c1 = new TCanvas("c1", "c1", 1);
 c1->cd();
 TLatex* latex = new TLatex(0.8534483,0.965035,"CMS Phase II Simulation                             #sqrt{s}=14 TeV, PU=140 ");
 latex->SetNDC();
 latex->SetTextAngle(0);
 latex->SetTextColor(kBlack);  
 latex->SetTextFont(42);
 latex->SetTextAlign(31); 
 latex->SetTextSize(0.03);    
  
 TFile* f_noFit = TFile::Open("output_noFit_.root");
 TEfficiency* eff_noFit= (TEfficiency*) f_noFit->Get("efficiency");
 TFile* f_noFit02mm = TFile::Open("output_noFit02mm_.root");
 TEfficiency* eff_noFit02mm= (TEfficiency*) f_noFit02mm->Get("efficiency");
 TFile* f_yesFit = TFile::Open("output_yesFit_.root");
 TEfficiency* eff_yesFit= (TEfficiency*) f_yesFit->Get("efficiency");
 TLatex* latexY = new TLatex(0.06, 0.95,"Tracker Efficiency");
 latexY->SetNDC();
 latexY->SetTextAngle(90);
 latexY->SetTextColor(kBlack);  
 latexY->SetTextFont(42);
 latexY->SetTextAlign(31); 
 latexY->SetTextSize(0.05);    
 TLatex* latexX = new TLatex(0.8491379,0.04020979,"VBF Quark |#eta|");
 latexX->SetNDC();
 latexX->SetTextAngle(0);
 latexX->SetTextColor(kBlack);  
 latexX->SetTextFont(42);
 latexX->SetTextAlign(31); 
 latexX->SetTextSize(0.05);  
 
 TLegend* legmc;
 legmc = new TLegend(0.502299,0.6377622,0.8405172,0.9388112, "", "brNDC"); 
 legmc->SetTextFont(42);
 legmc->SetBorderSize(0);
 legmc->SetFillStyle(0);
 legmc->AddEntry(eff_noFit, "All tracks with #Delta Z < 1 cm", "PL");
 legmc->AddEntry(eff_noFit02mm, "All tracks with #Delta Z < 0.2 cm", "PL");
 legmc->AddEntry(eff_yesFit, "All tracks from fit to vertices", "PL");
 
 eff_noFit->SetLineColor(kRed);
 eff_noFit->SetMarkerColor(kRed);
 eff_noFit02mm->SetLineColor(kBlue);
 eff_noFit02mm->SetMarkerColor(kBlue);
 eff_yesFit->SetLineColor(kBlack);
 eff_yesFit->SetMarkerColor(kBlack);

 eff_noFit->Draw("AP");
 eff_noFit02mm->Draw("Psame");
 eff_yesFit->Draw("Psame");
 legmc->Draw("same");
 latex->Draw("same"); 
 latexY->Draw("same"); 
 latexX->Draw("same"); 

 c1->SetLogy(0);
 c1->SaveAs("~/www/upgrade/VFPix-trkJets/recoEff_comparison.png");
 c1->SaveAs("~/www/upgrade/VFPix-trkJets/recoEff_comparison.pdf");
 }





void PlotResponse() {


 TCanvas* c1 = new TCanvas("c1", "c1", 1);
 c1->cd();
 TLatex* latex = new TLatex(0.8534483,0.965035,"CMS Phase II Simulation                             #sqrt{s}=14 TeV, PU=140 ");
 latex->SetNDC();
 latex->SetTextAngle(0);
 latex->SetTextColor(kBlack);  
 latex->SetTextFont(42);
 latex->SetTextAlign(31); 
 latex->SetTextSize(0.03);    
  
 TFile* f_noFit = TFile::Open("output_noFit_.root");
 TGraphErrors* resp_noFit= (TGraphErrors*) f_noFit->Get("meanResponse");
 TFile* f_noFit02mm = TFile::Open("output_noFit02mm_.root");
 TGraphErrors* resp_noFit02mm= (TGraphErrors*) f_noFit02mm->Get("meanResponse");
 TFile* f_yesFit = TFile::Open("output_yesFit_.root");
 TGraphErrors* resp_yesFit= (TGraphErrors*) f_yesFit->Get("meanResponse");
 TLatex* latexY = new TLatex(0.06, 0.95,"Tracker Response");
 latexY->SetNDC();
 latexY->SetTextAngle(90);
 latexY->SetTextColor(kBlack);  
 latexY->SetTextFont(42);
 latexY->SetTextAlign(31); 
 latexY->SetTextSize(0.05);    
 TLatex* latexX = new TLatex(0.8491379,0.04020979,"Trk Jet |#eta|");
 latexX->SetNDC();
 latexX->SetTextAngle(0);
 latexX->SetTextColor(kBlack);  
 latexX->SetTextFont(42);
 latexX->SetTextAlign(31); 
 latexX->SetTextSize(0.05);  
 
 TLegend* legmc;
 legmc = new TLegend(0.402299,0.577622,0.8405172,0.9388112, "", "brNDC"); 
 legmc->SetTextFont(42);
 legmc->SetBorderSize(0);
 legmc->SetFillStyle(0);
 legmc->AddEntry(resp_noFit, "All tracks with #Delta Z < 1 cm", "PL");
 legmc->AddEntry(resp_noFit02mm, "All tracks with #Delta Z < 0.2 cm", "PL");
 legmc->AddEntry(resp_yesFit, "All tracks from fit to vertices", "PL");
 
 resp_noFit->SetLineColor(kRed);
 resp_noFit->SetMarkerColor(kRed);
 resp_noFit02mm->SetLineColor(kBlue);
 resp_noFit02mm->SetMarkerColor(kBlue);
 resp_yesFit->SetLineColor(kBlack);
 resp_yesFit->SetMarkerColor(kBlack);
 resp_noFit->Draw("AP");
 resp_noFit02mm->Draw("Psame");
 resp_yesFit->Draw("Psame");
 legmc->Draw("same");
 latex->Draw("same"); 
 // latexY->Draw("same"); 
 latexX->Draw("same"); 

 c1->SetLogy(0);
 c1->SaveAs("~/www/upgrade/VFPix-trkJets/recoResp_comparison.png");
 c1->SaveAs("~/www/upgrade/VFPix-trkJets/recoResp_comparison.pdf");
 }






void PlotGGHVBF(std::string var, std::string channel) {


 TCanvas* c1 = new TCanvas("c1", "c1", 1);
 c1->cd();
 TLatex* latex = new TLatex(0.8534483,0.965035,"CMS Phase II Simulation                             #sqrt{s}=14 TeV, PU=140 ");
 latex->SetNDC();
 latex->SetTextAngle(0);
 latex->SetTextColor(kBlack);  
 latex->SetTextFont(42);
 latex->SetTextAlign(31); 
 latex->SetTextSize(0.03);    
  
 TFile* f_VBF = TFile::Open(("output_VBF_"+channel+".root").c_str());
 TFile* f_GGH = TFile::Open(("output_GGH_"+channel+".root").c_str());
 TH1F* hvbf = (TH1F*)f_VBF->Get(("h_"+var).c_str());
 TH1F* hggh = (TH1F*)f_GGH->Get(("h_"+var).c_str());
 TLatex* latexY = new TLatex(0.06, 0.95,"a.u.");
 latexY->SetNDC();
 latexY->SetTextAngle(90);
 latexY->SetTextColor(kBlack);  
 latexY->SetTextFont(42);
 latexY->SetTextAlign(31); 
 latexY->SetTextSize(0.05);    
 TLatex* latexX = new TLatex(0.8491379,0.04020979,var.c_str());
 latexX->SetNDC();
 latexX->SetTextAngle(0);
 latexX->SetTextColor(kBlack);  
 latexX->SetTextFont(42);
 latexX->SetTextAlign(31); 
 latexX->SetTextSize(0.05);  
 
 TLegend* legmc;
 legmc = new TLegend(0.402299,0.577622,0.8405172,0.9388112, "", "brNDC"); 
 legmc->SetTextFont(42);
 legmc->SetBorderSize(0);
 legmc->SetFillStyle(0);
 legmc->AddEntry(hggh, "Gluon fusion H #rightarrow 4L", "PL");
 legmc->AddEntry(hvbf, "VBF H #rightarrow 4L", "PL");
 
 hggh->SetLineColor(kRed);
 hggh->SetMarkerColor(kRed);
 hvbf->SetLineColor(kBlue);
 hvbf->SetMarkerColor(kBlue);


 if(var=="nJetsNoOverl")hggh->GetXaxis()->SetTitle("# Non-Overlap. Trk Jets");
 if(var=="deltaEtajj") hggh->GetXaxis()->SetTitle("#Delta #eta_{jj}");
 if(var=="Mjj")hggh->GetXaxis()->SetTitle("M_{jj}");
 if(var=="Djj")hggh->GetXaxis()->SetTitle("D_{jj}");
 hggh->DrawNormalized("hist");
 hvbf->DrawNormalized("histsame");
 legmc->Draw("same");
 latex->Draw("same"); 
 // latexY->Draw("same"); 
 //latexX->Draw("same"); 

 c1->SetLogy(0);
 c1->SaveAs(("plots/GGHVBF_"+var+"_"+channel+".pdf").c_str());
 c1->SaveAs(("plots/GGHVBF_"+var+"_"+channel+".png").c_str());
 c1->SetLogy(1);
 c1->SaveAs(("plots/GGHVBF_"+var+"_"+channel+"_LOG.pdf").c_str());
 c1->SaveAs(("plots/GGHVBF_"+var+"_"+channel+"_LOG.png").c_str());

 }



void PlotAllDelphesFullSim(){

  PlotDelphesFullSim("nJetsNoOverl", "VBF");   
  PlotDelphesFullSim("Mjj", "VBF");	       
  PlotDelphesFullSim("deltaEtajj", "VBF");     
  PlotDelphesFullSim("Jet1Pt", "VBF"); 	       
  PlotDelphesFullSim("Jet1Eta", "VBF");	       
  PlotDelphesFullSim("Jet1Phi", "VBF");	       
  PlotDelphesFullSim("Jet2Pt", "VBF"); 	       
  PlotDelphesFullSim("Jet2Eta", "VBF");	       
  PlotDelphesFullSim("Jet2Phi", "VBF");        
  PlotDelphesFullSim("nGenJet", "VBF");        
  //  PlotDelphesFullSim("nMatchJet", "VBF");        
  PlotDelphesFullSim("nTrkJet", "VBF");        
  PlotDelphesFullSim("ntrk", "VBF");        
  PlotDelphesFullSim("trkPt", "VBF");        
  PlotDelphesFullSim("trkEta", "VBF");        
  PlotDelphesFullSim("trkPhi", "VBF");        
  PlotDelphesFullSim("ntrkInJet", "VBF");        
  PlotDelphesFullSim("trkInJetPt", "VBF");        
  PlotDelphesFullSim("trkInJetEta", "VBF");        
  PlotDelphesFullSim("trkInJetPhi", "VBF");        
  PlotDelphesFullSim("genJetPt", "VBF");        
  PlotDelphesFullSim("genJetEta", "VBF");        
  PlotDelphesFullSim("genJetPhi", "VBF");        
  
  PlotDelphesFullSim("nJetsNoOverl", "GGH");   
  PlotDelphesFullSim("Mjj", "GGH");	       
  PlotDelphesFullSim("deltaEtajj", "GGH");     
  PlotDelphesFullSim("Jet1Pt", "GGH"); 	       
  PlotDelphesFullSim("Jet1Eta", "GGH");	       
  PlotDelphesFullSim("Jet1Phi", "GGH");	       
  PlotDelphesFullSim("Jet2Pt", "GGH"); 	       
  PlotDelphesFullSim("Jet2Eta", "GGH");	       
  PlotDelphesFullSim("Jet2Phi", "GGH");        
  PlotDelphesFullSim("nGenJet", "GGH");        
  //  PlotDelphesFullSim("nMatchJet", "GGH");        
  PlotDelphesFullSim("nTrkJet", "GGH");        
  PlotDelphesFullSim("ntrk", "GGH");        
  PlotDelphesFullSim("trkPt", "GGH");        
  PlotDelphesFullSim("trkEta", "GGH");        
  PlotDelphesFullSim("trkPhi", "GGH");  
  PlotDelphesFullSim("ntrkInJet", "GGH");        
  PlotDelphesFullSim("trkInJetPt", "GGH");        
  PlotDelphesFullSim("trkInJetEta", "GGH");        
  PlotDelphesFullSim("trkInJetPhi", "GGH");         
  PlotDelphesFullSim("genJetPt", "GGH");        
  PlotDelphesFullSim("genJetEta", "GGH");        
  PlotDelphesFullSim("genJetPhi", "GGH");        





}

void PlotDelphesFullSim(std::string var, std::string channel) {


 TCanvas* c1 = new TCanvas("c1", "c1", 1);
 c1->cd();
 TLatex* latex = new TLatex(0.8534483,0.965035,"CMS Phase II Simulation                             #sqrt{s}=14 TeV, PU=140 ");
 latex->SetNDC();
 latex->SetTextAngle(0);
 latex->SetTextColor(kBlack);  
 latex->SetTextFont(42);
 latex->SetTextAlign(31); 
 latex->SetTextSize(0.03);    
  
 TFile* f_d = TFile::Open(("output_"+channel+"_true.root").c_str());
 TFile* f_f = TFile::Open(("output_"+channel+"_false.root").c_str());
 TH1F* hd = (TH1F*)f_d->Get(("h_"+var).c_str());
 TH1F* hf = (TH1F*)f_f->Get(("h_"+var).c_str());
 TLatex* latexY = new TLatex(0.06, 0.95,"a.u.");
 latexY->SetNDC();
 latexY->SetTextAngle(90);
 latexY->SetTextColor(kBlack);  
 latexY->SetTextFont(42);
 latexY->SetTextAlign(31); 
 latexY->SetTextSize(0.05);    
 TLatex* latexX = new TLatex(0.8491379,0.04020979,var.c_str());
 latexX->SetNDC();
 latexX->SetTextAngle(0);
 latexX->SetTextColor(kBlack);  
 latexX->SetTextFont(42);
 latexX->SetTextAlign(31); 
 latexX->SetTextSize(0.05);  
 
 TLegend* legmc;
 legmc = new TLegend(0.402299,0.577622,0.8405172,0.9388112, "", "brNDC"); 
 legmc->SetTextFont(42);
 legmc->SetBorderSize(0);
 legmc->SetFillStyle(0);
 if(channel=="VBF"){
   legmc->AddEntry(hd, "VBF Delphes", "L");
   legmc->AddEntry(hf, "VBF Fullsim", "L");
 }else{
   legmc->AddEntry(hd, "GGH Delphes", "L");
   legmc->AddEntry(hf, "GGH Fullsim", "L");
 }
 TColor col;
 if(channel=="VBF"){ 
   hf->SetLineColor(kBlue);
   hf->SetMarkerColor(kBlue);
   hd->SetLineColor(kBlue);
   hd->SetMarkerColor(kBlue);
 }else{ 
   hf->SetLineColor(kRed);
   hf->SetMarkerColor(kRed);
   hd->SetLineColor(kRed);
   hd->SetMarkerColor(kRed);
 }
 hd->SetLineStyle(kDashed);
 hf->Scale(1./hf->Integral());
 hd->Scale(1./hd->Integral());
 if(var=="Mjj"){
   hf->Rebin(4);
   hd->Rebin(4);
   hf->GetYaxis()->SetRangeUser(0.001,0.25);
 }
if(var=="deltaEtajj"){
   hf->GetYaxis()->SetRangeUser(0.001,0.25);
 }
 if(var=="nJetsNoOverl")hf->GetXaxis()->SetTitle("# Non-Overlap. Trk Jets");
 if(var=="deltaEtajj") hf->GetXaxis()->SetTitle("#Delta #eta_{jj}");
 if(var=="Mjj")hf->GetXaxis()->SetTitle("M_{jj}");
 if(var=="Djet")hf->GetXaxis()->SetTitle("D_{jj}");
 hd->Draw("hist");
 hf->Draw("histsame");
 legmc->Draw("same");
 latex->Draw("same"); 
 // latexY->Draw("same"); 
 //latexX->Draw("same"); 

 c1->SetLogy(0);
 c1->SaveAs(("plots/DELPHESvsFullSIM_"+var+"_"+channel+".pdf").c_str());
 c1->SaveAs(("plots/DELPHESvsFullSIM_"+var+"_"+channel+".png").c_str());
 c1->SetLogy(1);
 c1->SaveAs(("plots/DELPHESvsFullSIM_"+var+"_"+channel+"_LOG.pdf").c_str());
 c1->SaveAs(("plots/DELPHESvsFullSIM_"+var+"_"+channel+"_LOG.png").c_str());

 }




