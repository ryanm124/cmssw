#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

void SetPlotStyle();
void mySmallText(Double_t x,Double_t y,Color_t color,char *text); 


// ----------------------------------------------------------------------------------------------------------------
// Main script
void overlay(TString what) {

  SetPlotStyle();

  TFile* tree1 = new TFile("output_Combined_D49_DispMu_Dev_NoTrunc_MaxR10.root");
  TFile* tree2 = new TFile("output_Combined_D49_DispMu_CMSSW_NoTrunc_MaxR10.root");
  
  TH1F* h1 = (TH1F*) tree1->Get(what);
  TH1F* h2 = (TH1F*) tree2->Get(what);

  if (what.Contains("ntrk")) {
    h1->Rebin(10);
    h2->Rebin(10);
  }

  TCanvas c;
  h1->SetLineColor(1);
  h1->SetMarkerColor(1);
  h1->SetMarkerStyle(8);
  h1->Draw("p");
  h2->SetLineColor(2);
  h2->SetMarkerColor(2);
  h2->SetMarkerStyle(24);
  h2->Draw("p,same");
  
  TLegend* l;
  if (what.Contains("eff")) {
    l = new TLegend(0.55,0.22,0.85,0.40);
    mySmallText(0.25,0.22,1,"DisplacedMu, PU=0");
  }
  else {
    l = new TLegend(0.2,0.72,0.5,0.9);
    mySmallText(0.6,0.82,1,"DisplacedMu, PU=0");
  }
  l->SetFillColor(0);
  l->SetLineColor(0);
  l->SetTextSize(0.04);
  l->AddEntry(h1,"Dev","lep");
  l->AddEntry(h2,"CMSSW","lep");
  l->SetTextFont(42);
  l->Draw();	

  gPad->SetGridy();

  c.SaveAs("Overlay/D49_DispMu_Dev_vs_CMSSW_NoApproxBoth_NoTrunc_"+what+".pdf");


}


void SetPlotStyle() {

  // from ATLAS plot style macro
  gStyle->SetErrorX(0.);
  // use plain black on white colors
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetHistLineColor(1);

  gStyle->SetPalette(1);

  // set the paper & margin sizes
  gStyle->SetPaperSize(20,26);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  gStyle->SetTitleXOffset(1.4);
  gStyle->SetTitleYOffset(1.4);

  // use large fonts
  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.05);
  gStyle->SetLabelFont(42,"x");
  gStyle->SetTitleFont(42,"x");
  gStyle->SetLabelFont(42,"y");
  gStyle->SetTitleFont(42,"y");
  gStyle->SetLabelFont(42,"z");
  gStyle->SetTitleFont(42,"z");
  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetTitleSize(0.05,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetTitleSize(0.05,"y");
  gStyle->SetLabelSize(0.05,"z");
  gStyle->SetTitleSize(0.05,"z");

  // use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetHistLineWidth(2.);
  gStyle->SetLineStyleString(2,"[12 12]");

  // get rid of error bar caps
  gStyle->SetEndErrorSize(0);

  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

}


void mySmallText(Double_t x,Double_t y,Color_t color,char *text) {
  Double_t tsize=0.044;
  TLatex l;
  l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}


