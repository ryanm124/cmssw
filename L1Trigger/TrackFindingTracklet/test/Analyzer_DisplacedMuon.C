
// ----------------------------------------------------------------------------------------------------------------
// Feasibility study of using L1 Tracks to identify Displaced Vertex
//
// By Bharadwaj Harikrishnan, May 2021
// Edited by Ryan McCarthy, Sept 2021
// ----------------------------------------------------------------------------------------------------------------

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TVector3.h>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include <TCanvas.h>
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include <math.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <valarray>
#include <deque>
#include <THStack.h>
#include <TF2.h>
#include <TEllipse.h>
#include <TMarker.h>
#include <chrono>
#include <boost/variant.hpp>
using namespace std;

bool VERBOSE[]={false,false,false,false};
bool detailedPlots = false;

//float d0_res = 0.0554; //cm
float d0_res = 0.02152; //cm
float CUTOFF = 1.0; // Common z cut off in cm for Displaced Vertex and verification
float bendChi2Max = 14.0;
float chi2RPhiMax = 6.0;
//float chi2RZMax = 2.5;
float cosTMin = 0.96;
//float dTMax = 1.0;
float deltaZMax = 0.5;
float deltaEtaMax = 2.0;
//float deltaEtaMin = 0.05;
//float d0Min = 0.025;
float delxyMax = 0.1;
float delzMax = 0.5;

void SetPlotStyle();
void mySmallText(Double_t x, Double_t y, Color_t color, char *text);
void removeFlows(TH1F* h);
void removeFlows(TH2F* h);

class Track_Parameters
{
public:
  float pt;
  float d0;
  float dxy = -99999;
  float z0;
  float eta;
  float phi;
  float charge;
  float rho;
  int index;
  int pdgid = -99999;
  float vx;
  float vy;
  float vz;
  Track_Parameters* tp;
  float x0;
  float y0;
  int nstubs;
  float chi2rphi;
  float chi2rz;
  float bendchi2;
  float MVA1;
  float MVA2;

  float dist_calc(float x_dv, float y_dv, float x, float y){
    dxy = TMath::Sqrt((x_dv-x)*(x_dv-x) + (y_dv-y)*(y_dv-y));
    return dxy;
  }
  float x(float phi_T=0){
    return (-charge * rho * TMath::Sin(phi - charge*phi_T) + (d0 + charge * rho) * TMath::Sin(phi));
  }
  float y(float phi_T=0){
    return ( charge * rho * TMath::Cos(phi - charge*phi_T) - (d0 + charge * rho) * TMath::Cos(phi));
  }
  float z(float phi_T=0){
    float theta = 2 * TMath::ATan(TMath::Exp(-eta));
    return (z0 + rho*phi_T/TMath::Tan(theta));
  }
  float deltaPhi_T(Double_t phi1, Double_t phi2)
  {
    Double_t dPhi = phi1 - phi2;
    if (dPhi >= TMath::Pi())
      dPhi -= 2. * TMath::Pi();
    if (dPhi <= -TMath::Pi())
      dPhi += 2. * TMath::Pi();
    return dPhi;
  }
  float phi_T(float x, float y){
    float num = x - (d0 + charge * rho) * TMath::Sin(phi);
    float den = y + (d0 + charge * rho) * TMath::Cos(phi);
    return ((phi-TMath::ATan2(num,-den))/charge);
  }
  float z(float x, float y){
    float t = std::sinh(eta);
    float r = TMath::Sqrt(pow(x,2)+pow(y,2));
    return (z0+(t*r*(1+(pow(d0,2)/pow(r,2))+(1.0/6.0)*pow(r/(2*rho),2)))); // can do higher order terms if necessary from displaced math
  }
  Track_Parameters(float pt_in, float d0_in, float z0_in, float eta_in, float phi_in, int pdgid_in, float vx_in, float vy_in, float vz_in, float charge_in=0, int index_in=-1, Track_Parameters* tp_in=nullptr, int nstubs_in=0, float chi2rphi_in=0, float chi2rz_in=0, float bendchi2_in=0, float MVA1_in=0, float MVA2_in=0)
  {
    pt = pt_in;
    d0 = d0_in;
    z0 = z0_in;
    eta = eta_in;
    phi = phi_in;
    if(charge_in > 0){
      charge = 1;
    }
    else if (charge_in < 0){
      charge = -1;
    }
    else{
      charge = 0;
    }
    index = index_in;
    pdgid = pdgid_in;
    vx = vx_in;
    vy = vy_in;
    vz = vz_in;
    tp = tp_in;
    rho = fabs(1/charge_in);
    x0 = (rho+charge*d0)*TMath::Cos(phi-(charge*TMath::Pi()/2));
    y0 = (rho+charge*d0)*TMath::Sin(phi-(charge*TMath::Pi()/2));
    nstubs = nstubs_in;
    chi2rphi = chi2rphi_in;
    chi2rz = chi2rz_in;
    bendchi2 = bendchi2_in;
    MVA1 = MVA1_in;
    MVA2 = MVA2_in;
  }
  Track_Parameters(){};
  ~Track_Parameters(){};
};

constexpr bool operator==(const Track_Parameters* lhs, const Track_Parameters& rhs)
{
  return (lhs->pt==rhs.pt && lhs->d0==rhs.d0 && lhs->z0==rhs.z0 && lhs->eta==rhs.eta && lhs->phi==rhs.phi);
}
constexpr bool operator==(const Track_Parameters& lhs, const Track_Parameters* rhs)
{
  return (lhs.pt==rhs->pt && lhs.d0==rhs->d0 && lhs.z0==rhs->z0 && lhs.eta==rhs->eta && lhs.phi==rhs->phi);
}
constexpr bool operator==(const Track_Parameters& lhs, const Track_Parameters& rhs)
{
  return (lhs.pt==rhs.pt && lhs.d0==rhs.d0 && lhs.z0==rhs.z0 && lhs.eta==rhs.eta && lhs.phi==rhs.phi);
}

std::valarray<float> calcPVec(Track_Parameters a, double_t v_x, double_t v_y)
{
  std::valarray<float> r_vec = {float(v_x)-a.x0,float(v_y)-a.y0};
  std::valarray<float> p_vec = {-r_vec[1],r_vec[0]};
  if(a.charge>0){
    p_vec *= -1;
  }
  p_vec /= TMath::Sqrt(pow(p_vec[0],2)+pow(p_vec[1],2));
  p_vec *= a.pt;
  return p_vec;
}

class Vertex_Parameters
{
public:
  Double_t x_dv;
  Double_t y_dv;
  Double_t z_dv;
  float score;
  Track_Parameters a;
  Track_Parameters b;
  int inTraj;
  bool matched = false;
  std::vector<Track_Parameters> tracks = {};
  float p_mag;
  float p2_mag;
  float openingAngle;
  float R_T;
  float cos_T;
  float alpha_T;
  float d_T;
  float chi2rphidofSum;
  float chi2rzdofSum;
  float bendchi2Sum;
  float MVA1Sum;
  float MVA2Sum;
  int numStubsSum;
  float delta_z;
  float delta_eta;
  float phi;
  Vertex_Parameters(Double_t x_dv_in, Double_t y_dv_in, Double_t z_dv_in, Track_Parameters a_in, Track_Parameters b_in, float score_in=-1, int inTraj_in=4):   
    a(a_in),
    b(b_in)
  {
    x_dv = x_dv_in;
    y_dv = y_dv_in;
    z_dv = z_dv_in;
    score = score_in;
    tracks.push_back(a_in);
    tracks.push_back(b_in);
    inTraj = inTraj_in;
    std::valarray<float> p_trk_1 = calcPVec(a_in,x_dv_in,y_dv_in);
    std::valarray<float> p_trk_2 = calcPVec(b_in,x_dv_in,y_dv_in);
    std::valarray<float> p_tot = p_trk_1+p_trk_2;
    p_mag = TMath::Sqrt(pow(p_tot[0],2)+pow(p_tot[1],2));
    openingAngle = (p_trk_1[0]*p_trk_2[0]+p_trk_1[1]*p_trk_2[1]) / (TMath::Sqrt(pow(p_trk_1[0],2)+pow(p_trk_1[1],2))*TMath::Sqrt(pow(p_trk_2[0],2)+pow(p_trk_2[1],2)));
    R_T = TMath::Sqrt(pow(x_dv_in,2)+pow(y_dv_in,2));
    cos_T = (p_tot[0]*x_dv_in+p_tot[1]*y_dv_in)/(R_T*TMath::Sqrt(pow(p_tot[0],2)+pow(p_tot[1],2)));
    alpha_T = acos(cos_T);
    phi = atan2(p_tot[1],p_tot[0]);
    d_T = fabs(cos(phi)*y_dv_in-sin(phi)*x_dv_in);
    int ndof_1 = 2 * a_in.nstubs - 5;
    float chi2rphidof_1 = a_in.chi2rphi / ndof_1;
    float chi2rzdof_1 = a_in.chi2rz / ndof_1;
    float bendchi2_1 = a_in.bendchi2;
    int ndof_2 = 2 * b_in.nstubs - 5;
    float chi2rphidof_2 = b_in.chi2rphi / ndof_2;
    float chi2rzdof_2 = b_in.chi2rz / ndof_2;
    float bendchi2_2 = b_in.bendchi2;
    chi2rphidofSum = chi2rphidof_1 + chi2rphidof_2;
    chi2rzdofSum = chi2rzdof_1 + chi2rzdof_2;
    bendchi2Sum = bendchi2_1 + bendchi2_2;
    MVA1Sum = a_in.MVA1 + b_in.MVA1;
    MVA2Sum = a_in.MVA2 + b_in.MVA2;
    numStubsSum = a_in.nstubs + b_in.nstubs;
    p2_mag = pow(a_in.pt,2)+pow(b_in.pt,2);
    delta_z = fabs(a_in.z(x_dv_in,y_dv_in)-b_in.z(x_dv_in,y_dv_in));
    delta_eta = fabs(a_in.eta-b_in.eta);
  }

  void addTrack(Track_Parameters trk){
    tracks.push_back(trk);
    std::valarray<float> p_tot = {0,0};
    for(auto track : tracks){
      p_tot+= calcPVec(track,x_dv,y_dv);
    }
    p_mag = TMath::Sqrt(pow(p_tot[0],2)+pow(p_tot[1],2));
    cos_T = (p_tot[0]*x_dv+p_tot[1]*y_dv)/(R_T*TMath::Sqrt(pow(p_tot[0],2)+pow(p_tot[1],2)));
    alpha_T = acos(cos_T);
    phi = atan2(p_tot[1],p_tot[0]);
    d_T = fabs(cos(phi)*y_dv-sin(phi)*x_dv);
    int ndof = 2 * trk.nstubs - 5;
    float chi2rphidof = trk.chi2rphi / ndof;
    float chi2rzdof = trk.chi2rz / ndof;
    float bendchi2 = trk.bendchi2;
    chi2rphidofSum+= chi2rphidof;
    chi2rzdofSum+= chi2rzdof;
    bendchi2Sum+= bendchi2;
    numStubsSum+= trk.nstubs;
    p2_mag+= pow(trk.pt,2);
    MVA1Sum+= trk.MVA1;
    MVA2Sum+= trk.MVA2;
  }

  Vertex_Parameters(){};
  ~Vertex_Parameters(){};
};

constexpr bool operator==(const Vertex_Parameters& lhs, const Vertex_Parameters& rhs)
{
  return (lhs.x_dv==rhs.x_dv && lhs.y_dv==rhs.y_dv && lhs.z_dv==rhs.z_dv);
}

void displayProgress(long current, long max)
{
  using std::cerr;
  if (max < 2500)
    return;
  if (current % (max / 2500) != 0 && current < max - 1)
    return;

  int width = 52; // Hope the terminal is at least that wide.
  int barWidth = width - 2;
  cerr << "\x1B[2K";    // Clear line
  cerr << "\x1B[2000D"; // Cursor left
  cerr << '[';
  for (int i = 0; i < barWidth; ++i)
    {
      if (i < barWidth * current / max)
	{
	  cerr << '=';
	}
      else
	{
	  cerr << ' ';
	}
    }
  cerr << ']';
  cerr << " " << Form("%8d/%8d (%5.2f%%)", (int)current, (int)max, 100.0 * current / max);
  cerr.flush();
}

template <typename T, typename S = TH1F>
void raiseMax(T *hist1, S *hist2=nullptr, T *hist3=nullptr, T *hist4=nullptr)
{
  Double_t max = hist1->GetMaximum();
  if(hist2!=nullptr){
    Double_t max2 = hist2->GetMaximum();
    if(max2>max) max = max2;
  }
  if(hist3!=nullptr){
    Double_t max3 = hist3->GetMaximum();
    if(max3>max) max = max3;
  }
  if(hist4!=nullptr){
    Double_t max4 = hist4->GetMaximum();
    if(max4>max) max = max4;
  }
  if(max>0.0){
    hist1->GetYaxis()->SetRangeUser(0.,1.2*max);
    if(hist2!=nullptr) hist2->GetYaxis()->SetRangeUser(0.,1.2*max);
    if(hist3!=nullptr) hist3->GetYaxis()->SetRangeUser(0.,1.2*max);
    if(hist4!=nullptr) hist4->GetYaxis()->SetRangeUser(0.,1.2*max);
  }
}

template <typename T, typename S>
void drawSame(T *hist1, S *hist2, T *hist3=nullptr, T *hist4=nullptr)
{
  if(hist1->GetMaximum()!=0.0){
    hist1->Draw("HIST");
    hist2->Draw("HIST,SAME");
    if(hist3!=nullptr) hist3->Draw("HIST,SAME");
    if(hist4!=nullptr) hist4->Draw("HIST,SAME");
  }
  else if(hist2->GetMaximum()!=0.0){
    hist2->Draw("HIST");
    if(hist3!=nullptr) hist3->Draw("HIST,SAME");
    if(hist4!=nullptr) hist4->Draw("HIST,SAME");
  }
  else if(hist3!=nullptr){
    if(hist3->GetMaximum()!=0.0){
      hist3->Draw("HIST");
      if(hist4!=nullptr) hist4->Draw("HIST,SAME");
    }
  }
  else if(hist4!=nullptr){
    if(hist4->GetMaximum()!=0.0){
      hist4->Draw("HIST");
    }
  }
  else{
    hist1->Draw("HIST");
  }
}

bool ComparePtTrack(Track_Parameters a, Track_Parameters b) { return a.pt > b.pt; }
bool CompareZ0Track(Track_Parameters a, Track_Parameters b) { return a.z0 > b.z0; }
bool CompareD0Track(Track_Parameters a, Track_Parameters b) { return a.d0 > b.d0; }
bool ComparePtVert(Vertex_Parameters v1, Vertex_Parameters v2) {return v1.a.pt > v2.a.pt; }
bool CompareDelzVert(Vertex_Parameters v1, Vertex_Parameters v2) {return v1.delta_z > v2.delta_z; }
bool CompareDtVert(Vertex_Parameters v1, Vertex_Parameters v2) {return v1.d_T > v2.d_T; }
bool CompareChi2rphidofSumVert(Vertex_Parameters v1, Vertex_Parameters v2) {return v1.chi2rphidofSum > v2.chi2rphidofSum; }
bool CompareRtVert(Vertex_Parameters v1, Vertex_Parameters v2) {return v1.R_T > v2.R_T; }

template<typename T>
std::vector<T> linspace(T start, T end, int num){
  std::vector<T> out;
  T delta = (end - start) / (num-1);
  for(int i=0; i<num-1; i++){
    out.push_back(start+delta*i);
  }
  out.push_back(end);
  return out;
}

Double_t deltaPhi(Double_t phi1, Double_t phi2)
{
  Double_t dPhi = phi1 - phi2;
  if (dPhi > TMath::Pi())
    dPhi -= 2. * TMath::Pi();
  if (dPhi < -TMath::Pi())
    dPhi += 2. * TMath::Pi();
  return dPhi;
}

Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2)
{
  Double_t dEta, dPhi;
  dEta = eta1 - eta2;
  dPhi = deltaPhi(phi1, phi2);
  return sqrt(dEta * dEta + dPhi * dPhi);
}

Double_t dist(Double_t x1, Double_t y1 , Double_t x2=0, Double_t y2=0){ // Distance between 2 points
  return (TMath::Sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)));
}

Double_t dist_Vertex(Double_t x_vtx, Double_t y_vtx, Track_Parameters a){ // Distance between track and displaced vertex
  float R = dist(x_vtx,y_vtx,a.x0,a.y0);
  return (fabs(R-(a.rho)));
}

Double_t dist_TPs(Track_Parameters* a, Track_Parameters* b); // Closest distance between 2 tracks
Double_t dist_TPs(Track_Parameters a, Track_Parameters b); // Closest distance between 2 tracks
bool CompareDeltaXY(Vertex_Parameters v1, Vertex_Parameters v2) {return dist_TPs(v1.a,v1.b) < dist_TPs(v2.a,v2.b); }

Int_t calcVertex(Track_Parameters a, Track_Parameters b, Double_t &x_vtx, Double_t &y_vtx, Double_t &z_vtx); 
Int_t Vertex(Track_Parameters a, Track_Parameters b, Double_t &x_vtx, Double_t &y_vtx, Double_t &z_vtx); 
// Identify the displaced vertex (x_vtx,y_vtx,z_vtx) and return the status 
//-2 = Circles with same center. No Intersection
//-1 = Circles don't Intersect. A point on the line connecting the centers is chosen.
// 0 = Only 1 Intersection not satisfying Z cutoff
// 1 = Only 1 Intersection satisfying Z cutoff
// 2 = Only 1 Intersection detectable dist(x,y)<20
// 3 = 2 Intersections 

void Analyzer_DisplacedMuon(TString inputFilePath,
			    TString outputDir,
			    float TP_maxD0 = 1.9,
			    float TP_minD0 = 0.0004196)
{
  TChain *tree = new TChain("L1TrackNtuple/eventTree"); 
  tree->Add(inputFilePath);
  //TChain *vertTree = new TChain("L1TrackNtuple/dispVertTree");
  //vertTree->Add(inputFilePath);
  std::string inputFileString(inputFilePath.Data());
  inputFileString = inputFileString.substr(inputFileString.find_last_of("/")+1);
  TString inputFile(inputFileString);
  std::cout<<"input: "<<inputFile<<std::endl;
  gROOT->SetBatch();
  gErrorIgnoreLevel = kWarning;

  SetPlotStyle();

  //float TP_minPt = 9.0;
  //TP_maxD0 = 1.5;
  //TP_minD0 = 0.08;
  //float TP_maxZ0 = 14.0;
  //float maxChi2rzdof = 5.0;
  //float maxChi2rphidof = 5.00;
  float barrelEta = 0.95;

  //general preselection cuts
  float TP_maxEta = 2.4;
  //float trk_maxBendChi2 = 9.0;
  float maxChi2rzdof = 3.0;
  float minMVA2 = 0.2;
  //float minMVA2_D = 0.95;
  float minMVA1 = 0.2;
  float minMVA1_D = 0.8;
  float TP_minPt = 3.0;
  //float TP_maxZ0 = 20.0;
  //TP_maxD0 = 2.0;
  //barrel preselection cuts
  float TP_minD0_barrel = 0.06;
  float maxChi2rphidof_barrel = 5.0;
  //disk preselection cuts
  float TP_minD0_disk = 0.08;
  //float minChi2rphidof_disk = 10000.0;

  vector<float>   *trk_pt;
  vector<float>   *trk_eta;
  vector<float>   *trk_phi;
  vector<float>   *trk_d0;
  vector<float>   *trk_rinv;
  vector<float>   *trk_z0;
  vector<float>   *trk_chi2;
  vector<float>   *trk_chi2rphi;
  vector<float>   *trk_chi2rz;
  vector<float>   *trk_bendchi2;
  vector<int>     *trk_nstub;
  vector<int>     *trk_lhits;
  vector<int>     *trk_dhits;
  vector<int>     *trk_seed;
  vector<int>     *trk_hitpattern;
  vector<unsigned int> *trk_phiSector;
  vector<int>     *trk_genuine;
  vector<int>     *trk_loose;
  vector<int>     *trk_unknown;
  vector<int>     *trk_combinatoric;
  vector<int>     *trk_fake;
  vector<float>     *trk_MVA1;
  vector<float>     *trk_MVA2;
  vector<int>     *trk_matchtp_pdgid;
  vector<bool>     *trk_matchtp_isHToMu;
  vector<bool>     *trk_matchtp_isHToB;
  vector<float>   *trk_matchtp_pt;
  vector<float>   *trk_matchtp_eta;
  vector<float>   *trk_matchtp_phi;
  vector<float>   *trk_matchtp_z0;
  vector<float>   *trk_matchtp_dxy;
  vector<float>   *trk_matchtp_d0;
  vector<float>   *trk_matchtp_x;
  vector<float>   *trk_matchtp_y;
  vector<float>   *trk_matchtp_z;
  vector<float>   *tp_pt;
  vector<float>   *tp_eta;
  vector<float>   *tp_phi;
  vector<float>   *tp_dxy;
  vector<float>   *tp_d0;
  vector<float>   *tp_z0;
  vector<float>   *tp_x;
  vector<float>   *tp_y;
  vector<float>   *tp_z;
  vector<float>   *tp_d0_prod;
  vector<float>   *tp_z0_prod;
  vector<int>     *tp_pdgid;
  vector<bool>    *tp_isHToMu;
  vector<bool>    *tp_isHToB;
  vector<int>     *tp_nmatch;
  vector<int>     *tp_nstub;
  vector<int>     *tp_eventid;
  vector<int>     *tp_charge;
  vector<float>   *matchtrk_pt;
  vector<float>   *matchtrk_eta;
  vector<float>   *matchtrk_phi;
  vector<float>   *matchtrk_z0;
  vector<float>   *matchtrk_d0;
  vector<float>   *matchtrk_rinv;
  vector<float>   *matchtrk_chi2;
  vector<float>   *matchtrk_chi2rphi;
  vector<float>   *matchtrk_chi2rz;
  vector<float>   *matchtrk_bendchi2;
  vector<float>   *matchtrk_MVA1;
  vector<float>   *matchtrk_MVA2;
  vector<int>     *matchtrk_nstub;
  vector<int>     *matchtrk_lhits;
  vector<int>     *matchtrk_dhits;
  vector<int>     *matchtrk_seed;
  vector<int>     *matchtrk_hitpattern;
#if 0
  vector<float>* vert_pt;
  vector<float>* vert_z0;
  vector<float>* vert_d0;
  vector<float>* vert_eta;
  vector<float>* vert_phi;
  vector<float>* vert_delta_z;
  vector<float>* vert_R_T;
  vector<float>* vert_cos_T;
  vector<float>* vert_d_T;
  vector<float>* vert_chi2rzdofSum;
  vector<int>* vert_numStubsSum;
  vector<float>* vert_chi2rphidofSum;
  vector<float>* vert_minD0;
  vector<float>* vert_sumPt;
  vector<float>* vert_score;
  vector<float>* vert_x;
  vector<float>* vert_y;
  vector<float>* vert_z;
#endif
  TBranch        *b_trk_pt;
  TBranch        *b_trk_eta;
  TBranch        *b_trk_phi;
  TBranch        *b_trk_d0;
  TBranch        *b_trk_rinv;
  TBranch        *b_trk_z0;
  TBranch        *b_trk_chi2;
  TBranch        *b_trk_chi2rphi;
  TBranch        *b_trk_chi2rz;
  TBranch        *b_trk_bendchi2;
  TBranch        *b_trk_nstub;
  TBranch        *b_trk_lhits;
  TBranch        *b_trk_dhits;
  TBranch        *b_trk_seed;
  TBranch        *b_trk_hitpattern;
  TBranch        *b_trk_phiSector;
  TBranch        *b_trk_genuine;
  TBranch        *b_trk_loose;
  TBranch        *b_trk_unknown;
  TBranch        *b_trk_combinatoric;
  TBranch        *b_trk_fake;
  TBranch        *b_trk_MVA1;
  TBranch        *b_trk_MVA2;
  TBranch        *b_trk_matchtp_pdgid;
  TBranch        *b_trk_matchtp_isHToMu;
  TBranch        *b_trk_matchtp_isHToB;
  TBranch        *b_trk_matchtp_pt;
  TBranch        *b_trk_matchtp_eta;
  TBranch        *b_trk_matchtp_phi;
  TBranch        *b_trk_matchtp_z0;
  TBranch        *b_trk_matchtp_dxy;
  TBranch        *b_trk_matchtp_d0;
  TBranch        *b_trk_matchtp_x;
  TBranch        *b_trk_matchtp_y;
  TBranch        *b_trk_matchtp_z;
  TBranch        *b_tp_pt;
  TBranch        *b_tp_eta;
  TBranch        *b_tp_phi;
  TBranch        *b_tp_dxy;
  TBranch        *b_tp_d0;
  TBranch        *b_tp_z0;
  TBranch        *b_tp_x;
  TBranch        *b_tp_y;
  TBranch        *b_tp_z;
  TBranch        *b_tp_d0_prod;
  TBranch        *b_tp_z0_prod;
  TBranch        *b_tp_pdgid;
  TBranch        *b_tp_isHToMu;
  TBranch        *b_tp_isHToB;
  TBranch        *b_tp_nmatch;
  TBranch        *b_tp_nstub;
  TBranch        *b_tp_eventid;
  TBranch        *b_tp_charge;
  TBranch        *b_matchtrk_pt;
  TBranch        *b_matchtrk_eta;
  TBranch        *b_matchtrk_phi;
  TBranch        *b_matchtrk_z0;
  TBranch        *b_matchtrk_d0;
  TBranch        *b_matchtrk_rinv;
  TBranch        *b_matchtrk_chi2;
  TBranch        *b_matchtrk_chi2rphi;
  TBranch        *b_matchtrk_chi2rz;
  TBranch        *b_matchtrk_bendchi2;
  TBranch        *b_matchtrk_MVA1;
  TBranch        *b_matchtrk_MVA2;
  TBranch        *b_matchtrk_nstub;
  TBranch        *b_matchtrk_lhits;
  TBranch        *b_matchtrk_dhits;
  TBranch        *b_matchtrk_seed;
  TBranch        *b_matchtrk_hitpattern;
#if 0
  TBranch *b_vert_pt;
  TBranch *b_vert_z0;
  TBranch *b_vert_d0;
  TBranch *b_vert_eta;
  TBranch *b_vert_phi;
  TBranch *b_vert_delta_z;
  TBranch *b_vert_R_T;
  TBranch *b_vert_cos_T;
  TBranch *b_vert_d_T;
  TBranch *b_vert_chi2rzdofSum;
  TBranch *b_vert_numStubsSum;
  TBranch *b_vert_chi2rphidofSum;
  TBranch *b_vert_minD0;
  TBranch *b_vert_sumPt;
  TBranch *b_vert_score;
  TBranch *b_vert_x;
  TBranch *b_vert_y;
  TBranch *b_vert_z;
#endif
  trk_pt = 0;
  trk_eta = 0;
  trk_phi = 0;
  trk_d0 = 0;
  trk_rinv = 0;
  trk_z0 = 0;
  trk_chi2 = 0;
  trk_chi2rphi = 0;
  trk_chi2rz = 0;
  trk_bendchi2 = 0;
  trk_nstub = 0;
  trk_lhits = 0;
  trk_dhits = 0;
  trk_seed = 0;
  trk_hitpattern = 0;
  trk_phiSector = 0;
  trk_genuine = 0;
  trk_loose = 0;
  trk_unknown = 0;
  trk_combinatoric = 0;
  trk_fake = 0;
  trk_MVA1 = 0;
  trk_MVA2 = 0;
  trk_matchtp_pdgid = 0;
  trk_matchtp_isHToMu = 0;
  trk_matchtp_isHToB = 0;
  trk_matchtp_pt = 0;
  trk_matchtp_eta = 0;
  trk_matchtp_phi = 0;
  trk_matchtp_z0 = 0;
  trk_matchtp_dxy = 0;
  trk_matchtp_d0 = 0;
  trk_matchtp_x = 0;
  trk_matchtp_y = 0;
  trk_matchtp_z = 0;
  tp_pt = 0;
  tp_eta = 0;
  tp_phi = 0;
  tp_dxy = 0;
  tp_d0 = 0;
  tp_z0 = 0;
  tp_x = 0;
  tp_y = 0;
  tp_z = 0;
  tp_d0_prod = 0;
  tp_z0_prod = 0;
  tp_pdgid = 0;
  tp_isHToMu = 0;
  tp_isHToB = 0;
  tp_nmatch = 0;
  tp_nstub = 0;
  tp_eventid = 0;
  tp_charge = 0;
  matchtrk_pt = 0;
  matchtrk_eta = 0;
  matchtrk_phi = 0;
  matchtrk_z0 = 0;
  matchtrk_d0 = 0;
  matchtrk_rinv = 0;
  matchtrk_chi2 = 0;
  matchtrk_chi2rphi = 0;
  matchtrk_chi2rz = 0;
  matchtrk_bendchi2 = 0;
  matchtrk_MVA1 = 0;
  matchtrk_MVA2 = 0;
  matchtrk_nstub = 0;
  matchtrk_lhits = 0;
  matchtrk_dhits = 0;
  matchtrk_seed = 0;
  matchtrk_hitpattern = 0;
#if 0
  vert_pt = 0;
  vert_z0 = 0;
  vert_d0 = 0;
  vert_eta = 0;
  vert_phi = 0;
  vert_delta_z = 0;
  vert_R_T = 0;
  vert_cos_T = 0;
  vert_d_T = 0;
  vert_chi2rzdofSum = 0;
  vert_numStubsSum = 0;
  vert_chi2rphidofSum = 0;
  vert_minD0 = 0;
  vert_sumPt = 0;
  vert_score = 0;
  vert_x = 0;
  vert_y = 0;
  vert_z = 0;
#endif
  tree->SetMakeClass(1);
  tree->SetBranchAddress("trk_pt", &trk_pt, &b_trk_pt);
  tree->SetBranchAddress("trk_eta", &trk_eta, &b_trk_eta);
  tree->SetBranchAddress("trk_phi", &trk_phi, &b_trk_phi);
  tree->SetBranchAddress("trk_d0", &trk_d0, &b_trk_d0);
  tree->SetBranchAddress("trk_rinv", &trk_rinv, &b_trk_rinv);
  tree->SetBranchAddress("trk_z0", &trk_z0, &b_trk_z0);
  tree->SetBranchAddress("trk_chi2", &trk_chi2, &b_trk_chi2);
  tree->SetBranchAddress("trk_chi2rphi", &trk_chi2rphi, &b_trk_chi2rphi);
  tree->SetBranchAddress("trk_chi2rz", &trk_chi2rz, &b_trk_chi2rz);
  tree->SetBranchAddress("trk_bendchi2", &trk_bendchi2, &b_trk_bendchi2);
  tree->SetBranchAddress("trk_nstub", &trk_nstub, &b_trk_nstub);
  tree->SetBranchAddress("trk_lhits", &trk_lhits, &b_trk_lhits);
  tree->SetBranchAddress("trk_dhits", &trk_dhits, &b_trk_dhits);
  tree->SetBranchAddress("trk_seed", &trk_seed, &b_trk_seed);
  tree->SetBranchAddress("trk_hitpattern", &trk_hitpattern, &b_trk_hitpattern);
  tree->SetBranchAddress("trk_phiSector", &trk_phiSector, &b_trk_phiSector);
  tree->SetBranchAddress("trk_genuine", &trk_genuine, &b_trk_genuine);
  tree->SetBranchAddress("trk_loose", &trk_loose, &b_trk_loose);
  tree->SetBranchAddress("trk_unknown", &trk_unknown, &b_trk_unknown);
  tree->SetBranchAddress("trk_combinatoric", &trk_combinatoric, &b_trk_combinatoric);
  tree->SetBranchAddress("trk_fake", &trk_fake, &b_trk_fake);
  tree->SetBranchAddress("trk_MVA1", &trk_MVA1, &b_trk_MVA1);
  tree->SetBranchAddress("trk_MVA2", &trk_MVA2, &b_trk_MVA2);
  tree->SetBranchAddress("trk_matchtp_pdgid", &trk_matchtp_pdgid, &b_trk_matchtp_pdgid);
  tree->SetBranchAddress("trk_matchtp_isHToMu", &trk_matchtp_isHToMu, &b_trk_matchtp_isHToMu);
  tree->SetBranchAddress("trk_matchtp_isHToB", &trk_matchtp_isHToB, &b_trk_matchtp_isHToB);
  tree->SetBranchAddress("trk_matchtp_pt", &trk_matchtp_pt, &b_trk_matchtp_pt);
  tree->SetBranchAddress("trk_matchtp_eta", &trk_matchtp_eta, &b_trk_matchtp_eta);
  tree->SetBranchAddress("trk_matchtp_phi", &trk_matchtp_phi, &b_trk_matchtp_phi);
  tree->SetBranchAddress("trk_matchtp_z0", &trk_matchtp_z0, &b_trk_matchtp_z0);
  tree->SetBranchAddress("trk_matchtp_dxy", &trk_matchtp_dxy, &b_trk_matchtp_dxy);
  tree->SetBranchAddress("trk_matchtp_d0", &trk_matchtp_d0, &b_trk_matchtp_d0);
  tree->SetBranchAddress("trk_matchtp_x", &trk_matchtp_x, &b_trk_matchtp_x);
  tree->SetBranchAddress("trk_matchtp_y", &trk_matchtp_y, &b_trk_matchtp_y);
  tree->SetBranchAddress("trk_matchtp_z", &trk_matchtp_z, &b_trk_matchtp_z);
  tree->SetBranchAddress("tp_pt", &tp_pt, &b_tp_pt);
  tree->SetBranchAddress("tp_eta", &tp_eta, &b_tp_eta);
  tree->SetBranchAddress("tp_phi", &tp_phi, &b_tp_phi);
  tree->SetBranchAddress("tp_dxy", &tp_dxy, &b_tp_dxy);
  tree->SetBranchAddress("tp_d0", &tp_d0, &b_tp_d0);
  tree->SetBranchAddress("tp_z0", &tp_z0, &b_tp_z0);
  tree->SetBranchAddress("tp_x", &tp_x, &b_tp_x);
  tree->SetBranchAddress("tp_y", &tp_y, &b_tp_y);
  tree->SetBranchAddress("tp_z", &tp_z, &b_tp_z);
  tree->SetBranchAddress("tp_d0_prod", &tp_d0_prod, &b_tp_d0_prod);
  tree->SetBranchAddress("tp_z0_prod", &tp_z0_prod, &b_tp_z0_prod);
  tree->SetBranchAddress("tp_pdgid", &tp_pdgid, &b_tp_pdgid);
  tree->SetBranchAddress("tp_isHToMu", &tp_isHToMu, &b_tp_isHToMu);
  tree->SetBranchAddress("tp_isHToB", &tp_isHToB, &b_tp_isHToB);
  tree->SetBranchAddress("tp_nmatch", &tp_nmatch, &b_tp_nmatch);
  tree->SetBranchAddress("tp_nstub", &tp_nstub, &b_tp_nstub);
  tree->SetBranchAddress("tp_eventid", &tp_eventid, &b_tp_eventid);
  tree->SetBranchAddress("tp_charge", &tp_charge, &b_tp_charge);
  tree->SetBranchAddress("matchtrk_pt", &matchtrk_pt, &b_matchtrk_pt);
  tree->SetBranchAddress("matchtrk_eta", &matchtrk_eta, &b_matchtrk_eta);
  tree->SetBranchAddress("matchtrk_phi", &matchtrk_phi, &b_matchtrk_phi);
  tree->SetBranchAddress("matchtrk_z0", &matchtrk_z0, &b_matchtrk_z0);
  tree->SetBranchAddress("matchtrk_d0", &matchtrk_d0, &b_matchtrk_d0);
  tree->SetBranchAddress("matchtrk_rinv", &matchtrk_rinv, &b_matchtrk_rinv);
  tree->SetBranchAddress("matchtrk_chi2", &matchtrk_chi2, &b_matchtrk_chi2);
  tree->SetBranchAddress("matchtrk_chi2rphi", &matchtrk_chi2rphi, &b_matchtrk_chi2rphi);
  tree->SetBranchAddress("matchtrk_chi2rz", &matchtrk_chi2rz, &b_matchtrk_chi2rz);
  tree->SetBranchAddress("matchtrk_bendchi2", &matchtrk_bendchi2, &b_matchtrk_bendchi2);
  tree->SetBranchAddress("matchtrk_MVA1", &matchtrk_MVA1, &b_matchtrk_MVA1);
  tree->SetBranchAddress("matchtrk_MVA2", &matchtrk_MVA2, &b_matchtrk_MVA2);
  tree->SetBranchAddress("matchtrk_nstub", &matchtrk_nstub, &b_matchtrk_nstub);
  tree->SetBranchAddress("matchtrk_lhits", &matchtrk_lhits, &b_matchtrk_lhits);
  tree->SetBranchAddress("matchtrk_dhits", &matchtrk_dhits, &b_matchtrk_dhits);
  tree->SetBranchAddress("matchtrk_seed", &matchtrk_seed, &b_matchtrk_seed);
  tree->SetBranchAddress("matchtrk_hitpattern", &matchtrk_hitpattern, &b_matchtrk_hitpattern);
#if 0
  vertTree->SetBranchAddress("vert_pt", &vert_pt, &b_vert_pt);
  vertTree->SetBranchAddress("vert_z0", &vert_z0, &b_vert_z0);
  vertTree->SetBranchAddress("vert_d0", &vert_d0, &b_vert_d0);
  vertTree->SetBranchAddress("vert_eta", &vert_eta, &b_vert_eta);
  vertTree->SetBranchAddress("vert_phi", &vert_phi, &b_vert_phi);
  vertTree->SetBranchAddress("vert_delta_z", &vert_delta_z, &b_vert_delta_z);
  vertTree->SetBranchAddress("vert_R_T", &vert_R_T, &b_vert_R_T);
  vertTree->SetBranchAddress("vert_cos_T", &vert_cos_T, &b_vert_cos_T);
  vertTree->SetBranchAddress("vert_d_T", &vert_d_T, &b_vert_d_T);
  vertTree->SetBranchAddress("vert_chi2rzdofSum", &vert_chi2rzdofSum, &b_vert_chi2rzdofSum);
  vertTree->SetBranchAddress("vert_numStubsSum", &vert_numStubsSum, &b_vert_numStubsSum);
  vertTree->SetBranchAddress("vert_chi2rphidofSum", &vert_chi2rphidofSum, &b_vert_chi2rphidofSum);
  vertTree->SetBranchAddress("vert_minD0", &vert_minD0, &b_vert_minD0);
  vertTree->SetBranchAddress("vert_sumPt", &vert_sumPt, &b_vert_sumPt);
  vertTree->SetBranchAddress("vert_score", &vert_score, &b_vert_score);
  vertTree->SetBranchAddress("vert_x", &vert_x, &b_vert_x);
  vertTree->SetBranchAddress("vert_y", &vert_y, &b_vert_y);
  vertTree->SetBranchAddress("vert_z", &vert_z, &b_vert_z);
  tree->AddFriend(vertTree);
#endif  
  using milli = std::chrono::milliseconds;
  typedef boost::variant<vector<float>**,vector<int>**> boostVector;
  std::map<TString, std::pair<boostVector, float> > preselCuts{{{"maxEta",{&trk_eta,2.4}}, {"maxChi2rzdof",{&trk_chi2rz,3.0}}, {"maxChi2rphidof_barrel",{&trk_chi2rphi, 5.0}}, {"minMVA2",{&trk_MVA2,0.2}},{"minMVA1",{&trk_MVA1,0.2}}, {"minMVA1_D",{&trk_MVA1,0.8}}, {"minNumStub_overlap",{&trk_nstub,5}}, {"minPt",{&trk_pt,3.0}}, {"minD0_barrel",{&trk_d0,0.06}}, {"minD0_disk",{&trk_d0,0.08}}}};
  std::map<TString, std::pair<boostVector, float> > preselCutsTP{{{"maxEta",{&tp_eta,2.4}}, {"minPt",{&tp_pt,3.0}}, {"minD0_barrel",{&tp_d0,0.06}}, {"minD0_disk",{&tp_d0,0.08}}}};

  std::map<std::vector<TString>,std::pair<boostVector,std::vector<float>> > varCutFlows = { {{"d0","cm"},{&trk_d0,{200,-2.0,2.0}}},
													{{"pt","GeV"},{&trk_pt,{200,0.0,100.0}}},
													{{"eta",""},{&trk_eta,{50,-2.5,2.5}}},
													{{"z0","cm"},{&trk_z0,{100,-20.0,20.0}}},
													{{"phi",""},{&trk_phi,{100,-2*TMath::Pi(),2*TMath::Pi()}}},
													{{"sectorPhi",""},{&trk_phi,{100,-1.0,1.0}}},
													{{"MVA1",""},{&trk_MVA1,{100,0.0,1.0}}},
													{{"MVA2",""},{&trk_MVA2,{100,0.0,1.0}}},
													{{"chi2rphidof",""},{&trk_chi2rphi,{100,0.0,6.0}}},
													{{"chi2rzdof",""},{&trk_chi2rz,{100,0.0,6.0}}},
													{{"bendchi2",""},{&trk_bendchi2,{100,0.0,10.0}}}
  };
  
  std::map<std::vector<TString>,std::pair<boostVector,std::vector<float>> > varCutFlowsTP = { {{"d0","cm"},{&tp_d0,{200,-2.0,2.0}}},
													  {{"pt","GeV"},{&tp_pt,{200,0.0,100.0}}},
													  {{"eta",""},{&tp_eta,{50,-2.5,2.5}}},
													  {{"z0","cm"},{&tp_z0,{100,-20.0,20.0}}},
													  {{"phi",""},{&tp_phi,{100,-2*TMath::Pi(),2*TMath::Pi()}}},
													  {{"sectorPhi",""},{&tp_phi,{100,-1.0,1.0}}},
													  {{"dxy","cm"},{&tp_dxy,{50,-2.0,2.0}}}
  };

  std::map<std::vector<TString>,std::pair<std::vector<boostVector>,std::vector<float>> > varCutFlows2D = { {{"d0","cm","pt","GeV"},{{&trk_d0,&trk_pt},{200,-2.0,2.0,200,0.0,30.0}}},
														       {{"eta","","pt","GeV"},{{&trk_eta,&trk_pt},{200,-2.4,2.4,200,0.0,30.0}}},
														       {{"d0","cm","eta",""},{{&trk_d0,&trk_eta},{200,-2.0,2.0,200,-2.4,2.4}}},
														       {{"eta","","nstub",""},{{&trk_eta,&trk_nstub},{200,-2.4,2.4,7,0.0,7.0}}}
  };

  std::map<std::vector<TString>,std::pair<std::vector<boostVector>,std::vector<float>> > varCutFlowsTP2D = { {{"d0","cm","pt","GeV"},{{&tp_d0,&tp_pt},{200,-2.0,2.0,200,0.0,30.0}}},
															 {{"eta","","pt","GeV"},{{&tp_eta,&tp_pt},{200,-2.4,2.4,200,0.0,30.0}}},
															 {{"d0","cm","eta",""},{{&tp_d0,&tp_eta},{200,-2.0,2.0,200,-2.4,2.4}}},
															 {{"eta","","nstub",""},{{&tp_eta,&tp_nstub},{200,-2.4,2.4,7,0.0,7.0}}}
  };
    
  std::vector<TString> trackType = {"primary","np","fake","PU","notHiggs"};
  std::vector<TString> tpType = {"primary","np","PU","notHiggs","match",""};
  std::vector<TString> plotModifiers = {"","_H","_L","_P","_D","_barrel","_disk"};
  if(!detailedPlots) plotModifiers = {""};
  TH1F* preselCutFlows[varCutFlows.size()][trackType.size()][preselCuts.size()][plotModifiers.size()];
  TH2F* preselCutFlows2D[varCutFlows2D.size()][trackType.size()][preselCuts.size()][plotModifiers.size()];
  TH1F* preselCutFlowsTP[varCutFlowsTP.size()][tpType.size()][preselCutsTP.size()][plotModifiers.size()];
  TH2F* preselCutFlowsTP2D[varCutFlowsTP2D.size()][tpType.size()][preselCutsTP.size()][plotModifiers.size()];
  std::map<string,int> numPartCutFlows[trackType.size()][preselCuts.size()];
  std::map<string,int> numPartCutFlowsTP[tpType.size()][preselCutsTP.size()];
  
  int it_counter = 0;
  for(auto it=varCutFlows.cbegin(); it!=varCutFlows.cend(); ++it){
    for(uint i=0; i<trackType.size(); ++i){
      int jt_counter = 0;
      for(auto jt=preselCuts.cbegin(); jt!=preselCuts.cend(); ++jt){
	for(uint j=0; j<plotModifiers.size(); ++j){
	  TString name = "h_trk_"+it->first.at(0)+"_"+trackType[i]+"_"+jt->first+"Cut"+plotModifiers[j];
	  float binWidth = (it->second.second.at(2) - it->second.second.at(1)) / it->second.second.at(0);
	  TString binLabel = std::to_string(binWidth);
	  TString labels = name+"; Track "+it->first.at(0)+" ("+it->first.at(1)+") ; Events / "+binLabel+" "+it->first.at(1);
	  TH1F* hist = new TH1F(name,labels,it->second.second.at(0),it->second.second.at(1),it->second.second.at(2));
	  preselCutFlows[it_counter][i][jt_counter][j] = hist;
	}
	jt_counter++;
      }
    }
    it_counter++;
  }
  
  it_counter = 0;
  for(auto it=varCutFlows2D.cbegin(); it!=varCutFlows2D.cend(); ++it){
    for(uint i=0; i<trackType.size(); ++i){
      int jt_counter = 0;
      for(auto jt=preselCuts.cbegin(); jt!=preselCuts.cend(); ++jt){
	for(uint j=0; j<plotModifiers.size(); ++j){
	  TString name = "h_trk_"+it->first.at(2)+"_vs_"+it->first.at(0)+"_"+trackType[i]+"_"+jt->first+"Cut"+plotModifiers[j];
	  TString labels = name+"; Track "+it->first.at(0)+" ("+it->first.at(1)+") ; Track "+it->first.at(2)+" ("+it->first.at(3)+")";
	  TH2F* hist = new TH2F(name,labels,it->second.second.at(0),it->second.second.at(1),it->second.second.at(2),it->second.second.at(3),it->second.second.at(4),it->second.second.at(5));
	  preselCutFlows2D[it_counter][i][jt_counter][j] = hist;	  
	}
	jt_counter++;
      }
    }
    it_counter++;
  }
  
  it_counter = 0;
  for(auto it=varCutFlowsTP.cbegin(); it!=varCutFlowsTP.cend(); ++it){
    for(uint i=0; i<tpType.size(); ++i){
      int jt_counter = 0;
      for(auto jt=preselCutsTP.cbegin(); jt!=preselCutsTP.cend(); ++jt){
	TString cutName = jt->first;
	for(uint j=0; j<plotModifiers.size(); ++j){
	  TString name = "h_tp_"+it->first.at(0)+"_"+tpType[i]+"_"+jt->first+"Cut"+plotModifiers[j];
	  float binWidth = (it->second.second.at(2) - it->second.second.at(1)) / it->second.second.at(0);
	  TString binLabel = std::to_string(binWidth);
	  TString labels = name+"; Tp "+it->first.at(0)+" ("+it->first.at(1)+") ; Events / "+binLabel+" "+it->first.at(1);
	  TH1F* hist = new TH1F(name,labels,it->second.second.at(0),it->second.second.at(1),it->second.second.at(2));
	  preselCutFlowsTP[it_counter][i][jt_counter][j] = hist;
	}
	jt_counter++;
      }
    }
    it_counter++;
  }
  
  it_counter = 0;
  for(auto it=varCutFlowsTP2D.cbegin(); it!=varCutFlowsTP2D.cend(); ++it){
    for(uint i=0; i<tpType.size(); ++i){
      int jt_counter = 0;
      for(auto jt=preselCutsTP.cbegin(); jt!=preselCutsTP.cend(); ++jt){
	TString cutName = jt->first;
	for(uint j=0; j<plotModifiers.size(); ++j){
	  TString name = "h_tp_"+it->first.at(2)+"_vs_"+it->first.at(0)+"_"+tpType[i]+"_"+jt->first+"Cut"+plotModifiers[j];
	  TString labels = name+"; Tp "+it->first.at(0)+" ("+it->first.at(1)+") ; Tp "+it->first.at(2)+" ("+it->first.at(3)+")";
	  TH2F* hist = new TH2F(name,labels,it->second.second.at(0),it->second.second.at(1),it->second.second.at(2),it->second.second.at(3),it->second.second.at(4),it->second.second.at(5));
	  preselCutFlowsTP2D[it_counter][i][jt_counter][j] = hist;
	}
	jt_counter++;
      }
    }
    it_counter++;
  }
  
  TH1F *h_trk_d0 = new TH1F("h_trk_d0","h_trk_d0; Track d_{0} Distribution (cm) ; Events / 0.02 cm",200,-2,2);
  TH1F *h_trk_d0_primary = new TH1F("h_trk_d0_primary","h_trk_d0_primary; Track d_{0} Distribution (cm) ; Events / 0.02 cm",200,-2,2);
  TH1F *h_trk_d0_primary_noCuts = new TH1F("h_trk_d0_primary_noCuts","h_trk_d0_primary_noCuts; Track d_{0} Distribution (cm) ; Events / 0.02 cm",200,-2,2);
  TH1F *h_trk_d0_primary_noCuts_zoomOut = new TH1F("h_trk_d0_primary_noCuts_zoomOut","h_trk_d0_primary_noCuts_zoomOut; Track d_{0} Distribution (cm) ; Events / 0.2 cm",200,-20,20);
  TH1F *h_trk_d0_primary_noCuts_barrel = new TH1F("h_trk_d0_primary_noCuts_barrel","h_trk_d0_primary_noCuts_barrel; Track d_{0} Distribution (cm) (#eta<=0.95); Events / 0.02 cm",200,-2,2);
  TH1F *h_trk_d0_primary_noCuts_disk = new TH1F("h_trk_d0_primary_noCuts_disk","h_trk_d0_primary_noCuts_disk; Track d_{0} Distribution (cm) (#eta>0.95); Events / 0.02 cm",200,-2,2);
  TH1F *h_trk_d0_primary_noCuts_H = new TH1F("h_trk_d0_primary_noCuts_H","h_trk_d0_primary_noCuts_H; Track d_{0} Distribution (cm) (p_{T}>10 GeV); Events / 0.02 cm",200,-2,2);
  TH1F *h_trk_d0_primary_noCuts_L = new TH1F("h_trk_d0_primary_noCuts_L","h_trk_d0_primary_noCuts_L; Track d_{0} Distribution (cm) (p_{T}<10 GeV); Events / 0.02 cm",200,-2,2);
  TH1F *h_trk_d0_primary_qualCuts = new TH1F("h_trk_d0_primary_qualCuts","h_trk_d0_primary_qualCuts; Track d_{0} Distribution (cm) ; Events / 0.02 cm",200,-2,2);
  TH1F *h_trk_d0_primary_allCuts = new TH1F("h_trk_d0_primary_allCuts","h_trk_d0_primary_allCuts; Track d_{0} Distribution (cm) ; Events / 0.02 cm",200,-2,2);
  TH1F *h_trk_d0_primary_allCuts_zoomOut = new TH1F("h_trk_d0_primary_allCuts_zoomOut","h_trk_d0_primary_allCuts_zoomOut; Track d_{0} Distribution (cm) ; Events / 0.2 cm",200,-20,20);
  TH1F *h_trk_d0_primary_allCuts_barrel = new TH1F("h_trk_d0_primary_allCuts_barrel","h_trk_d0_primary_allCuts_barrel; Track d_{0} Distribution (cm) (#eta<=0.95); Events / 0.02 cm",200,-2,2);
  TH1F *h_trk_d0_primary_allCuts_disk = new TH1F("h_trk_d0_primary_allCuts_disk","h_trk_d0_primary_allCuts_disk; Track d_{0} Distribution (cm) (#eta>0.95); Events / 0.02 cm",200,-2,2);
  TH1F *h_trk_ptWeightedD0_primary_qualCuts = new TH1F("h_trk_ptWeightedD0_primary_qualCuts","h_trk_ptWeightedD0_primary_qualCuts; Track p_{T}*d_{0} Distribution (GeV cm) ; Events / 0.02 GeV cm",200,-2,2);
  // np = not primary
  TH1F *h_trk_d0_np = new TH1F("h_trk_d0_np","h_trk_d0_np; Track d_{0} Distribution (cm) ; Events / 0.02 cm",200,-2,2);
  TH1F *h_trk_d0_np_noCuts = new TH1F("h_trk_d0_np_noCuts","h_trk_d0_np_noCuts; Track d_{0} Distribution (cm) ; Events / 0.02 cm",200,-2,2);
  TH1F *h_trk_d0_np_noCuts_zoomOut = new TH1F("h_trk_d0_np_noCuts_zoomOut","h_trk_d0_np_noCuts_zoomOut; Track d_{0} Distribution (cm) ; Events / 0.2 cm",200,-20,20);
  TH1F *h_trk_d0_np_noCuts_barrel = new TH1F("h_trk_d0_np_noCuts_barrel","h_trk_d0_np_noCuts_barrel; Track d_{0} Distribution (cm) (#eta<=0.95); Events / 0.02 cm",200,-2,2);
  TH1F *h_trk_d0_np_noCuts_disk = new TH1F("h_trk_d0_np_noCuts_disk","h_trk_d0_np_noCuts_disk; Track d_{0} Distribution (cm) (#eta>0.95); Events / 0.02 cm",200,-2,2);
  TH1F *h_trk_d0_np_noCuts_H = new TH1F("h_trk_d0_np_noCuts_H","h_trk_d0_np_noCuts_H; Track d_{0} Distribution (cm) (p_{T}>10 GeV); Events / 0.02 cm",200,-2,2);
  TH1F *h_trk_d0_np_noCuts_L = new TH1F("h_trk_d0_np_noCuts_L","h_trk_d0_np_noCuts_L; Track d_{0} Distribution (cm) (p_{T}<10 GeV); Events / 0.02 cm",200,-2,2);
  TH1F *h_trk_d0_fake_noCuts = new TH1F("h_trk_d0_fake_noCuts","h_trk_d0_fake_noCuts; Track d_{0} Distribution (cm) ; Events / 0.02 cm",200,-2,2);
  TH1F *h_trk_d0_PU_noCuts = new TH1F("h_trk_d0_PU_noCuts","h_trk_d0_PU_noCuts; Track d_{0} Distribution (cm) ; Events / 0.02 cm",200,-2,2);
  TH1F *h_trk_d0_notHiggs_noCuts = new TH1F("h_trk_d0_notHiggs_noCuts","h_trk_d0_notHiggs_noCuts; Track d_{0} Distribution (cm) ; Events / 0.02 cm",200,-2,2);
  TH1F *h_trk_d0_np_qualCuts = new TH1F("h_trk_d0_np_qualCuts","h_trk_d0_np_qualCuts; Track d_{0} Distribution (cm) ; Events / 0.02 cm",200,-2,2);
  TH1F *h_trk_d0_np_allCuts = new TH1F("h_trk_d0_np_allCuts","h_trk_d0_np_allCuts; Track d_{0} Distribution (cm) ; Events / 0.02 cm",200,-2,2);
  TH1F *h_trk_d0_np_allCuts_zoomOut = new TH1F("h_trk_d0_np_allCuts_zoomOut","h_trk_d0_np_allCuts_zoomOut; Track d_{0} Distribution (cm) ; Events / 0.2 cm",200,-20,20);
  TH1F *h_trk_d0_np_allCuts_barrel = new TH1F("h_trk_d0_np_allCuts_barrel","h_trk_d0_np_allCuts_barrel; Track d_{0} Distribution (cm) (#eta<=0.95); Events / 0.02 cm",200,-2,2);
  TH1F *h_trk_d0_np_allCuts_disk = new TH1F("h_trk_d0_np_allCuts_disk","h_trk_d0_np_allCuts_disk; Track d_{0} Distribution (cm) (#eta>0.95); Events / 0.02 cm",200,-2,2);
  TH1F *h_trk_ptWeightedD0_np_qualCuts = new TH1F("h_trk_ptWeightedD0_np_qualCuts","h_trk_ptWeightedD0_np_qualCuts; Track p_{T}*d_{0} Distribution (GeV cm) ; Events / 0.02 GeV cm",200,-2,2);
  TH1F *h_trk_d0_fake_qualCuts = new TH1F("h_trk_d0_fake_qualCuts","h_trk_d0_fake_qualCuts; Track d_{0} Distribution (cm) ; Events / 0.02 cm",200,-2,2);
  TH1F *h_trk_d0_PU_qualCuts = new TH1F("h_trk_d0_PU_qualCuts","h_trk_d0_PU_qualCuts; Track d_{0} Distribution (cm) ; Events / 0.02 cm",200,-2,2);
  TH1F *h_trk_d0_notHiggs_qualCuts = new TH1F("h_trk_d0_notHiggs_qualCuts","h_trk_d0_notHiggs_qualCuts; Track d_{0} Distribution (cm) ; Events / 0.02 cm",200,-2,2);
  TH1F *h_trk_pt = new TH1F("h_trk_pt","h_trk_pt; Track p_{T} Distribution (GeV); Events / 0.6 GeV", 100, 0, 60.0);
  TH1F *h_trk_pt_noCuts = new TH1F("h_trk_pt_noCuts","h_trk_pt_noCuts; Track p_{T} Distribution (GeV); Events / 0.6 GeV", 100, 0, 60.0);
  TH1F *h_trk_pt_oldCuts = new TH1F("h_trk_pt_oldCuts","h_trk_pt_oldCuts; Track p_{T} Distribution (GeV); Events / 0.6 GeV", 100, 0, 60.0);
  TH1F *h_trk_pt_allCuts = new TH1F("h_trk_pt_allCuts","h_trk_pt_allCuts; Track p_{T} Distribution (GeV); Events / 0.6 GeV", 100, 0, 60.0);
  TH1F *h_trk_pt_primary = new TH1F("h_trk_pt_primary","h_trk_pt_primary; Track p_{T} Distribution (GeV); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_primary_noCuts = new TH1F("h_trk_pt_primary_noCuts","h_trk_pt_primary_noCuts; Track p_{T} Distribution (GeV); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_primary_noCuts_barrel = new TH1F("h_trk_pt_primary_noCuts_barrel","h_trk_pt_primary_noCuts_barrel; Track p_{T} Distribution (GeV) (#eta<=0.95); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_primary_noCuts_disk = new TH1F("h_trk_pt_primary_noCuts_disk","h_trk_pt_primary_noCuts_disk; Track p_{T} Distribution (GeV) (#eta>0.95); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_primary_qualCuts = new TH1F("h_trk_pt_primary_qualCuts","h_trk_pt_primary_qualCuts; Track p_{T} Distribution (GeV); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_primary_ptCuts = new TH1F("h_trk_pt_primary_ptCuts","h_trk_pt_primary_ptCuts; Track p_{T} Distribution (GeV); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_primary_d0Cuts = new TH1F("h_trk_pt_primary_d0Cuts","h_trk_pt_primary_d0Cuts; Track p_{T} Distribution (GeV); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_primary_chi2rzdofCuts = new TH1F("h_trk_pt_primary_chi2rzdofCuts","h_trk_pt_primary_chi2rzdofCuts; Track p_{T} Distribution (GeV); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_primary_bendchi2Cuts = new TH1F("h_trk_pt_primary_bendchi2Cuts","h_trk_pt_primary_bendchi2Cuts; Track p_{T} Distribution (GeV); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_primary_chi2rphidofCuts = new TH1F("h_trk_pt_primary_chi2rphidofCuts","h_trk_pt_primary_chi2rphidofCuts; Track p_{T} Distribution (GeV); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_primary_nstubCuts = new TH1F("h_trk_pt_primary_nstubCuts","h_trk_pt_primary_nstubCuts; Track p_{T} Distribution (GeV); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_primary_z0Cuts = new TH1F("h_trk_pt_primary_z0Cuts","h_trk_pt_primary_z0Cuts; Track p_{T} Distribution (GeV); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_primary_allCuts = new TH1F("h_trk_pt_primary_allCuts","h_trk_pt_primary_allCuts; Track p_{T} Distribution (GeV); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_primary_allCuts_barrel = new TH1F("h_trk_pt_primary_allCuts_barrel","h_trk_pt_primary_allCuts_barrel; Track p_{T} Distribution (GeV) (#eta<=0.95); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_primary_allCuts_disk = new TH1F("h_trk_pt_primary_allCuts_disk","h_trk_pt_primary_allCuts_disk; Track p_{T} Distribution (GeV) (#eta>0.95); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_ptIso4_primary_allCuts = new TH1F("h_trk_ptIso4_primary_allCuts","h_trk_ptIso4_primary_allCuts; Track p_{T} Isolation Distribution; Events / 0.02", 1000, 0, 20.0);
  TH1F *h_trk_ptIso8_primary_allCuts = new TH1F("h_trk_ptIso8_primary_allCuts","h_trk_ptIso8_primary_allCuts; Track p_{T} Isolation Distribution; Events / 0.02", 1000, 0, 20.0);
  TH1F *h_trk_pt_np = new TH1F("h_trk_pt_np","h_trk_pt_np; Track p_{T} Distribution (GeV); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_np_noCuts = new TH1F("h_trk_pt_np_noCuts","h_trk_pt_np_noCuts; Track p_{T} Distribution (GeV); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_np_noCuts_barrel = new TH1F("h_trk_pt_np_noCuts_barrel","h_trk_pt_np_noCuts_barrel; Track p_{T} Distribution (GeV) (#eta<=0.95); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_np_noCuts_disk = new TH1F("h_trk_pt_np_noCuts_disk","h_trk_pt_np_noCuts_disk; Track p_{T} Distribution (GeV) (#eta>0.95); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_fake_noCuts = new TH1F("h_trk_pt_fake_noCuts","h_trk_pt_fake_noCuts; Track p_{T} Distribution (GeV); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_PU_noCuts = new TH1F("h_trk_pt_PU_noCuts","h_trk_pt_PU_noCuts; Track p_{T} Distribution (GeV); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_notHiggs_noCuts = new TH1F("h_trk_pt_notHiggs_noCuts","h_trk_pt_notHiggs_noCuts; Track p_{T} Distribution (GeV); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_np_qualCuts = new TH1F("h_trk_pt_np_qualCuts","h_trk_pt_np_qualCuts; Track p_{T} Distribution (GeV); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_fake_qualCuts = new TH1F("h_trk_pt_fake_qualCuts","h_trk_pt_fake_qualCuts; Track p_{T} Distribution (GeV); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_PU_qualCuts = new TH1F("h_trk_pt_PU_qualCuts","h_trk_pt_PU_qualCuts; Track p_{T} Distribution (GeV); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_notHiggs_qualCuts = new TH1F("h_trk_pt_notHiggs_qualCuts","h_trk_pt_notHiggs_qualCuts; Track p_{T} Distribution (GeV); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_np_ptCuts = new TH1F("h_trk_pt_np_ptCuts","h_trk_pt_np_ptCuts; Track p_{T} Distribution (GeV); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_np_d0Cuts = new TH1F("h_trk_pt_np_d0Cuts","h_trk_pt_np_d0Cuts; Track p_{T} Distribution (GeV); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_np_chi2rzdofCuts = new TH1F("h_trk_pt_np_chi2rzdofCuts","h_trk_pt_np_chi2rzdofCuts; Track p_{T} Distribution (GeV); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_np_bendchi2Cuts = new TH1F("h_trk_pt_np_bendchi2Cuts","h_trk_pt_np_bendchi2Cuts; Track p_{T} Distribution (GeV); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_np_chi2rphidofCuts = new TH1F("h_trk_pt_np_chi2rphidofCuts","h_trk_pt_np_chi2rphidofCuts; Track p_{T} Distribution (GeV); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_np_nstubCuts = new TH1F("h_trk_pt_np_nstubCuts","h_trk_pt_np_nstubCuts; Track p_{T} Distribution (GeV); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_np_z0Cuts = new TH1F("h_trk_pt_np_z0Cuts","h_trk_pt_np_z0Cuts; Track p_{T} Distribution (GeV); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_np_allCuts = new TH1F("h_trk_pt_np_allCuts","h_trk_pt_np_allCuts; Track p_{T} Distribution (GeV); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_np_allCuts_barrel = new TH1F("h_trk_pt_np_allCuts_barrel","h_trk_pt_np_allCuts_barrel; Track p_{T} Distribution (GeV) (#eta<=0.95); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_pt_np_allCuts_disk = new TH1F("h_trk_pt_np_allCuts_disk","h_trk_pt_np_allCuts_disk; Track p_{T} Distribution (GeV) (#eta>0.95); Events / 0.5 GeV", 200, 0, 100.0);
  TH1F *h_trk_ptIso4_np_allCuts = new TH1F("h_trk_ptIso4_np_allCuts","h_trk_ptIso4_np_allCuts; Track p_{T} Isolation Distribution; Events / 0.02", 1000, 0, 20.0);
  TH1F *h_trk_ptIso8_np_allCuts = new TH1F("h_trk_ptIso8_np_allCuts","h_trk_ptIso8_np_allCuts; Track p_{T} Isolation Distribution; Events / 0.02", 1000, 0, 20.0);
  TH1F *h_trk_eta = new TH1F("h_trk_eta","h_trk_eta; Track #eta Distribution; Events / 0.1", 50, -2.5, 2.5);
  TH1F *h_trk_eta_primary = new TH1F("h_trk_eta_primary","h_trk_eta_primary; Track #eta Distribution; Events / 0.1", 50, -2.5, 2.5);
  TH1F *h_trk_eta_primary_noCuts = new TH1F("h_trk_eta_primary_noCuts","h_trk_eta_primary_noCuts; Track #eta Distribution; Events / 0.1", 50, -2.5, 2.5);
  TH1F *h_trk_eta_primary_noCuts_H = new TH1F("h_trk_eta_primary_noCuts_H","h_trk_eta_primary_noCuts_H; Track #eta Distribution (p_{T}>10 GeV); Events / 0.1", 50, -2.5, 2.5);
  TH1F *h_trk_eta_primary_noCuts_L = new TH1F("h_trk_eta_primary_noCuts_L","h_trk_eta_primary_noCuts_L; Track #eta Distribution (p_{T}<10 GeV); Events / 0.1", 50, -2.5, 2.5);
  TH1F *h_trk_eta_primary_qualCuts = new TH1F("h_trk_eta_primary_qualCuts","h_trk_eta_primary_qualCuts; Track #eta Distribution; Events / 0.1", 50, -2.5, 2.5);
  TH1F *h_trk_eta_primary_allCuts = new TH1F("h_trk_eta_primary_allCuts","h_trk_eta_primary_allCuts; Track #eta Distribution; Events / 0.1", 50, -2.5, 2.5);
  TH1F *h_trk_eta_np = new TH1F("h_trk_eta_np","h_trk_eta_np; Track #eta Distribution; Events / 0.1", 50, -2.5, 2.5);
  TH1F *h_trk_eta_np_noCuts = new TH1F("h_trk_eta_np_noCuts","h_trk_eta_np_noCuts; Track #eta Distribution; Events / 0.1", 50, -2.5, 2.5);
  TH1F *h_trk_eta_np_noCuts_H = new TH1F("h_trk_eta_np_noCuts_H","h_trk_eta_np_noCuts_H; Track #eta Distribution (p_{T}>10 GeV); Events / 0.1", 50, -2.5, 2.5);
  TH1F *h_trk_eta_np_noCuts_L = new TH1F("h_trk_eta_np_noCuts_L","h_trk_eta_np_noCuts_L; Track #eta Distribution (p_{T}<10 GeV); Events / 0.1", 50, -2.5, 2.5);
  TH1F *h_trk_eta_fake_noCuts = new TH1F("h_trk_eta_fake_noCuts","h_trk_eta_fake_noCuts; Track #eta Distribution; Events / 0.1", 50, -2.5, 2.5);
  TH1F *h_trk_eta_PU_noCuts = new TH1F("h_trk_eta_PU_noCuts","h_trk_eta_PU_noCuts; Track #eta Distribution; Events / 0.1", 50, -2.5, 2.5);
  TH1F *h_trk_eta_notHiggs_noCuts = new TH1F("h_trk_eta_notHiggs_noCuts","h_trk_eta_notHiggs_noCuts; Track #eta Distribution; Events / 0.1", 50, -2.5, 2.5);
  TH1F *h_trk_eta_np_qualCuts = new TH1F("h_trk_eta_np_qualCuts","h_trk_eta_np_qualCuts; Track #eta Distribution; Events / 0.1", 50, -2.5, 2.5);
  TH1F *h_trk_eta_fake_qualCuts = new TH1F("h_trk_eta_fake_qualCuts","h_trk_eta_fake_qualCuts; Track #eta Distribution; Events / 0.1", 50, -2.5, 2.5);
  TH1F *h_trk_eta_PU_qualCuts = new TH1F("h_trk_eta_PU_qualCuts","h_trk_eta_PU_qualCuts; Track #eta Distribution; Events / 0.1", 50, -2.5, 2.5);
  TH1F *h_trk_eta_notHiggs_qualCuts = new TH1F("h_trk_eta_notHiggs_qualCuts","h_trk_eta_notHiggs_qualCuts; Track #eta Distribution; Events / 0.1", 50, -2.5, 2.5);
  TH1F *h_trk_eta_np_allCuts = new TH1F("h_trk_eta_np_allCuts","h_trk_eta_np_allCuts; Track #eta Distribution; Events / 0.1", 50, -2.5, 2.5);
  TH1F* h_trk_z0 = new TH1F("h_trk_z0","h_trk_z0; Track z_{0} Distribution (cm); Events / 0.1 cm", 200, 0, 20);
  TH1F* h_trk_z0_primary = new TH1F("h_trk_z0_primary","h_trk_z0_primary; Track z_{0} Distribution (cm); Events / 0.1 cm", 200, 0, 20);
  TH1F* h_trk_z0_primary_noCuts = new TH1F("h_trk_z0_primary_noCuts","h_trk_z0_primary_noCuts; Track z_{0} Distribution (cm); Events / 0.4 cm", 100, -20, 20);
  TH1F* h_trk_z0_primary_noCuts_zoomOut = new TH1F("h_trk_z0_primary_noCuts_zoomOut","h_trk_z0_primary_noCuts_zoomOut; Track z_{0} Distribution (cm); Events / 1.0 cm", 100, -50, 50);
  TH1F* h_trk_z0_primary_noCuts_barrel = new TH1F("h_trk_z0_primary_noCuts_barrel","h_trk_z0_primary_noCuts_barrel; Track z_{0} Distribution (cm) (#eta<=0.95); Events / 0.2 cm", 100, -20, 20);
  TH1F* h_trk_z0_primary_noCuts_disk = new TH1F("h_trk_z0_primary_noCuts_disk","h_trk_z0_primary_noCuts_disk; Track z_{0} Distribution (cm) (#eta>0.95); Events / 0.2 cm", 100, -20, 20);
  TH1F* h_trk_z0_primary_noCuts_H = new TH1F("h_trk_z0_primary_noCuts_H","h_trk_z0_primary_noCuts_H; Track z_{0} Distribution (cm) (p_{T}>10 GeV); Events / 0.1 cm", 100, 0, 20);
  TH1F* h_trk_z0_primary_noCuts_L = new TH1F("h_trk_z0_primary_noCuts_L","h_trk_z0_primary_noCuts_L; Track z_{0} Distribution (cm) (p_{T}<10 GeV); Events / 0.1 cm", 100, 0, 20);
  TH1F* h_trk_z0_primary_qualCuts = new TH1F("h_trk_z0_primary_qualCuts","h_trk_z0_primary_qualCuts; Track z_{0} Distribution (cm); Events / 0.1 cm", 200, 0, 20);
  TH1F* h_trk_z0_primary_allCuts = new TH1F("h_trk_z0_primary_allCuts","h_trk_z0_primary_allCuts; Track z_{0} Distribution (cm); Events / 0.2 cm", 100, -20, 20);
  TH1F* h_trk_z0_primary_allCuts_zoomOut = new TH1F("h_trk_z0_primary_allCuts_zoomOut","h_trk_z0_primary_allCuts_zoomOut; Track z_{0} Distribution (cm); Events / 1.0 cm", 100, -50, 50);
  TH1F* h_trk_z0_primary_allCuts_barrel = new TH1F("h_trk_z0_primary_allCuts_barrel","h_trk_z0_primary_allCuts_barrel; Track z_{0} Distribution (cm) (#eta<=0.95); Events / 0.2 cm", 100, -20, 20);
  TH1F* h_trk_z0_primary_allCuts_disk = new TH1F("h_trk_z0_primary_allCuts_disk","h_trk_z0_primary_allCuts_disk; Track z_{0} Distribution (cm) (#eta>0.95); Events / 0.2 cm", 100, -20, 20);
  TH1F* h_trk_z0_np = new TH1F("h_trk_z0_np","h_trk_z0_np; Track z_{0} Distribution (cm); Events / 0.1 cm", 200, 0, 20);
  TH1F* h_trk_z0_np_noCuts = new TH1F("h_trk_z0_np_noCuts","h_trk_z0_np_noCuts; Track z_{0} Distribution (cm); Events / 0.2 cm", 100, -20, 20);
  TH1F* h_trk_z0_np_noCuts_zoomOut = new TH1F("h_trk_z0_np_noCuts_zoomOut","h_trk_z0_np_noCuts_zoomOut; Track z_{0} Distribution (cm); Events / 1.0 cm", 100, -50, 50);
  TH1F* h_trk_z0_np_noCuts_barrel = new TH1F("h_trk_z0_np_noCuts_barrel","h_trk_z0_np_noCuts_barrel; Track z_{0} Distribution (cm) (#eta<=0.95); Events / 0.2 cm", 100, -20, 20);
  TH1F* h_trk_z0_np_noCuts_disk = new TH1F("h_trk_z0_np_noCuts_disk","h_trk_z0_np_noCuts_disk; Track z_{0} Distribution (cm) (#eta>0.95); Events / 0.2 cm", 100, -20, 20);
  TH1F* h_trk_z0_np_noCuts_H = new TH1F("h_trk_z0_np_noCuts_H","h_trk_z0_np_noCuts_H; Track z_{0} Distribution (cm) (p_{T}>10 GeV); Events / 0.1 cm", 200, 0, 20);
  TH1F* h_trk_z0_np_noCuts_L = new TH1F("h_trk_z0_np_noCuts_L","h_trk_z0_np_noCuts_L; Track z_{0} Distribution (cm) (p_{T}<10 GeV); Events / 0.1 cm", 200, 0, 20);
  TH1F* h_trk_z0_fake_noCuts = new TH1F("h_trk_z0_fake_noCuts","h_trk_z0_fake_noCuts; Track z_{0} Distribution (cm); Events / 0.2 cm", 100, -20, 20);
  TH1F* h_trk_z0_PU_noCuts = new TH1F("h_trk_z0_PU_noCuts","h_trk_z0_PU_noCuts; Track z_{0} Distribution (cm); Events / 0.2 cm", 100, -20, 20);
  TH1F* h_trk_z0_notHiggs_noCuts = new TH1F("h_trk_z0_notHiggs_noCuts","h_trk_z0_notHiggs_noCuts; Track z_{0} Distribution (cm); Events / 0.2 cm", 100, -20, 20);
  TH1F* h_trk_z0_np_qualCuts = new TH1F("h_trk_z0_np_qualCuts","h_trk_z0_np_qualCuts; Track z_{0} Distribution (cm); Events / 0.1 cm", 200, 0, 20);
  TH1F* h_trk_z0_fake_qualCuts = new TH1F("h_trk_z0_fake_qualCuts","h_trk_z0_fake_qualCuts; Track z_{0} Distribution (cm); Events / 0.1 cm", 200, 0, 20);
  TH1F* h_trk_z0_PU_qualCuts = new TH1F("h_trk_z0_PU_qualCuts","h_trk_z0_PU_qualCuts; Track z_{0} Distribution (cm); Events / 0.1 cm", 200, 0, 20);
  TH1F* h_trk_z0_notHiggs_qualCuts = new TH1F("h_trk_z0_notHiggs_qualCuts","h_trk_z0_notHiggs_qualCuts; Track z_{0} Distribution (cm); Events / 0.1 cm", 200, 0, 20);
  TH1F* h_trk_z0_np_allCuts = new TH1F("h_trk_z0_np_allCuts","h_trk_z0_np_allCuts; Track z_{0} Distribution (cm); Events / 0.2 cm", 100, -20, 20);
  TH1F* h_trk_z0_np_allCuts_zoomOut = new TH1F("h_trk_z0_np_allCuts_zoomOut","h_trk_z0_np_allCuts_zoomOut; Track z_{0} Distribution (cm); Events / 1.0 cm", 100, -50, 50);
  TH1F* h_trk_z0_np_allCuts_barrel = new TH1F("h_trk_z0_np_allCuts_barrel","h_trk_z0_np_allCuts_barrel; Track z_{0} Distribution (cm) (#eta<=0.95); Events / 0.2 cm", 100, -20, 20);
  TH1F* h_trk_z0_np_allCuts_disk = new TH1F("h_trk_z0_np_allCuts_disk","h_trk_z0_np_allCuts_disk; Track z_{0} Distribution (cm) (#eta>0.95); Events / 0.2 cm", 100, -20, 20);
  TH1F *h_trk_phi = new TH1F("h_trk_phi","h_trk_phi; Track #phi_{0} Distribution; Events / 0.1256",100,-2*TMath::Pi(),2*TMath::Pi());
  TH1F *h_trk_phi_primary = new TH1F("h_trk_phi_primary","h_trk_phi_primary; Track #phi_{0} Distribution; Events / 0.1256",100,-2*TMath::Pi(),2*TMath::Pi());
  TH1F *h_trk_phi_primary_noCuts = new TH1F("h_trk_phi_primary_noCuts","h_trk_phi_primary_noCuts; Track #phi_{0} Distribution; Events / 0.1256",100,-2*TMath::Pi(),2*TMath::Pi());
  TH1F *h_trk_phi_primary_noCuts_barrel = new TH1F("h_trk_phi_primary_noCuts_barrel","h_trk_phi_primary_noCuts_barrel; Track #phi_{0} Distribution (#eta<=0.95); Events / 0.1256",100,-2*TMath::Pi(),2*TMath::Pi());
  TH1F *h_trk_phi_primary_noCuts_disk = new TH1F("h_trk_phi_primary_noCuts_disk","h_trk_phi_primary_noCuts_disk; Track #phi_{0} Distribution (#eta>0.95); Events / 0.1256",100,-2*TMath::Pi(),2*TMath::Pi());
  TH1F *h_trk_phi_primary_noCuts_H = new TH1F("h_trk_phi_primary_noCuts_H","h_trk_phi_primary_noCuts_H; Track #phi_{0} Distribution (p_{T}>10 GeV); Events / 0.1256",100,-2*TMath::Pi(),2*TMath::Pi());
  TH1F *h_trk_phi_primary_noCuts_L = new TH1F("h_trk_phi_primary_noCuts_L","h_trk_phi_primary_noCuts_L; Track #phi_{0} Distribution (p_{T}<10 GeV); Events / 0.1256",100,-2*TMath::Pi(),2*TMath::Pi());
  TH1F *h_trk_phi_primary_qualCuts = new TH1F("h_trk_phi_primary_qualCuts","h_trk_phi_primary_qualCuts; Track #phi_{0} Distribution; Events / 0.1256",100,-2*TMath::Pi(),2*TMath::Pi());
  TH1F *h_trk_phi_primary_allCuts = new TH1F("h_trk_phi_primary_allCuts","h_trk_phi_primary_allCuts; Track #phi_{0} Distribution; Events / 0.1256",100,-2*TMath::Pi(),2*TMath::Pi());
  TH1F *h_trk_phi_primary_allCuts_barrel = new TH1F("h_trk_phi_primary_allCuts_barrel","h_trk_phi_primary_allCuts_barrel; Track #phi_{0} Distribution (#eta<=0.95); Events / 0.1256",100,-2*TMath::Pi(),2*TMath::Pi());
  TH1F *h_trk_phi_primary_allCuts_disk = new TH1F("h_trk_phi_primary_allCuts_disk","h_trk_phi_primary_allCuts_disk; Track #phi_{0} Distribution (#eta>0.95); Events / 0.1256",100,-2*TMath::Pi(),2*TMath::Pi());
  TH1F *h_trk_phi_np = new TH1F("h_trk_phi_np","h_trk_phi_np; Track #phi_{0} Distribution; Events / 0.1256",100,-2*TMath::Pi(),2*TMath::Pi());
  TH1F *h_trk_phi_np_noCuts = new TH1F("h_trk_phi_np_noCuts","h_trk_phi_np_noCuts; Track #phi_{0} Distribution; Events / 0.1256",100,-2*TMath::Pi(),2*TMath::Pi());
  TH1F *h_trk_phi_np_noCuts_barrel = new TH1F("h_trk_phi_np_noCuts_barrel","h_trk_phi_np_noCuts_barrel; Track #phi_{0} Distribution (#eta<=0.95); Events / 0.1256",100,-2*TMath::Pi(),2*TMath::Pi());
  TH1F *h_trk_phi_np_noCuts_disk = new TH1F("h_trk_phi_np_noCuts_disk","h_trk_phi_np_noCuts_disk; Track #phi_{0} Distribution (#eta>0.95); Events / 0.1256",100,-2*TMath::Pi(),2*TMath::Pi());
  TH1F *h_trk_phi_np_noCuts_H = new TH1F("h_trk_phi_np_noCuts_H","h_trk_phi_np_noCuts_H; Track #phi_{0} Distribution (p_{T}>10 GeV); Events / 0.1256",100,-2*TMath::Pi(),2*TMath::Pi());
  TH1F *h_trk_phi_np_noCuts_L = new TH1F("h_trk_phi_np_noCuts_L","h_trk_phi_np_noCuts_L; Track #phi_{0} Distribution (p_{T}<10 GeV); Events / 0.1256",100,-2*TMath::Pi(),2*TMath::Pi());
  TH1F *h_trk_phi_fake_noCuts = new TH1F("h_trk_phi_fake_noCuts","h_trk_phi_fake_noCuts; Track #phi_{0} Distribution; Events / 0.1256",100,-2*TMath::Pi(),2*TMath::Pi());
  TH1F *h_trk_phi_PU_noCuts = new TH1F("h_trk_phi_PU_noCuts","h_trk_phi_PU_noCuts; Track #phi_{0} Distribution; Events / 0.1256",100,-2*TMath::Pi(),2*TMath::Pi());
  TH1F *h_trk_phi_notHiggs_noCuts = new TH1F("h_trk_phi_notHiggs_noCuts","h_trk_phi_notHiggs_noCuts; Track #phi_{0} Distribution; Events / 0.1256",100,-2*TMath::Pi(),2*TMath::Pi());
  TH1F *h_trk_phi_np_qualCuts = new TH1F("h_trk_phi_np_qualCuts","h_trk_phi_np_qualCuts; Track #phi_{0} Distribution; Events / 0.1256",100,-2*TMath::Pi(),2*TMath::Pi());
  TH1F *h_trk_phi_fake_qualCuts = new TH1F("h_trk_phi_fake_qualCuts","h_trk_phi_fake_qualCuts; Track #phi_{0} Distribution; Events / 0.1256",100,-2*TMath::Pi(),2*TMath::Pi());
  TH1F *h_trk_phi_PU_qualCuts = new TH1F("h_trk_phi_PU_qualCuts","h_trk_phi_PU_qualCuts; Track #phi_{0} Distribution; Events / 0.1256",100,-2*TMath::Pi(),2*TMath::Pi());
  TH1F *h_trk_phi_notHiggs_qualCuts = new TH1F("h_trk_phi_notHiggs_qualCuts","h_trk_phi_notHiggs_qualCuts; Track #phi_{0} Distribution; Events / 0.1256",100,-2*TMath::Pi(),2*TMath::Pi());
  TH1F *h_trk_phi_np_allCuts = new TH1F("h_trk_phi_np_allCuts","h_trk_phi_np_allCuts; Track #phi_{0} Distribution; Events / 0.1256",100,-2*TMath::Pi(),2*TMath::Pi());
  TH1F *h_trk_phi_np_allCuts_barrel = new TH1F("h_trk_phi_np_allCuts_barrel","h_trk_phi_np_allCuts_barrel; Track #phi_{0} Distribution (#eta<=0.95); Events / 0.1256",100,-2*TMath::Pi(),2*TMath::Pi());
  TH1F *h_trk_phi_np_allCuts_disk = new TH1F("h_trk_phi_np_allCuts_disk","h_trk_phi_np_allCuts_disk; Track #phi_{0} Distribution (#eta>0.95); Events / 0.1256",100,-2*TMath::Pi(),2*TMath::Pi());
  TH1F *h_trk_sectorPhi = new TH1F("h_trk_sectorPhi","h_trk_sectorPhi; Track #sectorPhi_{0} Distribution; Events / 0.02",100,-1,1);
  TH1F *h_trk_sectorPhi_primary = new TH1F("h_trk_sectorPhi_primary","h_trk_sectorPhi_primary; Track #sectorPhi_{0} Distribution; Events / 0.02",100,-1,1);
  TH1F *h_trk_sectorPhi_primary_noCuts = new TH1F("h_trk_sectorPhi_primary_noCuts","h_trk_sectorPhi_primary_noCuts; Track #sectorPhi_{0} Distribution; Events / 0.02",100,-1,1);
  TH1F *h_trk_sectorPhi_primary_noCuts_barrel = new TH1F("h_trk_sectorPhi_primary_noCuts_barrel","h_trk_sectorPhi_primary_noCuts_barrel; Track #sectorPhi_{0} Distribution (#eta<=0.95); Events / 0.02",100,-1,1);
  TH1F *h_trk_sectorPhi_primary_noCuts_disk = new TH1F("h_trk_sectorPhi_primary_noCuts_disk","h_trk_sectorPhi_primary_noCuts_disk; Track #sectorPhi_{0} Distribution (#eta>0.95); Events / 0.02",100,-1,1);
  TH1F *h_trk_sectorPhi_primary_noCuts_H = new TH1F("h_trk_sectorPhi_primary_noCuts_H","h_trk_sectorPhi_primary_noCuts_H; Track #sectorPhi_{0} Distribution (p_{T}>10 GeV); Events / 0.02",100,-1,1);
  TH1F *h_trk_sectorPhi_primary_noCuts_L = new TH1F("h_trk_sectorPhi_primary_noCuts_L","h_trk_sectorPhi_primary_noCuts_L; Track #sectorPhi_{0} Distribution (p_{T}<10 GeV); Events / 0.02",100,-1,1);
  TH1F *h_trk_sectorPhi_primary_qualCuts = new TH1F("h_trk_sectorPhi_primary_qualCuts","h_trk_sectorPhi_primary_qualCuts; Track #sectorPhi_{0} Distribution; Events / 0.02",100,-1,1);
  TH1F *h_trk_sectorPhi_primary_allCuts = new TH1F("h_trk_sectorPhi_primary_allCuts","h_trk_sectorPhi_primary_allCuts; Track #sectorPhi_{0} Distribution; Events / 0.02",100,-1,1);
  TH1F *h_trk_sectorPhi_primary_allCuts_barrel = new TH1F("h_trk_sectorPhi_primary_allCuts_barrel","h_trk_sectorPhi_primary_allCuts_barrel; Track #sectorPhi_{0} Distribution (#eta<=0.95); Events / 0.02",100,-1,1);
  TH1F *h_trk_sectorPhi_primary_allCuts_disk = new TH1F("h_trk_sectorPhi_primary_allCuts_disk","h_trk_sectorPhi_primary_allCuts_disk; Track #sectorPhi_{0} Distribution (#eta>0.95); Events / 0.02",100,-1,1);
  TH1F *h_trk_sectorPhi_np = new TH1F("h_trk_sectorPhi_np","h_trk_sectorPhi_np; Track #sectorPhi_{0} Distribution; Events / 0.02",100,-1,1);
  TH1F *h_trk_sectorPhi_np_noCuts = new TH1F("h_trk_sectorPhi_np_noCuts","h_trk_sectorPhi_np_noCuts; Track #sectorPhi_{0} Distribution; Events / 0.02",100,-1,1);
  TH1F *h_trk_sectorPhi_np_noCuts_barrel = new TH1F("h_trk_sectorPhi_np_noCuts_barrel","h_trk_sectorPhi_np_noCuts_barrel; Track #sectorPhi_{0} Distribution (#eta<=0.95); Events / 0.02",100,-1,1);
  TH1F *h_trk_sectorPhi_np_noCuts_disk = new TH1F("h_trk_sectorPhi_np_noCuts_disk","h_trk_sectorPhi_np_noCuts_disk; Track #sectorPhi_{0} Distribution (#eta>0.95); Events / 0.02",100,-1,1);
  TH1F *h_trk_sectorPhi_np_noCuts_H = new TH1F("h_trk_sectorPhi_np_noCuts_H","h_trk_sectorPhi_np_noCuts_H; Track #sectorPhi_{0} Distribution (p_{T}>10 GeV); Events / 0.02",100,-1,1);
  TH1F *h_trk_sectorPhi_np_noCuts_L = new TH1F("h_trk_sectorPhi_np_noCuts_L","h_trk_sectorPhi_np_noCuts_L; Track #sectorPhi_{0} Distribution (p_{T}<10 GeV); Events / 0.02",100,-1,1);
  TH1F *h_trk_sectorPhi_fake_noCuts = new TH1F("h_trk_sectorPhi_fake_noCuts","h_trk_sectorPhi_fake_noCuts; Track #sectorPhi_{0} Distribution; Events / 0.02",100,-1,1);
  TH1F *h_trk_sectorPhi_PU_noCuts = new TH1F("h_trk_sectorPhi_PU_noCuts","h_trk_sectorPhi_PU_noCuts; Track #sectorPhi_{0} Distribution; Events / 0.02",100,-1,1);
  TH1F *h_trk_sectorPhi_notHiggs_noCuts = new TH1F("h_trk_sectorPhi_notHiggs_noCuts","h_trk_sectorPhi_notHiggs_noCuts; Track #sectorPhi_{0} Distribution; Events / 0.02",100,-1,1);
  TH1F *h_trk_sectorPhi_np_qualCuts = new TH1F("h_trk_sectorPhi_np_qualCuts","h_trk_sectorPhi_np_qualCuts; Track #sectorPhi_{0} Distribution; Events / 0.02",100,-1,1);
  TH1F *h_trk_sectorPhi_fake_qualCuts = new TH1F("h_trk_sectorPhi_fake_qualCuts","h_trk_sectorPhi_fake_qualCuts; Track #sectorPhi_{0} Distribution; Events / 0.02",100,-1,1);
  TH1F *h_trk_sectorPhi_PU_qualCuts = new TH1F("h_trk_sectorPhi_PU_qualCuts","h_trk_sectorPhi_PU_qualCuts; Track #sectorPhi_{0} Distribution; Events / 0.02",100,-1,1);
  TH1F *h_trk_sectorPhi_notHiggs_qualCuts = new TH1F("h_trk_sectorPhi_notHiggs_qualCuts","h_trk_sectorPhi_notHiggs_qualCuts; Track #sectorPhi_{0} Distribution; Events / 0.02",100,-1,1);
  TH1F *h_trk_sectorPhi_np_allCuts = new TH1F("h_trk_sectorPhi_np_allCuts","h_trk_sectorPhi_np_allCuts; Track #sectorPhi_{0} Distribution; Events / 0.02",100,-1,1);
  TH1F *h_trk_sectorPhi_np_allCuts_barrel = new TH1F("h_trk_sectorPhi_np_allCuts_barrel","h_trk_sectorPhi_np_allCuts_barrel; Track #sectorPhi_{0} Distribution (#eta<=0.95); Events / 0.02",100,-1,1); 
  TH1F *h_trk_sectorPhi_np_allCuts_disk = new TH1F("h_trk_sectorPhi_np_allCuts_disk","h_trk_sectorPhi_np_allCuts_disk; Track #sectorPhi_{0} Distribution (#eta>0.95); Events / 0.02",100,-1,1); 

  TH1F *h_tp_pt = new TH1F("h_tp_pt", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
  TH1F *h_tp_eta = new TH1F("h_tp_eta", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
  TH1F *h_tp_d0 = new TH1F("h_tp_d0", ";Tracking particle d_{0} [cm]; Tracking particles / 0.01 cm", 50, -2, 2);
  TH1F *h_numSelectedTrks = new TH1F("h_numSelectedTrks","h_numSelectedTrks; Number of Selected Tracks; Events / 1.0",100,0,100);
  TH1F *h_numSelectedTrks_zoomOut = new TH1F("h_numSelectedTrks_zoomOut","h_numSelectedTrks_zoomOut; Number of Selected Tracks; Events / 10.0",100,0,1000);
  TH1F *h_trk_H_T = new TH1F("h_trk_H_T","h_trk_H_T; Event Track Scalar p_{T} Sum [GeV]; Events / 10.0",100,0,1000);
  TH1F *h_trk_MET = new TH1F("h_trk_MET","h_trk_MET; Event Track Missing E_{T} [GeV]; Events / 4.0",100,0,400);
  TH1F *h_trk_oneMatch_H_T = new TH1F("h_trk_oneMatch_H_T","h_trk_oneMatch_H_T; Event Track Scalar p_{T} Sum [GeV]; Events / 10.0",100,0,1000);
  TH1F *h_trk_oneMatch_MET = new TH1F("h_trk_oneMatch_MET","h_trk_oneMatch_MET; Event Track Missing E_{T} [GeV]; Events / 4.0",100,0,400);
  TH1F *h_trk_oneVert_H_T = new TH1F("h_trk_oneVert_H_T","h_trk_oneVert_H_T; Event Track Scalar p_{T} Sum [GeV]; Events / 10.0",100,0,1000);
  TH1F *h_trk_oneVert_MET = new TH1F("h_trk_oneVert_MET","h_trk_oneVert_MET; Event Track Missing E_{T} [GeV]; Events / 4.0",100,0,400);
  TH1F *h_tp_H_T = new TH1F("h_tp_H_T","h_tp_H_T; Event TP Scalar p_{T} Sum [GeV]; Events / 10.0",100,0,1000);
  TH1F *h_tp_MET = new TH1F("h_tp_MET","h_tp_MET; Event TP Missing E_{T} [GeV]; Events / 4.0",100,0,400);

  TH1F *h_trueVertex_numAllCuts = new TH1F("h_trueVertex_numAllCuts","h_trueVertex_numAllCuts; TP Vertices; Events / 1.0",40,0,40);
  TH1F *h_trueVertex_numNoCuts = new TH1F("h_trueVertex_numNoCuts","h_trueVertex_numNoCuts; TP Vertices; Events / 1.0",40,0,40);
  TH1F *h_trueVertex_numTPs = new TH1F("h_trueVertex_numTPs","h_trueVertex_numTPs; TPs Associated with Vertex; Events / 1.0",6,0,6);
  TH1F *h_trueVertex_numTracks = new TH1F("h_trueVertex_numTracks","h_trueVertex_numTracks; Tracks Associated with Vertex; Events / 1.0",6,0,6);
  TH1F *h_trueVertex_inTraj = new TH1F("h_trueVertex_inTraj","h_trueVertex_inTraj; Calc Vertex Return Code; Events / 1.0",6,0,6);
  TH1F *h_trueVertex_x = new TH1F("h_trueVertex_x","h_trueVertex_x; TP Vertex X Position (cm); Events / 0.1 cm",100,-5.0,5.0);
  TH1F *h_trueVertex_y = new TH1F("h_trueVertex_y","h_trueVertex_y; TP Vertex Y Position (cm); Events / 0.1 cm",100,-5.0,5.0);
  TH1F *h_trueVertex_z = new TH1F("h_trueVertex_z","h_trueVertex_z; TP Vertex Z Position (cm); Events / 0.025 cm",4000,-50,50);
  TH1F *h_trueVertex_sumPt = new TH1F("h_trueVertex_sumPt","h_trueVertex_sumPt; TP Vertex Momentum (GeV); Events / 1.0 GeV",200,0,200.0);
  TH1F *h_trueVertex_highPt = new TH1F("h_trueVertex_highPt","h_trueVertex_highPt; p_{T} of Leading p_{T} TP from vertex (GeV); Events / 1.0 GeV",200,0,200.0);
  TH1F *h_trueVertex_lowPt = new TH1F("h_trueVertex_lowPt","h_trueVertex_lowPt; p_{T} of Lower p_{T} TP from vertex (GeV); Events / 1.0 GeV",200,0,200.0);
  TH1F *h_trueVertexCuts_x = new TH1F("h_trueVertexCuts_x","h_trueVertexCuts_x; TP Vertex X Position (cm); Events / 0.1 cm",100,-5.0,5.0);
  TH1F *h_trueVertexCuts_y = new TH1F("h_trueVertexCuts_y","h_trueVertexCuts_y; TP Vertex Y Position (cm); Events / 0.1 cm",100,-5.0,5.0);
  TH1F *h_trueVertexCuts_z = new TH1F("h_trueVertexCuts_z","h_trueVertexCuts_z; TP Vertex Z Position (cm); Events / 0.025 cm",4000,-50,50);
  TH1F *h_trueVertexCuts_sumPt = new TH1F("h_trueVertexCuts_sumPt","h_trueVertexCuts_sumPt; TP Vertex Momentum (GeV); Events / 1.0 GeV",200,0,200.0);
  TH1F *h_trueVertexCuts_indexPt = new TH1F("h_trueVertexCuts_indexPt","h_trueVertexCuts_indexPt; Pt Ranking of TP; Events / 1.0",100,0,100.0);
  TH1F *h_trackVertex_numAllCuts = new TH1F("h_trackVertex_numAllCuts","h_trackVertex_numAllCuts; Track Vertices; Events / 1.0",40,0,40);
  TH1F *h_trackVertex_numNoCuts = new TH1F("h_trackVertex_numNoCuts","h_trackVertex_numNoCuts; Track Vertices; Events / 1.0",40,0,40);
  TH1F *h_trackVertex_numTracks = new TH1F("h_trackVertex_numTracks","h_trackVertex_numTracks; Tracks Associated with Vertex; Events / 1.0",20,0,20);
  TH1F *h_trackVertex_numVertPerTrack = new TH1F("h_trackVertex_numVertPerTrack","h_trackVertex_numVertPerTrack; Vertices Associated with Track; Events / 1.0",20,0,20);
  TH1F *h_trackVertex_inTraj = new TH1F("h_trackVertex_inTraj","h_trackVertex_inTraj; Calc Vertex Return Code; Events / 1.0",6,0,6);
  TH1F *h_trackVertex_x = new TH1F("h_trackVertex_x","h_trackVertex_x; Track Vertex X Position (cm); Events / 0.1 cm",100,-5.0,5.0);
  TH1F *h_trackVertex_y = new TH1F("h_trackVertex_y","h_trackVertex_y; Track Vertex Y Position (cm); Events / 0.1 cm",100,-5.0,5.0);
  TH1F *h_trackVertex_z = new TH1F("h_trackVertex_z","h_trackVertex_z; Track Vertex Z Position (cm); Events / 0.025 cm",4000,-50,50);
  TH1F *h_trackVertex_sumPt = new TH1F("h_trackVertex_sumPt","h_trackVertex_sumPt; Track Vertex Momentum (GeV); Events / 1.0 GeV",200,0,200.0);
  TH1F *h_trackVertexCuts_x = new TH1F("h_trackVertexCuts_x","h_trackVertexCuts_x; Track Vertex X Position (cm); Events / 0.1 cm",100,-5.0,5.0);
  TH1F *h_trackVertexCuts_y = new TH1F("h_trackVertexCuts_y","h_trackVertexCuts_y; Track Vertex Y Position (cm); Events / 0.1 cm",100,-5.0,5.0);
  TH1F *h_trackVertexCuts_z = new TH1F("h_trackVertexCuts_z","h_trackVertexCuts_z; Track Vertex Z Position (cm); Events / 0.025 cm",4000,-50,50);
  TH1F *h_trackVertexCuts_sumPt = new TH1F("h_trackVertexCuts_sumPt","h_trackVertexCuts_sumPt; Track Vertex Momentum (GeV); Events / 1.0 GeV",200,0,200.0);
  TH1F *h_trackVertexCuts_indexPt = new TH1F("h_trackVertexCuts_indexPt","h_trackVertexCuts_indexPt; Pt Ranking of Track; Events / 1.0",50,0,50.0);
  TH1F *h_correct_trackVertex_indexPt = new TH1F("h_correct_trackVertex_indexPt","h_correct_trackVertex_indexPt; Pt Ranking of Track; Events / 1.0",100,0,100.0);

  TH1F *h_trk_MVA1_primary_noCuts = new TH1F("h_trk_MVA1_primary_noCuts","h_trk_MVA1_primary_noCuts; Track MVA Score ; Events / 0.01",100,0,1);
  TH1F *h_trk_MVA1_primary_noCuts_H = new TH1F("h_trk_MVA1_primary_noCuts_H","h_trk_MVA1_primary_noCuts_H; Track MVA Score (p_{T}>10 GeV); Events / 0.01",100,0,1);
  TH1F *h_trk_MVA1_primary_noCuts_L = new TH1F("h_trk_MVA1_primary_noCuts_L","h_trk_MVA1_primary_noCuts_L; Track MVA Score (p_{T}<10 GeV); Events / 0.01",100,0,1);
  TH1F *h_trk_MVA1_primary_qualCuts = new TH1F("h_trk_MVA1_primary_qualCuts","h_trk_MVA1_primary_qualCuts; Track MVA Score ; Events / 0.01",100,0,1);
  TH1F *h_trk_MVA1_primary_allCuts = new TH1F("h_trk_MVA1_primary_allCuts","h_trk_MVA1_primary_allCuts; Track MVA Score ; Events / 0.01",100,0,1);
  TH1F *h_trk_MVA1_primary_allCuts_P = new TH1F("h_trk_MVA1_primary_allCuts_P","h_trk_MVA1_primary_allCuts_P; Track MVA Score (d_{0}<1 cm); Events / 0.01",100,0,1);
  TH1F *h_trk_MVA1_primary_allCuts_D = new TH1F("h_trk_MVA1_primary_allCuts_D","h_trk_MVA1_primary_allCuts_D; Track MVA Score (d_{0}>1 cm); Events / 0.01",100,0,1);
  TH1F *h_trk_MVA1_np_noCuts = new TH1F("h_trk_MVA1_np_noCuts","h_trk_MVA1_np_noCuts; Track MVA Score ; Events / 0.01",100,0,1);
  TH1F *h_trk_MVA1_np_noCuts_H = new TH1F("h_trk_MVA1_np_noCuts_H","h_trk_MVA1_np_noCuts_H; Track MVA Score (p_{T}>10 GeV); Events / 0.01",100,0,1);
  TH1F *h_trk_MVA1_np_noCuts_L = new TH1F("h_trk_MVA1_np_noCuts_L","h_trk_MVA1_np_noCuts_L; Track MVA Score (p_{T}<10 GeV); Events / 0.01",100,0,1);
  TH1F *h_trk_MVA1_fake_noCuts = new TH1F("h_trk_MVA1_fake_noCuts","h_trk_MVA1_fake_noCuts; Track MVA Score ; Events / 0.01",100,0,1);
  TH1F *h_trk_MVA1_PU_noCuts = new TH1F("h_trk_MVA1_PU_noCuts","h_trk_MVA1_PU_noCuts; Track MVA Score ; Events / 0.01",100,0,1);
  TH1F *h_trk_MVA1_notHiggs_noCuts = new TH1F("h_trk_MVA1_notHiggs_noCuts","h_trk_MVA1_notHiggs_noCuts; Track MVA Score ; Events / 0.01",100,0,1);
  TH1F *h_trk_MVA1_np_qualCuts = new TH1F("h_trk_MVA1_np_qualCuts","h_trk_MVA1_np_qualCuts; Track MVA Score ; Events / 0.01",100,0,1);
  TH1F *h_trk_MVA1_fake_qualCuts = new TH1F("h_trk_MVA1_fake_qualCuts","h_trk_MVA1_fake_qualCuts; Track MVA Score ; Events / 0.01",100,0,1);
  TH1F *h_trk_MVA1_PU_qualCuts = new TH1F("h_trk_MVA1_PU_qualCuts","h_trk_MVA1_PU_qualCuts; Track MVA Score ; Events / 0.01",100,0,1);
  TH1F *h_trk_MVA1_notHiggs_qualCuts = new TH1F("h_trk_MVA1_notHiggs_qualCuts","h_trk_MVA1_notHiggs_qualCuts; Track MVA Score ; Events / 0.01",100,0,1);
  TH1F *h_trk_MVA1_np_allCuts = new TH1F("h_trk_MVA1_np_allCuts","h_trk_MVA1_np_allCuts; Track MVA Score ; Events / 0.01",100,0,1);
  TH1F *h_trk_MVA1_np_allCuts_P = new TH1F("h_trk_MVA1_np_allCuts_P","h_trk_MVA1_np_allCuts_P; Track MVA Score (d_{0}<1 cm); Events / 0.01",100,0,1);
  TH1F *h_trk_MVA1_np_allCuts_D = new TH1F("h_trk_MVA1_np_allCuts_D","h_trk_MVA1_np_allCuts_D; Track MVA Score (d_{0}>1 cm); Events / 0.01",100,0,1);
  TH1F *h_trk_MVA2_primary_noCuts = new TH1F("h_trk_MVA2_primary_noCuts","h_trk_MVA2_primary_noCuts; Track MVA Score ; Events / 0.01",100,0,1);
  TH1F *h_trk_MVA2_primary_allCuts = new TH1F("h_trk_MVA2_primary_allCuts","h_trk_MVA2_primary_allCuts; Track MVA Score ; Events / 0.01",100,0,1);
  TH1F *h_trk_MVA2_primary_allCuts_P = new TH1F("h_trk_MVA2_primary_allCuts_P","h_trk_MVA2_primary_allCuts_P; Track MVA Score (d_{0}<1 cm); Events / 0.01",100,0,1);
  TH1F *h_trk_MVA2_primary_allCuts_D = new TH1F("h_trk_MVA2_primary_allCuts_D","h_trk_MVA2_primary_allCuts_D; Track MVA Score (d_{0}>1 cm); Events / 0.01",100,0,1);
  TH1F *h_trk_MVA2_np_noCuts = new TH1F("h_trk_MVA2_np_noCuts","h_trk_MVA2_np_noCuts; Track MVA Score ; Events / 0.01",100,0,1);
  TH1F *h_trk_MVA2_np_allCuts = new TH1F("h_trk_MVA2_np_allCuts","h_trk_MVA2_np_allCuts; Track MVA Score ; Events / 0.01",100,0,1);
  TH1F *h_trk_MVA2_np_allCuts_P = new TH1F("h_trk_MVA2_np_allCuts_P","h_trk_MVA2_np_allCuts_P; Track MVA Score (d_{0}<1 cm); Events / 0.01",100,0,1);
  TH1F *h_trk_MVA2_np_allCuts_D = new TH1F("h_trk_MVA2_np_allCuts_D","h_trk_MVA2_np_allCuts_D; Track MVA Score (d_{0}>1 cm); Events / 0.01",100,0,1);
  // Chi2 plots
  TH1F *h_trk_chi2rphidof = new TH1F("h_trk_chi2rphidof","h_trk_chi2rphidof; Track #chi^{2}_{r#phi}/d.o.f ; Events / 0.1",50,0,5);
  TH1F *h_trk_chi2rphidof_primary = new TH1F("h_trk_chi2rphidof_primary","h_trk_chi2rphidof_primary; Track #chi^{2}_{r#phi}/d.o.f ; Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rphidof_primary_noCuts = new TH1F("h_trk_chi2rphidof_primary_noCuts","h_trk_chi2rphidof_primary_noCuts; Track #chi^{2}_{r#phi}/d.o.f ; Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rphidof_primary_noCuts_zoomOut = new TH1F("h_trk_chi2rphidof_primary_noCuts_zoomOut","h_trk_chi2rphidof_primary_noCuts_zoomOut; Track #chi^{2}_{r#phi}/d.o.f ; Events / 0.6",100,0,60);
  TH1F *h_trk_chi2rphidof_primary_noCuts_barrel = new TH1F("h_trk_chi2rphidof_primary_noCuts_barrel","h_trk_chi2rphidof_primary_noCuts_barrel; Track #chi^{2}_{r#phi}/d.o.f (#eta<=0.95); Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rphidof_primary_noCuts_disk = new TH1F("h_trk_chi2rphidof_primary_noCuts_disk","h_trk_chi2rphidof_primary_noCuts_disk; Track #chi^{2}_{r#phi}/d.o.f (#eta>0.95); Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rphidof_primary_noCuts_H = new TH1F("h_trk_chi2rphidof_primary_noCuts_H","h_trk_chi2rphidof_primary_noCuts_H; Track #chi^{2}_{r#phi}/d.o.f (p_{T}>10 GeV); Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rphidof_primary_noCuts_L = new TH1F("h_trk_chi2rphidof_primary_noCuts_L","h_trk_chi2rphidof_primary_noCuts_L; Track #chi^{2}_{r#phi}/d.o.f (p_{T}<10 GeV); Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rphidof_primary_qualCuts = new TH1F("h_trk_chi2rphidof_primary_qualCuts","h_trk_chi2rphidof_primary_qualCuts; Track #chi^{2}_{r#phi}/d.o.f ; Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rphidof_primary_allCuts = new TH1F("h_trk_chi2rphidof_primary_allCuts","h_trk_chi2rphidof_primary_allCuts; Track #chi^{2}_{r#phi}/d.o.f ; Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rphidof_primary_allCuts_zoomOut = new TH1F("h_trk_chi2rphidof_primary_allCuts_zoomOut","h_trk_chi2rphidof_primary_allCuts_zoomOut; Track #chi^{2}_{r#phi}/d.o.f ; Events / 0.6",100,0,60);
  TH1F *h_trk_chi2rphidof_primary_allCuts_barrel = new TH1F("h_trk_chi2rphidof_primary_allCuts_barrel","h_trk_chi2rphidof_primary_allCuts_barrel; Track #chi^{2}_{r#phi}/d.o.f (#eta<=0.95); Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rphidof_primary_allCuts_disk = new TH1F("h_trk_chi2rphidof_primary_allCuts_disk","h_trk_chi2rphidof_primary_allCuts_disk; Track #chi^{2}_{r#phi}/d.o.f (#eta>0.95); Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rphidof_np = new TH1F("h_trk_chi2rphidof_np","h_trk_chi2rphidof_np; Track #chi^{2}_{r#phi}/d.o.f ; Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rphidof_np_noCuts = new TH1F("h_trk_chi2rphidof_np_noCuts","h_trk_chi2rphidof_np_noCuts; Track #chi^{2}_{r#phi}/d.o.f ; Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rphidof_np_noCuts_zoomOut = new TH1F("h_trk_chi2rphidof_np_noCuts_zoomOut","h_trk_chi2rphidof_np_noCuts_zoomOut; Track #chi^{2}_{r#phi}/d.o.f ; Events / 0.6",100,0,60);
  TH1F *h_trk_chi2rphidof_np_noCuts_barrel = new TH1F("h_trk_chi2rphidof_np_noCuts_barrel","h_trk_chi2rphidof_np_noCuts_barrel; Track #chi^{2}_{r#phi}/d.o.f (#eta<=0.95); Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rphidof_np_noCuts_disk = new TH1F("h_trk_chi2rphidof_np_noCuts_disk","h_trk_chi2rphidof_np_noCuts_disk; Track #chi^{2}_{r#phi}/d.o.f (#eta>0.95); Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rphidof_np_noCuts_H = new TH1F("h_trk_chi2rphidof_np_noCuts_H","h_trk_chi2rphidof_np_noCuts_H; Track #chi^{2}_{r#phi}/d.o.f (p_{T}>10 GeV); Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rphidof_np_noCuts_L = new TH1F("h_trk_chi2rphidof_np_noCuts_L","h_trk_chi2rphidof_np_noCuts_L; Track #chi^{2}_{r#phi}/d.o.f (p_{T}<10 GeV); Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rphidof_fake_noCuts = new TH1F("h_trk_chi2rphidof_fake_noCuts","h_trk_chi2rphidof_fake_noCuts; Track #chi^{2}_{r#phi}/d.o.f ; Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rphidof_PU_noCuts = new TH1F("h_trk_chi2rphidof_PU_noCuts","h_trk_chi2rphidof_PU_noCuts; Track #chi^{2}_{r#phi}/d.o.f ; Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rphidof_notHiggs_noCuts = new TH1F("h_trk_chi2rphidof_notHiggs_noCuts","h_trk_chi2rphidof_notHiggs_noCuts; Track #chi^{2}_{r#phi}/d.o.f ; Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rphidof_np_qualCuts = new TH1F("h_trk_chi2rphidof_np_qualCuts","h_trk_chi2rphidof_np_qualCuts; Track #chi^{2}_{r#phi}/d.o.f ; Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rphidof_fake_qualCuts = new TH1F("h_trk_chi2rphidof_fake_qualCuts","h_trk_chi2rphidof_fake_qualCuts; Track #chi^{2}_{r#phi}/d.o.f ; Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rphidof_PU_qualCuts = new TH1F("h_trk_chi2rphidof_PU_qualCuts","h_trk_chi2rphidof_PU_qualCuts; Track #chi^{2}_{r#phi}/d.o.f ; Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rphidof_notHiggs_qualCuts = new TH1F("h_trk_chi2rphidof_notHiggs_qualCuts","h_trk_chi2rphidof_notHiggs_qualCuts; Track #chi^{2}_{r#phi}/d.o.f ; Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rphidof_np_allCuts = new TH1F("h_trk_chi2rphidof_np_allCuts","h_trk_chi2rphidof_np_allCuts; Track #chi^{2}_{r#phi}/d.o.f ; Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rphidof_np_allCuts_zoomOut = new TH1F("h_trk_chi2rphidof_np_allCuts_zoomOut","h_trk_chi2rphidof_np_allCuts_zoomOut; Track #chi^{2}_{r#phi}/d.o.f ; Events / 0.6",100,0,60);
  TH1F *h_trk_chi2rphidof_np_allCuts_barrel = new TH1F("h_trk_chi2rphidof_np_allCuts_barrel","h_trk_chi2rphidof_np_allCuts_barrel; Track #chi^{2}_{r#phi}/d.o.f (#eta<=0.95); Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rphidof_np_allCuts_disk = new TH1F("h_trk_chi2rphidof_np_allCuts_disk","h_trk_chi2rphidof_np_allCuts_disk; Track #chi^{2}_{r#phi}/d.o.f (#eta>0.95); Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rphidof_H = new TH1F("h_trk_chi2rphidof_H","h_trk_chi2rphidof_H; Track #chi^{2}_{r#phi}/d.o.f (p_{T}>8 GeV); Events / 0.1",50,0,5);   
  TH1F *h_trk_chi2rphidof_L = new TH1F("h_trk_chi2rphidof_L","h_trk_chi2rphidof_L; Track #chi^{2}_{r#phi}/d.o.f (p_{T}<8 GeV); Events / 0.1",50,0,5);   
  TH1F *h_trk_chi2rphidof_C = new TH1F("h_trk_chi2rphidof_C","h_trk_chi2rphidof_C; Track #chi^{2}_{r#phi}/d.o.f (|#eta|<0.8); Events / 0.1",50,0,5);   
  TH1F *h_trk_chi2rphidof_I = new TH1F("h_trk_chi2rphidof_I","h_trk_chi2rphidof_I; Track #chi^{2}_{r#phi}/d.o.f (0.8#leq|#eta|#leq1.6); Events / 0.1",50,0,5);   
  TH1F *h_trk_chi2rphidof_F = new TH1F("h_trk_chi2rphidof_F","h_trk_chi2rphidof_F; Track #chi^{2}_{r#phi}/d.o.f (|#eta|>1.6); Events / 0.1",50,0,5);   
  TH1F *h_trk_chi2rphidof_P = new TH1F("h_trk_chi2rphidof_P","h_trk_chi2rphidof_P; Track #chi^{2}_{r#phi}/d.o.f (d_{0}#leq 1 cm); Events / 0.1",50,0,5);
  TH1F *h_trk_chi2rphidof_D = new TH1F("h_trk_chi2rphidof_D","h_trk_chi2rphidof_D; Track #chi^{2}_{r#phi}/d.o.f (d_{0}>1 cm); Events / 0.1",50,0,5);
  TH1F *h_trk_chi2rzdof = new TH1F("h_trk_chi2rzdof","h_trk_chi2rzdof; Track #chi^{2}_{rz}/d.o.f ; Events / 0.1",50,0,5);
  TH1F *h_trk_chi2rzdof_primary = new TH1F("h_trk_chi2rzdof_primary","h_trk_chi2rzdof_primary; Track #chi^{2}_{rz}/d.o.f ; Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rzdof_primary_noCuts = new TH1F("h_trk_chi2rzdof_primary_noCuts","h_trk_chi2rzdof_primary_noCuts; Track #chi^{2}_{rz}/d.o.f ; Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rzdof_primary_noCuts_zoomOut = new TH1F("h_trk_chi2rzdof_primary_noCuts_zoomOut","h_trk_chi2rzdof_primary_noCuts_zoomOut; Track #chi^{2}_{rz}/d.o.f ; Events / 0.3",100,0,30);
  TH1F *h_trk_chi2rzdof_primary_noCuts_barrel = new TH1F("h_trk_chi2rzdof_primary_noCuts_barrel","h_trk_chi2rzdof_primary_noCuts_barrel; Track #chi^{2}_{rz}/d.o.f (#eta<=0.95); Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rzdof_primary_noCuts_disk = new TH1F("h_trk_chi2rzdof_primary_noCuts_disk","h_trk_chi2rzdof_primary_noCuts_disk; Track #chi^{2}_{rz}/d.o.f (#eta>0.95); Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rzdof_primary_noCuts_H = new TH1F("h_trk_chi2rzdof_primary_noCuts_H","h_trk_chi2rzdof_primary_noCuts_H; Track #chi^{2}_{rz}/d.o.f (p_{T}>10 GeV); Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rzdof_primary_noCuts_L = new TH1F("h_trk_chi2rzdof_primary_noCuts_L","h_trk_chi2rzdof_primary_noCuts_L; Track #chi^{2}_{rz}/d.o.f (p_{T}<10 GeV); Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rzdof_primary_qualCuts = new TH1F("h_trk_chi2rzdof_primary_qualCuts","h_trk_chi2rzdof_primary_qualCuts; Track #chi^{2}_{rz}/d.o.f ; Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rzdof_primary_allCuts = new TH1F("h_trk_chi2rzdof_primary_allCuts","h_trk_chi2rzdof_primary_allCuts; Track #chi^{2}_{rz}/d.o.f ; Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rzdof_primary_allCuts_zoomOut = new TH1F("h_trk_chi2rzdof_primary_allCuts_zoomOut","h_trk_chi2rzdof_primary_allCuts_zoomOut; Track #chi^{2}_{rz}/d.o.f ; Events / 0.3",100,0,30);
  TH1F *h_trk_chi2rzdof_primary_allCuts_barrel = new TH1F("h_trk_chi2rzdof_primary_allCuts_barrel","h_trk_chi2rzdof_primary_allCuts_barrel; Track #chi^{2}_{rz}/d.o.f (#eta<=0.95); Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rzdof_primary_allCuts_disk = new TH1F("h_trk_chi2rzdof_primary_allCuts_disk","h_trk_chi2rzdof_primary_allCuts_disk; Track #chi^{2}_{rz}/d.o.f (#eta>0.95); Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rzdof_np = new TH1F("h_trk_chi2rzdof_np","h_trk_chi2rzdof_np; Track #chi^{2}_{rz}/d.o.f ; Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rzdof_np_noCuts = new TH1F("h_trk_chi2rzdof_np_noCuts","h_trk_chi2rzdof_np_noCuts; Track #chi^{2}_{rz}/d.o.f ; Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rzdof_np_noCuts_zoomOut = new TH1F("h_trk_chi2rzdof_np_noCuts_zoomOut","h_trk_chi2rzdof_np_noCuts_zoomOut; Track #chi^{2}_{rz}/d.o.f ; Events / 0.3",100,0,30);
  TH1F *h_trk_chi2rzdof_np_noCuts_barrel = new TH1F("h_trk_chi2rzdof_np_noCuts_barrel","h_trk_chi2rzdof_np_noCuts_barrel; Track #chi^{2}_{rz}/d.o.f (#eta<=0.95); Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rzdof_np_noCuts_disk = new TH1F("h_trk_chi2rzdof_np_noCuts_disk","h_trk_chi2rzdof_np_noCuts_disk; Track #chi^{2}_{rz}/d.o.f (#eta>0.95); Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rzdof_np_noCuts_H = new TH1F("h_trk_chi2rzdof_np_noCuts_H","h_trk_chi2rzdof_np_noCuts_H; Track #chi^{2}_{rz}/d.o.f (p_{T}>10 GeV); Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rzdof_np_noCuts_L = new TH1F("h_trk_chi2rzdof_np_noCuts_L","h_trk_chi2rzdof_np_noCuts_L; Track #chi^{2}_{rz}/d.o.f (p_{T}<10 GeV); Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rzdof_fake_noCuts = new TH1F("h_trk_chi2rzdof_fake_noCuts","h_trk_chi2rzdof_fake_noCuts; Track #chi^{2}_{rz}/d.o.f ; Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rzdof_PU_noCuts = new TH1F("h_trk_chi2rzdof_PU_noCuts","h_trk_chi2rzdof_PU_noCuts; Track #chi^{2}_{rz}/d.o.f ; Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rzdof_notHiggs_noCuts = new TH1F("h_trk_chi2rzdof_notHiggs_noCuts","h_trk_chi2rzdof_notHiggs_noCuts; Track #chi^{2}_{rz}/d.o.f ; Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rzdof_np_qualCuts = new TH1F("h_trk_chi2rzdof_np_qualCuts","h_trk_chi2rzdof_np_qualCuts; Track #chi^{2}_{rz}/d.o.f ; Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rzdof_fake_qualCuts = new TH1F("h_trk_chi2rzdof_fake_qualCuts","h_trk_chi2rzdof_fake_qualCuts; Track #chi^{2}_{rz}/d.o.f ; Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rzdof_PU_qualCuts = new TH1F("h_trk_chi2rzdof_PU_qualCuts","h_trk_chi2rzdof_PU_qualCuts; Track #chi^{2}_{rz}/d.o.f ; Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rzdof_notHiggs_qualCuts = new TH1F("h_trk_chi2rzdof_notHiggs_qualCuts","h_trk_chi2rzdof_notHiggs_qualCuts; Track #chi^{2}_{rz}/d.o.f ; Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rzdof_np_allCuts = new TH1F("h_trk_chi2rzdof_np_allCuts","h_trk_chi2rzdof_np_allCuts; Track #chi^{2}_{rz}/d.o.f ; Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rzdof_np_allCuts_zoomOut = new TH1F("h_trk_chi2rzdof_np_allCuts_zoomOut","h_trk_chi2rzdof_np_allCuts_zoomOut; Track #chi^{2}_{rz}/d.o.f ; Events / 0.3",100,0,30);
  TH1F *h_trk_chi2rzdof_np_allCuts_barrel = new TH1F("h_trk_chi2rzdof_np_allCuts_barrel","h_trk_chi2rzdof_np_allCuts_barrel; Track #chi^{2}_{rz}/d.o.f (#eta<=0.95); Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rzdof_np_allCuts_disk = new TH1F("h_trk_chi2rzdof_np_allCuts_disk","h_trk_chi2rzdof_np_allCuts_disk; Track #chi^{2}_{rz}/d.o.f (#eta>0.95); Events / 0.06",100,0,6);
  TH1F *h_trk_chi2rzdof_H = new TH1F("h_trk_chi2rzdof_H","h_trk_chi2rzdof_H; Track #chi^{2}_{rz}/d.o.f (p_{T}>8 GeV); Events / 0.1",50,0,5);
  TH1F *h_trk_chi2rzdof_L = new TH1F("h_trk_chi2rzdof_L","h_trk_chi2rzdof_L; Track #chi^{2}_{rz}/d.o.f (p_{T}<8 GeV); Events / 0.1",50,0,5);
  TH1F *h_trk_chi2rzdof_C = new TH1F("h_trk_chi2rzdof_C","h_trk_chi2rzdof_C; Track #chi^{2}_{rz}/d.o.f (|#eta|<0.8); Events / 0.1",50,0,5);
  TH1F *h_trk_chi2rzdof_I = new TH1F("h_trk_chi2rzdof_I","h_trk_chi2rzdof_I; Track #chi^{2}_{rz}/d.o.f (0.8#leq|#eta|#leq1.6); Events / 0.1",50,0,5);
  TH1F *h_trk_chi2rzdof_F = new TH1F("h_trk_chi2rzdof_F","h_trk_chi2rzdof_F; Track #chi^{2}_{rz}/d.o.f (|#eta|>1.6); Events / 0.1",50,0,5);
  TH1F *h_trk_chi2rzdof_P = new TH1F("h_trk_chi2rzdof_P","h_trk_chi2rzdof_P; Track #chi^{2}_{rz}/d.o.f (d_{0}#leq 1 cm); Events / 0.1",50,0,5);
  TH1F *h_trk_chi2rzdof_D = new TH1F("h_trk_chi2rzdof_D","h_trk_chi2rzdof_D; Track #chi^{2}_{rz}/d.o.f (d_{0}>1 cm); Events / 0.1",50,0,5);
  TH1F *h_trk_bendchi2 = new TH1F("h_trk_bendchi2","h_trk_bendchi2; Track bend #chi^{2} ; Events / 0.1",50,0,5);
  TH1F *h_trk_bendchi2_primary = new TH1F("h_trk_bendchi2_primary","h_trk_bendchi2_primary; Track bend #chi^{2} ; Events / 0.1",100,0,10);
  TH1F *h_trk_bendchi2_primary_noCuts = new TH1F("h_trk_bendchi2_primary_noCuts","h_trk_bendchi2_primary_noCuts; Track bend #chi^{2} ; Events / 0.1",100,0,10);
  TH1F *h_trk_bendchi2_primary_noCuts_zoomOut = new TH1F("h_trk_bendchi2_primary_noCuts_zoomOut","h_trk_bendchi2_primary_noCuts_zoomOut; Track bend #chi^{2} ; Events / 0.7",100,0,70);
  TH1F *h_trk_bendchi2_primary_noCuts_barrel = new TH1F("h_trk_bendchi2_primary_noCuts_barrel","h_trk_bendchi2_primary_noCuts_barrel; Track bend #chi^{2} (#eta<=0.95); Events / 0.1",100,0,10);
  TH1F *h_trk_bendchi2_primary_noCuts_disk = new TH1F("h_trk_bendchi2_primary_noCuts_disk","h_trk_bendchi2_primary_noCuts_disk; Track bend #chi^{2} (#eta>0.95); Events / 0.1",100,0,10);
  TH1F *h_trk_bendchi2_primary_noCuts_H = new TH1F("h_trk_bendchi2_primary_noCuts_H","h_trk_bendchi2_primary_noCuts_H; Track bend #chi^{2} (p_{T}>10 GeV); Events / 0.1",100,0,10);
  TH1F *h_trk_bendchi2_primary_noCuts_L = new TH1F("h_trk_bendchi2_primary_noCuts_L","h_trk_bendchi2_primary_noCuts_L; Track bend #chi^{2} (p_{T}<10 GeV); Events / 0.1",100,0,10);
  TH1F *h_trk_bendchi2_primary_qualCuts = new TH1F("h_trk_bendchi2_primary_qualCuts","h_trk_bendchi2_primary_qualCuts; Track bend #chi^{2} ; Events / 0.1",100,0,10);
  TH1F *h_trk_bendchi2_primary_allCuts = new TH1F("h_trk_bendchi2_primary_allCuts","h_trk_bendchi2_primary_allCuts; Track bend #chi^{2} ; Events / 0.1",100,0,10);
  TH1F *h_trk_bendchi2_primary_allCuts_zoomOut = new TH1F("h_trk_bendchi2_primary_allCuts_zoomOut","h_trk_bendchi2_primary_allCuts_zoomOut; Track bend #chi^{2} ; Events / 0.7",100,0,70);
  TH1F *h_trk_bendchi2_primary_allCuts_barrel = new TH1F("h_trk_bendchi2_primary_allCuts_barrel","h_trk_bendchi2_primary_allCuts_barrel; Track bend #chi^{2} (#eta<=0.95); Events / 0.1",100,0,10);
  TH1F *h_trk_bendchi2_primary_allCuts_disk = new TH1F("h_trk_bendchi2_primary_allCuts_disk","h_trk_bendchi2_primary_allCuts_disk; Track bend #chi^{2} (#eta>0.95); Events / 0.1",100,0,10);
  TH1F *h_trk_bendchi2_np = new TH1F("h_trk_bendchi2_np","h_trk_bendchi2_np; Track bend #chi^{2} ; Events / 0.1",100,0,10);
  TH1F *h_trk_bendchi2_np_noCuts = new TH1F("h_trk_bendchi2_np_noCuts","h_trk_bendchi2_np_noCuts; Track bend #chi^{2} ; Events / 0.1",100,0,10);
  TH1F *h_trk_bendchi2_np_noCuts_zoomOut = new TH1F("h_trk_bendchi2_np_noCuts_zoomOut","h_trk_bendchi2_np_noCuts_zoomOut; Track bend #chi^{2} ; Events / 0.7",100,0,70);
  TH1F *h_trk_bendchi2_np_noCuts_barrel = new TH1F("h_trk_bendchi2_np_noCuts_barrel","h_trk_bendchi2_np_noCuts_barrel; Track bend #chi^{2} (#eta<=0.95); Events / 0.1",100,0,10);
  TH1F *h_trk_bendchi2_np_noCuts_disk = new TH1F("h_trk_bendchi2_np_noCuts_disk","h_trk_bendchi2_np_noCuts_disk; Track bend #chi^{2} (#eta>0.95); Events / 0.1",100,0,10);
  TH1F *h_trk_bendchi2_np_noCuts_H = new TH1F("h_trk_bendchi2_np_noCuts_H","h_trk_bendchi2_np_noCuts_H; Track bend #chi^{2} (p_{T}>10 GeV); Events / 0.1",100,0,10);
  TH1F *h_trk_bendchi2_np_noCuts_L = new TH1F("h_trk_bendchi2_np_noCuts_L","h_trk_bendchi2_np_noCuts_L; Track bend #chi^{2} (p_{T}<10 GeV); Events / 0.1",100,0,10);
  TH1F *h_trk_bendchi2_fake_noCuts = new TH1F("h_trk_bendchi2_fake_noCuts","h_trk_bendchi2_fake_noCuts; Track bend #chi^{2} ; Events / 0.1",100,0,10);
  TH1F *h_trk_bendchi2_PU_noCuts = new TH1F("h_trk_bendchi2_PU_noCuts","h_trk_bendchi2_PU_noCuts; Track bend #chi^{2} ; Events / 0.1",100,0,10);
  TH1F *h_trk_bendchi2_notHiggs_noCuts = new TH1F("h_trk_bendchi2_notHiggs_noCuts","h_trk_bendchi2_notHiggs_noCuts; Track bend #chi^{2} ; Events / 0.1",100,0,10);
  TH1F *h_trk_bendchi2_np_qualCuts = new TH1F("h_trk_bendchi2_np_qualCuts","h_trk_bendchi2_np_qualCuts; Track bend #chi^{2} ; Events / 0.1",100,0,10);
  TH1F *h_trk_bendchi2_fake_qualCuts = new TH1F("h_trk_bendchi2_fake_qualCuts","h_trk_bendchi2_fake_qualCuts; Track bend #chi^{2} ; Events / 0.1",100,0,10);
  TH1F *h_trk_bendchi2_PU_qualCuts = new TH1F("h_trk_bendchi2_PU_qualCuts","h_trk_bendchi2_PU_qualCuts; Track bend #chi^{2} ; Events / 0.1",100,0,10);
  TH1F *h_trk_bendchi2_notHiggs_qualCuts = new TH1F("h_trk_bendchi2_notHiggs_qualCuts","h_trk_bendchi2_notHiggs_qualCuts; Track bend #chi^{2} ; Events / 0.1",100,0,10);
  TH1F *h_trk_bendchi2_np_allCuts = new TH1F("h_trk_bendchi2_np_allCuts","h_trk_bendchi2_np_allCuts; Track bend #chi^{2} ; Events / 0.1",100,0,10);
  TH1F *h_trk_bendchi2_np_allCuts_zoomOut = new TH1F("h_trk_bendchi2_np_allCuts_zoomOut","h_trk_bendchi2_np_allCuts_zoomOut; Track bend #chi^{2} ; Events / 0.7",100,0,70);
  TH1F *h_trk_bendchi2_np_allCuts_barrel = new TH1F("h_trk_bendchi2_np_allCuts_barrel","h_trk_bendchi2_np_allCuts_barrel; Track bend #chi^{2} (#eta<=0.95); Events / 0.1",100,0,10);
  TH1F *h_trk_bendchi2_np_allCuts_disk = new TH1F("h_trk_bendchi2_np_allCuts_disk","h_trk_bendchi2_np_allCuts_disk; Track bend #chi^{2} (#eta>0.95); Events / 0.1",100,0,10);
  TH1F *h_trk_bendchi2_H = new TH1F("h_trk_bendchi2_H","h_trk_bendchi2_H; Track bend #chi^{2} (p_{T}>8 GeV); Events / 0.1",50,0,5);
  TH1F *h_trk_bendchi2_L = new TH1F("h_trk_bendchi2_L","h_trk_bendchi2_L; Track bend #chi^{2} (p_{T}<8 GeV); Events / 0.1",50,0,5);
  TH1F *h_trk_bendchi2_C = new TH1F("h_trk_bendchi2_C","h_trk_bendchi2_C; Track bend #chi^{2} (|#eta|<0.8); Events / 0.1",50,0,5);
  TH1F *h_trk_bendchi2_I = new TH1F("h_trk_bendchi2_I","h_trk_bendchi2_I; Track bend #chi^{2} (0.8#leq|#eta|#leq1.6); Events / 0.1",50,0,5);
  TH1F *h_trk_bendchi2_F = new TH1F("h_trk_bendchi2_F","h_trk_bendchi2_F; Track bend #chi^{2} (|#eta|>1.6); Events / 0.1",50,0,5);
  TH1F *h_trk_bendchi2_P = new TH1F("h_trk_bendchi2_P","h_trk_bendchi2_P; Track bend #chi^{2} (d_{0}#leq 1 cm); Events / 0.1",50,0,5);
  TH1F *h_trk_bendchi2_D = new TH1F("h_trk_bendchi2_D","h_trk_bendchi2_D; Track bend #chi^{2} (d_{0}>1 cm); Events / 0.1",50,0,5);

  // Efficiency of Identifying Tracks Plots
  TH1F* h_tp_pt_noCuts = new TH1F("h_tp_pt_noCuts", ";Tracking particle p_{T} [GeV]; Tracking particles / 0.4 GeV", 100, 0, 100.0);
  TH1F* h_tp_pt_noCuts_primary = new TH1F("h_tp_pt_noCuts_primary", ";Tracking particle p_{T} [GeV]; Tracking particles / 0.4 GeV", 100, 0, 100.0);
  TH1F* h_tp_pt_noCuts_PU = new TH1F("h_tp_pt_noCuts_PU", ";Tracking particle p_{T} [GeV]; Tracking particles / 0.4 GeV", 100, 0, 100.0);
  TH1F* h_match_tp_pt_noCuts = new TH1F("h_match_tp_pt_noCuts", ";Tracking particle p_{T} [GeV]; Tracking particles / 0.4 GeV", 100, 0, 100.0);
  TH1F* h_match_tp_pt_noCuts_primary = new TH1F("h_match_tp_pt_noCuts_primary", ";Tracking particle p_{T} [GeV]; Tracking particles / 0.4 GeV", 100, 0, 100.0);
  TH1F* h_tp_eta_noCuts = new TH1F("h_tp_eta_noCuts", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
  TH1F* h_tp_eta_noCuts_primary = new TH1F("h_tp_eta_noCuts_primary", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
  TH1F* h_tp_eta_noCuts_PU = new TH1F("h_tp_eta_noCuts_PU", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
  TH1F* h_match_tp_eta_noCuts = new TH1F("h_match_tp_eta_noCuts", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
  TH1F* h_tp_d0_noCuts = new TH1F("h_tp_d0_noCuts", ";Tracking particle d_{0} [cm]; Tracking particles / 0.0005 cm", 120, -0.03, 0.03);
  TH1F* h_tp_d0_noCuts_primary = new TH1F("h_tp_d0_noCuts_primary", ";Tracking particle d_{0} [cm]; Tracking particles / 0.0005 cm", 120, -0.03, 0.03);
  TH1F* h_tp_d0_noCuts_PU = new TH1F("h_tp_d0_noCuts_PU", ";Tracking particle d_{0} [cm]; Tracking particles / 0.0005 cm", 120, -0.03, 0.03);
  TH1F* h_match_tp_d0_noCuts = new TH1F("h_match_tp_d0_noCuts", ";Tracking particle d_{0} [cm]; Tracking particles / 0.0005 cm", 120, -0.03, 0.03);
  TH1F* h_tp_z0_noCuts = new TH1F("h_tp_z0_noCuts", ";Tracking particle z_{0} [cm]; Tracking particles / 0.1 cm", 220, -2, 20);
  TH1F* h_tp_z0_noCuts_primary = new TH1F("h_tp_z0_noCuts_primary", ";Tracking particle z_{0} [cm]; Tracking particles / 0.1 cm", 220, -2, 20);
  TH1F* h_tp_z0_noCuts_PU = new TH1F("h_tp_z0_noCuts_PU", ";Tracking particle z_{0} [cm]; Tracking particles / 0.1 cm", 220, -2, 20);

  TH1F* h_tp_pt_maxD0Cut = new TH1F("h_tp_pt_maxD0Cut", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
  TH1F* h_match_tp_pt_maxD0Cut = new TH1F("h_match_tp_pt_maxD0Cut", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
  TH1F* h_tp_eta_maxD0Cut = new TH1F("h_tp_eta_maxD0Cut", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
  TH1F* h_match_tp_eta_maxD0Cut = new TH1F("h_match_tp_eta_maxD0Cut", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
  TH1F* h_tp_d0_maxD0Cut = new TH1F("h_tp_d0_maxD0Cut", ";Tracking particle d_{0} [cm]; Tracking particles / 0.08 cm", 50, -2, 2);
  TH1F* h_match_tp_d0_maxD0Cut = new TH1F("h_match_tp_d0_maxD0Cut", ";Tracking particle d_{0} [cm]; Tracking particles / 0.08 cm", 50, -2, 2);
   
  TH1F* h_tp_pt_minD0Cut = new TH1F("h_tp_pt_minD0Cut", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
  TH1F* h_match_tp_pt_minD0Cut = new TH1F("h_match_tp_pt_minD0Cut", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
  TH1F* h_tp_eta_minD0Cut = new TH1F("h_tp_eta_minD0Cut", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
  TH1F* h_match_tp_eta_minD0Cut = new TH1F("h_match_tp_eta_minD0Cut", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
  TH1F* h_tp_d0_minD0Cut = new TH1F("h_tp_d0_minD0Cut", ";Tracking particle d_{0} [cm]; Tracking particles / 0.08 cm", 50, -2, 2);
  TH1F* h_match_tp_d0_minD0Cut = new TH1F("h_match_tp_d0_minD0Cut", ";Tracking particle d_{0} [cm]; Tracking particles / 0.08 cm", 50, -2, 2);

  TH1F* h_tp_pt_minPtCut = new TH1F("h_tp_pt_minPtCut", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
  TH1F* h_match_tp_pt_minPtCut = new TH1F("h_match_tp_pt_minPtCut", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
  TH1F* h_tp_eta_minPtCut = new TH1F("h_tp_eta_minPtCut", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
  TH1F* h_match_tp_eta_minPtCut = new TH1F("h_match_tp_eta_minPtCut", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
  TH1F* h_tp_d0_minPtCut = new TH1F("h_tp_d0_minPtCut", ";Tracking particle d_{0} [cm]; Tracking particles / 0.08 cm", 50, -2, 2);
  TH1F* h_match_tp_d0_minPtCut = new TH1F("h_match_tp_d0_minPtCut", ";Tracking particle d_{0} [cm]; Tracking particles / 0.08 cm", 50, -2, 2);

  TH1F* h_tp_pt_maxPtCut = new TH1F("h_tp_pt_maxPtCut", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
  TH1F* h_match_tp_pt_maxPtCut = new TH1F("h_match_tp_pt_maxPtCut", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
  TH1F* h_tp_eta_maxPtCut = new TH1F("h_tp_eta_maxPtCut", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
  TH1F* h_match_tp_eta_maxPtCut = new TH1F("h_match_tp_eta_maxPtCut", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
  TH1F* h_tp_d0_maxPtCut = new TH1F("h_tp_d0_maxPtCut", ";Tracking particle d_{0} [cm]; Tracking particles / 0.08 cm", 50, -2, 2);
  TH1F* h_match_tp_d0_maxPtCut = new TH1F("h_match_tp_d0_maxPtCut", ";Tracking particle d_{0} [cm]; Tracking particles / 0.08 cm", 50, -2, 2);

  TH1F* h_tp_pt_maxEtaCut = new TH1F("h_tp_pt_maxEtaCut", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
  TH1F* h_match_tp_pt_maxEtaCut = new TH1F("h_match_tp_pt_maxEtaCut", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
  TH1F* h_tp_eta_maxEtaCut = new TH1F("h_tp_eta_maxEtaCut", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
  TH1F* h_match_tp_eta_maxEtaCut = new TH1F("h_match_tp_eta_maxEtaCut", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
  TH1F* h_tp_d0_maxEtaCut = new TH1F("h_tp_d0_maxEtaCut", ";Tracking particle d_{0} [cm]; Tracking particles / 0.08 cm", 50, -2, 2);
  TH1F* h_match_tp_d0_maxEtaCut = new TH1F("h_match_tp_d0_maxEtaCut", ";Tracking particle d_{0} [cm]; Tracking particles / 0.08 cm", 50, -2, 2);

  TH1F* h_tp_pt_allCuts = new TH1F("h_tp_pt_allCuts", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
  TH1F* h_match_tp_pt_allCuts = new TH1F("h_match_tp_pt_allCuts", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
  TH1F* h_match_tp_pt_allCuts_trkCuts = new TH1F("h_match_tp_pt_allCuts_trkCuts", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
  TH1F* h_tp_eta_allCuts = new TH1F("h_tp_eta_allCuts", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
  TH1F* h_match_tp_eta_allCuts = new TH1F("h_match_tp_eta_allCuts", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
  TH1F* h_tp_d0_allCuts = new TH1F("h_tp_d0_allCuts", ";Tracking particle d_{0} [cm]; Tracking particles / 0.08 cm", 50, -2, 2);
  TH1F* h_match_tp_d0_allCuts = new TH1F("h_match_tp_d0_allCuts", ";Tracking particle d_{0} [cm]; Tracking particles / 0.08 cm", 50, -2, 2);

  TH1F* h_trk_matchtp_pt_allCuts = new TH1F("h_trk_matchtp_pt_allCuts", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
  TH1F* h_trk_matchtp_eta_allCuts = new TH1F("h_trk_matchtp_eta_allCuts", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
  TH1F* h_tp_dxy_allCuts = new TH1F("h_tp_dxy_allCuts", ";Tracking particle d_{xy} [cm]; Tracking particles / 0.08 cm", 50, -2, 2);
  TH1F* h_trk_matchtp_dxy_allCuts = new TH1F("h_trk_matchtp_dxy_allCuts", ";Tracking particle d_{xy} [cm]; Tracking particles / 0.08 cm", 50, -2, 2);

  TH1F* h_trk_matchtp_pt_noCuts = new TH1F("h_trk_matchtp_pt_noCuts", ";Tracking particle p_{T} [GeV]; Tracking particles / 0.4 GeV", 100, 0, 100.0);
  TH1F* h_trk_matchtp_eta_noCuts = new TH1F("h_trk_matchtp_eta_noCuts", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
  TH1F* h_tp_dxy_noCuts = new TH1F("h_tp_dxy_noCuts", ";Tracking particle d_{xy} [cm]; Tracking particles / 0.08 cm", 50, -2, 2);
  TH1F* h_trk_matchtp_dxy_noCuts = new TH1F("h_trk_matchtp_dxy_noCuts", ";Tracking particle d_{xy} [cm]; Tracking particles / 0.08 cm", 50, -2, 2);

  TH1F* h_trk_matchtp_pt_oldCuts = new TH1F("h_trk_matchtp_pt_oldCuts", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
  TH1F* h_trk_matchtp_eta_oldCuts = new TH1F("h_trk_matchtp_eta_oldCuts", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
  TH1F* h_tp_dxy_oldCuts = new TH1F("h_tp_dxy_oldCuts", ";Tracking particle d_{xy} [cm]; Tracking particles / 0.08 cm", 50, -2, 2);
  TH1F* h_trk_matchtp_dxy_oldCuts = new TH1F("h_trk_matchtp_dxy_oldCuts", ";Tracking particle d_{xy} [cm]; Tracking particles / 0.08 cm", 50, -2, 2);

  // Displaced Vertex Plots
  TH1F *h_delta_dist_xy = new TH1F("h_delta_dist_xy","h_delta_dist_xy; (a) Distance between chosen TPs in x-y (cm) ; Events / 0.005 cm",100,0,0.5);
  TH1F *h_error_delta_x = new TH1F("h_error_delta_x","h_error_delta_x; (b) #Delta x_{displaced} with chosen TP true x (cm) ; Events / 0.02 cm",100,0,2);
  TH1F *h_delta_dist_z = new TH1F("h_delta_dist_z","h_delta_dist_z; (c) Distance between chosen TPs in z (cm) ; Events / 0.05 cm",40,0,2);
  TH1F *h_error_delta_z = new TH1F("h_error_delta_z","h_error_delta_z; (d) Chosen TP error in z (cm) ; Events / 0.1 cm",100,0,10);
  TH1F *h_delta_x = new TH1F("h_delta_x","h_delta_x; #Delta x between chosen TPs true x (cm) ; Events / 0.05 cm",40,0,2);
  TH1F *h_trueVertexAssoc_delxy = new TH1F("h_trueVertexAssoc_delxy","h_trueVertexAssoc_delxy; Distance between TP and trueVertex in x-y (cm) ; Events / 0.04 cm",100,0,4);
  TH1F *h_trueVertexAssoc_delz = new TH1F("h_trueVertexAssoc_delz","h_trueVertexAssoc_delz; Distance between TP and trueVertex in z (cm) ; Events / 0.04 cm",100,0,4);
  TH1F *h_trueVertexAssoc_calcVertDelxy = new TH1F("h_trueVertexAssoc_calcVertDelxy","h_trueVertexAssoc_calcVertDelxy; Distance between TP and trueVertex in x-y (cm) ; Events / 0.04 cm",100,0,4);
  TH1F *h_trueVertexAssoc_calcVertDelz = new TH1F("h_trueVertexAssoc_calcVertDelz","h_trueVertexAssoc_calcVertDelz; Distance between TP and trueVertex in z (cm) ; Events / 0.04 cm",100,0,4);
  TH1F *h_trackVertexAssoc_delxy = new TH1F("h_trackVertexAssoc_delxy","h_trackVertexAssoc_delxy; Distance between track and trackVertex in x-y (cm) ; Events / 0.04 cm",100,0,4);
  TH1F *h_trackVertexAssoc_delz = new TH1F("h_trackVertexAssoc_delz","h_trackVertexAssoc_delz; Distance between track and trackVertex in z (cm) ; Events / 0.04 cm",100,0,4);
  TH1F* h_tpAssoc_pt = new TH1F("h_tpAssoc_pt", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
  TH1F* h_tpAssoc_pt_matched = new TH1F("h_tpAssoc_pt_matched", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
  TH1F* h_trkAssoc_pt = new TH1F("h_trkAssoc_pt", ";Track p_{T} [GeV]; Tracks / 1.0 GeV", 100, 0, 100.0);
  TH1F* h_trkAssoc_pt_noMatch = new TH1F("h_trkAssoc_pt_noMatch", ";Track p_{T} [GeV]; Tracks / 1.0 GeV", 100, 0, 100.0);
  TH1F* h_correct_trkAssoc_delPhi = new TH1F("h_correct_trkAssoc_delPhi", ";Delta Phi Between Vertex Momentum and Associated Track; Tracks / 0.006",100,0,0.6);
  TH1F* h_false_trkAssoc_delPhi = new TH1F("h_false_trkAssoc_delPhi", ";Delta Phi Between Vertex Momentum and Associated Track; Tracks / 0.006",100,0,0.6);
  TH1F* h_correct_trkAssoc_delPhiProp = new TH1F("h_correct_trkAssoc_delPhiProp", ";Delta Phi Between Vertex Momentum and Track's at Vertex; Tracks / 0.006",100,0,0.6);
  TH1F* h_false_trkAssoc_delPhiProp = new TH1F("h_false_trkAssoc_delPhiProp", ";Delta Phi Between Vertex Momentum and Track's at Vertex; Tracks / 0.006",100,0,0.6);
  TH1F* h_correct_trkAssoc_delxy = new TH1F("h_correct_trkAssoc_delxy", ";Distance between track and trackVertex in x-y (cm); Tracks / 0.001 cm",100,0,0.1);
  TH1F* h_false_trkAssoc_delxy = new TH1F("h_false_trkAssoc_delxy", ";Distance between track and trackVertex in x-y (cm); Tracks / 0.001 cm",100,0,0.1);
  TH1F* h_correct_trkAssoc_delz = new TH1F("h_correct_trkAssoc_delz", ";Distance between track and trackVertex in z (cm); Tracks / 0.005 cm",100,0,0.5);
  TH1F* h_false_trkAssoc_delz = new TH1F("h_false_trkAssoc_delz", ";Distance between track and trackVertex in z (cm); Tracks / 0.005 cm",100,0,0.5);
  TH1F* h_correct_trkAssoc_pt = new TH1F("h_correct_trkAssoc_pt", ";Track p_{T} [GeV]; Tracks / 1.0 GeV", 100, 0, 100.0);
  TH1F* h_false_trkAssoc_pt = new TH1F("h_false_trkAssoc_pt", ";Track p_{T} [GeV]; Tracks / 1.0 GeV", 100, 0, 100.0);
  TH1F* h_correct_trkAssoc_delPt = new TH1F("h_correct_trkAssoc_delPt", ";Delta p_{T} Between Track and Leading Track from Vertex [GeV]; Tracks / 1.0 GeV", 100, 0, 100.0);
  TH1F* h_false_trkAssoc_delPt = new TH1F("h_false_trkAssoc_delPt", ";Delta p_{T} Between Track and Leading Track from Vertex [GeV]; Tracks / 1.0 GeV", 100, 0, 100.0);

  TH1I *h_trk_Counter_TPcombination = new TH1I("h_trk_Counter_TPcombination","h_trk_Counter_TPcombination; Track combination chosen; Events / 1.0",6,0,6);
  TH1F *h_trk_delta_dist_xy = new TH1F("h_trk_delta_dist_xy","h_trk_delta_dist_xy; Distance between chosen Tracks in x-y (cm) ; Events / 2E-8 cm",50,0,0.000001);
  TH1F *h_trk_delta_dist_z = new TH1F("h_trk_delta_dist_z","h_trk_delta_dist_z; Distance between chosen Tracks in z (cm) ; Events / 0.05 cm",40,0,2);
  TH1F *h_trackVertex_cos_T = new TH1F("h_trackVertex_cos_T","h_trackVertex_cos_T; Cos(angle): parent momentum and vertex position; Events / 0.05",40,-1,1);
  TH1F *h_trackVertex_alpha_T = new TH1F("h_trackVertex_alpha_T","h_trackVertex_alpha_T; angle btwn parent momentum and vertex position; Events / 0.315",40,-6.3,6.3);
  TH1F *h_trackVertex_openingAngle = new TH1F("h_trackVertex_openingAngle","h_trackVertex_openingAngle; angle btwn tracks; Events / 0.157",40,-3.14,3.14);
  TH1F *h_trackVertex_parentPt = new TH1F("h_trackVertex_parentPt","h_trackVertex_parentPt; Pt magnitude of parent; Events / 1.0",200,0,200.0);
  TH1F *h_trackVertex_d_T = new TH1F("h_trackVertex_d_T","h_trackVertex_d_T; Impact parameter of parent (cm) ; Events / 0.005 cm",40,0,0.2);
  TH1F *h_trackVertex_R_T = new TH1F("h_trackVertex_R_T","h_trackVertex_R_T; Tranverse distance of vertex (cm) ; Events / 0.25 cm",80,0,20);
  TH1F *h_trueVertex_cos_T = new TH1F("h_trueVertex_cos_T","h_trueVertex_cos_T; Cos(angle): parent momentum and vertex position; Events / 0.05",40,-1,1);
  TH1F *h_trueVertex_alpha_T = new TH1F("h_trueVertex_alpha_T","h_trueVertex_alpha_T; angle btwn parent momentum and vertex position; Events / 0.315",40,-6.3,6.3);
  TH1F *h_trueVertex_openingAngle = new TH1F("h_trueVertex_openingAngle","h_trueVertex_openingAngle; angle btwn TPs; Events / 0.157",40,-3.14,3.14);
  TH1F *h_trueVertex_parentPt = new TH1F("h_trueVertex_parentPt","h_trueVertex_parentPt; Pt magnitude of parent; Events / 1.0",200,0,200.0);
  TH1F *h_trueVertex_d_T = new TH1F("h_trueVertex_d_T","h_trueVertex_d_T; Impact parameter of parent (cm) ; Events / 0.005 cm",40,0,0.2);
  TH1F *h_trueVertex_R_T = new TH1F("h_trueVertex_R_T","h_trueVertex_R_T; Tranverse distance of vertex (cm) ; Events / 0.25 cm",80,0,20);
  TH1F *h_trueVertex_delta_dist_z0 = new TH1F("h_trueVertex_delta_dist_z0","h_trueVertex_delta_dist_z0; Distance in z0 btwn two hardest TPs from vertex (cm) ; Events / 0.1 cm",100,0,10);
  TH1F *h_trueVertex_delta_dist_d0 = new TH1F("h_trueVertex_delta_dist_d0","h_trueVertex_delta_dist_d0; Distance in d0 btwn two hardest TPs from vertex (cm) ; Events / 0.1 cm",100,0,10);
  TH1F *h_trueVertex_delta_dist_eta = new TH1F("h_trueVertex_delta_dist_eta","h_trueVertex_delta_dist_eta; Distance in eta btwn two hardest TPs from vertex ; Events / 0.024",100,0,2.4);
  TH1F *h_trueVertex_delta_dist_phi = new TH1F("h_trueVertex_delta_dist_phi","h_trueVertex_delta_dist_phi; Distance in phi btwn two hardest TPs from vertex ; Events / 0.063",100,0,6.3);
  TH1F *h_trueVertex_delta_dist_R = new TH1F("h_trueVertex_delta_dist_R","h_trueVertex_delta_dist_R; Distance in R btwn two hardest TPs from vertex ; Events / 0.063",100,0,6.3);
  TH1F *h_trueVertex_delta_dist_indexPt = new TH1F("h_trueVertex_delta_dist_indexPt","h_trueVertex_delta_dist_indexPt; Distance in pt index btwn two hardest TPs from vertex ; Events / 1.0",20,0,20);
  TH1F *h_trueVertex_delta_dist_pt = new TH1F("h_trueVertex_delta_dist_pt","h_trueVertex_delta_dist_pt; Distance in pt btwn two hardest TPs from vertex (GeV); Events / 1.0",50,0,50);
  
  TH1F *h_res_tp_trk_x = new TH1F("h_res_tp_trk_x","h_res_tp_trk_x; x residual of vertex (cm) ; Events / 0.02 cm",100,-1,1);
  TH1F *h_res_tp_trk_y = new TH1F("h_res_tp_trk_y","h_res_tp_trk_y; y residual of vertex (cm) ; Events / 0.02 cm",100,-1,1);
  TH1F *h_res_tp_trk_x_zoomOut = new TH1F("h_res_tp_trk_x_zoomOut","h_res_tp_trk_x_zoomOut; x residual of vertex (cm) ; Events / 0.022 cm",500,-1,10);
  TH1F *h_res_tp_trk_y_zoomOut = new TH1F("h_res_tp_trk_y_zoomOut","h_res_tp_trk_y_zoomOut; y residual of vertex (cm) ; Events / 0.022 cm",500,-1,10);
  TH1F *h_res_tp_trk_x_findVert = new TH1F("h_res_tp_trk_x_findVert","h_res_tp_trk_x_findVert; x residual of vertex (cm) ; Events / 0.022 cm",500,-1,10);
  TH1F *h_res_tp_trk_y_findVert = new TH1F("h_res_tp_trk_y_findVert","h_res_tp_trk_y_findVert; y residual of vertex (cm) ; Events / 0.022 cm",500,-1,10);
  TH1F *h_res_tp_trk_z = new TH1F("h_res_tp_trk_z","h_res_tp_trk_z; z residual of vertex (cm) ; Events / 0.02 cm",1000,-10,10);
  TH1F *h_res_tp_trk_r = new TH1F("h_res_tp_trk_r","h_res_tp_trk_r; r residual of vertex (cm) ; Events / 0.02 cm",100,-1,1);
  TH1F *h_res_tp_trk_phi = new TH1F("h_res_tp_trk_phi","h_res_tp_trk_phi; phi residual of vertex ; Events / 0.02",100,-1,1);

  // Efficiency of Identifying Displaced Vertex Plots
  TH1F *h_all_trueVertex_pt = new TH1F("h_all_trueVertex_pt","h_all_trueVertex_pt; p_{T} of Leading p_{T} TP from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_all_trueVertex_minD0 = new TH1F("h_all_trueVertex_minD0","h_all_trueVertex_minD0; minimum d_{0} of two hardest TPs from vertex (cm); Events / 0.05 cm",80,-2,2);
  TH1F *h_all_trueVertex_minD0_allTPs = new TH1F("h_all_trueVertex_minD0_allTPs","h_all_trueVertex_minD0_allTPs; minimum d_{0} of TPs from vertex (cm); Events / 0.02 cm",100,0,2);
  TH1F *h_all_trueVertex_maxD0 = new TH1F("h_all_trueVertex_maxD0","h_all_trueVertex_maxD0; maximum d_{0} of two hardest TPs from vertex (cm); Events / 0.05 cm",80,-2,2);
  TH1F *h_all_trueVertex_maxD0_allTPs = new TH1F("h_all_trueVertex_maxD0_allTPs","h_all_trueVertex_maxD0_allTPs; maximum d_{0} of TPs from vertex (cm); Events / 0.025 cm",100,0,2.5);
  TH1F *h_all_trueVertex_lowPt = new TH1F("h_all_trueVertex_lowPt","h_all_trueVertex_lowPt; Lower p_{T} TP of two hardest TPs from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_all_trueVertex_lowPt_allTPs = new TH1F("h_all_trueVertex_lowPt_allTPs","h_all_trueVertex_lowPt_allTPs; Minimum p_{T} TP from vertex (GeV); Events / 0.25 GeV",100,0,25);
  TH1F *h_all_oneMatch_trueVertex_pt = new TH1F("h_all_oneMatch_trueVertex_pt","h_all_oneMatch_trueVertex_pt; p_{T} of Leading p_{T} TP from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_all_oneMatch_trueVertex_lowPt = new TH1F("h_all_oneMatch_trueVertex_lowPt","h_all_oneMatch_trueVertex_lowPt; p_{T} of Lower p_{T} TP from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_all_oneMatchAlt_trueVertex_pt = new TH1F("h_all_oneMatchAlt_trueVertex_pt","h_all_oneMatchAlt_trueVertex_pt; p_{T} of Leading p_{T} TP from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_all_trackVertex_pt = new TH1F("h_all_trackVertex_pt","h_all_trackVertex_pt; p_{T} of Leading p_{T} trk from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_correct_trackVertex_deltaPos = new TH1F("h_correct_trackVertex_deltaPos","h_correct_trackVertex_deltaPos; Distance Between Vertices (cm); Events / 0.2 cm",100,0,20);
  TH1F *h_false_trackVertex_deltaPos = new TH1F("h_false_trackVertex_deltaPos","h_false_trackVertex_deltaPos; Distance Between Vertices (cm); Events / 0.2 cm",100,0,20);
  TH1F *h_correct_trackVertex_pt = new TH1F("h_correct_trackVertex_pt","h_correct_trackVertex_pt; p_{T} of Leading p_{T} trk from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_false_trackVertex_maxPt = new TH1F("h_false_trackVertex_maxPt","h_false_trackVertex_maxPt; p_{T} of Leading p_{T} trk from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_all_trackVertex_minD0 = new TH1F("h_all_trackVertex_minD0","h_all_trackVertex_minD0; minimum d_{0} of trks from vertex (cm); Events / 0.05 cm",80,-2,2);
  TH1F *h_correct_trackVertex_minD0 = new TH1F("h_correct_trackVertex_minD0","h_correct_trackVertex_minD0; minimum d_{0} of trks from vertex (cm); Events / 0.025 cm",80,0,2);
  TH1F *h_false_trackVertex_minD0 = new TH1F("h_false_trackVertex_minD0","h_false_trackVertex_minD0; minimum d_{0} of trks from vertex (cm); Events / 0.025 cm",80,0,2);
  TH1F *h_correct_trackVertex_maxD0 = new TH1F("h_correct_trackVertex_maxD0","h_correct_trackVertex_maxD0; maximum d_{0} of trks from vertex (cm); Events / 0.05 cm",80,-2,2);
  TH1F *h_false_trackVertex_maxD0 = new TH1F("h_false_trackVertex_maxD0","h_false_trackVertex_maxD0; maximum d_{0} of trks from vertex (cm); Events / 0.05 cm",80,-2,2);
  TH1F *h_all_trackVertex_lowPt = new TH1F("h_all_trackVertex_lowPt","h_all_trackVertex_lowPt; p_{T} of Lower p_{T} trk from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_correct_trueVertex_pt = new TH1F("h_correct_trueVertex_pt","h_correct_trueVertex_pt; p_{T} of Leading p_{T} TP from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_findable_trueVertex_pt = new TH1F("h_findable_trueVertex_pt","h_findable_trueVertex_pt; p_{T} of Leading p_{T} TP from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_findable_trueVertexBinned_pt = new TH1F("h_findable_trueVertexBinned_pt","h_findable_trueVertexBinned_pt; p_{T} of Leading p_{T} TP from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_findableIntersect_trueVertexBinned_pt = new TH1F("h_findableIntersect_trueVertexBinned_pt","h_findableIntersect_trueVertexBinned_pt; p_{T} of Leading p_{T} TP from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_findableNonPrompt_trueVertexBinned_pt = new TH1F("h_findableNonPrompt_trueVertexBinned_pt","h_findableNonPrompt_trueVertexBinned_pt; p_{T} of Leading p_{T} TP from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_findableBeforeTracker_trueVertexBinned_pt = new TH1F("h_findableBeforeTracker_trueVertexBinned_pt","h_findableBeforeTracker_trueVertexBinned_pt; p_{T} of Leading p_{T} TP from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_findableBendChi2Cut_trueVertexBinned_pt = new TH1F("h_findableBendChi2Cut_trueVertexBinned_pt","h_findableBendChi2Cut_trueVertexBinned_pt; p_{T} of Leading p_{T} TP from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_findableChi2RPhiCut_trueVertexBinned_pt = new TH1F("h_findableChi2RPhiCut_trueVertexBinned_pt","h_findableChi2RPhiCut_trueVertexBinned_pt; p_{T} of Leading p_{T} TP from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_findableCosTCut_trueVertexBinned_pt = new TH1F("h_findableCosTCut_trueVertexBinned_pt","h_findableCosTCut_trueVertexBinned_pt; p_{T} of Leading p_{T} TP from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_findableDeltaZCut_trueVertexBinned_pt = new TH1F("h_findableDeltaZCut_trueVertexBinned_pt","h_findableDeltaZCut_trueVertexBinned_pt; p_{T} of Leading p_{T} TP from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_findFake_trackVertex_pt = new TH1F("h_findFake_trackVertex_pt","h_findFake_trackVertex_pt; p_{T} of Leading p_{T} track from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_findFake_trackVertexBinned_pt = new TH1F("h_findFake_trackVertexBinned_pt","h_findFake_trackVertexBinned_pt; p_{T} of Leading p_{T} track from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_findFakeIntersect_trackVertexBinned_pt = new TH1F("h_findFakeIntersect_trackVertexBinned_pt","h_findFakeIntersect_trackVertexBinned_pt; p_{T} of Leading p_{T} track from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_findFakeNonPrompt_trackVertexBinned_pt = new TH1F("h_findFakeNonPrompt_trackVertexBinned_pt","h_findFakeNonPrompt_trackVertexBinned_pt; p_{T} of Leading p_{T} track from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_findFakeBeforeTracker_trackVertexBinned_pt = new TH1F("h_findFakeBeforeTracker_trackVertexBinned_pt","h_findFakeBeforeTracker_trackVertexBinned_pt; p_{T} of Leading p_{T} track from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_findFakeBendChi2Cut_trackVertexBinned_pt = new TH1F("h_findFakeBendChi2Cut_trackVertexBinned_pt","h_findFakeBendChi2Cut_trackVertexBinned_pt; p_{T} of Leading p_{T} track from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_findFakeChi2RPhiCut_trackVertexBinned_pt = new TH1F("h_findFakeChi2RPhiCut_trackVertexBinned_pt","h_findFakeChi2RPhiCut_trackVertexBinned_pt; p_{T} of Leading p_{T} track from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_findFakeCosTCut_trackVertexBinned_pt = new TH1F("h_findFakeCosTCut_trackVertexBinned_pt","h_findFakeCosTCut_trackVertexBinned_pt; p_{T} of Leading p_{T} track from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_findFakeDeltaZCut_trackVertexBinned_pt = new TH1F("h_findFakeDeltaZCut_trackVertexBinned_pt","h_findFakeDeltaZCut_trackVertexBinned_pt; p_{T} of Leading p_{T} track from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_correct_trueVertex_lowPt = new TH1F("h_correct_trueVertex_lowPt","h_correct_trueVertex_lowPt; p_{T} of Lower p_{T} TP from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_correct_oneMatch_trueVertex_pt = new TH1F("h_correct_oneMatch_trueVertex_pt","h_correct_oneMatch_trueVertex_pt; p_{T} of Leading p_{T} TP from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_correct_oneMatch_trueVertex_lowPt = new TH1F("h_correct_oneMatch_trueVertex_lowPt","h_correct_oneMatch_trueVertex_lowPt; p_{T} of Lower p_{T} TP from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_correct_oneMatchAlt_trueVertex_pt = new TH1F("h_correct_oneMatchAlt_trueVertex_pt","h_correct_oneMatchAlt_trueVertex_pt; p_{T} of Leading p_{T} TP from vertex (GeV); Events / 1.0 GeV",100,0,100);
  TH1F *h_false_trackVertex_pt = new TH1F("h_false_trackVertex_pt","h_false_trackVertex_pt; p_{T} of Leading p_{T} trk from vertex (GeV); Events / 1.0 GeV",100,0,100);

  TH1F *h_false_trackVertex_delta_dist_z = new TH1F("h_false_trackVertex_delta_dist_z","h_false_trackVertex_delta_dist_z; Distance in z btwn trks from vertex (cm); Events / 0.015 cm",100,0,1.5);
  TH1F *h_correct_trackVertex_delta_dist_z = new TH1F("h_correct_trackVertex_delta_dist_z","h_correct_trackVertex_delta_dist_z; Distance in z btwn trks from vertex (cm); Events / 0.015 cm",100,0,1.5);
  TH1F *h_false_trackVertex_delta_dist_z_inBothTraj = new TH1F("h_false_trackVertex_delta_dist_z_inBothTraj","h_false_trackVertex_delta_dist_z_inBothTraj; Distance in z btwn trks from vertex (cm); Events / 0.015 cm",100,0,1.5);
  TH1F *h_correct_trackVertex_delta_dist_z_inBothTraj = new TH1F("h_correct_trackVertex_delta_dist_z_inBothTraj","h_correct_trackVertex_delta_dist_z_inBothTraj; Distance in z btwn trks from vertex (cm); Events / 0.015 cm",100,0,1.5);
  TH1F *h_false_trackVertex_delta_dist_z_inOneTraj = new TH1F("h_false_trackVertex_delta_dist_z_inOneTraj","h_false_trackVertex_delta_dist_z_inOneTraj; Distance in z btwn trks from vertex (cm); Events / 0.015 cm",100,0,1.5);
  TH1F *h_correct_trackVertex_delta_dist_z_inOneTraj = new TH1F("h_correct_trackVertex_delta_dist_z_inOneTraj","h_correct_trackVertex_delta_dist_z_inOneTraj; Distance in z btwn trks from vertex (cm); Events / 0.015 cm",100,0,1.5);
  TH1F *h_false_trackVertex_delta_dist_z_inNoTraj = new TH1F("h_false_trackVertex_delta_dist_z_inNoTraj","h_false_trackVertex_delta_dist_z_inNoTraj; Distance in z btwn trks from vertex (cm); Events / 0.015 cm",100,0,1.5);
  TH1F *h_correct_trackVertex_delta_dist_z_inNoTraj = new TH1F("h_correct_trackVertex_delta_dist_z_inNoTraj","h_correct_trackVertex_delta_dist_z_inNoTraj; Distance in z btwn trks from vertex (cm); Events / 0.015 cm",100,0,1.5);
  TH1F *h_correct_trackVertex_delta_dist_z0 = new TH1F("h_correct_trackVertex_delta_dist_z0","h_correct_trackVertex_delta_dist_z0; Distance in z0 btwn trks from vertex (cm); Events / 0.1 cm",100,0,10);
  TH1F *h_correct_trackVertex_delta_dist_d0 = new TH1F("h_correct_trackVertex_delta_dist_d0","h_correct_trackVertex_delta_dist_d0; Distance in d0 btwn trks from vertex (cm); Events / 0.1 cm",100,0,10);
  TH1F *h_correct_trackVertex_delta_dist_eta = new TH1F("h_correct_trackVertex_delta_dist_eta","h_correct_trackVertex_delta_dist_eta; Distance in eta btwn trks from vertex ; Events / 0.024",100,0,2.4);
  TH1F *h_false_trackVertex_delta_dist_eta = new TH1F("h_false_trackVertex_delta_dist_eta","h_false_trackVertex_delta_dist_eta; Distance in eta btwn trks from vertex ; Events / 0.024",100,0,2.4);
  TH1F *h_correct_trackVertex_delta_dist_phi = new TH1F("h_correct_trackVertex_delta_dist_phi","h_correct_trackVertex_delta_dist_phi; Distance in phi btwn trks from vertex ; Events / 0.063",100,0,6.3);
  TH1F *h_false_trackVertex_delta_dist_phi = new TH1F("h_false_trackVertex_delta_dist_phi","h_false_trackVertex_delta_dist_phi; Distance in phi btwn trks from vertex ; Events / 0.063",100,0,6.3);
  TH1F *h_correct_trackVertex_delta_dist_indexPt = new TH1F("h_correct_trackVertex_delta_dist_indexPt","h_correct_trackVertex_delta_dist_indexPt; Distance in pt index btwn trks from vertex ; Events / 1.0",20,0,20);
  TH1F *h_all_trackVertex_numStubs = new TH1F("h_all_trackVertex_numStubs","h_all_trackVertex_numStubs; Number of stubs of trks from vertex; Events / 1.0",12,0,12);
  TH1F *h_all_trackVertex_fakeId = new TH1F("h_all_trackVertex_fakeId","h_all_trackVertex_fakeId; Fake Id of Tracks; Events / 1.0",4,0,4);
  TH1F *h_correct_trackVertex_numStubs = new TH1F("h_correct_trackVertex_numStubs","h_correct_trackVertex_numStubs; Number of stubs of trks from vertex; Events / 1.0",12,0,12);
  TH1F *h_correct_trackVertex_numStubsSum = new TH1F("h_correct_trackVertex_numStubsSum","h_correct_trackVertex_numStubsSum; Number of stubs of trks from vertex; Events / 1.0",12,0,12);
  TH1F *h_correct_trackVertex_numTracks = new TH1F("h_correct_trackVertex_numTracks","h_correct_trackVertex_numTracks; Tracks Associated with Vertex; Events / 1.0",20,0,20);
  TH1F *h_false_trackVertex_numStubs = new TH1F("h_false_trackVertex_numStubs","h_false_trackVertex_numStubs; Number of stubs of trks from vertex; Events / 1.0",12,0,12);
  TH1F *h_false_trackVertex_numStubsSum = new TH1F("h_false_trackVertex_numStubsSum","h_false_trackVertex_numStubsSum; Number of stubs of trks from vertex; Events / 1.0",12,0,12);
  TH1F *h_false_trackVertex_numTracks = new TH1F("h_false_trackVertex_numTracks","h_false_trackVertex_numTracks; Tracks Associated with Vertex; Events / 1.0",20,0,20);
  TH1F *h_false_trackVertex_d0 = new TH1F("h_false_trackVertex_d0","h_false_trackVertex_d0; d_{0} of trks from vertex; Events / 0.05 cm",100,0,5);
  TH1F *h_false_trackVertex_d_T = new TH1F("h_false_trackVertex_d_T","h_false_trackVertex_d_T; Impact parameter of parent (cm) ; Events / 0.025 cm",40,0,1);
  TH1F *h_false_trackVertex_chi2rphidofSum = new TH1F("h_false_trackVertex_chi2rphidofSum","h_false_trackVertex_chi2rphidofSum; #Sigma Track #chi^{2}_{r#phi}/d.o.f ; Events / 0.04",200,0,8);
  TH1F *h_false_trackVertex_MVA1Sum = new TH1F("h_false_trackVertex_MVA1Sum","h_false_trackVertex_MVA1Sum; #Sigma Track MVA1 ; Events / 0.02",100,0,2);
  TH1F *h_false_trackVertex_MVA2Sum = new TH1F("h_false_trackVertex_MVA2Sum","h_false_trackVertex_MVA2Sum; #Sigma Track MVA2 ; Events / 0.02",100,0,2);
  TH1F *h_false_trackVertex_chi2rzdofSum = new TH1F("h_false_trackVertex_chi2rzdofSum","h_false_trackVertex_chi2rzdofSum; #Sigma Track #chi^{2}_{rz}/d.o.f ; Events / 0.015",200,0,3);
  TH1F *h_false_trackVertex_bendchi2Sum = new TH1F("h_false_trackVertex_bendchi2Sum","h_false_trackVertex_bendchi2Sum; #Sigma Track bend #chi^{2} ; Events / 0.07",200,0,14);
  TH1F *h_false_trackVertex_score = new TH1F("h_false_trackVertex_score","h_false_trackVertex_score; Vertex NN score ; Events / 0.1",100,0,1);
  TH1F *h_false_trackVertex_R_T = new TH1F("h_false_trackVertex_R_T","h_false_trackVertex_R_T; Tranverse distance of vertex (cm) ; Events / 0.25 cm",80,0,20);
  TH1F *h_false_trackVertex_cos_T = new TH1F("h_false_trackVertex_cos_T","h_false_trackVertex_cos_T; Cos(angle): parent momentum and vertex position; Events / 0.0014",50,0.93,1);
  TH1F *h_false_trackVertex_p2_mag = new TH1F("h_false_trackVertex_p2_mag","h_false_trackVertex_p2_mag; #Sigma p_{T}^2 of trks from vertex (GeV^{2}); Events / 10.0 GeV^{2}",1000,0,10000);
  TH1F *h_false_trackVertex_lowPt = new TH1F("h_false_trackVertex_lowPt","h_false_trackVertex_lowPt; p_{T} of lower p_{T} trk from vertex (GeV); Events / 1 GeV",100,0,100);
  TH1F *h_false_trackVertex_inTraj = new TH1F("h_false_trackVertex_inTraj","h_false_trackVertex_inTraj; Calc Vertex Return Code; Events / 1.0",6,0,6);
  TH1F *h_false_trackVertex_numMatched = new TH1F("h_false_trackVertex_numMatched","h_false_trackVertex_numMatched; Number of Matched Tracks; Events / 1.0",4,0,4);
  TH1F *h_false_trackVertex_fakeId = new TH1F("h_false_trackVertex_fakeId","h_false_trackVertex_fakeId; Fake Id of Tracks; Events / 1.0",4,0,4);
  TH1F *h_correct_trackVertex_chi2rphidofSum = new TH1F("h_correct_trackVertex_chi2rphidofSum","h_correct_trackVertex_chi2rphidofSum; #Sigma Track #chi^{2}_{r#phi}/d.o.f ; Events / 0.04",200,0,8);
  TH1F *h_correct_trackVertex_MVA1Sum = new TH1F("h_correct_trackVertex_MVA1Sum","h_correct_trackVertex_MVA1Sum; #Sigma Track MVA1 ; Events / 0.02",100,0,2);
  TH1F *h_correct_trackVertex_MVA2Sum = new TH1F("h_correct_trackVertex_MVA2Sum","h_correct_trackVertex_MVA2Sum; #Sigma Track MVA2 ; Events / 0.02",100,0,2);
  TH1F *h_correct_trackVertex_chi2rzdofSum = new TH1F("h_correct_trackVertex_chi2rzdofSum","h_correct_trackVertex_chi2rzdofSum; #Sigma Track #chi^{2}_{rz}/d.o.f ; Events / 0.015",200,0,3);
  TH1F *h_correct_trackVertex_bendchi2Sum = new TH1F("h_correct_trackVertex_bendchi2Sum","h_correct_trackVertex_bendchi2Sum; #Sigma Track bend #chi^{2} ; Events / 0.07",200,0,14);
  TH1F *h_correct_trackVertex_score = new TH1F("h_correct_trackVertex_score","h_correct_trackVertex_score; Vertex NN score ; Events / 0.1",100,0,1);
  TH1F *h_correct_trackVertex_d_T = new TH1F("h_correct_trackVertex_d_T","h_correct_trackVertex_d_T; Impact parameter of parent (cm) ; Events / 0.025 cm",40,0,1);
  TH1F *h_correct_trackVertex_R_T = new TH1F("h_correct_trackVertex_R_T","h_correct_trackVertex_R_T; Tranverse distance of vertex (cm) ; Events / 0.25 cm",80,0,20);
  TH1F *h_correct_trackVertex_cos_T = new TH1F("h_correct_trackVertex_cos_T","h_correct_trackVertex_cos_T; Cos(angle): parent momentum and vertex position; Events / 0.0014",50,0.93,1);
  TH1F *h_correct_trackVertex_p2_mag = new TH1F("h_correct_trackVertex_p2_mag","h_correct_trackVertex_p2_mag; #Sigma p_{T}^2 of trks from vertex (GeV^{2}); Events / 10.0 GeV^{2}",1000,0,10000);
  TH1F *h_correct_trackVertex_lowPt = new TH1F("h_correct_trackVertex_lowPt","h_correct_trackVertex_lowPt; p_{T} of lower p_{T} trk from vertex (GeV); Events / 1 GeV",100,0,100);
  TH1F *h_correct_trackVertex_inTraj = new TH1F("h_correct_trackVertex_inTraj","h_correct_trackVertex_inTraj; Calc Vertex Return Code; Events / 1.0",6,0,6);

  TH1F *h_all_trueVertex_eta = new TH1F("h_all_trueVertex_eta","h_all_trueVertex_eta; #eta of Leading p_{T} TP from vertex; Events / 0.096",50,-2.4,2.4);
  TH1F *h_all_oneMatch_trueVertex_eta = new TH1F("h_all_oneMatch_trueVertex_eta","h_all_oneMatch_trueVertex_eta; #eta of Leading p_{T} TP from vertex ; Events / 0.096",50,-2.4,2.4);
  TH1F *h_all_oneMatchAlt_trueVertex_eta = new TH1F("h_all_oneMatchAlt_trueVertex_eta","h_all_oneMatchAlt_trueVertex_eta; #eta of Leading p_{T} TP from vertex ; Events / 0.096",50,-2.4,2.4);
  TH1F *h_all_trackVertex_eta = new TH1F("h_all_trackVertex_eta","h_all_trackVertex_eta; #eta of Leading p_{T} trk from vertex ; Events / 0.096",50,-2.4,2.4);
  TH1F *h_correct_trueVertex_eta = new TH1F("h_correct_trueVertex_eta","h_correct_trueVertex_eta; #eta of Leading p_{T} TP from vertex ; Events / 0.096",50,-2.4,2.4);
  TH1F *h_correct_oneMatch_trueVertex_eta = new TH1F("h_correct_oneMatch_trueVertex_eta","h_correct_oneMatch_trueVertex_eta; #eta of Leading p_{T} TP from vertex ; Events / 0.096",50,-2.4,2.4);
  TH1F *h_correct_oneMatchAlt_trueVertex_eta = new TH1F("h_correct_oneMatchAlt_trueVertex_eta","h_correct_oneMatchAlt_trueVertex_eta; #eta of Leading p_{T} TP from vertex ; Events / 0.096",50,-2.4,2.4);
  TH1F *h_false_trackVertex_eta = new TH1F("h_false_trackVertex_eta","h_false_trackVertex_eta; #eta of Leading p_{T} trk from vertex ; Events / 0.096",50,-2.4,2.4);
   
  TH1F *h_all_trueVertex_dxy = new TH1F("h_all_trueVertex_dxy","h_all_trueVertex_dxy; dxy of vertex (cm); Events / 0.1 cm",20,0,2);
  TH1F *h_all_oneMatch_trueVertex_dxy = new TH1F("h_all_oneMatch_trueVertex_dxy","h_all_oneMatch_trueVertex_dxy; dxy of vertex (cm); Events / 0.4 cm",50,0,20);
  TH1F *h_all_oneMatchAlt_trueVertex_dxy = new TH1F("h_all_oneMatchAlt_trueVertex_dxy","h_all_oneMatchAlt_trueVertex_dxy; dxy of vertex (cm); Events / 0.1 cm",20,0,2);
  TH1F *h_all_trackVertex_dxy = new TH1F("h_all_trackVertex_dxy","h_all_trackVertex_dxy; dxy of vertex (cm); Events / 0.1 cm",20,0,2);
  TH1F *h_correct_trueVertex_dxy = new TH1F("h_correct_trueVertex_dxy","h_correct_trueVertex_dxy; dxy of vertex (cm) ; Events / 0.1 cm",20,0,2);
  TH1F *h_correct_oneMatch_trueVertex_dxy = new TH1F("h_correct_oneMatch_trueVertex_dxy","h_correct_oneMatch_trueVertex_dxy; dxy of vertex (cm) ; Events / 0.4 cm",50,0,20);
  TH1F *h_correct_oneMatchAlt_trueVertex_dxy = new TH1F("h_correct_oneMatchAlt_trueVertex_dxy","h_correct_oneMatchAlt_trueVertex_dxy; dxy of vertex (cm) ; Events / 0.1 cm",20,0,2);
  TH1F *h_false_trackVertex_dxy = new TH1F("h_false_trackVertex_dxy","h_false_trackVertex_dxy; dxy of vertex (cm) ; Events / 0.1 cm",20,0,2);

  // Trigger Rates Study 
  TH2F *h_Count_trk_pt_d0 = new TH2F("h_Count_trk_pt_d0","h_Count_trk_pt_d0; Transverse Impact Parameter d_{0} (cm); Transverse Momentum p_{T}(GeV)",5,0.1,0.6,5,2,12); // Count of Selected tracks
  TH2F *h_Count_trk_pt_d0_dv = new TH2F("h_Count_trk_pt_d0_dv","h_Count_trk_pt_d0_dv;(DV selection) Transverse Impact Parameter d_{0} (cm); Transverse Momentum p_{T}(GeV)",5,0.1,0.6,5,2,12); // Count of Selected tracks including the displaced vertex selection
  TH2F *h_trk_pt_vs_d0_primary_qualCuts = new TH2F("h_trk_pt_vs_d0_primary_qualCuts","h_trk_pt_vs_d0_primary_qualCuts; Transverse Impact Parameter d_{0} (cm); Transverse Momentum p_{T} (GeV)",200,-2,2,200,0,30);
  TH2F *h_trk_pt_vs_d0_np_qualCuts = new TH2F("h_trk_pt_vs_d0_np_qualCuts","h_trk_pt_vs_d0_np_qualCuts; Transverse Impact Parameter d_{0} (cm); Transverse Momentum p_{T} (GeV)",200,-2,2,200,0,30);
  TH2F *h_trk_pt_vs_eta_primary_qualCuts = new TH2F("h_trk_pt_vs_eta_primary_qualCuts","h_trk_pt_vs_eta_primary_qualCuts; Eta; Transverse Momentum p_{T} (GeV)",200,-2.4,2.4,200,0,30);
  TH2F *h_trk_pt_vs_eta_np_qualCuts = new TH2F("h_trk_pt_vs_eta_np_qualCuts","h_trk_pt_vs_eta_np_qualCuts; Eta; Transverse Momentum p_{T} (GeV)",200,-2.4,2.4,200,0,30);
  TH2F *h_trk_eta_vs_d0_primary_qualCuts = new TH2F("h_trk_eta_vs_d0_primary_qualCuts","h_trk_eta_vs_d0_primary_qualCuts; Transverse Impact Parameter d_{0} (cm); Eta",200,-2,2,200,-2.4,2.4);
  TH2F *h_trk_eta_vs_d0_np_qualCuts = new TH2F("h_trk_eta_vs_d0_np_qualCuts","h_trk_eta_vs_d0_np_qualCuts; Transverse Impact Parameter d_{0} (cm); Eta",200,-2,2,200,-2.4,2.4);
  
  TH2F *h_trk_pt_vs_d0_primary_allCuts = new TH2F("h_trk_pt_vs_d0_primary_allCuts","h_trk_pt_vs_d0_primary_allCuts; Transverse Impact Parameter d_{0} (cm); Transverse Momentum p_{T} (GeV)",200,-2,2,200,0,30);
  TH2F *h_trk_pt_vs_d0_np_allCuts = new TH2F("h_trk_pt_vs_d0_np_allCuts","h_trk_pt_vs_d0_np_allCuts; Transverse Impact Parameter d_{0} (cm); Transverse Momentum p_{T} (GeV)",200,-2,2,200,0,30);
  TH2F *h_trk_pt_vs_eta_primary_allCuts = new TH2F("h_trk_pt_vs_eta_primary_allCuts","h_trk_pt_vs_eta_primary_allCuts; Eta; Transverse Momentum p_{T} (GeV)",200,-2.4,2.4,200,0,30);
  TH2F *h_trk_pt_vs_eta_np_allCuts = new TH2F("h_trk_pt_vs_eta_np_allCuts","h_trk_pt_vs_eta_np_allCuts; Eta; Transverse Momentum p_{T} (GeV)",200,-2.4,2.4,200,0,30);
  TH2F *h_trk_eta_vs_d0_primary_allCuts = new TH2F("h_trk_eta_vs_d0_primary_allCuts","h_trk_eta_vs_d0_primary_allCuts; Transverse Impact Parameter d_{0} (cm); Eta",200,-2,2,200,-2.4,2.4);
  TH2F *h_trk_eta_vs_d0_np_allCuts = new TH2F("h_trk_eta_vs_d0_np_allCuts","h_trk_eta_vs_d0_np_allCuts; Transverse Impact Parameter d_{0} (cm); Eta",200,-2,2,200,-2.4,2.4);
  TH2F *h_trk_numStubs_vs_eta_primary_allCuts = new TH2F("h_trk_numStubs_vs_eta_primary_allCuts","h_trk_numStubs_vs_eta_primary_allCuts; Eta; Track Stubs",200,-2.4,2.4,7,0,7);
  TH2F *h_trk_numStubs_vs_eta_np_allCuts = new TH2F("h_trk_numStubs_vs_eta_np_allCuts","h_trk_numStubs_vs_eta_np_allCuts; Eta; Track Stubs",200,-2.4,2.4,7,0,7);
  TH2F *h_trk_numStubs_vs_eta_fake_allCuts = new TH2F("h_trk_numStubs_vs_eta_fake_allCuts","h_trk_numStubs_vs_eta_fake_allCuts; Eta; Track Stubs",200,-2.4,2.4,7,0,7);
  TH2F *h_trk_numStubs_vs_eta_PU_allCuts = new TH2F("h_trk_numStubs_vs_eta_PU_allCuts","h_trk_numStubs_vs_eta_PU_allCuts; Eta; Track Stubs",200,-2.4,2.4,7,0,7);
  TH2F *h_trk_numStubs_vs_eta_notHiggs_allCuts = new TH2F("h_trk_numStubs_vs_eta_notHiggs_allCuts","h_trk_numStubs_vs_eta_notHiggs_allCuts; Eta; Track Stubs",200,-2.4,2.4,7,0,7);
  TH2F *h_trueVertex_charge_vs_numTPs = new TH2F("h_trueVertex_charge_vs_numTPs","h_trueVertex_charge_vs_numTPs; TPs Associated with Vertex; Net Charge",6,0,6,12,-6,6);
  TH2F *h_trueVertex_charge_vs_numTracks = new TH2F("h_trueVertex_charge_vs_numTracks","h_trueVertex_charge_vs_numTracks; Tracks Associated with Vertex; Net Charge",6,0,6,12,-6,6);
  TH2F *h_correct_trackVertex_charge_vs_numTracks = new TH2F("h_correct_trackVertex_charge_vs_numTracks","h_correct_trackVertex_charge_vs_numTracks; Tracks Associated with Vertex; Net Charge",20,0,20,40,-20,20);
  TH2F *h_false_trackVertex_charge_vs_numTracks = new TH2F("h_false_trackVertex_charge_vs_numTracks","h_false_trackVertex_charge_vs_numTracks; Tracks Associated with Vertex; Net Charge",20,0,20,40,-20,20);
  TH2F *h_trackVertex_rankPt_vs_numVertPerTrack = new TH2F("h_trackVertex_rankPt_vs_numVertPerTrack","h_trackVertex_rankPt_vs_numVertPerTrack; Vertices Associated with Track; Ranking of Lower-p_{T} Track / 1.0",20,0,20,20,0,20);
  TH2F *h_trackVertex_rankDelz_vs_numVertPerTrack = new TH2F("h_trackVertex_rankDelz_vs_numVertPerTrack","h_trackVertex_rankDelz_vs_numVertPerTrack; Vertices Associated with Track; Ranking of Vertex #Deltaz / 1.0",20,0,20,20,0,20);
  TH2F *h_trackVertex_rankDt_vs_numVertPerTrack = new TH2F("h_trackVertex_rankDt_vs_numVertPerTrack","h_trackVertex_rankDt_vs_numVertPerTrack; Vertices Associated with Track; Ranking of Vertex d_{T} / 1.0",20,0,20,20,0,20);
  TH2F *h_trackVertex_rankChi2rphidofSum_vs_numVertPerTrack = new TH2F("h_trackVertex_rankChi2rphidofSum_vs_numVertPerTrack","h_trackVertex_rankChi2rphidofSum_vs_numVertPerTrack; Vertices Associated with Track; Ranking of Vertex #Sigma #chi^{2}_{r#phi}/d.o.f / 1.0",20,0,20,20,0,20);
  TH2F *h_trackVertex_rankRt_vs_numVertPerTrack = new TH2F("h_trackVertex_rankRt_vs_numVertPerTrack","h_trackVertex_rankRt_vs_numVertPerTrack; Vertices Associated with Track; Ranking of Vertex R_{T} / 1.0",20,0,20,20,0,20);
  
  std::string binVariable = "";
  std::vector<std::vector<double>> track_bins = {{-1,1}};
  std::vector<std::vector<double>> z0_bins;
  double z0_bin_min = -20.0;
  double z0_bin_width = 4.0;
  double z0_bin_max = z0_bin_min+z0_bin_width;
  while(z0_bin_max<20.0){
    z0_bins.push_back({z0_bin_min,z0_bin_max});
    z0_bin_min += (z0_bin_width/2);
    z0_bin_max += (z0_bin_width/2);
  }
  std::vector<std::vector<double>> phi_bins;
  double phi_bin_min = -TMath::Pi();
  double phi_bin_width = 0.6;
  double phi_bin_max = phi_bin_min+phi_bin_width;
  while(phi_bin_max<TMath::Pi()){
    phi_bins.push_back({phi_bin_min,phi_bin_max});
    phi_bin_min += (phi_bin_width/2);
    phi_bin_max += (phi_bin_width/2);
  }
  if(binVariable=="z0") track_bins = z0_bins;
  if(binVariable=="phi") track_bins = phi_bins;

  std::map<string,int> numPart_primary_noCuts{};
  std::map<string,int> numPart_primary_chi2rzdofCuts{};
  std::map<string,int> numPart_primary_bendchi2Cuts{};
  std::map<string,int> numPart_primary_chi2rphidofCuts{};
  std::map<string,int> numPart_primary_nstubCuts{};
  std::map<string,int> numPart_primary_ptCuts{};
  std::map<string,int> numPart_primary_d0Cuts{};
  std::map<string,int> numPart_primary_z0Cuts{};
  std::map<string,int> numPart_np_noCuts{};
  std::map<string,int> numPart_np_chi2rzdofCuts{};
  std::map<string,int> numPart_np_bendchi2Cuts{};
  std::map<string,int> numPart_np_chi2rphidofCuts{};
  std::map<string,int> numPart_np_nstubCuts{};
  std::map<string,int> numPart_np_ptCuts{};
  std::map<string,int> numPart_np_d0Cuts{};
  std::map<string,int> numPart_np_z0Cuts{};

  
  float pt_cuts[5] = {2.0,4.0,6.0,8.0,10.0};     // Cuts to control event rate
  float d0_cuts[5] = {0.1,0.2,0.3,0.4,0.5}; 
  std::vector<double> dxy_cuts = linspace(0.0,0.0000001,5); //cut <dxy
  std::vector<double> dz_cuts = linspace(0.0,3.0,5); //cut <dz
  std::vector<double> cos_T_cuts = linspace(.95,.999,5); //cut >cos
  std::vector<double> d_T_cuts = linspace(0.0,0.1,5); //cut <d_T
  std::vector<double> R_T_cuts = linspace(0.0,2.0,5); //cut >R_T
  std::vector<double> chi2rz_cuts = linspace(0.0,10.0,5); //cut <chi2rz
  std::vector<double> minD0_cuts = linspace(0.0,1.0,5); //cut >minD0
  std::vector<double> stubSum_cuts = linspace(8.0,12.0,5); //cut >stubSum
  
  
  double true_vertices = 0.0;
  
  std::vector<double> correct_vert_dxy_cut(dxy_cuts.size()); 
  std::vector<double> false_vert_dxy_cut(dxy_cuts.size());
  std::vector<double> all_vert_dxy_cut(dxy_cuts.size());
  
  std::vector<double> correct_vert_dz_cut(dz_cuts.size());
  std::vector<double> false_vert_dz_cut(dz_cuts.size());
  std::vector<double> all_vert_dz_cut(dz_cuts.size());
  std::vector<double> prev_dz_cut(dz_cuts.size());
  
  std::vector<double> correct_vert_cos_T_cut(cos_T_cuts.size());
  std::vector<double> false_vert_cos_T_cut(cos_T_cuts.size());
  std::vector<double> all_vert_cos_T_cut(cos_T_cuts.size());
  std::vector<double> prev_cos_T_cut(cos_T_cuts.size());
  
  std::vector<double> correct_vert_d_T_cut(d_T_cuts.size());
  std::vector<double> false_vert_d_T_cut(d_T_cuts.size());
  std::vector<double> all_vert_d_T_cut(d_T_cuts.size());
  std::vector<double> prev_d_T_cut(d_T_cuts.size());
  
  std::vector<double> correct_vert_R_T_cut(R_T_cuts.size());
  std::vector<double> false_vert_R_T_cut(R_T_cuts.size());
  std::vector<double> all_vert_R_T_cut(R_T_cuts.size());
  std::vector<double> prev_R_T_cut(R_T_cuts.size());
  
  std::vector<double> correct_vert_chi2rz_cut(chi2rz_cuts.size());
  std::vector<double> false_vert_chi2rz_cut(chi2rz_cuts.size());
  std::vector<double> all_vert_chi2rz_cut(chi2rz_cuts.size());
  std::vector<double> prev_chi2rz_cut(chi2rz_cuts.size());
  
  std::vector<double> correct_vert_minD0_cut(minD0_cuts.size());
  std::vector<double> false_vert_minD0_cut(minD0_cuts.size());
  std::vector<double> all_vert_minD0_cut(minD0_cuts.size());
  std::vector<double> prev_minD0_cut(minD0_cuts.size());
  
  std::vector<double> correct_vert_stubSum_cut(stubSum_cuts.size());
  std::vector<double> false_vert_stubSum_cut(stubSum_cuts.size());
  std::vector<double> all_vert_stubSum_cut(stubSum_cuts.size());
  std::vector<double> prev_stubSum_cut(stubSum_cuts.size());
  
  double correct_vert_array_cut[dz_cuts.size()][cos_T_cuts.size()][d_T_cuts.size()][R_T_cuts.size()][chi2rz_cuts.size()][minD0_cuts.size()][stubSum_cuts.size()];
  double false_vert_array_cut[dz_cuts.size()][cos_T_cuts.size()][d_T_cuts.size()][R_T_cuts.size()][chi2rz_cuts.size()][minD0_cuts.size()][stubSum_cuts.size()];
  double all_vert_array_cut[dz_cuts.size()][cos_T_cuts.size()][d_T_cuts.size()][R_T_cuts.size()][chi2rz_cuts.size()][minD0_cuts.size()][stubSum_cuts.size()];
  double prev_array_cut[dz_cuts.size()][cos_T_cuts.size()][d_T_cuts.size()][R_T_cuts.size()][chi2rz_cuts.size()][minD0_cuts.size()][stubSum_cuts.size()];
  
  if(detailedPlots){
    for(uint k_a=0;k_a<dz_cuts.size();k_a++){
      for(uint k_b=0;k_b<cos_T_cuts.size();k_b++){
	for(uint k_c=0;k_c<d_T_cuts.size();k_c++){
	  for(uint k_d=0;k_d<R_T_cuts.size();k_d++){
	    for(uint k_e=0;k_e<chi2rz_cuts.size();k_e++){
	      for(uint k_f=0;k_f<minD0_cuts.size();k_f++){
		for(uint k_g=0;k_g<stubSum_cuts.size();k_g++){
		  correct_vert_array_cut[k_a][k_b][k_c][k_d][k_e][k_f][k_g] = 0;
		  false_vert_array_cut[k_a][k_b][k_c][k_d][k_e][k_f][k_g] = 0;
		  all_vert_array_cut[k_a][k_b][k_c][k_d][k_e][k_f][k_g] = 0;
		  prev_array_cut[k_a][k_b][k_c][k_d][k_e][k_f][k_g] = 0;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  if (tree == 0) return;
  Long64_t nevt = tree->GetEntries();
  nevt = 1;
  Vertex_Parameters geomTrackVertex;
  Vertex_Parameters geomTrueVertex;

  for (Long64_t i_evnt=0; i_evnt<nevt; i_evnt++) {
    //std::cout<<"event number: "<<i_evnt<<std::endl;
    tree->GetEntry(i_evnt);
    displayProgress(i_evnt, nevt);
    std::vector<std::deque<Track_Parameters>> binnedSelectedTracks;
    for(uint i=0; i<track_bins.size(); i++){
      binnedSelectedTracks.push_back({});
    }
    std::deque<Track_Parameters> selectedTracks;      // Tracks 
    std::deque<Track_Parameters> selectedTPs;         // Tracking particles
    std::vector<Vertex_Parameters> trueVertices;
    std::vector<Vertex_Parameters> trackVertices;
    std::vector<Vertex_Parameters> tempVertices;
    int maxPT_i = 0;
    int firstMatch_j = 0;
    bool oneMatch = false;
    std::valarray<float> trkMET = {0.0,0.0};
    float trkH_T = 0.0;
    std::vector<float> matchtp_pt_noCuts_vec;
    std::vector<float> matchtp_pt_oldCuts_vec;
    std::vector<float> matchtp_pt_allCuts_vec;
    std::vector<float> tp_pt_noCuts_primary_vec;
    std::vector<float> tp_pt_oldCuts_vec;
    std::vector<float> tp_pt_allCuts_vec;
    // ----------------------------------------------------------------------------------------------------------------
    // track loop
    for (int it = 0; it < (int)trk_pt->size(); it++){
      bool isPrimary = true;
      if(inputFile.Contains("DarkPhoton")) isPrimary = trk_matchtp_isHToMu->at(it);
      if(inputFile.Contains("DisplacedTrackJet")) isPrimary = trk_matchtp_isHToB->at(it);
      //std::cout<<"track pt: "<<trk_pt->at(it)<<" eta: "<<trk_eta->at(it)<<" d0 : "<<trk_d0->at(it)<<" phi: "<<trk_phi->at(it)<<" z0: "<<trk_z0->at(it)<<" nstub: "<<trk_nstub<<std::endl;
      uint icut_counter = 0;
      for(auto icut=preselCuts.cbegin(); icut!=preselCuts.cend(); ++icut){
	bool mods = true;
	float param;
	if(icut->second.first.type() == typeid(std::vector<float>**)){
	  std::vector<float>** paramVector = boost::get<std::vector<float>**>(icut->second.first);
	  param = (*paramVector)->at(it);
	}
	else{
	  std::vector<int>** paramVector = boost::get<std::vector<int>**>(icut->second.first);
	  int tmp_param = (*paramVector)->at(it);
	  param = float(tmp_param);
	}
	TString cutName = icut->first;
	float cutValue = icut->second.second;
	if(cutName.Contains("D0") || cutName.Contains("Eta")) param = fabs(param);
	//std::cout<<"cutName: "<<cutName<<" cutValue: "<<cutValue<<" param: "<<param<<std::endl;
	if(cutName.Contains("barrel") && fabs(trk_eta->at(it))>barrelEta) mods = false;
	if(cutName.Contains("disk") && fabs(trk_eta->at(it))<=barrelEta) mods = false;
	if(cutName.Contains("_D") && fabs(trk_d0->at(it))<=1 ) mods = false;
	if(cutName.Contains("overlap") && (fabs(trk_eta->at(it))<=1.1 || fabs(trk_eta->at(it))>=1.7)) mods = false;
	if(cutName.Contains("dof")) param /= 2*trk_nstub->at(it) - 5;
	if(mods){
	  if(cutName.Contains("max") && param>cutValue) break;
	  if(cutName.Contains("min") && param<cutValue) break;
	}
	//std::cout<<"passed cut"<<std::endl;
	for(uint i=0; i<trackType.size(); ++i){
	  bool primary = trk_fake->at(it)==1 && isPrimary;
	  if(trackType[i]=="primary" && !primary) continue;
	  if(trackType[i]=="np" && primary) continue;
	  if(trackType[i]=="fake" && trk_fake->at(it)!=0) continue;
	  if(trackType[i]=="PU" && trk_fake->at(it)!=2) continue;
	  if(trackType[i]=="notHiggs" && !(trk_fake->at(it)==1 && !isPrimary)) continue;
	  string partId = to_string(trk_matchtp_pdgid->at(it));
	  numPartCutFlows[i][icut_counter][partId]++;
	  for(uint j=0; j<plotModifiers.size(); ++j){
	    if(plotModifiers[j]=="_H" && trk_pt->at(it)<=10) continue;
	    if(plotModifiers[j]=="_L" && trk_pt->at(it)>10) continue;
	    if(plotModifiers[j]=="_P" && fabs(trk_d0->at(it))>1) continue;
	    if(plotModifiers[j]=="_D" && fabs(trk_d0->at(it))<=1) continue;
	    if(plotModifiers[j]=="_barrel" && fabs(trk_eta->at(it))>barrelEta) continue;
	    if(plotModifiers[j]=="_disk" && fabs(trk_eta->at(it))<=barrelEta) continue;
	    int ivar_counter = 0;
	    for(auto ivar=varCutFlows.cbegin(); ivar!=varCutFlows.cend(); ++ivar){
	      if(ivar->second.first.type() == typeid(std::vector<float>**)){
		std::vector<float>** paramVector = boost::get<std::vector<float>**>(ivar->second.first);
		param = (*paramVector)->at(it);
	      }
	      else{
		std::vector<int>** paramVector = boost::get<std::vector<int>**>(ivar->second.first);
		int tmp_param = (*paramVector)->at(it);
		param = float(tmp_param);
	      }
	      TString varName = ivar->first.at(0);
	      if(varName.Contains("dof")) param /= 2*trk_nstub->at(it) - 5;
	      if(varName.Contains("sector")){
		while (param < -TMath::Pi()/9 ) param += 2*TMath::Pi();
		while (param > TMath::Pi()*2 ) param -= 2*TMath::Pi();
		while (param > TMath::Pi()/9) param -= 2*TMath::Pi()/9;
	      }
	      preselCutFlows[ivar_counter][i][icut_counter][j]->Fill(param);
	      ivar_counter++;
	    }
	    int ivar2D_counter = 0;
	    for(auto ivar2D=varCutFlows2D.cbegin(); ivar2D!=varCutFlows2D.cend(); ++ivar2D){
	      std::vector<float> params;
	      for(int k=0; k<2; ++k){
		if(ivar2D->second.first.at(k).type() == typeid(std::vector<float>**)){
		  std::vector<float>** paramVector = boost::get<std::vector<float>**>(ivar2D->second.first.at(k));
		  param = (*paramVector)->at(it);
		}
		else{
		  std::vector<int>** paramVector = boost::get<std::vector<int>**>(ivar2D->second.first.at(k));
		  int tmp_param = (*paramVector)->at(it);
		  param = float(tmp_param);
		}
		TString varName = ivar2D->first.at(2*k);
		if(varName.Contains("dof")) param /= 2*trk_nstub->at(it) - 5;
		if(varName.Contains("sector")){
		  while (param < -TMath::Pi()/9 ) param += 2*TMath::Pi();
		  while (param > TMath::Pi()*2 ) param -= 2*TMath::Pi();
		  while (param > TMath::Pi()/9) param -= 2*TMath::Pi()/9;
		}
		params.push_back(param);
	      }
	      preselCutFlows2D[ivar2D_counter][i][icut_counter][j]->Fill(params[0],params[1]);
	      ivar2D_counter++;
	    }
	  }
	}
	icut_counter++;
      }
      if(icut_counter==preselCuts.size()){
	Track_Parameters* tp_params = new Track_Parameters(trk_matchtp_pt->at(it), -1*trk_matchtp_d0->at(it), trk_matchtp_z0->at(it), trk_matchtp_eta->at(it), trk_matchtp_phi->at(it), trk_matchtp_pdgid->at(it), trk_matchtp_x->at(it), trk_matchtp_y->at(it), trk_matchtp_z->at(it));
	for(uint i=0; i<track_bins.size(); i++){
	  float trkVariable = 0.0;
	  if(binVariable=="phi") trkVariable = trk_phi->at(it);
	  if (binVariable=="z0") trkVariable = fabs(trk_z0->at(it));
	  if(trkVariable<track_bins[i][1] && trkVariable>track_bins[i][0] ){
	    binnedSelectedTracks[i].push_back(Track_Parameters(trk_pt->at(it), trk_d0->at(it), trk_z0->at(it), trk_eta->at(it), trk_phi->at(it), -99999, -999., -999., -999., trk_rinv->at(it), it, tp_params, trk_nstub->at(it), trk_chi2rphi->at(it), trk_chi2rz->at(it), trk_bendchi2->at(it), trk_MVA1->at(it), trk_MVA2->at(it)));
	  }
	}
	
	//std::cout<<"pushed track icut_counter: "<<icut_counter<<std::endl;
	selectedTracks.push_back(Track_Parameters(trk_pt->at(it), trk_d0->at(it), trk_z0->at(it), trk_eta->at(it), trk_phi->at(it), -99999, -999., -999., -999., trk_rinv->at(it), it, tp_params, trk_nstub->at(it), trk_chi2rphi->at(it), trk_chi2rz->at(it), trk_bendchi2->at(it), trk_MVA1->at(it), trk_MVA2->at(it) ));
	trkH_T += trk_pt->at(it);
	std::valarray<float> trackPtVec = {trk_pt->at(it)*cos(trk_phi->at(it)),trk_pt->at(it)*sin(trk_phi->at(it))};
	trkMET -= trackPtVec;
      }
    }
    
    h_trk_H_T->Fill(trkH_T);
    h_trk_MET->Fill(TMath::Sqrt(pow(trkMET[0],2)+pow(trkMET[1],2)));
#if 0
    //pt isolation calculation
    std::vector<int> deleteTracks;
    for(uint i=0; i<selectedTracks.size(); i++){
      float ptSum4 = 0.0;
      float ptSum8 = 0.0;
      for(uint j=0; j<selectedTracks.size(); j++){
	if(i==j) continue;
	float deltaR = TMath::Sqrt(pow(fabs(selectedTracks[i].phi-selectedTracks[j].phi),2)+fabs(selectedTracks[i].eta-selectedTracks[j].eta));
	if(deltaR<0.8){
	  ptSum8 += selectedTracks[j].pt;
	  if(deltaR<0.4){
	    ptSum4 += selectedTracks[j].pt;
	  }
	}
      }
      float ptIso4 = ptSum4 / selectedTracks[i].pt;
      float ptIso8 = ptSum8 / selectedTracks[i].pt;
      bool isPrimary = true;
      if(inputFile.Contains("DarkPhoton")) isPrimary = trk_matchtp_isHToMu->at(selectedTracks[i].index);
      if(inputFile.Contains("DisplacedTrackJet")) isPrimary = trk_matchtp_isHToB->at(selectedTracks[i].index);
      if(trk_fake->at(selectedTracks[i].index)==1 && isPrimary){
	h_trk_ptIso4_primary_allCuts->Fill(ptIso4);
	h_trk_ptIso8_primary_allCuts->Fill(ptIso8);
      }
      else{
	h_trk_ptIso4_np_allCuts->Fill(ptIso4);
	h_trk_ptIso8_np_allCuts->Fill(ptIso8);
      }
      //if(ptIso4<0.05) deleteTracks.push_back(selectedTracks[i].index);
    }

    for(uint i=0; i<deleteTracks.size(); i++){
      for(uint j=0; j<selectedTracks.size();){
	if(selectedTracks[j].index==deleteTracks[i]){
	  selectedTracks.erase(selectedTracks.begin()+j); 
	}
	else{
	  j++;
	}
      }
      for(uint j=0; j<binnedSelectedTracks.size(); j++){
	for(uint k=0; k<binnedSelectedTracks[j].size();){
	  if(binnedSelectedTracks[j][k].index==deleteTracks[i]){
	    binnedSelectedTracks[j].erase(binnedSelectedTracks[j].begin()+k);
	  }
	  else{
	    k++;
	  }
	}
      }
    }
#endif
    h_numSelectedTrks->Fill(selectedTracks.size());
    h_numSelectedTrks_zoomOut->Fill(selectedTracks.size());
    
    // ----------------------------------------------------------------------------------------------------------------
    // tracking particle loop
    float tpH_T = 0.0;
    std::valarray<float> tpMET = {0.0,0.0};
    //std::cout<<"tp_pt size: "<<tp_pt->size()<<std::endl;
    for (int it = 0; it < (int)tp_pt->size(); it++){
      
      /*if(tp_isHToMu->at(it)==true && tp_eventid->at(it)==0){
	std::cout<<"tp pos: "<<tp_x->at(it)<<" "<<tp_y->at(it)<<" "<<tp_z->at(it)<<" pdgid: "<<tp_pdgid->at(it)<<" pt: "<<tp_pt->at(it)<<std::endl;
	}*/
      
      float tmp_d0 = -tp_d0->at(it);	// Sign difference in the NTupleMaker
      float tmp_z0 = tp_z0->at(it);
	
      bool isPrimary = true;
      if(inputFile.Contains("DarkPhoton")) isPrimary = tp_isHToMu->at(it);
      if(inputFile.Contains("DisplacedTrackJet")) isPrimary = tp_isHToB->at(it);
      //std::cout<<"tp pt: "<<tp_pt->at(it)<<" eta: "<<tp_eta->at(it)<<" d0 : "<<tp_d0->at(it)<<" phi: "<<tp_phi->at(it)<<" z0: "<<tp_z0->at(it)<<" nstub: "<<tp_nstub<<std::endl;
      uint icut_counter = 0;
      for(auto icut=preselCutsTP.cbegin(); icut!=preselCutsTP.cend(); ++icut){
	bool mods = true;
	float param;
	if(icut->second.first.type() == typeid(std::vector<float>**)){
	  std::vector<float>** paramVector = boost::get<std::vector<float>**>(icut->second.first);
	  param = (*paramVector)->at(it);
	}
	else{
	  std::vector<int>** paramVector = boost::get<std::vector<int>**>(icut->second.first);
	  int tmp_param = (*paramVector)->at(it);
	  param = float(tmp_param);
	}
	TString cutName = icut->first;
	if(cutName.Contains("D0") || cutName.Contains("Eta")) param = fabs(param);
	float cutValue = icut->second.second;
	//std::cout<<"cutName: "<<cutName<<" cutValue: "<<cutValue<<" param: "<<param<<std::endl;
	if(cutName.Contains("barrel") && fabs(tp_eta->at(it))>barrelEta) mods = false;
	if(cutName.Contains("disk") && fabs(tp_eta->at(it))<=barrelEta) mods = false;
	if(cutName.Contains("_D") && fabs(tp_d0->at(it))<=1 ) mods = false;
	if(cutName.Contains("overlap") && (fabs(tp_eta->at(it))<=1.1 || fabs(tp_eta->at(it))>=1.7)) mods = false;
	if(cutName.Contains("dof")) param /= 2*tp_nstub->at(it) - 5;
	if(mods){
	  if(cutName.Contains("max") && param>cutValue) break;
	  if(cutName.Contains("min") && param<cutValue) break;
	}
	//std::cout<<"passed cut"<<std::endl;
	for(uint i=0; i<tpType.size(); ++i){
	  bool primary = tp_eventid->at(it)==0 && isPrimary;
	  if(tpType[i]=="primary" && !primary) continue;
	  if(tpType[i]=="np" && primary) continue;
	  if(tpType[i]=="PU" && tp_eventid->at(it)==0) continue;
	  if(tpType[i]=="notHiggs" && !(tp_eventid->at(it)==0 && !isPrimary)) continue;
	  if(tpType[i]=="match" && tp_nmatch->at(it)==0) continue;
	  string partId = to_string(tp_pdgid->at(it));
	  numPartCutFlowsTP[i][icut_counter][partId]++;
	  for(uint j=0; j<plotModifiers.size(); ++j){
	    if(plotModifiers[j]=="_H" && tp_pt->at(it)<=10) continue;
	    if(plotModifiers[j]=="_L" && tp_pt->at(it)>10) continue;
	    if(plotModifiers[j]=="_P" && fabs(tp_d0->at(it))>1) continue;
	    if(plotModifiers[j]=="_D" && fabs(tp_d0->at(it))<=1) continue;
	    if(plotModifiers[j]=="_barrel" && fabs(tp_eta->at(it))>barrelEta) continue;
	    if(plotModifiers[j]=="_disk" && fabs(tp_eta->at(it))<=barrelEta) continue;
	    int ivar_counter = 0;
	    for(auto ivar=varCutFlowsTP.cbegin(); ivar!=varCutFlowsTP.cend(); ++ivar){
	      if(ivar->second.first.type() == typeid(std::vector<float>**)){
		std::vector<float>** paramVector = boost::get<std::vector<float>**>(ivar->second.first);
		param = (*paramVector)->at(it);
	      }
	      else{
		std::vector<int>** paramVector = boost::get<std::vector<int>**>(ivar->second.first);
		int tmp_param = (*paramVector)->at(it);
		param = float(tmp_param);
	      }
	      TString varName = ivar->first.at(0);
	      if(varName.Contains("dof")) param /= 2*tp_nstub->at(it) - 5;
	      if(varName.Contains("sector")){
		while (param < -TMath::Pi()/9 ) param += 2*TMath::Pi();
		while (param > TMath::Pi()*2 ) param -= 2*TMath::Pi();
		while (param > TMath::Pi()/9) param -= 2*TMath::Pi()/9;
	      }
	      preselCutFlowsTP[ivar_counter][i][icut_counter][j]->Fill(param);
	      ivar_counter++;
	    }
	    int ivar2D_counter = 0;
	    for(auto ivar2D=varCutFlowsTP2D.cbegin(); ivar2D!=varCutFlowsTP2D.cend(); ++ivar2D){
	      std::vector<float> params;
	      for(int k=0; k<2; ++k){
		if(ivar2D->second.first.at(k).type() == typeid(std::vector<float>**)){
		  std::vector<float>** paramVector = boost::get<std::vector<float>**>(ivar2D->second.first.at(k));
		  param = (*paramVector)->at(it);
		}
		else{
		  std::vector<int>** paramVector = boost::get<std::vector<int>**>(ivar2D->second.first.at(k));
		  int tmp_param = (*paramVector)->at(it);
		  param = float(tmp_param);
		}
		TString varName = ivar2D->first.at(2*k);
		if(varName.Contains("dof")) param /= 2*tp_nstub->at(it) - 5;
		if(varName.Contains("sector")){
		  while (param < -TMath::Pi()/9 ) param += 2*TMath::Pi();
		  while (param > TMath::Pi()*2 ) param -= 2*TMath::Pi();
		  while (param > TMath::Pi()/9) param -= 2*TMath::Pi()/9;
		}
		params.push_back(param);
	      }
	      preselCutFlowsTP2D[ivar2D_counter][i][icut_counter][j]->Fill(params[0],params[1]);
	      ivar2D_counter++;
	    }
	  }
	}
	icut_counter++;
      }
      if(icut_counter==preselCutsTP.size()){
	//std::cout<<"pushed tp icut_counter: "<<icut_counter<<std::endl;
	selectedTPs.push_back(Track_Parameters(tp_pt->at(it), tmp_d0, tmp_z0, tp_eta->at(it), tp_phi->at(it), tp_pdgid->at(it), tp_x->at(it), tp_y->at(it), tp_z->at(it), tp_charge->at(it), it));
	if (tp_eventid->at(it)>0){
	  tpH_T += tp_pt->at(it);
	  std::valarray<float> tpPtVec = {tp_pt->at(it)*cos(tp_phi->at(it)),tp_pt->at(it)*sin(tp_phi->at(it))};
	  tpMET -= tpPtVec;
	}
      }
    }
    h_tp_H_T->Fill(tpH_T);
    h_tp_MET->Fill(TMath::Sqrt(pow(tpMET[0],2)+pow(tpMET[1],2)));
    
    // --------------------------------------------------------------------------------------------
    //         Vertex finding in Tracking Particles
    // --------------------------------------------------------------------------------------------
    if (!(selectedTracks.size() >= 2)) continue;
    bool true_DV = false;
    if(VERBOSE[0])
      std::cout<<"End of z-vertex Finding"<<endl;
      
    double_t x_dv = -9999.0;// (tp_x->at((*selectedTPs)[0]->index));//+tp_x->at((*selectedTPs)[1]->index))/2.0;
    double_t y_dv = -9999.0;// (tp_y->at((*selectedTPs)[0]->index));//+tp_y->at((*selectedTPs)[1]->index))/2.0;
    double_t z_dv = -9999.0;// (tp_z->at((*selectedTPs)[0]->index));//+tp_z->at((*selectedTPs)[1]->index))/2.0;

    double_t x_tmp = x_dv;
    double_t y_tmp = y_dv;
    double_t z_tmp = z_dv;
    Int_t Vertex_check_tmp = -1;
    int noCuts_trueVertices = 0;
    
    if(selectedTPs.size()>=2){
      //std::cout<<"selectedTPs size: "<<selectedTPs.size()<<std::endl;
      sort(selectedTPs.begin(), selectedTPs.end(), ComparePtTrack);
      std::vector<Track_Parameters> copyTPs;
      if(detailedPlots){
	for (uint i=0; i<selectedTPs.size(); i++){
	  copyTPs.push_back(selectedTPs[i]);
	}
      }
      /*std::cout<<"Selected TPs"<<std::endl;
      for( uint i=0; i<selectedTPs.size(); i++ ){
	std::cout<<"tp pos: "<<tp_x->at(selectedTPs[i].index)<<" "<<tp_y->at(selectedTPs[i].index)<<" "<<tp_z->at(selectedTPs[i].index)<<" pdgid: "<<tp_pdgid->at(selectedTPs[i].index)<<" pt: "<<tp_pt->at(selectedTPs[i].index)<<std::endl;
	}*/
      
      while(selectedTPs.size()>1){
	bool foundTrueVertex = false;
	for( uint i=1; i<selectedTPs.size();){
	  int index0 = selectedTPs[0].index;
	  int index1 = selectedTPs[i].index;
	  if( fabs(tp_x->at(index0)-tp_x->at(index1))<0.0001 && fabs(tp_y->at(index0)-tp_y->at(index1))<0.0001 && fabs(tp_z->at(index0)-tp_z->at(index1))<0.0001 ){
	    //std::cout<<"tp1 pos: "<<tp_x->at(index0)<<" "<<tp_y->at(index0)<<" "<<tp_z->at(index0)<<" pdgid: "<<tp_pdgid->at(index0)<<" tp2 pos: "<<tp_x->at(index1)<<" "<<tp_y->at(index1)<<" "<<tp_z->at(index1)<<" pdgid: "<<tp_pdgid->at(index1)<<std::endl;
	    x_dv = tp_x->at(index0);
	    y_dv = tp_y->at(index0);
	    z_dv = tp_z->at(index0);
	    noCuts_trueVertices++;
	    if(dist(x_dv,y_dv)>d0_res && dist(x_dv,y_dv)<20){
	      //std::cout<<"true vertex: "<<x_dv<<" "<<y_dv<<" "<<z_dv<<" tp_pt: "<<selectedTPs[0].pt<<" "<<selectedTPs[i].pt<<" eventid's: "<<tp_eventid->at(selectedTPs[0].index)<<" "<<tp_eventid->at(selectedTPs[i].index)<<std::endl;
	      if(!foundTrueVertex){
		h_trueVertex_x->Fill(x_dv);
		h_trueVertex_y->Fill(y_dv);
		h_trueVertex_z->Fill(z_dv);
		h_trueVertex_sumPt->Fill(selectedTPs[0].pt+selectedTPs[i].pt);
		h_trueVertex_lowPt->Fill(selectedTPs[i].pt);
		h_trueVertex_highPt->Fill(selectedTPs[0].pt);
		trueVertices.push_back(Vertex_Parameters(x_dv, y_dv, z_dv, selectedTPs[0], selectedTPs[i]) );
		foundTrueVertex = true;
	      }
	      else{
		trueVertices.back().tracks.push_back(selectedTPs[i]);
	      }
	      selectedTPs.erase(selectedTPs.begin()+i);
	    }
	    else{
	      i++;
	    }
	  }
	  else{
	    i++;
	  }
	}
	selectedTPs.pop_front();
      }
      
      h_trueVertex_numNoCuts->Fill(noCuts_trueVertices);
      h_trueVertex_numAllCuts->Fill(trueVertices.size());
      true_vertices += trueVertices.size();
      float maxPT = 0.0;
      for(uint i=0; i<trueVertices.size(); i++){
	if(trueVertices[i].a.pt>maxPT){
	  maxPT = trueVertices[i].a.pt;
	  maxPT_i = i;
	}
	
	/*for(uint j=0; j<trueVertices[i].tracks.size(); j++){
	  std::cout<<"vertex "<<i<<" tp pos: "<<tp_x->at(trueVertices[i].tracks[j].index)<<" "<<tp_y->at(trueVertices[i].tracks[j].index)<<" "<<tp_z->at(trueVertices[i].tracks[j].index)<<" pdgid: "<<tp_pdgid->at(trueVertices[i].tracks[j].index)<<" pt: "<<tp_pt->at(trueVertices[i].tracks[j].index)<<std::endl;
	  }*/
	
	h_trueVertex_d_T->Fill(trueVertices[i].d_T);
	h_trueVertex_cos_T->Fill(trueVertices[i].cos_T);
	h_trueVertex_alpha_T->Fill(trueVertices[i].alpha_T);
	h_trueVertex_R_T->Fill(trueVertices[i].R_T);
	h_trueVertex_openingAngle->Fill(trueVertices[i].openingAngle);
	h_trueVertex_parentPt->Fill(trueVertices[i].p_mag);

	int i_pt = 9999;
	int j_pt = 9999;
	if(detailedPlots){
	  for( uint j=0; j<copyTPs.size();j++ ){
	    if(trueVertices[i].a.index == copyTPs[j].index){
	      i_pt = j;
	    }
	    if(trueVertices[i].b.index == copyTPs[j].index){
	      j_pt = j;
	    }
	  }	
	  
	  h_trueVertexCuts_indexPt->Fill(i_pt);
	  h_trueVertexCuts_indexPt->Fill(j_pt);
	  h_trueVertex_delta_dist_indexPt->Fill(abs(i_pt-j_pt));
	}
	h_delta_dist_xy->Fill(dist_TPs(trueVertices[i].a,trueVertices[i].b));
	Double_t x_dv_tp = -9999.0;
	Double_t y_dv_tp = -9999.0;
	Double_t z_dv_tp = -9999.0;
	int inTraj = calcVertex(trueVertices[i].a,trueVertices[i].b,x_dv_tp,y_dv_tp,z_dv_tp);
	h_trueVertex_inTraj->Fill(inTraj);
	h_delta_dist_z->Fill(fabs(trueVertices[i].a.z(x_dv_tp,y_dv_tp)-trueVertices[i].b.z(x_dv_tp,y_dv_tp)));
	h_trueVertex_delta_dist_z0->Fill(fabs(tp_z0->at(trueVertices[i].a.index)-tp_z0->at(trueVertices[i].b.index)));
	h_trueVertex_delta_dist_d0->Fill(fabs(tp_d0->at(trueVertices[i].a.index)-tp_d0->at(trueVertices[i].b.index)));
	float deltaEta = fabs(tp_eta->at(trueVertices[i].a.index)-tp_eta->at(trueVertices[i].b.index));
	h_trueVertex_delta_dist_eta->Fill(deltaEta);
	float deltaPhi = fabs(tp_phi->at(trueVertices[i].a.index)-tp_phi->at(trueVertices[i].b.index));
	h_trueVertex_delta_dist_phi->Fill(deltaPhi);
	float R = TMath::Sqrt(pow(deltaEta,2)+pow(deltaPhi,2));
	h_trueVertex_delta_dist_R->Fill(R);
	h_trueVertex_delta_dist_pt->Fill(fabs(tp_pt->at(trueVertices[i].a.index)-tp_pt->at(trueVertices[i].b.index)));
	
	if(trueVertices[i].x_dv!=-9999.0){
	  h_trueVertexCuts_x->Fill(trueVertices[i].x_dv);
	  h_trueVertexCuts_y->Fill(trueVertices[i].y_dv);
	  h_trueVertexCuts_z->Fill(trueVertices[i].z_dv);
	  h_trueVertexCuts_sumPt->Fill(trueVertices[i].a.pt+trueVertices[i].b.pt);
	}

	h_delta_x->Fill(fabs((tp_x->at(trueVertices[i].a.index) - (tp_x->at(trueVertices[i].b.index)))));
	true_DV = true;
	   
	h_all_trueVertex_pt->Fill(trueVertices[i].a.pt);
	if(fabs(trueVertices[i].a.d0) < fabs(trueVertices[i].b.d0)){
	  h_all_trueVertex_minD0->Fill(trueVertices[i].a.d0);
	  h_all_trueVertex_maxD0->Fill(trueVertices[i].b.d0);
	}
	else{
	  h_all_trueVertex_minD0->Fill(trueVertices[i].b.d0);
	  h_all_trueVertex_maxD0->Fill(trueVertices[i].a.d0);
	}
	
	float minD0 = fabs(trueVertices[i].tracks[0].d0);
	float lowPt = trueVertices[i].tracks[0].pt;
	float maxD0 = minD0;
	float netCharge = 0;
	for(uint itrack=0; itrack<trueVertices[i].tracks.size(); itrack++){
	  if(fabs(trueVertices[i].tracks[itrack].d0)<minD0) minD0 = fabs(trueVertices[i].tracks[itrack].d0);
	  if(fabs(trueVertices[i].tracks[itrack].d0)>maxD0) maxD0 = fabs(trueVertices[i].tracks[itrack].d0);
	  if(trueVertices[i].tracks[itrack].pt<lowPt) lowPt = trueVertices[i].tracks[itrack].pt;
	  netCharge+=trueVertices[i].tracks[itrack].charge;
	  float delz = fabs(trueVertices[i].tracks[itrack].z(trueVertices[i].x_dv,trueVertices[i].y_dv) - trueVertices[i].z_dv);
	  float delxy = dist_Vertex(trueVertices[i].x_dv,trueVertices[i].y_dv,trueVertices[i].tracks[itrack]);
	  h_trueVertexAssoc_delxy->Fill(delxy);
	  h_trueVertexAssoc_delz->Fill(delz);
	  if(itrack>1){
	    delz = fabs(trueVertices[i].tracks[itrack].z(x_dv_tp,y_dv_tp) - z_dv_tp);
	    delxy = dist_Vertex(x_dv_tp,y_dv_tp,trueVertices[i].tracks[itrack]);
	    h_trueVertexAssoc_calcVertDelxy->Fill(delxy);
	    h_trueVertexAssoc_calcVertDelz->Fill(delz);
	  }
	}

	h_all_trueVertex_minD0_allTPs->Fill(minD0);
	h_all_trueVertex_maxD0_allTPs->Fill(maxD0);
	h_all_trueVertex_lowPt->Fill(trueVertices[i].b.pt);
	h_all_trueVertex_lowPt_allTPs->Fill(lowPt);
	h_all_trueVertex_eta->Fill(trueVertices[i].a.eta);
	h_all_trueVertex_dxy->Fill(dist(trueVertices[i].x_dv,trueVertices[i].y_dv));
	h_trueVertex_numTPs->Fill(trueVertices[i].tracks.size());
	h_trueVertex_charge_vs_numTPs->Fill(trueVertices[i].tracks.size(),netCharge);
      }
      if(trueVertices.size()>0){
	h_all_oneMatch_trueVertex_pt->Fill(trueVertices[maxPT_i].a.pt);
	h_all_oneMatch_trueVertex_lowPt->Fill(trueVertices[maxPT_i].b.pt);
	h_all_oneMatch_trueVertex_eta->Fill(trueVertices[maxPT_i].a.eta);
	h_all_oneMatch_trueVertex_dxy->Fill(dist(trueVertices[maxPT_i].x_dv,trueVertices[maxPT_i].y_dv));
      }
    }
    
    // --------------------------------------------------------------------------------------------
    //                Vertex finding in Tracks
    // --------------------------------------------------------------------------------------------
    sort(selectedTracks.begin(), selectedTracks.end(), ComparePtTrack);
    std::vector<Track_Parameters> copyTracks;
    if(detailedPlots){
      for (uint i=0; i<selectedTracks.size(); i++){
	copyTracks.push_back(selectedTracks[i]);
	//std::cout<<"selected tracks tp pt: "<<selectedTracks[i].tp->pt<<std::endl;
      }
    }
    for(auto trackBin : binnedSelectedTracks){
      sort(trackBin.begin(), trackBin.end(), ComparePtTrack);
    }
    if(detailedPlots){
      for( uint i=0; i<trueVertices.size(); i++ ){
	int numMatched = 0;
	int netCharge = 0;
	std::vector<Track_Parameters> matchedTracks;
	for ( uint j=0; j<trueVertices[i].tracks.size(); j++ ){
	  for ( uint k=0; k<selectedTracks.size(); k++ ){
	    if(trueVertices[i].tracks[j]==selectedTracks[k].tp){
	      numMatched++;
	      netCharge+= selectedTracks[k].charge;
	      matchedTracks.push_back(selectedTracks[k]);
	    }
	  }
	}
	h_trueVertex_numTracks->Fill(numMatched);
	h_trueVertex_charge_vs_numTracks->Fill(numMatched,netCharge);
	if(numMatched>=2){
	  h_findable_trueVertex_pt->Fill(trueVertices[i].a.pt);
	  if(numMatched>=3){
	    sort(matchedTracks.begin(),matchedTracks.end(),ComparePtTrack);
	    uint itrack = 0;
	    bool foundVert = false;
	    while(itrack<matchedTracks.size()-1 && !foundVert){
	      for(uint jtrack=itrack+1; jtrack<matchedTracks.size(); jtrack++){
		if(dist_TPs(matchedTracks[itrack],matchedTracks[jtrack])==0){
		  Double_t x_dv_trk = -9999.0;
		  Double_t y_dv_trk = -9999.0;
		  Double_t z_dv_trk = -9999.0;
		  calcVertex(matchedTracks[itrack],matchedTracks[jtrack],x_dv_trk,y_dv_trk,z_dv_trk);
		  Vertex_Parameters matchedVert(x_dv_trk,y_dv_trk,z_dv_trk,matchedTracks[itrack],matchedTracks[jtrack]);
		  if(dist(x_dv_trk,y_dv_trk)>d0_res && dist(x_dv_trk,y_dv_trk)<20){
		    foundVert = true;
		    for(uint ktrack=0; ktrack<matchedTracks.size(); ktrack++){
		      if(ktrack==itrack || ktrack==jtrack) continue;
		      float delxy = dist_Vertex(x_dv_trk,y_dv_trk,matchedTracks[ktrack]);
		      float delz = fabs(matchedTracks[ktrack].z(x_dv_trk,y_dv_trk)-z_dv_trk);
		      h_trackVertexAssoc_delxy->Fill(delxy);
		      h_trackVertexAssoc_delz->Fill(delz);
		    }
		    break;
		  }
		}
	      }
	      itrack++; 
	    }
	  }
	}
      }
      
      for( uint i=0; i<trueVertices.size(); i++ ){
	//std::cout<<"trueVertex find eff loop i: "<<i<<std::endl;
	bool findable = false;
	bool intersect = false;
	bool nonPrompt = false;
	bool beforeTracker = false;
	bool bendChi2Cut = false;
	bool chi2RPhiCut = false;
	bool cosTCut = false;
	bool deltaZCut = false;
	std::vector<Vertex_Parameters> matchedVertices;
	for( auto trackBin : binnedSelectedTracks ){
	  int numMatched = 0;
	  std::vector<Track_Parameters> matchedTracks;
	  for( uint j=0; j<trueVertices[i].tracks.size(); j++ ){
	    for( uint k=0; k<trackBin.size(); k++){
	      if(trueVertices[i].tracks[j]==trackBin[k].tp){ 
		numMatched++;
		matchedTracks.push_back(trackBin[k]);
	      }
	    } 
	  }
	  if(numMatched>=2) findable = true;
	  if(matchedTracks.size()>=2){
	    for(uint i_trk=0; i_trk<matchedTracks.size()-1; i_trk++){
	      for(uint j_trk=i_trk+1; j_trk<matchedTracks.size(); j_trk++){
		if(dist_TPs(matchedTracks[i_trk],matchedTracks[j_trk])==0){
		  intersect = true;
		  Double_t x_dv_trk = -9999.0;
		  Double_t y_dv_trk = -9999.0;
		  Double_t z_dv_trk = -9999.0;
		  calcVertex(matchedTracks[i_trk],matchedTracks[j_trk],x_dv_trk,y_dv_trk,z_dv_trk);
		  Vertex_Parameters matchedVert(x_dv_trk,y_dv_trk,z_dv_trk,matchedTracks[i_trk],matchedTracks[j_trk]);
		  matchedVertices.push_back(matchedVert);
		  if(dist(x_dv_trk,y_dv_trk)>d0_res) nonPrompt = true;
		  if(dist(x_dv_trk,y_dv_trk)<20) beforeTracker = true;
		  if(matchedVert.bendchi2Sum<bendChi2Max) bendChi2Cut = true;
		  if(matchedVert.chi2rphidofSum<chi2RPhiMax) chi2RPhiCut = true;
		  if(matchedVert.cos_T>cosTMin) cosTCut = true;
		  if(matchedVert.delta_z<deltaZMax) deltaZCut = true;
		}
	      }
	    }
	  }
	}
	if(findable) h_findable_trueVertexBinned_pt->Fill(trueVertices[i].a.pt);
	if(findable&&intersect){
	  h_findableIntersect_trueVertexBinned_pt->Fill(trueVertices[i].a.pt);
	  for(uint i_vert=0; i_vert<matchedVertices.size(); i_vert++){
	    if(matchedVertices[i_vert]==matchedVertices[i_vert+1]) continue;
	    h_res_tp_trk_x_findVert->Fill(trueVertices[i].x_dv-matchedVertices[i_vert].x_dv);
	    h_res_tp_trk_y_findVert->Fill(trueVertices[i].y_dv-matchedVertices[i_vert].y_dv);
	  }
	}
	if(findable&&intersect&&nonPrompt) h_findableNonPrompt_trueVertexBinned_pt->Fill(trueVertices[i].a.pt);
	if(findable&&intersect&&nonPrompt&&beforeTracker) h_findableBeforeTracker_trueVertexBinned_pt->Fill(trueVertices[i].a.pt);
	if(findable&&intersect&&nonPrompt&&beforeTracker&&bendChi2Cut) h_findableBendChi2Cut_trueVertexBinned_pt->Fill(trueVertices[i].a.pt);
	if(findable&&intersect&&nonPrompt&&beforeTracker&&bendChi2Cut&&chi2RPhiCut) h_findableChi2RPhiCut_trueVertexBinned_pt->Fill(trueVertices[i].a.pt);
	if(findable&&intersect&&nonPrompt&&beforeTracker&&bendChi2Cut&&chi2RPhiCut&&cosTCut) h_findableCosTCut_trueVertexBinned_pt->Fill(trueVertices[i].a.pt);
	if(findable&&intersect&&nonPrompt&&beforeTracker&&bendChi2Cut&&chi2RPhiCut&&cosTCut&&deltaZCut) h_findableDeltaZCut_trueVertexBinned_pt->Fill(trueVertices[i].a.pt);
      }
      
      for(uint i=0; i<selectedTracks.size()-1; i++){
	for(uint j=i+1; j<selectedTracks.size(); j++){
	  bool isReal = false;
	  for(uint k=0; k<trueVertices.size(); k++){
	    int numMatched = 0;
	    for(uint l=0; l<trueVertices[k].tracks.size(); l++){
	      if(selectedTracks[i].tp==trueVertices[k].tracks[l] || selectedTracks[j].tp==trueVertices[k].tracks[l]) numMatched++;
	    }
	    if(numMatched>=2) isReal = true;
	  }
	  if(!isReal) h_findFake_trackVertex_pt->Fill(selectedTracks[i].pt);
	}
      }
      
      for(auto trackBin: binnedSelectedTracks){
	if(trackBin.size()<2) continue;
	for(uint i=0; i<trackBin.size()-1; i++){
	  for(uint j=i+1; j<trackBin.size(); j++){
	    bool isReal = false;
	    for(uint k=0; k<trueVertices.size(); k++){
	      int numMatched = 0;
	      for(uint l=0; l<trueVertices[k].tracks.size(); l++){
		if(trackBin[i].tp==trueVertices[k].tracks[l] || trackBin[j].tp==trueVertices[k].tracks[l]) numMatched++;
	      }
	      if(numMatched>=2) isReal = true;
	    }
	    if(!isReal){
	      h_findFake_trackVertexBinned_pt->Fill(trackBin[i].pt);
	      if(dist_TPs(trackBin[i],trackBin[j])==0){
		h_findFakeIntersect_trackVertexBinned_pt->Fill(trackBin[i].pt);
		Double_t x_dv_trk = -9999.0;
		Double_t y_dv_trk = -9999.0;
		Double_t z_dv_trk = -9999.0;
		calcVertex(trackBin[i],trackBin[j],x_dv_trk,y_dv_trk,z_dv_trk);
		Vertex_Parameters fakeVertex(x_dv_trk,y_dv_trk,z_dv_trk,trackBin[i],trackBin[j]);
		if(dist(x_dv_trk,y_dv_trk)>d0_res){
		  h_findFakeNonPrompt_trackVertexBinned_pt->Fill(trackBin[i].pt);
		  if(dist(x_dv_trk,y_dv_trk)<20){
		    h_findFakeBeforeTracker_trackVertexBinned_pt->Fill(trackBin[i].pt);
		    if(fakeVertex.bendchi2Sum<bendChi2Max){
		      h_findFakeBendChi2Cut_trackVertexBinned_pt->Fill(trackBin[i].pt);
		      if(fakeVertex.chi2rphidofSum<chi2RPhiMax){
			h_findFakeChi2RPhiCut_trackVertexBinned_pt->Fill(trackBin[i].pt);
			if(fakeVertex.cos_T>cosTMin){
			  h_findFakeCosTCut_trackVertexBinned_pt->Fill(trackBin[i].pt);
			  if(fakeVertex.delta_z<deltaZMax){
			    h_findFakeDeltaZCut_trackVertexBinned_pt->Fill(trackBin[i].pt);
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	} 
      }

    }//detailed plots
    int noCuts_trackVertices = 0;
    for(auto trackBin : binnedSelectedTracks){
      Double_t x_dv_trk = -9999.0;
      Double_t y_dv_trk = -9999.0;
      Double_t z_dv_trk = -9999.0;
      
      if(trackBin.size()>100){
	trackBin.erase(trackBin.begin()+100,trackBin.end());
      }
    
      while (trackBin.size()>1){
	std::vector<Vertex_Parameters> vertexCandidates;
	for( uint j=1; j<trackBin.size(); j++){
	  //std::cout<<"trying tp_pt pair: "<<trackBin[0].tp_pt<<" "<<trackBin[j].tp_pt<<std::endl;
	  if( dist_TPs( trackBin[0], trackBin[j] ) == 0 ){  
	    int inTraj = calcVertex(trackBin[0],trackBin[j],x_dv_trk,y_dv_trk,z_dv_trk);
	    //std::cout<<"no cut track vertex: "<<x_dv_trk<<" "<<y_dv_trk<<" "<<z_dv_trk<<" indices: "<<trackBin[0].index<<" "<<trackBin[j].index<<std::endl;
	    h_trackVertex_x->Fill(x_dv_trk);
	    h_trackVertex_y->Fill(y_dv_trk);
	    h_trackVertex_z->Fill(z_dv_trk);
	    h_trackVertex_sumPt->Fill(trackBin[0].pt+trackBin[j].pt);
	    noCuts_trackVertices++;
	    
	    float delta_z = fabs(trackBin[0].z(x_dv_trk,y_dv_trk)-trackBin[j].z(x_dv_trk,y_dv_trk));
	    float sumPt = trackBin[0].pt+trackBin[j].pt;
	    float minD0;
	    if(trackBin[0].d0<trackBin[j].d0){
	      minD0 = trackBin[0].d0;
	    }
	    else{
	      minD0 = trackBin[j].d0;
	    }
	    
	    float score = -1;
#if 0
	    for (int i_vert = 0; i_vert < (int)vert_score->size(); i_vert++){
	      if(vert_chi2rzdofSum->at(i_vert)==chi2rzdofSum && vert_chi2rphidofSum->at(i_vert)==chi2rphidofSum && vert_sumPt->at(i_vert)==sumPt){
		score = vert_score->at(i_vert);
	      }
	    }
#endif
	  //if(fabs(trackBin[0].z(x_dv_trk,y_dv_trk)-trackBin[j].z(x_dv_trk,y_dv_trk))< 1.25 && dist(x_dv_trk,y_dv_trk)>d0_res && dist(x_dv_trk,y_dv_trk)<20 && trackVertices.size()<2 && cos_T > 0.995 && d_T < .05 && (trackBin[0].charge+trackBin[j].charge)==0 && chi2rzdofSum<5 && numStubsSum>8 && chi2rphidofSum<20 && (fabs(trackBin[0].d0)>0.05 || fabs(trackBin[j].d0)>0.05) && (trackBin[0].pt+trackBin[j].pt)>20 && TMath::Sqrt(pow(trkMET[0],2)+pow(trkMET[1],2))>20 && trkH_T>160){ //vertex cuts go here (chi2rphi, MET, sumPT, and H_T cuts are not from 0 PU case)
	  //if(fabs(trackBin[0].z(x_dv_trk,y_dv_trk)-trackBin[j].z(x_dv_trk,y_dv_trk))< 1.05 && dist(x_dv_trk,y_dv_trk)>d0_res && dist(x_dv_trk,y_dv_trk)<20 && trackVertices.size()<2 && cos_T > 0.97 && d_T < .02 && R_T > 0.5 && (trackBin[0].charge+trackBin[j].charge)==0 && chi2rzdofSum<4 && numStubsSum>9 && (fabs(trackBin[0].d0)>0.2 || fabs(trackBin[j].d0)>0.2) ){ //vertex cuts from grid search
	    Vertex_Parameters trackVert = Vertex_Parameters(x_dv_trk, y_dv_trk, z_dv_trk, trackBin[0], trackBin[j], score, inTraj);
	    //if( dist(x_dv_trk,y_dv_trk)>d0_res && dist(x_dv_trk,y_dv_trk)<20 && trackVert.bendchi2Sum<bendChi2Max && trackVert.chi2rphidofSum<chi2RPhiMax && trackVert.cos_T>cosTMin && trackVert.delta_z<deltaZMax && trackVert.R_T>2 && trackBin[0].charge+trackBin[j].charge==0 ){ //vertex cuts go here
	    //if( dist(x_dv_trk,y_dv_trk)>d0_res && dist(x_dv_trk,y_dv_trk)<20 && trackBin[0].charge+trackBin[j].charge==0 && trackVert.chi2rphidofSum<chi2RPhiMax && trackVert.chi2rzdofSum<chi2RZMax && trackVert.cos_T>cosTMin && trackVert.delta_z<deltaZMax && trackVert.delta_eta<deltaEtaMax && trackVert.delta_eta>deltaEtaMin && trackVert.d_T<dTMax){ //vertex cuts go here
	    //if( dist(x_dv_trk,y_dv_trk)>d0_res && dist(x_dv_trk,y_dv_trk)<20 && trackBin[0].charge+trackBin[j].charge==0 ){ //vertex cuts go here
	    if( dist(x_dv_trk,y_dv_trk)>d0_res && dist(x_dv_trk,y_dv_trk)<20 && trackBin[0].charge+trackBin[j].charge==0 && trackVert.bendchi2Sum<bendChi2Max && trackVert.cos_T>cosTMin && trackVert.delta_eta<deltaEtaMax && trackVert.delta_z<deltaZMax && trackVert.R_T>0.25 && trackVert.a.pt>13.0 && trackVert.chi2rphidofSum<chi2RPhiMax && trackVert.numStubsSum>10 && (fabs(trackBin[0].d0)>0.3 || fabs(trackBin[j].d0)>0.3) ){ //vertex cuts go here
	      for(uint k=1; k<selectedTracks.size(); k++){
		if(selectedTracks[k]==trackBin[0] || selectedTracks[k]==trackBin[j]) continue;
		float delxy = dist_Vertex(x_dv_trk,y_dv_trk,selectedTracks[k]);
		float delz = fabs(selectedTracks[k].z(x_dv_trk,y_dv_trk)-z_dv_trk);
		if(delxy<delxyMax && delz<delzMax){
		  trackVert.tracks.push_back(selectedTracks[k]);
		}
	      }

	      h_trackVertex_inTraj->Fill(inTraj);
	      trackVertices.push_back(trackVert);
	      vertexCandidates.push_back(trackVert);
	      //std::cout<<"track vertex candidate: "<<x_dv_trk<<" "<<y_dv_trk<<" "<<z_dv_trk<<" tp_pt: "<<trackBin[0].tp->pt<<" "<<trackBin[j].tp->pt<<std::endl;
	    }
	  }
	}
	trackBin.pop_front();
	h_trackVertex_numVertPerTrack->Fill(vertexCandidates.size());
#if 0
	uint numCandidates = vertexCandidates.size();
	for(uint i=0; i<(trackVertices.size()-numCandidates); i++){
	  for(uint j=0; j<vertexCandidates.size();){
	    if(trackVertices[i]==vertexCandidates[j]){ 
	      vertexCandidates.erase(vertexCandidates.begin()+j);
	    }
	    else{
	      j++;
	    }
	  }
	}
#endif
	for(uint i=0; i<vertexCandidates.size(); i++){
	  for(uint j=0; j<trueVertices.size(); j++){
	    int numMatched = 0;
	    for(uint k=0; k<trueVertices[j].tracks.size(); k++){
	      if(vertexCandidates[i].a.tp==trueVertices[j].tracks[k] || vertexCandidates[i].b.tp==trueVertices[j].tracks[k]) numMatched++;
	    }
	    if(numMatched>=2) h_trackVertex_rankPt_vs_numVertPerTrack->Fill(vertexCandidates.size(),i);
	  }
	}
	sort(vertexCandidates.begin(),vertexCandidates.end(),CompareDelzVert);
	for(uint i=0; i<vertexCandidates.size(); i++){
	  for(uint j=0; j<trueVertices.size(); j++){
	    int numMatched = 0;
	    for(uint k=0; k<trueVertices[j].tracks.size(); k++){
	      if(vertexCandidates[i].a.tp==trueVertices[j].tracks[k] || vertexCandidates[i].b.tp==trueVertices[j].tracks[k]) numMatched++;
	    }
	    if(numMatched>=2) h_trackVertex_rankDelz_vs_numVertPerTrack->Fill(vertexCandidates.size(),i);
	  }
	}

	sort(vertexCandidates.begin(),vertexCandidates.end(),CompareDtVert);
	for(uint i=0; i<vertexCandidates.size(); i++){
	  for(uint j=0; j<trueVertices.size(); j++){
	    int numMatched = 0;
	    for(uint k=0; k<trueVertices[j].tracks.size(); k++){
	      if(vertexCandidates[i].a.tp==trueVertices[j].tracks[k] || vertexCandidates[i].b.tp==trueVertices[j].tracks[k]) numMatched++;
	    }
	    if(numMatched>=2) h_trackVertex_rankDt_vs_numVertPerTrack->Fill(vertexCandidates.size(),i);
	  }
	}

	sort(vertexCandidates.begin(),vertexCandidates.end(),CompareChi2rphidofSumVert);
	for(uint i=0; i<vertexCandidates.size(); i++){
	  for(uint j=0; j<trueVertices.size(); j++){
	    int numMatched = 0;
	    for(uint k=0; k<trueVertices[j].tracks.size(); k++){
	      if(vertexCandidates[i].a.tp==trueVertices[j].tracks[k] || vertexCandidates[i].b.tp==trueVertices[j].tracks[k]) numMatched++;
	    }
	    if(numMatched>=2) h_trackVertex_rankChi2rphidofSum_vs_numVertPerTrack->Fill(vertexCandidates.size(),i);
	  }
	}
	sort(vertexCandidates.begin(),vertexCandidates.end(),CompareRtVert);
	for(uint i=0; i<vertexCandidates.size(); i++){
	  for(uint j=0; j<trueVertices.size(); j++){
	    int numMatched = 0;
	    for(uint k=0; k<trueVertices[j].tracks.size(); k++){
	      if(vertexCandidates[i].a.tp==trueVertices[j].tracks[k] || vertexCandidates[i].b.tp==trueVertices[j].tracks[k]) numMatched++;
	    }
	    if(numMatched>=2) h_trackVertex_rankRt_vs_numVertPerTrack->Fill(vertexCandidates.size(),i);
	  }
	}
#if 0
	if(vertexCandidates.size()>0){
	  trackVertices.push_back(vertexCandidates[0]);
	  for(uint i=0; i<trackBin.size();){
	    if(trackBin[i]==vertexCandidates[0].b){
	      trackBin.erase(trackBin.begin()+i);
	    }
	    else{
	      i++;
	    }
	  }
	}
#endif	
      }
    }
    if(binVariable!=""){
      if(trackVertices.size()>1){
	for(uint i=0; i<trackVertices.size()-1; i++){
	  for(uint j=i+1; j<trackVertices.size();){
	    if(trackVertices[i]==trackVertices[j]){
	      trackVertices.erase(trackVertices.begin()+j);
	    }
	    else{
	      j++;
	    }
	  }
	}
      }
    }
    h_trackVertex_numNoCuts->Fill(noCuts_trackVertices);
    h_trackVertex_numAllCuts->Fill(trackVertices.size());
    std::vector<bool> check_dz(dz_cuts.size(),true);
    std::vector<bool> check_cos_T(cos_T_cuts.size(),true);
    std::vector<bool> check_d_T(d_T_cuts.size(),true);
    std::vector<bool> check_R_T(R_T_cuts.size(),true);
    std::vector<bool> check_chi2rz(chi2rz_cuts.size(),true);
    std::vector<bool> check_minD0(minD0_cuts.size(),true);
    std::vector<bool> check_stubSum(stubSum_cuts.size(),true);
    bool check_array[dz_cuts.size()][cos_T_cuts.size()][d_T_cuts.size()][R_T_cuts.size()][chi2rz_cuts.size()][minD0_cuts.size()][stubSum_cuts.size()];
    if(detailedPlots){
      for(uint k_a=0;k_a<dz_cuts.size();k_a++){
	for(uint k_b=0;k_b<cos_T_cuts.size();k_b++){
	  for(uint k_c=0;k_c<d_T_cuts.size();k_c++){
	    for(uint k_d=0;k_d<R_T_cuts.size();k_d++){
	      for(uint k_e=0;k_e<chi2rz_cuts.size();k_e++){
		for(uint k_f=0;k_f<minD0_cuts.size();k_f++){
		  for(uint k_g=0;k_g<stubSum_cuts.size();k_g++){
		    check_array[k_a][k_b][k_c][k_d][k_e][k_f][k_g] = true;
		  }
		}
	      }
	    }
	  }
	}
      }
    }

    int geomCounter = 0;
    for( uint i=0; i<trackVertices.size(); i++ ){
      //std::cout<<"trackVertex pos: "<<trackVertices[i].x_dv<<" "<<trackVertices[i].y_dv<<" "<<trackVertices[i].z_dv<<std::endl;
      int i_pt = 9999;
      int j_pt = 9999;
      if(detailedPlots){
	for( uint j=0; j<copyTracks.size();j++ ){
	  if(trackVertices[i].a.index == copyTracks[j].index){
	    i_pt = j;
	  }
	  if(trackVertices[i].b.index == copyTracks[j].index){
	    j_pt = j;
	  }
	}
      
      h_trackVertexCuts_indexPt->Fill(i_pt);
      h_trackVertexCuts_indexPt->Fill(j_pt);
      }
      h_trk_delta_dist_xy->Fill(dist_TPs(trackVertices[i].a,trackVertices[i].b));
      float z_dv_trk_1 = trackVertices[i].a.z(trackVertices[i].x_dv,trackVertices[i].y_dv);
      float z_dv_trk_2 = trackVertices[i].b.z(trackVertices[i].x_dv,trackVertices[i].y_dv);
      h_trk_delta_dist_z->Fill(fabs(z_dv_trk_1-z_dv_trk_2)); //* Use this condition
      h_trackVertex_d_T->Fill(trackVertices[i].d_T);
      h_trackVertex_cos_T->Fill(trackVertices[i].cos_T);
      h_trackVertex_alpha_T->Fill(trackVertices[i].alpha_T);
      h_trackVertex_R_T->Fill(trackVertices[i].R_T);
      h_trackVertex_openingAngle->Fill(trackVertices[i].openingAngle);
      h_trackVertex_parentPt->Fill(trackVertices[i].p_mag);

      if(trackVertices[i].x_dv!=-9999.0){
	h_trackVertexCuts_x->Fill(trackVertices[i].x_dv);
	h_trackVertexCuts_y->Fill(trackVertices[i].y_dv);
	h_trackVertexCuts_z->Fill(trackVertices[i].z_dv);
	h_trackVertexCuts_sumPt->Fill(trackVertices[i].a.pt+trackVertices[i].b.pt);


	h_all_trackVertex_pt->Fill(trackVertices[i].a.pt);
	if(fabs(trackVertices[i].a.d0) < fabs(trackVertices[i].b.d0)){
	  h_all_trackVertex_minD0->Fill(trackVertices[i].a.d0);
	}
	else{
	  h_all_trackVertex_minD0->Fill(trackVertices[i].b.d0);
	}
	h_all_trackVertex_lowPt->Fill(trackVertices[i].b.pt);
	h_all_trackVertex_eta->Fill(trackVertices[i].a.eta);
	h_all_trackVertex_dxy->Fill(dist(trackVertices[i].x_dv,trackVertices[i].y_dv));
	h_all_trackVertex_numStubs->Fill(trk_nstub->at(trackVertices[i].a.index));
	h_all_trackVertex_numStubs->Fill(trk_nstub->at(trackVertices[i].b.index));
	h_all_trackVertex_fakeId->Fill(trk_fake->at(trackVertices[i].a.index));
	h_all_trackVertex_fakeId->Fill(trk_fake->at(trackVertices[i].b.index));
	h_trackVertex_numTracks->Fill(trackVertices[i].tracks.size());

	if(detailedPlots){
	  for(uint k=0;k<dxy_cuts.size();k++){
	    if(dist_TPs(trackVertices[i].a,trackVertices[i].b) < dxy_cuts[k]){
	      all_vert_dxy_cut[k]++;
	    }
	  }
	  for(uint k=0;k<dz_cuts.size();k++){
	    if(fabs(z_dv_trk_1-z_dv_trk_2) < dz_cuts[k] && check_dz[k]==true){
	      all_vert_dz_cut[k]++;
	    }
	  }
	  for(uint k=0;k<cos_T_cuts.size();k++){
	    if(trackVertices[i].cos_T > cos_T_cuts[k] && check_cos_T[k]==true){
	      all_vert_cos_T_cut[k]++;
	    }
	  }
	  for(uint k=0;k<d_T_cuts.size();k++){
	    if(trackVertices[i].d_T < d_T_cuts[k] && check_d_T[k]==true){
	      all_vert_d_T_cut[k]++;
	    }
	  }
	  for(uint k=0;k<R_T_cuts.size();k++){
	    if(trackVertices[i].R_T > R_T_cuts[k] && check_R_T[k]==true){
	      all_vert_R_T_cut[k]++;
	    }
	  }
	  
	  for(uint k=0;k<chi2rz_cuts.size();k++){
	    if(trackVertices[i].chi2rzdofSum < chi2rz_cuts[k] && check_chi2rz[k]==true){
	      all_vert_chi2rz_cut[k]++;
	    }
	  }
	  for(uint k=0;k<minD0_cuts.size();k++){
	    if( ((fabs(trackVertices[i].a.d0) > minD0_cuts[k]) || (fabs(trackVertices[i].b.d0) > minD0_cuts[k])) && check_minD0[k]==true){
	      all_vert_minD0_cut[k]++;
	    }
	  }
	  for(uint k=0;k<stubSum_cuts.size();k++){
	    if(trackVertices[i].numStubsSum > stubSum_cuts[k] && check_stubSum[k]==true){
	      all_vert_stubSum_cut[k]++;
	    }
	  }
	  for(uint k_a=0;k_a<dz_cuts.size();k_a++){
	    for(uint k_b=0;k_b<cos_T_cuts.size();k_b++){
	      for(uint k_c=0;k_c<d_T_cuts.size();k_c++){
		for(uint k_d=0;k_d<R_T_cuts.size();k_d++){
		  for(uint k_e=0;k_e<chi2rz_cuts.size();k_e++){
		    for(uint k_f=0;k_f<minD0_cuts.size();k_f++){
		      for(uint k_g=0;k_g<stubSum_cuts.size();k_g++){
			if( (fabs(z_dv_trk_1-z_dv_trk_2) < dz_cuts[k_a]) && (trackVertices[i].cos_T > cos_T_cuts[k_b]) && (trackVertices[i].d_T < d_T_cuts[k_c]) && (trackVertices[i].R_T > R_T_cuts[k_d]) && (trackVertices[i].chi2rzdofSum < chi2rz_cuts[k_e]) && ((fabs(trackVertices[i].a.d0) > minD0_cuts[k_f]) || (fabs(trackVertices[i].b.d0) > minD0_cuts[k_f])) && (trackVertices[i].numStubsSum > stubSum_cuts[k_g]) && check_array[k_a][k_b][k_c][k_d][k_e][k_f][k_g]==true ){
			  all_vert_array_cut[k_a][k_b][k_c][k_d][k_e][k_f][k_g]++;
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }

	}
      }
      uint matched_j = 0;
      bool foundMatch = false;
      int aMatched = 0;
      int bMatched = 0;
      for(uint j=0; j<trueVertices.size(); j++){
	int numMatched = 0;
	for(uint k=0; k<trueVertices[j].tracks.size(); k++){
	  if( (trackVertices[i].a.tp == trueVertices[j].tracks[k]) && trueVertices[j].matched==false ){
	    numMatched++;
	    aMatched = 1;
	  }
	  if( (trackVertices[i].b.tp == trueVertices[j].tracks[k]) && trueVertices[j].matched==false ){ 
	    numMatched++;
	    bMatched = 1;
	  }
	}
	if(numMatched>=2){
	  //trueVertices[j].matched = true;
	  trueVertices[j].matched = false;
	  matched_j = j;
	  foundMatch = true;
	  if(oneMatch==false){
	    oneMatch = true;
	    firstMatch_j = j;
	  }
	  break;	    
	}
      }
    
      if(true_DV && foundMatch){
	geomCounter++;
	if(geomCounter==1&&i_evnt==276){
	  Double_t dummyVert_x;
	  Double_t dummyVert_y;
	  Double_t dummyVert_z;
	  calcVertex(trackVertices[i].a,trackVertices[i].b,dummyVert_x,dummyVert_y,dummyVert_z);
	  geomTrackVertex = trackVertices[i];
	  geomTrueVertex = trueVertices[matched_j];
	}

	for(uint k=0; k<trueVertices[matched_j].tracks.size(); k++){
	  //std::cout<<"trueVertex tp_pt: "<<trueVertices[matched_j].tracks[k].pt<<std::endl;
	  if(trackVertices[i].a.tp == trueVertices[matched_j].tracks[k] || trackVertices[i].b.tp == trueVertices[matched_j].tracks[k]) continue; 
	  h_tpAssoc_pt->Fill(trueVertices[matched_j].tracks[k].pt);   
	  bool isAssoc = false;
	  for(uint itrack=2; itrack<trackVertices[i].tracks.size(); itrack++){
	    if(trackVertices[i].tracks[itrack].tp==trueVertices[matched_j].tracks[k]) isAssoc=true;
	  }
	  if(isAssoc) h_tpAssoc_pt_matched->Fill(trueVertices[matched_j].tracks[k].pt);
	}

	for(uint itrack=2; itrack<trackVertices[i].tracks.size(); itrack++){
	  h_trkAssoc_pt->Fill(trackVertices[i].tracks[itrack].pt);
	  bool isAssoc = false;
	  for(uint k=0; k<trueVertices[matched_j].tracks.size(); k++){
	    if(trackVertices[i].tracks[itrack].tp==trueVertices[matched_j].tracks[k]) isAssoc = true;
	  }
	  float delxy = dist_Vertex(trackVertices[i].x_dv,trackVertices[i].y_dv,trackVertices[i].tracks[itrack]);
	  float delz = fabs(trackVertices[i].tracks[itrack].z(trackVertices[i].x_dv,trackVertices[i].y_dv)-trackVertices[i].z_dv);
	  std::valarray<float> p_proj = calcPVec(trackVertices[i].tracks[itrack],trackVertices[i].x_dv,trackVertices[i].y_dv);
	  float delPhiProj = fabs(trackVertices[i].phi-atan2(p_proj[1],p_proj[0]));
	  float delPt = fabs(trackVertices[i].tracks[itrack].pt - trackVertices[i].a.pt);
	  if(!isAssoc){
	    h_trkAssoc_pt_noMatch->Fill(trackVertices[i].tracks[itrack].pt);
	    h_false_trkAssoc_delPhi->Fill(fabs(trackVertices[i].phi-trackVertices[i].tracks[itrack].phi));
	    h_false_trkAssoc_delxy->Fill(delxy);
	    h_false_trkAssoc_delz->Fill(delz);
	    h_false_trkAssoc_delPhiProp->Fill(delPhiProj);
	    h_false_trkAssoc_pt->Fill(trackVertices[i].tracks[itrack].pt);
	    h_false_trkAssoc_delPt->Fill(delPt);
	  }
	  else{
	    h_correct_trkAssoc_delPhi->Fill(fabs(trackVertices[i].phi-trackVertices[i].tracks[itrack].phi));
	    h_correct_trkAssoc_delxy->Fill(delxy);
	    h_correct_trkAssoc_delz->Fill(delz);
	    h_correct_trkAssoc_delPhiProp->Fill(delPhiProj);
	    h_correct_trkAssoc_pt->Fill(trackVertices[i].tracks[itrack].pt);
	    h_correct_trkAssoc_delPt->Fill(delPt);
	  }
	}

	/*for(uint itrack=0; itrack<trackVertices[i].tracks.size(); itrack++){
	  std::cout<<"trackVertex track tp_pt: "<<trackVertices[i].tracks[itrack].tp->pt<<std::endl;
	  }*/
      
	h_correct_trueVertex_pt->Fill(trueVertices[matched_j].a.pt);
	h_correct_trueVertex_lowPt->Fill(trueVertices[matched_j].b.pt);
	h_correct_trueVertex_eta->Fill(trueVertices[matched_j].a.eta);
	h_correct_trueVertex_dxy->Fill(dist(trueVertices[matched_j].x_dv,trueVertices[matched_j].y_dv));
	h_correct_trackVertex_numStubs->Fill(trk_nstub->at(trackVertices[i].a.index));
	h_correct_trackVertex_numStubs->Fill(trk_nstub->at(trackVertices[i].b.index));
	h_correct_trackVertex_numStubsSum->Fill(trk_nstub->at(trackVertices[i].a.index)+trk_nstub->at(trackVertices[i].b.index));
	h_correct_trackVertex_chi2rphidofSum->Fill(trackVertices[i].chi2rphidofSum);
	h_correct_trackVertex_MVA1Sum->Fill(trackVertices[i].MVA1Sum);
	h_correct_trackVertex_MVA2Sum->Fill(trackVertices[i].MVA2Sum);
	h_correct_trackVertex_chi2rzdofSum->Fill(trackVertices[i].chi2rzdofSum);
	h_correct_trackVertex_bendchi2Sum->Fill(trackVertices[i].bendchi2Sum);
	h_correct_trackVertex_score->Fill(trackVertices[i].score);
	h_correct_trackVertex_delta_dist_z->Fill(fabs(z_dv_trk_1-z_dv_trk_2));
	if(trackVertices[i].inTraj==0){
	  h_correct_trackVertex_delta_dist_z_inBothTraj->Fill(fabs(z_dv_trk_1-z_dv_trk_2));
	}
	else if(trackVertices[i].inTraj==1 || trackVertices[i].inTraj==2){
	  h_correct_trackVertex_delta_dist_z_inOneTraj->Fill(fabs(z_dv_trk_1-z_dv_trk_2));
	}
	else if(trackVertices[i].inTraj==3){
	  h_correct_trackVertex_delta_dist_z_inNoTraj->Fill(fabs(z_dv_trk_1-z_dv_trk_2));
	}
	h_correct_trackVertex_delta_dist_z0->Fill(fabs(trk_z0->at(trackVertices[i].a.index)-trk_z0->at(trackVertices[i].b.index)));
	h_correct_trackVertex_delta_dist_d0->Fill(fabs(trk_d0->at(trackVertices[i].a.index)-trk_d0->at(trackVertices[i].b.index)));
	h_correct_trackVertex_delta_dist_eta->Fill(fabs(trk_eta->at(trackVertices[i].a.index)-trk_eta->at(trackVertices[i].b.index)));
	h_correct_trackVertex_delta_dist_phi->Fill(fabs(trk_phi->at(trackVertices[i].a.index)-trk_phi->at(trackVertices[i].b.index)));
	if(detailedPlots){
	  h_correct_trackVertex_delta_dist_indexPt->Fill(abs(i_pt-j_pt));
	  h_correct_trackVertex_indexPt->Fill(i_pt);
	  h_correct_trackVertex_indexPt->Fill(j_pt);
	}
	h_correct_trackVertex_R_T->Fill(trackVertices[i].R_T);
	h_correct_trackVertex_cos_T->Fill(trackVertices[i].cos_T);
	h_correct_trackVertex_p2_mag->Fill(trackVertices[i].p2_mag);
	h_correct_trackVertex_lowPt->Fill(trackVertices[i].b.pt);
	h_correct_trackVertex_pt->Fill(trackVertices[i].a.pt);
      
	for(uint j=i+1; j<trackVertices.size(); j++){
	  float deltaPos = TMath::Sqrt(pow(trackVertices[i].x_dv-trackVertices[j].x_dv,2)+pow(trackVertices[i].y_dv-trackVertices[j].y_dv,2)+pow(trackVertices[i].z_dv-trackVertices[j].z_dv,2));
	  h_correct_trackVertex_deltaPos->Fill(deltaPos);
	}
      
	h_correct_trackVertex_d_T->Fill(trackVertices[i].d_T);
	h_correct_trackVertex_numTracks->Fill(trackVertices[i].tracks.size());
	h_correct_trackVertex_inTraj->Fill(trackVertices[i].inTraj);
	if(fabs(trackVertices[i].a.d0) < fabs(trackVertices[i].b.d0)){
	  h_correct_trackVertex_minD0->Fill(fabs(trackVertices[i].a.d0));
	  h_correct_trackVertex_maxD0->Fill(fabs(trackVertices[i].b.d0));
	}
	else{
	  h_correct_trackVertex_minD0->Fill(fabs(trackVertices[i].b.d0));
	  h_correct_trackVertex_maxD0->Fill(fabs(trackVertices[i].a.d0));
	}
	int netCharge = 0;
	for(uint itrack=0; itrack<trackVertices[i].tracks.size(); itrack++){
	  netCharge+=trackVertices[i].tracks[itrack].charge;
	}
	h_correct_trackVertex_charge_vs_numTracks->Fill(trackVertices[i].tracks.size(),netCharge);
	if(detailedPlots){ 
	  for(uint k=0;k<dxy_cuts.size();k++){
	    if(dist_TPs(trackVertices[i].a,trackVertices[i].b) < dxy_cuts[k]){
	      correct_vert_dxy_cut[k]++;
	    }
	  }
	  for(uint k=0;k<dz_cuts.size();k++){
	    if(fabs(z_dv_trk_1-z_dv_trk_2) < dz_cuts[k] && check_dz[k]==true){
	      correct_vert_dz_cut[k]++;
	    }
	  }
	  for(uint k=0;k<cos_T_cuts.size();k++){
	    if(trackVertices[i].cos_T > cos_T_cuts[k] && check_cos_T[k]==true){
	      correct_vert_cos_T_cut[k]++;
	    }
	  }
	  for(uint k=0;k<d_T_cuts.size();k++){
	    if(trackVertices[i].d_T < d_T_cuts[k] && check_d_T[k]==true){
	      correct_vert_d_T_cut[k]++;
	    }
	  }
	  for(uint k=0;k<R_T_cuts.size();k++){
	    if(trackVertices[i].R_T > R_T_cuts[k] && check_R_T[k]==true){
	      correct_vert_R_T_cut[k]++;
	    }
	  }
	  for(uint k=0;k<chi2rz_cuts.size();k++){
	    if(trackVertices[i].chi2rzdofSum < chi2rz_cuts[k] && check_chi2rz[k]==true){
	      correct_vert_chi2rz_cut[k]++;
	    }
	  }
	  for(uint k=0;k<minD0_cuts.size();k++){
	    if( ((fabs(trackVertices[i].a.d0) > minD0_cuts[k]) || (fabs(trackVertices[i].b.d0) > minD0_cuts[k])) && check_minD0[k]==true){
	      correct_vert_minD0_cut[k]++;
	    }
	  }
	  for(uint k=0;k<stubSum_cuts.size();k++){
	    if(trackVertices[i].numStubsSum > stubSum_cuts[k] && check_stubSum[k]==true){
	      correct_vert_stubSum_cut[k]++;
	    }
	  }
	  
	  for(uint k_a=0;k_a<dz_cuts.size();k_a++){
	    for(uint k_b=0;k_b<cos_T_cuts.size();k_b++){
	      for(uint k_c=0;k_c<d_T_cuts.size();k_c++){
		for(uint k_d=0;k_d<R_T_cuts.size();k_d++){
		  for(uint k_e=0;k_e<chi2rz_cuts.size();k_e++){
		    for(uint k_f=0;k_f<minD0_cuts.size();k_f++){
		      for(uint k_g=0;k_g<stubSum_cuts.size();k_g++){
			if( (fabs(z_dv_trk_1-z_dv_trk_2) < dz_cuts[k_a]) && (trackVertices[i].cos_T > cos_T_cuts[k_b]) && (trackVertices[i].d_T < d_T_cuts[k_c]) && (trackVertices[i].R_T > R_T_cuts[k_d]) && (trackVertices[i].chi2rzdofSum < chi2rz_cuts[k_e]) && ((fabs(trackVertices[i].a.d0) > minD0_cuts[k_f]) || (fabs(trackVertices[i].b.d0) > minD0_cuts[k_f])) && (trackVertices[i].numStubsSum > stubSum_cuts[k_g]) && check_array[k_a][k_b][k_c][k_d][k_e][k_f][k_g]==true ){
			  correct_vert_array_cut[k_a][k_b][k_c][k_d][k_e][k_f][k_g]++;
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }

	}
	if(fabs(trueVertices[matched_j].x_dv-trackVertices[i].x_dv)>5){
	  //std::cout<<"matched trueVertex pos: "<<trueVertices[matched_j].x_dv<<" "<<trueVertices[matched_j].y_dv<<" "<<trueVertices[matched_j].z_dv<<" trackVertex pos: "<<trackVertices[i].x_dv<<" "<<trackVertices[i].y_dv<<" "<<trackVertices[i].z_dv<<std::endl;
	  //std::cout<<"res_x: "<<trueVertices[matched_j].x_dv-trackVertices[i].x_dv<<std::endl;
	  //std::cout<<"track vertex charges: "<<trackVertices[i].a.charge<<" "<<trackVertices[i].b.charge<<" trueVertex charges: "<<trueVertices[matched_j].a.charge<<" "<<trueVertices[matched_j].b.charge<<std::endl;
	}

	h_res_tp_trk_x->Fill(trueVertices[matched_j].x_dv-trackVertices[i].x_dv);
	h_res_tp_trk_x_zoomOut->Fill(trueVertices[matched_j].x_dv-trackVertices[i].x_dv);
	h_res_tp_trk_y->Fill(trueVertices[matched_j].y_dv-trackVertices[i].y_dv);
	h_res_tp_trk_y_zoomOut->Fill(trueVertices[matched_j].y_dv-trackVertices[i].y_dv);
	h_res_tp_trk_z->Fill(trueVertices[matched_j].z_dv-trackVertices[i].z_dv);
	float r_j = dist(trueVertices[matched_j].x_dv,trueVertices[matched_j].y_dv);
	float r_i = dist(trackVertices[i].x_dv,trackVertices[i].y_dv);
	float phi_j = atan2(trueVertices[matched_j].y_dv,trueVertices[matched_j].x_dv);
	float phi_i = atan2(trackVertices[i].y_dv,trackVertices[i].x_dv);
	h_res_tp_trk_r->Fill(r_j-r_i);
	h_res_tp_trk_phi->Fill(phi_j-phi_i);
       
      }
      else{
	if(trackVertices[i].x_dv!=-9999.0){
	  h_false_trackVertex_pt->Fill(trackVertices[i].a.pt);
	  h_false_trackVertex_eta->Fill(trackVertices[i].a.eta);
	  h_false_trackVertex_dxy->Fill(dist(trackVertices[i].x_dv,trackVertices[i].y_dv));
	  h_false_trackVertex_delta_dist_z->Fill(fabs(z_dv_trk_1-z_dv_trk_2));
	  if(trackVertices[i].inTraj==0){
	    h_false_trackVertex_delta_dist_z_inBothTraj->Fill(fabs(z_dv_trk_1-z_dv_trk_2));
	  }
	  else if(trackVertices[i].inTraj==1 || trackVertices[i].inTraj==2){
	    h_false_trackVertex_delta_dist_z_inOneTraj->Fill(fabs(z_dv_trk_1-z_dv_trk_2));
	  }
	  else if(trackVertices[i].inTraj==3){
	    h_false_trackVertex_delta_dist_z_inNoTraj->Fill(fabs(z_dv_trk_1-z_dv_trk_2));
	  }
	  h_false_trackVertex_numStubs->Fill(trk_nstub->at(trackVertices[i].a.index));
	  h_false_trackVertex_numStubs->Fill(trk_nstub->at(trackVertices[i].b.index));
	  h_false_trackVertex_numStubsSum->Fill(trk_nstub->at(trackVertices[i].a.index)+trk_nstub->at(trackVertices[i].b.index));
	  h_false_trackVertex_chi2rphidofSum->Fill(trackVertices[i].chi2rphidofSum);
	  h_false_trackVertex_MVA1Sum->Fill(trackVertices[i].MVA1Sum);
	  h_false_trackVertex_MVA2Sum->Fill(trackVertices[i].MVA2Sum);
	  h_false_trackVertex_chi2rzdofSum->Fill(trackVertices[i].chi2rzdofSum);
	  h_false_trackVertex_bendchi2Sum->Fill(trackVertices[i].bendchi2Sum);
	  h_false_trackVertex_score->Fill(trackVertices[i].score);
	  h_false_trackVertex_d0->Fill(trackVertices[i].a.d0);
	  h_false_trackVertex_d_T->Fill(trackVertices[i].d_T);
	  h_false_trackVertex_R_T->Fill(trackVertices[i].R_T);
	  h_false_trackVertex_cos_T->Fill(trackVertices[i].cos_T);
	  h_false_trackVertex_p2_mag->Fill(trackVertices[i].p2_mag);
	  h_false_trackVertex_lowPt->Fill(trackVertices[i].b.pt);
	  h_false_trackVertex_maxPt->Fill(trackVertices[i].a.pt);

	  for(uint j=i+1; j<trackVertices.size(); j++){
	    float deltaPos = TMath::Sqrt(pow(trackVertices[i].x_dv-trackVertices[j].x_dv,2)+pow(trackVertices[i].y_dv-trackVertices[j].y_dv,2)+pow(trackVertices[i].z_dv-trackVertices[j].z_dv,2));
	    h_false_trackVertex_deltaPos->Fill(deltaPos);
	  }

	  h_false_trackVertex_numTracks->Fill(trackVertices[i].tracks.size());
	  h_false_trackVertex_inTraj->Fill(trackVertices[i].inTraj);
	  h_false_trackVertex_numMatched->Fill(aMatched+bMatched);
	  h_false_trackVertex_fakeId->Fill(trk_fake->at(trackVertices[i].a.index));
	  h_false_trackVertex_fakeId->Fill(trk_fake->at(trackVertices[i].b.index));
	  //std::cout<<"fakeVertex track pdgid's: "<<trackVertices[i].a.tp->pdgid<<" "<<trackVertices[i].b.tp->pdgid<<" isHToB: "<<trk_matchtp_isHToB->at(trackVertices[i].a.index)<<" "<<trk_matchtp_isHToB->at(trackVertices[i].b.index)<<" vertex pos x y z: ("<<trackVertices[i].a.tp->vx<<" "<<trackVertices[i].a.tp->vy<<" "<<trackVertices[i].a.tp->vz<<") ("<<trackVertices[i].b.tp->vx<<" "<<trackVertices[i].b.tp->vy<<" "<<trackVertices[i].b.tp->vz<<")"<<std::endl;
	  h_false_trackVertex_delta_dist_eta->Fill(fabs(trk_eta->at(trackVertices[i].a.index)-trk_eta->at(trackVertices[i].b.index)));
	  h_false_trackVertex_delta_dist_phi->Fill(fabs(trk_phi->at(trackVertices[i].a.index)-trk_phi->at(trackVertices[i].b.index)));
	  if(fabs(trackVertices[i].a.d0) < fabs(trackVertices[i].b.d0)){
	    h_false_trackVertex_minD0->Fill(fabs(trackVertices[i].a.d0));
	    h_false_trackVertex_maxD0->Fill(fabs(trackVertices[i].b.d0));
	  }
	  else{
	    h_false_trackVertex_minD0->Fill(fabs(trackVertices[i].b.d0));
	    h_false_trackVertex_maxD0->Fill(fabs(trackVertices[i].a.d0));
	  }
	  int netCharge = 0;
	  for(uint itrack=0; itrack<trackVertices[i].tracks.size(); itrack++){
	    netCharge+=trackVertices[i].tracks[itrack].charge;
	  }
	  h_false_trackVertex_charge_vs_numTracks->Fill(trackVertices[i].tracks.size(),netCharge);
	  if(detailedPlots){
	    for(uint k=0;k<dxy_cuts.size();k++){
	      if(dist_TPs(trackVertices[i].a,trackVertices[i].b) < dxy_cuts[k]){
		false_vert_dxy_cut[k]++;
	      }
	    }
	    for(uint k=0;k<dz_cuts.size();k++){
	      if(fabs(z_dv_trk_1-z_dv_trk_2) < dz_cuts[k] && check_dz[k]==true){
		false_vert_dz_cut[k]++;
	      }
	    }
	    for(uint k=0;k<cos_T_cuts.size();k++){
	      if(trackVertices[i].cos_T > cos_T_cuts[k] && check_cos_T[k]==true){
		false_vert_cos_T_cut[k]++;
	      }
	    }
	    for(uint k=0;k<d_T_cuts.size();k++){
	      if(trackVertices[i].d_T < d_T_cuts[k] && check_d_T[k]==true){
		false_vert_d_T_cut[k]++;
	      }
	    }
	    for(uint k=0;k<R_T_cuts.size();k++){
	      if(trackVertices[i].R_T > R_T_cuts[k] && check_R_T[k]==true){
		false_vert_R_T_cut[k]++;
	      }
	    }
	    for(uint k=0;k<chi2rz_cuts.size();k++){
	      if(trackVertices[i].chi2rzdofSum < chi2rz_cuts[k] && check_chi2rz[k]==true){
		false_vert_chi2rz_cut[k]++;
	      }
	  }
	    for(uint k=0;k<minD0_cuts.size();k++){
	      if( ((fabs(trackVertices[i].a.d0) > minD0_cuts[k]) || (fabs(trackVertices[i].b.d0) > minD0_cuts[k])) && check_minD0[k]==true){
		false_vert_minD0_cut[k]++;
	    }
	    }
	    for(uint k=0;k<stubSum_cuts.size();k++){
	      if(trackVertices[i].numStubsSum > stubSum_cuts[k] && check_stubSum[k]==true){
		false_vert_stubSum_cut[k]++;
	      }
	    }
	    
	    for(uint k_a=0;k_a<dz_cuts.size();k_a++){
	      for(uint k_b=0;k_b<cos_T_cuts.size();k_b++){
		for(uint k_c=0;k_c<d_T_cuts.size();k_c++){
		  for(uint k_d=0;k_d<R_T_cuts.size();k_d++){
		    for(uint k_e=0;k_e<chi2rz_cuts.size();k_e++){
		      for(uint k_f=0;k_f<minD0_cuts.size();k_f++){
			for(uint k_g=0;k_g<stubSum_cuts.size();k_g++){
			  if( (fabs(z_dv_trk_1-z_dv_trk_2) < dz_cuts[k_a]) && (trackVertices[i].cos_T > cos_T_cuts[k_b]) && (trackVertices[i].d_T < d_T_cuts[k_c]) && (trackVertices[i].R_T > R_T_cuts[k_d]) && (trackVertices[i].chi2rzdofSum < chi2rz_cuts[k_e]) && ((fabs(trackVertices[i].a.d0) > minD0_cuts[k_f]) || (fabs(trackVertices[i].b.d0) > minD0_cuts[k_f])) && (trackVertices[i].numStubsSum > stubSum_cuts[k_g]) && check_array[k_a][k_b][k_c][k_d][k_e][k_f][k_g]==true ){
			    false_vert_array_cut[k_a][k_b][k_c][k_d][k_e][k_f][k_g]++;
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }

	  }
	}
      }
      if(detailedPlots){
	for(uint k=0;k<dz_cuts.size();k++){
	  if(all_vert_dz_cut[k]==prev_dz_cut[k]+2){
	    //check_dz[k] = false;
	  }
	}
	for(uint k=0;k<cos_T_cuts.size();k++){
	  if(all_vert_cos_T_cut[k]==prev_cos_T_cut[k]+2){
	    //check_cos_T[k] = false;
	  }
	}
	for(uint k=0;k<d_T_cuts.size();k++){
	  if(all_vert_d_T_cut[k]==prev_d_T_cut[k]+2){
	    //check_d_T[k] = false;
	  }
	}
	for(uint k=0;k<R_T_cuts.size();k++){
	  if(all_vert_R_T_cut[k]==prev_R_T_cut[k]+2){
	    //check_R_T[k] = false;
	  }
	}
	for(uint k=0;k<chi2rz_cuts.size();k++){
	  if(all_vert_chi2rz_cut[k]==prev_chi2rz_cut[k]+2){
	    //check_chi2rz[k] = false;
	  }
	}
	for(uint k=0;k<minD0_cuts.size();k++){
	  if(all_vert_minD0_cut[k]==prev_minD0_cut[k]+2){
	    //check_minD0[k] = false;
	  }
	}
	for(uint k=0;k<stubSum_cuts.size();k++){
	  if(all_vert_stubSum_cut[k]==prev_stubSum_cut[k]+2){
	    //check_stubSum[k] = false;
	  }
	}
	for(uint k_a=0;k_a<dz_cuts.size();k_a++){
	  for(uint k_b=0;k_b<cos_T_cuts.size();k_b++){
	    for(uint k_c=0;k_c<d_T_cuts.size();k_c++){
	      for(uint k_d=0;k_d<R_T_cuts.size();k_d++){
		for(uint k_e=0;k_e<chi2rz_cuts.size();k_e++){
		  for(uint k_f=0;k_f<minD0_cuts.size();k_f++){
		    for(uint k_g=0;k_g<stubSum_cuts.size();k_g++){
		      if( all_vert_array_cut[k_a][k_b][k_c][k_d][k_e][k_f][k_g]==prev_array_cut[k_a][k_b][k_c][k_d][k_e][k_f][k_g]+2){
			check_array[k_a][k_b][k_c][k_d][k_e][k_f][k_g] = false;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}

      }
    } // End of Loop of TrackVertices
    if(detailedPlots){
      prev_dz_cut = all_vert_dz_cut;
      prev_cos_T_cut = all_vert_cos_T_cut;
      prev_d_T_cut = all_vert_d_T_cut;
      prev_R_T_cut = all_vert_R_T_cut;
      prev_chi2rz_cut = all_vert_chi2rz_cut;
      prev_minD0_cut = all_vert_minD0_cut;
      prev_stubSum_cut = all_vert_stubSum_cut;
      
      for(uint k_a=0;k_a<dz_cuts.size();k_a++){
	for(uint k_b=0;k_b<cos_T_cuts.size();k_b++){
	  for(uint k_c=0;k_c<d_T_cuts.size();k_c++){
	    for(uint k_d=0;k_d<R_T_cuts.size();k_d++){
	      for(uint k_e=0;k_e<chi2rz_cuts.size();k_e++){
		for(uint k_f=0;k_f<minD0_cuts.size();k_f++){
		  for(uint k_g=0;k_g<stubSum_cuts.size();k_g++){
		    prev_array_cut[k_a][k_b][k_c][k_d][k_e][k_f][k_g] = all_vert_array_cut[k_a][k_b][k_c][k_d][k_e][k_f][k_g];
		  }
		}
	      }
	    }
	  }
	}
      }

    }

    if(trackVertices.size()>0){
      h_trk_oneVert_MET->Fill(TMath::Sqrt(pow(trkMET[0],2)+pow(trkMET[1],2)));
      h_trk_oneVert_H_T->Fill(trkH_T);
    }
      
    if(oneMatch){
      h_correct_oneMatch_trueVertex_pt->Fill(trueVertices[maxPT_i].a.pt);
      h_correct_oneMatch_trueVertex_lowPt->Fill(trueVertices[maxPT_i].b.pt);
      h_correct_oneMatch_trueVertex_eta->Fill(trueVertices[maxPT_i].a.eta);
      h_correct_oneMatch_trueVertex_dxy->Fill(dist(trueVertices[maxPT_i].x_dv,trueVertices[maxPT_i].y_dv));
      h_correct_oneMatchAlt_trueVertex_pt->Fill(trueVertices[firstMatch_j].a.pt);
      h_correct_oneMatchAlt_trueVertex_eta->Fill(trueVertices[firstMatch_j].a.eta);
      h_correct_oneMatchAlt_trueVertex_dxy->Fill(dist(trueVertices[firstMatch_j].x_dv,trueVertices[firstMatch_j].y_dv));
      h_all_oneMatchAlt_trueVertex_pt->Fill(trueVertices[firstMatch_j].a.pt);
      h_all_oneMatchAlt_trueVertex_eta->Fill(trueVertices[firstMatch_j].a.eta);
      h_all_oneMatchAlt_trueVertex_dxy->Fill(dist(trueVertices[firstMatch_j].x_dv,trueVertices[firstMatch_j].y_dv));
      h_trk_oneMatch_MET->Fill(TMath::Sqrt(pow(trkMET[0],2)+pow(trkMET[1],2)));
      h_trk_oneMatch_H_T->Fill(trkH_T);
    }
    else{
      if(trueVertices.size()>0){
	h_all_oneMatchAlt_trueVertex_pt->Fill(trueVertices[maxPT_i].a.pt);
	h_all_oneMatchAlt_trueVertex_eta->Fill(trueVertices[maxPT_i].a.eta);
	h_all_oneMatchAlt_trueVertex_dxy->Fill(dist(trueVertices[maxPT_i].x_dv,trueVertices[maxPT_i].y_dv));
      }
    }
  } // End of Event Loop

  // ---------------------------------------------------------------------------------------------------------
  //some Histograms

  char ctxt[500];
  if(inputFile.Contains("DarkPhoton")){
    if(inputFile.Contains("cT0")){
      sprintf(ctxt, "Dark Photon, PU=0, #tau=0mm");
    }
    else if(inputFile.Contains("cT10000")){
      sprintf(ctxt, "Dark Photon, PU=0, #tau=10000mm");
    }
    else if(inputFile.Contains("cT5000")){
      sprintf(ctxt, "Dark Photon, PU=0, #tau=5000mm");
    }   
    else if(inputFile.Contains("cT100")){
      sprintf(ctxt, "Dark Photon, PU=200, #tau=100mm");
    }
    else if(inputFile.Contains("cT10")){
      if(inputFile.Contains("PU200")){
	sprintf(ctxt, "Dark Photon, PU=200, #tau=10mm");
      }
      else{
	sprintf(ctxt, "Dark Photon, PU=200, #tau=10mm");
      }
    }
  }
  else if(inputFile.Contains("DisplacedTrackJet")){
    if(inputFile.Contains("cT10")){
      if(inputFile.Contains("PU200")){
	sprintf(ctxt, "DispTrkJet, PU=200, #tau=10mm");
      }
      else{
	sprintf(ctxt, "DispTrkJet, PU=200, #tau=10mm");
      }
    }
  }
  else if(inputFile.Contains("NeutrinoGun")){
    sprintf(ctxt, "Neutrino Gun, PU=200");
  }
  else if(inputFile.Contains("DispMu")){
    if(inputFile.Contains("PU200")){
      sprintf(ctxt, "Displaced Mu, PU=200");
    }
    else{
      sprintf(ctxt, "Displaced Mu, PU=0");
    }
  }
  else if(inputFile.Contains("TTbar")){
    if(inputFile.Contains("PU200")){
      sprintf(ctxt, "TTbar, PU=200");
    }
    else{
      sprintf(ctxt, "TTbar, PU=0");
    }
  }
  else{
    sprintf(ctxt, " ");
  }
  TCanvas c;

  TString DIR = outputDir + "AnalyzerTrkPlots/";
  TString makedir = "mkdir -p " + DIR;
  const char *mkDIR = makedir.Data();
  gSystem->Exec(mkDIR);
  TString PRESELDIR = DIR + "PreselectionPlots/";
  TString makedirPreSel = "mkdir -p " + PRESELDIR;
  const char *mkDIRPRESEL = makedirPreSel.Data();
  gSystem->Exec(mkDIRPRESEL);

  TFile *fout;
  fout = new TFile(outputDir + "output_" + inputFile + ".root", "recreate");
  
  TLegend* l = new TLegend(0.86,0.35,0.98,0.7);
  l->SetFillColor(0);
  l->SetLineColor(0);
  l->SetTextSize(0.04);
  l->SetTextFont(42);
  //std::cout<<"plots"<<std::endl;

  int i = 0;
  for(auto ivar=varCutFlows.cbegin(); ivar!=varCutFlows.cend(); ++ivar){
    for(uint j=0; j<trackType.size(); ++j){
      for(uint k=0; k<plotModifiers.size(); ++k){
	auto start = std::chrono::high_resolution_clock::now();
	l->Clear();
	TH1F* h_trkEff[preselCuts.size()];
	int m = 0;
	for(auto mcut=preselCuts.cbegin(); mcut!=preselCuts.cend(); ++mcut){
	  if(m!=0){
	    //std::cout<<"trkEffOverlay i j m k: "<<i<<" "<<j<<" "<<m<<" "<<k<<std::endl;
	    h_trkEff[m] = (TH1F*)preselCutFlows[i][j][0][k]->Clone();
	    h_trkEff[m]->GetYaxis()->SetNoExponent(kTRUE);
	    removeFlows(h_trkEff[m]);
	    h_trkEff[m]->SetStats(0);
	    removeFlows(preselCutFlows[i][j][m][k]);
	    h_trkEff[m]->Divide(preselCutFlows[i][j][m][k],h_trkEff[m],1.0,1.0,"B");
	    h_trkEff[m]->SetLineColor(m);
	    h_trkEff[m]->SetMarkerColor(m);
	    TString cutName = mcut->first;
	    //std::cout<<"cutName: "<<cutName<<std::endl;
	    l->AddEntry(h_trkEff[m],cutName,"lp");
	    if(m==1){
	      raiseMax(h_trkEff[m]);
	      h_trkEff[m]->Draw();
	    }
	    else{
	      h_trkEff[m]->Draw("SAME");
	    }
	  }
	  m++;
	}
	mySmallText(0.3, 0.9, 1, ctxt);
	l->Draw();
	c.SaveAs(PRESELDIR + "/h_trkEffOverlay_"+ivar->first.at(0)+"_"+trackType[j]+plotModifiers[k]+".pdf");
	auto finish = std::chrono::high_resolution_clock::now();
	/*std::cout << "trkEffOverlay took "
		  << std::chrono::duration_cast<milli>(finish - start).count()
		  << " milliseconds\n";*/
      }
    }
    i++;
  }

  i = 0;
  int m_primary = distance(trackType.begin(), find(trackType.begin(), trackType.end(), "primary"));
  int m_np = distance(trackType.begin(), find(trackType.begin(), trackType.end(), "np"));
  int m_fake = distance(trackType.begin(), find(trackType.begin(), trackType.end(), "fake"));
  int m_PU = distance(trackType.begin(), find(trackType.begin(), trackType.end(), "PU"));
  int m_notHiggs = distance(trackType.begin(), find(trackType.begin(), trackType.end(), "notHiggs"));

  for(auto icut=preselCuts.cbegin(); icut!=preselCuts.cend(); ++icut){
    for(uint j=0; j<plotModifiers.size(); ++j){
      int k = 0;
      for(auto kvar=varCutFlows.cbegin(); kvar!=varCutFlows.cend(); ++kvar){
	TString cutName = icut->first;
	auto h_stack = new THStack("hs_"+kvar->first.at(0)+"_"+cutName+"Cut"+plotModifiers[j],"Stacked BG histograms");
	float integralSum = 0;
	l->Clear();
	for(uint m=0; m<trackType.size(); ++m){
	  auto start = std::chrono::high_resolution_clock::now();
	  preselCutFlows[k][m][i][j]->GetYaxis()->SetNoExponent(kTRUE);
	  removeFlows(preselCutFlows[k][m][i][j]);
	  if(detailedPlots){
	    raiseMax(preselCutFlows[k][m][i][j]);
	    preselCutFlows[k][m][i][j]->Draw();
	    mySmallText(0.3, 0.9, 1, ctxt);
	    preselCutFlows[k][m][i][j]->Write("", TObject::kOverwrite);
	    c.SaveAs(PRESELDIR + "/"+ preselCutFlows[k][m][i][j]->GetName() + ".pdf");
	  }
	  preselCutFlows[k][m][i][j]->SetLineColor(m+1);
	  preselCutFlows[k][m][i][j]->SetMarkerColor(m+1);
	  if(trackType[m]=="fake" || trackType[m]=="PU" || trackType[m]=="notHiggs"){
	    integralSum+=preselCutFlows[k][m][i][j]->Integral();
	  }
	  auto finish = std::chrono::high_resolution_clock::now();
	  /*std::cout << "preselCutFlows took "
		    << std::chrono::duration_cast<milli>(finish - start).count()
		    << " milliseconds\n";*/
	}
	for(uint m=0; m<trackType.size(); ++m){
	  if(trackType[m]=="fake" || trackType[m]=="PU" || trackType[m]=="notHiggs"){
	    preselCutFlows[k][m][i][j]->Scale(1./integralSum);
	    h_stack->Add(preselCutFlows[k][m][i][j]);
	  }
	}
	auto start = std::chrono::high_resolution_clock::now();
	h_stack->Draw("HIST");
	preselCutFlows[k][m_primary][i][j]->Scale(1./preselCutFlows[k][m_primary][i][j]->Integral());
	raiseMax(preselCutFlows[k][m_primary][i][j],h_stack);
	drawSame(preselCutFlows[k][m_primary][i][j],h_stack);
	mySmallText(0.3, 0.9, 1, ctxt);
	l->Clear();
	l->AddEntry(preselCutFlows[k][m_primary][i][j],"Primary","l");
	l->AddEntry(preselCutFlows[k][m_fake][i][j],"Fake","l");
	l->AddEntry(preselCutFlows[k][m_PU][i][j],"PU","l");
	l->AddEntry(preselCutFlows[k][m_notHiggs][i][j],"notHiggs","l");
	l->Draw();
	c.SaveAs(PRESELDIR + "/h_signalVsBGStack_"+kvar->first.at(0)+"_"+icut->first+"Cut"+plotModifiers[j]+".pdf");
	for(uint m=0; m<trackType.size(); ++m){
	  preselCutFlows[k][m][i][j]->Scale(1./preselCutFlows[k][m][i][j]->Integral());
	  preselCutFlows[k][m][i][j]->SetStats(0);
	}
	raiseMax(preselCutFlows[k][m_primary][i][j],preselCutFlows[k][m_fake][i][j],preselCutFlows[k][m_PU][i][j],preselCutFlows[k][m_notHiggs][i][j]);
	drawSame(preselCutFlows[k][m_primary][i][j],preselCutFlows[k][m_fake][i][j],preselCutFlows[k][m_PU][i][j],preselCutFlows[k][m_notHiggs][i][j]);
	mySmallText(0.3, 0.9, 1, ctxt);
	l->Draw();
	//std::cout<<"signalvsBGOverlay primary fake PU notHiggs: "<<m_primary<<" "<<m_fake<<" "<<m_PU<<" "<<m_notHiggs<<std::endl;
	//std::cout<<"signalvsBGOverlay k i j: "<<k<<" "<<i<<" "<<j<<std::endl;
	c.SaveAs(PRESELDIR + "/h_signalVsBGOverlay_"+kvar->first.at(0)+"_"+icut->first+"Cut"+plotModifiers[j]+".pdf");
	auto finish = std::chrono::high_resolution_clock::now();
	/*std::cout << "signalVsBGOverlay took "
		  << std::chrono::duration_cast<milli>(finish - start).count()
		  << " milliseconds\n";*/
	raiseMax(preselCutFlows[k][m_primary][i][j],preselCutFlows[k][m_np][i][j]);
	drawSame(preselCutFlows[k][m_primary][i][j],preselCutFlows[k][m_np][i][j]);
	mySmallText(0.3, 0.9, 1, ctxt);
	l->Clear();
	l->AddEntry(preselCutFlows[k][m_primary][i][j],"Primary","l");
	l->AddEntry(preselCutFlows[k][m_np][i][j],"NP","l");
	l->Draw();
	c.SaveAs(PRESELDIR + "/h_signalVsBG_"+kvar->first.at(0)+"_"+icut->first+"Cut"+plotModifiers[j]+".pdf");
	
	k++;
      }
    }
    i++;
  }

  i=0;
  int m_match = distance(tpType.begin(), find(tpType.begin(), tpType.end(), "match"));
  int m_tp = distance(tpType.begin(), find(tpType.begin(), tpType.end(), ""));
  for(auto icut=preselCutsTP.cbegin(); icut!=preselCutsTP.cend(); ++icut){
    TString cutName = icut->first;
    for(uint j=0; j<plotModifiers.size(); ++j){
      int k = 0;
      for(auto kvar=varCutFlowsTP.cbegin(); kvar!=varCutFlowsTP.cend(); ++kvar){
	for(uint m=0; m<tpType.size(); ++m){
	  preselCutFlowsTP[k][m][i][j]->GetYaxis()->SetNoExponent(kTRUE);
	  removeFlows(preselCutFlowsTP[k][m][i][j]);
	  if(detailedPlots){
	    raiseMax(preselCutFlowsTP[k][m][i][j]);
	    preselCutFlowsTP[k][m][i][j]->Draw();
	    mySmallText(0.3, 0.9, 1, ctxt);
	    preselCutFlowsTP[k][m][i][j]->Write("", TObject::kOverwrite);
	    c.SaveAs(PRESELDIR + "/"+ preselCutFlowsTP[k][m][i][j]->GetName() + ".pdf");
	  }
	}
	preselCutFlowsTP[k][m_match][i][j]->Divide(preselCutFlowsTP[k][m_match][i][j],preselCutFlowsTP[k][m_tp][i][j]);
	raiseMax(preselCutFlowsTP[k][m_match][i][j]);
	mySmallText(0.3, 0.9, 1, ctxt);
	c.SaveAs(PRESELDIR + "/h_trackFindingEff_"+kvar->first.at(0)+"_"+cutName+"Cut"+plotModifiers[j]+".pdf");
	k++;
      }
    }
    i++;
  }

  i = 0;
  for(auto icut=preselCuts.cbegin(); icut!=preselCuts.cend(); ++icut){
    TString cutName = icut->first;
    for(uint j=0; j<plotModifiers.size(); ++j){
      int k = 0;
      for(auto kvar=varCutFlows2D.cbegin(); kvar!=varCutFlows2D.cend(); ++kvar){
	for(uint m=0; m<trackType.size(); ++m){
	  removeFlows(preselCutFlows2D[k][m][i][j]);
	  preselCutFlows2D[k][m][i][j]->Draw();
	  mySmallText(0.3, 0.9, 1, ctxt);
	  c.SaveAs(PRESELDIR + "/"+ preselCutFlows2D[k][m][i][j]->GetName() + ".pdf");
	}
	k++;
      }
    }
    i++;
  }

  i = 0;
  if(detailedPlots){
    for(auto icut=preselCutsTP.cbegin(); icut!=preselCutsTP.cend(); ++icut){
      TString cutName = icut->first;
      for(uint j=0; j<plotModifiers.size(); ++j){
	int k = 0;
	for(auto kvar=varCutFlowsTP2D.cbegin(); kvar!=varCutFlowsTP2D.cend(); ++kvar){
	  for(uint m=0; m<tpType.size(); ++m){
	    removeFlows(preselCutFlowsTP2D[k][m][i][j]);
	    preselCutFlowsTP2D[k][m][i][j]->Draw();
	    mySmallText(0.3, 0.9, 1, ctxt);
	    c.SaveAs(PRESELDIR + "/"+ preselCutFlowsTP2D[k][m][i][j]->GetName() + ".pdf");
	  }
	  k++;
	}
      }
      i++;
    }
  }

  h_trk_d0->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0);
  h_trk_d0->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_d0->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_d0->GetName() + ".pdf");
  delete h_trk_d0;

  h_trk_d0_primary->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_primary);
  h_trk_d0_primary->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_d0_primary->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_d0_primary->GetName() + ".pdf");
  delete h_trk_d0_primary;

  h_trk_d0_primary_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_primary_noCuts);
  h_trk_d0_primary_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_d0_primary_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_d0_primary_noCuts->GetName() + ".pdf");

  h_trk_d0_np->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_np);
  h_trk_d0_np->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_d0_np->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_d0_np->GetName() + ".pdf");
  delete h_trk_d0_np;

  h_trk_d0_np_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_np_noCuts);
  h_trk_d0_np_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_d0_np_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_d0_np_noCuts->GetName() + ".pdf");
  
  h_trk_d0_primary_noCuts->Scale(1./h_trk_d0_primary_noCuts->Integral());
  h_trk_d0_primary_noCuts->SetLineColor(1);
  h_trk_d0_primary_noCuts->SetMarkerColor(1);
  TH1F *h_trk_d0_np_noCutsNorm = (TH1F*)h_trk_d0_np_noCuts->Clone(); 
  h_trk_d0_np_noCutsNorm->Scale(1./h_trk_d0_np_noCutsNorm->Integral());
  h_trk_d0_np_noCutsNorm->SetLineColor(2);
  h_trk_d0_np_noCutsNorm->SetMarkerColor(2);
  raiseMax(h_trk_d0_np_noCutsNorm,h_trk_d0_primary_noCuts);
  h_trk_d0_np_noCutsNorm->SetStats(0);
  h_trk_d0_primary_noCuts->SetStats(0);
  h_trk_d0_np_noCutsNorm->Draw("HIST");
  h_trk_d0_primary_noCuts->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->SetFillColor(0);
  l->SetLineColor(0);
  l->SetTextSize(0.04);
  l->AddEntry(h_trk_d0_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_d0_np_noCutsNorm,"NP","l");
  l->SetTextFont(42);
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_d0_noCuts.pdf");
  delete h_trk_d0_np_noCutsNorm;

  h_trk_d0_primary_noCuts_zoomOut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_primary_noCuts_zoomOut);
  h_trk_d0_np_noCuts_zoomOut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_np_noCuts_zoomOut);
  h_trk_d0_primary_noCuts_zoomOut->Scale(1./h_trk_d0_primary_noCuts_zoomOut->Integral());
  h_trk_d0_primary_noCuts_zoomOut->SetLineColor(1);
  h_trk_d0_primary_noCuts_zoomOut->SetMarkerColor(1);
  h_trk_d0_np_noCuts_zoomOut->Scale(1./h_trk_d0_np_noCuts_zoomOut->Integral());
  h_trk_d0_np_noCuts_zoomOut->SetLineColor(2);
  h_trk_d0_np_noCuts_zoomOut->SetMarkerColor(2);
  raiseMax(h_trk_d0_np_noCuts_zoomOut,h_trk_d0_primary_noCuts_zoomOut);
  h_trk_d0_np_noCuts_zoomOut->SetStats(0);
  h_trk_d0_primary_noCuts_zoomOut->SetStats(0);
  h_trk_d0_np_noCuts_zoomOut->Draw("HIST");
  h_trk_d0_primary_noCuts_zoomOut->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->SetFillColor(0);
  l->SetLineColor(0);
  l->SetTextSize(0.04);
  l->AddEntry(h_trk_d0_primary_noCuts_zoomOut,"Primary","l");
  l->AddEntry(h_trk_d0_np_noCuts_zoomOut,"NP","l");
  l->SetTextFont(42);
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_d0_noCuts_zoomOut.pdf");
  delete h_trk_d0_primary_noCuts_zoomOut;
  delete h_trk_d0_np_noCuts_zoomOut;

  TH1F *h_trk_d0_fake_noCutsNorm = (TH1F*)h_trk_d0_fake_noCuts->Clone(); 
  h_trk_d0_fake_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_fake_noCutsNorm);
  TH1F *h_trk_d0_PU_noCutsNorm = (TH1F*)h_trk_d0_PU_noCuts->Clone(); 
  h_trk_d0_PU_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_PU_noCutsNorm);
  TH1F *h_trk_d0_notHiggs_noCutsNorm = (TH1F*)h_trk_d0_notHiggs_noCuts->Clone(); 
  h_trk_d0_notHiggs_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_notHiggs_noCutsNorm);
  h_trk_d0_fake_noCutsNorm->Scale(1./h_trk_d0_np_noCuts->Integral());
  h_trk_d0_fake_noCutsNorm->SetLineColor(2);
  h_trk_d0_fake_noCutsNorm->SetMarkerColor(2);
  h_trk_d0_PU_noCutsNorm->Scale(1./h_trk_d0_np_noCuts->Integral());
  h_trk_d0_PU_noCutsNorm->SetLineColor(3);
  h_trk_d0_PU_noCutsNorm->SetMarkerColor(3);
  h_trk_d0_notHiggs_noCutsNorm->Scale(1./h_trk_d0_np_noCuts->Integral());
  h_trk_d0_notHiggs_noCutsNorm->SetLineColor(4);
  h_trk_d0_notHiggs_noCutsNorm->SetMarkerColor(4);
  auto hs_d0_noCuts = new THStack("hs_d0_noCuts","Stacked BG histograms");
  hs_d0_noCuts->Add(h_trk_d0_fake_noCutsNorm);
  hs_d0_noCuts->Add(h_trk_d0_PU_noCutsNorm);
  hs_d0_noCuts->Add(h_trk_d0_notHiggs_noCutsNorm);
  hs_d0_noCuts->Draw("HIST");
  h_trk_d0_primary_noCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_d0_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_d0_fake_noCutsNorm,"Fake","l");
  l->AddEntry(h_trk_d0_PU_noCutsNorm,"PU","l");
  l->AddEntry(h_trk_d0_notHiggs_noCutsNorm,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_d0_noCutsNorm.pdf");
  delete h_trk_d0_fake_noCutsNorm;
  delete h_trk_d0_PU_noCutsNorm;
  delete h_trk_d0_notHiggs_noCutsNorm;
  delete h_trk_d0_np_noCuts;
  delete hs_d0_noCuts;

  h_trk_d0_fake_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_fake_noCuts);
  h_trk_d0_PU_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_PU_noCuts);
  h_trk_d0_notHiggs_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_notHiggs_noCuts);
  h_trk_d0_fake_noCuts->Scale(1./h_trk_d0_fake_noCuts->Integral());
  h_trk_d0_fake_noCuts->SetLineColor(2);
  h_trk_d0_fake_noCuts->SetMarkerColor(2);
  h_trk_d0_PU_noCuts->Scale(1./h_trk_d0_PU_noCuts->Integral());
  h_trk_d0_PU_noCuts->SetLineColor(3);
  h_trk_d0_PU_noCuts->SetMarkerColor(3);
  h_trk_d0_notHiggs_noCuts->Scale(1./h_trk_d0_notHiggs_noCuts->Integral());
  h_trk_d0_notHiggs_noCuts->SetLineColor(4);
  h_trk_d0_notHiggs_noCuts->SetMarkerColor(4);
  h_trk_d0_fake_noCuts->Draw("HIST");
  h_trk_d0_primary_noCuts->Draw("HIST,SAME");
  h_trk_d0_PU_noCuts->Draw("HIST,SAME");
  h_trk_d0_notHiggs_noCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_d0_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_d0_fake_noCuts,"Fake","l");
  l->AddEntry(h_trk_d0_PU_noCuts,"PU","l");
  l->AddEntry(h_trk_d0_notHiggs_noCuts,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_d0_noCuts.pdf");
  delete h_trk_d0_primary_noCuts;
  delete h_trk_d0_fake_noCuts;
  delete h_trk_d0_PU_noCuts;
  delete h_trk_d0_notHiggs_noCuts;

  h_trk_d0_primary_noCuts_H->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_primary_noCuts_H);
  h_trk_d0_np_noCuts_H->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_np_noCuts_H);
  h_trk_d0_primary_noCuts_H->Scale(1./h_trk_d0_primary_noCuts_H->Integral());
  h_trk_d0_primary_noCuts_H->SetLineColor(1);
  h_trk_d0_primary_noCuts_H->SetMarkerColor(1);
  h_trk_d0_np_noCuts_H->Scale(1./h_trk_d0_np_noCuts_H->Integral());
  h_trk_d0_np_noCuts_H->SetLineColor(2);
  h_trk_d0_np_noCuts_H->SetMarkerColor(2);
  h_trk_d0_np_noCuts_H->Draw("HIST");
  h_trk_d0_primary_noCuts_H->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_d0_primary_noCuts_H,"Primary","l");
  l->AddEntry(h_trk_d0_np_noCuts_H,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_d0_noCuts_H.pdf");
  delete h_trk_d0_primary_noCuts_H;
  delete h_trk_d0_np_noCuts_H;

  h_trk_d0_primary_noCuts_L->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_primary_noCuts_L);
  h_trk_d0_np_noCuts_L->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_np_noCuts_L);
  h_trk_d0_primary_noCuts_L->Scale(1./h_trk_d0_primary_noCuts_L->Integral());
  h_trk_d0_primary_noCuts_L->SetLineColor(1);
  h_trk_d0_primary_noCuts_L->SetMarkerColor(1);
  h_trk_d0_np_noCuts_L->Scale(1./h_trk_d0_np_noCuts_L->Integral());
  h_trk_d0_np_noCuts_L->SetLineColor(2);
  h_trk_d0_np_noCuts_L->SetMarkerColor(2);
  h_trk_d0_np_noCuts_L->Draw("HIST");
  h_trk_d0_primary_noCuts_L->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_d0_primary_noCuts_L,"Primary","l");
  l->AddEntry(h_trk_d0_np_noCuts_L,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_d0_noCuts_L.pdf");
  delete h_trk_d0_primary_noCuts_L;
  delete h_trk_d0_np_noCuts_L;

  h_trk_d0_primary_noCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_primary_noCuts_barrel);
  h_trk_d0_np_noCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_np_noCuts_barrel);
  h_trk_d0_primary_noCuts_barrel->Scale(1./h_trk_d0_primary_noCuts_barrel->Integral());
  h_trk_d0_primary_noCuts_barrel->SetLineColor(1);
  h_trk_d0_primary_noCuts_barrel->SetMarkerColor(1);
  h_trk_d0_np_noCuts_barrel->Scale(1./h_trk_d0_np_noCuts_barrel->Integral());
  h_trk_d0_np_noCuts_barrel->SetLineColor(2);
  h_trk_d0_np_noCuts_barrel->SetMarkerColor(2);
  h_trk_d0_np_noCuts_barrel->Draw("HIST");
  h_trk_d0_primary_noCuts_barrel->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_d0_primary_noCuts_barrel,"Primary","l");
  l->AddEntry(h_trk_d0_np_noCuts_barrel,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_d0_noCuts_barrel.pdf");
  delete h_trk_d0_primary_noCuts_barrel;
  delete h_trk_d0_np_noCuts_barrel;

  h_trk_d0_primary_noCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_primary_noCuts_disk);
  h_trk_d0_np_noCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_np_noCuts_disk);
  h_trk_d0_primary_noCuts_disk->Scale(1./h_trk_d0_primary_noCuts_disk->Integral());
  h_trk_d0_primary_noCuts_disk->SetLineColor(1);
  h_trk_d0_primary_noCuts_disk->SetMarkerColor(1);
  h_trk_d0_np_noCuts_disk->Scale(1./h_trk_d0_np_noCuts_disk->Integral());
  h_trk_d0_np_noCuts_disk->SetLineColor(2);
  h_trk_d0_np_noCuts_disk->SetMarkerColor(2);
  h_trk_d0_np_noCuts_disk->Draw("HIST");
  h_trk_d0_primary_noCuts_disk->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_d0_primary_noCuts_disk,"Primary","l");
  l->AddEntry(h_trk_d0_np_noCuts_disk,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_d0_noCuts_disk.pdf");
  delete h_trk_d0_primary_noCuts_disk;
  delete h_trk_d0_np_noCuts_disk;

  h_trk_d0_primary_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_primary_qualCuts);
  h_trk_d0_np_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_np_qualCuts);
  h_trk_d0_primary_qualCuts->Scale(1./h_trk_d0_primary_qualCuts->Integral());
  h_trk_d0_primary_qualCuts->SetLineColor(1);
  h_trk_d0_primary_qualCuts->SetMarkerColor(1);
  TH1F *h_trk_d0_np_qualCutsNorm = (TH1F*)h_trk_d0_np_qualCuts->Clone(); 
  h_trk_d0_np_qualCutsNorm->Scale(1./h_trk_d0_np_qualCutsNorm->Integral());
  h_trk_d0_np_qualCutsNorm->SetLineColor(2);
  h_trk_d0_np_qualCutsNorm->SetMarkerColor(2);
  h_trk_d0_np_qualCutsNorm->Draw("HIST");
  h_trk_d0_primary_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_d0_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_d0_np_qualCutsNorm,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_d0_qualCuts.pdf");
  delete h_trk_d0_np_qualCutsNorm;

  h_trk_d0_primary_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_primary_allCuts);
  h_trk_d0_np_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_np_allCuts);
  h_trk_d0_primary_allCuts->Scale(1./h_trk_d0_primary_allCuts->Integral());
  h_trk_d0_primary_allCuts->SetLineColor(1);
  h_trk_d0_primary_allCuts->SetMarkerColor(1);
  TH1F *h_trk_d0_np_allCutsNorm = (TH1F*)h_trk_d0_np_allCuts->Clone(); 
  h_trk_d0_np_allCutsNorm->Scale(1./h_trk_d0_np_allCutsNorm->Integral());
  h_trk_d0_np_allCutsNorm->SetLineColor(2);
  h_trk_d0_np_allCutsNorm->SetMarkerColor(2);
  raiseMax(h_trk_d0_primary_allCuts,h_trk_d0_np_allCutsNorm);
  h_trk_d0_np_allCutsNorm->SetStats(0);
  h_trk_d0_primary_allCuts->SetStats(0);
  h_trk_d0_primary_allCuts->Draw("HIST");
  h_trk_d0_np_allCutsNorm->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_d0_primary_allCuts,"Primary","l");
  l->AddEntry(h_trk_d0_np_allCutsNorm,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_d0_allCuts.pdf");
  delete h_trk_d0_np_allCutsNorm;

  h_trk_d0_primary_allCuts_zoomOut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_primary_allCuts_zoomOut);
  h_trk_d0_np_allCuts_zoomOut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_np_allCuts_zoomOut);
  h_trk_d0_primary_allCuts_zoomOut->Scale(1./h_trk_d0_primary_allCuts_zoomOut->Integral());
  h_trk_d0_primary_allCuts_zoomOut->SetLineColor(1);
  h_trk_d0_primary_allCuts_zoomOut->SetMarkerColor(1);
  h_trk_d0_np_allCuts_zoomOut->Scale(1./h_trk_d0_np_allCuts_zoomOut->Integral());
  h_trk_d0_np_allCuts_zoomOut->SetLineColor(2);
  h_trk_d0_np_allCuts_zoomOut->SetMarkerColor(2);
  raiseMax(h_trk_d0_np_allCuts_zoomOut,h_trk_d0_primary_allCuts_zoomOut);
  h_trk_d0_np_allCuts_zoomOut->SetStats(0);
  h_trk_d0_primary_allCuts_zoomOut->SetStats(0);
  h_trk_d0_np_allCuts_zoomOut->Draw("HIST");
  h_trk_d0_primary_allCuts_zoomOut->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->SetFillColor(0);
  l->SetLineColor(0);
  l->SetTextSize(0.04);
  l->AddEntry(h_trk_d0_primary_allCuts_zoomOut,"Primary","l");
  l->AddEntry(h_trk_d0_np_allCuts_zoomOut,"NP","l");
  l->SetTextFont(42);
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_d0_allCuts_zoomOut.pdf");
  delete h_trk_d0_primary_allCuts_zoomOut;
  delete h_trk_d0_np_allCuts_zoomOut;

  h_trk_d0_primary_allCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_primary_allCuts_barrel);
  h_trk_d0_np_allCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_np_allCuts_barrel);
  h_trk_d0_primary_allCuts_barrel->Scale(1./h_trk_d0_primary_allCuts_barrel->Integral());
  h_trk_d0_primary_allCuts_barrel->SetLineColor(1);
  h_trk_d0_primary_allCuts_barrel->SetMarkerColor(1);
  h_trk_d0_primary_allCuts_barrel->SetStats(0);
  h_trk_d0_np_allCuts_barrel->Scale(1./h_trk_d0_np_allCuts_barrel->Integral());
  h_trk_d0_np_allCuts_barrel->SetLineColor(2);
  h_trk_d0_np_allCuts_barrel->SetMarkerColor(2);
  h_trk_d0_np_allCuts_barrel->SetStats(0);
  h_trk_d0_np_allCuts_barrel->Draw("HIST");
  h_trk_d0_primary_allCuts_barrel->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_d0_primary_allCuts_barrel,"Primary","l");
  l->AddEntry(h_trk_d0_np_allCuts_barrel,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_d0_allCuts_barrel.pdf");
  delete h_trk_d0_primary_allCuts_barrel;
  delete h_trk_d0_np_allCuts_barrel;

  h_trk_d0_primary_allCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_primary_allCuts_disk);
  h_trk_d0_np_allCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_np_allCuts_disk);
  h_trk_d0_primary_allCuts_disk->Scale(1./h_trk_d0_primary_allCuts_disk->Integral());
  h_trk_d0_primary_allCuts_disk->SetLineColor(1);
  h_trk_d0_primary_allCuts_disk->SetMarkerColor(1);
  h_trk_d0_primary_allCuts_disk->SetStats(0);
  h_trk_d0_np_allCuts_disk->Scale(1./h_trk_d0_np_allCuts_disk->Integral());
  h_trk_d0_np_allCuts_disk->SetLineColor(2);
  h_trk_d0_np_allCuts_disk->SetMarkerColor(2);
  h_trk_d0_np_allCuts_disk->SetStats(0);
  h_trk_d0_np_allCuts_disk->Draw("HIST");
  h_trk_d0_primary_allCuts_disk->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_d0_primary_allCuts_disk,"Primary","l");
  l->AddEntry(h_trk_d0_np_allCuts_disk,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_d0_allCuts_disk.pdf");
  delete h_trk_d0_primary_allCuts_disk;
  delete h_trk_d0_np_allCuts_disk;

  h_trk_pt_vs_d0_primary_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_vs_d0_primary_qualCuts);
  h_trk_pt_vs_d0_np_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_vs_d0_np_qualCuts);
  h_trk_pt_vs_d0_primary_qualCuts->SetLineColor(1);
  h_trk_pt_vs_d0_primary_qualCuts->SetMarkerColor(1);
  h_trk_pt_vs_d0_np_qualCuts->SetLineColor(2);
  h_trk_pt_vs_d0_np_qualCuts->SetMarkerColor(2);
  h_trk_pt_vs_d0_np_qualCuts->SetStats(false);
  h_trk_pt_vs_d0_np_qualCuts->Draw("COLZ");
  auto func = new TF2("fit func","abs(x*pow(y,1.8))<0.7",-2,2,0,30);
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_trk_pt_vs_d0_np_qualCuts.pdf");
  func->Draw("SAME");
  c.SaveAs(DIR + "/h_trk_pt_vs_d0_np_qualCutsFit.pdf");
  h_trk_pt_vs_d0_primary_qualCuts->SetStats(false);
  h_trk_pt_vs_d0_primary_qualCuts->Draw("COLZ");
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_trk_pt_vs_d0_primary_qualCuts.pdf");
  func->Draw("SAME");
  c.SaveAs(DIR + "/h_trk_pt_vs_d0_primary_qualCutsFit.pdf");
  delete h_trk_pt_vs_d0_primary_qualCuts;
  delete h_trk_pt_vs_d0_np_qualCuts;

  h_trk_pt_vs_d0_primary_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_vs_d0_primary_allCuts);
  h_trk_pt_vs_d0_np_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_vs_d0_np_allCuts);
  h_trk_pt_vs_d0_primary_allCuts->SetLineColor(1);
  h_trk_pt_vs_d0_primary_allCuts->SetMarkerColor(1);
  h_trk_pt_vs_d0_np_allCuts->SetLineColor(2);
  h_trk_pt_vs_d0_np_allCuts->SetMarkerColor(2);
  h_trk_pt_vs_d0_np_allCuts->Draw("COLZ");
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_trk_pt_vs_d0_np_allCuts.pdf");
  h_trk_pt_vs_d0_primary_allCuts->Draw("COLZ");
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_trk_pt_vs_d0_primary_allCuts.pdf");
  delete h_trk_pt_vs_d0_primary_allCuts;
  delete h_trk_pt_vs_d0_np_allCuts;

  h_trueVertex_charge_vs_numTPs->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertex_charge_vs_numTPs);
  h_trueVertex_charge_vs_numTPs->SetStats(0);
  c.SetLogz();
  h_trueVertex_charge_vs_numTPs->Draw("COLZ");
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_trueVertex_charge_vs_numTPs.pdf");
  delete h_trueVertex_charge_vs_numTPs;
  c.SetLogz(0);

  h_trueVertex_charge_vs_numTracks->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertex_charge_vs_numTracks);
  h_trueVertex_charge_vs_numTracks->Draw("COLZ");
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_trueVertex_charge_vs_numTracks.pdf");
  delete h_trueVertex_charge_vs_numTracks;

  h_trackVertex_numVertPerTrack->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trackVertex_numVertPerTrack);
  c.SetLogy();
  h_trackVertex_numVertPerTrack->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_trackVertex_numVertPerTrack.pdf");
  c.SetLogy(0);
  delete h_trackVertex_numVertPerTrack;

  h_trackVertex_rankPt_vs_numVertPerTrack->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trackVertex_rankPt_vs_numVertPerTrack);
  h_trackVertex_rankPt_vs_numVertPerTrack->SetStats(0);
  h_trackVertex_rankPt_vs_numVertPerTrack->Draw("COLZ");
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_trackVertex_rankPt_vs_numVertPerTrack.pdf");
  delete h_trackVertex_rankPt_vs_numVertPerTrack;

  h_trackVertex_rankDelz_vs_numVertPerTrack->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trackVertex_rankDelz_vs_numVertPerTrack);
  h_trackVertex_rankDelz_vs_numVertPerTrack->SetStats(0);
  h_trackVertex_rankDelz_vs_numVertPerTrack->Draw("COLZ");
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_trackVertex_rankDelz_vs_numVertPerTrack.pdf");
  delete h_trackVertex_rankDelz_vs_numVertPerTrack;

  h_trackVertex_rankDt_vs_numVertPerTrack->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trackVertex_rankDt_vs_numVertPerTrack);
  h_trackVertex_rankDt_vs_numVertPerTrack->SetStats(0);
  h_trackVertex_rankDt_vs_numVertPerTrack->Draw("COLZ");
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_trackVertex_rankDt_vs_numVertPerTrack.pdf");
  delete h_trackVertex_rankDt_vs_numVertPerTrack;

  h_trackVertex_rankChi2rphidofSum_vs_numVertPerTrack->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trackVertex_rankChi2rphidofSum_vs_numVertPerTrack);
  h_trackVertex_rankChi2rphidofSum_vs_numVertPerTrack->SetStats(0);
  h_trackVertex_rankChi2rphidofSum_vs_numVertPerTrack->Draw("COLZ");
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_trackVertex_rankChi2rphidofSum_vs_numVertPerTrack.pdf");
  delete h_trackVertex_rankChi2rphidofSum_vs_numVertPerTrack;

  h_trackVertex_rankRt_vs_numVertPerTrack->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trackVertex_rankRt_vs_numVertPerTrack);
  h_trackVertex_rankRt_vs_numVertPerTrack->SetStats(0);
  h_trackVertex_rankRt_vs_numVertPerTrack->Draw("COLZ");
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_trackVertex_rankRt_vs_numVertPerTrack.pdf");
  delete h_trackVertex_rankRt_vs_numVertPerTrack;

  h_trueVertex_inTraj->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertex_inTraj);
  h_trueVertex_inTraj->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_trueVertex_inTraj.pdf");
  delete h_trueVertex_inTraj;

  h_trackVertex_inTraj->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trackVertex_inTraj);
  h_trackVertex_inTraj->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_trackVertex_inTraj.pdf");
  delete h_trackVertex_inTraj;

  h_correct_trackVertex_charge_vs_numTracks->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trackVertex_charge_vs_numTracks);
  h_correct_trackVertex_charge_vs_numTracks->Draw("COLZ");
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_correct_trackVertex_charge_vs_numTracks.pdf");
  delete h_correct_trackVertex_charge_vs_numTracks;

  h_false_trackVertex_charge_vs_numTracks->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trackVertex_charge_vs_numTracks);
  h_false_trackVertex_charge_vs_numTracks->Draw("COLZ");
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_false_trackVertex_charge_vs_numTracks.pdf");
  delete h_false_trackVertex_charge_vs_numTracks;

  h_trk_pt_vs_eta_primary_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_vs_eta_primary_qualCuts);
  h_trk_pt_vs_eta_np_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_vs_eta_np_qualCuts);
  h_trk_pt_vs_eta_primary_qualCuts->SetLineColor(1);
  h_trk_pt_vs_eta_primary_qualCuts->SetMarkerColor(1);
  h_trk_pt_vs_eta_np_qualCuts->SetLineColor(2);
  h_trk_pt_vs_eta_np_qualCuts->SetMarkerColor(2);
  h_trk_pt_vs_eta_np_qualCuts->Draw("COLZ");
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_trk_pt_vs_eta_np_qualCuts.pdf");
  h_trk_pt_vs_eta_primary_qualCuts->Draw("COLZ");
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_trk_pt_vs_eta_primary_qualCuts.pdf");
  delete h_trk_pt_vs_eta_primary_qualCuts;
  delete h_trk_pt_vs_eta_np_qualCuts;

  h_trk_pt_vs_eta_primary_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_vs_eta_primary_allCuts);
  h_trk_pt_vs_eta_np_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_vs_eta_np_allCuts);
  h_trk_pt_vs_eta_primary_allCuts->SetLineColor(1);
  h_trk_pt_vs_eta_primary_allCuts->SetMarkerColor(1);
  h_trk_pt_vs_eta_np_allCuts->SetLineColor(2);
  h_trk_pt_vs_eta_np_allCuts->SetMarkerColor(2);
  h_trk_pt_vs_eta_np_allCuts->Draw("COLZ");
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_trk_pt_vs_eta_np_allCuts.pdf");
  h_trk_pt_vs_eta_primary_allCuts->Draw("COLZ");
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_trk_pt_vs_eta_primary_allCuts.pdf");
  delete h_trk_pt_vs_eta_primary_allCuts;
  delete h_trk_pt_vs_eta_np_allCuts;

  h_trk_numStubs_vs_eta_primary_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_numStubs_vs_eta_primary_allCuts);
  h_trk_numStubs_vs_eta_np_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_numStubs_vs_eta_np_allCuts);
  h_trk_numStubs_vs_eta_fake_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_numStubs_vs_eta_fake_allCuts);
  h_trk_numStubs_vs_eta_PU_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_numStubs_vs_eta_PU_allCuts);
  h_trk_numStubs_vs_eta_notHiggs_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_numStubs_vs_eta_notHiggs_allCuts);
  h_trk_numStubs_vs_eta_primary_allCuts->SetLineColor(1);
  h_trk_numStubs_vs_eta_primary_allCuts->SetMarkerColor(1);
  h_trk_numStubs_vs_eta_np_allCuts->SetLineColor(2);
  h_trk_numStubs_vs_eta_np_allCuts->SetMarkerColor(2);
  h_trk_numStubs_vs_eta_np_allCuts->Draw("COLZ");
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_trk_numStubs_vs_eta_np_allCuts.pdf");
  h_trk_numStubs_vs_eta_primary_allCuts->Draw("COLZ");
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_trk_numStubs_vs_eta_primary_allCuts.pdf");
  h_trk_numStubs_vs_eta_fake_allCuts->SetStats(0);
  h_trk_numStubs_vs_eta_fake_allCuts->Draw("COLZ");
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_trk_numStubs_vs_eta_fake_allCuts.pdf");
  h_trk_numStubs_vs_eta_PU_allCuts->Draw("COLZ");
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_trk_numStubs_vs_eta_PU_allCuts.pdf");
  h_trk_numStubs_vs_eta_notHiggs_allCuts->Draw("COLZ");
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_trk_numStubs_vs_eta_notHiggs_allCuts.pdf");
  delete h_trk_numStubs_vs_eta_primary_allCuts;
  delete h_trk_numStubs_vs_eta_np_allCuts;
  delete h_trk_numStubs_vs_eta_fake_allCuts;
  delete h_trk_numStubs_vs_eta_PU_allCuts;
  delete h_trk_numStubs_vs_eta_notHiggs_allCuts;

  h_trk_eta_vs_d0_primary_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta_vs_d0_primary_qualCuts);
  h_trk_eta_vs_d0_np_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta_vs_d0_np_qualCuts);
  h_trk_eta_vs_d0_primary_qualCuts->SetLineColor(1);
  h_trk_eta_vs_d0_primary_qualCuts->SetMarkerColor(1);
  h_trk_eta_vs_d0_np_qualCuts->SetLineColor(2);
  h_trk_eta_vs_d0_np_qualCuts->SetMarkerColor(2);
  h_trk_eta_vs_d0_np_qualCuts->Draw("COLZ");
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_trk_eta_vs_d0_np_qualCuts.pdf");
  h_trk_eta_vs_d0_primary_qualCuts->Draw("COLZ");
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_trk_eta_vs_d0_primary_qualCuts.pdf");
  delete h_trk_eta_vs_d0_primary_qualCuts;
  delete h_trk_eta_vs_d0_np_qualCuts;

  h_trk_eta_vs_d0_primary_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta_vs_d0_primary_allCuts);
  h_trk_eta_vs_d0_np_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta_vs_d0_np_allCuts);
  h_trk_eta_vs_d0_primary_allCuts->SetLineColor(1);
  h_trk_eta_vs_d0_primary_allCuts->SetMarkerColor(1);
  h_trk_eta_vs_d0_np_allCuts->SetLineColor(2);
  h_trk_eta_vs_d0_np_allCuts->SetMarkerColor(2);
  h_trk_eta_vs_d0_np_allCuts->Draw("COLZ");
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_trk_eta_vs_d0_np_allCuts.pdf");
  h_trk_eta_vs_d0_primary_allCuts->Draw("COLZ");
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/h_trk_eta_vs_d0_primary_allCuts.pdf");
  delete h_trk_eta_vs_d0_primary_allCuts;
  delete h_trk_eta_vs_d0_np_allCuts;

  h_trk_ptWeightedD0_primary_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_ptWeightedD0_primary_qualCuts);
  h_trk_ptWeightedD0_np_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_ptWeightedD0_np_qualCuts);
  h_trk_ptWeightedD0_primary_qualCuts->Scale(1./h_trk_ptWeightedD0_primary_qualCuts->Integral());
  h_trk_ptWeightedD0_primary_qualCuts->SetLineColor(1);
  h_trk_ptWeightedD0_primary_qualCuts->SetMarkerColor(1);
  TH1F *h_trk_ptWeightedD0_np_qualCutsNorm = (TH1F*)h_trk_ptWeightedD0_np_qualCuts->Clone(); 
  h_trk_ptWeightedD0_np_qualCutsNorm->Scale(1./h_trk_ptWeightedD0_np_qualCutsNorm->Integral());
  h_trk_ptWeightedD0_np_qualCutsNorm->SetLineColor(2);
  h_trk_ptWeightedD0_np_qualCutsNorm->SetMarkerColor(2);
  h_trk_ptWeightedD0_np_qualCutsNorm->Draw("HIST");
  h_trk_ptWeightedD0_primary_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_ptWeightedD0_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_ptWeightedD0_np_qualCutsNorm,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_ptWeightedD0_qualCuts.pdf");
  delete h_trk_ptWeightedD0_np_qualCutsNorm;

  h_correct_trkAssoc_delPhi->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trkAssoc_delPhi);
  h_false_trkAssoc_delPhi->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trkAssoc_delPhi);
  h_correct_trkAssoc_delPhi->Scale(1./h_correct_trkAssoc_delPhi->Integral());
  h_correct_trkAssoc_delPhi->SetLineColor(1);
  h_correct_trkAssoc_delPhi->SetMarkerColor(1);
  h_false_trkAssoc_delPhi->Scale(1./h_false_trkAssoc_delPhi->Integral());
  h_false_trkAssoc_delPhi->SetLineColor(2);
  h_false_trkAssoc_delPhi->SetMarkerColor(2);
  raiseMax(h_correct_trkAssoc_delPhi,h_false_trkAssoc_delPhi);
  h_correct_trkAssoc_delPhi->Draw("HIST");
  h_false_trkAssoc_delPhi->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_correct_trkAssoc_delPhi,"Correct","l");
  l->AddEntry(h_false_trkAssoc_delPhi,"False","l");
  l->Draw();
  c.SaveAs(DIR + "/h_correctVsFalse_trkAssoc_delPhi.pdf");
  delete h_false_trkAssoc_delPhi;
  delete h_correct_trkAssoc_delPhi;

  h_correct_trkAssoc_delPhiProp->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trkAssoc_delPhiProp);
  h_false_trkAssoc_delPhiProp->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trkAssoc_delPhiProp);
  h_correct_trkAssoc_delPhiProp->Scale(1./h_correct_trkAssoc_delPhiProp->Integral());
  h_correct_trkAssoc_delPhiProp->SetLineColor(1);
  h_correct_trkAssoc_delPhiProp->SetMarkerColor(1);
  h_false_trkAssoc_delPhiProp->Scale(1./h_false_trkAssoc_delPhiProp->Integral());
  h_false_trkAssoc_delPhiProp->SetLineColor(2);
  h_false_trkAssoc_delPhiProp->SetMarkerColor(2);
  raiseMax(h_correct_trkAssoc_delPhiProp,h_false_trkAssoc_delPhiProp);
  h_correct_trkAssoc_delPhiProp->Draw("HIST");
  h_false_trkAssoc_delPhiProp->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_correct_trkAssoc_delPhiProp,"Correct","l");
  l->AddEntry(h_false_trkAssoc_delPhiProp,"False","l");
  l->Draw();
  c.SaveAs(DIR + "/h_correctVsFalse_trkAssoc_delPhiProp.pdf");
  delete h_false_trkAssoc_delPhiProp;
  delete h_correct_trkAssoc_delPhiProp;

  h_correct_trkAssoc_pt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trkAssoc_pt);
  h_false_trkAssoc_pt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trkAssoc_pt);
  h_correct_trkAssoc_pt->Scale(1./h_correct_trkAssoc_pt->Integral());
  h_correct_trkAssoc_pt->SetLineColor(1);
  h_correct_trkAssoc_pt->SetMarkerColor(1);
  h_false_trkAssoc_pt->Scale(1./h_false_trkAssoc_pt->Integral());
  h_false_trkAssoc_pt->SetLineColor(2);
  h_false_trkAssoc_pt->SetMarkerColor(2);
  raiseMax(h_correct_trkAssoc_pt,h_false_trkAssoc_pt);
  h_correct_trkAssoc_pt->Draw("HIST");
  h_false_trkAssoc_pt->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_correct_trkAssoc_pt,"Correct","l");
  l->AddEntry(h_false_trkAssoc_pt,"False","l");
  l->Draw();
  c.SaveAs(DIR + "/h_correctVsFalse_trkAssoc_pt.pdf");
  delete h_false_trkAssoc_pt;
  delete h_correct_trkAssoc_pt;

  h_correct_trkAssoc_delPt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trkAssoc_delPt);
  h_false_trkAssoc_delPt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trkAssoc_delPt);
  h_correct_trkAssoc_delPt->Scale(1./h_correct_trkAssoc_delPt->Integral());
  h_correct_trkAssoc_delPt->SetLineColor(1);
  h_correct_trkAssoc_delPt->SetMarkerColor(1);
  h_false_trkAssoc_delPt->Scale(1./h_false_trkAssoc_delPt->Integral());
  h_false_trkAssoc_delPt->SetLineColor(2);
  h_false_trkAssoc_delPt->SetMarkerColor(2);
  raiseMax(h_correct_trkAssoc_delPt,h_false_trkAssoc_delPt);
  h_correct_trkAssoc_delPt->Draw("HIST");
  h_false_trkAssoc_delPt->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_correct_trkAssoc_delPt,"Correct","l");
  l->AddEntry(h_false_trkAssoc_delPt,"False","l");
  l->Draw();
  c.SaveAs(DIR + "/h_correctVsFalse_trkAssoc_delPt.pdf");
  delete h_false_trkAssoc_delPt;
  delete h_correct_trkAssoc_delPt;

  h_correct_trkAssoc_delxy->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trkAssoc_delxy);
  h_false_trkAssoc_delxy->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trkAssoc_delxy);
  h_correct_trkAssoc_delxy->Scale(1./h_correct_trkAssoc_delxy->Integral());
  h_correct_trkAssoc_delxy->SetLineColor(1);
  h_correct_trkAssoc_delxy->SetMarkerColor(1);
  h_false_trkAssoc_delxy->Scale(1./h_false_trkAssoc_delxy->Integral());
  h_false_trkAssoc_delxy->SetLineColor(2);
  h_false_trkAssoc_delxy->SetMarkerColor(2);
  raiseMax(h_correct_trkAssoc_delxy,h_false_trkAssoc_delxy);
  h_correct_trkAssoc_delxy->Draw("HIST");
  h_false_trkAssoc_delxy->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_correct_trkAssoc_delxy,"Correct","l");
  l->AddEntry(h_false_trkAssoc_delxy,"False","l");
  l->Draw();
  c.SaveAs(DIR + "/h_correctVsFalse_trkAssoc_delxy.pdf");
  delete h_false_trkAssoc_delxy;
  delete h_correct_trkAssoc_delxy;

  h_correct_trkAssoc_delz->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trkAssoc_delz);
  h_false_trkAssoc_delz->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trkAssoc_delz);
  h_correct_trkAssoc_delz->Scale(1./h_correct_trkAssoc_delz->Integral());
  h_correct_trkAssoc_delz->SetLineColor(1);
  h_correct_trkAssoc_delz->SetMarkerColor(1);
  h_false_trkAssoc_delz->Scale(1./h_false_trkAssoc_delz->Integral());
  h_false_trkAssoc_delz->SetLineColor(2);
  h_false_trkAssoc_delz->SetMarkerColor(2);
  raiseMax(h_correct_trkAssoc_delz,h_false_trkAssoc_delz);
  h_correct_trkAssoc_delz->Draw("HIST");
  h_false_trkAssoc_delz->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_correct_trkAssoc_delz,"Correct","l");
  l->AddEntry(h_false_trkAssoc_delz,"False","l");
  l->Draw();
  c.SaveAs(DIR + "/h_correctVsFalse_trkAssoc_delz.pdf");
  delete h_false_trkAssoc_delz;
  delete h_correct_trkAssoc_delz;

  TH1F *h_trk_d0_fake_qualCutsNorm = (TH1F*)h_trk_d0_fake_qualCuts->Clone(); 
  h_trk_d0_fake_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_fake_qualCutsNorm);
  TH1F *h_trk_d0_PU_qualCutsNorm = (TH1F*)h_trk_d0_PU_qualCuts->Clone(); 
  h_trk_d0_PU_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_PU_qualCutsNorm);
  TH1F *h_trk_d0_notHiggs_qualCutsNorm = (TH1F*)h_trk_d0_notHiggs_qualCuts->Clone(); 
  h_trk_d0_notHiggs_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_notHiggs_qualCutsNorm);
  h_trk_d0_fake_qualCutsNorm->Scale(1./h_trk_d0_np_qualCuts->Integral());
  h_trk_d0_fake_qualCutsNorm->SetLineColor(2);
  h_trk_d0_fake_qualCutsNorm->SetMarkerColor(2);
  h_trk_d0_PU_qualCutsNorm->Scale(1./h_trk_d0_np_qualCuts->Integral());
  h_trk_d0_PU_qualCutsNorm->SetLineColor(3);
  h_trk_d0_PU_qualCutsNorm->SetMarkerColor(3);
  h_trk_d0_notHiggs_qualCutsNorm->Scale(1./h_trk_d0_np_qualCuts->Integral());
  h_trk_d0_notHiggs_qualCutsNorm->SetLineColor(4);
  h_trk_d0_notHiggs_qualCutsNorm->SetMarkerColor(4);
  auto hs_d0_qualCuts = new THStack("hs_d0_qualCuts","Stacked BG histograms");
  hs_d0_qualCuts->Add(h_trk_d0_fake_qualCutsNorm);
  hs_d0_qualCuts->Add(h_trk_d0_PU_qualCutsNorm);
  hs_d0_qualCuts->Add(h_trk_d0_notHiggs_qualCutsNorm);
  hs_d0_qualCuts->Draw("HIST");
  h_trk_d0_primary_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_d0_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_d0_fake_qualCutsNorm,"Fake","l");
  l->AddEntry(h_trk_d0_PU_qualCutsNorm,"PU","l");
  l->AddEntry(h_trk_d0_notHiggs_qualCutsNorm,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_d0_qualCutsNorm.pdf");
  delete h_trk_d0_fake_qualCutsNorm;
  delete h_trk_d0_PU_qualCutsNorm;
  delete h_trk_d0_notHiggs_qualCutsNorm;
  delete h_trk_d0_np_qualCuts;
  delete hs_d0_qualCuts;

  h_trk_d0_fake_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_fake_qualCuts);
  h_trk_d0_PU_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_PU_qualCuts);
  h_trk_d0_notHiggs_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_d0_notHiggs_qualCuts);
  h_trk_d0_fake_qualCuts->Scale(1./h_trk_d0_fake_qualCuts->Integral());
  h_trk_d0_fake_qualCuts->SetLineColor(2);
  h_trk_d0_fake_qualCuts->SetMarkerColor(2);
  h_trk_d0_PU_qualCuts->Scale(1./h_trk_d0_PU_qualCuts->Integral());
  h_trk_d0_PU_qualCuts->SetLineColor(3);
  h_trk_d0_PU_qualCuts->SetMarkerColor(3);
  h_trk_d0_notHiggs_qualCuts->Scale(1./h_trk_d0_notHiggs_qualCuts->Integral());
  h_trk_d0_notHiggs_qualCuts->SetLineColor(4);
  h_trk_d0_notHiggs_qualCuts->SetMarkerColor(4);
  h_trk_d0_fake_qualCuts->Draw("HIST");
  h_trk_d0_primary_qualCuts->Draw("HIST,SAME");
  h_trk_d0_PU_qualCuts->Draw("HIST,SAME");
  h_trk_d0_notHiggs_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_d0_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_d0_fake_qualCuts,"Fake","l");
  l->AddEntry(h_trk_d0_PU_qualCuts,"PU","l");
  l->AddEntry(h_trk_d0_notHiggs_qualCuts,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_d0_qualCuts.pdf");
  delete h_trk_d0_primary_qualCuts;
  delete h_trk_d0_fake_qualCuts;
  delete h_trk_d0_PU_qualCuts;
  delete h_trk_d0_notHiggs_qualCuts;

  h_tp_d0->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_tp_d0);
  h_tp_d0->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_tp_d0->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_tp_d0->GetName() + ".pdf");
  delete h_tp_d0;
   
  h_tp_pt_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_tp_pt_noCuts);
  h_tp_pt_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_tp_pt_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_tp_pt_noCuts->GetName() + ".pdf");
  //delete h_tp_pt_noCuts;

  h_tp_eta_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_tp_eta_noCuts);
  h_tp_eta_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_tp_eta_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_tp_eta_noCuts->GetName() + ".pdf");
  //delete h_tp_eta_noCuts;

  h_tp_d0_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_tp_d0_noCuts);
  h_tp_d0_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_tp_d0_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_tp_d0_noCuts->GetName() + ".pdf");
  //delete h_tp_d0_noCuts;

  h_trk_z0->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0);
  h_trk_z0->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_z0->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_z0->GetName() + ".pdf");
  delete h_trk_z0;

  h_trk_z0_primary->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_primary);
  h_trk_z0_primary->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_z0_primary->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_z0_primary->GetName() + ".pdf");
  delete h_trk_z0_primary;

  h_trk_z0_primary_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_primary_noCuts);
  h_trk_z0_primary_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_z0_primary_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_z0_primary_noCuts->GetName() + ".pdf");

  h_trk_z0_np->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_np);
  h_trk_z0_np->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_z0_np->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_z0_np->GetName() + ".pdf");
  delete h_trk_z0_np;

  h_trk_z0_np_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_np_noCuts);
  h_trk_z0_np_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_z0_np_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_z0_np_noCuts->GetName() + ".pdf");

  h_trk_z0_primary_noCuts->Scale(1./h_trk_z0_primary_noCuts->Integral());
  h_trk_z0_primary_noCuts->SetLineColor(1);
  h_trk_z0_primary_noCuts->SetMarkerColor(1);
  TH1F *h_trk_z0_np_noCutsNorm = (TH1F*)h_trk_z0_np_noCuts->Clone(); 
  h_trk_z0_np_noCutsNorm->Scale(1./h_trk_z0_np_noCutsNorm->Integral());
  h_trk_z0_np_noCutsNorm->SetLineColor(2);
  h_trk_z0_np_noCutsNorm->SetMarkerColor(2);
  raiseMax(h_trk_z0_np_noCutsNorm,h_trk_z0_primary_noCuts);
  h_trk_z0_np_noCutsNorm->SetStats(0);
  h_trk_z0_primary_noCuts->SetStats(0);
  h_trk_z0_np_noCutsNorm->Draw("HIST");
  h_trk_z0_primary_noCuts->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_z0_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_z0_np_noCutsNorm,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_z0_noCuts.pdf");
  delete h_trk_z0_np_noCutsNorm;

  h_trk_z0_primary_noCuts_zoomOut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_primary_noCuts_zoomOut);
  h_trk_z0_np_noCuts_zoomOut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_np_noCuts_zoomOut);
  h_trk_z0_primary_noCuts_zoomOut->Scale(1./h_trk_z0_primary_noCuts_zoomOut->Integral());
  h_trk_z0_primary_noCuts_zoomOut->SetLineColor(1);
  h_trk_z0_primary_noCuts_zoomOut->SetMarkerColor(1);
  h_trk_z0_np_noCuts_zoomOut->Scale(1./h_trk_z0_np_noCuts_zoomOut->Integral());
  h_trk_z0_np_noCuts_zoomOut->SetLineColor(2);
  h_trk_z0_np_noCuts_zoomOut->SetMarkerColor(2);
  raiseMax(h_trk_z0_np_noCuts_zoomOut,h_trk_z0_primary_noCuts_zoomOut);
  h_trk_z0_np_noCuts_zoomOut->SetStats(0);
  h_trk_z0_primary_noCuts_zoomOut->SetStats(0);
  h_trk_z0_np_noCuts_zoomOut->Draw("HIST");
  h_trk_z0_primary_noCuts_zoomOut->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->SetFillColor(0);
  l->SetLineColor(0);
  l->SetTextSize(0.04);
  l->AddEntry(h_trk_z0_primary_noCuts_zoomOut,"Primary","l");
  l->AddEntry(h_trk_z0_np_noCuts_zoomOut,"NP","l");
  l->SetTextFont(42);
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_z0_noCuts_zoomOut.pdf");
  delete h_trk_z0_primary_noCuts_zoomOut;
  delete h_trk_z0_np_noCuts_zoomOut;

  TH1F *h_trk_z0_fake_noCutsNorm = (TH1F*)h_trk_z0_fake_noCuts->Clone(); 
  h_trk_z0_fake_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_fake_noCutsNorm);
  TH1F *h_trk_z0_PU_noCutsNorm = (TH1F*)h_trk_z0_PU_noCuts->Clone(); 
  h_trk_z0_PU_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_PU_noCutsNorm);
  TH1F *h_trk_z0_notHiggs_noCutsNorm = (TH1F*)h_trk_z0_notHiggs_noCuts->Clone(); 
  h_trk_z0_notHiggs_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_notHiggs_noCutsNorm);
  h_trk_z0_fake_noCutsNorm->Scale(1./h_trk_z0_np_noCuts->Integral());
  h_trk_z0_fake_noCutsNorm->SetLineColor(2);
  h_trk_z0_fake_noCutsNorm->SetMarkerColor(2);
  h_trk_z0_PU_noCutsNorm->Scale(1./h_trk_z0_np_noCuts->Integral());
  h_trk_z0_PU_noCutsNorm->SetLineColor(3);
  h_trk_z0_PU_noCutsNorm->SetMarkerColor(3);
  h_trk_z0_notHiggs_noCutsNorm->Scale(1./h_trk_z0_np_noCuts->Integral());
  h_trk_z0_notHiggs_noCutsNorm->SetLineColor(4);
  h_trk_z0_notHiggs_noCutsNorm->SetMarkerColor(4);
  auto hs_z0_noCuts = new THStack("hs_z0_noCuts","Stacked BG histograms");
  hs_z0_noCuts->Add(h_trk_z0_fake_noCutsNorm);
  hs_z0_noCuts->Add(h_trk_z0_PU_noCutsNorm);
  hs_z0_noCuts->Add(h_trk_z0_notHiggs_noCutsNorm);
  hs_z0_noCuts->Draw("HIST");
  h_trk_z0_primary_noCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_z0_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_z0_fake_noCutsNorm,"Fake","l");
  l->AddEntry(h_trk_z0_PU_noCutsNorm,"PU","l");
  l->AddEntry(h_trk_z0_notHiggs_noCutsNorm,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_z0_noCutsNorm.pdf");
  delete h_trk_z0_fake_noCutsNorm;
  delete h_trk_z0_PU_noCutsNorm;
  delete h_trk_z0_notHiggs_noCutsNorm;
  delete h_trk_z0_np_noCuts;
  delete hs_z0_noCuts;

  h_trk_z0_fake_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_fake_noCuts);
  h_trk_z0_PU_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_PU_noCuts);
  h_trk_z0_notHiggs_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_notHiggs_noCuts);
  h_trk_z0_fake_noCuts->Scale(1./h_trk_z0_fake_noCuts->Integral());
  h_trk_z0_fake_noCuts->SetLineColor(2);
  h_trk_z0_fake_noCuts->SetMarkerColor(2);
  h_trk_z0_PU_noCuts->Scale(1./h_trk_z0_PU_noCuts->Integral());
  h_trk_z0_PU_noCuts->SetLineColor(3);
  h_trk_z0_PU_noCuts->SetMarkerColor(3);
  h_trk_z0_notHiggs_noCuts->Scale(1./h_trk_z0_notHiggs_noCuts->Integral());
  h_trk_z0_notHiggs_noCuts->SetLineColor(4);
  h_trk_z0_notHiggs_noCuts->SetMarkerColor(4);
  h_trk_z0_fake_noCuts->Draw("HIST");
  h_trk_z0_primary_noCuts->Draw("HIST,SAME");
  h_trk_z0_PU_noCuts->Draw("HIST,SAME");
  h_trk_z0_notHiggs_noCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_z0_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_z0_fake_noCuts,"Fake","l");
  l->AddEntry(h_trk_z0_PU_noCuts,"PU","l");
  l->AddEntry(h_trk_z0_notHiggs_noCuts,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_z0_noCuts.pdf");
  delete h_trk_z0_primary_noCuts;
  delete h_trk_z0_fake_noCuts;
  delete h_trk_z0_PU_noCuts;
  delete h_trk_z0_notHiggs_noCuts;

  h_trk_z0_primary_noCuts_H->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_primary_noCuts_H);
  h_trk_z0_np_noCuts_H->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_np_noCuts_H);
  h_trk_z0_primary_noCuts_H->Scale(1./h_trk_z0_primary_noCuts_H->Integral());
  h_trk_z0_primary_noCuts_H->SetLineColor(1);
  h_trk_z0_primary_noCuts_H->SetMarkerColor(1);
  h_trk_z0_np_noCuts_H->Scale(1./h_trk_z0_np_noCuts_H->Integral());
  h_trk_z0_np_noCuts_H->SetLineColor(2);
  h_trk_z0_np_noCuts_H->SetMarkerColor(2);
  h_trk_z0_np_noCuts_H->Draw("HIST");
  h_trk_z0_primary_noCuts_H->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_z0_primary_noCuts_H,"Primary","l");
  l->AddEntry(h_trk_z0_np_noCuts_H,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_z0_noCuts_H.pdf");
  delete h_trk_z0_primary_noCuts_H;
  delete h_trk_z0_np_noCuts_H;

  h_trk_z0_primary_noCuts_L->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_primary_noCuts_L);
  h_trk_z0_np_noCuts_L->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_np_noCuts_L);
  h_trk_z0_primary_noCuts_L->Scale(1./h_trk_z0_primary_noCuts_L->Integral());
  h_trk_z0_primary_noCuts_L->SetLineColor(1);
  h_trk_z0_primary_noCuts_L->SetMarkerColor(1);
  h_trk_z0_np_noCuts_L->Scale(1./h_trk_z0_np_noCuts_L->Integral());
  h_trk_z0_np_noCuts_L->SetLineColor(2);
  h_trk_z0_np_noCuts_L->SetMarkerColor(2);
  h_trk_z0_np_noCuts_L->Draw("HIST");
  h_trk_z0_primary_noCuts_L->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_z0_primary_noCuts_L,"Primary","l");
  l->AddEntry(h_trk_z0_np_noCuts_L,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_z0_noCuts_L.pdf");
  delete h_trk_z0_primary_noCuts_L;
  delete h_trk_z0_np_noCuts_L;

  h_trk_z0_primary_noCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_primary_noCuts_barrel);
  h_trk_z0_np_noCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_np_noCuts_barrel);
  h_trk_z0_primary_noCuts_barrel->Scale(1./h_trk_z0_primary_noCuts_barrel->Integral());
  h_trk_z0_primary_noCuts_barrel->SetLineColor(1);
  h_trk_z0_primary_noCuts_barrel->SetMarkerColor(1);
  h_trk_z0_np_noCuts_barrel->Scale(1./h_trk_z0_np_noCuts_barrel->Integral());
  h_trk_z0_np_noCuts_barrel->SetLineColor(2);
  h_trk_z0_np_noCuts_barrel->SetMarkerColor(2);
  h_trk_z0_np_noCuts_barrel->Draw("HIST");
  h_trk_z0_primary_noCuts_barrel->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_z0_primary_noCuts_barrel,"Primary","l");
  l->AddEntry(h_trk_z0_np_noCuts_barrel,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_z0_noCuts_barrel.pdf");
  delete h_trk_z0_primary_noCuts_barrel;
  delete h_trk_z0_np_noCuts_barrel;

  h_trk_z0_primary_noCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_primary_noCuts_disk);
  h_trk_z0_np_noCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_np_noCuts_disk);
  h_trk_z0_primary_noCuts_disk->Scale(1./h_trk_z0_primary_noCuts_disk->Integral());
  h_trk_z0_primary_noCuts_disk->SetLineColor(1);
  h_trk_z0_primary_noCuts_disk->SetMarkerColor(1);
  h_trk_z0_np_noCuts_disk->Scale(1./h_trk_z0_np_noCuts_disk->Integral());
  h_trk_z0_np_noCuts_disk->SetLineColor(2);
  h_trk_z0_np_noCuts_disk->SetMarkerColor(2);
  h_trk_z0_np_noCuts_disk->Draw("HIST");
  h_trk_z0_primary_noCuts_disk->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_z0_primary_noCuts_disk,"Primary","l");
  l->AddEntry(h_trk_z0_np_noCuts_disk,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_z0_noCuts_disk.pdf");
  delete h_trk_z0_primary_noCuts_disk;
  delete h_trk_z0_np_noCuts_disk;

  h_trk_z0_primary_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_primary_qualCuts);
  h_trk_z0_np_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_np_qualCuts);
  h_trk_z0_primary_qualCuts->Scale(1./h_trk_z0_primary_qualCuts->Integral());
  h_trk_z0_primary_qualCuts->SetLineColor(1);
  h_trk_z0_primary_qualCuts->SetMarkerColor(1);
  TH1F *h_trk_z0_np_qualCutsNorm = (TH1F*)h_trk_z0_np_qualCuts->Clone(); 
  h_trk_z0_np_qualCutsNorm->Scale(1./h_trk_z0_np_qualCutsNorm->Integral());
  h_trk_z0_np_qualCutsNorm->SetLineColor(2);
  h_trk_z0_np_qualCutsNorm->SetMarkerColor(2);
  h_trk_z0_np_qualCutsNorm->Draw("HIST");
  h_trk_z0_primary_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_z0_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_z0_np_qualCutsNorm,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_z0_qualCuts.pdf");
  delete h_trk_z0_np_qualCutsNorm;

  h_trk_z0_primary_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_primary_allCuts);
  h_trk_z0_np_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_np_allCuts);
  h_trk_z0_primary_allCuts->Scale(1./h_trk_z0_primary_allCuts->Integral());
  h_trk_z0_primary_allCuts->SetLineColor(1);
  h_trk_z0_primary_allCuts->SetMarkerColor(1);
  TH1F *h_trk_z0_np_allCutsNorm = (TH1F*)h_trk_z0_np_allCuts->Clone(); 
  h_trk_z0_np_allCutsNorm->Scale(1./h_trk_z0_np_allCutsNorm->Integral());
  h_trk_z0_np_allCutsNorm->SetLineColor(2);
  h_trk_z0_np_allCutsNorm->SetMarkerColor(2);
  raiseMax(h_trk_z0_primary_allCuts,h_trk_z0_np_allCutsNorm);
  h_trk_z0_np_allCutsNorm->SetStats(0);
  h_trk_z0_primary_allCuts->SetStats(0);
  h_trk_z0_primary_allCuts->Draw("HIST");
  h_trk_z0_np_allCutsNorm->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_z0_primary_allCuts,"Primary","l");
  l->AddEntry(h_trk_z0_np_allCutsNorm,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_z0_allCuts.pdf");
  delete h_trk_z0_np_allCutsNorm;

  h_trk_z0_primary_allCuts_zoomOut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_primary_allCuts_zoomOut);
  h_trk_z0_np_allCuts_zoomOut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_np_allCuts_zoomOut);
  h_trk_z0_primary_allCuts_zoomOut->Scale(1./h_trk_z0_primary_allCuts_zoomOut->Integral());
  h_trk_z0_primary_allCuts_zoomOut->SetLineColor(1);
  h_trk_z0_primary_allCuts_zoomOut->SetMarkerColor(1);
  h_trk_z0_np_allCuts_zoomOut->Scale(1./h_trk_z0_np_allCuts_zoomOut->Integral());
  h_trk_z0_np_allCuts_zoomOut->SetLineColor(2);
  h_trk_z0_np_allCuts_zoomOut->SetMarkerColor(2);
  raiseMax(h_trk_z0_np_allCuts_zoomOut,h_trk_z0_primary_allCuts_zoomOut);
  h_trk_z0_np_allCuts_zoomOut->SetStats(0);
  h_trk_z0_primary_allCuts_zoomOut->SetStats(0);
  h_trk_z0_np_allCuts_zoomOut->Draw("HIST");
  h_trk_z0_primary_allCuts_zoomOut->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->SetFillColor(0);
  l->SetLineColor(0);
  l->SetTextSize(0.04);
  l->AddEntry(h_trk_z0_primary_allCuts_zoomOut,"Primary","l");
  l->AddEntry(h_trk_z0_np_allCuts_zoomOut,"NP","l");
  l->SetTextFont(42);
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_z0_allCuts_zoomOut.pdf");
  delete h_trk_z0_primary_allCuts_zoomOut;
  delete h_trk_z0_np_allCuts_zoomOut;

  h_trk_z0_primary_allCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_primary_allCuts_barrel);
  h_trk_z0_np_allCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_np_allCuts_barrel);
  h_trk_z0_primary_allCuts_barrel->Scale(1./h_trk_z0_primary_allCuts_barrel->Integral());
  h_trk_z0_primary_allCuts_barrel->SetLineColor(1);
  h_trk_z0_primary_allCuts_barrel->SetMarkerColor(1);
  h_trk_z0_np_allCuts_barrel->Scale(1./h_trk_z0_np_allCuts_barrel->Integral());
  h_trk_z0_np_allCuts_barrel->SetLineColor(2);
  h_trk_z0_np_allCuts_barrel->SetMarkerColor(2);
  h_trk_z0_np_allCuts_barrel->Draw("HIST");
  h_trk_z0_primary_allCuts_barrel->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_z0_primary_allCuts_barrel,"Primary","l");
  l->AddEntry(h_trk_z0_np_allCuts_barrel,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_z0_allCuts_barrel.pdf");
  delete h_trk_z0_primary_allCuts_barrel;
  delete h_trk_z0_np_allCuts_barrel;

  h_trk_z0_primary_allCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_primary_allCuts_disk);
  h_trk_z0_np_allCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_np_allCuts_disk);
  h_trk_z0_primary_allCuts_disk->Scale(1./h_trk_z0_primary_allCuts_disk->Integral());
  h_trk_z0_primary_allCuts_disk->SetLineColor(1);
  h_trk_z0_primary_allCuts_disk->SetMarkerColor(1);
  h_trk_z0_np_allCuts_disk->Scale(1./h_trk_z0_np_allCuts_disk->Integral());
  h_trk_z0_np_allCuts_disk->SetLineColor(2);
  h_trk_z0_np_allCuts_disk->SetMarkerColor(2);
  h_trk_z0_np_allCuts_disk->Draw("HIST");
  h_trk_z0_primary_allCuts_disk->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_z0_primary_allCuts_disk,"Primary","l");
  l->AddEntry(h_trk_z0_np_allCuts_disk,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_z0_allCuts_disk.pdf");
  delete h_trk_z0_primary_allCuts_disk;
  delete h_trk_z0_np_allCuts_disk;

  TH1F *h_trk_z0_fake_qualCutsNorm = (TH1F*)h_trk_z0_fake_qualCuts->Clone(); 
  h_trk_z0_fake_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_fake_qualCutsNorm);
  TH1F *h_trk_z0_PU_qualCutsNorm = (TH1F*)h_trk_z0_PU_qualCuts->Clone(); 
  h_trk_z0_PU_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_PU_qualCutsNorm);
  TH1F *h_trk_z0_notHiggs_qualCutsNorm = (TH1F*)h_trk_z0_notHiggs_qualCuts->Clone(); 
  h_trk_z0_notHiggs_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_notHiggs_qualCutsNorm);
  h_trk_z0_fake_qualCutsNorm->Scale(1./h_trk_z0_np_qualCuts->Integral());
  h_trk_z0_fake_qualCutsNorm->SetLineColor(2);
  h_trk_z0_fake_qualCutsNorm->SetMarkerColor(2);
  h_trk_z0_PU_qualCutsNorm->Scale(1./h_trk_z0_np_qualCuts->Integral());
  h_trk_z0_PU_qualCutsNorm->SetLineColor(3);
  h_trk_z0_PU_qualCutsNorm->SetMarkerColor(3);
  h_trk_z0_notHiggs_qualCutsNorm->Scale(1./h_trk_z0_np_qualCuts->Integral());
  h_trk_z0_notHiggs_qualCutsNorm->SetLineColor(4);
  h_trk_z0_notHiggs_qualCutsNorm->SetMarkerColor(4);
  auto hs_z0_qualCuts = new THStack("hs_z0_qualCuts","Stacked BG histograms");
  hs_z0_qualCuts->Add(h_trk_z0_fake_qualCutsNorm);
  hs_z0_qualCuts->Add(h_trk_z0_PU_qualCutsNorm);
  hs_z0_qualCuts->Add(h_trk_z0_notHiggs_qualCutsNorm);
  hs_z0_qualCuts->Draw("HIST");
  h_trk_z0_primary_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_z0_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_z0_fake_qualCutsNorm,"Fake","l");
  l->AddEntry(h_trk_z0_PU_qualCutsNorm,"PU","l");
  l->AddEntry(h_trk_z0_notHiggs_qualCutsNorm,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_z0_qualCutsNorm.pdf");
  delete h_trk_z0_fake_qualCutsNorm;
  delete h_trk_z0_PU_qualCutsNorm;
  delete h_trk_z0_notHiggs_qualCutsNorm;
  delete h_trk_z0_np_qualCuts;
  delete hs_z0_qualCuts;

  h_trk_z0_fake_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_fake_qualCuts);
  h_trk_z0_PU_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_PU_qualCuts);
  h_trk_z0_notHiggs_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_z0_notHiggs_qualCuts);
  h_trk_z0_fake_qualCuts->Scale(1./h_trk_z0_fake_qualCuts->Integral());
  h_trk_z0_fake_qualCuts->SetLineColor(2);
  h_trk_z0_fake_qualCuts->SetMarkerColor(2);
  h_trk_z0_PU_qualCuts->Scale(1./h_trk_z0_PU_qualCuts->Integral());
  h_trk_z0_PU_qualCuts->SetLineColor(3);
  h_trk_z0_PU_qualCuts->SetMarkerColor(3);
  h_trk_z0_notHiggs_qualCuts->Scale(1./h_trk_z0_notHiggs_qualCuts->Integral());
  h_trk_z0_notHiggs_qualCuts->SetLineColor(4);
  h_trk_z0_notHiggs_qualCuts->SetMarkerColor(4);
  h_trk_z0_fake_qualCuts->Draw("HIST");
  h_trk_z0_primary_qualCuts->Draw("HIST,SAME");
  h_trk_z0_PU_qualCuts->Draw("HIST,SAME");
  h_trk_z0_notHiggs_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_z0_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_z0_fake_qualCuts,"Fake","l");
  l->AddEntry(h_trk_z0_PU_qualCuts,"PU","l");
  l->AddEntry(h_trk_z0_notHiggs_qualCuts,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_z0_qualCuts.pdf");
  delete h_trk_z0_primary_qualCuts;
  delete h_trk_z0_fake_qualCuts;
  delete h_trk_z0_PU_qualCuts;
  delete h_trk_z0_notHiggs_qualCuts;

  h_trk_phi->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi);
  h_trk_phi->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_phi->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_phi->GetName() + ".pdf");
  delete h_trk_phi;

  h_trk_phi_primary->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_primary);
  h_trk_phi_primary->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_phi_primary->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_phi_primary->GetName() + ".pdf");
  delete h_trk_phi_primary;

  h_trk_phi_primary_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_primary_noCuts);
  h_trk_phi_primary_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_phi_primary_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_phi_primary_noCuts->GetName() + ".pdf");

  h_trk_phi_np->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_np);
  h_trk_phi_np->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_phi_np->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_phi_np->GetName() + ".pdf");
  delete h_trk_phi_np;

  h_trk_phi_np_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_np_noCuts);
  h_trk_phi_np_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_phi_np_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_phi_np_noCuts->GetName() + ".pdf");

  h_trk_phi_primary_noCuts->Scale(1./h_trk_phi_primary_noCuts->Integral());
  h_trk_phi_primary_noCuts->SetLineColor(1);
  h_trk_phi_primary_noCuts->SetMarkerColor(1);
  TH1F *h_trk_phi_np_noCutsNorm = (TH1F*)h_trk_phi_np_noCuts->Clone(); 
  h_trk_phi_np_noCutsNorm->Scale(1./h_trk_phi_np_noCutsNorm->Integral());
  h_trk_phi_np_noCutsNorm->SetLineColor(2);
  h_trk_phi_np_noCutsNorm->SetMarkerColor(2);
  raiseMax(h_trk_phi_primary_noCuts,h_trk_phi_np_noCutsNorm);
  h_trk_phi_primary_noCuts->SetStats(0);
  h_trk_phi_np_noCutsNorm->SetStats(0);
  h_trk_phi_primary_noCuts->Draw("HIST");
  h_trk_phi_np_noCutsNorm->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_phi_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_phi_np_noCutsNorm,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_phi_noCuts.pdf");
  delete h_trk_phi_np_noCutsNorm;

  TH1F *h_trk_phi_fake_noCutsNorm = (TH1F*)h_trk_phi_fake_noCuts->Clone(); 
  h_trk_phi_fake_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_fake_noCutsNorm);
  TH1F *h_trk_phi_PU_noCutsNorm = (TH1F*)h_trk_phi_PU_noCuts->Clone(); 
  h_trk_phi_PU_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_PU_noCutsNorm);
  TH1F *h_trk_phi_notHiggs_noCutsNorm = (TH1F*)h_trk_phi_notHiggs_noCuts->Clone(); 
  h_trk_phi_notHiggs_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_notHiggs_noCutsNorm);
  h_trk_phi_fake_noCutsNorm->Scale(1./h_trk_phi_np_noCuts->Integral());
  h_trk_phi_fake_noCutsNorm->SetLineColor(2);
  h_trk_phi_fake_noCutsNorm->SetMarkerColor(2);
  h_trk_phi_PU_noCutsNorm->Scale(1./h_trk_phi_np_noCuts->Integral());
  h_trk_phi_PU_noCutsNorm->SetLineColor(3);
  h_trk_phi_PU_noCutsNorm->SetMarkerColor(3);
  h_trk_phi_notHiggs_noCutsNorm->Scale(1./h_trk_phi_np_noCuts->Integral());
  h_trk_phi_notHiggs_noCutsNorm->SetLineColor(4);
  h_trk_phi_notHiggs_noCutsNorm->SetMarkerColor(4);
  auto hs_phi_noCuts = new THStack("hs_phi_noCuts","Stacked BG histograms");
  hs_phi_noCuts->Add(h_trk_phi_fake_noCutsNorm);
  hs_phi_noCuts->Add(h_trk_phi_PU_noCutsNorm);
  hs_phi_noCuts->Add(h_trk_phi_notHiggs_noCutsNorm);
  hs_phi_noCuts->Draw("HIST");
  h_trk_phi_primary_noCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_phi_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_phi_fake_noCutsNorm,"Fake","l");
  l->AddEntry(h_trk_phi_PU_noCutsNorm,"PU","l");
  l->AddEntry(h_trk_phi_notHiggs_noCutsNorm,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_phi_noCutsNorm.pdf");
  delete h_trk_phi_fake_noCutsNorm;
  delete h_trk_phi_PU_noCutsNorm;
  delete h_trk_phi_notHiggs_noCutsNorm;
  delete h_trk_phi_np_noCuts;
  delete hs_phi_noCuts;

  h_trk_phi_fake_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_fake_noCuts);
  h_trk_phi_PU_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_PU_noCuts);
  h_trk_phi_notHiggs_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_notHiggs_noCuts);
  h_trk_phi_fake_noCuts->Scale(1./h_trk_phi_fake_noCuts->Integral());
  h_trk_phi_fake_noCuts->SetLineColor(2);
  h_trk_phi_fake_noCuts->SetMarkerColor(2);
  h_trk_phi_PU_noCuts->Scale(1./h_trk_phi_PU_noCuts->Integral());
  h_trk_phi_PU_noCuts->SetLineColor(3);
  h_trk_phi_PU_noCuts->SetMarkerColor(3);
  h_trk_phi_notHiggs_noCuts->Scale(1./h_trk_phi_notHiggs_noCuts->Integral());
  h_trk_phi_notHiggs_noCuts->SetLineColor(4);
  h_trk_phi_notHiggs_noCuts->SetMarkerColor(4);
  h_trk_phi_fake_noCuts->Draw("HIST");
  h_trk_phi_primary_noCuts->Draw("HIST,SAME");
  h_trk_phi_PU_noCuts->Draw("HIST,SAME");
  h_trk_phi_notHiggs_noCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_phi_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_phi_fake_noCuts,"Fake","l");
  l->AddEntry(h_trk_phi_PU_noCuts,"PU","l");
  l->AddEntry(h_trk_phi_notHiggs_noCuts,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_phi_noCuts.pdf");
  delete h_trk_phi_primary_noCuts;
  delete h_trk_phi_fake_noCuts;
  delete h_trk_phi_PU_noCuts;
  delete h_trk_phi_notHiggs_noCuts;

  h_trk_phi_primary_noCuts_H->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_primary_noCuts_H);
  h_trk_phi_np_noCuts_H->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_np_noCuts_H);
  h_trk_phi_primary_noCuts_H->Scale(1./h_trk_phi_primary_noCuts_H->Integral());
  h_trk_phi_primary_noCuts_H->SetLineColor(1);
  h_trk_phi_primary_noCuts_H->SetMarkerColor(1);
  h_trk_phi_np_noCuts_H->Scale(1./h_trk_phi_np_noCuts_H->Integral());
  h_trk_phi_np_noCuts_H->SetLineColor(2);
  h_trk_phi_np_noCuts_H->SetMarkerColor(2);
  h_trk_phi_np_noCuts_H->Draw("HIST");
  h_trk_phi_primary_noCuts_H->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_phi_primary_noCuts_H,"Primary","l");
  l->AddEntry(h_trk_phi_np_noCuts_H,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_phi_noCuts_H.pdf");
  delete h_trk_phi_primary_noCuts_H;
  delete h_trk_phi_np_noCuts_H;

  h_trk_phi_primary_noCuts_L->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_primary_noCuts_L);
  h_trk_phi_np_noCuts_L->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_np_noCuts_L);
  h_trk_phi_primary_noCuts_L->Scale(1./h_trk_phi_primary_noCuts_L->Integral());
  h_trk_phi_primary_noCuts_L->SetLineColor(1);
  h_trk_phi_primary_noCuts_L->SetMarkerColor(1);
  h_trk_phi_np_noCuts_L->Scale(1./h_trk_phi_np_noCuts_L->Integral());
  h_trk_phi_np_noCuts_L->SetLineColor(2);
  h_trk_phi_np_noCuts_L->SetMarkerColor(2);
  h_trk_phi_np_noCuts_L->Draw("HIST");
  h_trk_phi_primary_noCuts_L->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_phi_primary_noCuts_L,"Primary","l");
  l->AddEntry(h_trk_phi_np_noCuts_L,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_phi_noCuts_L.pdf");
  delete h_trk_phi_primary_noCuts_L;
  delete h_trk_phi_np_noCuts_L;

  h_trk_phi_primary_noCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_primary_noCuts_barrel);
  h_trk_phi_np_noCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_np_noCuts_barrel);
  h_trk_phi_primary_noCuts_barrel->Scale(1./h_trk_phi_primary_noCuts_barrel->Integral());
  h_trk_phi_primary_noCuts_barrel->SetLineColor(1);
  h_trk_phi_primary_noCuts_barrel->SetMarkerColor(1);
  h_trk_phi_np_noCuts_barrel->Scale(1./h_trk_phi_np_noCuts_barrel->Integral());
  h_trk_phi_np_noCuts_barrel->SetLineColor(2);
  h_trk_phi_np_noCuts_barrel->SetMarkerColor(2);
  h_trk_phi_np_noCuts_barrel->Draw("HIST");
  h_trk_phi_primary_noCuts_barrel->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_phi_primary_noCuts_barrel,"Primary","l");
  l->AddEntry(h_trk_phi_np_noCuts_barrel,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_phi_noCuts_barrel.pdf");
  delete h_trk_phi_primary_noCuts_barrel;
  delete h_trk_phi_np_noCuts_barrel;

  h_trk_phi_primary_noCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_primary_noCuts_disk);
  h_trk_phi_np_noCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_np_noCuts_disk);
  h_trk_phi_primary_noCuts_disk->Scale(1./h_trk_phi_primary_noCuts_disk->Integral());
  h_trk_phi_primary_noCuts_disk->SetLineColor(1);
  h_trk_phi_primary_noCuts_disk->SetMarkerColor(1);
  h_trk_phi_np_noCuts_disk->Scale(1./h_trk_phi_np_noCuts_disk->Integral());
  h_trk_phi_np_noCuts_disk->SetLineColor(2);
  h_trk_phi_np_noCuts_disk->SetMarkerColor(2);
  h_trk_phi_np_noCuts_disk->Draw("HIST");
  h_trk_phi_primary_noCuts_disk->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_phi_primary_noCuts_disk,"Primary","l");
  l->AddEntry(h_trk_phi_np_noCuts_disk,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_phi_noCuts_disk.pdf");
  delete h_trk_phi_primary_noCuts_disk;
  delete h_trk_phi_np_noCuts_disk;

  h_trk_phi_primary_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_primary_qualCuts);
  h_trk_phi_np_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_np_qualCuts);
  h_trk_phi_primary_qualCuts->Scale(1./h_trk_phi_primary_qualCuts->Integral());
  h_trk_phi_primary_qualCuts->SetLineColor(1);
  h_trk_phi_primary_qualCuts->SetMarkerColor(1);
  TH1F *h_trk_phi_np_qualCutsNorm = (TH1F*)h_trk_phi_np_qualCuts->Clone(); 
  h_trk_phi_np_qualCutsNorm->Scale(1./h_trk_phi_np_qualCutsNorm->Integral());
  h_trk_phi_np_qualCutsNorm->SetLineColor(2);
  h_trk_phi_np_qualCutsNorm->SetMarkerColor(2);
  h_trk_phi_np_qualCutsNorm->Draw("HIST");
  h_trk_phi_primary_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_phi_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_phi_np_qualCutsNorm,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_phi_qualCuts.pdf");
  delete h_trk_phi_np_qualCutsNorm;

  h_trk_phi_primary_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_primary_allCuts);
  h_trk_phi_np_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_np_allCuts);
  h_trk_phi_primary_allCuts->Scale(1./h_trk_phi_primary_allCuts->Integral());
  h_trk_phi_primary_allCuts->SetLineColor(1);
  h_trk_phi_primary_allCuts->SetMarkerColor(1);
  TH1F *h_trk_phi_np_allCutsNorm = (TH1F*)h_trk_phi_np_allCuts->Clone(); 
  h_trk_phi_np_allCutsNorm->Scale(1./h_trk_phi_np_allCutsNorm->Integral());
  h_trk_phi_np_allCutsNorm->SetLineColor(2);
  h_trk_phi_np_allCutsNorm->SetMarkerColor(2);
  raiseMax(h_trk_phi_primary_allCuts,h_trk_phi_np_allCutsNorm);
  h_trk_phi_np_allCutsNorm->SetStats(0);
  h_trk_phi_primary_allCuts->SetStats(0);
  h_trk_phi_primary_allCuts->Draw("HIST");
  h_trk_phi_np_allCutsNorm->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_phi_primary_allCuts,"Primary","l");
  l->AddEntry(h_trk_phi_np_allCutsNorm,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_phi_allCuts.pdf");
  delete h_trk_phi_np_allCutsNorm;

  h_trk_phi_primary_allCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_primary_allCuts_barrel);
  h_trk_phi_np_allCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_np_allCuts_barrel);
  h_trk_phi_primary_allCuts_barrel->Scale(1./h_trk_phi_primary_allCuts_barrel->Integral());
  h_trk_phi_primary_allCuts_barrel->SetLineColor(1);
  h_trk_phi_primary_allCuts_barrel->SetMarkerColor(1);
  h_trk_phi_np_allCuts_barrel->Scale(1./h_trk_phi_np_allCuts_barrel->Integral());
  h_trk_phi_np_allCuts_barrel->SetLineColor(2);
  h_trk_phi_np_allCuts_barrel->SetMarkerColor(2);
  h_trk_phi_np_allCuts_barrel->Draw("HIST");
  h_trk_phi_primary_allCuts_barrel->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_phi_primary_allCuts_barrel,"Primary","l");
  l->AddEntry(h_trk_phi_np_allCuts_barrel,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_phi_allCuts_barrel.pdf");
  delete h_trk_phi_primary_allCuts_barrel;
  delete h_trk_phi_np_allCuts_barrel;

  h_trk_phi_primary_allCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_primary_allCuts_disk);
  h_trk_phi_np_allCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_np_allCuts_disk);
  h_trk_phi_primary_allCuts_disk->Scale(1./h_trk_phi_primary_allCuts_disk->Integral());
  h_trk_phi_primary_allCuts_disk->SetLineColor(1);
  h_trk_phi_primary_allCuts_disk->SetMarkerColor(1);
  h_trk_phi_np_allCuts_disk->Scale(1./h_trk_phi_np_allCuts_disk->Integral());
  h_trk_phi_np_allCuts_disk->SetLineColor(2);
  h_trk_phi_np_allCuts_disk->SetMarkerColor(2);
  h_trk_phi_np_allCuts_disk->Draw("HIST");
  h_trk_phi_primary_allCuts_disk->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_phi_primary_allCuts_disk,"Primary","l");
  l->AddEntry(h_trk_phi_np_allCuts_disk,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_phi_allCuts_disk.pdf");
  delete h_trk_phi_primary_allCuts_disk;
  delete h_trk_phi_np_allCuts_disk;

  TH1F *h_trk_phi_fake_qualCutsNorm = (TH1F*)h_trk_phi_fake_qualCuts->Clone(); 
  h_trk_phi_fake_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_fake_qualCutsNorm);
  TH1F *h_trk_phi_PU_qualCutsNorm = (TH1F*)h_trk_phi_PU_qualCuts->Clone(); 
  h_trk_phi_PU_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_PU_qualCutsNorm);
  TH1F *h_trk_phi_notHiggs_qualCutsNorm = (TH1F*)h_trk_phi_notHiggs_qualCuts->Clone(); 
  h_trk_phi_notHiggs_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_notHiggs_qualCutsNorm);
  h_trk_phi_fake_qualCutsNorm->Scale(1./h_trk_phi_np_qualCuts->Integral());
  h_trk_phi_fake_qualCutsNorm->SetLineColor(2);
  h_trk_phi_fake_qualCutsNorm->SetMarkerColor(2);
  h_trk_phi_PU_qualCutsNorm->Scale(1./h_trk_phi_np_qualCuts->Integral());
  h_trk_phi_PU_qualCutsNorm->SetLineColor(3);
  h_trk_phi_PU_qualCutsNorm->SetMarkerColor(3);
  h_trk_phi_notHiggs_qualCutsNorm->Scale(1./h_trk_phi_np_qualCuts->Integral());
  h_trk_phi_notHiggs_qualCutsNorm->SetLineColor(4);
  h_trk_phi_notHiggs_qualCutsNorm->SetMarkerColor(4);
  auto hs_phi_qualCuts = new THStack("hs_phi_qualCuts","Stacked BG histograms");
  hs_phi_qualCuts->Add(h_trk_phi_fake_qualCutsNorm);
  hs_phi_qualCuts->Add(h_trk_phi_PU_qualCutsNorm);
  hs_phi_qualCuts->Add(h_trk_phi_notHiggs_qualCutsNorm);
  hs_phi_qualCuts->Draw("HIST");
  h_trk_phi_primary_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_phi_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_phi_fake_qualCutsNorm,"Fake","l");
  l->AddEntry(h_trk_phi_PU_qualCutsNorm,"PU","l");
  l->AddEntry(h_trk_phi_notHiggs_qualCutsNorm,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_phi_qualCutsNorm.pdf");
  delete h_trk_phi_fake_qualCutsNorm;
  delete h_trk_phi_PU_qualCutsNorm;
  delete h_trk_phi_notHiggs_qualCutsNorm;
  delete h_trk_phi_np_qualCuts;
  delete hs_phi_qualCuts;

  h_trk_phi_fake_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_fake_qualCuts);
  h_trk_phi_PU_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_PU_qualCuts);
  h_trk_phi_notHiggs_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_phi_notHiggs_qualCuts);
  h_trk_phi_fake_qualCuts->Scale(1./h_trk_phi_fake_qualCuts->Integral());
  h_trk_phi_fake_qualCuts->SetLineColor(2);
  h_trk_phi_fake_qualCuts->SetMarkerColor(2);
  h_trk_phi_PU_qualCuts->Scale(1./h_trk_phi_PU_qualCuts->Integral());
  h_trk_phi_PU_qualCuts->SetLineColor(3);
  h_trk_phi_PU_qualCuts->SetMarkerColor(3);
  h_trk_phi_notHiggs_qualCuts->Scale(1./h_trk_phi_notHiggs_qualCuts->Integral());
  h_trk_phi_notHiggs_qualCuts->SetLineColor(4);
  h_trk_phi_notHiggs_qualCuts->SetMarkerColor(4);
  h_trk_phi_fake_qualCuts->Draw("HIST");
  h_trk_phi_primary_qualCuts->Draw("HIST,SAME");
  h_trk_phi_PU_qualCuts->Draw("HIST,SAME");
  h_trk_phi_notHiggs_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_phi_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_phi_fake_qualCuts,"Fake","l");
  l->AddEntry(h_trk_phi_PU_qualCuts,"PU","l");
  l->AddEntry(h_trk_phi_notHiggs_qualCuts,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_phi_qualCuts.pdf");
  delete h_trk_phi_primary_qualCuts;
  delete h_trk_phi_fake_qualCuts;
  delete h_trk_phi_PU_qualCuts;
  delete h_trk_phi_notHiggs_qualCuts;

  h_trk_sectorPhi->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi);
  h_trk_sectorPhi->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_sectorPhi->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_sectorPhi->GetName() + ".pdf");
  delete h_trk_sectorPhi;

  h_trk_sectorPhi_primary->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_primary);
  h_trk_sectorPhi_primary->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_sectorPhi_primary->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_sectorPhi_primary->GetName() + ".pdf");
  delete h_trk_sectorPhi_primary;

  h_trk_sectorPhi_primary_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_primary_noCuts);
  h_trk_sectorPhi_primary_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_sectorPhi_primary_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_sectorPhi_primary_noCuts->GetName() + ".pdf");

  h_trk_sectorPhi_np->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_np);
  h_trk_sectorPhi_np->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_sectorPhi_np->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_sectorPhi_np->GetName() + ".pdf");
  delete h_trk_sectorPhi_np;

  h_trk_sectorPhi_np_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_np_noCuts);
  h_trk_sectorPhi_np_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_sectorPhi_np_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_sectorPhi_np_noCuts->GetName() + ".pdf");

  h_trk_sectorPhi_primary_noCuts->Scale(1./h_trk_sectorPhi_primary_noCuts->Integral());
  h_trk_sectorPhi_primary_noCuts->SetLineColor(1);
  h_trk_sectorPhi_primary_noCuts->SetMarkerColor(1);
  TH1F *h_trk_sectorPhi_np_noCutsNorm = (TH1F*)h_trk_sectorPhi_np_noCuts->Clone();
  h_trk_sectorPhi_np_noCutsNorm->Scale(1./h_trk_sectorPhi_np_noCutsNorm->Integral());
  h_trk_sectorPhi_np_noCutsNorm->SetLineColor(2);
  h_trk_sectorPhi_np_noCutsNorm->SetMarkerColor(2);
  raiseMax(h_trk_sectorPhi_primary_noCuts,h_trk_sectorPhi_np_noCutsNorm);
  h_trk_sectorPhi_primary_noCuts->SetStats(0);
  h_trk_sectorPhi_np_noCutsNorm->SetStats(0);
  h_trk_sectorPhi_primary_noCuts->Draw("HIST");
  h_trk_sectorPhi_np_noCutsNorm->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_sectorPhi_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_sectorPhi_np_noCutsNorm,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_sectorPhi_noCuts.pdf");
  delete h_trk_sectorPhi_np_noCutsNorm;

  TH1F *h_trk_sectorPhi_fake_noCutsNorm = (TH1F*)h_trk_sectorPhi_fake_noCuts->Clone(); 
  h_trk_sectorPhi_fake_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_fake_noCutsNorm);
  TH1F *h_trk_sectorPhi_PU_noCutsNorm = (TH1F*)h_trk_sectorPhi_PU_noCuts->Clone(); 
  h_trk_sectorPhi_PU_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_PU_noCutsNorm);
  TH1F *h_trk_sectorPhi_notHiggs_noCutsNorm = (TH1F*)h_trk_sectorPhi_notHiggs_noCuts->Clone(); 
  h_trk_sectorPhi_notHiggs_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_notHiggs_noCutsNorm);
  h_trk_sectorPhi_fake_noCutsNorm->Scale(1./h_trk_sectorPhi_np_noCuts->Integral());
  h_trk_sectorPhi_fake_noCutsNorm->SetLineColor(2);
  h_trk_sectorPhi_fake_noCutsNorm->SetMarkerColor(2);
  h_trk_sectorPhi_PU_noCutsNorm->Scale(1./h_trk_sectorPhi_np_noCuts->Integral());
  h_trk_sectorPhi_PU_noCutsNorm->SetLineColor(3);
  h_trk_sectorPhi_PU_noCutsNorm->SetMarkerColor(3);
  h_trk_sectorPhi_notHiggs_noCutsNorm->Scale(1./h_trk_sectorPhi_np_noCuts->Integral());
  h_trk_sectorPhi_notHiggs_noCutsNorm->SetLineColor(4);
  h_trk_sectorPhi_notHiggs_noCutsNorm->SetMarkerColor(4);
  auto hs_sectorPhi_noCuts = new THStack("hs_sectorPhi_noCuts","Stacked BG histograms");
  hs_sectorPhi_noCuts->Add(h_trk_sectorPhi_fake_noCutsNorm);
  hs_sectorPhi_noCuts->Add(h_trk_sectorPhi_PU_noCutsNorm);
  hs_sectorPhi_noCuts->Add(h_trk_sectorPhi_notHiggs_noCutsNorm);
  hs_sectorPhi_noCuts->Draw("HIST");
  h_trk_sectorPhi_primary_noCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_sectorPhi_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_sectorPhi_fake_noCutsNorm,"Fake","l");
  l->AddEntry(h_trk_sectorPhi_PU_noCutsNorm,"PU","l");
  l->AddEntry(h_trk_sectorPhi_notHiggs_noCutsNorm,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_sectorPhi_noCutsNorm.pdf");
  delete h_trk_sectorPhi_fake_noCutsNorm;
  delete h_trk_sectorPhi_PU_noCutsNorm;
  delete h_trk_sectorPhi_notHiggs_noCutsNorm;
  delete h_trk_sectorPhi_np_noCuts;
  delete hs_sectorPhi_noCuts;

  h_trk_sectorPhi_fake_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_fake_noCuts);
  h_trk_sectorPhi_PU_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_PU_noCuts);
  h_trk_sectorPhi_notHiggs_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_notHiggs_noCuts);
  h_trk_sectorPhi_fake_noCuts->Scale(1./h_trk_sectorPhi_fake_noCuts->Integral());
  h_trk_sectorPhi_fake_noCuts->SetLineColor(2);
  h_trk_sectorPhi_fake_noCuts->SetMarkerColor(2);
  h_trk_sectorPhi_PU_noCuts->Scale(1./h_trk_sectorPhi_PU_noCuts->Integral());
  h_trk_sectorPhi_PU_noCuts->SetLineColor(3);
  h_trk_sectorPhi_PU_noCuts->SetMarkerColor(3);
  h_trk_sectorPhi_notHiggs_noCuts->Scale(1./h_trk_sectorPhi_notHiggs_noCuts->Integral());
  h_trk_sectorPhi_notHiggs_noCuts->SetLineColor(4);
  h_trk_sectorPhi_notHiggs_noCuts->SetMarkerColor(4);
  h_trk_sectorPhi_fake_noCuts->Draw("HIST");
  h_trk_sectorPhi_primary_noCuts->Draw("HIST,SAME");
  h_trk_sectorPhi_PU_noCuts->Draw("HIST,SAME");
  h_trk_sectorPhi_notHiggs_noCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_sectorPhi_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_sectorPhi_fake_noCuts,"Fake","l");
  l->AddEntry(h_trk_sectorPhi_PU_noCuts,"PU","l");
  l->AddEntry(h_trk_sectorPhi_notHiggs_noCuts,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_sectorPhi_noCuts.pdf");
  delete h_trk_sectorPhi_primary_noCuts;
  delete h_trk_sectorPhi_fake_noCuts;
  delete h_trk_sectorPhi_PU_noCuts;
  delete h_trk_sectorPhi_notHiggs_noCuts;

  h_trk_sectorPhi_primary_noCuts_H->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_primary_noCuts_H);
  h_trk_sectorPhi_np_noCuts_H->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_np_noCuts_H);
  h_trk_sectorPhi_primary_noCuts_H->Scale(1./h_trk_sectorPhi_primary_noCuts_H->Integral());
  h_trk_sectorPhi_primary_noCuts_H->SetLineColor(1);
  h_trk_sectorPhi_primary_noCuts_H->SetMarkerColor(1);
  h_trk_sectorPhi_np_noCuts_H->Scale(1./h_trk_sectorPhi_np_noCuts_H->Integral());
  h_trk_sectorPhi_np_noCuts_H->SetLineColor(2);
  h_trk_sectorPhi_np_noCuts_H->SetMarkerColor(2);
  h_trk_sectorPhi_np_noCuts_H->Draw("HIST");
  h_trk_sectorPhi_primary_noCuts_H->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_sectorPhi_primary_noCuts_H,"Primary","l");
  l->AddEntry(h_trk_sectorPhi_np_noCuts_H,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_sectorPhi_noCuts_H.pdf");
  delete h_trk_sectorPhi_primary_noCuts_H;
  delete h_trk_sectorPhi_np_noCuts_H;

  h_trk_sectorPhi_primary_noCuts_L->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_primary_noCuts_L);
  h_trk_sectorPhi_np_noCuts_L->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_np_noCuts_L);
  h_trk_sectorPhi_primary_noCuts_L->Scale(1./h_trk_sectorPhi_primary_noCuts_L->Integral());
  h_trk_sectorPhi_primary_noCuts_L->SetLineColor(1);
  h_trk_sectorPhi_primary_noCuts_L->SetMarkerColor(1);
  h_trk_sectorPhi_np_noCuts_L->Scale(1./h_trk_sectorPhi_np_noCuts_L->Integral());
  h_trk_sectorPhi_np_noCuts_L->SetLineColor(2);
  h_trk_sectorPhi_np_noCuts_L->SetMarkerColor(2);
  h_trk_sectorPhi_np_noCuts_L->Draw("HIST");
  h_trk_sectorPhi_primary_noCuts_L->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_sectorPhi_primary_noCuts_L,"Primary","l");
  l->AddEntry(h_trk_sectorPhi_np_noCuts_L,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_sectorPhi_noCuts_L.pdf");
  delete h_trk_sectorPhi_primary_noCuts_L;
  delete h_trk_sectorPhi_np_noCuts_L;

  h_trk_sectorPhi_primary_noCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_primary_noCuts_barrel);
  h_trk_sectorPhi_np_noCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_np_noCuts_barrel);
  h_trk_sectorPhi_primary_noCuts_barrel->Scale(1./h_trk_sectorPhi_primary_noCuts_barrel->Integral());
  h_trk_sectorPhi_primary_noCuts_barrel->SetLineColor(1);
  h_trk_sectorPhi_primary_noCuts_barrel->SetMarkerColor(1);
  h_trk_sectorPhi_np_noCuts_barrel->Scale(1./h_trk_sectorPhi_np_noCuts_barrel->Integral());
  h_trk_sectorPhi_np_noCuts_barrel->SetLineColor(2);
  h_trk_sectorPhi_np_noCuts_barrel->SetMarkerColor(2);
  h_trk_sectorPhi_np_noCuts_barrel->Draw("HIST");
  h_trk_sectorPhi_primary_noCuts_barrel->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_sectorPhi_primary_noCuts_barrel,"Primary","l");
  l->AddEntry(h_trk_sectorPhi_np_noCuts_barrel,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_sectorPhi_noCuts_barrel.pdf");
  delete h_trk_sectorPhi_primary_noCuts_barrel;
  delete h_trk_sectorPhi_np_noCuts_barrel;

  h_trk_sectorPhi_primary_noCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_primary_noCuts_disk);
  h_trk_sectorPhi_np_noCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_np_noCuts_disk);
  h_trk_sectorPhi_primary_noCuts_disk->Scale(1./h_trk_sectorPhi_primary_noCuts_disk->Integral());
  h_trk_sectorPhi_primary_noCuts_disk->SetLineColor(1);
  h_trk_sectorPhi_primary_noCuts_disk->SetMarkerColor(1);
  h_trk_sectorPhi_np_noCuts_disk->Scale(1./h_trk_sectorPhi_np_noCuts_disk->Integral());
  h_trk_sectorPhi_np_noCuts_disk->SetLineColor(2);
  h_trk_sectorPhi_np_noCuts_disk->SetMarkerColor(2);
  h_trk_sectorPhi_np_noCuts_disk->Draw("HIST");
  h_trk_sectorPhi_primary_noCuts_disk->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_sectorPhi_primary_noCuts_disk,"Primary","l");
  l->AddEntry(h_trk_sectorPhi_np_noCuts_disk,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_sectorPhi_noCuts_disk.pdf");
  delete h_trk_sectorPhi_primary_noCuts_disk;
  delete h_trk_sectorPhi_np_noCuts_disk;

  h_trk_sectorPhi_primary_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_primary_qualCuts);
  h_trk_sectorPhi_np_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_np_qualCuts);
  h_trk_sectorPhi_primary_qualCuts->Scale(1./h_trk_sectorPhi_primary_qualCuts->Integral());
  h_trk_sectorPhi_primary_qualCuts->SetLineColor(1);
  h_trk_sectorPhi_primary_qualCuts->SetMarkerColor(1);
  TH1F *h_trk_sectorPhi_np_qualCutsNorm = (TH1F*)h_trk_sectorPhi_np_qualCuts->Clone(); 
  h_trk_sectorPhi_np_qualCutsNorm->Scale(1./h_trk_sectorPhi_np_qualCutsNorm->Integral());
  h_trk_sectorPhi_np_qualCutsNorm->SetLineColor(2);
  h_trk_sectorPhi_np_qualCutsNorm->SetMarkerColor(2);
  h_trk_sectorPhi_np_qualCutsNorm->Draw("HIST");
  h_trk_sectorPhi_primary_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_sectorPhi_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_sectorPhi_np_qualCutsNorm,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_sectorPhi_qualCuts.pdf");
  delete h_trk_sectorPhi_np_qualCutsNorm;

  h_trk_sectorPhi_primary_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_primary_allCuts);
  h_trk_sectorPhi_np_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_np_allCuts);
  h_trk_sectorPhi_primary_allCuts->Scale(1./h_trk_sectorPhi_primary_allCuts->Integral());
  h_trk_sectorPhi_primary_allCuts->SetLineColor(1);
  h_trk_sectorPhi_primary_allCuts->SetMarkerColor(1);
  TH1F *h_trk_sectorPhi_np_allCutsNorm = (TH1F*)h_trk_sectorPhi_np_allCuts->Clone(); 
  h_trk_sectorPhi_np_allCutsNorm->Scale(1./h_trk_sectorPhi_np_allCutsNorm->Integral());
  h_trk_sectorPhi_np_allCutsNorm->SetLineColor(2);
  h_trk_sectorPhi_np_allCutsNorm->SetMarkerColor(2);
  raiseMax(h_trk_sectorPhi_primary_allCuts,h_trk_sectorPhi_np_allCutsNorm);
  h_trk_sectorPhi_np_allCutsNorm->SetStats(0);
  h_trk_sectorPhi_primary_allCuts->SetStats(0);
  h_trk_sectorPhi_primary_allCuts->Draw("HIST");
  h_trk_sectorPhi_np_allCutsNorm->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_sectorPhi_primary_allCuts,"Primary","l");
  l->AddEntry(h_trk_sectorPhi_np_allCutsNorm,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_sectorPhi_allCuts.pdf");
  delete h_trk_sectorPhi_np_allCutsNorm;

  h_trk_sectorPhi_primary_allCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_primary_allCuts_barrel);
  h_trk_sectorPhi_np_allCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_np_allCuts_barrel);
  h_trk_sectorPhi_primary_allCuts_barrel->Scale(1./h_trk_sectorPhi_primary_allCuts_barrel->Integral());
  h_trk_sectorPhi_primary_allCuts_barrel->SetLineColor(1);
  h_trk_sectorPhi_primary_allCuts_barrel->SetMarkerColor(1);
  h_trk_sectorPhi_np_allCuts_barrel->Scale(1./h_trk_sectorPhi_np_allCuts_barrel->Integral());
  h_trk_sectorPhi_np_allCuts_barrel->SetLineColor(2);
  h_trk_sectorPhi_np_allCuts_barrel->SetMarkerColor(2);
  h_trk_sectorPhi_np_allCuts_barrel->Draw("HIST");
  h_trk_sectorPhi_primary_allCuts_barrel->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_sectorPhi_primary_allCuts_barrel,"Primary","l");
  l->AddEntry(h_trk_sectorPhi_np_allCuts_barrel,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_sectorPhi_allCuts_barrel.pdf");
  delete h_trk_sectorPhi_primary_allCuts_barrel;
  delete h_trk_sectorPhi_np_allCuts_barrel;

  h_trk_sectorPhi_primary_allCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_primary_allCuts_disk);
  h_trk_sectorPhi_np_allCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_np_allCuts_disk);
  h_trk_sectorPhi_primary_allCuts_disk->Scale(1./h_trk_sectorPhi_primary_allCuts_disk->Integral());
  h_trk_sectorPhi_primary_allCuts_disk->SetLineColor(1);
  h_trk_sectorPhi_primary_allCuts_disk->SetMarkerColor(1);
  h_trk_sectorPhi_np_allCuts_disk->Scale(1./h_trk_sectorPhi_np_allCuts_disk->Integral());
  h_trk_sectorPhi_np_allCuts_disk->SetLineColor(2);
  h_trk_sectorPhi_np_allCuts_disk->SetMarkerColor(2);
  h_trk_sectorPhi_np_allCuts_disk->Draw("HIST");
  h_trk_sectorPhi_primary_allCuts_disk->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_sectorPhi_primary_allCuts_disk,"Primary","l");
  l->AddEntry(h_trk_sectorPhi_np_allCuts_disk,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_sectorPhi_allCuts_disk.pdf");
  delete h_trk_sectorPhi_primary_allCuts_disk;
  delete h_trk_sectorPhi_np_allCuts_disk;

  TH1F *h_trk_sectorPhi_fake_qualCutsNorm = (TH1F*)h_trk_sectorPhi_fake_qualCuts->Clone(); 
  h_trk_sectorPhi_fake_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_fake_qualCutsNorm);
  TH1F *h_trk_sectorPhi_PU_qualCutsNorm = (TH1F*)h_trk_sectorPhi_PU_qualCuts->Clone(); 
  h_trk_sectorPhi_PU_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_PU_qualCutsNorm);
  TH1F *h_trk_sectorPhi_notHiggs_qualCutsNorm = (TH1F*)h_trk_sectorPhi_notHiggs_qualCuts->Clone(); 
  h_trk_sectorPhi_notHiggs_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_notHiggs_qualCutsNorm);
  h_trk_sectorPhi_fake_qualCutsNorm->Scale(1./h_trk_sectorPhi_np_qualCuts->Integral());
  h_trk_sectorPhi_fake_qualCutsNorm->SetLineColor(2);
  h_trk_sectorPhi_fake_qualCutsNorm->SetMarkerColor(2);
  h_trk_sectorPhi_PU_qualCutsNorm->Scale(1./h_trk_sectorPhi_np_qualCuts->Integral());
  h_trk_sectorPhi_PU_qualCutsNorm->SetLineColor(3);
  h_trk_sectorPhi_PU_qualCutsNorm->SetMarkerColor(3);
  h_trk_sectorPhi_notHiggs_qualCutsNorm->Scale(1./h_trk_sectorPhi_np_qualCuts->Integral());
  h_trk_sectorPhi_notHiggs_qualCutsNorm->SetLineColor(4);
  h_trk_sectorPhi_notHiggs_qualCutsNorm->SetMarkerColor(4);
  auto hs_sectorPhi_qualCuts = new THStack("hs_sectorPhi_qualCuts","Stacked BG histograms");
  hs_sectorPhi_qualCuts->Add(h_trk_sectorPhi_fake_qualCutsNorm);
  hs_sectorPhi_qualCuts->Add(h_trk_sectorPhi_PU_qualCutsNorm);
  hs_sectorPhi_qualCuts->Add(h_trk_sectorPhi_notHiggs_qualCutsNorm);
  hs_sectorPhi_qualCuts->Draw("HIST");
  h_trk_sectorPhi_primary_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_sectorPhi_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_sectorPhi_fake_qualCutsNorm,"Fake","l");
  l->AddEntry(h_trk_sectorPhi_PU_qualCutsNorm,"PU","l");
  l->AddEntry(h_trk_sectorPhi_notHiggs_qualCutsNorm,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_sectorPhi_qualCutsNorm.pdf");
  delete h_trk_sectorPhi_fake_qualCutsNorm;
  delete h_trk_sectorPhi_PU_qualCutsNorm;
  delete h_trk_sectorPhi_notHiggs_qualCutsNorm;
  delete h_trk_sectorPhi_np_qualCuts;
  delete hs_sectorPhi_qualCuts;

  h_trk_sectorPhi_fake_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_fake_qualCuts);
  h_trk_sectorPhi_PU_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_PU_qualCuts);
  h_trk_sectorPhi_notHiggs_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_sectorPhi_notHiggs_qualCuts);
  h_trk_sectorPhi_fake_qualCuts->Scale(1./h_trk_sectorPhi_fake_qualCuts->Integral());
  h_trk_sectorPhi_fake_qualCuts->SetLineColor(2);
  h_trk_sectorPhi_fake_qualCuts->SetMarkerColor(2);
  h_trk_sectorPhi_PU_qualCuts->Scale(1./h_trk_sectorPhi_PU_qualCuts->Integral());
  h_trk_sectorPhi_PU_qualCuts->SetLineColor(3);
  h_trk_sectorPhi_PU_qualCuts->SetMarkerColor(3);
  h_trk_sectorPhi_notHiggs_qualCuts->Scale(1./h_trk_sectorPhi_notHiggs_qualCuts->Integral());
  h_trk_sectorPhi_notHiggs_qualCuts->SetLineColor(4);
  h_trk_sectorPhi_notHiggs_qualCuts->SetMarkerColor(4);
  h_trk_sectorPhi_fake_qualCuts->Draw("HIST");
  h_trk_sectorPhi_primary_qualCuts->Draw("HIST,SAME");
  h_trk_sectorPhi_PU_qualCuts->Draw("HIST,SAME");
  h_trk_sectorPhi_notHiggs_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_sectorPhi_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_sectorPhi_fake_qualCuts,"Fake","l");
  l->AddEntry(h_trk_sectorPhi_PU_qualCuts,"PU","l");
  l->AddEntry(h_trk_sectorPhi_notHiggs_qualCuts,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_sectorPhi_qualCuts.pdf");
  delete h_trk_sectorPhi_primary_qualCuts;
  delete h_trk_sectorPhi_fake_qualCuts;
  delete h_trk_sectorPhi_PU_qualCuts;
  delete h_trk_sectorPhi_notHiggs_qualCuts;

  h_tp_z0_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_tp_z0_noCuts);
  h_tp_z0_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_tp_z0_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_tp_z0_noCuts->GetName() + ".pdf");
  delete h_tp_z0_noCuts;
   
  h_tp_pt_noCuts_primary->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_tp_pt_noCuts_primary);
  h_tp_pt_noCuts_primary->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_tp_pt_noCuts_primary->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_tp_pt_noCuts_primary->GetName() + ".pdf");
  //delete h_tp_pt_noCuts_primary;

  h_tp_eta_noCuts_primary->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_tp_eta_noCuts_primary);
  h_tp_eta_noCuts_primary->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_tp_eta_noCuts_primary->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_tp_eta_noCuts_primary->GetName() + ".pdf");
  delete h_tp_eta_noCuts_primary;

  h_tp_d0_noCuts_primary->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_tp_d0_noCuts_primary);
  h_tp_d0_noCuts_primary->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_tp_d0_noCuts_primary->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_tp_d0_noCuts_primary->GetName() + ".pdf");
  delete h_tp_d0_noCuts_primary;

  h_tp_z0_noCuts_primary->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_tp_z0_noCuts_primary);
  h_tp_z0_noCuts_primary->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_tp_z0_noCuts_primary->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_tp_z0_noCuts_primary->GetName() + ".pdf");
  delete h_tp_z0_noCuts_primary;
   
  h_tp_pt_noCuts_PU->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_tp_pt_noCuts_PU);
  h_tp_pt_noCuts_PU->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_tp_pt_noCuts_PU->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_tp_pt_noCuts_PU->GetName() + ".pdf");
  delete h_tp_pt_noCuts_PU;

  h_tp_eta_noCuts_PU->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_tp_eta_noCuts_PU);
  h_tp_eta_noCuts_PU->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_tp_eta_noCuts_PU->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_tp_eta_noCuts_PU->GetName() + ".pdf");
  delete h_tp_eta_noCuts_PU;

  h_tp_d0_noCuts_PU->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_tp_d0_noCuts_PU);
  h_tp_d0_noCuts_PU->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_tp_d0_noCuts_PU->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_tp_d0_noCuts_PU->GetName() + ".pdf");
  delete h_tp_d0_noCuts_PU;

  h_tp_z0_noCuts_PU->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_tp_z0_noCuts_PU);
  h_tp_z0_noCuts_PU->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_tp_z0_noCuts_PU->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_tp_z0_noCuts_PU->GetName() + ".pdf");
  delete h_tp_z0_noCuts_PU;
   
  h_numSelectedTrks->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_numSelectedTrks);
  h_numSelectedTrks->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_numSelectedTrks->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_numSelectedTrks->GetName() + ".pdf");
  h_numSelectedTrks->SetStats(0);
  h_numSelectedTrks->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/"+ h_numSelectedTrks->GetName() + "_noStatBox.pdf");
  delete h_numSelectedTrks;

  h_numSelectedTrks_zoomOut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_numSelectedTrks_zoomOut);
  h_numSelectedTrks_zoomOut->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_numSelectedTrks_zoomOut->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_numSelectedTrks_zoomOut->GetName() + ".pdf");
  h_numSelectedTrks_zoomOut->SetStats(0);
  h_numSelectedTrks_zoomOut->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  c.SaveAs(DIR + "/"+ h_numSelectedTrks_zoomOut->GetName() + "_noStatBox.pdf");
  delete h_numSelectedTrks_zoomOut;

  h_trk_H_T->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_H_T);
  h_trk_H_T->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_H_T->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_H_T->GetName() + ".pdf");
  delete h_trk_H_T;

  h_trk_oneMatch_H_T->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_oneMatch_H_T);
  h_trk_oneMatch_H_T->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_oneMatch_H_T->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_oneMatch_H_T->GetName() + ".pdf");
  delete h_trk_oneMatch_H_T;

  h_trk_oneVert_H_T->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_oneVert_H_T);
  h_trk_oneVert_H_T->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_oneVert_H_T->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_oneVert_H_T->GetName() + ".pdf");
  delete h_trk_oneVert_H_T;

  h_trk_MET->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_MET);
  h_trk_MET->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_MET->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_MET->GetName() + ".pdf");
  delete h_trk_MET;

  h_trk_oneMatch_MET->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_oneMatch_MET);
  h_trk_oneMatch_MET->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_oneMatch_MET->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_oneMatch_MET->GetName() + ".pdf");
  delete h_trk_oneMatch_MET;

  h_trk_oneVert_MET->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_oneVert_MET);
  h_trk_oneVert_MET->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_oneVert_MET->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_oneVert_MET->GetName() + ".pdf");
  delete h_trk_oneVert_MET;

  h_tp_H_T->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_tp_H_T);
  h_tp_H_T->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_tp_H_T->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_tp_H_T->GetName() + ".pdf");
  delete h_tp_H_T;

  h_tp_MET->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_tp_MET);
  h_tp_MET->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_tp_MET->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_tp_MET->GetName() + ".pdf");
  delete h_tp_MET;

  h_trk_pt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt);
  h_trk_pt->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_pt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_pt->GetName() + ".pdf");
  delete h_trk_pt;

  h_trk_pt_primary->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_primary);
  h_trk_pt_primary->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_pt_primary->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_pt_primary->GetName() + ".pdf");

  h_trk_pt_primary_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_primary_noCuts);
  h_trk_pt_primary_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_pt_primary_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_pt_primary_noCuts->GetName() + ".pdf");

  h_trk_pt_primary_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_primary_allCuts);
  h_trk_pt_primary_allCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_pt_primary_allCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_pt_primary_allCuts->GetName() + ".pdf");

  h_trk_pt_np->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_np);
  h_trk_pt_np->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_pt_np->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_pt_np->GetName() + ".pdf");

  h_trk_pt_np_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_np_noCuts);
  h_trk_pt_np_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_pt_np_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_pt_np_noCuts->GetName() + ".pdf");

  h_trk_pt_np_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_np_allCuts);
  h_trk_pt_np_allCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_pt_np_allCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_pt_np_allCuts->GetName() + ".pdf");

  TH1F *h_trk_pt_primary_noCutsNorm = (TH1F*)h_trk_pt_primary_noCuts->Clone(); 
  h_trk_pt_primary_noCutsNorm->Scale(1./h_trk_pt_primary_noCutsNorm->Integral());
  h_trk_pt_primary_noCutsNorm->SetLineColor(1);
  h_trk_pt_primary_noCutsNorm->SetMarkerColor(1);
  TH1F *h_trk_pt_np_noCutsNorm = (TH1F*)h_trk_pt_np_noCuts->Clone(); 
  h_trk_pt_np_noCutsNorm->Scale(1./h_trk_pt_np_noCutsNorm->Integral());
  h_trk_pt_np_noCutsNorm->SetLineColor(2);
  h_trk_pt_np_noCutsNorm->SetMarkerColor(2);
  raiseMax(h_trk_pt_np_noCutsNorm,h_trk_pt_primary_noCutsNorm);
  h_trk_pt_np_noCutsNorm->SetStats(0);
  h_trk_pt_primary_noCutsNorm->SetStats(0);
  h_trk_pt_np_noCutsNorm->Draw("HIST");
  h_trk_pt_primary_noCutsNorm->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_pt_primary_noCutsNorm,"Primary","l");
  l->AddEntry(h_trk_pt_np_noCutsNorm,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_pt_noCuts.pdf");
  delete h_trk_pt_np_noCutsNorm;

  h_trk_pt_primary_noCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_primary_noCuts_barrel);
  h_trk_pt_np_noCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_np_noCuts_barrel);
  h_trk_pt_primary_noCuts_barrel->Scale(1./h_trk_pt_primary_noCuts_barrel->Integral());
  h_trk_pt_primary_noCuts_barrel->SetLineColor(1);
  h_trk_pt_primary_noCuts_barrel->SetMarkerColor(1);
  h_trk_pt_np_noCuts_barrel->Scale(1./h_trk_pt_np_noCuts_barrel->Integral());
  h_trk_pt_np_noCuts_barrel->SetLineColor(2);
  h_trk_pt_np_noCuts_barrel->SetMarkerColor(2);
  h_trk_pt_np_noCuts_barrel->Draw("HIST");
  h_trk_pt_primary_noCuts_barrel->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_pt_primary_noCuts_barrel,"Primary","l");
  l->AddEntry(h_trk_pt_np_noCuts_barrel,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_pt_noCuts_barrel.pdf");
  delete h_trk_pt_primary_noCuts_barrel;
  delete h_trk_pt_np_noCuts_barrel;

  h_trk_pt_primary_noCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_primary_noCuts_disk);
  h_trk_pt_np_noCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_np_noCuts_disk);
  h_trk_pt_primary_noCuts_disk->Scale(1./h_trk_pt_primary_noCuts_disk->Integral());
  h_trk_pt_primary_noCuts_disk->SetLineColor(1);
  h_trk_pt_primary_noCuts_disk->SetMarkerColor(1);
  h_trk_pt_np_noCuts_disk->Scale(1./h_trk_pt_np_noCuts_disk->Integral());
  h_trk_pt_np_noCuts_disk->SetLineColor(2);
  h_trk_pt_np_noCuts_disk->SetMarkerColor(2);
  h_trk_pt_np_noCuts_disk->Draw("HIST");
  h_trk_pt_primary_noCuts_disk->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_pt_primary_noCuts_disk,"Primary","l");
  l->AddEntry(h_trk_pt_np_noCuts_disk,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_pt_noCuts_disk.pdf");
  delete h_trk_pt_primary_noCuts_disk;
  delete h_trk_pt_np_noCuts_disk;

  TH1F *h_trk_pt_fake_noCutsNorm = (TH1F*)h_trk_pt_fake_noCuts->Clone(); 
  h_trk_pt_fake_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_fake_noCutsNorm);
  TH1F *h_trk_pt_PU_noCutsNorm = (TH1F*)h_trk_pt_PU_noCuts->Clone(); 
  h_trk_pt_PU_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_PU_noCutsNorm);
  TH1F *h_trk_pt_notHiggs_noCutsNorm = (TH1F*)h_trk_pt_notHiggs_noCuts->Clone(); 
  h_trk_pt_notHiggs_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_notHiggs_noCutsNorm);
  h_trk_pt_fake_noCutsNorm->Scale(1./h_trk_pt_np_noCuts->Integral());
  h_trk_pt_fake_noCutsNorm->SetLineColor(2);
  h_trk_pt_fake_noCutsNorm->SetMarkerColor(2);
  h_trk_pt_PU_noCutsNorm->Scale(1./h_trk_pt_np_noCuts->Integral());
  h_trk_pt_PU_noCutsNorm->SetLineColor(3);
  h_trk_pt_PU_noCutsNorm->SetMarkerColor(3);
  h_trk_pt_notHiggs_noCutsNorm->Scale(1./h_trk_pt_np_noCuts->Integral());
  h_trk_pt_notHiggs_noCutsNorm->SetLineColor(4);
  h_trk_pt_notHiggs_noCutsNorm->SetMarkerColor(4);
  auto hs_pt_noCuts = new THStack("hs_pt_noCuts","Stacked BG histograms");
  hs_pt_noCuts->Add(h_trk_pt_fake_noCutsNorm);
  hs_pt_noCuts->Add(h_trk_pt_PU_noCutsNorm);
  hs_pt_noCuts->Add(h_trk_pt_notHiggs_noCutsNorm);
  hs_pt_noCuts->Draw("HIST");
  h_trk_pt_primary_noCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_pt_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_pt_fake_noCutsNorm,"Fake","l");
  l->AddEntry(h_trk_pt_PU_noCutsNorm,"PU","l");
  l->AddEntry(h_trk_pt_notHiggs_noCutsNorm,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_pt_noCutsNorm.pdf");
  delete h_trk_pt_fake_noCutsNorm;
  delete h_trk_pt_PU_noCutsNorm;
  delete h_trk_pt_notHiggs_noCutsNorm;
  delete hs_pt_noCuts;

  h_trk_pt_fake_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_fake_noCuts);
  h_trk_pt_PU_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_PU_noCuts);
  h_trk_pt_notHiggs_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_notHiggs_noCuts);
  h_trk_pt_fake_noCuts->Scale(1./h_trk_pt_fake_noCuts->Integral());
  h_trk_pt_fake_noCuts->SetLineColor(2);
  h_trk_pt_fake_noCuts->SetMarkerColor(2);
  h_trk_pt_PU_noCuts->Scale(1./h_trk_pt_PU_noCuts->Integral());
  h_trk_pt_PU_noCuts->SetLineColor(3);
  h_trk_pt_PU_noCuts->SetMarkerColor(3);
  h_trk_pt_notHiggs_noCuts->Scale(1./h_trk_pt_notHiggs_noCuts->Integral());
  h_trk_pt_notHiggs_noCuts->SetLineColor(4);
  h_trk_pt_notHiggs_noCuts->SetMarkerColor(4);
  h_trk_pt_fake_noCuts->Draw("HIST");
  h_trk_pt_primary_noCuts->Draw("HIST,SAME");
  h_trk_pt_PU_noCuts->Draw("HIST,SAME");
  h_trk_pt_notHiggs_noCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_pt_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_pt_fake_noCuts,"Fake","l");
  l->AddEntry(h_trk_pt_PU_noCuts,"PU","l");
  l->AddEntry(h_trk_pt_notHiggs_noCuts,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_pt_noCuts.pdf");
  delete h_trk_pt_fake_noCuts;
  delete h_trk_pt_PU_noCuts;
  delete h_trk_pt_notHiggs_noCuts;

  h_trk_pt_primary_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_primary_qualCuts);
  h_trk_pt_np_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_np_qualCuts);
  h_trk_pt_primary_qualCuts->Scale(1./h_trk_pt_primary_qualCuts->Integral());
  h_trk_pt_primary_qualCuts->SetLineColor(1);
  h_trk_pt_primary_qualCuts->SetMarkerColor(1);
  TH1F *h_trk_pt_np_qualCutsNorm = (TH1F*)h_trk_pt_np_qualCuts->Clone(); 
  h_trk_pt_np_qualCutsNorm->Scale(1./h_trk_pt_np_qualCutsNorm->Integral());
  h_trk_pt_np_qualCutsNorm->SetLineColor(2);
  h_trk_pt_np_qualCutsNorm->SetMarkerColor(2);
  h_trk_pt_np_qualCutsNorm->Draw("HIST");
  h_trk_pt_primary_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_pt_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_pt_np_qualCutsNorm,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_pt_qualCuts.pdf");
  delete h_trk_pt_np_qualCutsNorm;

  TH1F *h_trk_pt_fake_qualCutsNorm = (TH1F*)h_trk_pt_fake_qualCuts->Clone(); 
  h_trk_pt_fake_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_fake_qualCutsNorm);
  TH1F *h_trk_pt_PU_qualCutsNorm = (TH1F*)h_trk_pt_PU_qualCuts->Clone(); 
  h_trk_pt_PU_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_PU_qualCutsNorm);
  TH1F *h_trk_pt_notHiggs_qualCutsNorm = (TH1F*)h_trk_pt_notHiggs_qualCuts->Clone(); 
  h_trk_pt_notHiggs_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_notHiggs_qualCutsNorm);
  h_trk_pt_fake_qualCutsNorm->Scale(1./h_trk_pt_np_qualCuts->Integral());
  h_trk_pt_fake_qualCutsNorm->SetLineColor(2);
  h_trk_pt_fake_qualCutsNorm->SetMarkerColor(2);
  h_trk_pt_PU_qualCutsNorm->Scale(1./h_trk_pt_np_qualCuts->Integral());
  h_trk_pt_PU_qualCutsNorm->SetLineColor(3);
  h_trk_pt_PU_qualCutsNorm->SetMarkerColor(3);
  h_trk_pt_notHiggs_qualCutsNorm->Scale(1./h_trk_pt_np_qualCuts->Integral());
  h_trk_pt_notHiggs_qualCutsNorm->SetLineColor(4);
  h_trk_pt_notHiggs_qualCutsNorm->SetMarkerColor(4);
  auto hs_pt_qualCuts = new THStack("hs_pt_qualCuts","Stacked BG histograms");
  hs_pt_qualCuts->Add(h_trk_pt_fake_qualCutsNorm);
  hs_pt_qualCuts->Add(h_trk_pt_PU_qualCutsNorm);
  hs_pt_qualCuts->Add(h_trk_pt_notHiggs_qualCutsNorm);
  hs_pt_qualCuts->Draw("HIST");
  h_trk_pt_primary_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_pt_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_pt_fake_qualCutsNorm,"Fake","l");
  l->AddEntry(h_trk_pt_PU_qualCutsNorm,"PU","l");
  l->AddEntry(h_trk_pt_notHiggs_qualCutsNorm,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_pt_qualCutsNorm.pdf");
  delete h_trk_pt_fake_qualCutsNorm;
  delete h_trk_pt_PU_qualCutsNorm;
  delete h_trk_pt_notHiggs_qualCutsNorm;
  delete h_trk_pt_np_qualCuts;
  delete hs_pt_qualCuts;

  h_trk_pt_fake_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_fake_qualCuts);
  h_trk_pt_PU_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_PU_qualCuts);
  h_trk_pt_notHiggs_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_notHiggs_qualCuts);
  h_trk_pt_fake_qualCuts->Scale(1./h_trk_pt_fake_qualCuts->Integral());
  h_trk_pt_fake_qualCuts->SetLineColor(2);
  h_trk_pt_fake_qualCuts->SetMarkerColor(2);
  h_trk_pt_PU_qualCuts->Scale(1./h_trk_pt_PU_qualCuts->Integral());
  h_trk_pt_PU_qualCuts->SetLineColor(3);
  h_trk_pt_PU_qualCuts->SetMarkerColor(3);
  h_trk_pt_notHiggs_qualCuts->Scale(1./h_trk_pt_notHiggs_qualCuts->Integral());
  h_trk_pt_notHiggs_qualCuts->SetLineColor(4);
  h_trk_pt_notHiggs_qualCuts->SetMarkerColor(4);
  h_trk_pt_fake_qualCuts->Draw("HIST");
  h_trk_pt_primary_qualCuts->Draw("HIST,SAME");
  h_trk_pt_PU_qualCuts->Draw("HIST,SAME");
  h_trk_pt_notHiggs_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_pt_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_pt_fake_qualCuts,"Fake","l");
  l->AddEntry(h_trk_pt_PU_qualCuts,"PU","l");
  l->AddEntry(h_trk_pt_notHiggs_qualCuts,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_pt_qualCuts.pdf");
  delete h_trk_pt_primary_qualCuts;
  delete h_trk_pt_fake_qualCuts;
  delete h_trk_pt_PU_qualCuts;
  delete h_trk_pt_notHiggs_qualCuts;
  
  h_trk_pt_primary_ptCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_primary_ptCuts);
  h_trk_pt_primary_d0Cuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_primary_d0Cuts);
  h_trk_pt_primary_chi2rzdofCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_primary_chi2rzdofCuts);
  h_trk_pt_primary_bendchi2Cuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_primary_bendchi2Cuts);
  h_trk_pt_primary_chi2rphidofCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_primary_chi2rphidofCuts);
  h_trk_pt_primary_nstubCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_primary_nstubCuts);
  h_trk_pt_primary_z0Cuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_primary_z0Cuts);
  h_trk_pt_primary_chi2rzdofCuts->SetName("trkEff_pt_primary");
  h_trk_pt_primary_chi2rzdofCuts->GetXaxis()->SetTitle("Track p_{T} (GeV)");
  h_trk_pt_primary_chi2rzdofCuts->GetYaxis()->SetTitle("Cut Efficiency");
  h_trk_pt_primary_chi2rzdofCuts->Divide(h_trk_pt_primary_chi2rzdofCuts,h_trk_pt_primary_noCuts,1.0,1.0,"B");
  h_trk_pt_primary_chi2rzdofCuts->SetLineColor(1);
  h_trk_pt_primary_chi2rzdofCuts->SetMarkerColor(1);
  h_trk_pt_primary_chi2rzdofCuts->SetStats(0);
  h_trk_pt_primary_chi2rzdofCuts->Draw();
  h_trk_pt_primary_bendchi2Cuts->Divide(h_trk_pt_primary_bendchi2Cuts,h_trk_pt_primary_noCuts,1.0,1.0,"B");
  h_trk_pt_primary_bendchi2Cuts->SetLineColor(2);
  h_trk_pt_primary_bendchi2Cuts->SetMarkerColor(2);
  h_trk_pt_primary_bendchi2Cuts->SetStats(0);
  h_trk_pt_primary_bendchi2Cuts->Draw("SAME");
  h_trk_pt_primary_chi2rphidofCuts->Divide(h_trk_pt_primary_chi2rphidofCuts,h_trk_pt_primary_noCuts,1.0,1.0,"B");
  h_trk_pt_primary_chi2rphidofCuts->SetLineColor(3);
  h_trk_pt_primary_chi2rphidofCuts->SetMarkerColor(3);
  h_trk_pt_primary_chi2rphidofCuts->SetStats(0);
  h_trk_pt_primary_chi2rphidofCuts->Draw("SAME");
  h_trk_pt_primary_nstubCuts->Divide(h_trk_pt_primary_nstubCuts,h_trk_pt_primary_noCuts,1.0,1.0,"B");
  h_trk_pt_primary_nstubCuts->SetLineColor(4);
  h_trk_pt_primary_nstubCuts->SetMarkerColor(4);
  h_trk_pt_primary_nstubCuts->SetStats(0);
  h_trk_pt_primary_nstubCuts->Draw("SAME");
  h_trk_pt_primary_ptCuts->Divide(h_trk_pt_primary_ptCuts,h_trk_pt_primary_noCuts,1.0,1.0,"B");
  h_trk_pt_primary_ptCuts->SetLineColor(5);
  h_trk_pt_primary_ptCuts->SetMarkerColor(5);
  h_trk_pt_primary_ptCuts->SetStats(0);
  h_trk_pt_primary_ptCuts->Draw("SAME");
  h_trk_pt_primary_d0Cuts->Divide(h_trk_pt_primary_d0Cuts,h_trk_pt_primary_noCuts,1.0,1.0,"B");
  h_trk_pt_primary_d0Cuts->SetLineColor(6);
  h_trk_pt_primary_d0Cuts->SetMarkerColor(6);
  h_trk_pt_primary_d0Cuts->SetStats(0);
  h_trk_pt_primary_d0Cuts->Draw("SAME");
  h_trk_pt_primary_z0Cuts->Divide(h_trk_pt_primary_z0Cuts,h_trk_pt_primary_noCuts,1.0,1.0,"B");
  h_trk_pt_primary_z0Cuts->SetLineColor(7);
  h_trk_pt_primary_z0Cuts->SetMarkerColor(7);
  h_trk_pt_primary_z0Cuts->SetStats(0);
  h_trk_pt_primary_z0Cuts->Draw("SAME");
  TH1F *h_trk_pt_primary_allCutsEff = (TH1F*)h_trk_pt_primary_allCuts->Clone(); 
  h_trk_pt_primary_allCutsEff->Divide(h_trk_pt_primary_allCutsEff,h_trk_pt_primary_noCuts,1.0,1.0,"B");
  h_trk_pt_primary_allCutsEff->SetLineColor(28);
  h_trk_pt_primary_allCutsEff->SetMarkerColor(28);
  h_trk_pt_primary_allCutsEff->SetStats(0);
  h_trk_pt_primary_allCutsEff->Draw("SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_pt_primary_chi2rzdofCuts,"#chi^{2}_{rz}/d.o.f Cut","lp");
  l->AddEntry(h_trk_pt_primary_bendchi2Cuts,"#chi^{2}_{bend} Cut","lp");
  l->AddEntry(h_trk_pt_primary_chi2rphidofCuts,"#chi^{2}_{r#phi}/d.o.f Cut","lp");
  l->AddEntry(h_trk_pt_primary_nstubCuts,"n_{stub} Cut","lp");
  l->AddEntry(h_trk_pt_primary_ptCuts,"p_{T} Cut","lp");
  l->AddEntry(h_trk_pt_primary_d0Cuts,"d_{0} Cut","lp");
  l->AddEntry(h_trk_pt_primary_z0Cuts,"z_{0} Cut","lp");
  l->AddEntry(h_trk_pt_primary_allCutsEff,"All Cuts","lp");
  l->Draw();
  c.SaveAs(DIR + "/h_trkEffOverlay_pt_primary.pdf");
  delete h_trk_pt_primary;
  delete h_trk_pt_primary_noCuts;
  delete h_trk_pt_primary_ptCuts;
  delete h_trk_pt_primary_d0Cuts;
  delete h_trk_pt_primary_chi2rzdofCuts;
  delete h_trk_pt_primary_bendchi2Cuts;
  delete h_trk_pt_primary_chi2rphidofCuts;
  delete h_trk_pt_primary_nstubCuts;
  delete h_trk_pt_primary_z0Cuts;
  delete h_trk_pt_primary_allCutsEff;

  int numPart = numPart_primary_noCuts.size();
  TH1F *h_numPart_primary_noCuts = new TH1F("h_numPart_primary_noCuts","h_numPart_primary_noCuts; pdgid; Number of Particles",numPart,0,numPart);
  TH1F *h_numPart_primary_chi2rzdofCuts = new TH1F("h_numPart_primary_chi2rzdofCuts","h_numPart_primary_chi2rzdofCuts; pdgid; Number of Particles",numPart,0,numPart);
  TH1F *h_numPart_primary_bendchi2Cuts = new TH1F("h_numPart_primary_bendchi2Cuts","h_numPart_primary_bendchi2Cuts; pdgid; Number of Particles",numPart,0,numPart);
  TH1F *h_numPart_primary_chi2rphidofCuts = new TH1F("h_numPart_primary_chi2rphidofCuts","h_numPart_primary_chi2rphidofCuts; pdgid; Number of Particles",numPart,0,numPart);
  TH1F *h_numPart_primary_nstubCuts = new TH1F("h_numPart_primary_nstubCuts","h_numPart_primary_nstubCuts; pdgid; Number of Particles",numPart,0,numPart);
  TH1F *h_numPart_primary_ptCuts = new TH1F("h_numPart_primary_ptCuts","h_numPart_primary_ptCuts; pdgid; Number of Particles",numPart,0,numPart);
  TH1F *h_numPart_primary_d0Cuts = new TH1F("h_numPart_primary_d0Cuts","h_numPart_primary_d0Cuts; pdgid; Number of Particles",numPart,0,numPart);
  TH1F *h_numPart_primary_z0Cuts = new TH1F("h_numPart_primary_z0Cuts","h_numPart_primary_z0Cuts; pdgid; Number of Particles",numPart,0,numPart);

  int binNum = 1;
  for(const auto & [key, value] : numPart_primary_noCuts){
    h_numPart_primary_noCuts->SetBinContent(binNum,value);
    h_numPart_primary_noCuts->GetXaxis()->SetBinLabel(binNum,key.c_str());
    binNum++;
  }
  h_numPart_primary_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_numPart_primary_noCuts);
  binNum = 1;
  for(const auto & [key, value] : numPart_primary_ptCuts){
    h_numPart_primary_ptCuts->SetBinContent(binNum,value);
    h_numPart_primary_ptCuts->GetXaxis()->SetBinLabel(binNum,key.c_str());
    binNum++;
  }
  h_numPart_primary_ptCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_numPart_primary_ptCuts);
  binNum = 1;
  for(const auto & [key, value] : numPart_primary_d0Cuts){
    h_numPart_primary_d0Cuts->SetBinContent(binNum,value);
    h_numPart_primary_d0Cuts->GetXaxis()->SetBinLabel(binNum,key.c_str());
    binNum++;
  }
  h_numPart_primary_d0Cuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_numPart_primary_d0Cuts);
  binNum = 1;
  for(const auto & [key, value] : numPart_primary_chi2rzdofCuts){
    h_numPart_primary_chi2rzdofCuts->SetBinContent(binNum,value);
    h_numPart_primary_chi2rzdofCuts->GetXaxis()->SetBinLabel(binNum,key.c_str());
    binNum++;
  }
  h_numPart_primary_chi2rzdofCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_numPart_primary_chi2rzdofCuts);
  binNum = 1;
  for(const auto & [key, value] : numPart_primary_bendchi2Cuts){
    h_numPart_primary_bendchi2Cuts->SetBinContent(binNum,value);
    h_numPart_primary_bendchi2Cuts->GetXaxis()->SetBinLabel(binNum,key.c_str());
    binNum++;
  }
  h_numPart_primary_bendchi2Cuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_numPart_primary_bendchi2Cuts);
  binNum = 1;
  for(const auto & [key, value] : numPart_primary_chi2rphidofCuts){
    h_numPart_primary_chi2rphidofCuts->SetBinContent(binNum,value);
    h_numPart_primary_chi2rphidofCuts->GetXaxis()->SetBinLabel(binNum,key.c_str());
    binNum++;
  }
  h_numPart_primary_chi2rphidofCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_numPart_primary_chi2rphidofCuts);
  binNum = 1;
  for(const auto & [key, value] : numPart_primary_nstubCuts){
    h_numPart_primary_nstubCuts->SetBinContent(binNum,value);
    h_numPart_primary_nstubCuts->GetXaxis()->SetBinLabel(binNum,key.c_str());
    binNum++;
  }
  h_numPart_primary_nstubCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_numPart_primary_nstubCuts);
  binNum = 1;
  for(const auto & [key, value] : numPart_primary_z0Cuts){
    h_numPart_primary_z0Cuts->SetBinContent(binNum,value);
    h_numPart_primary_z0Cuts->GetXaxis()->SetBinLabel(binNum,key.c_str());
    binNum++;
  }
  h_numPart_primary_z0Cuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_numPart_primary_z0Cuts);
  h_numPart_primary_chi2rzdofCuts->SetName("partEff_pt_primary");
  h_numPart_primary_chi2rzdofCuts->GetXaxis()->SetTitle("pdgid");
  h_numPart_primary_chi2rzdofCuts->GetYaxis()->SetTitle("Cut Efficiency");
  h_numPart_primary_chi2rzdofCuts->Divide(h_numPart_primary_chi2rzdofCuts,h_numPart_primary_noCuts,1.0,1.0,"B");
  h_numPart_primary_chi2rzdofCuts->SetLineColor(1);
  h_numPart_primary_chi2rzdofCuts->SetMarkerColor(1);
  h_numPart_primary_chi2rzdofCuts->SetStats(0);
  h_numPart_primary_chi2rzdofCuts->GetYaxis()->SetRangeUser(0,1);
  h_numPart_primary_chi2rzdofCuts->GetXaxis()->SetRangeUser(0,numPart);
  h_numPart_primary_chi2rzdofCuts->Draw();
  h_numPart_primary_bendchi2Cuts->Divide(h_numPart_primary_bendchi2Cuts,h_numPart_primary_noCuts,1.0,1.0,"B");
  h_numPart_primary_bendchi2Cuts->SetLineColor(2);
  h_numPart_primary_bendchi2Cuts->SetMarkerColor(2);
  h_numPart_primary_bendchi2Cuts->SetStats(0);
  h_numPart_primary_bendchi2Cuts->Draw("SAME");
  h_numPart_primary_chi2rphidofCuts->Divide(h_numPart_primary_chi2rphidofCuts,h_numPart_primary_noCuts,1.0,1.0,"B");
  h_numPart_primary_chi2rphidofCuts->SetLineColor(3);
  h_numPart_primary_chi2rphidofCuts->SetMarkerColor(3);
  h_numPart_primary_chi2rphidofCuts->SetStats(0);
  h_numPart_primary_chi2rphidofCuts->Draw("SAME");
  h_numPart_primary_nstubCuts->Divide(h_numPart_primary_nstubCuts,h_numPart_primary_noCuts,1.0,1.0,"B");
  h_numPart_primary_nstubCuts->SetLineColor(4);
  h_numPart_primary_nstubCuts->SetMarkerColor(4);
  h_numPart_primary_nstubCuts->SetStats(0);
  h_numPart_primary_nstubCuts->Draw("SAME");
  h_numPart_primary_ptCuts->Divide(h_numPart_primary_ptCuts,h_numPart_primary_noCuts,1.0,1.0,"B");
  h_numPart_primary_ptCuts->SetLineColor(5);
  h_numPart_primary_ptCuts->SetMarkerColor(5);
  h_numPart_primary_ptCuts->SetStats(0);
  h_numPart_primary_ptCuts->Draw("SAME");
  h_numPart_primary_d0Cuts->Divide(h_numPart_primary_d0Cuts,h_numPart_primary_noCuts,1.0,1.0,"B");
  h_numPart_primary_d0Cuts->SetLineColor(6);
  h_numPart_primary_d0Cuts->SetMarkerColor(6);
  h_numPart_primary_d0Cuts->SetStats(0);
  h_numPart_primary_d0Cuts->Draw("SAME");
  h_numPart_primary_z0Cuts->Divide(h_numPart_primary_z0Cuts,h_numPart_primary_noCuts,1.0,1.0,"B");
  h_numPart_primary_z0Cuts->SetLineColor(7);
  h_numPart_primary_z0Cuts->SetMarkerColor(7);
  h_numPart_primary_z0Cuts->SetStats(0);
  h_numPart_primary_z0Cuts->Draw("SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_numPart_primary_chi2rzdofCuts,"#chi^{2}_{rz}/d.o.f Cut","lp");
  l->AddEntry(h_numPart_primary_bendchi2Cuts,"#chi^{2}_{bend} Cut","lp");
  l->AddEntry(h_numPart_primary_chi2rphidofCuts,"#chi^{2}_{r#phi}/d.o.f Cut","lp");
  l->AddEntry(h_numPart_primary_nstubCuts,"n_{stub} Cut","lp");
  l->AddEntry(h_numPart_primary_ptCuts,"p_{T} Cut","lp");
  l->AddEntry(h_numPart_primary_d0Cuts,"d_{0} Cut","lp");
  l->AddEntry(h_numPart_primary_z0Cuts,"z_{0} Cut","lp");
  l->Draw();
  c.SaveAs(DIR + "/h_partEffOverlay_pt_primary.pdf");
  delete h_numPart_primary_noCuts;
  delete h_numPart_primary_ptCuts;
  delete h_numPart_primary_d0Cuts;
  delete h_numPart_primary_chi2rzdofCuts;
  delete h_numPart_primary_bendchi2Cuts;
  delete h_numPart_primary_chi2rphidofCuts;
  delete h_numPart_primary_nstubCuts;
  delete h_numPart_primary_z0Cuts;
  
  h_trk_pt_np_ptCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_np_ptCuts);
  h_trk_pt_np_ptCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_pt_np_ptCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_pt_np_ptCuts->GetName() + ".pdf");

  h_trk_pt_np_d0Cuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_np_d0Cuts);
  h_trk_pt_np_d0Cuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_pt_np_d0Cuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_pt_np_d0Cuts->GetName() + ".pdf");

  h_trk_pt_np_chi2rzdofCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_np_chi2rzdofCuts);
  h_trk_pt_np_chi2rzdofCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_pt_np_chi2rzdofCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_pt_np_chi2rzdofCuts->GetName() + ".pdf");
  
  h_trk_pt_np_bendchi2Cuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_np_bendchi2Cuts);
  h_trk_pt_np_bendchi2Cuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_pt_np_bendchi2Cuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_pt_np_bendchi2Cuts->GetName() + ".pdf");
  
  h_trk_pt_np_z0Cuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_np_z0Cuts);
  h_trk_pt_np_z0Cuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_pt_np_z0Cuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_pt_np_z0Cuts->GetName() + ".pdf");

  h_trk_pt_np_chi2rphidofCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_np_chi2rphidofCuts);
  h_trk_pt_np_chi2rphidofCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_pt_np_chi2rphidofCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_pt_np_chi2rphidofCuts->GetName() + ".pdf");

  h_trk_pt_np_nstubCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_np_nstubCuts);
  h_trk_pt_np_chi2rzdofCuts->SetName("trkEff_pt_np");
  h_trk_pt_np_chi2rzdofCuts->GetXaxis()->SetTitle("Track p_{T} (GeV)");
  h_trk_pt_np_chi2rzdofCuts->GetYaxis()->SetTitle("Cut Efficiency");
  h_trk_pt_np_chi2rzdofCuts->Divide(h_trk_pt_np_chi2rzdofCuts,h_trk_pt_np_noCuts,1.0,1.0,"B");
  h_trk_pt_np_chi2rzdofCuts->SetLineColor(1);
  h_trk_pt_np_chi2rzdofCuts->SetMarkerColor(1);
  h_trk_pt_np_chi2rzdofCuts->SetStats(0);
  h_trk_pt_np_chi2rzdofCuts->Draw();
  h_trk_pt_np_bendchi2Cuts->Divide(h_trk_pt_np_bendchi2Cuts,h_trk_pt_np_noCuts,1.0,1.0,"B");
  h_trk_pt_np_bendchi2Cuts->SetLineColor(2);
  h_trk_pt_np_bendchi2Cuts->SetMarkerColor(2);
  h_trk_pt_np_bendchi2Cuts->SetStats(0);
  h_trk_pt_np_bendchi2Cuts->Draw("SAME");
  h_trk_pt_np_chi2rphidofCuts->Divide(h_trk_pt_np_chi2rphidofCuts,h_trk_pt_np_noCuts,1.0,1.0,"B");
  h_trk_pt_np_chi2rphidofCuts->SetLineColor(3);
  h_trk_pt_np_chi2rphidofCuts->SetMarkerColor(3);
  h_trk_pt_np_chi2rphidofCuts->SetStats(0);
  h_trk_pt_np_chi2rphidofCuts->Draw("SAME");
  h_trk_pt_np_nstubCuts->Divide(h_trk_pt_np_nstubCuts,h_trk_pt_np_noCuts,1.0,1.0,"B");
  h_trk_pt_np_nstubCuts->SetLineColor(4);
  h_trk_pt_np_nstubCuts->SetMarkerColor(4);
  h_trk_pt_np_nstubCuts->SetStats(0);
  h_trk_pt_np_nstubCuts->Draw("SAME");
  h_trk_pt_np_ptCuts->Divide(h_trk_pt_np_ptCuts,h_trk_pt_np_noCuts,1.0,1.0,"B");
  h_trk_pt_np_ptCuts->SetLineColor(5);
  h_trk_pt_np_ptCuts->SetMarkerColor(5);
  h_trk_pt_np_ptCuts->SetStats(0);
  h_trk_pt_np_ptCuts->Draw("SAME");
  h_trk_pt_np_d0Cuts->Divide(h_trk_pt_np_d0Cuts,h_trk_pt_np_noCuts,1.0,1.0,"B");
  h_trk_pt_np_d0Cuts->SetLineColor(6);
  h_trk_pt_np_d0Cuts->SetMarkerColor(6);
  h_trk_pt_np_d0Cuts->SetStats(0);
  h_trk_pt_np_d0Cuts->Draw("SAME");
  h_trk_pt_np_z0Cuts->Divide(h_trk_pt_np_z0Cuts,h_trk_pt_np_noCuts,1.0,1.0,"B");
  h_trk_pt_np_z0Cuts->SetLineColor(7);
  h_trk_pt_np_z0Cuts->SetMarkerColor(7);
  h_trk_pt_np_z0Cuts->SetStats(0);
  h_trk_pt_np_z0Cuts->Draw("SAME");
  TH1F *h_trk_pt_np_allCutsEff = (TH1F*)h_trk_pt_np_allCuts->Clone();
  h_trk_pt_np_allCutsEff->Divide(h_trk_pt_np_allCutsEff,h_trk_pt_np_noCuts,1.0,1.0,"B");
  h_trk_pt_np_allCutsEff->SetLineColor(28);
  h_trk_pt_np_allCutsEff->SetMarkerColor(28);
  h_trk_pt_np_allCutsEff->SetStats(0);
  h_trk_pt_np_allCutsEff->Draw("SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_pt_np_chi2rzdofCuts,"#chi^{2}_{rz}/d.o.f Cut","lp");
  l->AddEntry(h_trk_pt_np_bendchi2Cuts,"#chi^{2}_{bend} Cut","lp");
  l->AddEntry(h_trk_pt_np_chi2rphidofCuts,"#chi^{2}_{r#phi}/d.o.f Cut","lp");
  l->AddEntry(h_trk_pt_np_nstubCuts,"n_{stub} Cut","lp");
  l->AddEntry(h_trk_pt_np_ptCuts,"p_{T} Cut","lp");
  l->AddEntry(h_trk_pt_np_d0Cuts,"d_{0} Cut","lp");
  l->AddEntry(h_trk_pt_np_z0Cuts,"z_{0} Cut","lp");
  l->AddEntry(h_trk_pt_np_allCutsEff,"All Cuts","lp");
  l->Draw();
  c.SaveAs(DIR + "/h_trkEffOverlay_pt_np.pdf");
  delete h_trk_pt_np;
  delete h_trk_pt_np_noCuts;
  delete h_trk_pt_np_ptCuts;
  delete h_trk_pt_np_d0Cuts;
  delete h_trk_pt_np_chi2rzdofCuts;
  delete h_trk_pt_np_bendchi2Cuts;
  delete h_trk_pt_np_chi2rphidofCuts;
  delete h_trk_pt_np_nstubCuts;
  delete h_trk_pt_np_z0Cuts;
  delete h_trk_pt_np_allCutsEff;
  
  
  numPart = numPart_np_noCuts.size();
  TH1F *h_numPart_np_noCuts = new TH1F("h_numPart_np_noCuts","h_numPart_np_noCuts; pdgid; Number of Particles",numPart,0,numPart);
  TH1F *h_numPart_np_chi2rzdofCuts = new TH1F("h_numPart_np_chi2rzdofCuts","h_numPart_np_chi2rzdofCuts; pdgid; Number of Particles",numPart,0,numPart);
  TH1F *h_numPart_np_bendchi2Cuts = new TH1F("h_numPart_np_bendchi2Cuts","h_numPart_np_bendchi2Cuts; pdgid; Number of Particles",numPart,0,numPart);
  TH1F *h_numPart_np_chi2rphidofCuts = new TH1F("h_numPart_np_chi2rphidofCuts","h_numPart_np_chi2rphidofCuts; pdgid; Number of Particles",numPart,0,numPart);
  TH1F *h_numPart_np_nstubCuts = new TH1F("h_numPart_np_nstubCuts","h_numPart_np_nstubCuts; pdgid; Number of Particles",numPart,0,numPart);
  TH1F *h_numPart_np_ptCuts = new TH1F("h_numPart_np_ptCuts","h_numPart_np_ptCuts; pdgid; Number of Particles",numPart,0,numPart);
  TH1F *h_numPart_np_d0Cuts = new TH1F("h_numPart_np_d0Cuts","h_numPart_np_d0Cuts; pdgid; Number of Particles",numPart,0,numPart);
  TH1F *h_numPart_np_z0Cuts = new TH1F("h_numPart_np_z0Cuts","h_numPart_np_z0Cuts; pdgid; Number of Particles",numPart,0,numPart);

  binNum = 1;
  for(const auto & [key, value] : numPart_np_noCuts){
    h_numPart_np_noCuts->SetBinContent(binNum,value);
    h_numPart_np_noCuts->GetXaxis()->SetBinLabel(binNum,key.c_str());
    binNum++;
  }
  h_numPart_np_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_numPart_np_noCuts);
  binNum = 1;
  for(const auto & [key, value] : numPart_np_ptCuts){
    h_numPart_np_ptCuts->SetBinContent(binNum,value);
    h_numPart_np_ptCuts->GetXaxis()->SetBinLabel(binNum,key.c_str());
    binNum++;
  }
  h_numPart_np_ptCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_numPart_np_ptCuts);
  binNum = 1;
  for(const auto & [key, value] : numPart_np_d0Cuts){
    h_numPart_np_d0Cuts->SetBinContent(binNum,value);
    h_numPart_np_d0Cuts->GetXaxis()->SetBinLabel(binNum,key.c_str());
    binNum++;
  }
  h_numPart_np_d0Cuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_numPart_np_d0Cuts);
  binNum = 1;
  for(const auto & [key, value] : numPart_np_chi2rzdofCuts){
    h_numPart_np_chi2rzdofCuts->SetBinContent(binNum,value);
    h_numPart_np_chi2rzdofCuts->GetXaxis()->SetBinLabel(binNum,key.c_str());
    binNum++;
  }
  h_numPart_np_chi2rzdofCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_numPart_np_chi2rzdofCuts);
  binNum = 1;
  for(const auto & [key, value] : numPart_np_bendchi2Cuts){
    h_numPart_np_bendchi2Cuts->SetBinContent(binNum,value);
    h_numPart_np_bendchi2Cuts->GetXaxis()->SetBinLabel(binNum,key.c_str());
    binNum++;
  }
  h_numPart_np_bendchi2Cuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_numPart_np_bendchi2Cuts);
  binNum = 1;
  for(const auto & [key, value] : numPart_np_chi2rphidofCuts){
    h_numPart_np_chi2rphidofCuts->SetBinContent(binNum,value);
    h_numPart_np_chi2rphidofCuts->GetXaxis()->SetBinLabel(binNum,key.c_str());
    binNum++;
  }
  h_numPart_np_chi2rphidofCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_numPart_np_chi2rphidofCuts);
  binNum = 1;
  for(const auto & [key, value] : numPart_np_nstubCuts){
    h_numPart_np_nstubCuts->SetBinContent(binNum,value);
    h_numPart_np_nstubCuts->GetXaxis()->SetBinLabel(binNum,key.c_str());
    binNum++;
  }
  h_numPart_np_nstubCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_numPart_np_nstubCuts);
  binNum = 1;
  for(const auto & [key, value] : numPart_np_z0Cuts){
    h_numPart_np_z0Cuts->SetBinContent(binNum,value);
    h_numPart_np_z0Cuts->GetXaxis()->SetBinLabel(binNum,key.c_str());
    binNum++;
  }
  h_numPart_np_z0Cuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_numPart_np_z0Cuts);
  h_numPart_np_chi2rzdofCuts->SetName("partEff_pt_np");
  h_numPart_np_chi2rzdofCuts->GetXaxis()->SetTitle("pdgid");
  h_numPart_np_chi2rzdofCuts->GetYaxis()->SetTitle("Cut Efficiency");
  h_numPart_np_chi2rzdofCuts->Divide(h_numPart_np_chi2rzdofCuts,h_numPart_np_noCuts,1.0,1.0,"B");
  h_numPart_np_chi2rzdofCuts->SetLineColor(1);
  h_numPart_np_chi2rzdofCuts->SetMarkerColor(1);
  h_numPart_np_chi2rzdofCuts->SetStats(0);
  h_numPart_np_chi2rzdofCuts->GetYaxis()->SetRangeUser(0,1);
  h_numPart_np_chi2rzdofCuts->GetXaxis()->SetRangeUser(0,numPart);
  h_numPart_np_chi2rzdofCuts->Draw();
  h_numPart_np_bendchi2Cuts->Divide(h_numPart_np_bendchi2Cuts,h_numPart_np_noCuts,1.0,1.0,"B");
  h_numPart_np_bendchi2Cuts->SetLineColor(2);
  h_numPart_np_bendchi2Cuts->SetMarkerColor(2);
  h_numPart_np_bendchi2Cuts->SetStats(0);
  h_numPart_np_bendchi2Cuts->Draw("SAME");
  h_numPart_np_chi2rphidofCuts->Divide(h_numPart_np_chi2rphidofCuts,h_numPart_np_noCuts,1.0,1.0,"B");
  h_numPart_np_chi2rphidofCuts->SetLineColor(3);
  h_numPart_np_chi2rphidofCuts->SetMarkerColor(3);
  h_numPart_np_chi2rphidofCuts->SetStats(0);
  h_numPart_np_chi2rphidofCuts->Draw("SAME");
  h_numPart_np_nstubCuts->Divide(h_numPart_np_nstubCuts,h_numPart_np_noCuts,1.0,1.0,"B");
  h_numPart_np_nstubCuts->SetLineColor(4);
  h_numPart_np_nstubCuts->SetMarkerColor(4);
  h_numPart_np_nstubCuts->SetStats(0);
  h_numPart_np_nstubCuts->Draw("SAME");
  h_numPart_np_ptCuts->Divide(h_numPart_np_ptCuts,h_numPart_np_noCuts,1.0,1.0,"B");
  h_numPart_np_ptCuts->SetLineColor(5);
  h_numPart_np_ptCuts->SetMarkerColor(5);
  h_numPart_np_ptCuts->SetStats(0);
  h_numPart_np_ptCuts->Draw("SAME");
  h_numPart_np_d0Cuts->Divide(h_numPart_np_d0Cuts,h_numPart_np_noCuts,1.0,1.0,"B");
  h_numPart_np_d0Cuts->SetLineColor(6);
  h_numPart_np_d0Cuts->SetMarkerColor(6);
  h_numPart_np_d0Cuts->SetStats(0);
  h_numPart_np_d0Cuts->Draw("SAME");
  h_numPart_np_z0Cuts->Divide(h_numPart_np_z0Cuts,h_numPart_np_noCuts,1.0,1.0,"B");
  h_numPart_np_z0Cuts->SetLineColor(7);
  h_numPart_np_z0Cuts->SetMarkerColor(7);
  h_numPart_np_z0Cuts->SetStats(0);
  h_numPart_np_z0Cuts->Draw("SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_numPart_np_chi2rzdofCuts,"#chi^{2}_{rz}/d.o.f Cut","lp");
  l->AddEntry(h_numPart_np_bendchi2Cuts,"#chi^{2}_{bend} Cut","lp");
  l->AddEntry(h_numPart_np_chi2rphidofCuts,"#chi^{2}_{r#phi}/d.o.f Cut","lp");
  l->AddEntry(h_numPart_np_nstubCuts,"n_{stub} Cut","lp");
  l->AddEntry(h_numPart_np_ptCuts,"p_{T} Cut","lp");
  l->AddEntry(h_numPart_np_d0Cuts,"d_{0} Cut","lp");
  l->AddEntry(h_numPart_np_z0Cuts,"z_{0} Cut","lp");
  l->Draw();
  c.SaveAs(DIR + "/h_partEffOverlay_pt_np.pdf");
  delete h_numPart_np_noCuts;
  delete h_numPart_np_ptCuts;
  delete h_numPart_np_d0Cuts;
  delete h_numPart_np_chi2rzdofCuts;
  delete h_numPart_np_bendchi2Cuts;
  delete h_numPart_np_chi2rphidofCuts;
  delete h_numPart_np_nstubCuts;
  delete h_numPart_np_z0Cuts;

  h_trk_pt_primary_allCuts->Scale(1./h_trk_pt_primary_allCuts->Integral());
  h_trk_pt_primary_allCuts->SetLineColor(1);
  h_trk_pt_primary_allCuts->SetMarkerColor(1);
  h_trk_pt_np_allCuts->Scale(1./h_trk_pt_np_allCuts->Integral());
  h_trk_pt_np_allCuts->SetLineColor(2);
  h_trk_pt_np_allCuts->SetMarkerColor(2);
  raiseMax(h_trk_pt_primary_allCuts,h_trk_pt_np_allCuts);
  h_trk_pt_np_allCuts->SetStats(0);
  h_trk_pt_primary_allCuts->SetStats(0);
  h_trk_pt_primary_allCuts->Draw("HIST");
  h_trk_pt_np_allCuts->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_pt_primary_allCuts,"Primary","l");
  l->AddEntry(h_trk_pt_np_allCuts,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_pt_allCuts.pdf");
  delete h_trk_pt_primary_allCuts;
  delete h_trk_pt_np_allCuts;

  h_trk_pt_primary_allCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_primary_allCuts_barrel);
  h_trk_pt_np_allCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_np_allCuts_barrel);
  h_trk_pt_primary_allCuts_barrel->Scale(1./h_trk_pt_primary_allCuts_barrel->Integral());
  h_trk_pt_primary_allCuts_barrel->SetLineColor(1);
  h_trk_pt_primary_allCuts_barrel->SetMarkerColor(1);
  h_trk_pt_np_allCuts_barrel->Scale(1./h_trk_pt_np_allCuts_barrel->Integral());
  h_trk_pt_np_allCuts_barrel->SetLineColor(2);
  h_trk_pt_np_allCuts_barrel->SetMarkerColor(2);
  h_trk_pt_np_allCuts_barrel->Draw("HIST");
  h_trk_pt_primary_allCuts_barrel->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_pt_primary_allCuts_barrel,"Primary","l");
  l->AddEntry(h_trk_pt_np_allCuts_barrel,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_pt_allCuts_barrel.pdf");
  delete h_trk_pt_primary_allCuts_barrel;
  delete h_trk_pt_np_allCuts_barrel;

  h_trk_pt_primary_allCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_primary_allCuts_disk);
  h_trk_pt_np_allCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_pt_np_allCuts_disk);
  h_trk_pt_primary_allCuts_disk->Scale(1./h_trk_pt_primary_allCuts_disk->Integral());
  h_trk_pt_primary_allCuts_disk->SetLineColor(1);
  h_trk_pt_primary_allCuts_disk->SetMarkerColor(1);
  h_trk_pt_np_allCuts_disk->Scale(1./h_trk_pt_np_allCuts_disk->Integral());
  h_trk_pt_np_allCuts_disk->SetLineColor(2);
  h_trk_pt_np_allCuts_disk->SetMarkerColor(2);
  h_trk_pt_np_allCuts_disk->Draw("HIST");
  h_trk_pt_primary_allCuts_disk->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_pt_primary_allCuts_disk,"Primary","l");
  l->AddEntry(h_trk_pt_np_allCuts_disk,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_pt_allCuts_disk.pdf");
  delete h_trk_pt_primary_allCuts_disk;
  delete h_trk_pt_np_allCuts_disk;

  h_trk_ptIso4_primary_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_ptIso4_primary_allCuts);
  h_trk_ptIso4_np_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_ptIso4_np_allCuts);
  h_trk_ptIso4_primary_allCuts->Scale(1./h_trk_ptIso4_primary_allCuts->Integral());
  h_trk_ptIso4_primary_allCuts->SetLineColor(1);
  h_trk_ptIso4_primary_allCuts->SetMarkerColor(1);
  h_trk_ptIso4_np_allCuts->Scale(1./h_trk_ptIso4_np_allCuts->Integral());
  h_trk_ptIso4_np_allCuts->SetLineColor(2);
  h_trk_ptIso4_np_allCuts->SetMarkerColor(2);
  raiseMax(h_trk_ptIso4_np_allCuts,h_trk_ptIso4_primary_allCuts);
  h_trk_ptIso4_np_allCuts->Draw("HIST");
  h_trk_ptIso4_primary_allCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_ptIso4_primary_allCuts,"Primary","l");
  l->AddEntry(h_trk_ptIso4_np_allCuts,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_ptIso4_allCuts.pdf");
  delete h_trk_ptIso4_primary_allCuts;
  delete h_trk_ptIso4_np_allCuts;

  h_trk_ptIso8_primary_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_ptIso8_primary_allCuts);
  h_trk_ptIso8_np_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_ptIso8_np_allCuts);
  h_trk_ptIso8_primary_allCuts->Scale(1./h_trk_ptIso8_primary_allCuts->Integral());
  h_trk_ptIso8_primary_allCuts->SetLineColor(1);
  h_trk_ptIso8_primary_allCuts->SetMarkerColor(1);
  h_trk_ptIso8_np_allCuts->Scale(1./h_trk_ptIso8_np_allCuts->Integral());
  h_trk_ptIso8_np_allCuts->SetLineColor(2);
  h_trk_ptIso8_np_allCuts->SetMarkerColor(2);
  raiseMax(h_trk_ptIso8_np_allCuts,h_trk_ptIso8_primary_allCuts);
  h_trk_ptIso8_np_allCuts->Draw("HIST");
  h_trk_ptIso8_primary_allCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_ptIso8_primary_allCuts,"Primary","l");
  l->AddEntry(h_trk_ptIso8_np_allCuts,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_ptIso8_allCuts.pdf");
  delete h_trk_ptIso8_primary_allCuts;
  delete h_trk_ptIso8_np_allCuts;

  h_tp_pt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_tp_pt);
  h_tp_pt->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_tp_pt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_tp_pt->GetName() + ".pdf");
  delete h_tp_pt;

  h_trk_eta->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta);
  h_trk_eta->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_eta->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_eta->GetName() + ".pdf");
  delete h_trk_eta;

  h_trk_eta_primary->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta_primary);
  h_trk_eta_primary->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_eta_primary->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_eta_primary->GetName() + ".pdf");
  delete h_trk_eta_primary;

  h_trk_eta_primary_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta_primary_noCuts);
  h_trk_eta_primary_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_eta_primary_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_eta_primary_noCuts->GetName() + ".pdf");

  h_trk_eta_primary_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta_primary_allCuts);
  h_trk_eta_primary_allCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_eta_primary_allCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_eta_primary_allCuts->GetName() + ".pdf");

  h_trk_eta_np->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta_np);
  h_trk_eta_np->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_eta_np->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_eta_np->GetName() + ".pdf");
  delete h_trk_eta_np;

  h_trk_eta_np_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta_np_noCuts);
  h_trk_eta_np_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_eta_np_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_eta_np_noCuts->GetName() + ".pdf");

  h_trk_eta_np_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta_np_allCuts);
  h_trk_eta_np_allCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_eta_np_allCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_eta_np_allCuts->GetName() + ".pdf");

  h_trk_eta_primary_noCuts->Scale(1./h_trk_eta_primary_noCuts->Integral());
  h_trk_eta_primary_noCuts->SetLineColor(1);
  h_trk_eta_primary_noCuts->SetMarkerColor(1);
  TH1F *h_trk_eta_np_noCutsNorm = (TH1F*)h_trk_eta_np_noCuts->Clone(); 
  h_trk_eta_np_noCutsNorm->Scale(1./h_trk_eta_np_noCutsNorm->Integral());
  h_trk_eta_np_noCutsNorm->SetLineColor(2);
  h_trk_eta_np_noCutsNorm->SetMarkerColor(2);
  raiseMax(h_trk_eta_np_noCutsNorm,h_trk_eta_primary_noCuts);
  h_trk_eta_np_noCutsNorm->SetStats(0);
  h_trk_eta_primary_noCuts->SetStats(0);
  h_trk_eta_np_noCutsNorm->Draw("HIST");
  h_trk_eta_primary_noCuts->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_eta_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_eta_np_noCutsNorm,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_eta_noCuts.pdf");
  delete h_trk_eta_np_noCutsNorm;
  
  TH1F *h_trk_eta_fake_noCutsNorm = (TH1F*)h_trk_eta_fake_noCuts->Clone(); 
  h_trk_eta_fake_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta_fake_noCutsNorm);
  TH1F *h_trk_eta_PU_noCutsNorm = (TH1F*)h_trk_eta_PU_noCuts->Clone(); 
  h_trk_eta_PU_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta_PU_noCutsNorm);
  TH1F *h_trk_eta_notHiggs_noCutsNorm = (TH1F*)h_trk_eta_notHiggs_noCuts->Clone(); 
  h_trk_eta_notHiggs_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta_notHiggs_noCutsNorm);
  h_trk_eta_fake_noCutsNorm->Scale(1./h_trk_eta_np_noCuts->Integral());
  h_trk_eta_fake_noCutsNorm->SetLineColor(2);
  h_trk_eta_fake_noCutsNorm->SetMarkerColor(2);
  h_trk_eta_PU_noCutsNorm->Scale(1./h_trk_eta_np_noCuts->Integral());
  h_trk_eta_PU_noCutsNorm->SetLineColor(3);
  h_trk_eta_PU_noCutsNorm->SetMarkerColor(3);
  h_trk_eta_notHiggs_noCutsNorm->Scale(1./h_trk_eta_np_noCuts->Integral());
  h_trk_eta_notHiggs_noCutsNorm->SetLineColor(4);
  h_trk_eta_notHiggs_noCutsNorm->SetMarkerColor(4);
  auto hs_eta_noCuts = new THStack("hs_eta_noCuts","Stacked BG histograms");
  hs_eta_noCuts->Add(h_trk_eta_fake_noCutsNorm);
  hs_eta_noCuts->Add(h_trk_eta_PU_noCutsNorm);
  hs_eta_noCuts->Add(h_trk_eta_notHiggs_noCutsNorm);
  hs_eta_noCuts->Draw("HIST");
  h_trk_eta_primary_noCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_eta_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_eta_fake_noCutsNorm,"Fake","l");
  l->AddEntry(h_trk_eta_PU_noCutsNorm,"PU","l");
  l->AddEntry(h_trk_eta_notHiggs_noCutsNorm,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_eta_noCutsNorm.pdf");
  delete h_trk_eta_fake_noCutsNorm;
  delete h_trk_eta_PU_noCutsNorm;
  delete h_trk_eta_notHiggs_noCutsNorm;
  delete h_trk_eta_np_noCuts;
  delete hs_eta_noCuts;

  h_trk_eta_fake_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta_fake_noCuts);
  h_trk_eta_PU_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta_PU_noCuts);
  h_trk_eta_notHiggs_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta_notHiggs_noCuts);
  h_trk_eta_fake_noCuts->Scale(1./h_trk_eta_fake_noCuts->Integral());
  h_trk_eta_fake_noCuts->SetLineColor(2);
  h_trk_eta_fake_noCuts->SetMarkerColor(2);
  h_trk_eta_PU_noCuts->Scale(1./h_trk_eta_PU_noCuts->Integral());
  h_trk_eta_PU_noCuts->SetLineColor(3);
  h_trk_eta_PU_noCuts->SetMarkerColor(3);
  h_trk_eta_notHiggs_noCuts->Scale(1./h_trk_eta_notHiggs_noCuts->Integral());
  h_trk_eta_notHiggs_noCuts->SetLineColor(4);
  h_trk_eta_notHiggs_noCuts->SetMarkerColor(4);
  h_trk_eta_fake_noCuts->Draw("HIST");
  h_trk_eta_primary_noCuts->Draw("HIST,SAME");
  h_trk_eta_PU_noCuts->Draw("HIST,SAME");
  h_trk_eta_notHiggs_noCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_eta_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_eta_fake_noCuts,"Fake","l");
  l->AddEntry(h_trk_eta_PU_noCuts,"PU","l");
  l->AddEntry(h_trk_eta_notHiggs_noCuts,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_eta_noCuts.pdf");
  delete h_trk_eta_primary_noCuts;
  delete h_trk_eta_fake_noCuts;
  delete h_trk_eta_PU_noCuts;
  delete h_trk_eta_notHiggs_noCuts;

  h_trk_eta_primary_noCuts_H->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta_primary_noCuts_H);
  h_trk_eta_np_noCuts_H->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta_np_noCuts_H);
  h_trk_eta_primary_noCuts_H->Scale(1./h_trk_eta_primary_noCuts_H->Integral());
  h_trk_eta_primary_noCuts_H->SetLineColor(1);
  h_trk_eta_primary_noCuts_H->SetMarkerColor(1);
  h_trk_eta_np_noCuts_H->Scale(1./h_trk_eta_np_noCuts_H->Integral());
  h_trk_eta_np_noCuts_H->SetLineColor(2);
  h_trk_eta_np_noCuts_H->SetMarkerColor(2);
  h_trk_eta_np_noCuts_H->Draw("HIST");
  h_trk_eta_primary_noCuts_H->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_eta_primary_noCuts_H,"Primary","l");
  l->AddEntry(h_trk_eta_np_noCuts_H,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_eta_noCuts_H.pdf");
  delete h_trk_eta_primary_noCuts_H;
  delete h_trk_eta_np_noCuts_H;

  h_trk_eta_primary_noCuts_L->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta_primary_noCuts_L);
  h_trk_eta_np_noCuts_L->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta_np_noCuts_L);
  h_trk_eta_primary_noCuts_L->Scale(1./h_trk_eta_primary_noCuts_L->Integral());
  h_trk_eta_primary_noCuts_L->SetLineColor(1);
  h_trk_eta_primary_noCuts_L->SetMarkerColor(1);
  h_trk_eta_np_noCuts_L->Scale(1./h_trk_eta_np_noCuts_L->Integral());
  h_trk_eta_np_noCuts_L->SetLineColor(2);
  h_trk_eta_np_noCuts_L->SetMarkerColor(2);
  h_trk_eta_np_noCuts_L->Draw("HIST");
  h_trk_eta_primary_noCuts_L->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_eta_primary_noCuts_L,"Primary","l");
  l->AddEntry(h_trk_eta_np_noCuts_L,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_eta_noCuts_L.pdf");
  delete h_trk_eta_primary_noCuts_L;
  delete h_trk_eta_np_noCuts_L;

  h_trk_eta_primary_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta_primary_qualCuts);
  h_trk_eta_np_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta_np_qualCuts);
  h_trk_eta_primary_qualCuts->Scale(1./h_trk_eta_primary_qualCuts->Integral());
  h_trk_eta_primary_qualCuts->SetLineColor(1);
  h_trk_eta_primary_qualCuts->SetMarkerColor(1);
  TH1F *h_trk_eta_np_qualCutsNorm = (TH1F*)h_trk_eta_np_qualCuts->Clone(); 
  h_trk_eta_np_qualCutsNorm->Scale(1./h_trk_eta_np_qualCutsNorm->Integral());
  h_trk_eta_np_qualCutsNorm->SetLineColor(2);
  h_trk_eta_np_qualCutsNorm->SetMarkerColor(2);
  h_trk_eta_np_qualCutsNorm->Draw("HIST");
  h_trk_eta_primary_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_eta_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_eta_np_qualCutsNorm,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_eta_qualCuts.pdf");
  delete h_trk_eta_np_qualCutsNorm;

  TH1F *h_trk_eta_fake_qualCutsNorm = (TH1F*)h_trk_eta_fake_qualCuts->Clone(); 
  h_trk_eta_fake_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta_fake_qualCutsNorm);
  TH1F *h_trk_eta_PU_qualCutsNorm = (TH1F*)h_trk_eta_PU_qualCuts->Clone(); 
  h_trk_eta_PU_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta_PU_qualCutsNorm);
  TH1F *h_trk_eta_notHiggs_qualCutsNorm = (TH1F*)h_trk_eta_notHiggs_qualCuts->Clone(); 
  h_trk_eta_notHiggs_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta_notHiggs_qualCutsNorm);
  h_trk_eta_fake_qualCutsNorm->Scale(1./h_trk_eta_np_qualCuts->Integral());
  h_trk_eta_fake_qualCutsNorm->SetLineColor(2);
  h_trk_eta_fake_qualCutsNorm->SetMarkerColor(2);
  h_trk_eta_PU_qualCutsNorm->Scale(1./h_trk_eta_np_qualCuts->Integral());
  h_trk_eta_PU_qualCutsNorm->SetLineColor(3);
  h_trk_eta_PU_qualCutsNorm->SetMarkerColor(3);
  h_trk_eta_notHiggs_qualCutsNorm->Scale(1./h_trk_eta_np_qualCuts->Integral());
  h_trk_eta_notHiggs_qualCutsNorm->SetLineColor(4);
  h_trk_eta_notHiggs_qualCutsNorm->SetMarkerColor(4);
  auto hs_eta_qualCuts = new THStack("hs_eta_qualCuts","Stacked BG histograms");
  hs_eta_qualCuts->Add(h_trk_eta_fake_qualCutsNorm);
  hs_eta_qualCuts->Add(h_trk_eta_PU_qualCutsNorm);
  hs_eta_qualCuts->Add(h_trk_eta_notHiggs_qualCutsNorm);
  hs_eta_qualCuts->Draw("HIST");
  h_trk_eta_primary_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_eta_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_eta_fake_qualCutsNorm,"Fake","l");
  l->AddEntry(h_trk_eta_PU_qualCutsNorm,"PU","l");
  l->AddEntry(h_trk_eta_notHiggs_qualCutsNorm,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_eta_qualCutsNorm.pdf");
  delete h_trk_eta_fake_qualCutsNorm;
  delete h_trk_eta_PU_qualCutsNorm;
  delete h_trk_eta_notHiggs_qualCutsNorm;
  delete h_trk_eta_np_qualCuts;
  delete hs_eta_qualCuts;

  h_trk_eta_fake_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta_fake_qualCuts);
  h_trk_eta_PU_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta_PU_qualCuts);
  h_trk_eta_notHiggs_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_eta_notHiggs_qualCuts);
  h_trk_eta_fake_qualCuts->Scale(1./h_trk_eta_fake_qualCuts->Integral());
  h_trk_eta_fake_qualCuts->SetLineColor(2);
  h_trk_eta_fake_qualCuts->SetMarkerColor(2);
  h_trk_eta_PU_qualCuts->Scale(1./h_trk_eta_PU_qualCuts->Integral());
  h_trk_eta_PU_qualCuts->SetLineColor(3);
  h_trk_eta_PU_qualCuts->SetMarkerColor(3);
  h_trk_eta_notHiggs_qualCuts->Scale(1./h_trk_eta_notHiggs_qualCuts->Integral());
  h_trk_eta_notHiggs_qualCuts->SetLineColor(4);
  h_trk_eta_notHiggs_qualCuts->SetMarkerColor(4);
  h_trk_eta_fake_qualCuts->Draw("HIST");
  h_trk_eta_primary_qualCuts->Draw("HIST,SAME");
  h_trk_eta_PU_qualCuts->Draw("HIST,SAME");
  h_trk_eta_notHiggs_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_eta_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_eta_fake_qualCuts,"Fake","l");
  l->AddEntry(h_trk_eta_PU_qualCuts,"PU","l");
  l->AddEntry(h_trk_eta_notHiggs_qualCuts,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_eta_qualCuts.pdf");
  delete h_trk_eta_primary_qualCuts;
  delete h_trk_eta_fake_qualCuts;
  delete h_trk_eta_PU_qualCuts;
  delete h_trk_eta_notHiggs_qualCuts;

  h_trk_eta_primary_allCuts->Scale(1./h_trk_eta_primary_allCuts->Integral());
  h_trk_eta_primary_allCuts->SetLineColor(1);
  h_trk_eta_primary_allCuts->SetMarkerColor(1);
  h_trk_eta_np_allCuts->Scale(1./h_trk_eta_np_allCuts->Integral());
  h_trk_eta_np_allCuts->SetLineColor(2);
  h_trk_eta_np_allCuts->SetMarkerColor(2);
  raiseMax(h_trk_eta_primary_allCuts,h_trk_eta_np_allCuts);
  h_trk_eta_primary_allCuts->SetStats(0);
  h_trk_eta_np_allCuts->SetStats(0);
  h_trk_eta_primary_allCuts->Draw("HIST");
  h_trk_eta_np_allCuts->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_eta_primary_allCuts,"Primary","l");
  l->AddEntry(h_trk_eta_np_allCuts,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_eta_allCuts.pdf");
  delete h_trk_eta_primary_allCuts;
  delete h_trk_eta_np_allCuts;

  h_tp_eta->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_tp_eta);
  h_tp_eta->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_tp_eta->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_tp_eta->GetName() + ".pdf");
  delete h_tp_eta;

  h_trueVertex_numAllCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertex_numAllCuts);
  h_trueVertex_numAllCuts->SetStats(0);
  h_trueVertex_numAllCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertex_numAllCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertex_numAllCuts->GetName() + ".pdf");
  delete h_trueVertex_numAllCuts;
   
  h_trueVertex_numNoCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertex_numNoCuts);
  h_trueVertex_numNoCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertex_numNoCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertex_numNoCuts->GetName() + ".pdf");
  delete h_trueVertex_numNoCuts;

  h_trackVertex_numAllCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trackVertex_numAllCuts);
  c.SetLogy();
  h_trackVertex_numAllCuts->SetStats(0);
  h_trackVertex_numAllCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trackVertex_numAllCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trackVertex_numAllCuts->GetName() + ".pdf");
  delete h_trackVertex_numAllCuts;
  c.SetLogy(0);

  h_trackVertex_numNoCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trackVertex_numNoCuts);
  h_trackVertex_numNoCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trackVertex_numNoCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trackVertex_numNoCuts->GetName() + ".pdf");
  delete h_trackVertex_numNoCuts;

  h_trueVertex_x->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertex_x);
  h_trueVertex_x->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertex_x->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertex_x->GetName() + ".pdf");
  delete h_trueVertex_x;

  h_trueVertex_y->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertex_y);
  h_trueVertex_y->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertex_y->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertex_y->GetName() + ".pdf");
  delete h_trueVertex_y;
   
  h_trueVertex_z->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertex_z);
  h_trueVertex_z->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertex_z->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertex_z->GetName() + ".pdf");
  delete h_trueVertex_z; 

  h_trueVertex_sumPt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertex_sumPt);
  h_trueVertex_sumPt->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertex_sumPt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertex_sumPt->GetName() + ".pdf");
  delete h_trueVertex_sumPt;

  h_trueVertex_lowPt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertex_lowPt);
  h_trueVertex_lowPt->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertex_lowPt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertex_lowPt->GetName() + ".pdf");
  delete h_trueVertex_lowPt;

  h_trueVertex_highPt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertex_highPt);
  h_trueVertex_highPt->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertex_highPt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertex_highPt->GetName() + ".pdf");
  delete h_trueVertex_highPt;

  h_false_trackVertex_delta_dist_z->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trackVertex_delta_dist_z);
  h_false_trackVertex_delta_dist_z->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_false_trackVertex_delta_dist_z->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_false_trackVertex_delta_dist_z->GetName() + ".pdf");

  h_correct_trackVertex_delta_dist_z->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trackVertex_delta_dist_z);
  h_correct_trackVertex_delta_dist_z->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_correct_trackVertex_delta_dist_z->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_correct_trackVertex_delta_dist_z->GetName() + ".pdf");

  h_correct_trackVertex_delta_dist_z0->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trackVertex_delta_dist_z0);
  h_correct_trackVertex_delta_dist_z0->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_correct_trackVertex_delta_dist_z0->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_correct_trackVertex_delta_dist_z0->GetName() + ".pdf");
  delete h_correct_trackVertex_delta_dist_z0;

  h_correct_trackVertex_delta_dist_d0->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trackVertex_delta_dist_d0);
  h_correct_trackVertex_delta_dist_d0->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_correct_trackVertex_delta_dist_d0->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_correct_trackVertex_delta_dist_d0->GetName() + ".pdf");
  delete h_correct_trackVertex_delta_dist_d0;

  h_correct_trackVertex_delta_dist_eta->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trackVertex_delta_dist_eta);
  h_correct_trackVertex_delta_dist_eta->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_correct_trackVertex_delta_dist_eta->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_correct_trackVertex_delta_dist_eta->GetName() + ".pdf");

  h_false_trackVertex_delta_dist_eta->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trackVertex_delta_dist_eta);
  h_false_trackVertex_delta_dist_eta->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_false_trackVertex_delta_dist_eta->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_false_trackVertex_delta_dist_eta->GetName() + ".pdf");

  h_correct_trackVertex_delta_dist_phi->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trackVertex_delta_dist_phi);
  h_correct_trackVertex_delta_dist_phi->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_correct_trackVertex_delta_dist_phi->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_correct_trackVertex_delta_dist_phi->GetName() + ".pdf");
  delete h_correct_trackVertex_delta_dist_phi;

  h_false_trackVertex_delta_dist_phi->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trackVertex_delta_dist_phi);
  h_false_trackVertex_delta_dist_phi->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_false_trackVertex_delta_dist_phi->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_false_trackVertex_delta_dist_phi->GetName() + ".pdf");
  delete h_false_trackVertex_delta_dist_phi;

  h_correct_trackVertex_delta_dist_indexPt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trackVertex_delta_dist_indexPt);
  h_correct_trackVertex_delta_dist_indexPt->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_correct_trackVertex_delta_dist_indexPt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_correct_trackVertex_delta_dist_indexPt->GetName() + ".pdf");
  delete h_correct_trackVertex_delta_dist_indexPt;
   
  h_false_trackVertex_numStubs->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trackVertex_numStubs);
  h_false_trackVertex_numStubs->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_false_trackVertex_numStubs->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_false_trackVertex_numStubs->GetName() + ".pdf");
  delete h_false_trackVertex_numStubs;

  h_false_trackVertex_numStubsSum->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trackVertex_numStubsSum);
  h_false_trackVertex_numStubsSum->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_false_trackVertex_numStubsSum->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_false_trackVertex_numStubsSum->GetName() + ".pdf");

  h_false_trackVertex_chi2rphidofSum->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trackVertex_chi2rphidofSum);
  h_false_trackVertex_chi2rphidofSum->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_false_trackVertex_chi2rphidofSum->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_false_trackVertex_chi2rphidofSum->GetName() + ".pdf");

  h_false_trackVertex_MVA1Sum->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trackVertex_MVA1Sum);
  h_false_trackVertex_MVA1Sum->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_false_trackVertex_MVA1Sum->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_false_trackVertex_MVA1Sum->GetName() + ".pdf");
  
  h_false_trackVertex_MVA2Sum->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trackVertex_MVA2Sum);
  h_false_trackVertex_MVA2Sum->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_false_trackVertex_MVA2Sum->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_false_trackVertex_MVA2Sum->GetName() + ".pdf");

  h_false_trackVertex_chi2rzdofSum->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trackVertex_chi2rzdofSum);
  h_false_trackVertex_chi2rzdofSum->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_false_trackVertex_chi2rzdofSum->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_false_trackVertex_chi2rzdofSum->GetName() + ".pdf");

  h_false_trackVertex_bendchi2Sum->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trackVertex_bendchi2Sum);
  h_false_trackVertex_bendchi2Sum->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_false_trackVertex_bendchi2Sum->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_false_trackVertex_bendchi2Sum->GetName() + ".pdf");

  h_false_trackVertex_score->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trackVertex_score);
  h_false_trackVertex_score->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_false_trackVertex_score->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_false_trackVertex_score->GetName() + ".pdf");
  delete h_false_trackVertex_score;

  h_correct_trackVertex_chi2rphidofSum->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trackVertex_chi2rphidofSum);
  h_correct_trackVertex_chi2rphidofSum->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_correct_trackVertex_chi2rphidofSum->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_correct_trackVertex_chi2rphidofSum->GetName() + ".pdf");

  h_correct_trackVertex_MVA1Sum->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trackVertex_MVA1Sum);
  h_correct_trackVertex_MVA1Sum->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_correct_trackVertex_MVA1Sum->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_correct_trackVertex_MVA1Sum->GetName() + ".pdf");

  h_correct_trackVertex_MVA2Sum->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trackVertex_MVA2Sum);
  h_correct_trackVertex_MVA2Sum->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_correct_trackVertex_MVA2Sum->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_correct_trackVertex_MVA2Sum->GetName() + ".pdf");

  h_correct_trackVertex_chi2rzdofSum->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trackVertex_chi2rzdofSum);
  h_correct_trackVertex_chi2rzdofSum->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_correct_trackVertex_chi2rzdofSum->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_correct_trackVertex_chi2rzdofSum->GetName() + ".pdf");
   
  h_correct_trackVertex_bendchi2Sum->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trackVertex_bendchi2Sum);
  h_correct_trackVertex_bendchi2Sum->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_correct_trackVertex_bendchi2Sum->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_correct_trackVertex_bendchi2Sum->GetName() + ".pdf");

  h_correct_trackVertex_score->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trackVertex_score);
  h_correct_trackVertex_score->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_correct_trackVertex_score->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_correct_trackVertex_score->GetName() + ".pdf");
  delete h_correct_trackVertex_score;

  h_all_trackVertex_numStubs->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_all_trackVertex_numStubs);
  h_all_trackVertex_numStubs->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_all_trackVertex_numStubs->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_all_trackVertex_numStubs->GetName() + ".pdf");
  delete h_all_trackVertex_numStubs;

  h_all_trackVertex_fakeId->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_all_trackVertex_fakeId);
  h_all_trackVertex_fakeId->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_all_trackVertex_fakeId->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_all_trackVertex_fakeId->GetName() + ".pdf");
  delete h_all_trackVertex_fakeId;

  h_false_trackVertex_fakeId->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trackVertex_fakeId);
  h_false_trackVertex_fakeId->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_false_trackVertex_fakeId->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_false_trackVertex_fakeId->GetName() + ".pdf");
  delete h_false_trackVertex_fakeId;

  h_correct_trackVertex_numStubs->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trackVertex_numStubs);
  h_correct_trackVertex_numStubs->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_correct_trackVertex_numStubs->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_correct_trackVertex_numStubs->GetName() + ".pdf");
  delete h_correct_trackVertex_numStubs;

  h_correct_trackVertex_numStubsSum->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trackVertex_numStubsSum);
  h_correct_trackVertex_numStubsSum->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_correct_trackVertex_numStubsSum->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_correct_trackVertex_numStubsSum->GetName() + ".pdf");

  h_false_trackVertex_d0->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trackVertex_d0);
  h_false_trackVertex_d0->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_false_trackVertex_d0->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_false_trackVertex_d0->GetName() + ".pdf");
  delete h_false_trackVertex_d0;

  h_false_trackVertex_d_T->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trackVertex_d_T);
  h_false_trackVertex_d_T->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_false_trackVertex_d_T->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_false_trackVertex_d_T->GetName() + ".pdf");

  h_trackVertex_x->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trackVertex_x);
  h_trackVertex_x->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trackVertex_x->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trackVertex_x->GetName() + ".pdf");
  delete h_trackVertex_x;
   
  h_trackVertex_y->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trackVertex_y);
  h_trackVertex_y->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trackVertex_y->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trackVertex_y->GetName() + ".pdf");
  delete h_trackVertex_y;

  h_trackVertex_z->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trackVertex_z);
  h_trackVertex_z->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trackVertex_z->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trackVertex_z->GetName() + ".pdf");
  delete h_trackVertex_z;

  h_trackVertex_sumPt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trackVertex_sumPt);
  h_trackVertex_sumPt->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trackVertex_sumPt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trackVertex_sumPt->GetName() + ".pdf");
  delete h_trackVertex_sumPt;

  h_trackVertexCuts_x->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trackVertexCuts_x);
  h_trackVertexCuts_x->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trackVertexCuts_x->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trackVertexCuts_x->GetName() + ".pdf");
  delete h_trackVertexCuts_x;

  h_trackVertexCuts_y->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trackVertexCuts_y);
  h_trackVertexCuts_y->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trackVertexCuts_y->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trackVertexCuts_y->GetName() + ".pdf");
  delete h_trackVertexCuts_y;

  h_trackVertexCuts_z->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trackVertexCuts_z);
  h_trackVertexCuts_z->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trackVertexCuts_z->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trackVertexCuts_z->GetName() + ".pdf");
  delete h_trackVertexCuts_z;

  h_trackVertexCuts_sumPt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trackVertexCuts_sumPt);
  h_trackVertexCuts_sumPt->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trackVertexCuts_sumPt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trackVertexCuts_sumPt->GetName() + ".pdf");
  delete h_trackVertexCuts_sumPt;

  h_trackVertexCuts_indexPt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trackVertexCuts_indexPt);
  h_trackVertexCuts_indexPt->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trackVertexCuts_indexPt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trackVertexCuts_indexPt->GetName() + ".pdf");
  delete h_trackVertexCuts_indexPt;

  h_correct_trackVertex_indexPt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trackVertex_indexPt);
  h_correct_trackVertex_indexPt->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_correct_trackVertex_indexPt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_correct_trackVertex_indexPt->GetName() + ".pdf");
  delete h_correct_trackVertex_indexPt;

  h_trueVertexCuts_x->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertexCuts_x);
  h_trueVertexCuts_x->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertexCuts_x->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertexCuts_x->GetName() + ".pdf");
  delete h_trueVertexCuts_x;
   
  h_trueVertexCuts_y->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertexCuts_y);
  h_trueVertexCuts_y->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertexCuts_y->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertexCuts_y->GetName() + ".pdf");
  delete h_trueVertexCuts_y;

  h_trueVertexCuts_z->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertexCuts_z);
  h_trueVertexCuts_z->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertexCuts_z->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertexCuts_z->GetName() + ".pdf");
  delete h_trueVertexCuts_z;

  h_trueVertexCuts_sumPt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertexCuts_sumPt);
  h_trueVertexCuts_sumPt->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertexCuts_sumPt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertexCuts_sumPt->GetName() + ".pdf");
  delete h_trueVertexCuts_sumPt;

  h_trueVertexCuts_indexPt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertexCuts_indexPt);
  h_trueVertexCuts_indexPt->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertexCuts_indexPt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertexCuts_indexPt->GetName() + ".pdf");
  delete h_trueVertexCuts_indexPt;

  h_trk_chi2rphidof->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof);
  h_trk_chi2rphidof->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_chi2rphidof->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_chi2rphidof->GetName() + ".pdf");
  delete h_trk_chi2rphidof;

  h_trk_chi2rphidof_primary->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_primary);
  h_trk_chi2rphidof_primary->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_chi2rphidof_primary->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_chi2rphidof_primary->GetName() + ".pdf");
  delete h_trk_chi2rphidof_primary;

  h_trk_chi2rphidof_primary_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_primary_noCuts);
  h_trk_chi2rphidof_primary_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_chi2rphidof_primary_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_chi2rphidof_primary_noCuts->GetName() + ".pdf");

  h_trk_MVA1_primary_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_MVA1_primary_noCuts);
  h_trk_MVA1_primary_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_MVA1_primary_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_MVA1_primary_noCuts->GetName() + ".pdf");

  h_trk_chi2rphidof_np->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_np);
  h_trk_chi2rphidof_np->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_chi2rphidof_np->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_chi2rphidof_np->GetName() + ".pdf");
  delete h_trk_chi2rphidof_np;

  h_trk_chi2rphidof_np_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_np_noCuts);
  h_trk_chi2rphidof_np_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_chi2rphidof_np_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_chi2rphidof_np_noCuts->GetName() + ".pdf");

  h_trk_MVA1_np_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_MVA1_np_noCuts);
  h_trk_MVA1_np_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_MVA1_np_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_MVA1_np_noCuts->GetName() + ".pdf");

  h_correct_trackVertex_d_T->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trackVertex_d_T);
  h_correct_trackVertex_d_T->Scale(1./h_correct_trackVertex_d_T->Integral());
  h_correct_trackVertex_d_T->SetLineColor(1);
  h_correct_trackVertex_d_T->SetMarkerColor(1);
  h_correct_trackVertex_d_T->SetStats(0);
  h_false_trackVertex_d_T->Scale(1./h_false_trackVertex_d_T->Integral());
  h_false_trackVertex_d_T->SetLineColor(2);
  h_false_trackVertex_d_T->SetMarkerColor(2);
  h_false_trackVertex_d_T->SetStats(0);
  raiseMax(h_false_trackVertex_d_T,h_correct_trackVertex_d_T);
  h_false_trackVertex_d_T->Draw("HIST");
  h_correct_trackVertex_d_T->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_correct_trackVertex_d_T,"Correct","l");
  l->AddEntry(h_false_trackVertex_d_T,"False","l");
  l->Draw();
  c.SaveAs(DIR + "/h_correctVsFalse_d_T.pdf");
  delete h_correct_trackVertex_d_T;
  delete h_false_trackVertex_d_T;

  h_correct_trackVertex_deltaPos->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trackVertex_deltaPos);
  h_correct_trackVertex_deltaPos->Scale(1./h_correct_trackVertex_deltaPos->Integral());
  h_correct_trackVertex_deltaPos->SetLineColor(1);
  h_correct_trackVertex_deltaPos->SetMarkerColor(1);
  h_correct_trackVertex_deltaPos->SetStats(0);
  h_false_trackVertex_deltaPos->Scale(1./h_false_trackVertex_deltaPos->Integral());
  h_false_trackVertex_deltaPos->SetLineColor(2);
  h_false_trackVertex_deltaPos->SetMarkerColor(2);
  h_false_trackVertex_deltaPos->SetStats(0);
  raiseMax(h_false_trackVertex_deltaPos,h_correct_trackVertex_deltaPos);
  h_false_trackVertex_deltaPos->Draw("HIST");
  h_correct_trackVertex_deltaPos->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->AddEntry(h_correct_trackVertex_deltaPos,"Correct","l");
  l->AddEntry(h_false_trackVertex_deltaPos,"False","l");
  l->Draw();
  c.SaveAs(DIR + "/h_correctVsFalse_deltaPos.pdf");
  delete h_correct_trackVertex_deltaPos;
  delete h_false_trackVertex_deltaPos;

  h_correct_trackVertex_minD0->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trackVertex_minD0);
  h_false_trackVertex_minD0->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trackVertex_minD0);
  h_correct_trackVertex_minD0->Scale(1./h_correct_trackVertex_minD0->Integral());
  h_correct_trackVertex_minD0->SetLineColor(1);
  h_correct_trackVertex_minD0->SetMarkerColor(1);
  h_correct_trackVertex_minD0->SetStats(0);
  h_false_trackVertex_minD0->Scale(1./h_false_trackVertex_minD0->Integral());
  h_false_trackVertex_minD0->SetLineColor(2);
  h_false_trackVertex_minD0->SetMarkerColor(2);
  h_false_trackVertex_minD0->SetStats(0);
  raiseMax(h_false_trackVertex_minD0,h_correct_trackVertex_minD0);
  h_false_trackVertex_minD0->Draw("HIST");
  h_correct_trackVertex_minD0->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->AddEntry(h_correct_trackVertex_minD0,"Correct","l");
  l->AddEntry(h_false_trackVertex_minD0,"False","l");
  l->Draw();
  c.SaveAs(DIR + "/h_correctVsFalse_minD0.pdf");
  delete h_correct_trackVertex_minD0;
  delete h_false_trackVertex_minD0;

  h_correct_trackVertex_maxD0->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trackVertex_maxD0);
  h_false_trackVertex_maxD0->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trackVertex_maxD0);
  h_correct_trackVertex_maxD0->Scale(1./h_correct_trackVertex_maxD0->Integral());
  h_correct_trackVertex_maxD0->SetLineColor(1);
  h_correct_trackVertex_maxD0->SetMarkerColor(1);
  h_correct_trackVertex_maxD0->SetStats(0);
  h_false_trackVertex_maxD0->Scale(1./h_false_trackVertex_maxD0->Integral());
  h_false_trackVertex_maxD0->SetLineColor(2);
  h_false_trackVertex_maxD0->SetMarkerColor(2);
  h_false_trackVertex_maxD0->SetStats(0);
  raiseMax(h_false_trackVertex_maxD0,h_correct_trackVertex_maxD0);
  h_false_trackVertex_maxD0->Draw("HIST");
  h_correct_trackVertex_maxD0->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->AddEntry(h_correct_trackVertex_maxD0,"Correct","l");
  l->AddEntry(h_false_trackVertex_maxD0,"False","l");
  l->Draw();
  c.SaveAs(DIR + "/h_correctVsFalse_maxD0.pdf");
  delete h_correct_trackVertex_maxD0;
  delete h_false_trackVertex_maxD0;

  h_correct_trackVertex_pt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trackVertex_pt);
  h_false_trackVertex_maxPt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trackVertex_maxPt);
  h_correct_trackVertex_pt->Scale(1./h_correct_trackVertex_pt->Integral());
  h_correct_trackVertex_pt->SetLineColor(1);
  h_correct_trackVertex_pt->SetMarkerColor(1);
  h_correct_trackVertex_pt->SetStats(0);
  h_false_trackVertex_maxPt->Scale(1./h_false_trackVertex_maxPt->Integral());
  h_false_trackVertex_maxPt->SetLineColor(2);
  h_false_trackVertex_maxPt->SetMarkerColor(2);
  h_false_trackVertex_maxPt->SetStats(0);
  raiseMax(h_false_trackVertex_maxPt,h_correct_trackVertex_pt);
  h_false_trackVertex_maxPt->Draw("HIST");
  h_correct_trackVertex_pt->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->AddEntry(h_correct_trackVertex_pt,"Correct","l");
  l->AddEntry(h_false_trackVertex_maxPt,"False","l");
  l->Draw();
  c.SaveAs(DIR + "/h_correctVsFalse_pt.pdf");
  delete h_correct_trackVertex_pt;
  delete h_false_trackVertex_maxPt;

  h_correct_trackVertex_delta_dist_eta->Scale(1./h_correct_trackVertex_delta_dist_eta->Integral());
  h_correct_trackVertex_delta_dist_eta->SetLineColor(1);
  h_correct_trackVertex_delta_dist_eta->SetMarkerColor(1);
  h_correct_trackVertex_delta_dist_eta->SetStats(0);
  h_false_trackVertex_delta_dist_eta->Scale(1./h_false_trackVertex_delta_dist_eta->Integral());
  h_false_trackVertex_delta_dist_eta->SetLineColor(2);
  h_false_trackVertex_delta_dist_eta->SetMarkerColor(2);
  h_false_trackVertex_delta_dist_eta->SetStats(0);
  raiseMax(h_false_trackVertex_delta_dist_eta,h_correct_trackVertex_delta_dist_eta);
  h_false_trackVertex_delta_dist_eta->Draw("HIST");
  h_correct_trackVertex_delta_dist_eta->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_correct_trackVertex_delta_dist_eta,"Correct","l");
  l->AddEntry(h_false_trackVertex_delta_dist_eta,"False","l");
  l->Draw();
  c.SaveAs(DIR + "/h_correctVsFalse_delta_dist_eta.pdf");
  delete h_correct_trackVertex_delta_dist_eta;
  delete h_false_trackVertex_delta_dist_eta;

  h_correct_trackVertex_R_T->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trackVertex_R_T);
  h_false_trackVertex_R_T->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trackVertex_R_T);
  h_correct_trackVertex_R_T->Scale(1./h_correct_trackVertex_R_T->Integral());
  h_correct_trackVertex_R_T->SetLineColor(1);
  h_correct_trackVertex_R_T->SetMarkerColor(1);
  h_correct_trackVertex_R_T->SetStats(0);
  h_false_trackVertex_R_T->Scale(1./h_false_trackVertex_R_T->Integral());
  h_false_trackVertex_R_T->SetLineColor(2);
  h_false_trackVertex_R_T->SetMarkerColor(2);
  h_false_trackVertex_R_T->SetStats(0);
  raiseMax(h_false_trackVertex_R_T,h_correct_trackVertex_R_T);
  h_false_trackVertex_R_T->Draw("HIST");
  h_correct_trackVertex_R_T->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_correct_trackVertex_R_T,"Correct","l");
  l->AddEntry(h_false_trackVertex_R_T,"False","l");
  l->Draw();
  c.SaveAs(DIR + "/h_correctVsFalse_R_T.pdf");
  delete h_correct_trackVertex_R_T;
  delete h_false_trackVertex_R_T;

  h_correct_trackVertex_cos_T->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trackVertex_cos_T);
  h_false_trackVertex_cos_T->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trackVertex_cos_T);
  h_correct_trackVertex_cos_T->Scale(1./h_correct_trackVertex_cos_T->Integral());
  h_correct_trackVertex_cos_T->SetLineColor(1);
  h_correct_trackVertex_cos_T->SetMarkerColor(1);
  h_correct_trackVertex_cos_T->SetStats(0);
  h_false_trackVertex_cos_T->Scale(1./h_false_trackVertex_cos_T->Integral());
  h_false_trackVertex_cos_T->SetLineColor(2);
  h_false_trackVertex_cos_T->SetMarkerColor(2);
  h_false_trackVertex_cos_T->SetStats(0);
  raiseMax(h_false_trackVertex_cos_T, h_correct_trackVertex_cos_T);
  h_false_trackVertex_cos_T->Draw("HIST");
  h_correct_trackVertex_cos_T->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_correct_trackVertex_cos_T,"Correct","l");
  l->AddEntry(h_false_trackVertex_cos_T,"False","l");
  l->Draw();
  c.SaveAs(DIR + "/h_correctVsFalse_cos_T.pdf");
  delete h_correct_trackVertex_cos_T;
  delete h_false_trackVertex_cos_T;

  h_correct_trackVertex_p2_mag->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trackVertex_p2_mag);
  h_false_trackVertex_p2_mag->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trackVertex_p2_mag);
  h_correct_trackVertex_p2_mag->Scale(1./h_correct_trackVertex_p2_mag->Integral());
  h_correct_trackVertex_p2_mag->SetLineColor(1);
  h_correct_trackVertex_p2_mag->SetMarkerColor(1);
  h_false_trackVertex_p2_mag->Scale(1./h_false_trackVertex_p2_mag->Integral());
  h_false_trackVertex_p2_mag->SetLineColor(2);
  h_false_trackVertex_p2_mag->SetMarkerColor(2);
  raiseMax(h_correct_trackVertex_p2_mag,h_false_trackVertex_p2_mag);
  h_correct_trackVertex_p2_mag->Draw("HIST");
  h_false_trackVertex_p2_mag->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_correct_trackVertex_p2_mag,"Correct","l");
  l->AddEntry(h_false_trackVertex_p2_mag,"False","l");
  l->Draw();
  c.SaveAs(DIR + "/h_correctVsFalse_p2_mag.pdf");
  delete h_correct_trackVertex_p2_mag;
  delete h_false_trackVertex_p2_mag;

  h_correct_trackVertex_lowPt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trackVertex_lowPt);
  h_false_trackVertex_lowPt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trackVertex_lowPt);
  h_correct_trackVertex_lowPt->Scale(1./h_correct_trackVertex_lowPt->Integral());
  h_correct_trackVertex_lowPt->SetLineColor(1);
  h_correct_trackVertex_lowPt->SetMarkerColor(1);
  h_false_trackVertex_lowPt->Scale(1./h_false_trackVertex_lowPt->Integral());
  h_false_trackVertex_lowPt->SetLineColor(2);
  h_false_trackVertex_lowPt->SetMarkerColor(2);
  raiseMax(h_correct_trackVertex_lowPt,h_false_trackVertex_lowPt);
  h_correct_trackVertex_lowPt->Draw("HIST");
  h_false_trackVertex_lowPt->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_correct_trackVertex_lowPt,"Correct","l");
  l->AddEntry(h_false_trackVertex_lowPt,"False","l");
  l->Draw();
  c.SaveAs(DIR + "/h_correctVsFalse_lowPt.pdf");
  delete h_correct_trackVertex_lowPt;
  delete h_false_trackVertex_lowPt;

  h_correct_trackVertex_delta_dist_z->Scale(1./h_correct_trackVertex_delta_dist_z->Integral());
  h_correct_trackVertex_delta_dist_z->SetLineColor(1);
  h_correct_trackVertex_delta_dist_z->SetMarkerColor(1);
  h_correct_trackVertex_delta_dist_z->SetStats(0);
  h_false_trackVertex_delta_dist_z->Scale(1./h_false_trackVertex_delta_dist_z->Integral());
  h_false_trackVertex_delta_dist_z->SetLineColor(2);
  h_false_trackVertex_delta_dist_z->SetMarkerColor(2);
  h_false_trackVertex_delta_dist_z->SetStats(0);
  raiseMax(h_false_trackVertex_delta_dist_z,h_correct_trackVertex_delta_dist_z);
  h_false_trackVertex_delta_dist_z->Draw("HIST");
  h_correct_trackVertex_delta_dist_z->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_correct_trackVertex_delta_dist_z,"Correct","l");
  l->AddEntry(h_false_trackVertex_delta_dist_z,"False","l");
  l->Draw();
  c.SaveAs(DIR + "/h_correctVsFalse_delta_dist_z.pdf");
  delete h_correct_trackVertex_delta_dist_z;
  delete h_false_trackVertex_delta_dist_z;

  h_correct_trackVertex_delta_dist_z_inBothTraj->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trackVertex_delta_dist_z_inBothTraj);
  h_false_trackVertex_delta_dist_z_inBothTraj->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trackVertex_delta_dist_z_inBothTraj);
  h_correct_trackVertex_delta_dist_z_inBothTraj->Scale(1./h_correct_trackVertex_delta_dist_z_inBothTraj->Integral());
  h_correct_trackVertex_delta_dist_z_inBothTraj->SetLineColor(1);
  h_correct_trackVertex_delta_dist_z_inBothTraj->SetMarkerColor(1);
  h_false_trackVertex_delta_dist_z_inBothTraj->Scale(1./h_false_trackVertex_delta_dist_z_inBothTraj->Integral());
  h_false_trackVertex_delta_dist_z_inBothTraj->SetLineColor(2);
  h_false_trackVertex_delta_dist_z_inBothTraj->SetMarkerColor(2);
  raiseMax(h_correct_trackVertex_delta_dist_z_inBothTraj,h_false_trackVertex_delta_dist_z_inBothTraj);
  h_correct_trackVertex_delta_dist_z_inBothTraj->Draw("HIST");
  h_false_trackVertex_delta_dist_z_inBothTraj->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_correct_trackVertex_delta_dist_z_inBothTraj,"Correct","l");
  l->AddEntry(h_false_trackVertex_delta_dist_z_inBothTraj,"False","l");
  l->Draw();
  c.SaveAs(DIR + "/h_correctVsFalse_delta_dist_z_inBothTraj.pdf");
  delete h_correct_trackVertex_delta_dist_z_inBothTraj;
  delete h_false_trackVertex_delta_dist_z_inBothTraj;

  h_correct_trackVertex_delta_dist_z_inOneTraj->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trackVertex_delta_dist_z_inOneTraj);
  h_false_trackVertex_delta_dist_z_inOneTraj->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trackVertex_delta_dist_z_inOneTraj);
  h_correct_trackVertex_delta_dist_z_inOneTraj->Scale(1./h_correct_trackVertex_delta_dist_z_inOneTraj->Integral());
  h_correct_trackVertex_delta_dist_z_inOneTraj->SetLineColor(1);
  h_correct_trackVertex_delta_dist_z_inOneTraj->SetMarkerColor(1);
  h_false_trackVertex_delta_dist_z_inOneTraj->Scale(1./h_false_trackVertex_delta_dist_z_inOneTraj->Integral());
  h_false_trackVertex_delta_dist_z_inOneTraj->SetLineColor(2);
  h_false_trackVertex_delta_dist_z_inOneTraj->SetMarkerColor(2);
  raiseMax(h_correct_trackVertex_delta_dist_z_inOneTraj,h_false_trackVertex_delta_dist_z_inOneTraj);
  h_correct_trackVertex_delta_dist_z_inOneTraj->Draw("HIST");
  h_false_trackVertex_delta_dist_z_inOneTraj->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_correct_trackVertex_delta_dist_z_inOneTraj,"Correct","l");
  l->AddEntry(h_false_trackVertex_delta_dist_z_inOneTraj,"False","l");
  l->Draw();
  c.SaveAs(DIR + "/h_correctVsFalse_delta_dist_z_inOneTraj.pdf");
  delete h_correct_trackVertex_delta_dist_z_inOneTraj;
  delete h_false_trackVertex_delta_dist_z_inOneTraj;

  h_correct_trackVertex_delta_dist_z_inNoTraj->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_correct_trackVertex_delta_dist_z_inNoTraj);
  h_false_trackVertex_delta_dist_z_inNoTraj->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_false_trackVertex_delta_dist_z_inNoTraj);
  h_correct_trackVertex_delta_dist_z_inNoTraj->Scale(1./h_correct_trackVertex_delta_dist_z_inNoTraj->Integral());
  h_correct_trackVertex_delta_dist_z_inNoTraj->SetLineColor(1);
  h_correct_trackVertex_delta_dist_z_inNoTraj->SetMarkerColor(1);
  h_false_trackVertex_delta_dist_z_inNoTraj->Scale(1./h_false_trackVertex_delta_dist_z_inNoTraj->Integral());
  h_false_trackVertex_delta_dist_z_inNoTraj->SetLineColor(2);
  h_false_trackVertex_delta_dist_z_inNoTraj->SetMarkerColor(2);
  raiseMax(h_correct_trackVertex_delta_dist_z_inNoTraj,h_false_trackVertex_delta_dist_z_inNoTraj);
  h_correct_trackVertex_delta_dist_z_inNoTraj->Draw("HIST");
  h_false_trackVertex_delta_dist_z_inNoTraj->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_correct_trackVertex_delta_dist_z_inNoTraj,"Correct","l");
  l->AddEntry(h_false_trackVertex_delta_dist_z_inNoTraj,"False","l");
  l->Draw();
  c.SaveAs(DIR + "/h_correctVsFalse_delta_dist_z_inNoTraj.pdf");
  delete h_correct_trackVertex_delta_dist_z_inNoTraj;
  delete h_false_trackVertex_delta_dist_z_inNoTraj;
  
  h_correct_trackVertex_chi2rphidofSum->Scale(1./h_correct_trackVertex_chi2rphidofSum->Integral());
  h_correct_trackVertex_chi2rphidofSum->SetLineColor(1);
  h_correct_trackVertex_chi2rphidofSum->SetMarkerColor(1);
  h_correct_trackVertex_chi2rphidofSum->SetStats(0);
  h_false_trackVertex_chi2rphidofSum->Scale(1./h_false_trackVertex_chi2rphidofSum->Integral());
  h_false_trackVertex_chi2rphidofSum->SetLineColor(2);
  h_false_trackVertex_chi2rphidofSum->SetMarkerColor(2);
  h_false_trackVertex_chi2rphidofSum->SetStats(0);
  raiseMax(h_false_trackVertex_chi2rphidofSum,h_correct_trackVertex_chi2rphidofSum);
  h_false_trackVertex_chi2rphidofSum->Draw("HIST");
  h_correct_trackVertex_chi2rphidofSum->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_correct_trackVertex_chi2rphidofSum,"Correct","l");
  l->AddEntry(h_false_trackVertex_chi2rphidofSum,"False","l");
  l->Draw();
  c.SaveAs(DIR + "/h_correctVsFalse_chi2rphidofSum.pdf");
  delete h_correct_trackVertex_chi2rphidofSum;
  delete h_false_trackVertex_chi2rphidofSum;

  h_correct_trackVertex_MVA1Sum->Scale(1./h_correct_trackVertex_MVA1Sum->Integral());
  h_correct_trackVertex_MVA1Sum->SetLineColor(1);
  h_correct_trackVertex_MVA1Sum->SetMarkerColor(1);
  h_correct_trackVertex_MVA1Sum->SetStats(0);
  h_false_trackVertex_MVA1Sum->Scale(1./h_false_trackVertex_MVA1Sum->Integral());
  h_false_trackVertex_MVA1Sum->SetLineColor(2);
  h_false_trackVertex_MVA1Sum->SetMarkerColor(2);
  h_false_trackVertex_MVA1Sum->SetStats(0);
  raiseMax(h_false_trackVertex_MVA1Sum,h_correct_trackVertex_MVA1Sum);
  h_false_trackVertex_MVA1Sum->Draw("HIST");
  h_correct_trackVertex_MVA1Sum->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_correct_trackVertex_MVA1Sum,"Correct","l");
  l->AddEntry(h_false_trackVertex_MVA1Sum,"False","l");
  l->Draw();
  c.SaveAs(DIR + "/h_correctVsFalse_MVA1Sum.pdf");
  delete h_correct_trackVertex_MVA1Sum;
  delete h_false_trackVertex_MVA1Sum;

  h_correct_trackVertex_MVA2Sum->Scale(1./h_correct_trackVertex_MVA2Sum->Integral());
  h_correct_trackVertex_MVA2Sum->SetLineColor(1);
  h_correct_trackVertex_MVA2Sum->SetMarkerColor(1);
  h_correct_trackVertex_MVA2Sum->SetStats(0);
  h_false_trackVertex_MVA2Sum->Scale(1./h_false_trackVertex_MVA2Sum->Integral());
  h_false_trackVertex_MVA2Sum->SetLineColor(2);
  h_false_trackVertex_MVA2Sum->SetMarkerColor(2);
  h_false_trackVertex_MVA2Sum->SetStats(0);
  raiseMax(h_false_trackVertex_MVA2Sum,h_correct_trackVertex_MVA2Sum);
  h_false_trackVertex_MVA2Sum->Draw("HIST");
  h_correct_trackVertex_MVA2Sum->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_correct_trackVertex_MVA2Sum,"Correct","l");
  l->AddEntry(h_false_trackVertex_MVA2Sum,"False","l");
  l->Draw();
  c.SaveAs(DIR + "/h_correctVsFalse_MVA2Sum.pdf");
  delete h_correct_trackVertex_MVA2Sum;
  delete h_false_trackVertex_MVA2Sum;

  h_correct_trackVertex_chi2rzdofSum->Scale(1./h_correct_trackVertex_chi2rzdofSum->Integral());
  h_correct_trackVertex_chi2rzdofSum->SetLineColor(1);
  h_correct_trackVertex_chi2rzdofSum->SetMarkerColor(1);
  h_correct_trackVertex_chi2rzdofSum->SetStats(0);
  h_false_trackVertex_chi2rzdofSum->Scale(1./h_false_trackVertex_chi2rzdofSum->Integral());
  h_false_trackVertex_chi2rzdofSum->SetLineColor(2);
  h_false_trackVertex_chi2rzdofSum->SetMarkerColor(2);
  h_false_trackVertex_chi2rzdofSum->SetStats(0);
  raiseMax(h_false_trackVertex_chi2rzdofSum,h_correct_trackVertex_chi2rzdofSum);
  h_false_trackVertex_chi2rzdofSum->Draw("HIST");
  h_correct_trackVertex_chi2rzdofSum->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_correct_trackVertex_chi2rzdofSum,"Correct","l");
  l->AddEntry(h_false_trackVertex_chi2rzdofSum,"False","l");
  l->Draw();
  c.SaveAs(DIR + "/h_correctVsFalse_chi2rzdofSum.pdf");
  delete h_correct_trackVertex_chi2rzdofSum;
  delete h_false_trackVertex_chi2rzdofSum;

  h_correct_trackVertex_bendchi2Sum->Scale(1./h_correct_trackVertex_bendchi2Sum->Integral());
  h_correct_trackVertex_bendchi2Sum->SetLineColor(1);
  h_correct_trackVertex_bendchi2Sum->SetMarkerColor(1);
  h_correct_trackVertex_bendchi2Sum->SetStats(0);
  h_false_trackVertex_bendchi2Sum->Scale(1./h_false_trackVertex_bendchi2Sum->Integral());
  h_false_trackVertex_bendchi2Sum->SetLineColor(2);
  h_false_trackVertex_bendchi2Sum->SetMarkerColor(2);
  h_false_trackVertex_bendchi2Sum->SetStats(0);
  raiseMax(h_false_trackVertex_bendchi2Sum,h_correct_trackVertex_bendchi2Sum);
  h_false_trackVertex_bendchi2Sum->Draw("HIST");
  h_correct_trackVertex_bendchi2Sum->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_correct_trackVertex_bendchi2Sum,"Correct","l");
  l->AddEntry(h_false_trackVertex_bendchi2Sum,"False","l");
  l->Draw();
  c.SaveAs(DIR + "/h_correctVsFalse_bendchi2Sum.pdf");
  delete h_correct_trackVertex_bendchi2Sum;
  delete h_false_trackVertex_bendchi2Sum;

  h_correct_trackVertex_numStubsSum->Scale(1./h_correct_trackVertex_numStubsSum->Integral());
  h_correct_trackVertex_numStubsSum->SetLineColor(1);
  h_correct_trackVertex_numStubsSum->SetMarkerColor(1);
  h_false_trackVertex_numStubsSum->Scale(1./h_false_trackVertex_numStubsSum->Integral());
  h_false_trackVertex_numStubsSum->SetLineColor(2);
  h_false_trackVertex_numStubsSum->SetMarkerColor(2);
  raiseMax(h_false_trackVertex_numStubsSum,h_correct_trackVertex_numStubsSum);
  h_false_trackVertex_numStubsSum->Draw("HIST");
  h_correct_trackVertex_numStubsSum->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_correct_trackVertex_numStubsSum,"Correct","l");
  l->AddEntry(h_false_trackVertex_numStubsSum,"False","l");
  l->Draw();
  c.SaveAs(DIR + "/h_correctVsFalse_numStubsSum.pdf");
  delete h_correct_trackVertex_numStubsSum;
  delete h_false_trackVertex_numStubsSum;

  h_trk_chi2rphidof_primary_noCuts->Scale(1./h_trk_chi2rphidof_primary_noCuts->Integral());
  h_trk_chi2rphidof_primary_noCuts->SetLineColor(1);
  h_trk_chi2rphidof_primary_noCuts->SetMarkerColor(1);
  TH1F *h_trk_chi2rphidof_np_noCutsNorm = (TH1F*)h_trk_chi2rphidof_np_noCuts->Clone(); 
  h_trk_chi2rphidof_np_noCutsNorm->Scale(1./h_trk_chi2rphidof_np_noCutsNorm->Integral());
  h_trk_chi2rphidof_np_noCutsNorm->SetLineColor(2);
  h_trk_chi2rphidof_np_noCutsNorm->SetMarkerColor(2);
  raiseMax(h_trk_chi2rphidof_np_noCutsNorm,h_trk_chi2rphidof_primary_noCuts);
  h_trk_chi2rphidof_np_noCutsNorm->SetStats(0);
  h_trk_chi2rphidof_primary_noCuts->SetStats(0);
  h_trk_chi2rphidof_np_noCutsNorm->Draw("HIST");
  h_trk_chi2rphidof_primary_noCuts->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_chi2rphidof_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_chi2rphidof_np_noCutsNorm,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_chi2rphidof_noCuts.pdf");
  delete h_trk_chi2rphidof_np_noCutsNorm;

  h_trk_chi2rphidof_primary_noCuts_zoomOut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_primary_noCuts_zoomOut);
  h_trk_chi2rphidof_np_noCuts_zoomOut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_np_noCuts_zoomOut);
  h_trk_chi2rphidof_primary_noCuts_zoomOut->Scale(1./h_trk_chi2rphidof_primary_noCuts_zoomOut->Integral());
  h_trk_chi2rphidof_primary_noCuts_zoomOut->SetLineColor(1);
  h_trk_chi2rphidof_primary_noCuts_zoomOut->SetMarkerColor(1);
  h_trk_chi2rphidof_np_noCuts_zoomOut->Scale(1./h_trk_chi2rphidof_np_noCuts_zoomOut->Integral());
  h_trk_chi2rphidof_np_noCuts_zoomOut->SetLineColor(2);
  h_trk_chi2rphidof_np_noCuts_zoomOut->SetMarkerColor(2);
  raiseMax(h_trk_chi2rphidof_np_noCuts_zoomOut,h_trk_chi2rphidof_primary_noCuts_zoomOut);
  h_trk_chi2rphidof_np_noCuts_zoomOut->SetStats(0);
  h_trk_chi2rphidof_primary_noCuts_zoomOut->SetStats(0);
  h_trk_chi2rphidof_np_noCuts_zoomOut->Draw("HIST");
  h_trk_chi2rphidof_primary_noCuts_zoomOut->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->SetFillColor(0);
  l->SetLineColor(0);
  l->SetTextSize(0.04);
  l->AddEntry(h_trk_chi2rphidof_primary_noCuts_zoomOut,"Primary","l");
  l->AddEntry(h_trk_chi2rphidof_np_noCuts_zoomOut,"NP","l");
  l->SetTextFont(42);
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_chi2rphidof_noCuts_zoomOut.pdf");
  delete h_trk_chi2rphidof_primary_noCuts_zoomOut;
  delete h_trk_chi2rphidof_np_noCuts_zoomOut;

  TH1F *h_trk_chi2rphidof_fake_noCutsNorm = (TH1F*)h_trk_chi2rphidof_fake_noCuts->Clone(); 
  h_trk_chi2rphidof_fake_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_fake_noCutsNorm);
  TH1F *h_trk_chi2rphidof_PU_noCutsNorm = (TH1F*)h_trk_chi2rphidof_PU_noCuts->Clone(); 
  h_trk_chi2rphidof_PU_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_PU_noCutsNorm);
  TH1F *h_trk_chi2rphidof_notHiggs_noCutsNorm = (TH1F*)h_trk_chi2rphidof_notHiggs_noCuts->Clone(); 
  h_trk_chi2rphidof_notHiggs_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_notHiggs_noCutsNorm);
  h_trk_chi2rphidof_fake_noCutsNorm->Scale(1./h_trk_chi2rphidof_np_noCuts->Integral());
  h_trk_chi2rphidof_fake_noCutsNorm->SetLineColor(2);
  h_trk_chi2rphidof_fake_noCutsNorm->SetMarkerColor(2);
  h_trk_chi2rphidof_PU_noCutsNorm->Scale(1./h_trk_chi2rphidof_np_noCuts->Integral());
  h_trk_chi2rphidof_PU_noCutsNorm->SetLineColor(3);
  h_trk_chi2rphidof_PU_noCutsNorm->SetMarkerColor(3);
  h_trk_chi2rphidof_notHiggs_noCutsNorm->Scale(1./h_trk_chi2rphidof_np_noCuts->Integral());
  h_trk_chi2rphidof_notHiggs_noCutsNorm->SetLineColor(4);
  h_trk_chi2rphidof_notHiggs_noCutsNorm->SetMarkerColor(4);
  auto hs_chi2rphidof_noCuts = new THStack("hs_chi2rphidof_noCuts","Stacked BG histograms");
  hs_chi2rphidof_noCuts->Add(h_trk_chi2rphidof_fake_noCutsNorm);
  hs_chi2rphidof_noCuts->Add(h_trk_chi2rphidof_PU_noCutsNorm);
  hs_chi2rphidof_noCuts->Add(h_trk_chi2rphidof_notHiggs_noCutsNorm);
  hs_chi2rphidof_noCuts->Draw("HIST");
  h_trk_chi2rphidof_primary_noCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_chi2rphidof_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_chi2rphidof_fake_noCutsNorm,"Fake","l");
  l->AddEntry(h_trk_chi2rphidof_PU_noCutsNorm,"PU","l");
  l->AddEntry(h_trk_chi2rphidof_notHiggs_noCutsNorm,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_chi2rphidof_noCutsNorm.pdf");
  delete h_trk_chi2rphidof_fake_noCutsNorm;
  delete h_trk_chi2rphidof_PU_noCutsNorm;
  delete h_trk_chi2rphidof_notHiggs_noCutsNorm;
  delete h_trk_chi2rphidof_np_noCuts;
  delete hs_chi2rphidof_noCuts;

  h_trk_chi2rphidof_fake_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_fake_noCuts);
  h_trk_chi2rphidof_PU_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_PU_noCuts);
  h_trk_chi2rphidof_notHiggs_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_notHiggs_noCuts);
  h_trk_chi2rphidof_fake_noCuts->Scale(1./h_trk_chi2rphidof_fake_noCuts->Integral());
  h_trk_chi2rphidof_fake_noCuts->SetLineColor(2);
  h_trk_chi2rphidof_fake_noCuts->SetMarkerColor(2);
  h_trk_chi2rphidof_PU_noCuts->Scale(1./h_trk_chi2rphidof_PU_noCuts->Integral());
  h_trk_chi2rphidof_PU_noCuts->SetLineColor(3);
  h_trk_chi2rphidof_PU_noCuts->SetMarkerColor(3);
  h_trk_chi2rphidof_notHiggs_noCuts->Scale(1./h_trk_chi2rphidof_notHiggs_noCuts->Integral());
  h_trk_chi2rphidof_notHiggs_noCuts->SetLineColor(4);
  h_trk_chi2rphidof_notHiggs_noCuts->SetMarkerColor(4);
  h_trk_chi2rphidof_fake_noCuts->Draw("HIST");
  h_trk_chi2rphidof_primary_noCuts->Draw("HIST,SAME");
  h_trk_chi2rphidof_PU_noCuts->Draw("HIST,SAME");
  h_trk_chi2rphidof_notHiggs_noCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_chi2rphidof_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_chi2rphidof_fake_noCuts,"Fake","l");
  l->AddEntry(h_trk_chi2rphidof_PU_noCuts,"PU","l");
  l->AddEntry(h_trk_chi2rphidof_notHiggs_noCuts,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_chi2rphidof_noCuts.pdf");
  delete h_trk_chi2rphidof_primary_noCuts;
  delete h_trk_chi2rphidof_fake_noCuts;
  delete h_trk_chi2rphidof_PU_noCuts;
  delete h_trk_chi2rphidof_notHiggs_noCuts;

  h_trk_chi2rphidof_primary_noCuts_H->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_primary_noCuts_H);
  h_trk_chi2rphidof_np_noCuts_H->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_np_noCuts_H);
  h_trk_chi2rphidof_primary_noCuts_H->Scale(1./h_trk_chi2rphidof_primary_noCuts_H->Integral());
  h_trk_chi2rphidof_primary_noCuts_H->SetLineColor(1);
  h_trk_chi2rphidof_primary_noCuts_H->SetMarkerColor(1);
  h_trk_chi2rphidof_np_noCuts_H->Scale(1./h_trk_chi2rphidof_np_noCuts_H->Integral());
  h_trk_chi2rphidof_np_noCuts_H->SetLineColor(2);
  h_trk_chi2rphidof_np_noCuts_H->SetMarkerColor(2);
  h_trk_chi2rphidof_np_noCuts_H->Draw("HIST");
  h_trk_chi2rphidof_primary_noCuts_H->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_chi2rphidof_primary_noCuts_H,"Primary","l");
  l->AddEntry(h_trk_chi2rphidof_np_noCuts_H,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_chi2rphidof_noCuts_H.pdf");
  delete h_trk_chi2rphidof_primary_noCuts_H;
  delete h_trk_chi2rphidof_np_noCuts_H;

  h_trk_chi2rphidof_primary_noCuts_L->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_primary_noCuts_L);
  h_trk_chi2rphidof_np_noCuts_L->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_np_noCuts_L);
  h_trk_chi2rphidof_primary_noCuts_L->Scale(1./h_trk_chi2rphidof_primary_noCuts_L->Integral());
  h_trk_chi2rphidof_primary_noCuts_L->SetLineColor(1);
  h_trk_chi2rphidof_primary_noCuts_L->SetMarkerColor(1);
  h_trk_chi2rphidof_np_noCuts_L->Scale(1./h_trk_chi2rphidof_np_noCuts_L->Integral());
  h_trk_chi2rphidof_np_noCuts_L->SetLineColor(2);
  h_trk_chi2rphidof_np_noCuts_L->SetMarkerColor(2);
  h_trk_chi2rphidof_np_noCuts_L->Draw("HIST");
  h_trk_chi2rphidof_primary_noCuts_L->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_chi2rphidof_primary_noCuts_L,"Primary","l");
  l->AddEntry(h_trk_chi2rphidof_np_noCuts_L,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_chi2rphidof_noCuts_L.pdf");
  delete h_trk_chi2rphidof_primary_noCuts_L;
  delete h_trk_chi2rphidof_np_noCuts_L;

  h_trk_chi2rphidof_primary_noCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_primary_noCuts_barrel);
  h_trk_chi2rphidof_np_noCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_np_noCuts_barrel);
  h_trk_chi2rphidof_primary_noCuts_barrel->Scale(1./h_trk_chi2rphidof_primary_noCuts_barrel->Integral());
  h_trk_chi2rphidof_primary_noCuts_barrel->SetLineColor(1);
  h_trk_chi2rphidof_primary_noCuts_barrel->SetMarkerColor(1);
  h_trk_chi2rphidof_np_noCuts_barrel->Scale(1./h_trk_chi2rphidof_np_noCuts_barrel->Integral());
  h_trk_chi2rphidof_np_noCuts_barrel->SetLineColor(2);
  h_trk_chi2rphidof_np_noCuts_barrel->SetMarkerColor(2);
  h_trk_chi2rphidof_np_noCuts_barrel->Draw("HIST");
  h_trk_chi2rphidof_primary_noCuts_barrel->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_chi2rphidof_primary_noCuts_barrel,"Primary","l");
  l->AddEntry(h_trk_chi2rphidof_np_noCuts_barrel,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_chi2rphidof_noCuts_barrel.pdf");
  delete h_trk_chi2rphidof_primary_noCuts_barrel;
  delete h_trk_chi2rphidof_np_noCuts_barrel;

  h_trk_chi2rphidof_primary_noCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_primary_noCuts_disk);
  h_trk_chi2rphidof_np_noCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_np_noCuts_disk);
  h_trk_chi2rphidof_primary_noCuts_disk->Scale(1./h_trk_chi2rphidof_primary_noCuts_disk->Integral());
  h_trk_chi2rphidof_primary_noCuts_disk->SetLineColor(1);
  h_trk_chi2rphidof_primary_noCuts_disk->SetMarkerColor(1);
  h_trk_chi2rphidof_np_noCuts_disk->Scale(1./h_trk_chi2rphidof_np_noCuts_disk->Integral());
  h_trk_chi2rphidof_np_noCuts_disk->SetLineColor(2);
  h_trk_chi2rphidof_np_noCuts_disk->SetMarkerColor(2);
  h_trk_chi2rphidof_np_noCuts_disk->Draw("HIST");
  h_trk_chi2rphidof_primary_noCuts_disk->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_chi2rphidof_primary_noCuts_disk,"Primary","l");
  l->AddEntry(h_trk_chi2rphidof_np_noCuts_disk,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_chi2rphidof_noCuts_disk.pdf");
  delete h_trk_chi2rphidof_primary_noCuts_disk;
  delete h_trk_chi2rphidof_np_noCuts_disk;

  h_trk_chi2rphidof_primary_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_primary_qualCuts);
  h_trk_chi2rphidof_np_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_np_qualCuts);
  h_trk_chi2rphidof_primary_qualCuts->Scale(1./h_trk_chi2rphidof_primary_qualCuts->Integral());
  h_trk_chi2rphidof_primary_qualCuts->SetLineColor(1);
  h_trk_chi2rphidof_primary_qualCuts->SetMarkerColor(1);
  TH1F *h_trk_chi2rphidof_np_qualCutsNorm = (TH1F*)h_trk_chi2rphidof_np_qualCuts->Clone(); 
  h_trk_chi2rphidof_np_qualCutsNorm->Scale(1./h_trk_chi2rphidof_np_qualCutsNorm->Integral());
  h_trk_chi2rphidof_np_qualCutsNorm->SetLineColor(2);
  h_trk_chi2rphidof_np_qualCutsNorm->SetMarkerColor(2);
  h_trk_chi2rphidof_np_qualCutsNorm->Draw("HIST");
  h_trk_chi2rphidof_primary_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_chi2rphidof_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_chi2rphidof_np_qualCutsNorm,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_chi2rphidof_qualCuts.pdf");
  delete h_trk_chi2rphidof_np_qualCutsNorm;

  h_trk_chi2rphidof_primary_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_primary_allCuts);
  h_trk_chi2rphidof_np_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_np_allCuts);
  h_trk_chi2rphidof_primary_allCuts->Scale(1./h_trk_chi2rphidof_primary_allCuts->Integral());
  h_trk_chi2rphidof_primary_allCuts->SetLineColor(1);
  h_trk_chi2rphidof_primary_allCuts->SetMarkerColor(1);
  TH1F *h_trk_chi2rphidof_np_allCutsNorm = (TH1F*)h_trk_chi2rphidof_np_allCuts->Clone(); 
  h_trk_chi2rphidof_np_allCutsNorm->Scale(1./h_trk_chi2rphidof_np_allCutsNorm->Integral());
  h_trk_chi2rphidof_np_allCutsNorm->SetLineColor(2);
  h_trk_chi2rphidof_np_allCutsNorm->SetMarkerColor(2);
  raiseMax(h_trk_chi2rphidof_primary_allCuts,h_trk_chi2rphidof_np_allCutsNorm);
  h_trk_chi2rphidof_np_allCutsNorm->SetStats(0);
  h_trk_chi2rphidof_primary_allCuts->SetStats(0);
  h_trk_chi2rphidof_primary_allCuts->Draw("HIST");
  h_trk_chi2rphidof_np_allCutsNorm->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_chi2rphidof_primary_allCuts,"Primary","l");
  l->AddEntry(h_trk_chi2rphidof_np_allCutsNorm,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_chi2rphidof_allCuts.pdf");
  delete h_trk_chi2rphidof_np_allCutsNorm;

  h_trk_chi2rphidof_primary_allCuts_zoomOut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_primary_allCuts_zoomOut);
  h_trk_chi2rphidof_np_allCuts_zoomOut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_np_allCuts_zoomOut);
  h_trk_chi2rphidof_primary_allCuts_zoomOut->Scale(1./h_trk_chi2rphidof_primary_allCuts_zoomOut->Integral());
  h_trk_chi2rphidof_primary_allCuts_zoomOut->SetLineColor(1);
  h_trk_chi2rphidof_primary_allCuts_zoomOut->SetMarkerColor(1);
  h_trk_chi2rphidof_np_allCuts_zoomOut->Scale(1./h_trk_chi2rphidof_np_allCuts_zoomOut->Integral());
  h_trk_chi2rphidof_np_allCuts_zoomOut->SetLineColor(2);
  h_trk_chi2rphidof_np_allCuts_zoomOut->SetMarkerColor(2);
  raiseMax(h_trk_chi2rphidof_np_allCuts_zoomOut,h_trk_chi2rphidof_primary_allCuts_zoomOut);
  h_trk_chi2rphidof_np_allCuts_zoomOut->SetStats(0);
  h_trk_chi2rphidof_primary_allCuts_zoomOut->SetStats(0);
  h_trk_chi2rphidof_np_allCuts_zoomOut->Draw("HIST");
  h_trk_chi2rphidof_primary_allCuts_zoomOut->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->SetFillColor(0);
  l->SetLineColor(0);
  l->SetTextSize(0.04);
  l->AddEntry(h_trk_chi2rphidof_primary_allCuts_zoomOut,"Primary","l");
  l->AddEntry(h_trk_chi2rphidof_np_allCuts_zoomOut,"NP","l");
  l->SetTextFont(42);
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_chi2rphidof_allCuts_zoomOut.pdf");
  delete h_trk_chi2rphidof_primary_allCuts_zoomOut;
  delete h_trk_chi2rphidof_np_allCuts_zoomOut;

  h_trk_chi2rphidof_primary_allCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_primary_allCuts_barrel);
  h_trk_chi2rphidof_np_allCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_np_allCuts_barrel);
  h_trk_chi2rphidof_primary_allCuts_barrel->Scale(1./h_trk_chi2rphidof_primary_allCuts_barrel->Integral());
  h_trk_chi2rphidof_primary_allCuts_barrel->SetLineColor(1);
  h_trk_chi2rphidof_primary_allCuts_barrel->SetMarkerColor(1);
  h_trk_chi2rphidof_primary_allCuts_barrel->SetStats(0);
  h_trk_chi2rphidof_np_allCuts_barrel->Scale(1./h_trk_chi2rphidof_np_allCuts_barrel->Integral());
  h_trk_chi2rphidof_np_allCuts_barrel->SetLineColor(2);
  h_trk_chi2rphidof_np_allCuts_barrel->SetMarkerColor(2);
  h_trk_chi2rphidof_np_allCuts_barrel->SetStats(0);
  h_trk_chi2rphidof_np_allCuts_barrel->Draw("HIST");
  h_trk_chi2rphidof_primary_allCuts_barrel->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_chi2rphidof_primary_allCuts_barrel,"Primary","l");
  l->AddEntry(h_trk_chi2rphidof_np_allCuts_barrel,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_chi2rphidof_allCuts_barrel.pdf");
  delete h_trk_chi2rphidof_primary_allCuts_barrel;
  delete h_trk_chi2rphidof_np_allCuts_barrel;

  h_trk_chi2rphidof_primary_allCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_primary_allCuts_disk);
  h_trk_chi2rphidof_np_allCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_np_allCuts_disk);
  h_trk_chi2rphidof_primary_allCuts_disk->Scale(1./h_trk_chi2rphidof_primary_allCuts_disk->Integral());
  h_trk_chi2rphidof_primary_allCuts_disk->SetLineColor(1);
  h_trk_chi2rphidof_primary_allCuts_disk->SetMarkerColor(1);
  h_trk_chi2rphidof_primary_allCuts_disk->SetStats(0);
  h_trk_chi2rphidof_np_allCuts_disk->Scale(1./h_trk_chi2rphidof_np_allCuts_disk->Integral());
  h_trk_chi2rphidof_np_allCuts_disk->SetLineColor(2);
  h_trk_chi2rphidof_np_allCuts_disk->SetMarkerColor(2);
  h_trk_chi2rphidof_np_allCuts_disk->SetStats(0);
  h_trk_chi2rphidof_np_allCuts_disk->Draw("HIST");
  h_trk_chi2rphidof_primary_allCuts_disk->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_chi2rphidof_primary_allCuts_disk,"Primary","l");
  l->AddEntry(h_trk_chi2rphidof_np_allCuts_disk,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_chi2rphidof_allCuts_disk.pdf");
  delete h_trk_chi2rphidof_primary_allCuts_disk;
  delete h_trk_chi2rphidof_np_allCuts_disk;

  TH1F *h_trk_chi2rphidof_fake_qualCutsNorm = (TH1F*)h_trk_chi2rphidof_fake_qualCuts->Clone(); 
  h_trk_chi2rphidof_fake_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_fake_qualCutsNorm);
  TH1F *h_trk_chi2rphidof_PU_qualCutsNorm = (TH1F*)h_trk_chi2rphidof_PU_qualCuts->Clone(); 
  h_trk_chi2rphidof_PU_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_PU_qualCutsNorm);
  TH1F *h_trk_chi2rphidof_notHiggs_qualCutsNorm = (TH1F*)h_trk_chi2rphidof_notHiggs_qualCuts->Clone(); 
  h_trk_chi2rphidof_notHiggs_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_notHiggs_qualCutsNorm);
  h_trk_chi2rphidof_fake_qualCutsNorm->Scale(1./h_trk_chi2rphidof_np_qualCuts->Integral());
  h_trk_chi2rphidof_fake_qualCutsNorm->SetLineColor(2);
  h_trk_chi2rphidof_fake_qualCutsNorm->SetMarkerColor(2);
  h_trk_chi2rphidof_PU_qualCutsNorm->Scale(1./h_trk_chi2rphidof_np_qualCuts->Integral());
  h_trk_chi2rphidof_PU_qualCutsNorm->SetLineColor(3);
  h_trk_chi2rphidof_PU_qualCutsNorm->SetMarkerColor(3);
  h_trk_chi2rphidof_notHiggs_qualCutsNorm->Scale(1./h_trk_chi2rphidof_np_qualCuts->Integral());
  h_trk_chi2rphidof_notHiggs_qualCutsNorm->SetLineColor(4);
  h_trk_chi2rphidof_notHiggs_qualCutsNorm->SetMarkerColor(4);
  auto hs_chi2rphidof_qualCuts = new THStack("hs_chi2rphidof_qualCuts","Stacked BG histograms");
  hs_chi2rphidof_qualCuts->Add(h_trk_chi2rphidof_fake_qualCutsNorm);
  hs_chi2rphidof_qualCuts->Add(h_trk_chi2rphidof_PU_qualCutsNorm);
  hs_chi2rphidof_qualCuts->Add(h_trk_chi2rphidof_notHiggs_qualCutsNorm);
  hs_chi2rphidof_qualCuts->Draw("HIST");
  h_trk_chi2rphidof_primary_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_chi2rphidof_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_chi2rphidof_fake_qualCutsNorm,"Fake","l");
  l->AddEntry(h_trk_chi2rphidof_PU_qualCutsNorm,"PU","l");
  l->AddEntry(h_trk_chi2rphidof_notHiggs_qualCutsNorm,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_chi2rphidof_qualCutsNorm.pdf");
  delete h_trk_chi2rphidof_fake_qualCutsNorm;
  delete h_trk_chi2rphidof_PU_qualCutsNorm;
  delete h_trk_chi2rphidof_notHiggs_qualCutsNorm;
  delete h_trk_chi2rphidof_np_qualCuts;
  delete hs_chi2rphidof_qualCuts;

  h_trk_chi2rphidof_fake_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_fake_qualCuts);
  h_trk_chi2rphidof_PU_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_PU_qualCuts);
  h_trk_chi2rphidof_notHiggs_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_notHiggs_qualCuts);
  h_trk_chi2rphidof_fake_qualCuts->Scale(1./h_trk_chi2rphidof_fake_qualCuts->Integral());
  h_trk_chi2rphidof_fake_qualCuts->SetLineColor(2);
  h_trk_chi2rphidof_fake_qualCuts->SetMarkerColor(2);
  h_trk_chi2rphidof_PU_qualCuts->Scale(1./h_trk_chi2rphidof_PU_qualCuts->Integral());
  h_trk_chi2rphidof_PU_qualCuts->SetLineColor(3);
  h_trk_chi2rphidof_PU_qualCuts->SetMarkerColor(3);
  h_trk_chi2rphidof_notHiggs_qualCuts->Scale(1./h_trk_chi2rphidof_notHiggs_qualCuts->Integral());
  h_trk_chi2rphidof_notHiggs_qualCuts->SetLineColor(4);
  h_trk_chi2rphidof_notHiggs_qualCuts->SetMarkerColor(4);
  h_trk_chi2rphidof_fake_qualCuts->Draw("HIST");
  h_trk_chi2rphidof_primary_qualCuts->Draw("HIST,SAME");
  h_trk_chi2rphidof_PU_qualCuts->Draw("HIST,SAME");
  h_trk_chi2rphidof_notHiggs_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_chi2rphidof_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_chi2rphidof_fake_qualCuts,"Fake","l");
  l->AddEntry(h_trk_chi2rphidof_PU_qualCuts,"PU","l");
  l->AddEntry(h_trk_chi2rphidof_notHiggs_qualCuts,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_chi2rphidof_qualCuts.pdf");
  delete h_trk_chi2rphidof_primary_qualCuts;
  delete h_trk_chi2rphidof_fake_qualCuts;
  delete h_trk_chi2rphidof_PU_qualCuts;
  delete h_trk_chi2rphidof_notHiggs_qualCuts;

  removeFlows(h_trk_MVA1_primary_noCuts);
  h_trk_MVA1_primary_noCuts->Scale(1./h_trk_MVA1_primary_noCuts->Integral());
  h_trk_MVA1_primary_noCuts->SetLineColor(1);
  h_trk_MVA1_primary_noCuts->SetMarkerColor(1);
  TH1F *h_trk_MVA1_np_noCutsNorm = (TH1F*)h_trk_MVA1_np_noCuts->Clone(); 
  h_trk_MVA1_np_noCutsNorm->Scale(1./h_trk_MVA1_np_noCutsNorm->Integral());
  h_trk_MVA1_np_noCutsNorm->SetLineColor(2);
  h_trk_MVA1_np_noCutsNorm->SetMarkerColor(2);
  h_trk_MVA1_np_noCutsNorm->Draw("HIST");
  h_trk_MVA1_primary_noCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_MVA1_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_MVA1_np_noCutsNorm,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_MVA1_noCuts.pdf");
  delete h_trk_MVA1_np_noCutsNorm;

  removeFlows(h_trk_MVA2_primary_noCuts);
  h_trk_MVA2_primary_noCuts->Scale(1./h_trk_MVA2_primary_noCuts->Integral());
  h_trk_MVA2_primary_noCuts->SetLineColor(1);
  h_trk_MVA2_primary_noCuts->SetMarkerColor(1);
  h_trk_MVA2_primary_noCuts->SetStats(0);
  removeFlows(h_trk_MVA2_np_noCuts);
  h_trk_MVA2_np_noCuts->Scale(1./h_trk_MVA2_np_noCuts->Integral());
  h_trk_MVA2_np_noCuts->SetLineColor(2);
  h_trk_MVA2_np_noCuts->SetMarkerColor(2);
  h_trk_MVA2_np_noCuts->SetStats(0);
  raiseMax(h_trk_MVA2_np_noCuts,h_trk_MVA2_primary_noCuts);
  h_trk_MVA2_np_noCuts->Draw("HIST");
  h_trk_MVA2_primary_noCuts->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_MVA2_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_MVA2_np_noCuts,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_MVA2_noCuts.pdf");
  delete h_trk_MVA2_primary_noCuts;
  delete h_trk_MVA2_np_noCuts;

  TH1F *h_trk_MVA1_fake_noCutsNorm = (TH1F*)h_trk_MVA1_fake_noCuts->Clone(); 
  h_trk_MVA1_fake_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_MVA1_fake_noCutsNorm);
  TH1F *h_trk_MVA1_PU_noCutsNorm = (TH1F*)h_trk_MVA1_PU_noCuts->Clone(); 
  h_trk_MVA1_PU_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_MVA1_PU_noCutsNorm);
  TH1F *h_trk_MVA1_notHiggs_noCutsNorm = (TH1F*)h_trk_MVA1_notHiggs_noCuts->Clone(); 
  h_trk_MVA1_notHiggs_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_MVA1_notHiggs_noCutsNorm);
  h_trk_MVA1_fake_noCutsNorm->Scale(1./h_trk_MVA1_np_noCuts->Integral());
  h_trk_MVA1_fake_noCutsNorm->SetLineColor(2);
  h_trk_MVA1_fake_noCutsNorm->SetMarkerColor(2);
  h_trk_MVA1_PU_noCutsNorm->Scale(1./h_trk_MVA1_np_noCuts->Integral());
  h_trk_MVA1_PU_noCutsNorm->SetLineColor(3);
  h_trk_MVA1_PU_noCutsNorm->SetMarkerColor(3);
  h_trk_MVA1_notHiggs_noCutsNorm->Scale(1./h_trk_MVA1_np_noCuts->Integral());
  h_trk_MVA1_notHiggs_noCutsNorm->SetLineColor(4);
  h_trk_MVA1_notHiggs_noCutsNorm->SetMarkerColor(4);
  auto hs_MVA1_noCuts = new THStack("hs_MVA1_noCuts","Stacked BG histograms");
  hs_MVA1_noCuts->Add(h_trk_MVA1_fake_noCutsNorm);
  hs_MVA1_noCuts->Add(h_trk_MVA1_PU_noCutsNorm);
  hs_MVA1_noCuts->Add(h_trk_MVA1_notHiggs_noCutsNorm);
  hs_MVA1_noCuts->Draw("HIST");
  h_trk_MVA1_primary_noCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_MVA1_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_MVA1_fake_noCutsNorm,"Fake","l");
  l->AddEntry(h_trk_MVA1_PU_noCutsNorm,"PU","l");
  l->AddEntry(h_trk_MVA1_notHiggs_noCutsNorm,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_MVA1_noCutsNorm.pdf");
  delete h_trk_MVA1_fake_noCutsNorm;
  delete h_trk_MVA1_PU_noCutsNorm;
  delete h_trk_MVA1_notHiggs_noCutsNorm;
  delete h_trk_MVA1_np_noCuts;
  delete hs_MVA1_noCuts;

  h_trk_MVA1_fake_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_MVA1_fake_noCuts);
  h_trk_MVA1_PU_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_MVA1_PU_noCuts);
  h_trk_MVA1_notHiggs_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_MVA1_notHiggs_noCuts);
  h_trk_MVA1_fake_noCuts->Scale(1./h_trk_MVA1_fake_noCuts->Integral());
  h_trk_MVA1_fake_noCuts->SetLineColor(2);
  h_trk_MVA1_fake_noCuts->SetMarkerColor(2);
  h_trk_MVA1_PU_noCuts->Scale(1./h_trk_MVA1_PU_noCuts->Integral());
  h_trk_MVA1_PU_noCuts->SetLineColor(3);
  h_trk_MVA1_PU_noCuts->SetMarkerColor(3);
  h_trk_MVA1_notHiggs_noCuts->Scale(1./h_trk_MVA1_notHiggs_noCuts->Integral());
  h_trk_MVA1_notHiggs_noCuts->SetLineColor(4);
  h_trk_MVA1_notHiggs_noCuts->SetMarkerColor(4);
  h_trk_MVA1_fake_noCuts->Draw("HIST");
  h_trk_MVA1_primary_noCuts->Draw("HIST,SAME");
  h_trk_MVA1_PU_noCuts->Draw("HIST,SAME");
  h_trk_MVA1_notHiggs_noCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_MVA1_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_MVA1_fake_noCuts,"Fake","l");
  l->AddEntry(h_trk_MVA1_PU_noCuts,"PU","l");
  l->AddEntry(h_trk_MVA1_notHiggs_noCuts,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_MVA1_noCuts.pdf");
  delete h_trk_MVA1_primary_noCuts;
  delete h_trk_MVA1_fake_noCuts;
  delete h_trk_MVA1_PU_noCuts;
  delete h_trk_MVA1_notHiggs_noCuts;

  h_trk_MVA1_primary_noCuts_H->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_MVA1_primary_noCuts_H);
  h_trk_MVA1_np_noCuts_H->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_MVA1_np_noCuts_H);
  h_trk_MVA1_primary_noCuts_H->Scale(1./h_trk_MVA1_primary_noCuts_H->Integral());
  h_trk_MVA1_primary_noCuts_H->SetLineColor(1);
  h_trk_MVA1_primary_noCuts_H->SetMarkerColor(1);
  h_trk_MVA1_np_noCuts_H->Scale(1./h_trk_MVA1_np_noCuts_H->Integral());
  h_trk_MVA1_np_noCuts_H->SetLineColor(2);
  h_trk_MVA1_np_noCuts_H->SetMarkerColor(2);
  h_trk_MVA1_np_noCuts_H->Draw("HIST");
  h_trk_MVA1_primary_noCuts_H->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_MVA1_primary_noCuts_H,"Primary","l");
  l->AddEntry(h_trk_MVA1_np_noCuts_H,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_MVA1_noCuts_H.pdf");
  delete h_trk_MVA1_primary_noCuts_H;
  delete h_trk_MVA1_np_noCuts_H;

  h_trk_MVA1_primary_noCuts_L->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_MVA1_primary_noCuts_L);
  h_trk_MVA1_np_noCuts_L->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_MVA1_np_noCuts_L);
  h_trk_MVA1_primary_noCuts_L->Scale(1./h_trk_MVA1_primary_noCuts_L->Integral());
  h_trk_MVA1_primary_noCuts_L->SetLineColor(1);
  h_trk_MVA1_primary_noCuts_L->SetMarkerColor(1);
  h_trk_MVA1_np_noCuts_L->Scale(1./h_trk_MVA1_np_noCuts_L->Integral());
  h_trk_MVA1_np_noCuts_L->SetLineColor(2);
  h_trk_MVA1_np_noCuts_L->SetMarkerColor(2);
  h_trk_MVA1_np_noCuts_L->Draw("HIST");
  h_trk_MVA1_primary_noCuts_L->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_MVA1_primary_noCuts_L,"Primary","l");
  l->AddEntry(h_trk_MVA1_np_noCuts_L,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_MVA1_noCuts_L.pdf");
  delete h_trk_MVA1_primary_noCuts_L;
  delete h_trk_MVA1_np_noCuts_L;

  h_trk_MVA1_primary_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_MVA1_primary_qualCuts);
  h_trk_MVA1_np_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_MVA1_np_qualCuts);
  h_trk_MVA1_primary_qualCuts->Scale(1./h_trk_MVA1_primary_qualCuts->Integral());
  h_trk_MVA1_primary_qualCuts->SetLineColor(1);
  h_trk_MVA1_primary_qualCuts->SetMarkerColor(1);
  TH1F *h_trk_MVA1_np_qualCutsNorm = (TH1F*)h_trk_MVA1_np_qualCuts->Clone(); 
  h_trk_MVA1_np_qualCutsNorm->Scale(1./h_trk_MVA1_np_qualCutsNorm->Integral());
  h_trk_MVA1_np_qualCutsNorm->SetLineColor(2);
  h_trk_MVA1_np_qualCutsNorm->SetMarkerColor(2);
  h_trk_MVA1_np_qualCutsNorm->Draw("HIST");
  h_trk_MVA1_primary_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_MVA1_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_MVA1_np_qualCutsNorm,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_MVA1_qualCuts.pdf");
  delete h_trk_MVA1_np_qualCutsNorm;

  h_trk_MVA1_primary_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_MVA1_primary_allCuts);
  h_trk_MVA1_np_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_MVA1_np_allCuts);
  h_trk_MVA1_primary_allCuts->Scale(1./h_trk_MVA1_primary_allCuts->Integral());
  h_trk_MVA1_primary_allCuts->SetLineColor(1);
  h_trk_MVA1_primary_allCuts->SetMarkerColor(1);
  TH1F *h_trk_MVA1_np_allCutsNorm = (TH1F*)h_trk_MVA1_np_allCuts->Clone(); 
  h_trk_MVA1_np_allCutsNorm->Scale(1./h_trk_MVA1_np_allCutsNorm->Integral());
  h_trk_MVA1_np_allCutsNorm->SetLineColor(2);
  h_trk_MVA1_np_allCutsNorm->SetMarkerColor(2);
  raiseMax(h_trk_MVA1_primary_allCuts,h_trk_MVA1_np_allCutsNorm);
  h_trk_MVA1_primary_allCuts->Draw("HIST");
  h_trk_MVA1_np_allCutsNorm->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_MVA1_primary_allCuts,"Primary","l");
  l->AddEntry(h_trk_MVA1_np_allCutsNorm,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_MVA1_allCuts.pdf");
  delete h_trk_MVA1_np_allCutsNorm;

  removeFlows(h_trk_MVA2_primary_allCuts);
  h_trk_MVA2_primary_allCuts->Scale(1./h_trk_MVA2_primary_allCuts->Integral());
  h_trk_MVA2_primary_allCuts->SetLineColor(1);
  h_trk_MVA2_primary_allCuts->SetMarkerColor(1);
  h_trk_MVA2_primary_allCuts->SetStats(0);
  removeFlows(h_trk_MVA2_np_allCuts);
  h_trk_MVA2_np_allCuts->Scale(1./h_trk_MVA2_np_allCuts->Integral());
  h_trk_MVA2_np_allCuts->SetLineColor(2);
  h_trk_MVA2_np_allCuts->SetMarkerColor(2);
  h_trk_MVA2_np_allCuts->SetStats(0);
  raiseMax(h_trk_MVA2_np_allCuts,h_trk_MVA2_primary_allCuts);
  h_trk_MVA2_primary_allCuts->Draw("HIST");
  h_trk_MVA2_np_allCuts->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_MVA2_primary_allCuts,"Primary","l");
  l->AddEntry(h_trk_MVA2_np_allCuts,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_MVA2_allCuts.pdf");
  delete h_trk_MVA2_primary_allCuts;
  delete h_trk_MVA2_np_allCuts;

  removeFlows(h_trk_MVA1_primary_allCuts_P);
  h_trk_MVA1_primary_allCuts_P->Scale(1./h_trk_MVA1_primary_allCuts_P->Integral());
  h_trk_MVA1_primary_allCuts_P->SetLineColor(1);
  h_trk_MVA1_primary_allCuts_P->SetMarkerColor(1);
  h_trk_MVA1_primary_allCuts_P->SetStats(0);
  removeFlows(h_trk_MVA1_np_allCuts_P);
  h_trk_MVA1_np_allCuts_P->Scale(1./h_trk_MVA1_np_allCuts_P->Integral());
  h_trk_MVA1_np_allCuts_P->SetLineColor(2);
  h_trk_MVA1_np_allCuts_P->SetMarkerColor(2);
  h_trk_MVA1_np_allCuts_P->SetStats(0);
  raiseMax(h_trk_MVA1_np_allCuts_P,h_trk_MVA1_primary_allCuts_P);
  h_trk_MVA1_np_allCuts_P->Draw("HIST");
  h_trk_MVA1_primary_allCuts_P->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_MVA1_primary_allCuts_P,"Primary","l");
  l->AddEntry(h_trk_MVA1_np_allCuts_P,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_MVA1_allCuts_P.pdf");
  delete h_trk_MVA1_primary_allCuts_P;
  delete h_trk_MVA1_np_allCuts_P;

  removeFlows(h_trk_MVA1_primary_allCuts_D);
  h_trk_MVA1_primary_allCuts_D->Scale(1./h_trk_MVA1_primary_allCuts_D->Integral());
  h_trk_MVA1_primary_allCuts_D->SetLineColor(1);
  h_trk_MVA1_primary_allCuts_D->SetMarkerColor(1);
  h_trk_MVA1_primary_allCuts_D->SetStats(0);
  removeFlows(h_trk_MVA1_np_allCuts_D);
  h_trk_MVA1_np_allCuts_D->Scale(1./h_trk_MVA1_np_allCuts_D->Integral());
  h_trk_MVA1_np_allCuts_D->SetLineColor(2);
  h_trk_MVA1_np_allCuts_D->SetMarkerColor(2);
  h_trk_MVA1_np_allCuts_D->SetStats(0);
  raiseMax(h_trk_MVA1_np_allCuts_D,h_trk_MVA1_primary_allCuts_D);
  h_trk_MVA1_np_allCuts_D->Draw("HIST");
  h_trk_MVA1_primary_allCuts_D->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_MVA1_primary_allCuts_D,"Primary","l");
  l->AddEntry(h_trk_MVA1_np_allCuts_D,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_MVA1_allCuts_D.pdf");
  delete h_trk_MVA1_primary_allCuts_D;
  delete h_trk_MVA1_np_allCuts_D;

  removeFlows(h_trk_MVA2_primary_allCuts_P);
  h_trk_MVA2_primary_allCuts_P->Scale(1./h_trk_MVA2_primary_allCuts_P->Integral());
  h_trk_MVA2_primary_allCuts_P->SetLineColor(1);
  h_trk_MVA2_primary_allCuts_P->SetMarkerColor(1);
  h_trk_MVA2_primary_allCuts_P->SetStats(0);
  removeFlows(h_trk_MVA2_np_allCuts_P);
  h_trk_MVA2_np_allCuts_P->Scale(1./h_trk_MVA2_np_allCuts_P->Integral());
  h_trk_MVA2_np_allCuts_P->SetLineColor(2);
  h_trk_MVA2_np_allCuts_P->SetMarkerColor(2);
  h_trk_MVA2_np_allCuts_P->SetStats(0);
  raiseMax(h_trk_MVA2_np_allCuts_P,h_trk_MVA2_primary_allCuts_P);
  h_trk_MVA2_np_allCuts_P->Draw("HIST");
  h_trk_MVA2_primary_allCuts_P->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_MVA2_primary_allCuts_P,"Primary","l");
  l->AddEntry(h_trk_MVA2_np_allCuts_P,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_MVA2_allCuts_P.pdf");
  delete h_trk_MVA2_primary_allCuts_P;
  delete h_trk_MVA2_np_allCuts_P;

  removeFlows(h_trk_MVA2_primary_allCuts_D);
  h_trk_MVA2_primary_allCuts_D->Scale(1./h_trk_MVA2_primary_allCuts_D->Integral());
  h_trk_MVA2_primary_allCuts_D->SetLineColor(1);
  h_trk_MVA2_primary_allCuts_D->SetMarkerColor(1);
  h_trk_MVA2_primary_allCuts_D->SetStats(0);
  removeFlows(h_trk_MVA2_np_allCuts_D);
  h_trk_MVA2_np_allCuts_D->Scale(1./h_trk_MVA2_np_allCuts_D->Integral());
  h_trk_MVA2_np_allCuts_D->SetLineColor(2);
  h_trk_MVA2_np_allCuts_D->SetMarkerColor(2);
  h_trk_MVA2_np_allCuts_D->SetStats(0);
  raiseMax(h_trk_MVA2_np_allCuts_D,h_trk_MVA2_primary_allCuts_D);
  h_trk_MVA2_np_allCuts_D->Draw("HIST");
  h_trk_MVA2_primary_allCuts_D->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_MVA2_primary_allCuts_D,"Primary","l");
  l->AddEntry(h_trk_MVA2_np_allCuts_D,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_MVA2_allCuts_D.pdf");
  delete h_trk_MVA2_primary_allCuts_D;
  delete h_trk_MVA2_np_allCuts_D;

  TH1F *h_trk_MVA1_fake_qualCutsNorm = (TH1F*)h_trk_MVA1_fake_qualCuts->Clone(); 
  h_trk_MVA1_fake_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_MVA1_fake_qualCutsNorm);
  TH1F *h_trk_MVA1_PU_qualCutsNorm = (TH1F*)h_trk_MVA1_PU_qualCuts->Clone(); 
  h_trk_MVA1_PU_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_MVA1_PU_qualCutsNorm);
  TH1F *h_trk_MVA1_notHiggs_qualCutsNorm = (TH1F*)h_trk_MVA1_notHiggs_qualCuts->Clone(); 
  h_trk_MVA1_notHiggs_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_MVA1_notHiggs_qualCutsNorm);
  h_trk_MVA1_fake_qualCutsNorm->Scale(1./h_trk_MVA1_np_qualCuts->Integral());
  h_trk_MVA1_fake_qualCutsNorm->SetLineColor(2);
  h_trk_MVA1_fake_qualCutsNorm->SetMarkerColor(2);
  h_trk_MVA1_PU_qualCutsNorm->Scale(1./h_trk_MVA1_np_qualCuts->Integral());
  h_trk_MVA1_PU_qualCutsNorm->SetLineColor(3);
  h_trk_MVA1_PU_qualCutsNorm->SetMarkerColor(3);
  h_trk_MVA1_notHiggs_qualCutsNorm->Scale(1./h_trk_MVA1_np_qualCuts->Integral());
  h_trk_MVA1_notHiggs_qualCutsNorm->SetLineColor(4);
  h_trk_MVA1_notHiggs_qualCutsNorm->SetMarkerColor(4);
  auto hs_MVA1_qualCuts = new THStack("hs_MVA1_qualCuts","Stacked BG histograms");
  hs_MVA1_qualCuts->Add(h_trk_MVA1_fake_qualCutsNorm);
  hs_MVA1_qualCuts->Add(h_trk_MVA1_PU_qualCutsNorm);
  hs_MVA1_qualCuts->Add(h_trk_MVA1_notHiggs_qualCutsNorm);
  hs_MVA1_qualCuts->Draw("HIST");
  h_trk_MVA1_primary_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_MVA1_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_MVA1_fake_qualCutsNorm,"Fake","l");
  l->AddEntry(h_trk_MVA1_PU_qualCutsNorm,"PU","l");
  l->AddEntry(h_trk_MVA1_notHiggs_qualCutsNorm,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_MVA1_qualCutsNorm.pdf");
  delete h_trk_MVA1_fake_qualCutsNorm;
  delete h_trk_MVA1_PU_qualCutsNorm;
  delete h_trk_MVA1_notHiggs_qualCutsNorm;
  delete h_trk_MVA1_np_qualCuts;
  delete hs_MVA1_qualCuts;

  h_trk_MVA1_fake_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_MVA1_fake_qualCuts);
  h_trk_MVA1_PU_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_MVA1_PU_qualCuts);
  h_trk_MVA1_notHiggs_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_MVA1_notHiggs_qualCuts);
  h_trk_MVA1_fake_qualCuts->Scale(1./h_trk_MVA1_fake_qualCuts->Integral());
  h_trk_MVA1_fake_qualCuts->SetLineColor(2);
  h_trk_MVA1_fake_qualCuts->SetMarkerColor(2);
  h_trk_MVA1_PU_qualCuts->Scale(1./h_trk_MVA1_PU_qualCuts->Integral());
  h_trk_MVA1_PU_qualCuts->SetLineColor(3);
  h_trk_MVA1_PU_qualCuts->SetMarkerColor(3);
  h_trk_MVA1_notHiggs_qualCuts->Scale(1./h_trk_MVA1_notHiggs_qualCuts->Integral());
  h_trk_MVA1_notHiggs_qualCuts->SetLineColor(4);
  h_trk_MVA1_notHiggs_qualCuts->SetMarkerColor(4);
  h_trk_MVA1_fake_qualCuts->Draw("HIST");
  h_trk_MVA1_primary_qualCuts->Draw("HIST,SAME");
  h_trk_MVA1_PU_qualCuts->Draw("HIST,SAME");
  h_trk_MVA1_notHiggs_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_MVA1_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_MVA1_fake_qualCuts,"Fake","l");
  l->AddEntry(h_trk_MVA1_PU_qualCuts,"PU","l");
  l->AddEntry(h_trk_MVA1_notHiggs_qualCuts,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_MVA1_qualCuts.pdf");
  delete h_trk_MVA1_primary_qualCuts;
  delete h_trk_MVA1_fake_qualCuts;
  delete h_trk_MVA1_PU_qualCuts;
  delete h_trk_MVA1_notHiggs_qualCuts;

  h_trk_chi2rzdof->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof);
  h_trk_chi2rzdof->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_chi2rzdof->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_chi2rzdof->GetName() + ".pdf");
  delete h_trk_chi2rzdof;

  h_trk_chi2rzdof_primary->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_primary);
  h_trk_chi2rzdof_primary->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_chi2rzdof_primary->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_chi2rzdof_primary->GetName() + ".pdf");
  delete h_trk_chi2rzdof_primary;

  h_trk_chi2rzdof_primary_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_primary_noCuts);
  h_trk_chi2rzdof_primary_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_chi2rzdof_primary_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_chi2rzdof_primary_noCuts->GetName() + ".pdf");

  h_trk_chi2rzdof_np->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_np);
  h_trk_chi2rzdof_np->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_chi2rzdof_np->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_chi2rzdof_np->GetName() + ".pdf");
  delete h_trk_chi2rzdof_np;

  h_trk_chi2rzdof_np_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_np_noCuts);
  h_trk_chi2rzdof_np_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_chi2rzdof_np_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_chi2rzdof_np_noCuts->GetName() + ".pdf");
   
  h_trk_chi2rzdof_primary_noCuts->Scale(1./h_trk_chi2rzdof_primary_noCuts->Integral());
  h_trk_chi2rzdof_primary_noCuts->SetLineColor(1);
  h_trk_chi2rzdof_primary_noCuts->SetMarkerColor(1);
  TH1F *h_trk_chi2rzdof_np_noCutsNorm = (TH1F*)h_trk_chi2rzdof_np_noCuts->Clone(); 
  h_trk_chi2rzdof_np_noCutsNorm->Scale(1./h_trk_chi2rzdof_np_noCutsNorm->Integral());
  h_trk_chi2rzdof_np_noCutsNorm->SetLineColor(2);
  h_trk_chi2rzdof_np_noCutsNorm->SetMarkerColor(2);
  raiseMax(h_trk_chi2rzdof_primary_noCuts,h_trk_chi2rzdof_np_noCutsNorm);
  h_trk_chi2rzdof_primary_noCuts->SetStats(0);
  h_trk_chi2rzdof_np_noCutsNorm->SetStats(0);
  h_trk_chi2rzdof_primary_noCuts->Draw("HIST");
  h_trk_chi2rzdof_np_noCutsNorm->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_chi2rzdof_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_chi2rzdof_np_noCutsNorm,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_chi2rzdof_noCuts.pdf");
  delete h_trk_chi2rzdof_np_noCutsNorm;

  h_trk_chi2rzdof_primary_noCuts_zoomOut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_primary_noCuts_zoomOut);
  h_trk_chi2rzdof_np_noCuts_zoomOut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_np_noCuts_zoomOut);
  h_trk_chi2rzdof_primary_noCuts_zoomOut->Scale(1./h_trk_chi2rzdof_primary_noCuts_zoomOut->Integral());
  h_trk_chi2rzdof_primary_noCuts_zoomOut->SetLineColor(1);
  h_trk_chi2rzdof_primary_noCuts_zoomOut->SetMarkerColor(1);
  h_trk_chi2rzdof_np_noCuts_zoomOut->Scale(1./h_trk_chi2rzdof_np_noCuts_zoomOut->Integral());
  h_trk_chi2rzdof_np_noCuts_zoomOut->SetLineColor(2);
  h_trk_chi2rzdof_np_noCuts_zoomOut->SetMarkerColor(2);
  raiseMax(h_trk_chi2rzdof_np_noCuts_zoomOut,h_trk_chi2rzdof_primary_noCuts_zoomOut);
  h_trk_chi2rzdof_np_noCuts_zoomOut->SetStats(0);
  h_trk_chi2rzdof_primary_noCuts_zoomOut->SetStats(0);
  h_trk_chi2rzdof_np_noCuts_zoomOut->Draw("HIST");
  h_trk_chi2rzdof_primary_noCuts_zoomOut->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->SetFillColor(0);
  l->SetLineColor(0);
  l->SetTextSize(0.04);
  l->AddEntry(h_trk_chi2rzdof_primary_noCuts_zoomOut,"Primary","l");
  l->AddEntry(h_trk_chi2rzdof_np_noCuts_zoomOut,"NP","l");
  l->SetTextFont(42);
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_chi2rzdof_noCuts_zoomOut.pdf");
  delete h_trk_chi2rzdof_primary_noCuts_zoomOut;
  delete h_trk_chi2rzdof_np_noCuts_zoomOut;

  TH1F *h_trk_chi2rzdof_fake_noCutsNorm = (TH1F*)h_trk_chi2rzdof_fake_noCuts->Clone(); 
  h_trk_chi2rzdof_fake_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_fake_noCutsNorm);
  TH1F *h_trk_chi2rzdof_PU_noCutsNorm = (TH1F*)h_trk_chi2rzdof_PU_noCuts->Clone(); 
  h_trk_chi2rzdof_PU_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_PU_noCutsNorm);
  TH1F *h_trk_chi2rzdof_notHiggs_noCutsNorm = (TH1F*)h_trk_chi2rzdof_notHiggs_noCuts->Clone(); 
  h_trk_chi2rzdof_notHiggs_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_notHiggs_noCutsNorm);
  h_trk_chi2rzdof_fake_noCutsNorm->Scale(1./h_trk_chi2rzdof_np_noCuts->Integral());
  h_trk_chi2rzdof_fake_noCutsNorm->SetLineColor(2);
  h_trk_chi2rzdof_fake_noCutsNorm->SetMarkerColor(2);
  h_trk_chi2rzdof_PU_noCutsNorm->Scale(1./h_trk_chi2rzdof_np_noCuts->Integral());
  h_trk_chi2rzdof_PU_noCutsNorm->SetLineColor(3);
  h_trk_chi2rzdof_PU_noCutsNorm->SetMarkerColor(3);
  h_trk_chi2rzdof_notHiggs_noCutsNorm->Scale(1./h_trk_chi2rzdof_np_noCuts->Integral());
  h_trk_chi2rzdof_notHiggs_noCutsNorm->SetLineColor(4);
  h_trk_chi2rzdof_notHiggs_noCutsNorm->SetMarkerColor(4);
  auto hs_chi2rzdof_noCuts = new THStack("hs_chi2rzdof_noCuts","Stacked BG histograms");
  hs_chi2rzdof_noCuts->Add(h_trk_chi2rzdof_fake_noCutsNorm);
  hs_chi2rzdof_noCuts->Add(h_trk_chi2rzdof_PU_noCutsNorm);
  hs_chi2rzdof_noCuts->Add(h_trk_chi2rzdof_notHiggs_noCutsNorm);
  hs_chi2rzdof_noCuts->Draw("HIST");
  h_trk_chi2rzdof_primary_noCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_chi2rzdof_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_chi2rzdof_fake_noCutsNorm,"Fake","l");
  l->AddEntry(h_trk_chi2rzdof_PU_noCutsNorm,"PU","l");
  l->AddEntry(h_trk_chi2rzdof_notHiggs_noCutsNorm,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_chi2rzdof_noCutsNorm.pdf");
  delete h_trk_chi2rzdof_fake_noCutsNorm;
  delete h_trk_chi2rzdof_PU_noCutsNorm;
  delete h_trk_chi2rzdof_notHiggs_noCutsNorm;
  delete h_trk_chi2rzdof_np_noCuts;
  delete hs_chi2rzdof_noCuts;

  h_trk_chi2rzdof_fake_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_fake_noCuts);
  h_trk_chi2rzdof_PU_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_PU_noCuts);
  h_trk_chi2rzdof_notHiggs_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_notHiggs_noCuts);
  h_trk_chi2rzdof_fake_noCuts->Scale(1./h_trk_chi2rzdof_fake_noCuts->Integral());
  h_trk_chi2rzdof_fake_noCuts->SetLineColor(2);
  h_trk_chi2rzdof_fake_noCuts->SetMarkerColor(2);
  h_trk_chi2rzdof_PU_noCuts->Scale(1./h_trk_chi2rzdof_PU_noCuts->Integral());
  h_trk_chi2rzdof_PU_noCuts->SetLineColor(3);
  h_trk_chi2rzdof_PU_noCuts->SetMarkerColor(3);
  h_trk_chi2rzdof_notHiggs_noCuts->Scale(1./h_trk_chi2rzdof_notHiggs_noCuts->Integral());
  h_trk_chi2rzdof_notHiggs_noCuts->SetLineColor(4);
  h_trk_chi2rzdof_notHiggs_noCuts->SetMarkerColor(4);
  h_trk_chi2rzdof_fake_noCuts->Draw("HIST");
  h_trk_chi2rzdof_primary_noCuts->Draw("HIST,SAME");
  h_trk_chi2rzdof_PU_noCuts->Draw("HIST,SAME");
  h_trk_chi2rzdof_notHiggs_noCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_chi2rzdof_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_chi2rzdof_fake_noCuts,"Fake","l");
  l->AddEntry(h_trk_chi2rzdof_PU_noCuts,"PU","l");
  l->AddEntry(h_trk_chi2rzdof_notHiggs_noCuts,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_chi2rzdof_noCuts.pdf");
  delete h_trk_chi2rzdof_primary_noCuts;
  delete h_trk_chi2rzdof_fake_noCuts;
  delete h_trk_chi2rzdof_PU_noCuts;
  delete h_trk_chi2rzdof_notHiggs_noCuts;

  h_trk_chi2rzdof_primary_noCuts_H->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_primary_noCuts_H);
  h_trk_chi2rzdof_np_noCuts_H->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_np_noCuts_H);
  h_trk_chi2rzdof_primary_noCuts_H->Scale(1./h_trk_chi2rzdof_primary_noCuts_H->Integral());
  h_trk_chi2rzdof_primary_noCuts_H->SetLineColor(1);
  h_trk_chi2rzdof_primary_noCuts_H->SetMarkerColor(1);
  h_trk_chi2rzdof_np_noCuts_H->Scale(1./h_trk_chi2rzdof_np_noCuts_H->Integral());
  h_trk_chi2rzdof_np_noCuts_H->SetLineColor(2);
  h_trk_chi2rzdof_np_noCuts_H->SetMarkerColor(2);
  h_trk_chi2rzdof_np_noCuts_H->Draw("HIST");
  h_trk_chi2rzdof_primary_noCuts_H->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_chi2rzdof_primary_noCuts_H,"Primary","l");
  l->AddEntry(h_trk_chi2rzdof_np_noCuts_H,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_chi2rzdof_noCuts_H.pdf");
  delete h_trk_chi2rzdof_primary_noCuts_H;
  delete h_trk_chi2rzdof_np_noCuts_H;

  h_trk_chi2rzdof_primary_noCuts_L->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_primary_noCuts_L);
  h_trk_chi2rzdof_np_noCuts_L->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_np_noCuts_L);
  h_trk_chi2rzdof_primary_noCuts_L->Scale(1./h_trk_chi2rzdof_primary_noCuts_L->Integral());
  h_trk_chi2rzdof_primary_noCuts_L->SetLineColor(1);
  h_trk_chi2rzdof_primary_noCuts_L->SetMarkerColor(1);
  h_trk_chi2rzdof_np_noCuts_L->Scale(1./h_trk_chi2rzdof_np_noCuts_L->Integral());
  h_trk_chi2rzdof_np_noCuts_L->SetLineColor(2);
  h_trk_chi2rzdof_np_noCuts_L->SetMarkerColor(2);
  h_trk_chi2rzdof_np_noCuts_L->Draw("HIST");
  h_trk_chi2rzdof_primary_noCuts_L->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_chi2rzdof_primary_noCuts_L,"Primary","l");
  l->AddEntry(h_trk_chi2rzdof_np_noCuts_L,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_chi2rzdof_noCuts_L.pdf");
  delete h_trk_chi2rzdof_primary_noCuts_L;
  delete h_trk_chi2rzdof_np_noCuts_L;

  h_trk_chi2rzdof_primary_noCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_primary_noCuts_barrel);
  h_trk_chi2rzdof_np_noCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_np_noCuts_barrel);
  h_trk_chi2rzdof_primary_noCuts_barrel->Scale(1./h_trk_chi2rzdof_primary_noCuts_barrel->Integral());
  h_trk_chi2rzdof_primary_noCuts_barrel->SetLineColor(1);
  h_trk_chi2rzdof_primary_noCuts_barrel->SetMarkerColor(1);
  h_trk_chi2rzdof_np_noCuts_barrel->Scale(1./h_trk_chi2rzdof_np_noCuts_barrel->Integral());
  h_trk_chi2rzdof_np_noCuts_barrel->SetLineColor(2);
  h_trk_chi2rzdof_np_noCuts_barrel->SetMarkerColor(2);
  h_trk_chi2rzdof_np_noCuts_barrel->Draw("HIST");
  h_trk_chi2rzdof_primary_noCuts_barrel->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_chi2rzdof_primary_noCuts_barrel,"Primary","l");
  l->AddEntry(h_trk_chi2rzdof_np_noCuts_barrel,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_chi2rzdof_noCuts_barrel.pdf");
  delete h_trk_chi2rzdof_primary_noCuts_barrel;
  delete h_trk_chi2rzdof_np_noCuts_barrel;

  h_trk_chi2rzdof_primary_noCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_primary_noCuts_disk);
  h_trk_chi2rzdof_np_noCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_np_noCuts_disk);
  h_trk_chi2rzdof_primary_noCuts_disk->Scale(1./h_trk_chi2rzdof_primary_noCuts_disk->Integral());
  h_trk_chi2rzdof_primary_noCuts_disk->SetLineColor(1);
  h_trk_chi2rzdof_primary_noCuts_disk->SetMarkerColor(1);
  h_trk_chi2rzdof_np_noCuts_disk->Scale(1./h_trk_chi2rzdof_np_noCuts_disk->Integral());
  h_trk_chi2rzdof_np_noCuts_disk->SetLineColor(2);
  h_trk_chi2rzdof_np_noCuts_disk->SetMarkerColor(2);
  h_trk_chi2rzdof_np_noCuts_disk->Draw("HIST");
  h_trk_chi2rzdof_primary_noCuts_disk->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_chi2rzdof_primary_noCuts_disk,"Primary","l");
  l->AddEntry(h_trk_chi2rzdof_np_noCuts_disk,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_chi2rzdof_noCuts_disk.pdf");
  delete h_trk_chi2rzdof_primary_noCuts_disk;
  delete h_trk_chi2rzdof_np_noCuts_disk;

  h_trk_chi2rzdof_primary_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_primary_qualCuts);
  h_trk_chi2rzdof_np_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_np_qualCuts);
  h_trk_chi2rzdof_primary_qualCuts->Scale(1./h_trk_chi2rzdof_primary_qualCuts->Integral());
  h_trk_chi2rzdof_primary_qualCuts->SetLineColor(1);
  h_trk_chi2rzdof_primary_qualCuts->SetMarkerColor(1);
  TH1F *h_trk_chi2rzdof_np_qualCutsNorm = (TH1F*)h_trk_chi2rzdof_np_qualCuts->Clone(); 
  h_trk_chi2rzdof_np_qualCutsNorm->Scale(1./h_trk_chi2rzdof_np_qualCutsNorm->Integral());
  h_trk_chi2rzdof_np_qualCutsNorm->SetLineColor(2);
  h_trk_chi2rzdof_np_qualCutsNorm->SetMarkerColor(2);
  h_trk_chi2rzdof_np_qualCutsNorm->Draw("HIST");
  h_trk_chi2rzdof_primary_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_chi2rzdof_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_chi2rzdof_np_qualCutsNorm,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_chi2rzdof_qualCuts.pdf");
  delete h_trk_chi2rzdof_np_qualCutsNorm;

  h_trk_chi2rzdof_primary_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_primary_allCuts);
  h_trk_chi2rzdof_np_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_np_allCuts);
  h_trk_chi2rzdof_primary_allCuts->Scale(1./h_trk_chi2rzdof_primary_allCuts->Integral());
  h_trk_chi2rzdof_primary_allCuts->SetLineColor(1);
  h_trk_chi2rzdof_primary_allCuts->SetMarkerColor(1);
  TH1F *h_trk_chi2rzdof_np_allCutsNorm = (TH1F*)h_trk_chi2rzdof_np_allCuts->Clone(); 
  h_trk_chi2rzdof_np_allCutsNorm->Scale(1./h_trk_chi2rzdof_np_allCutsNorm->Integral());
  h_trk_chi2rzdof_np_allCutsNorm->SetLineColor(2);
  h_trk_chi2rzdof_np_allCutsNorm->SetMarkerColor(2);
  raiseMax(h_trk_chi2rzdof_primary_allCuts,h_trk_chi2rzdof_np_allCutsNorm);
  h_trk_chi2rzdof_np_allCutsNorm->SetStats(0);
  h_trk_chi2rzdof_primary_allCuts->SetStats(0);
  h_trk_chi2rzdof_primary_allCuts->Draw("HIST");
  h_trk_chi2rzdof_np_allCutsNorm->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_chi2rzdof_primary_allCuts,"Primary","l");
  l->AddEntry(h_trk_chi2rzdof_np_allCutsNorm,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_chi2rzdof_allCuts.pdf");
  delete h_trk_chi2rzdof_np_allCutsNorm;

  h_trk_chi2rzdof_primary_allCuts_zoomOut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_primary_allCuts_zoomOut);
  h_trk_chi2rzdof_np_allCuts_zoomOut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_np_allCuts_zoomOut);
  h_trk_chi2rzdof_primary_allCuts_zoomOut->Scale(1./h_trk_chi2rzdof_primary_allCuts_zoomOut->Integral());
  h_trk_chi2rzdof_primary_allCuts_zoomOut->SetLineColor(1);
  h_trk_chi2rzdof_primary_allCuts_zoomOut->SetMarkerColor(1);
  h_trk_chi2rzdof_np_allCuts_zoomOut->Scale(1./h_trk_chi2rzdof_np_allCuts_zoomOut->Integral());
  h_trk_chi2rzdof_np_allCuts_zoomOut->SetLineColor(2);
  h_trk_chi2rzdof_np_allCuts_zoomOut->SetMarkerColor(2);
  raiseMax(h_trk_chi2rzdof_np_allCuts_zoomOut,h_trk_chi2rzdof_primary_allCuts_zoomOut);
  h_trk_chi2rzdof_np_allCuts_zoomOut->SetStats(0);
  h_trk_chi2rzdof_primary_allCuts_zoomOut->SetStats(0);
  h_trk_chi2rzdof_np_allCuts_zoomOut->Draw("HIST");
  h_trk_chi2rzdof_primary_allCuts_zoomOut->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->SetFillColor(0);
  l->SetLineColor(0);
  l->SetTextSize(0.04);
  l->AddEntry(h_trk_chi2rzdof_primary_allCuts_zoomOut,"Primary","l");
  l->AddEntry(h_trk_chi2rzdof_np_allCuts_zoomOut,"NP","l");
  l->SetTextFont(42);
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_chi2rzdof_allCuts_zoomOut.pdf");
  delete h_trk_chi2rzdof_primary_allCuts_zoomOut;
  delete h_trk_chi2rzdof_np_allCuts_zoomOut;

  h_trk_chi2rzdof_primary_allCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_primary_allCuts_barrel);
  h_trk_chi2rzdof_np_allCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_np_allCuts_barrel);
  h_trk_chi2rzdof_primary_allCuts_barrel->Scale(1./h_trk_chi2rzdof_primary_allCuts_barrel->Integral());
  h_trk_chi2rzdof_primary_allCuts_barrel->SetLineColor(1);
  h_trk_chi2rzdof_primary_allCuts_barrel->SetMarkerColor(1);
  h_trk_chi2rzdof_np_allCuts_barrel->Scale(1./h_trk_chi2rzdof_np_allCuts_barrel->Integral());
  h_trk_chi2rzdof_np_allCuts_barrel->SetLineColor(2);
  h_trk_chi2rzdof_np_allCuts_barrel->SetMarkerColor(2);
  h_trk_chi2rzdof_np_allCuts_barrel->Draw("HIST");
  h_trk_chi2rzdof_primary_allCuts_barrel->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_chi2rzdof_primary_allCuts_barrel,"Primary","l");
  l->AddEntry(h_trk_chi2rzdof_np_allCuts_barrel,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_chi2rzdof_allCuts_barrel.pdf");
  delete h_trk_chi2rzdof_primary_allCuts_barrel;
  delete h_trk_chi2rzdof_np_allCuts_barrel;

  h_trk_chi2rzdof_primary_allCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_primary_allCuts_disk);
  h_trk_chi2rzdof_np_allCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_np_allCuts_disk);
  h_trk_chi2rzdof_primary_allCuts_disk->Scale(1./h_trk_chi2rzdof_primary_allCuts_disk->Integral());
  h_trk_chi2rzdof_primary_allCuts_disk->SetLineColor(1);
  h_trk_chi2rzdof_primary_allCuts_disk->SetMarkerColor(1);
  h_trk_chi2rzdof_np_allCuts_disk->Scale(1./h_trk_chi2rzdof_np_allCuts_disk->Integral());
  h_trk_chi2rzdof_np_allCuts_disk->SetLineColor(2);
  h_trk_chi2rzdof_np_allCuts_disk->SetMarkerColor(2);
  h_trk_chi2rzdof_np_allCuts_disk->Draw("HIST");
  h_trk_chi2rzdof_primary_allCuts_disk->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_chi2rzdof_primary_allCuts_disk,"Primary","l");
  l->AddEntry(h_trk_chi2rzdof_np_allCuts_disk,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_chi2rzdof_allCuts_disk.pdf");
  delete h_trk_chi2rzdof_primary_allCuts_disk;
  delete h_trk_chi2rzdof_np_allCuts_disk;

  TH1F *h_trk_chi2rzdof_fake_qualCutsNorm = (TH1F*)h_trk_chi2rzdof_fake_qualCuts->Clone(); 
  h_trk_chi2rzdof_fake_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_fake_qualCutsNorm);
  TH1F *h_trk_chi2rzdof_PU_qualCutsNorm = (TH1F*)h_trk_chi2rzdof_PU_qualCuts->Clone(); 
  h_trk_chi2rzdof_PU_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_PU_qualCutsNorm);
  TH1F *h_trk_chi2rzdof_notHiggs_qualCutsNorm = (TH1F*)h_trk_chi2rzdof_notHiggs_qualCuts->Clone(); 
  h_trk_chi2rzdof_notHiggs_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_notHiggs_qualCutsNorm);
  h_trk_chi2rzdof_fake_qualCutsNorm->Scale(1./h_trk_chi2rzdof_np_qualCuts->Integral());
  h_trk_chi2rzdof_fake_qualCutsNorm->SetLineColor(2);
  h_trk_chi2rzdof_fake_qualCutsNorm->SetMarkerColor(2);
  h_trk_chi2rzdof_PU_qualCutsNorm->Scale(1./h_trk_chi2rzdof_np_qualCuts->Integral());
  h_trk_chi2rzdof_PU_qualCutsNorm->SetLineColor(3);
  h_trk_chi2rzdof_PU_qualCutsNorm->SetMarkerColor(3);
  h_trk_chi2rzdof_notHiggs_qualCutsNorm->Scale(1./h_trk_chi2rzdof_np_qualCuts->Integral());
  h_trk_chi2rzdof_notHiggs_qualCutsNorm->SetLineColor(4);
  h_trk_chi2rzdof_notHiggs_qualCutsNorm->SetMarkerColor(4);
  auto hs_chi2rzdof_qualCuts = new THStack("hs_chi2rzdof_qualCuts","Stacked BG histograms");
  hs_chi2rzdof_qualCuts->Add(h_trk_chi2rzdof_fake_qualCutsNorm);
  hs_chi2rzdof_qualCuts->Add(h_trk_chi2rzdof_PU_qualCutsNorm);
  hs_chi2rzdof_qualCuts->Add(h_trk_chi2rzdof_notHiggs_qualCutsNorm);
  hs_chi2rzdof_qualCuts->Draw("HIST");
  h_trk_chi2rzdof_primary_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_chi2rzdof_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_chi2rzdof_fake_qualCutsNorm,"Fake","l");
  l->AddEntry(h_trk_chi2rzdof_PU_qualCutsNorm,"PU","l");
  l->AddEntry(h_trk_chi2rzdof_notHiggs_qualCutsNorm,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_chi2rzdof_qualCutsNorm.pdf");
  delete h_trk_chi2rzdof_fake_qualCutsNorm;
  delete h_trk_chi2rzdof_PU_qualCutsNorm;
  delete h_trk_chi2rzdof_notHiggs_qualCutsNorm;
  delete h_trk_chi2rzdof_np_qualCuts;
  delete hs_chi2rzdof_qualCuts;

  h_trk_chi2rzdof_fake_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_fake_qualCuts);
  h_trk_chi2rzdof_PU_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_PU_qualCuts);
  h_trk_chi2rzdof_notHiggs_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_notHiggs_qualCuts);
  h_trk_chi2rzdof_fake_qualCuts->Scale(1./h_trk_chi2rzdof_fake_qualCuts->Integral());
  h_trk_chi2rzdof_fake_qualCuts->SetLineColor(2);
  h_trk_chi2rzdof_fake_qualCuts->SetMarkerColor(2);
  h_trk_chi2rzdof_PU_qualCuts->Scale(1./h_trk_chi2rzdof_PU_qualCuts->Integral());
  h_trk_chi2rzdof_PU_qualCuts->SetLineColor(3);
  h_trk_chi2rzdof_PU_qualCuts->SetMarkerColor(3);
  h_trk_chi2rzdof_notHiggs_qualCuts->Scale(1./h_trk_chi2rzdof_notHiggs_qualCuts->Integral());
  h_trk_chi2rzdof_notHiggs_qualCuts->SetLineColor(4);
  h_trk_chi2rzdof_notHiggs_qualCuts->SetMarkerColor(4);
  h_trk_chi2rzdof_fake_qualCuts->Draw("HIST");
  h_trk_chi2rzdof_primary_qualCuts->Draw("HIST,SAME");
  h_trk_chi2rzdof_PU_qualCuts->Draw("HIST,SAME");
  h_trk_chi2rzdof_notHiggs_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_chi2rzdof_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_chi2rzdof_fake_qualCuts,"Fake","l");
  l->AddEntry(h_trk_chi2rzdof_PU_qualCuts,"PU","l");
  l->AddEntry(h_trk_chi2rzdof_notHiggs_qualCuts,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_chi2rzdof_qualCuts.pdf");
  delete h_trk_chi2rzdof_primary_qualCuts;
  delete h_trk_chi2rzdof_fake_qualCuts;
  delete h_trk_chi2rzdof_PU_qualCuts;
  delete h_trk_chi2rzdof_notHiggs_qualCuts;

  h_trk_bendchi2->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2);
  h_trk_bendchi2->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_bendchi2->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_bendchi2->GetName() + ".pdf");
  delete h_trk_bendchi2;

  h_trk_bendchi2_primary->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_primary);
  h_trk_bendchi2_primary->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_bendchi2_primary->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_bendchi2_primary->GetName() + ".pdf");
  delete h_trk_bendchi2_primary;

  h_trk_bendchi2_primary_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_primary_noCuts);
  h_trk_bendchi2_primary_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_bendchi2_primary_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_bendchi2_primary_noCuts->GetName() + ".pdf");

  h_trk_bendchi2_np->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_np);
  h_trk_bendchi2_np->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_bendchi2_np->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_bendchi2_np->GetName() + ".pdf");
  delete h_trk_bendchi2_np;

  h_trk_bendchi2_np_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_np_noCuts);
  h_trk_bendchi2_np_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_bendchi2_np_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_bendchi2_np_noCuts->GetName() + ".pdf");

  h_trk_bendchi2_primary_noCuts->Scale(1./h_trk_bendchi2_primary_noCuts->Integral());
  h_trk_bendchi2_primary_noCuts->SetLineColor(1);
  h_trk_bendchi2_primary_noCuts->SetMarkerColor(1);
  TH1F *h_trk_bendchi2_np_noCutsNorm = (TH1F*)h_trk_bendchi2_np_noCuts->Clone(); 
  h_trk_bendchi2_np_noCutsNorm->Scale(1./h_trk_bendchi2_np_noCutsNorm->Integral());
  h_trk_bendchi2_np_noCutsNorm->SetLineColor(2);
  h_trk_bendchi2_np_noCutsNorm->SetMarkerColor(2);
  raiseMax(h_trk_bendchi2_np_noCutsNorm,h_trk_bendchi2_primary_noCuts);
  h_trk_bendchi2_np_noCutsNorm->SetStats(0);
  h_trk_bendchi2_primary_noCuts->SetStats(0);
  h_trk_bendchi2_np_noCutsNorm->Draw("HIST");
  h_trk_bendchi2_primary_noCuts->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_bendchi2_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_bendchi2_np_noCutsNorm,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_bendchi2_noCuts.pdf");
  delete h_trk_bendchi2_np_noCutsNorm;

  h_trk_bendchi2_primary_noCuts_zoomOut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_primary_noCuts_zoomOut);
  h_trk_bendchi2_np_noCuts_zoomOut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_np_noCuts_zoomOut);
  h_trk_bendchi2_primary_noCuts_zoomOut->Scale(1./h_trk_bendchi2_primary_noCuts_zoomOut->Integral());
  h_trk_bendchi2_primary_noCuts_zoomOut->SetLineColor(1);
  h_trk_bendchi2_primary_noCuts_zoomOut->SetMarkerColor(1);
  h_trk_bendchi2_np_noCuts_zoomOut->Scale(1./h_trk_bendchi2_np_noCuts_zoomOut->Integral());
  h_trk_bendchi2_np_noCuts_zoomOut->SetLineColor(2);
  h_trk_bendchi2_np_noCuts_zoomOut->SetMarkerColor(2);
  raiseMax(h_trk_bendchi2_np_noCuts_zoomOut,h_trk_bendchi2_primary_noCuts_zoomOut);
  h_trk_bendchi2_np_noCuts_zoomOut->SetStats(0);
  h_trk_bendchi2_primary_noCuts_zoomOut->SetStats(0);
  h_trk_bendchi2_np_noCuts_zoomOut->Draw("HIST");
  h_trk_bendchi2_primary_noCuts_zoomOut->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->SetFillColor(0);
  l->SetLineColor(0);
  l->SetTextSize(0.04);
  l->AddEntry(h_trk_bendchi2_primary_noCuts_zoomOut,"Primary","l");
  l->AddEntry(h_trk_bendchi2_np_noCuts_zoomOut,"NP","l");
  l->SetTextFont(42);
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_bendchi2_noCuts_zoomOut.pdf");
  delete h_trk_bendchi2_primary_noCuts_zoomOut;
  delete h_trk_bendchi2_np_noCuts_zoomOut;

  TH1F *h_trk_bendchi2_fake_noCutsNorm = (TH1F*)h_trk_bendchi2_fake_noCuts->Clone(); 
  h_trk_bendchi2_fake_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_fake_noCutsNorm);
  TH1F *h_trk_bendchi2_PU_noCutsNorm = (TH1F*)h_trk_bendchi2_PU_noCuts->Clone(); 
  h_trk_bendchi2_PU_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_PU_noCutsNorm);
  TH1F *h_trk_bendchi2_notHiggs_noCutsNorm = (TH1F*)h_trk_bendchi2_notHiggs_noCuts->Clone(); 
  h_trk_bendchi2_notHiggs_noCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_notHiggs_noCutsNorm);
  h_trk_bendchi2_fake_noCutsNorm->Scale(1./h_trk_bendchi2_np_noCuts->Integral());
  h_trk_bendchi2_fake_noCutsNorm->SetLineColor(2);
  h_trk_bendchi2_fake_noCutsNorm->SetMarkerColor(2);
  h_trk_bendchi2_PU_noCutsNorm->Scale(1./h_trk_bendchi2_np_noCuts->Integral());
  h_trk_bendchi2_PU_noCutsNorm->SetLineColor(3);
  h_trk_bendchi2_PU_noCutsNorm->SetMarkerColor(3);
  h_trk_bendchi2_notHiggs_noCutsNorm->Scale(1./h_trk_bendchi2_np_noCuts->Integral());
  h_trk_bendchi2_notHiggs_noCutsNorm->SetLineColor(4);
  h_trk_bendchi2_notHiggs_noCutsNorm->SetMarkerColor(4);
  auto hs_bendchi2_noCuts = new THStack("hs_bendchi2_noCuts","Stacked BG histograms");
  hs_bendchi2_noCuts->Add(h_trk_bendchi2_fake_noCutsNorm);
  hs_bendchi2_noCuts->Add(h_trk_bendchi2_PU_noCutsNorm);
  hs_bendchi2_noCuts->Add(h_trk_bendchi2_notHiggs_noCutsNorm);
  hs_bendchi2_noCuts->Draw("HIST");
  h_trk_bendchi2_primary_noCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_bendchi2_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_bendchi2_fake_noCutsNorm,"Fake","l");
  l->AddEntry(h_trk_bendchi2_PU_noCutsNorm,"PU","l");
  l->AddEntry(h_trk_bendchi2_notHiggs_noCutsNorm,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_bendchi2_noCutsNorm.pdf");
  delete h_trk_bendchi2_fake_noCutsNorm;
  delete h_trk_bendchi2_PU_noCutsNorm;
  delete h_trk_bendchi2_notHiggs_noCutsNorm;
  delete h_trk_bendchi2_np_noCuts;
  delete hs_bendchi2_noCuts;

  h_trk_bendchi2_fake_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_fake_noCuts);
  h_trk_bendchi2_PU_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_PU_noCuts);
  h_trk_bendchi2_notHiggs_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_notHiggs_noCuts);
  h_trk_bendchi2_fake_noCuts->Scale(1./h_trk_bendchi2_fake_noCuts->Integral());
  h_trk_bendchi2_fake_noCuts->SetLineColor(2);
  h_trk_bendchi2_fake_noCuts->SetMarkerColor(2);
  h_trk_bendchi2_PU_noCuts->Scale(1./h_trk_bendchi2_PU_noCuts->Integral());
  h_trk_bendchi2_PU_noCuts->SetLineColor(3);
  h_trk_bendchi2_PU_noCuts->SetMarkerColor(3);
  h_trk_bendchi2_notHiggs_noCuts->Scale(1./h_trk_bendchi2_notHiggs_noCuts->Integral());
  h_trk_bendchi2_notHiggs_noCuts->SetLineColor(4);
  h_trk_bendchi2_notHiggs_noCuts->SetMarkerColor(4);
  h_trk_bendchi2_fake_noCuts->Draw("HIST");
  h_trk_bendchi2_primary_noCuts->Draw("HIST,SAME");
  h_trk_bendchi2_PU_noCuts->Draw("HIST,SAME");
  h_trk_bendchi2_notHiggs_noCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_bendchi2_primary_noCuts,"Primary","l");
  l->AddEntry(h_trk_bendchi2_fake_noCuts,"Fake","l");
  l->AddEntry(h_trk_bendchi2_PU_noCuts,"PU","l");
  l->AddEntry(h_trk_bendchi2_notHiggs_noCuts,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_bendchi2_noCuts.pdf");
  delete h_trk_bendchi2_primary_noCuts;
  delete h_trk_bendchi2_fake_noCuts;
  delete h_trk_bendchi2_PU_noCuts;
  delete h_trk_bendchi2_notHiggs_noCuts;

  h_trk_bendchi2_primary_noCuts_H->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_primary_noCuts_H);
  h_trk_bendchi2_np_noCuts_H->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_np_noCuts_H);
  h_trk_bendchi2_primary_noCuts_H->Scale(1./h_trk_bendchi2_primary_noCuts_H->Integral());
  h_trk_bendchi2_primary_noCuts_H->SetLineColor(1);
  h_trk_bendchi2_primary_noCuts_H->SetMarkerColor(1);
  h_trk_bendchi2_np_noCuts_H->Scale(1./h_trk_bendchi2_np_noCuts_H->Integral());
  h_trk_bendchi2_np_noCuts_H->SetLineColor(2);
  h_trk_bendchi2_np_noCuts_H->SetMarkerColor(2);
  h_trk_bendchi2_np_noCuts_H->Draw("HIST");
  h_trk_bendchi2_primary_noCuts_H->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_bendchi2_primary_noCuts_H,"Primary","l");
  l->AddEntry(h_trk_bendchi2_np_noCuts_H,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_bendchi2_noCuts_H.pdf");
  delete h_trk_bendchi2_primary_noCuts_H;
  delete h_trk_bendchi2_np_noCuts_H;

  h_trk_bendchi2_primary_noCuts_L->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_primary_noCuts_L);
  h_trk_bendchi2_np_noCuts_L->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_np_noCuts_L);
  h_trk_bendchi2_primary_noCuts_L->Scale(1./h_trk_bendchi2_primary_noCuts_L->Integral());
  h_trk_bendchi2_primary_noCuts_L->SetLineColor(1);
  h_trk_bendchi2_primary_noCuts_L->SetMarkerColor(1);
  h_trk_bendchi2_np_noCuts_L->Scale(1./h_trk_bendchi2_np_noCuts_L->Integral());
  h_trk_bendchi2_np_noCuts_L->SetLineColor(2);
  h_trk_bendchi2_np_noCuts_L->SetMarkerColor(2);
  h_trk_bendchi2_np_noCuts_L->Draw("HIST");
  h_trk_bendchi2_primary_noCuts_L->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_bendchi2_primary_noCuts_L,"Primary","l");
  l->AddEntry(h_trk_bendchi2_np_noCuts_L,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_bendchi2_noCuts_L.pdf");
  delete h_trk_bendchi2_primary_noCuts_L;
  delete h_trk_bendchi2_np_noCuts_L;

  h_trk_bendchi2_primary_noCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_primary_noCuts_barrel);
  h_trk_bendchi2_np_noCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_np_noCuts_barrel);
  h_trk_bendchi2_primary_noCuts_barrel->Scale(1./h_trk_bendchi2_primary_noCuts_barrel->Integral());
  h_trk_bendchi2_primary_noCuts_barrel->SetLineColor(1);
  h_trk_bendchi2_primary_noCuts_barrel->SetMarkerColor(1);
  h_trk_bendchi2_np_noCuts_barrel->Scale(1./h_trk_bendchi2_np_noCuts_barrel->Integral());
  h_trk_bendchi2_np_noCuts_barrel->SetLineColor(2);
  h_trk_bendchi2_np_noCuts_barrel->SetMarkerColor(2);
  h_trk_bendchi2_np_noCuts_barrel->Draw("HIST");
  h_trk_bendchi2_primary_noCuts_barrel->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_bendchi2_primary_noCuts_barrel,"Primary","l");
  l->AddEntry(h_trk_bendchi2_np_noCuts_barrel,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_bendchi2_noCuts_barrel.pdf");
  delete h_trk_bendchi2_primary_noCuts_barrel;
  delete h_trk_bendchi2_np_noCuts_barrel;

  h_trk_bendchi2_primary_noCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_primary_noCuts_disk);
  h_trk_bendchi2_np_noCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_np_noCuts_disk);
  h_trk_bendchi2_primary_noCuts_disk->Scale(1./h_trk_bendchi2_primary_noCuts_disk->Integral());
  h_trk_bendchi2_primary_noCuts_disk->SetLineColor(1);
  h_trk_bendchi2_primary_noCuts_disk->SetMarkerColor(1);
  h_trk_bendchi2_np_noCuts_disk->Scale(1./h_trk_bendchi2_np_noCuts_disk->Integral());
  h_trk_bendchi2_np_noCuts_disk->SetLineColor(2);
  h_trk_bendchi2_np_noCuts_disk->SetMarkerColor(2);
  h_trk_bendchi2_np_noCuts_disk->Draw("HIST");
  h_trk_bendchi2_primary_noCuts_disk->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_bendchi2_primary_noCuts_disk,"Primary","l");
  l->AddEntry(h_trk_bendchi2_np_noCuts_disk,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_bendchi2_noCuts_disk.pdf");
  delete h_trk_bendchi2_primary_noCuts_disk;
  delete h_trk_bendchi2_np_noCuts_disk;

  h_trk_bendchi2_primary_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_primary_qualCuts);
  h_trk_bendchi2_np_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_np_qualCuts);
  h_trk_bendchi2_primary_qualCuts->Scale(1./h_trk_bendchi2_primary_qualCuts->Integral());
  h_trk_bendchi2_primary_qualCuts->SetLineColor(1);
  h_trk_bendchi2_primary_qualCuts->SetMarkerColor(1);
  TH1F *h_trk_bendchi2_np_qualCutsNorm = (TH1F*)h_trk_bendchi2_np_qualCuts->Clone(); 
  h_trk_bendchi2_np_qualCutsNorm->Scale(1./h_trk_bendchi2_np_qualCutsNorm->Integral());
  h_trk_bendchi2_np_qualCutsNorm->SetLineColor(2);
  h_trk_bendchi2_np_qualCutsNorm->SetMarkerColor(2);
  h_trk_bendchi2_np_qualCutsNorm->Draw("HIST");
  h_trk_bendchi2_primary_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_bendchi2_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_bendchi2_np_qualCutsNorm,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_bendchi2_qualCuts.pdf");
  delete h_trk_bendchi2_np_qualCutsNorm;

  h_trk_bendchi2_primary_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_primary_allCuts);
  h_trk_bendchi2_np_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_np_allCuts);
  h_trk_bendchi2_primary_allCuts->Scale(1./h_trk_bendchi2_primary_allCuts->Integral());
  h_trk_bendchi2_primary_allCuts->SetLineColor(1);
  h_trk_bendchi2_primary_allCuts->SetMarkerColor(1);
  TH1F *h_trk_bendchi2_np_allCutsNorm = (TH1F*)h_trk_bendchi2_np_allCuts->Clone(); 
  h_trk_bendchi2_np_allCutsNorm->Scale(1./h_trk_bendchi2_np_allCutsNorm->Integral());
  h_trk_bendchi2_np_allCutsNorm->SetLineColor(2);
  h_trk_bendchi2_np_allCutsNorm->SetMarkerColor(2);
  raiseMax(h_trk_bendchi2_primary_allCuts,h_trk_bendchi2_np_allCutsNorm);
  h_trk_bendchi2_np_allCutsNorm->SetStats(0);
  h_trk_bendchi2_primary_allCuts->SetStats(0);
  h_trk_bendchi2_primary_allCuts->Draw("HIST");
  h_trk_bendchi2_np_allCutsNorm->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_bendchi2_primary_allCuts,"Primary","l");
  l->AddEntry(h_trk_bendchi2_np_allCutsNorm,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_bendchi2_allCuts.pdf");
  delete h_trk_bendchi2_np_allCutsNorm;

  h_trk_bendchi2_primary_allCuts_zoomOut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_primary_allCuts_zoomOut);
  h_trk_bendchi2_np_allCuts_zoomOut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_np_allCuts_zoomOut);
  h_trk_bendchi2_primary_allCuts_zoomOut->Scale(1./h_trk_bendchi2_primary_allCuts_zoomOut->Integral());
  h_trk_bendchi2_primary_allCuts_zoomOut->SetLineColor(1);
  h_trk_bendchi2_primary_allCuts_zoomOut->SetMarkerColor(1);
  h_trk_bendchi2_np_allCuts_zoomOut->Scale(1./h_trk_bendchi2_np_allCuts_zoomOut->Integral());
  h_trk_bendchi2_np_allCuts_zoomOut->SetLineColor(2);
  h_trk_bendchi2_np_allCuts_zoomOut->SetMarkerColor(2);
  raiseMax(h_trk_bendchi2_np_allCuts_zoomOut,h_trk_bendchi2_primary_allCuts_zoomOut);
  h_trk_bendchi2_np_allCuts_zoomOut->SetStats(0);
  h_trk_bendchi2_primary_allCuts_zoomOut->SetStats(0);
  h_trk_bendchi2_np_allCuts_zoomOut->Draw("HIST");
  h_trk_bendchi2_primary_allCuts_zoomOut->Draw("HIST,SAME");
  mySmallText(0.3, 0.9, 1, ctxt);
  l->Clear();
  l->SetFillColor(0);
  l->SetLineColor(0);
  l->SetTextSize(0.04);
  l->AddEntry(h_trk_bendchi2_primary_allCuts_zoomOut,"Primary","l");
  l->AddEntry(h_trk_bendchi2_np_allCuts_zoomOut,"NP","l");
  l->SetTextFont(42);
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_bendchi2_allCuts_zoomOut.pdf");
  delete h_trk_bendchi2_primary_allCuts_zoomOut;
  delete h_trk_bendchi2_np_allCuts_zoomOut;

  h_trk_bendchi2_primary_allCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_primary_allCuts_barrel);
  h_trk_bendchi2_np_allCuts_barrel->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_np_allCuts_barrel);
  h_trk_bendchi2_primary_allCuts_barrel->Scale(1./h_trk_bendchi2_primary_allCuts_barrel->Integral());
  h_trk_bendchi2_primary_allCuts_barrel->SetLineColor(1);
  h_trk_bendchi2_primary_allCuts_barrel->SetMarkerColor(1);
  h_trk_bendchi2_np_allCuts_barrel->Scale(1./h_trk_bendchi2_np_allCuts_barrel->Integral());
  h_trk_bendchi2_np_allCuts_barrel->SetLineColor(2);
  h_trk_bendchi2_np_allCuts_barrel->SetMarkerColor(2);
  h_trk_bendchi2_np_allCuts_barrel->Draw("HIST");
  h_trk_bendchi2_primary_allCuts_barrel->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_bendchi2_primary_allCuts_barrel,"Primary","l");
  l->AddEntry(h_trk_bendchi2_np_allCuts_barrel,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_bendchi2_allCuts_barrel.pdf");
  delete h_trk_bendchi2_primary_allCuts_barrel;
  delete h_trk_bendchi2_np_allCuts_barrel;

  h_trk_bendchi2_primary_allCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_primary_allCuts_disk);
  h_trk_bendchi2_np_allCuts_disk->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_np_allCuts_disk);
  h_trk_bendchi2_primary_allCuts_disk->Scale(1./h_trk_bendchi2_primary_allCuts_disk->Integral());
  h_trk_bendchi2_primary_allCuts_disk->SetLineColor(1);
  h_trk_bendchi2_primary_allCuts_disk->SetMarkerColor(1);
  h_trk_bendchi2_np_allCuts_disk->Scale(1./h_trk_bendchi2_np_allCuts_disk->Integral());
  h_trk_bendchi2_np_allCuts_disk->SetLineColor(2);
  h_trk_bendchi2_np_allCuts_disk->SetMarkerColor(2);
  h_trk_bendchi2_np_allCuts_disk->Draw("HIST");
  h_trk_bendchi2_primary_allCuts_disk->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_bendchi2_primary_allCuts_disk,"Primary","l");
  l->AddEntry(h_trk_bendchi2_np_allCuts_disk,"NP","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBG_bendchi2_allCuts_disk.pdf");
  delete h_trk_bendchi2_primary_allCuts_disk;
  delete h_trk_bendchi2_np_allCuts_disk;

  TH1F *h_trk_bendchi2_fake_qualCutsNorm = (TH1F*)h_trk_bendchi2_fake_qualCuts->Clone(); 
  h_trk_bendchi2_fake_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_fake_qualCutsNorm);
  TH1F *h_trk_bendchi2_PU_qualCutsNorm = (TH1F*)h_trk_bendchi2_PU_qualCuts->Clone(); 
  h_trk_bendchi2_PU_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_PU_qualCutsNorm);
  TH1F *h_trk_bendchi2_notHiggs_qualCutsNorm = (TH1F*)h_trk_bendchi2_notHiggs_qualCuts->Clone(); 
  h_trk_bendchi2_notHiggs_qualCutsNorm->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_notHiggs_qualCutsNorm);
  h_trk_bendchi2_fake_qualCutsNorm->Scale(1./h_trk_bendchi2_np_qualCuts->Integral());
  h_trk_bendchi2_fake_qualCutsNorm->SetLineColor(2);
  h_trk_bendchi2_fake_qualCutsNorm->SetMarkerColor(2);
  h_trk_bendchi2_PU_qualCutsNorm->Scale(1./h_trk_bendchi2_np_qualCuts->Integral());
  h_trk_bendchi2_PU_qualCutsNorm->SetLineColor(3);
  h_trk_bendchi2_PU_qualCutsNorm->SetMarkerColor(3);
  h_trk_bendchi2_notHiggs_qualCutsNorm->Scale(1./h_trk_bendchi2_np_qualCuts->Integral());
  h_trk_bendchi2_notHiggs_qualCutsNorm->SetLineColor(4);
  h_trk_bendchi2_notHiggs_qualCutsNorm->SetMarkerColor(4);
  auto hs_bendchi2_qualCuts = new THStack("hs_bendchi2_qualCuts","Stacked BG histograms");
  hs_bendchi2_qualCuts->Add(h_trk_bendchi2_fake_qualCutsNorm);
  hs_bendchi2_qualCuts->Add(h_trk_bendchi2_PU_qualCutsNorm);
  hs_bendchi2_qualCuts->Add(h_trk_bendchi2_notHiggs_qualCutsNorm);
  hs_bendchi2_qualCuts->Draw("HIST");
  h_trk_bendchi2_primary_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_bendchi2_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_bendchi2_fake_qualCutsNorm,"Fake","l");
  l->AddEntry(h_trk_bendchi2_PU_qualCutsNorm,"PU","l");
  l->AddEntry(h_trk_bendchi2_notHiggs_qualCutsNorm,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_bendchi2_qualCutsNorm.pdf");
  delete h_trk_bendchi2_fake_qualCutsNorm;
  delete h_trk_bendchi2_PU_qualCutsNorm;
  delete h_trk_bendchi2_notHiggs_qualCutsNorm;
  delete h_trk_bendchi2_np_qualCuts;
  delete hs_bendchi2_qualCuts;

  h_trk_bendchi2_fake_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_fake_qualCuts);
  h_trk_bendchi2_PU_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_PU_qualCuts);
  h_trk_bendchi2_notHiggs_qualCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_notHiggs_qualCuts);
  h_trk_bendchi2_fake_qualCuts->Scale(1./h_trk_bendchi2_fake_qualCuts->Integral());
  h_trk_bendchi2_fake_qualCuts->SetLineColor(2);
  h_trk_bendchi2_fake_qualCuts->SetMarkerColor(2);
  h_trk_bendchi2_PU_qualCuts->Scale(1./h_trk_bendchi2_PU_qualCuts->Integral());
  h_trk_bendchi2_PU_qualCuts->SetLineColor(3);
  h_trk_bendchi2_PU_qualCuts->SetMarkerColor(3);
  h_trk_bendchi2_notHiggs_qualCuts->Scale(1./h_trk_bendchi2_notHiggs_qualCuts->Integral());
  h_trk_bendchi2_notHiggs_qualCuts->SetLineColor(4);
  h_trk_bendchi2_notHiggs_qualCuts->SetMarkerColor(4);
  h_trk_bendchi2_notHiggs_qualCuts->Draw("HIST");
  h_trk_bendchi2_primary_qualCuts->Draw("HIST,SAME");
  h_trk_bendchi2_fake_qualCuts->Draw("HIST,SAME");
  h_trk_bendchi2_PU_qualCuts->Draw("HIST,SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_trk_bendchi2_primary_qualCuts,"Primary","l");
  l->AddEntry(h_trk_bendchi2_fake_qualCuts,"Fake","l");
  l->AddEntry(h_trk_bendchi2_PU_qualCuts,"PU","l");
  l->AddEntry(h_trk_bendchi2_notHiggs_qualCuts,"Not H->4b","l");
  l->Draw();
  c.SaveAs(DIR + "/h_signalVsBGOverlay_bendchi2_qualCuts.pdf");
  delete h_trk_bendchi2_primary_qualCuts;
  delete h_trk_bendchi2_fake_qualCuts;
  delete h_trk_bendchi2_PU_qualCuts;
  delete h_trk_bendchi2_notHiggs_qualCuts;

  h_trk_chi2rphidof_H->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_H);
  h_trk_chi2rphidof_H->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_chi2rphidof_H->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_chi2rphidof_H->GetName() + ".pdf");
  delete h_trk_chi2rphidof_H;

  h_trk_chi2rzdof_H->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_H);
  h_trk_chi2rzdof_H->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_chi2rzdof_H->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_chi2rzdof_H->GetName() + ".pdf");
  delete h_trk_chi2rzdof_H;
   
  h_trk_bendchi2_H->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_H);
  h_trk_bendchi2_H->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_bendchi2_H->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_bendchi2_H->GetName() + ".pdf");
  delete h_trk_bendchi2_H;

  h_trk_chi2rphidof_L->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_L);
  h_trk_chi2rphidof_L->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_chi2rphidof_L->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_chi2rphidof_L->GetName() + ".pdf");
  delete h_trk_chi2rphidof_L;

  h_trk_chi2rzdof_L->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_L);
  h_trk_chi2rzdof_L->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_chi2rzdof_L->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_chi2rzdof_L->GetName() + ".pdf");
  delete h_trk_chi2rzdof_L;

  h_trk_bendchi2_L->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_L);
  h_trk_bendchi2_L->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_bendchi2_L->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_bendchi2_L->GetName() + ".pdf");
  delete h_trk_bendchi2_L;

  h_trk_chi2rphidof_C->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_C);
  h_trk_chi2rphidof_C->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_chi2rphidof_C->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_chi2rphidof_C->GetName() + ".pdf");
  delete h_trk_chi2rphidof_C;

  h_trk_chi2rzdof_C->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_C);
  h_trk_chi2rzdof_C->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_chi2rzdof_C->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_chi2rzdof_C->GetName() + ".pdf");
  delete h_trk_chi2rzdof_C;

  h_trk_bendchi2_C->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_C);
  h_trk_bendchi2_C->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_bendchi2_C->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_bendchi2_C->GetName() + ".pdf");
  delete h_trk_bendchi2_C;

  h_trk_chi2rphidof_I->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_I);
  h_trk_chi2rphidof_I->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_chi2rphidof_I->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_chi2rphidof_I->GetName() + ".pdf");
  delete h_trk_chi2rphidof_I;
   
  h_trk_chi2rzdof_I->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_I);
  h_trk_chi2rzdof_I->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_chi2rzdof_I->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_chi2rzdof_I->GetName() + ".pdf");
  delete h_trk_chi2rzdof_I;
   
  h_trk_bendchi2_I->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_I);
  h_trk_bendchi2_I->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_bendchi2_I->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_bendchi2_I->GetName() + ".pdf");
  delete h_trk_bendchi2_I;
   
  h_trk_chi2rphidof_F->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_F);
  h_trk_chi2rphidof_F->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_chi2rphidof_F->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_chi2rphidof_F->GetName() + ".pdf");
  delete h_trk_chi2rphidof_F;
   
  h_trk_chi2rzdof_F->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_F);
  h_trk_chi2rzdof_F->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_chi2rzdof_F->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_chi2rzdof_F->GetName() + ".pdf");
  delete h_trk_chi2rzdof_F;
   
  h_trk_bendchi2_F->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_F);
  h_trk_bendchi2_F->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_bendchi2_F->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_bendchi2_F->GetName() + ".pdf");
  delete h_trk_bendchi2_F;
   
  h_trk_chi2rphidof_P->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_P);
  h_trk_chi2rphidof_P->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_chi2rphidof_P->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_chi2rphidof_P->GetName() + ".pdf");
  delete h_trk_chi2rphidof_P;
   
  h_trk_chi2rzdof_P->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_P);
  h_trk_chi2rzdof_P->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_chi2rzdof_P->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_chi2rzdof_P->GetName() + ".pdf");
  delete h_trk_chi2rzdof_P;
   
  h_trk_bendchi2_P->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_P);
  h_trk_bendchi2_P->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_bendchi2_P->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_bendchi2_P->GetName() + ".pdf");
  delete h_trk_bendchi2_P;
   
  h_trk_chi2rphidof_D->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rphidof_D);
  h_trk_chi2rphidof_D->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_chi2rphidof_D->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_chi2rphidof_D->GetName() + ".pdf");
  delete h_trk_chi2rphidof_D;
   
  h_trk_chi2rzdof_D->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_chi2rzdof_D);
  h_trk_chi2rzdof_D->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_chi2rzdof_D->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_chi2rzdof_D->GetName() + ".pdf");
  delete h_trk_chi2rzdof_D;
   
  h_trk_bendchi2_D->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_bendchi2_D);
  h_trk_bendchi2_D->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_bendchi2_D->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_bendchi2_D->GetName() + ".pdf");
  delete h_trk_bendchi2_D;
   
  h_match_tp_pt_noCuts->Sumw2();
  h_tp_pt_noCuts->Sumw2();
  TH1F* h_eff_pt_noCuts = (TH1F*)h_match_tp_pt_noCuts->Clone();
  h_eff_pt_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_pt_noCuts);
  h_eff_pt_noCuts->SetName("eff_pt_noCuts");
  h_eff_pt_noCuts->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
  h_eff_pt_noCuts->GetYaxis()->SetTitle("Efficiency");
  h_eff_pt_noCuts->Divide(h_match_tp_pt_noCuts, h_tp_pt_noCuts, 1.0, 1.0, "B");
  h_eff_pt_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_pt_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_pt_noCuts->GetName() + ".pdf");
  delete h_eff_pt_noCuts;
  delete h_match_tp_pt_noCuts;
  delete h_tp_pt_noCuts;

  h_tpAssoc_pt_matched->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_tpAssoc_pt_matched);
  h_tpAssoc_pt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_tpAssoc_pt);
  h_tpAssoc_pt_matched->Sumw2();
  h_tpAssoc_pt->Sumw2();
  TH1F* h_assocEff_pt = (TH1F*)h_tpAssoc_pt_matched->Clone();
  h_assocEff_pt->GetYaxis()->SetNoExponent(kTRUE);
  h_assocEff_pt->SetName("assocEff_pt");
  h_assocEff_pt->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
  h_assocEff_pt->GetYaxis()->SetTitle("Efficiency");
  h_assocEff_pt->Divide(h_tpAssoc_pt_matched, h_tpAssoc_pt, 1.0, 1.0, "B");
  h_assocEff_pt->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_assocEff_pt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_assocEff_pt->GetName() + ".pdf");
  delete h_assocEff_pt;
  delete h_tpAssoc_pt_matched;
  delete h_tpAssoc_pt;
   
  h_trkAssoc_pt_noMatch->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trkAssoc_pt_noMatch);
  h_trkAssoc_pt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trkAssoc_pt);
  h_trkAssoc_pt_noMatch->Sumw2();
  h_trkAssoc_pt->Sumw2();
  TH1F* h_assocFakeRate_pt = (TH1F*)h_trkAssoc_pt_noMatch->Clone();
  h_assocFakeRate_pt->GetYaxis()->SetNoExponent(kTRUE);
  h_assocFakeRate_pt->SetName("assocFakeRate_pt");
  h_assocFakeRate_pt->GetXaxis()->SetTitle("Track p_{T} (GeV)");
  h_assocFakeRate_pt->GetYaxis()->SetTitle("FakeRate");
  h_assocFakeRate_pt->Divide(h_trkAssoc_pt_noMatch, h_trkAssoc_pt, 1.0, 1.0, "B");
  h_assocFakeRate_pt->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_assocFakeRate_pt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_assocFakeRate_pt->GetName() + ".pdf");
  delete h_assocFakeRate_pt;
  delete h_trkAssoc_pt_noMatch;
  delete h_trkAssoc_pt;

  h_match_tp_eta_noCuts->Sumw2();
  h_tp_eta_noCuts->Sumw2();
  TH1F* h_eff_eta_noCuts = (TH1F*)h_match_tp_eta_noCuts->Clone();
  h_eff_eta_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_eta_noCuts);
  h_eff_eta_noCuts->SetName("eff_eta_noCuts");
  h_eff_eta_noCuts->GetXaxis()->SetTitle("Tracking Particle #eta");
  h_eff_eta_noCuts->GetYaxis()->SetTitle("Efficiency");
  h_eff_eta_noCuts->Divide(h_match_tp_eta_noCuts, h_tp_eta_noCuts, 1.0, 1.0, "B");
  h_eff_eta_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_eta_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_eta_noCuts->GetName() + ".pdf");
  delete h_eff_eta_noCuts;
  delete h_match_tp_eta_noCuts;
  

  h_match_tp_d0_noCuts->Sumw2();
  h_tp_d0_noCuts->Sumw2();
  TH1F* h_eff_d0_noCuts = (TH1F*)h_match_tp_d0_noCuts->Clone();
  h_eff_d0_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_d0_noCuts);
  h_eff_d0_noCuts->SetName("eff_d0_noCuts");
  h_eff_d0_noCuts->GetXaxis()->SetTitle("Tracking Particle d_{0} (cm)");
  h_eff_d0_noCuts->GetYaxis()->SetTitle("Efficiency");
  h_eff_d0_noCuts->Divide(h_match_tp_d0_noCuts, h_tp_d0_noCuts, 1.0, 1.0, "B");
  h_eff_d0_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_d0_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_d0_noCuts->GetName() + ".pdf");
  delete h_eff_d0_noCuts;
  delete h_match_tp_d0_noCuts;
  delete h_tp_d0_noCuts;

  h_match_tp_pt_maxD0Cut->Sumw2();
  h_tp_pt_maxD0Cut->Sumw2();
  TH1F* h_eff_pt_maxD0Cut = (TH1F*)h_match_tp_pt_maxD0Cut->Clone();
  h_eff_pt_maxD0Cut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_pt_maxD0Cut);
  h_eff_pt_maxD0Cut->SetName("eff_pt_maxD0Cut");
  h_eff_pt_maxD0Cut->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
  h_eff_pt_maxD0Cut->GetYaxis()->SetTitle("Efficiency");
  h_eff_pt_maxD0Cut->Divide(h_match_tp_pt_maxD0Cut, h_tp_pt_maxD0Cut, 1.0, 1.0, "B");
  h_eff_pt_maxD0Cut->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_pt_maxD0Cut->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_pt_maxD0Cut->GetName() + ".pdf");
  delete h_eff_pt_maxD0Cut;
  delete h_match_tp_pt_maxD0Cut;
  delete h_tp_pt_maxD0Cut;

  h_match_tp_eta_maxD0Cut->Sumw2();
  h_tp_eta_maxD0Cut->Sumw2();
  TH1F* h_eff_eta_maxD0Cut = (TH1F*)h_match_tp_eta_maxD0Cut->Clone();
  h_eff_eta_maxD0Cut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_eta_maxD0Cut);
  h_eff_eta_maxD0Cut->SetName("eff_eta_maxD0Cut");
  h_eff_eta_maxD0Cut->GetXaxis()->SetTitle("Tracking Particle #eta");
  h_eff_eta_maxD0Cut->GetYaxis()->SetTitle("Efficiency");
  h_eff_eta_maxD0Cut->Divide(h_match_tp_eta_maxD0Cut, h_tp_eta_maxD0Cut, 1.0, 1.0, "B");
  h_eff_eta_maxD0Cut->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_eta_maxD0Cut->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_eta_maxD0Cut->GetName() + ".pdf");
  delete h_eff_eta_maxD0Cut;
  delete h_match_tp_eta_maxD0Cut;
  delete h_tp_eta_maxD0Cut;

  h_match_tp_d0_maxD0Cut->Sumw2();
  h_tp_d0_maxD0Cut->Sumw2();
  TH1F* h_eff_d0_maxD0Cut = (TH1F*)h_match_tp_d0_maxD0Cut->Clone();
  h_eff_d0_maxD0Cut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_d0_maxD0Cut);
  h_eff_d0_maxD0Cut->SetName("eff_d0_maxD0Cut");
  h_eff_d0_maxD0Cut->GetXaxis()->SetTitle("Tracking Particle d_{0} (cm)");
  h_eff_d0_maxD0Cut->GetYaxis()->SetTitle("Efficiency");
  h_eff_d0_maxD0Cut->Divide(h_match_tp_d0_maxD0Cut, h_tp_d0_maxD0Cut, 1.0, 1.0, "B");
  h_eff_d0_maxD0Cut->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_d0_maxD0Cut->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_d0_maxD0Cut->GetName() + ".pdf");
  delete h_eff_d0_maxD0Cut;
  delete h_match_tp_d0_maxD0Cut;
  delete h_tp_d0_maxD0Cut;

  h_match_tp_pt_minD0Cut->Sumw2();
  h_tp_pt_minD0Cut->Sumw2();
  TH1F* h_eff_pt_minD0Cut = (TH1F*)h_match_tp_pt_minD0Cut->Clone();
  h_eff_pt_minD0Cut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_pt_minD0Cut);
  h_eff_pt_minD0Cut->SetName("eff_pt_minD0Cut");
  h_eff_pt_minD0Cut->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
  h_eff_pt_minD0Cut->GetYaxis()->SetTitle("Efficiency");
  h_eff_pt_minD0Cut->Divide(h_match_tp_pt_minD0Cut, h_tp_pt_minD0Cut, 1.0, 1.0, "B");
  h_eff_pt_minD0Cut->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_pt_minD0Cut->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_pt_minD0Cut->GetName() + ".pdf");
  delete h_eff_pt_minD0Cut;
  delete h_match_tp_pt_minD0Cut;
  delete h_tp_pt_minD0Cut;

  h_match_tp_eta_minD0Cut->Sumw2();
  h_tp_eta_minD0Cut->Sumw2();
  TH1F* h_eff_eta_minD0Cut = (TH1F*)h_match_tp_eta_minD0Cut->Clone();
  h_eff_eta_minD0Cut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_eta_minD0Cut);
  h_eff_eta_minD0Cut->SetName("eff_eta_minD0Cut");
  h_eff_eta_minD0Cut->GetXaxis()->SetTitle("Tracking Particle #eta");
  h_eff_eta_minD0Cut->GetYaxis()->SetTitle("Efficiency");
  h_eff_eta_minD0Cut->Divide(h_match_tp_eta_minD0Cut, h_tp_eta_minD0Cut, 1.0, 1.0, "B");
  h_eff_eta_minD0Cut->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_eta_minD0Cut->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_eta_minD0Cut->GetName() + ".pdf");
  delete h_eff_eta_minD0Cut;
  delete h_match_tp_eta_minD0Cut;
  delete h_tp_eta_minD0Cut;
   
  h_match_tp_d0_minD0Cut->Sumw2();
  h_tp_d0_minD0Cut->Sumw2();
  TH1F* h_eff_d0_minD0Cut = (TH1F*)h_match_tp_d0_minD0Cut->Clone();
  h_eff_d0_minD0Cut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_d0_minD0Cut);
  h_eff_d0_minD0Cut->SetName("eff_d0_minD0Cut");
  h_eff_d0_minD0Cut->GetXaxis()->SetTitle("Tracking Particle d_{0} (cm)");
  h_eff_d0_minD0Cut->GetYaxis()->SetTitle("Efficiency");
  h_eff_d0_minD0Cut->Divide(h_match_tp_d0_minD0Cut, h_tp_d0_minD0Cut, 1.0, 1.0, "B");
  h_eff_d0_minD0Cut->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_d0_minD0Cut->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_d0_minD0Cut->GetName() + ".pdf");
  delete h_eff_d0_minD0Cut;
  delete h_match_tp_d0_minD0Cut;
  delete h_tp_d0_minD0Cut;

  h_match_tp_pt_minPtCut->Sumw2();
  h_tp_pt_minPtCut->Sumw2();
  TH1F* h_eff_pt_minPtCut = (TH1F*)h_match_tp_pt_minPtCut->Clone();
  h_eff_pt_minPtCut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_pt_minPtCut);
  h_eff_pt_minPtCut->SetName("eff_pt_minPtCut");
  h_eff_pt_minPtCut->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
  h_eff_pt_minPtCut->GetYaxis()->SetTitle("Efficiency");
  h_eff_pt_minPtCut->Divide(h_match_tp_pt_minPtCut, h_tp_pt_minPtCut, 1.0, 1.0, "B");
  h_eff_pt_minPtCut->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_pt_minPtCut->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_pt_minPtCut->GetName() + ".pdf");
  delete h_eff_pt_minPtCut;
  delete h_match_tp_pt_minPtCut;
  delete h_tp_pt_minPtCut;

  h_match_tp_eta_minPtCut->Sumw2();
  h_tp_eta_minPtCut->Sumw2();
  TH1F* h_eff_eta_minPtCut = (TH1F*)h_match_tp_eta_minPtCut->Clone();
  h_eff_eta_minPtCut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_eta_minPtCut);
  h_eff_eta_minPtCut->SetName("eff_eta_minPtCut");
  h_eff_eta_minPtCut->GetXaxis()->SetTitle("Tracking Particle #eta");
  h_eff_eta_minPtCut->GetYaxis()->SetTitle("Efficiency");
  h_eff_eta_minPtCut->Divide(h_match_tp_eta_minPtCut, h_tp_eta_minPtCut, 1.0, 1.0, "B");
  h_eff_eta_minPtCut->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_eta_minPtCut->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_eta_minPtCut->GetName() + ".pdf");
  delete h_eff_eta_minPtCut;
  delete h_match_tp_eta_minPtCut;
  delete h_tp_eta_minPtCut;

  h_match_tp_d0_minPtCut->Sumw2();
  h_tp_d0_minPtCut->Sumw2();
  TH1F* h_eff_d0_minPtCut = (TH1F*)h_match_tp_d0_minPtCut->Clone();
  h_eff_d0_minPtCut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_d0_minPtCut);
  h_eff_d0_minPtCut->SetName("eff_d0_minPtCut");
  h_eff_d0_minPtCut->GetXaxis()->SetTitle("Tracking Particle d_{0} (cm)");
  h_eff_d0_minPtCut->GetYaxis()->SetTitle("Efficiency");
  h_eff_d0_minPtCut->Divide(h_match_tp_d0_minPtCut, h_tp_d0_minPtCut, 1.0, 1.0, "B");
  h_eff_d0_minPtCut->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_d0_minPtCut->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_d0_minPtCut->GetName() + ".pdf");
  delete h_eff_d0_minPtCut;
  delete h_match_tp_d0_minPtCut;
  delete h_tp_d0_minPtCut;

  h_match_tp_pt_maxPtCut->Sumw2();
  h_tp_pt_maxPtCut->Sumw2();
  TH1F* h_eff_pt_maxPtCut = (TH1F*)h_match_tp_pt_maxPtCut->Clone();
  h_eff_pt_maxPtCut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_pt_maxPtCut);
  h_eff_pt_maxPtCut->SetName("eff_pt_maxPtCut");
  h_eff_pt_maxPtCut->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
  h_eff_pt_maxPtCut->GetYaxis()->SetTitle("Efficiency");
  h_eff_pt_maxPtCut->Divide(h_match_tp_pt_maxPtCut, h_tp_pt_maxPtCut, 1.0, 1.0, "B");
  h_eff_pt_maxPtCut->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_pt_maxPtCut->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_pt_maxPtCut->GetName() + ".pdf");
  delete h_eff_pt_maxPtCut;
  delete h_match_tp_pt_maxPtCut;
  delete h_tp_pt_maxPtCut;

  h_match_tp_eta_maxPtCut->Sumw2();
  h_tp_eta_maxPtCut->Sumw2();
  TH1F* h_eff_eta_maxPtCut = (TH1F*)h_match_tp_eta_maxPtCut->Clone();
  h_eff_eta_maxPtCut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_eta_maxPtCut);
  h_eff_eta_maxPtCut->SetName("eff_eta_maxPtCut");
  h_eff_eta_maxPtCut->GetXaxis()->SetTitle("Tracking Particle #eta");
  h_eff_eta_maxPtCut->GetYaxis()->SetTitle("Efficiency");
  h_eff_eta_maxPtCut->Divide(h_match_tp_eta_maxPtCut, h_tp_eta_maxPtCut, 1.0, 1.0, "B");
  h_eff_eta_maxPtCut->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_eta_maxPtCut->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_eta_maxPtCut->GetName() + ".pdf");
  delete h_eff_eta_maxPtCut;
  delete h_match_tp_eta_maxPtCut;
  delete h_tp_eta_maxPtCut;

  h_match_tp_d0_maxPtCut->Sumw2();
  h_tp_d0_maxPtCut->Sumw2();
  TH1F* h_eff_d0_maxPtCut = (TH1F*)h_match_tp_d0_maxPtCut->Clone();
  h_eff_d0_maxPtCut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_d0_maxPtCut);
  h_eff_d0_maxPtCut->SetName("eff_d0_maxPtCut");
  h_eff_d0_maxPtCut->GetXaxis()->SetTitle("Tracking Particle d_{0} (cm)");
  h_eff_d0_maxPtCut->GetYaxis()->SetTitle("Efficiency");
  h_eff_d0_maxPtCut->Divide(h_match_tp_d0_maxPtCut, h_tp_d0_maxPtCut, 1.0, 1.0, "B");
  h_eff_d0_maxPtCut->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_d0_maxPtCut->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_d0_maxPtCut->GetName() + ".pdf");
  delete h_eff_d0_maxPtCut;
  delete h_match_tp_d0_maxPtCut;
  delete h_tp_d0_maxPtCut;

  h_match_tp_pt_maxEtaCut->Sumw2();
  h_tp_pt_maxEtaCut->Sumw2();
  TH1F* h_eff_pt_maxEtaCut = (TH1F*)h_match_tp_pt_maxEtaCut->Clone();
  h_eff_pt_maxEtaCut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_pt_maxEtaCut);
  h_eff_pt_maxEtaCut->SetName("eff_pt_maxEtaCut");
  h_eff_pt_maxEtaCut->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
  h_eff_pt_maxEtaCut->GetYaxis()->SetTitle("Efficiency");
  h_eff_pt_maxEtaCut->Divide(h_match_tp_pt_maxEtaCut, h_tp_pt_maxEtaCut, 1.0, 1.0, "B");
  h_eff_pt_maxEtaCut->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_pt_maxEtaCut->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_pt_maxEtaCut->GetName() + ".pdf");
  delete h_eff_pt_maxEtaCut;
  //delete h_match_tp_pt_maxEtaCut;  

  h_match_tp_eta_maxEtaCut->Sumw2();
  h_tp_eta_maxEtaCut->Sumw2();
  TH1F* h_eff_eta_maxEtaCut = (TH1F*)h_match_tp_eta_maxEtaCut->Clone();
  h_eff_eta_maxEtaCut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_eta_maxEtaCut);
  h_eff_eta_maxEtaCut->SetName("eff_eta_maxEtaCut");
  h_eff_eta_maxEtaCut->GetXaxis()->SetTitle("Tracking Particle #eta");
  h_eff_eta_maxEtaCut->GetYaxis()->SetTitle("Efficiency");
  h_eff_eta_maxEtaCut->Divide(h_match_tp_eta_maxEtaCut, h_tp_eta_maxEtaCut, 1.0, 1.0, "B");
  h_eff_eta_maxEtaCut->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_eta_maxEtaCut->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_eta_maxEtaCut->GetName() + ".pdf");
  delete h_eff_eta_maxEtaCut;
  delete h_match_tp_eta_maxEtaCut;
  

  h_match_tp_d0_maxEtaCut->Sumw2();
  h_tp_d0_maxEtaCut->Sumw2();
  TH1F* h_eff_d0_maxEtaCut = (TH1F*)h_match_tp_d0_maxEtaCut->Clone();
  h_eff_d0_maxEtaCut->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_d0_maxEtaCut);
  h_eff_d0_maxEtaCut->SetName("eff_d0_maxEtaCut");
  h_eff_d0_maxEtaCut->GetXaxis()->SetTitle("Tracking Particle d_{0} (cm)");
  h_eff_d0_maxEtaCut->GetYaxis()->SetTitle("Efficiency");
  h_eff_d0_maxEtaCut->Divide(h_match_tp_d0_maxEtaCut, h_tp_d0_maxEtaCut, 1.0, 1.0, "B");
  h_eff_d0_maxEtaCut->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_d0_maxEtaCut->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_d0_maxEtaCut->GetName() + ".pdf");
  delete h_eff_d0_maxEtaCut;
  delete h_match_tp_d0_maxEtaCut;
  delete h_tp_d0_maxEtaCut;

  h_match_tp_pt_allCuts->Sumw2();
  h_tp_pt_allCuts->Sumw2();
  TH1F* h_eff_pt_allCuts = (TH1F*)h_match_tp_pt_allCuts->Clone();
  h_eff_pt_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_pt_allCuts);
  h_eff_pt_allCuts->SetName("eff_pt_allCuts");
  h_eff_pt_allCuts->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
  h_eff_pt_allCuts->GetYaxis()->SetTitle("Efficiency");
  h_eff_pt_allCuts->Divide(h_match_tp_pt_allCuts, h_tp_pt_allCuts, 1.0, 1.0, "B");
  h_eff_pt_allCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_pt_allCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_pt_allCuts->GetName() + ".pdf");
  delete h_eff_pt_allCuts;
  //delete h_match_tp_pt_allCuts;
  

  h_match_tp_eta_allCuts->Sumw2();
  h_tp_eta_allCuts->Sumw2();
  TH1F* h_eff_eta_allCuts = (TH1F*)h_match_tp_eta_allCuts->Clone();
  h_eff_eta_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_eta_allCuts);
  h_eff_eta_allCuts->SetName("eff_eta_allCuts");
  h_eff_eta_allCuts->GetXaxis()->SetTitle("Tracking Particle #eta");
  h_eff_eta_allCuts->GetYaxis()->SetTitle("Efficiency");
  h_eff_eta_allCuts->Divide(h_match_tp_eta_allCuts, h_tp_eta_allCuts, 1.0, 1.0, "B");
  h_eff_eta_allCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_eta_allCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_eta_allCuts->GetName() + ".pdf");
  delete h_eff_eta_allCuts;
  delete h_match_tp_eta_allCuts;
  

  h_match_tp_d0_allCuts->Sumw2();
  h_tp_d0_allCuts->Sumw2();
  TH1F* h_eff_d0_allCuts = (TH1F*)h_match_tp_d0_allCuts->Clone();
  h_eff_d0_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_d0_allCuts);
  h_eff_d0_allCuts->SetName("eff_d0_allCuts");
  h_eff_d0_allCuts->GetXaxis()->SetTitle("Tracking Particle d_{0} (cm)");
  h_eff_d0_allCuts->GetYaxis()->SetTitle("Efficiency");
  h_eff_d0_allCuts->Divide(h_match_tp_d0_allCuts, h_tp_d0_allCuts, 1.0, 1.0, "B");
  h_eff_d0_allCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_d0_allCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_d0_allCuts->GetName() + ".pdf");
  delete h_eff_d0_allCuts;
  delete h_match_tp_d0_allCuts;
  delete h_tp_d0_allCuts;

  h_trk_matchtp_pt_allCuts->Sumw2();
  //h_tp_pt_allCuts->Sumw2();
  TH1F* h_psEff_pt_allCuts = (TH1F*)h_trk_matchtp_pt_allCuts->Clone();
  h_psEff_pt_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_psEff_pt_allCuts);
  h_psEff_pt_allCuts->SetName("psEff_pt_allCuts");
  h_psEff_pt_allCuts->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
  h_psEff_pt_allCuts->GetYaxis()->SetTitle("Efficiency");
  h_psEff_pt_allCuts->Divide(h_trk_matchtp_pt_allCuts, h_tp_pt_allCuts, 1.0, 1.0, "B");
  h_psEff_pt_allCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_psEff_pt_allCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_psEff_pt_allCuts->GetName() + ".pdf");
  delete h_psEff_pt_allCuts;
  delete h_trk_matchtp_pt_allCuts;
  delete h_tp_pt_allCuts;

  h_trk_matchtp_eta_allCuts->Sumw2();
  //h_tp_eta_allCuts->Sumw2();
  TH1F* h_psEff_eta_allCuts = (TH1F*)h_trk_matchtp_eta_allCuts->Clone();
  h_psEff_eta_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_psEff_eta_allCuts);
  h_psEff_eta_allCuts->SetName("psEff_eta_allCuts");
  h_psEff_eta_allCuts->GetXaxis()->SetTitle("Tracking Particle #eta (GeV)");
  h_psEff_eta_allCuts->GetYaxis()->SetTitle("Efficiency");
  h_psEff_eta_allCuts->Divide(h_trk_matchtp_eta_allCuts, h_tp_eta_allCuts, 1.0, 1.0, "B");
  h_psEff_eta_allCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_psEff_eta_allCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_psEff_eta_allCuts->GetName() + ".pdf");
  delete h_psEff_eta_allCuts;
  delete h_trk_matchtp_eta_allCuts;
  delete h_tp_eta_allCuts;

  h_trk_matchtp_dxy_allCuts->Sumw2();
  h_tp_dxy_allCuts->Sumw2();
  TH1F* h_psEff_dxy_allCuts = (TH1F*)h_trk_matchtp_dxy_allCuts->Clone();
  h_psEff_dxy_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_psEff_dxy_allCuts);
  h_psEff_dxy_allCuts->SetName("psEff_dxy_allCuts");
  h_psEff_dxy_allCuts->GetXaxis()->SetTitle("Tracking Particle #dxy (GeV)");
  h_psEff_dxy_allCuts->GetYaxis()->SetTitle("Efficiency");
  h_psEff_dxy_allCuts->Divide(h_trk_matchtp_dxy_allCuts, h_tp_dxy_allCuts, 1.0, 1.0, "B");
  h_psEff_dxy_allCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_psEff_dxy_allCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_psEff_dxy_allCuts->GetName() + ".pdf");
  delete h_psEff_dxy_allCuts;
  delete h_trk_matchtp_dxy_allCuts;
  delete h_tp_dxy_allCuts;

  h_trk_matchtp_pt_noCuts->Sumw2();
  h_tp_pt_noCuts_primary->Sumw2();
  TH1F* h_psEff_pt_noCuts = (TH1F*)h_trk_matchtp_pt_noCuts->Clone();
  h_psEff_pt_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_psEff_pt_noCuts);
  h_psEff_pt_noCuts->SetName("psEff_pt_noCuts");
  h_psEff_pt_noCuts->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
  h_psEff_pt_noCuts->GetYaxis()->SetTitle("Efficiency");
  h_psEff_pt_noCuts->Divide(h_trk_matchtp_pt_noCuts, h_tp_pt_noCuts_primary, 1.0, 1.0, "B");
  h_psEff_pt_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_psEff_pt_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_psEff_pt_noCuts->GetName() + ".pdf");
  delete h_psEff_pt_noCuts;
  delete h_trk_matchtp_pt_noCuts;

  h_trk_matchtp_eta_noCuts->Sumw2();
  //h_tp_eta_noCuts->Sumw2();
  TH1F* h_psEff_eta_noCuts = (TH1F*)h_trk_matchtp_eta_noCuts->Clone();
  h_psEff_eta_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_psEff_eta_noCuts);
  h_psEff_eta_noCuts->SetName("psEff_eta_noCuts");
  h_psEff_eta_noCuts->GetXaxis()->SetTitle("Tracking Particle #eta (GeV)");
  h_psEff_eta_noCuts->GetYaxis()->SetTitle("Efficiency");
  h_psEff_eta_noCuts->Divide(h_trk_matchtp_eta_noCuts, h_tp_eta_noCuts, 1.0, 1.0, "B");
  h_psEff_eta_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_psEff_eta_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_psEff_eta_noCuts->GetName() + ".pdf");
  delete h_psEff_eta_noCuts;
  delete h_trk_matchtp_eta_noCuts;
  delete h_tp_eta_noCuts;

  h_trk_matchtp_dxy_noCuts->Sumw2();
  h_tp_dxy_noCuts->Sumw2();
  TH1F* h_psEff_dxy_noCuts = (TH1F*)h_trk_matchtp_dxy_noCuts->Clone();
  h_psEff_dxy_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_psEff_dxy_noCuts);
  h_psEff_dxy_noCuts->SetName("psEff_dxy_noCuts");
  h_psEff_dxy_noCuts->GetXaxis()->SetTitle("Tracking Particle #dxy (GeV)");
  h_psEff_dxy_noCuts->GetYaxis()->SetTitle("Efficiency");
  h_psEff_dxy_noCuts->Divide(h_trk_matchtp_dxy_noCuts, h_tp_dxy_noCuts, 1.0, 1.0, "B");
  h_psEff_dxy_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_psEff_dxy_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_psEff_dxy_noCuts->GetName() + ".pdf");
  delete h_psEff_dxy_noCuts;
  delete h_trk_matchtp_dxy_noCuts;
  delete h_tp_dxy_noCuts;

  h_trk_matchtp_pt_oldCuts->Sumw2();
  //h_tp_pt_maxEtaCut->Sumw2();
  TH1F* h_psEff_pt_oldCuts = (TH1F*)h_trk_matchtp_pt_oldCuts->Clone();
  h_psEff_pt_oldCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_psEff_pt_oldCuts);
  h_psEff_pt_oldCuts->SetName("psEff_pt_oldCuts");
  h_psEff_pt_oldCuts->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
  h_psEff_pt_oldCuts->GetYaxis()->SetTitle("Efficiency");
  h_psEff_pt_oldCuts->Divide(h_trk_matchtp_pt_oldCuts, h_tp_pt_maxEtaCut, 1.0, 1.0, "B");
  h_psEff_pt_oldCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_psEff_pt_oldCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_psEff_pt_oldCuts->GetName() + ".pdf");
  delete h_psEff_pt_oldCuts;
  delete h_trk_matchtp_pt_oldCuts;
  delete h_tp_pt_maxEtaCut;

  h_trk_matchtp_eta_oldCuts->Sumw2();
  //h_tp_eta_maxEtaCut->Sumw2();
  TH1F* h_psEff_eta_oldCuts = (TH1F*)h_trk_matchtp_eta_oldCuts->Clone();
  h_psEff_eta_oldCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_psEff_eta_oldCuts);
  h_psEff_eta_oldCuts->SetName("psEff_eta_oldCuts");
  h_psEff_eta_oldCuts->GetXaxis()->SetTitle("Tracking Particle #eta (GeV)");
  h_psEff_eta_oldCuts->GetYaxis()->SetTitle("Efficiency");
  h_psEff_eta_oldCuts->Divide(h_trk_matchtp_eta_oldCuts, h_tp_eta_maxEtaCut, 1.0, 1.0, "B");
  h_psEff_eta_oldCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_psEff_eta_oldCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_psEff_eta_oldCuts->GetName() + ".pdf");
  delete h_psEff_eta_oldCuts;
  delete h_trk_matchtp_eta_oldCuts;
  delete h_tp_eta_maxEtaCut;

  h_trk_matchtp_dxy_oldCuts->Sumw2();
  h_tp_dxy_oldCuts->Sumw2();
  TH1F* h_psEff_dxy_oldCuts = (TH1F*)h_trk_matchtp_dxy_oldCuts->Clone();
  h_psEff_dxy_oldCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_psEff_dxy_oldCuts);
  h_psEff_dxy_oldCuts->SetName("psEff_dxy_oldCuts");
  h_psEff_dxy_oldCuts->GetXaxis()->SetTitle("Tracking Particle #dxy (GeV)");
  h_psEff_dxy_oldCuts->GetYaxis()->SetTitle("Efficiency");
  h_psEff_dxy_oldCuts->Divide(h_trk_matchtp_dxy_oldCuts, h_tp_dxy_oldCuts, 1.0, 1.0, "B");
  h_psEff_dxy_oldCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_psEff_dxy_oldCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_psEff_dxy_oldCuts->GetName() + ".pdf");
  delete h_psEff_dxy_oldCuts;
  delete h_trk_matchtp_dxy_oldCuts;
  delete h_tp_dxy_oldCuts;

  h_match_tp_pt_noCuts_primary->Sumw2();
  TH1F* h_primaryEff_pt_noCuts = (TH1F*)h_match_tp_pt_noCuts_primary->Clone();
  h_primaryEff_pt_noCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_primaryEff_pt_noCuts);
  h_primaryEff_pt_noCuts->SetName("primaryEff_pt_noCuts");
  h_primaryEff_pt_noCuts->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
  h_primaryEff_pt_noCuts->GetYaxis()->SetTitle("Efficiency");
  h_primaryEff_pt_noCuts->Divide(h_match_tp_pt_noCuts_primary, h_tp_pt_noCuts_primary, 1.0, 1.0, "B");
  h_primaryEff_pt_noCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_primaryEff_pt_noCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_primaryEff_pt_noCuts->GetName() + ".pdf");
  delete h_primaryEff_pt_noCuts;
  delete h_match_tp_pt_noCuts_primary;
  
  //h_match_tp_pt_allCuts->Sumw2();
  TH1F* h_primaryEff_pt_allCuts = (TH1F*)h_match_tp_pt_allCuts->Clone();
  h_primaryEff_pt_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_primaryEff_pt_allCuts);
  h_primaryEff_pt_allCuts->SetName("primaryEff_pt_allCuts");
  h_primaryEff_pt_allCuts->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
  h_primaryEff_pt_allCuts->GetYaxis()->SetTitle("Efficiency");
  h_primaryEff_pt_allCuts->Divide(h_match_tp_pt_allCuts, h_tp_pt_noCuts_primary, 1.0, 1.0, "B");
  h_primaryEff_pt_allCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_primaryEff_pt_allCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_primaryEff_pt_allCuts->GetName() + ".pdf");
  delete h_primaryEff_pt_allCuts;
  delete h_match_tp_pt_allCuts;

  h_match_tp_pt_allCuts_trkCuts->Sumw2();
  TH1F* h_primaryEff_pt_allCuts_trkCuts = (TH1F*)h_match_tp_pt_allCuts_trkCuts->Clone();
  h_primaryEff_pt_allCuts_trkCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_primaryEff_pt_allCuts_trkCuts);
  h_primaryEff_pt_allCuts_trkCuts->SetName("primaryEff_pt_allCuts_trkCuts");
  h_primaryEff_pt_allCuts_trkCuts->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
  h_primaryEff_pt_allCuts_trkCuts->GetYaxis()->SetTitle("Efficiency");
  h_primaryEff_pt_allCuts_trkCuts->Divide(h_match_tp_pt_allCuts_trkCuts, h_tp_pt_noCuts_primary, 1.0, 1.0, "B");
  h_primaryEff_pt_allCuts_trkCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_primaryEff_pt_allCuts_trkCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_primaryEff_pt_allCuts_trkCuts->GetName() + ".pdf");
  delete h_primaryEff_pt_allCuts_trkCuts;
  delete h_match_tp_pt_allCuts_trkCuts;
  
  //h_match_tp_pt_maxEtaCut->Sumw2();
  TH1F* h_primaryEff_pt_oldCuts = (TH1F*)h_match_tp_pt_maxEtaCut->Clone();
  h_primaryEff_pt_oldCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_primaryEff_pt_oldCuts);
  h_primaryEff_pt_oldCuts->SetName("primaryEff_pt_oldCuts");
  h_primaryEff_pt_oldCuts->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
  h_primaryEff_pt_oldCuts->GetYaxis()->SetTitle("Efficiency");
  h_primaryEff_pt_oldCuts->Divide(h_match_tp_pt_maxEtaCut, h_tp_pt_noCuts_primary, 1.0, 1.0, "B");
  h_primaryEff_pt_oldCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_primaryEff_pt_oldCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_primaryEff_pt_oldCuts->GetName() + ".pdf");
  delete h_primaryEff_pt_oldCuts;
  delete h_match_tp_pt_maxEtaCut;
  delete h_tp_pt_noCuts_primary;
  
  h_trk_pt_oldCuts->Sumw2();
  h_trk_pt_noCuts->Sumw2();
  TH1F* h_trkEff_pt_oldCuts = (TH1F*)h_trk_pt_oldCuts->Clone();
  h_trkEff_pt_oldCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trkEff_pt_oldCuts);
  h_trkEff_pt_oldCuts->SetName("trkEff_pt_oldCuts");
  h_trkEff_pt_oldCuts->GetXaxis()->SetTitle("Track p_{T} (GeV)");
  h_trkEff_pt_oldCuts->GetYaxis()->SetTitle("Efficiency");
  h_trkEff_pt_oldCuts->Divide(h_trk_pt_oldCuts, h_trk_pt_noCuts, 1.0, 1.0, "B");
  h_trkEff_pt_oldCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trkEff_pt_oldCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trkEff_pt_oldCuts->GetName() + ".pdf");
  delete h_trkEff_pt_oldCuts;
  delete h_trk_pt_oldCuts;

  h_trk_pt_allCuts->Sumw2();
  TH1F* h_trkEff_pt_allCuts = (TH1F*)h_trk_pt_allCuts->Clone();
  h_trkEff_pt_allCuts->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trkEff_pt_allCuts);
  h_trkEff_pt_allCuts->SetName("trkEff_pt_allCuts");
  h_trkEff_pt_allCuts->GetXaxis()->SetTitle("Track p_{T} (GeV)");
  h_trkEff_pt_allCuts->GetYaxis()->SetTitle("Efficiency");
  h_trkEff_pt_allCuts->Divide(h_trk_pt_allCuts, h_trk_pt_noCuts, 1.0, 1.0, "B");
  h_trkEff_pt_allCuts->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trkEff_pt_allCuts->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trkEff_pt_allCuts->GetName() + ".pdf");
  delete h_trkEff_pt_allCuts;
  delete h_trk_pt_allCuts;
  delete h_trk_pt_noCuts;

  h_delta_dist_xy->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_delta_dist_xy);
  h_delta_dist_xy->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_delta_dist_xy->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_delta_dist_xy->GetName() + ".pdf");
  delete h_delta_dist_xy;

  h_trueVertexAssoc_delxy->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertexAssoc_delxy);
  h_trueVertexAssoc_delxy->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertexAssoc_delxy->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertexAssoc_delxy->GetName() + ".pdf");
  delete h_trueVertexAssoc_delxy;

  h_trueVertexAssoc_delz->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertexAssoc_delz);
  h_trueVertexAssoc_delz->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertexAssoc_delz->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertexAssoc_delz->GetName() + ".pdf");
  delete h_trueVertexAssoc_delz;

  h_trueVertexAssoc_calcVertDelxy->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertexAssoc_calcVertDelxy);
  h_trueVertexAssoc_calcVertDelxy->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertexAssoc_calcVertDelxy->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertexAssoc_calcVertDelxy->GetName() + ".pdf");
  delete h_trueVertexAssoc_calcVertDelxy;

  h_trueVertexAssoc_calcVertDelz->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertexAssoc_calcVertDelz);
  h_trueVertexAssoc_calcVertDelz->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertexAssoc_calcVertDelz->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertexAssoc_calcVertDelz->GetName() + ".pdf");
  delete h_trueVertexAssoc_calcVertDelz;

  h_trackVertexAssoc_delxy->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trackVertexAssoc_delxy);
  h_trackVertexAssoc_delxy->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trackVertexAssoc_delxy->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trackVertexAssoc_delxy->GetName() + ".pdf");
  delete h_trackVertexAssoc_delxy;

  h_trackVertexAssoc_delz->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trackVertexAssoc_delz);
  h_trackVertexAssoc_delz->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trackVertexAssoc_delz->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trackVertexAssoc_delz->GetName() + ".pdf");
  delete h_trackVertexAssoc_delz;
   
  h_error_delta_x->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_error_delta_x);
  h_error_delta_x->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_error_delta_x->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_error_delta_x->GetName() + ".pdf");
  delete h_error_delta_x;
   
  h_delta_dist_z->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_delta_dist_z);
  h_delta_dist_z->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_delta_dist_z->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_delta_dist_z->GetName() + ".pdf");
  delete h_delta_dist_z;

  h_error_delta_z->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_error_delta_z);
  h_error_delta_z->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_error_delta_z->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_error_delta_z->GetName() + ".pdf");
  delete h_error_delta_z;
   
  h_delta_x->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_delta_x);
  h_delta_x->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_delta_x->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_delta_x->GetName() + ".pdf");
  delete h_delta_x;

  char binlabel[1000];
  sprintf(binlabel, "0 + 1,   0 + 2,   1 + 2,   0 + 3,   1 + 3,   2 + 3");

  h_trk_Counter_TPcombination->GetYaxis()->SetNoExponent(kTRUE);
  h_trk_Counter_TPcombination->GetXaxis()->SetRange(0, h_trk_Counter_TPcombination->GetNbinsX() + 1);
  h_trk_Counter_TPcombination->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  mySmallText(0.25, 0.95, 1, binlabel);
  h_trk_Counter_TPcombination->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_Counter_TPcombination->GetName() + ".pdf");
  delete h_trk_Counter_TPcombination;

  h_trk_delta_dist_xy->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_delta_dist_xy);
  h_trk_delta_dist_xy->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_delta_dist_xy->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_delta_dist_xy->GetName() + ".pdf");
  delete h_trk_delta_dist_xy;

  h_trackVertex_cos_T->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trackVertex_cos_T);
  h_trackVertex_cos_T->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trackVertex_cos_T->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trackVertex_cos_T->GetName() + ".pdf");
  delete h_trackVertex_cos_T;

  h_trackVertex_alpha_T->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trackVertex_alpha_T);
  h_trackVertex_alpha_T->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trackVertex_alpha_T->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trackVertex_alpha_T->GetName() + ".pdf");
  delete h_trackVertex_alpha_T;

  h_trackVertex_d_T->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trackVertex_d_T);
  h_trackVertex_d_T->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trackVertex_d_T->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trackVertex_d_T->GetName() + ".pdf");
  delete h_trackVertex_d_T;

  h_trackVertex_openingAngle->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trackVertex_openingAngle);
  h_trackVertex_openingAngle->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trackVertex_openingAngle->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trackVertex_openingAngle->GetName() + ".pdf");
  delete h_trackVertex_openingAngle;

  h_trackVertex_parentPt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trackVertex_parentPt);
  h_trackVertex_parentPt->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trackVertex_parentPt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trackVertex_parentPt->GetName() + ".pdf");
  delete h_trackVertex_parentPt;

  h_trackVertex_R_T->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trackVertex_R_T);
  h_trackVertex_R_T->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trackVertex_R_T->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trackVertex_R_T->GetName() + ".pdf");
  delete h_trackVertex_R_T;

  h_trueVertex_cos_T->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertex_cos_T);
  h_trueVertex_cos_T->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertex_cos_T->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertex_cos_T->GetName() + ".pdf");
  delete h_trueVertex_cos_T;

  h_trueVertex_alpha_T->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertex_alpha_T);
  h_trueVertex_alpha_T->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertex_alpha_T->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertex_alpha_T->GetName() + ".pdf");
  delete h_trueVertex_alpha_T;

  h_trueVertex_d_T->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertex_d_T);
  h_trueVertex_d_T->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertex_d_T->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertex_d_T->GetName() + ".pdf");
  delete h_trueVertex_d_T;

  h_trueVertex_openingAngle->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertex_openingAngle);
  h_trueVertex_openingAngle->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertex_openingAngle->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertex_openingAngle->GetName() + ".pdf");
  delete h_trueVertex_openingAngle;

  h_trueVertex_parentPt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertex_parentPt);
  h_trueVertex_parentPt->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertex_parentPt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertex_parentPt->GetName() + ".pdf");
  delete h_trueVertex_parentPt;
   
  h_trueVertex_R_T->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertex_R_T);
  h_trueVertex_R_T->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertex_R_T->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertex_R_T->GetName() + ".pdf");
  delete h_trueVertex_R_T;

  h_trueVertex_delta_dist_z0->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertex_delta_dist_z0);
  h_trueVertex_delta_dist_z0->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertex_delta_dist_z0->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertex_delta_dist_z0->GetName() + ".pdf");
  delete h_trueVertex_delta_dist_z0;

  h_trueVertex_delta_dist_d0->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertex_delta_dist_d0);
  h_trueVertex_delta_dist_d0->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertex_delta_dist_d0->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertex_delta_dist_d0->GetName() + ".pdf");
  delete h_trueVertex_delta_dist_d0;

  h_trueVertex_delta_dist_eta->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertex_delta_dist_eta);
  h_trueVertex_delta_dist_eta->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertex_delta_dist_eta->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertex_delta_dist_eta->GetName() + ".pdf");
  delete h_trueVertex_delta_dist_eta;

  h_trueVertex_delta_dist_phi->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertex_delta_dist_phi);
  h_trueVertex_delta_dist_phi->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertex_delta_dist_phi->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertex_delta_dist_phi->GetName() + ".pdf");
  delete h_trueVertex_delta_dist_phi;

  h_trueVertex_delta_dist_R->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertex_delta_dist_R);
  h_trueVertex_delta_dist_R->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertex_delta_dist_R->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertex_delta_dist_R->GetName() + ".pdf");
  delete h_trueVertex_delta_dist_R;

  h_trueVertex_delta_dist_indexPt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertex_delta_dist_indexPt);
  h_trueVertex_delta_dist_indexPt->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertex_delta_dist_indexPt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertex_delta_dist_indexPt->GetName() + ".pdf");
  delete h_trueVertex_delta_dist_indexPt;

  h_trueVertex_delta_dist_pt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trueVertex_delta_dist_pt);
  h_trueVertex_delta_dist_pt->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertex_delta_dist_pt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertex_delta_dist_pt->GetName() + ".pdf");
  delete h_trueVertex_delta_dist_pt;

  h_trk_delta_dist_z->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_trk_delta_dist_z);
  h_trk_delta_dist_z->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trk_delta_dist_z->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trk_delta_dist_z->GetName() + ".pdf");
  delete h_trk_delta_dist_z;

  char res[1000];
  float rms = 0;
  TF1* fit;
  fit = new TF1("fit", "gaus", -1, 1);
  removeFlows(h_res_tp_trk_x);
  h_res_tp_trk_x->Fit("fit","R");
  h_res_tp_trk_x->SetStats(0);
  h_res_tp_trk_x->Draw();
  rms = fit->GetParameter(2);
  sprintf(res, "RMS = %.4f", rms);
  mySmallText(0.22, 0.82, 1, res);
  mySmallText(0.4, 0.42, 1, ctxt);
  h_res_tp_trk_x->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_res_tp_trk_x->GetName() + ".pdf");
  delete h_res_tp_trk_x;
  delete fit;

  fit = new TF1("fit", "gaus", -1, 1);
  removeFlows(h_res_tp_trk_x_zoomOut);
  h_res_tp_trk_x_zoomOut->Fit("fit","R");
  h_res_tp_trk_x_zoomOut->Draw();
  rms = fit->GetParameter(2);
  sprintf(res, "RMS = %.4f", rms);
  mySmallText(0.22, 0.82, 1, res);
  mySmallText(0.4, 0.42, 1, ctxt);
  h_res_tp_trk_x_zoomOut->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_res_tp_trk_x_zoomOut->GetName() + ".pdf");
  delete h_res_tp_trk_x_zoomOut;
  delete fit;

  fit = new TF1("fit", "gaus", -1, 1);
  removeFlows(h_res_tp_trk_x_findVert);
  h_res_tp_trk_x_findVert->Fit("fit","R");
  h_res_tp_trk_x_findVert->Draw();
  rms = fit->GetParameter(2);
  sprintf(res, "RMS = %.4f", rms);
  mySmallText(0.22, 0.82, 1, res);
  mySmallText(0.4, 0.42, 1, ctxt);
  h_res_tp_trk_x_findVert->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_res_tp_trk_x_findVert->GetName() + ".pdf");
  delete h_res_tp_trk_x_findVert;
  delete fit;

  fit = new TF1("fit", "gaus", -1, 1);
  removeFlows(h_res_tp_trk_y);
  h_res_tp_trk_y->Fit("fit","R");
  h_res_tp_trk_y->SetStats(0);
  h_res_tp_trk_y->Draw();
  rms = fit->GetParameter(2);
  sprintf(res, "RMS = %.4f", rms);
  mySmallText(0.22, 0.82, 1, res);
  mySmallText(0.4, 0.42, 1, ctxt);
  h_res_tp_trk_y->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_res_tp_trk_y->GetName() + ".pdf");
  delete h_res_tp_trk_y;
  delete fit;

  fit = new TF1("fit", "gaus", -1, 1);
  removeFlows(h_res_tp_trk_y_zoomOut);
  h_res_tp_trk_y_zoomOut->Fit("fit","R");
  h_res_tp_trk_y_zoomOut->Draw();
  rms = fit->GetParameter(2);
  sprintf(res, "RMS = %.4f", rms);
  mySmallText(0.22, 0.82, 1, res);
  mySmallText(0.4, 0.42, 1, ctxt);
  h_res_tp_trk_y_zoomOut->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_res_tp_trk_y_zoomOut->GetName() + ".pdf");
  delete h_res_tp_trk_y_zoomOut;
  delete fit;

  fit = new TF1("fit", "gaus", -1, 1);
  removeFlows(h_res_tp_trk_y_findVert);
  h_res_tp_trk_y_findVert->Fit("fit","R");
  h_res_tp_trk_y_findVert->Draw();
  rms = fit->GetParameter(2);
  sprintf(res, "RMS = %.4f", rms);
  mySmallText(0.22, 0.82, 1, res);
  mySmallText(0.4, 0.42, 1, ctxt);
  h_res_tp_trk_y_findVert->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_res_tp_trk_y_findVert->GetName() + ".pdf");
  delete h_res_tp_trk_y_findVert;
  delete fit;

  fit = new TF1("fit", "gaus", -1, 1);
  removeFlows(h_res_tp_trk_r);
  h_res_tp_trk_r->Fit("fit");
  h_res_tp_trk_r->Draw();
  rms = fit->GetParameter(2);
  sprintf(res, "RMS = %.4f", rms);
  mySmallText(0.22, 0.82, 1, res);
  mySmallText(0.4, 0.42, 1, ctxt);
  h_res_tp_trk_r->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_res_tp_trk_r->GetName() + ".pdf");
  delete h_res_tp_trk_r;
  delete fit;

  fit = new TF1("fit", "gaus", -1, 1);
  removeFlows(h_res_tp_trk_phi);
  h_res_tp_trk_phi->Fit("fit");
  h_res_tp_trk_phi->Draw();
  rms = fit->GetParameter(2);
  sprintf(res, "RMS = %.4f", rms);
  mySmallText(0.22, 0.82, 1, res);
  mySmallText(0.4, 0.42, 1, ctxt);
  h_res_tp_trk_phi->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_res_tp_trk_phi->GetName() + ".pdf");
  delete h_res_tp_trk_phi;
  delete fit;

  fit = new TF1("fit", "gaus", -10, 10);
  removeFlows(h_res_tp_trk_z);
  h_res_tp_trk_z->Fit("fit");
  h_res_tp_trk_z->SetStats(0);
  h_res_tp_trk_z->Draw();
  rms = fit->GetParameter(2);
  sprintf(res, "RMS = %.4f", rms);
  mySmallText(0.22, 0.82, 1, res);
  mySmallText(0.4, 0.42, 1, ctxt);
  h_res_tp_trk_z->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_res_tp_trk_z->GetName() + ".pdf");
  delete h_res_tp_trk_z;
  delete fit;

  removeFlows(h_all_trueVertex_pt);
  h_all_trueVertex_pt->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_all_trueVertex_pt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_all_trueVertex_pt->GetName() + ".pdf");

  removeFlows(h_all_trueVertex_minD0);
  h_all_trueVertex_minD0->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_all_trueVertex_minD0->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_all_trueVertex_minD0->GetName() + ".pdf");

  removeFlows(h_all_trueVertex_maxD0);
  h_all_trueVertex_maxD0->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_all_trueVertex_maxD0->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_all_trueVertex_maxD0->GetName() + ".pdf");

  removeFlows(h_all_trueVertex_minD0_allTPs);
  h_all_trueVertex_minD0_allTPs->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_all_trueVertex_minD0_allTPs->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_all_trueVertex_minD0_allTPs->GetName() + ".pdf");

  removeFlows(h_all_trueVertex_maxD0_allTPs);
  h_all_trueVertex_maxD0_allTPs->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_all_trueVertex_maxD0_allTPs->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_all_trueVertex_maxD0_allTPs->GetName() + ".pdf");

  removeFlows(h_all_trueVertex_lowPt);
  h_all_trueVertex_lowPt->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_all_trueVertex_lowPt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_all_trueVertex_lowPt->GetName() + ".pdf");

  removeFlows(h_all_trueVertex_lowPt_allTPs);
  h_all_trueVertex_lowPt_allTPs->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_all_trueVertex_lowPt_allTPs->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_all_trueVertex_lowPt_allTPs->GetName() + ".pdf");

  removeFlows(h_trueVertex_numTPs);
  h_trueVertex_numTPs->SetStats(0);
  h_trueVertex_numTPs->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertex_numTPs->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertex_numTPs->GetName() + ".pdf");

  removeFlows(h_trueVertex_numTracks);
  h_trueVertex_numTracks->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trueVertex_numTracks->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trueVertex_numTracks->GetName() + ".pdf");

  removeFlows(h_trackVertex_numTracks);
  h_trackVertex_numTracks->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_trackVertex_numTracks->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_trackVertex_numTracks->GetName() + ".pdf");

  removeFlows(h_correct_trackVertex_numTracks);
  h_correct_trackVertex_numTracks->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_correct_trackVertex_numTracks->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_correct_trackVertex_numTracks->GetName() + ".pdf");

  removeFlows(h_false_trackVertex_numTracks);
  h_false_trackVertex_numTracks->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_false_trackVertex_numTracks->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_false_trackVertex_numTracks->GetName() + ".pdf");

  removeFlows(h_correct_trackVertex_inTraj);
  h_correct_trackVertex_inTraj->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_correct_trackVertex_inTraj->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_correct_trackVertex_inTraj->GetName() + ".pdf");

  removeFlows(h_false_trackVertex_inTraj);
  h_false_trackVertex_inTraj->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_false_trackVertex_inTraj->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_false_trackVertex_inTraj->GetName() + ".pdf");

  removeFlows(h_false_trackVertex_numMatched);
  h_false_trackVertex_numMatched->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_false_trackVertex_numMatched->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_false_trackVertex_numMatched->GetName() + ".pdf");

  removeFlows(h_all_trackVertex_pt);
  h_all_trackVertex_pt->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_all_trackVertex_pt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_all_trackVertex_pt->GetName() + ".pdf");

  removeFlows(h_all_trackVertex_minD0);
  h_all_trackVertex_minD0->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_all_trackVertex_minD0->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_all_trackVertex_minD0->GetName() + ".pdf");
   
  removeFlows(h_all_trackVertex_lowPt);
  h_all_trackVertex_lowPt->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_all_trackVertex_lowPt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_all_trackVertex_lowPt->GetName() + ".pdf");

  removeFlows(h_correct_trueVertex_pt);
  h_correct_trueVertex_pt->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_correct_trueVertex_pt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_correct_trueVertex_pt->GetName() + ".pdf");
   
  h_all_trueVertex_pt->Sumw2();
  h_correct_trueVertex_pt->Sumw2();
  TH1F* h_eff_trueVertex_pt = (TH1F*)h_correct_trueVertex_pt->Clone();
  h_eff_trueVertex_pt->GetYaxis()->SetNoExponent(kTRUE);
  h_eff_trueVertex_pt->SetName("eff_trueVertex_pt");
  h_eff_trueVertex_pt->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
  h_eff_trueVertex_pt->GetYaxis()->SetTitle("Efficiency");
  h_eff_trueVertex_pt->Divide(h_correct_trueVertex_pt,h_all_trueVertex_pt, 1.0, 1.0, "B");
  h_eff_trueVertex_pt->SetStats(0);
  h_eff_trueVertex_pt->Draw();
  h_eff_trueVertex_pt->SetAxisRange(0, 1.1, "Y");
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_trueVertex_pt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_trueVertex_pt->GetName() + ".pdf");
  delete h_eff_trueVertex_pt;
  delete h_correct_trueVertex_pt;

  h_findable_trueVertex_pt->Sumw2();
  h_findable_trueVertexBinned_pt->Sumw2();
  h_findableIntersect_trueVertexBinned_pt->Sumw2();
  h_findableNonPrompt_trueVertexBinned_pt->Sumw2();
  h_findableBeforeTracker_trueVertexBinned_pt->Sumw2();
  h_findableBendChi2Cut_trueVertexBinned_pt->Sumw2();
  h_findableChi2RPhiCut_trueVertexBinned_pt->Sumw2();
  h_findableCosTCut_trueVertexBinned_pt->Sumw2();
  h_findableDeltaZCut_trueVertexBinned_pt->Sumw2();
  TH1F* h_findEff_trueVertex_pt = (TH1F*)h_findable_trueVertex_pt->Clone();
  h_findEff_trueVertex_pt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_findEff_trueVertex_pt);
  h_findEff_trueVertex_pt->SetName("findEff_trueVertex_pt");
  h_findEff_trueVertex_pt->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
  h_findEff_trueVertex_pt->GetYaxis()->SetTitle("Findable Efficiency");
  h_findEff_trueVertex_pt->Divide(h_findable_trueVertex_pt,h_all_trueVertex_pt, 1.0, 1.0, "B");
  h_findEff_trueVertex_pt->SetLineColor(1);
  h_findEff_trueVertex_pt->SetMarkerColor(1);
  h_findEff_trueVertex_pt->Draw();
  h_findEff_trueVertex_pt->SetAxisRange(0, 1.1, "Y");
  h_findable_trueVertexBinned_pt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_findable_trueVertexBinned_pt);
  h_findable_trueVertexBinned_pt->Divide(h_findable_trueVertexBinned_pt,h_all_trueVertex_pt, 1.0, 1.0, "B");
  h_findable_trueVertexBinned_pt->SetLineColor(2);
  h_findable_trueVertexBinned_pt->SetMarkerColor(2);
  h_findable_trueVertexBinned_pt->Draw("SAME");
  h_findableIntersect_trueVertexBinned_pt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_findableIntersect_trueVertexBinned_pt);
  h_findableIntersect_trueVertexBinned_pt->Divide(h_findableIntersect_trueVertexBinned_pt,h_all_trueVertex_pt, 1.0, 1.0, "B");
  h_findableIntersect_trueVertexBinned_pt->SetLineColor(3);
  h_findableIntersect_trueVertexBinned_pt->SetMarkerColor(3);
  h_findableIntersect_trueVertexBinned_pt->Draw("SAME");
  h_findableNonPrompt_trueVertexBinned_pt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_findableNonPrompt_trueVertexBinned_pt);
  h_findableNonPrompt_trueVertexBinned_pt->Divide(h_findableNonPrompt_trueVertexBinned_pt,h_all_trueVertex_pt, 1.0, 1.0, "B");
  h_findableNonPrompt_trueVertexBinned_pt->SetLineColor(4);
  h_findableNonPrompt_trueVertexBinned_pt->SetMarkerColor(4);
  h_findableNonPrompt_trueVertexBinned_pt->Draw("SAME");
  h_findableBeforeTracker_trueVertexBinned_pt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_findableBeforeTracker_trueVertexBinned_pt);
  h_findableBeforeTracker_trueVertexBinned_pt->Divide(h_findableBeforeTracker_trueVertexBinned_pt,h_all_trueVertex_pt, 1.0, 1.0, "B");
  h_findableBeforeTracker_trueVertexBinned_pt->SetLineColor(5);
  h_findableBeforeTracker_trueVertexBinned_pt->SetMarkerColor(5);
  h_findableBeforeTracker_trueVertexBinned_pt->Draw("SAME");
  h_findableBendChi2Cut_trueVertexBinned_pt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_findableBendChi2Cut_trueVertexBinned_pt);
  h_findableBendChi2Cut_trueVertexBinned_pt->Divide(h_findableBendChi2Cut_trueVertexBinned_pt,h_all_trueVertex_pt, 1.0, 1.0, "B");
  h_findableBendChi2Cut_trueVertexBinned_pt->SetLineColor(6);
  h_findableBendChi2Cut_trueVertexBinned_pt->SetMarkerColor(6);
  h_findableBendChi2Cut_trueVertexBinned_pt->Draw("SAME");
  h_findableChi2RPhiCut_trueVertexBinned_pt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_findableChi2RPhiCut_trueVertexBinned_pt);
  h_findableChi2RPhiCut_trueVertexBinned_pt->Divide(h_findableChi2RPhiCut_trueVertexBinned_pt,h_all_trueVertex_pt, 1.0, 1.0, "B");
  h_findableChi2RPhiCut_trueVertexBinned_pt->SetLineColor(7);
  h_findableChi2RPhiCut_trueVertexBinned_pt->SetMarkerColor(7);
  h_findableChi2RPhiCut_trueVertexBinned_pt->Draw("SAME");
  h_findableCosTCut_trueVertexBinned_pt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_findableCosTCut_trueVertexBinned_pt);
  h_findableCosTCut_trueVertexBinned_pt->Divide(h_findableCosTCut_trueVertexBinned_pt,h_all_trueVertex_pt, 1.0, 1.0, "B");
  h_findableCosTCut_trueVertexBinned_pt->SetLineColor(8);
  h_findableCosTCut_trueVertexBinned_pt->SetMarkerColor(8);
  h_findableCosTCut_trueVertexBinned_pt->Draw("SAME");
  h_findableDeltaZCut_trueVertexBinned_pt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_findableDeltaZCut_trueVertexBinned_pt);
  h_findableDeltaZCut_trueVertexBinned_pt->Divide(h_findableDeltaZCut_trueVertexBinned_pt,h_all_trueVertex_pt, 1.0, 1.0, "B");
  h_findableDeltaZCut_trueVertexBinned_pt->SetLineColor(9);
  h_findableDeltaZCut_trueVertexBinned_pt->SetMarkerColor(9);
  h_findableDeltaZCut_trueVertexBinned_pt->Draw("SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_findEff_trueVertex_pt,"Pre-sel","lp");
  l->AddEntry(h_findable_trueVertexBinned_pt,"Binned","lp");
  l->AddEntry(h_findableIntersect_trueVertexBinned_pt,"#Delta_{xy}=0","lp");
  l->AddEntry(h_findableNonPrompt_trueVertexBinned_pt,"R>#sigma_{d0}","lp");
  l->AddEntry(h_findableBeforeTracker_trueVertexBinned_pt,"R<20cm","lp");
  l->AddEntry(h_findableBendChi2Cut_trueVertexBinned_pt,"#Sigma#chi^{2}_{bend}<12","lp");
  l->AddEntry(h_findableChi2RPhiCut_trueVertexBinned_pt,"#Sigma #chi^{2}_{r#phi}/d.o.f<6.0","lp");
  l->AddEntry(h_findableCosTCut_trueVertexBinned_pt,"cos_T>0","lp");
  l->AddEntry(h_findableDeltaZCut_trueVertexBinned_pt,"#Deltaz<1cm","lp");
  l->Draw();
  c.SaveAs(DIR + "/"+ h_findEff_trueVertex_pt->GetName() + ".pdf");
  delete h_findEff_trueVertex_pt;
  delete h_all_trueVertex_pt;
  delete h_findable_trueVertex_pt;
  delete h_findable_trueVertexBinned_pt;
  delete h_findableIntersect_trueVertexBinned_pt;
  delete h_findableNonPrompt_trueVertexBinned_pt;
  delete h_findableBeforeTracker_trueVertexBinned_pt;
  delete h_findableBendChi2Cut_trueVertexBinned_pt;
  delete h_findableChi2RPhiCut_trueVertexBinned_pt;
  delete h_findableCosTCut_trueVertexBinned_pt;
  delete h_findableDeltaZCut_trueVertexBinned_pt;

  h_findFake_trackVertex_pt->Sumw2();
  h_findFake_trackVertexBinned_pt->Sumw2();
  h_findFakeIntersect_trackVertexBinned_pt->Sumw2();
  h_findFakeNonPrompt_trackVertexBinned_pt->Sumw2();
  h_findFakeBeforeTracker_trackVertexBinned_pt->Sumw2();
  h_findFakeBendChi2Cut_trackVertexBinned_pt->Sumw2();
  h_findFakeChi2RPhiCut_trackVertexBinned_pt->Sumw2();
  h_findFakeCosTCut_trackVertexBinned_pt->Sumw2();
  h_findFakeDeltaZCut_trackVertexBinned_pt->Sumw2();
  TH1F* h_findFakeEff_trackVertex_pt = (TH1F*)h_findFake_trackVertex_pt->Clone();
  h_findFakeEff_trackVertex_pt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_findFakeEff_trackVertex_pt);
  h_findFakeEff_trackVertex_pt->SetName("findFakeEff_trackVertex_pt");
  h_findFakeEff_trackVertex_pt->GetXaxis()->SetTitle("Track p_{T} (GeV)");
  h_findFakeEff_trackVertex_pt->GetYaxis()->SetTitle("FindFake Efficiency");
  h_findFakeEff_trackVertex_pt->Divide(h_findFake_trackVertexBinned_pt,h_findFake_trackVertex_pt, 1.0, 1.0, "B");
  h_findFakeEff_trackVertex_pt->SetLineColor(2);
  h_findFakeEff_trackVertex_pt->SetMarkerColor(2);
  h_findFakeEff_trackVertex_pt->Draw();
  h_findFakeIntersect_trackVertexBinned_pt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_findFakeIntersect_trackVertexBinned_pt);
  h_findFakeIntersect_trackVertexBinned_pt->Divide(h_findFakeIntersect_trackVertexBinned_pt,h_findFake_trackVertex_pt, 1.0, 1.0, "B");
  h_findFakeIntersect_trackVertexBinned_pt->SetLineColor(3);
  h_findFakeIntersect_trackVertexBinned_pt->SetMarkerColor(3);
  h_findFakeIntersect_trackVertexBinned_pt->Draw("SAME");
  h_findFakeNonPrompt_trackVertexBinned_pt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_findFakeNonPrompt_trackVertexBinned_pt);
  h_findFakeNonPrompt_trackVertexBinned_pt->Divide(h_findFakeNonPrompt_trackVertexBinned_pt,h_findFake_trackVertex_pt, 1.0, 1.0, "B");
  h_findFakeNonPrompt_trackVertexBinned_pt->SetLineColor(4);
  h_findFakeNonPrompt_trackVertexBinned_pt->SetMarkerColor(4);
  h_findFakeNonPrompt_trackVertexBinned_pt->Draw("SAME");
  h_findFakeBeforeTracker_trackVertexBinned_pt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_findFakeBeforeTracker_trackVertexBinned_pt);
  h_findFakeBeforeTracker_trackVertexBinned_pt->Divide(h_findFakeBeforeTracker_trackVertexBinned_pt,h_findFake_trackVertex_pt, 1.0, 1.0, "B");
  h_findFakeBeforeTracker_trackVertexBinned_pt->SetLineColor(5);
  h_findFakeBeforeTracker_trackVertexBinned_pt->SetMarkerColor(5);
  h_findFakeBeforeTracker_trackVertexBinned_pt->Draw("SAME");
  h_findFakeBendChi2Cut_trackVertexBinned_pt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_findFakeBendChi2Cut_trackVertexBinned_pt);
  h_findFakeBendChi2Cut_trackVertexBinned_pt->Divide(h_findFakeBendChi2Cut_trackVertexBinned_pt,h_findFake_trackVertex_pt, 1.0, 1.0, "B");
  h_findFakeBendChi2Cut_trackVertexBinned_pt->SetLineColor(6);
  h_findFakeBendChi2Cut_trackVertexBinned_pt->SetMarkerColor(6);
  h_findFakeBendChi2Cut_trackVertexBinned_pt->Draw("SAME");
  h_findFakeChi2RPhiCut_trackVertexBinned_pt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_findFakeChi2RPhiCut_trackVertexBinned_pt);
  h_findFakeChi2RPhiCut_trackVertexBinned_pt->Divide(h_findFakeChi2RPhiCut_trackVertexBinned_pt,h_findFake_trackVertex_pt, 1.0, 1.0, "B");
  h_findFakeChi2RPhiCut_trackVertexBinned_pt->SetLineColor(7);
  h_findFakeChi2RPhiCut_trackVertexBinned_pt->SetMarkerColor(7);
  h_findFakeChi2RPhiCut_trackVertexBinned_pt->Draw("SAME");
  h_findFakeCosTCut_trackVertexBinned_pt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_findFakeCosTCut_trackVertexBinned_pt);
  h_findFakeCosTCut_trackVertexBinned_pt->Divide(h_findFakeCosTCut_trackVertexBinned_pt,h_findFake_trackVertex_pt, 1.0, 1.0, "B");
  h_findFakeCosTCut_trackVertexBinned_pt->SetLineColor(8);
  h_findFakeCosTCut_trackVertexBinned_pt->SetMarkerColor(8);
  h_findFakeCosTCut_trackVertexBinned_pt->Draw("SAME");
  h_findFakeDeltaZCut_trackVertexBinned_pt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_findFakeDeltaZCut_trackVertexBinned_pt);
  h_findFakeDeltaZCut_trackVertexBinned_pt->Divide(h_findFakeDeltaZCut_trackVertexBinned_pt,h_findFake_trackVertex_pt, 1.0, 1.0, "B");
  h_findFakeDeltaZCut_trackVertexBinned_pt->SetLineColor(9);
  h_findFakeDeltaZCut_trackVertexBinned_pt->SetMarkerColor(9);
  h_findFakeDeltaZCut_trackVertexBinned_pt->Draw("SAME");
  mySmallText(0.4, 0.82, 1, ctxt);
  l->Clear();
  l->AddEntry(h_findFakeEff_trackVertex_pt,"Binned","lp");
  l->AddEntry(h_findFakeIntersect_trackVertexBinned_pt,"#Delta_{xy}=0","lp");
  l->AddEntry(h_findFakeNonPrompt_trackVertexBinned_pt,"R>#sigma_{d0}","lp");
  l->AddEntry(h_findFakeBeforeTracker_trackVertexBinned_pt,"R<20cm","lp");
  l->AddEntry(h_findFakeBendChi2Cut_trackVertexBinned_pt,"#Sigma#chi^{2}_{bend}<12","lp");
  l->AddEntry(h_findFakeChi2RPhiCut_trackVertexBinned_pt,"#Sigma #chi^{2}_{r#phi}/d.o.f<6.0","lp");
  l->AddEntry(h_findFakeCosTCut_trackVertexBinned_pt,"cos_T>0","lp");
  l->AddEntry(h_findFakeDeltaZCut_trackVertexBinned_pt,"#Deltaz<1cm","lp");
  l->Draw();
  c.SaveAs(DIR + "/"+ h_findFakeEff_trackVertex_pt->GetName() + ".pdf");
  delete h_findFakeEff_trackVertex_pt;
  delete h_findFake_trackVertex_pt;
  delete h_findFake_trackVertexBinned_pt;
  delete h_findFakeIntersect_trackVertexBinned_pt;
  delete h_findFakeNonPrompt_trackVertexBinned_pt;
  delete h_findFakeBeforeTracker_trackVertexBinned_pt;
  delete h_findFakeBendChi2Cut_trackVertexBinned_pt;
  delete h_findFakeChi2RPhiCut_trackVertexBinned_pt;
  delete h_findFakeCosTCut_trackVertexBinned_pt;
  delete h_findFakeDeltaZCut_trackVertexBinned_pt;

  h_all_trueVertex_lowPt->Sumw2();
  h_correct_trueVertex_lowPt->Sumw2();
  TH1F* h_eff_trueVertex_lowPt = (TH1F*)h_correct_trueVertex_lowPt->Clone();
  h_eff_trueVertex_lowPt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_trueVertex_lowPt);
  h_eff_trueVertex_lowPt->SetName("eff_trueVertex_lowPt");
  h_eff_trueVertex_lowPt->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
  h_eff_trueVertex_lowPt->GetYaxis()->SetTitle("Efficiency");
  h_eff_trueVertex_lowPt->Divide(h_correct_trueVertex_lowPt,h_all_trueVertex_lowPt, 1.0, 1.0, "B");
  h_eff_trueVertex_lowPt->Draw();
  h_eff_trueVertex_lowPt->SetAxisRange(0, 1.1, "Y");
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_trueVertex_lowPt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_trueVertex_lowPt->GetName() + ".pdf");
  delete h_eff_trueVertex_lowPt;
  delete h_all_trueVertex_lowPt;
  delete h_correct_trueVertex_lowPt;

  h_all_oneMatch_trueVertex_pt->Sumw2();
  h_correct_oneMatch_trueVertex_pt->Sumw2();
  TH1F* h_eff_oneMatch_trueVertex_pt = (TH1F*)h_correct_oneMatch_trueVertex_pt->Clone();
  h_eff_oneMatch_trueVertex_pt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_oneMatch_trueVertex_pt);
  h_eff_oneMatch_trueVertex_pt->SetName("eff_oneMatch_trueVertex_pt");
  h_eff_oneMatch_trueVertex_pt->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
  h_eff_oneMatch_trueVertex_pt->GetYaxis()->SetTitle("Efficiency");
  h_eff_oneMatch_trueVertex_pt->Divide(h_correct_oneMatch_trueVertex_pt,h_all_oneMatch_trueVertex_pt, 1.0, 1.0, "B");
  h_eff_oneMatch_trueVertex_pt->SetStats(0);
  h_eff_oneMatch_trueVertex_pt->Draw();
  h_eff_oneMatch_trueVertex_pt->SetAxisRange(0, 1.1, "Y");
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_oneMatch_trueVertex_pt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_oneMatch_trueVertex_pt->GetName() + ".pdf");
  delete h_eff_oneMatch_trueVertex_pt;
  delete h_all_oneMatch_trueVertex_pt;
  delete h_correct_oneMatch_trueVertex_pt;

  h_all_oneMatch_trueVertex_lowPt->Sumw2();
  h_correct_oneMatch_trueVertex_lowPt->Sumw2();
  TH1F* h_eff_oneMatch_trueVertex_lowPt = (TH1F*)h_correct_oneMatch_trueVertex_lowPt->Clone();
  h_eff_oneMatch_trueVertex_lowPt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_oneMatch_trueVertex_lowPt);
  h_eff_oneMatch_trueVertex_lowPt->SetName("eff_oneMatch_trueVertex_lowPt");
  h_eff_oneMatch_trueVertex_lowPt->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
  h_eff_oneMatch_trueVertex_lowPt->GetYaxis()->SetTitle("Efficiency");
  h_eff_oneMatch_trueVertex_lowPt->Divide(h_correct_oneMatch_trueVertex_lowPt,h_all_oneMatch_trueVertex_lowPt, 1.0, 1.0, "B");
  h_eff_oneMatch_trueVertex_lowPt->Draw();
  h_eff_oneMatch_trueVertex_lowPt->SetAxisRange(0, 1.1, "Y");
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_oneMatch_trueVertex_lowPt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_oneMatch_trueVertex_lowPt->GetName() + ".pdf");
  delete h_eff_oneMatch_trueVertex_lowPt;
  delete h_all_oneMatch_trueVertex_lowPt;
  delete h_correct_oneMatch_trueVertex_lowPt;

  h_all_oneMatchAlt_trueVertex_pt->Sumw2();
  h_correct_oneMatchAlt_trueVertex_pt->Sumw2();
  TH1F* h_eff_oneMatchAlt_trueVertex_pt = (TH1F*)h_correct_oneMatchAlt_trueVertex_pt->Clone();
  h_eff_oneMatchAlt_trueVertex_pt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_oneMatchAlt_trueVertex_pt);
  h_eff_oneMatchAlt_trueVertex_pt->SetName("eff_oneMatchAlt_trueVertex_pt");
  h_eff_oneMatchAlt_trueVertex_pt->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
  h_eff_oneMatchAlt_trueVertex_pt->GetYaxis()->SetTitle("Efficiency");
  h_eff_oneMatchAlt_trueVertex_pt->Divide(h_correct_oneMatchAlt_trueVertex_pt,h_all_oneMatchAlt_trueVertex_pt, 1.0, 1.0, "B");
  h_eff_oneMatchAlt_trueVertex_pt->Draw();
  h_eff_oneMatchAlt_trueVertex_pt->SetAxisRange(0, 1.1, "Y");
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_oneMatchAlt_trueVertex_pt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_oneMatchAlt_trueVertex_pt->GetName() + ".pdf");
  delete h_eff_oneMatchAlt_trueVertex_pt;
  delete h_all_oneMatchAlt_trueVertex_pt;
  delete h_correct_oneMatchAlt_trueVertex_pt;
  
  if(detailedPlots){
    std::vector<double> eff_dxy_cuts;
    for(uint i=0;i<dxy_cuts.size();i++){
      eff_dxy_cuts.push_back(correct_vert_dxy_cut[i] / true_vertices);
    }
    TGraph* gr_eff_dxy_cuts = new TGraph(int(dxy_cuts.size()),dxy_cuts.data(),eff_dxy_cuts.data());
    gr_eff_dxy_cuts->SetTitle("Eff vs Dxy Cut; Dxy Cut (cm); Efficiency");
    gr_eff_dxy_cuts->SetName("gr_eff_dxy_cuts");
    gr_eff_dxy_cuts->Draw("AC*");
    c.SaveAs(DIR + "/"+ gr_eff_dxy_cuts->GetName() + ".pdf");
    delete gr_eff_dxy_cuts;
    
    std::vector<double> eff_dz_cuts;
    for(uint i=0;i<dz_cuts.size();i++){
      eff_dz_cuts.push_back(correct_vert_dz_cut[i] / true_vertices);
    }
    TGraph* gr_eff_dz_cuts = new TGraph(int(dz_cuts.size()),dz_cuts.data(),eff_dz_cuts.data());
    gr_eff_dz_cuts->SetTitle("Eff vs Dz Cut; Dz Cut (cm); Efficiency");
    gr_eff_dz_cuts->SetName("gr_eff_dz_cuts");
    gr_eff_dz_cuts->Draw("AC*");
    c.SaveAs(DIR + "/"+ gr_eff_dz_cuts->GetName() + ".pdf");
    delete gr_eff_dz_cuts;
    
    std::vector<double> eff_cos_T_cuts;
    for(uint i=0;i<cos_T_cuts.size();i++){
      eff_cos_T_cuts.push_back(correct_vert_cos_T_cut[i] / true_vertices);
    }
    TGraph* gr_eff_cos_T_cuts = new TGraph(int(cos_T_cuts.size()),cos_T_cuts.data(),eff_cos_T_cuts.data());
    gr_eff_cos_T_cuts->SetTitle("Eff vs Cos(alpha) Cut; Cos(alpha) Cut; Efficiency");
    gr_eff_cos_T_cuts->SetName("gr_eff_cos_T_cuts");
    gr_eff_cos_T_cuts->Draw("AC*");
    c.SaveAs(DIR + "/"+ gr_eff_cos_T_cuts->GetName() + ".pdf");
    delete gr_eff_cos_T_cuts;
    
    std::vector<double> eff_d_T_cuts;
    for(uint i=0;i<d_T_cuts.size();i++){
      eff_d_T_cuts.push_back(correct_vert_d_T_cut[i] / true_vertices);
    }
    TGraph* gr_eff_d_T_cuts = new TGraph(int(d_T_cuts.size()),d_T_cuts.data(),eff_d_T_cuts.data());
    gr_eff_d_T_cuts->SetTitle("Eff vs d_T Cut; d_T Cut (cm); Efficiency");
    gr_eff_d_T_cuts->SetName("gr_eff_d_T_cuts");
    gr_eff_d_T_cuts->Draw("AC*");
    c.SaveAs(DIR + "/"+ gr_eff_d_T_cuts->GetName() + ".pdf");
    delete gr_eff_d_T_cuts;
    
    std::vector<double> eff_R_T_cuts;
    for(uint i=0;i<R_T_cuts.size();i++){
      eff_R_T_cuts.push_back(correct_vert_R_T_cut[i] / true_vertices);
    }
    TGraph* gr_eff_R_T_cuts = new TGraph(int(R_T_cuts.size()),R_T_cuts.data(),eff_R_T_cuts.data());
    gr_eff_R_T_cuts->SetTitle("Eff vs R_T Cut; R_T Cut (cm); Efficiency");
    gr_eff_R_T_cuts->SetName("gr_eff_R_T_cuts");
    gr_eff_R_T_cuts->Draw("AC*");
    c.SaveAs(DIR + "/"+ gr_eff_R_T_cuts->GetName() + ".pdf");
    delete gr_eff_R_T_cuts;
    
    std::vector<double> eff_chi2rz_cuts;
    for(uint i=0;i<chi2rz_cuts.size();i++){
      eff_chi2rz_cuts.push_back(correct_vert_chi2rz_cut[i] / true_vertices);
    }
    TGraph* gr_eff_chi2rz_cuts = new TGraph(int(chi2rz_cuts.size()),chi2rz_cuts.data(),eff_chi2rz_cuts.data());
    gr_eff_chi2rz_cuts->SetTitle("Eff vs chi2rz Cut; chi2rz Cut; Efficiency");
    gr_eff_chi2rz_cuts->SetName("gr_eff_chi2rz_cuts");
    gr_eff_chi2rz_cuts->Draw("AC*");
    c.SaveAs(DIR + "/"+ gr_eff_chi2rz_cuts->GetName() + ".pdf");
    delete gr_eff_chi2rz_cuts;
    
    std::vector<double> eff_minD0_cuts;
    for(uint i=0;i<minD0_cuts.size();i++){
      eff_minD0_cuts.push_back(correct_vert_minD0_cut[i] / true_vertices);
    }
    TGraph* gr_eff_minD0_cuts = new TGraph(int(minD0_cuts.size()),minD0_cuts.data(),eff_minD0_cuts.data());
    gr_eff_minD0_cuts->SetTitle("Eff vs minD0 Cut; minD0 Cut (cm); Efficiency");
    gr_eff_minD0_cuts->SetName("gr_eff_minD0_cuts");
    gr_eff_minD0_cuts->Draw("AC*");
    c.SaveAs(DIR + "/"+ gr_eff_minD0_cuts->GetName() + ".pdf");
    delete gr_eff_minD0_cuts;
    
    std::vector<double> eff_stubSum_cuts;
    for(uint i=0;i<stubSum_cuts.size();i++){
      eff_stubSum_cuts.push_back(correct_vert_stubSum_cut[i] / true_vertices);
    }
    TGraph* gr_eff_stubSum_cuts = new TGraph(int(stubSum_cuts.size()),stubSum_cuts.data(),eff_stubSum_cuts.data());
    gr_eff_stubSum_cuts->SetTitle("Eff vs stubSum Cut; stubSum Cut; Efficiency");
    gr_eff_stubSum_cuts->SetName("gr_eff_stubSum_cuts");
    gr_eff_stubSum_cuts->Draw("AC*");
    c.SaveAs(DIR + "/"+ gr_eff_stubSum_cuts->GetName() + ".pdf");
    delete gr_eff_stubSum_cuts;
    
    std::vector<double> false_dxy_cuts;
    for(uint i=0;i<dxy_cuts.size();i++){
      if(all_vert_dxy_cut[i]!=0){
	false_dxy_cuts.push_back(false_vert_dxy_cut[i] / all_vert_dxy_cut[i]);
      }
      else{
	false_dxy_cuts.push_back(0);
      }
    }
    TGraph* gr_false_dxy_cuts = new TGraph(int(dxy_cuts.size()),dxy_cuts.data(),false_dxy_cuts.data());
    gr_false_dxy_cuts->SetTitle("False Rate vs Dxy Cut; Dxy Cut (cm); False Rate");
    gr_false_dxy_cuts->SetName("gr_false_dxy_cuts");
    gr_false_dxy_cuts->Draw("AC*");
    c.SaveAs(DIR + "/"+ gr_false_dxy_cuts->GetName() + ".pdf");
    delete gr_false_dxy_cuts;
    
    std::vector<double> false_dz_cuts;
    for(uint i=0;i<dz_cuts.size();i++){
      if(all_vert_dz_cut[i]!=0){
	false_dz_cuts.push_back(false_vert_dz_cut[i] / all_vert_dz_cut[i]);
      }
      else{
	false_dz_cuts.push_back(0);
      }
    }
    TGraph* gr_false_dz_cuts = new TGraph(int(dz_cuts.size()),dz_cuts.data(),false_dz_cuts.data());
    gr_false_dz_cuts->SetTitle("False Rate vs Dz Cut; Dz Cut (cm); False Rate");
    gr_false_dz_cuts->SetName("gr_false_dz_cuts");
    gr_false_dz_cuts->Draw("AC*");
    c.SaveAs(DIR + "/"+ gr_false_dz_cuts->GetName() + ".pdf");
    delete gr_false_dz_cuts;
    
    TGraph* gr_eff_vs_false_dz_cuts = new TGraph(int(dz_cuts.size()),eff_dz_cuts.data(),false_dz_cuts.data());
    gr_eff_vs_false_dz_cuts->SetTitle("False Rate vs Efficiency for Dz Cut; Efficiency; False Rate");
    gr_eff_vs_false_dz_cuts->SetName("gr_eff_vs_false_dz_cuts");
    gr_eff_vs_false_dz_cuts->Draw("AC*");
    c.SaveAs(DIR + "/"+ gr_eff_vs_false_dz_cuts->GetName() + ".pdf");
    delete gr_eff_vs_false_dz_cuts;
   
    std::vector<double> false_cos_T_cuts;
    for(uint i=0;i<cos_T_cuts.size();i++){
      if(all_vert_cos_T_cut[i]!=0){
	false_cos_T_cuts.push_back(false_vert_cos_T_cut[i] / all_vert_cos_T_cut[i]);
      }
      else{
	false_cos_T_cuts.push_back(0);
      }
    }
    TGraph* gr_false_cos_T_cuts = new TGraph(int(cos_T_cuts.size()),cos_T_cuts.data(),false_cos_T_cuts.data());
    gr_false_cos_T_cuts->SetTitle("False Rate vs Cos(alpha) Cut; Cos(alpha) Cut; False Rate");
    gr_false_cos_T_cuts->SetName("gr_false_cos_T_cuts");
    gr_false_cos_T_cuts->Draw("AC*");
    c.SaveAs(DIR + "/"+ gr_false_cos_T_cuts->GetName() + ".pdf");
    delete gr_false_cos_T_cuts;
    
    TGraph* gr_eff_vs_false_cos_T_cuts = new TGraph(int(cos_T_cuts.size()),eff_cos_T_cuts.data(),false_cos_T_cuts.data());
    gr_eff_vs_false_cos_T_cuts->SetTitle("False Rate vs Efficiency of Cos(alpha) Cut; Efficiency; False Rate");
    gr_eff_vs_false_cos_T_cuts->SetName("gr_eff_vs_false_cos_T_cuts");
    gr_eff_vs_false_cos_T_cuts->Draw("AC*");
    c.SaveAs(DIR + "/"+ gr_eff_vs_false_cos_T_cuts->GetName() + ".pdf");
    delete gr_eff_vs_false_cos_T_cuts;
    
    std::vector<double> false_d_T_cuts;
    for(uint i=0;i<d_T_cuts.size();i++){
      if(all_vert_d_T_cut[i]!=0){
	false_d_T_cuts.push_back(false_vert_d_T_cut[i] / all_vert_d_T_cut[i]);
      }
      else{
	false_d_T_cuts.push_back(0);
      }
    }
    TGraph* gr_false_d_T_cuts = new TGraph(int(d_T_cuts.size()),d_T_cuts.data(),false_d_T_cuts.data());
    gr_false_d_T_cuts->SetTitle("False Rate vs d_T Cut; d_T Cut (cm); False Rate");
    gr_false_d_T_cuts->SetName("gr_false_d_T_cuts");
    gr_false_d_T_cuts->Draw("AC*");
    c.SaveAs(DIR + "/"+ gr_false_d_T_cuts->GetName() + ".pdf");
    delete gr_false_d_T_cuts;
    
    TGraph* gr_eff_vs_false_d_T_cuts = new TGraph(int(d_T_cuts.size()),eff_d_T_cuts.data(),false_d_T_cuts.data());
    gr_eff_vs_false_d_T_cuts->SetTitle("False Rate vs Efficiency of d_T Cut; Efficiency; False Rate");
    gr_eff_vs_false_d_T_cuts->SetName("gr_eff_vs_false_d_T_cuts");
    gr_eff_vs_false_d_T_cuts->Draw("AC*");
    c.SaveAs(DIR + "/"+ gr_eff_vs_false_d_T_cuts->GetName() + ".pdf");
    delete gr_eff_vs_false_d_T_cuts;
    
    std::vector<double> false_R_T_cuts;
    for(uint i=0;i<R_T_cuts.size();i++){
      if(all_vert_R_T_cut[i]!=0){
	false_R_T_cuts.push_back(false_vert_R_T_cut[i] / all_vert_R_T_cut[i]);
      }
      else{
	false_R_T_cuts.push_back(0);
      }
    }
    TGraph* gr_false_R_T_cuts = new TGraph(int(R_T_cuts.size()),R_T_cuts.data(),false_R_T_cuts.data());
    gr_false_R_T_cuts->SetTitle("False Rate vs R_T Cut; R_T Cut (cm); False Rate");
    gr_false_R_T_cuts->SetName("gr_false_R_T_cuts");
    gr_false_R_T_cuts->Draw("AC*");
    c.SaveAs(DIR + "/"+ gr_false_R_T_cuts->GetName() + ".pdf");
    delete gr_false_R_T_cuts;
    
    TGraph* gr_eff_vs_false_R_T_cuts = new TGraph(int(R_T_cuts.size()),eff_R_T_cuts.data(),false_R_T_cuts.data());
    gr_eff_vs_false_R_T_cuts->SetTitle("False Rate vs Efficiency of R_T Cut; Efficiency; False Rate");
    gr_eff_vs_false_R_T_cuts->SetName("gr_eff_vs_false_R_T_cuts");
    gr_eff_vs_false_R_T_cuts->Draw("AC*");
    c.SaveAs(DIR + "/"+ gr_eff_vs_false_R_T_cuts->GetName() + ".pdf");
    delete gr_eff_vs_false_R_T_cuts;
    
    std::vector<double> false_chi2rz_cuts;
    for(uint i=0;i<chi2rz_cuts.size();i++){
      if(all_vert_chi2rz_cut[i]!=0){
	false_chi2rz_cuts.push_back(false_vert_chi2rz_cut[i] / all_vert_chi2rz_cut[i]);
      }
      else{
	false_chi2rz_cuts.push_back(0);
      }
    }
    TGraph* gr_false_chi2rz_cuts = new TGraph(int(chi2rz_cuts.size()),chi2rz_cuts.data(),false_chi2rz_cuts.data());
    gr_false_chi2rz_cuts->SetTitle("False Rate vs chi2rz Cut; chi2rz Cut; False Rate");
    gr_false_chi2rz_cuts->SetName("gr_false_chi2rz_cuts");
    gr_false_chi2rz_cuts->Draw("AC*");
    c.SaveAs(DIR + "/"+ gr_false_chi2rz_cuts->GetName() + ".pdf");
    delete gr_false_chi2rz_cuts;
    
    TGraph* gr_eff_vs_false_chi2rz_cuts = new TGraph(int(chi2rz_cuts.size()),eff_chi2rz_cuts.data(),false_chi2rz_cuts.data());
    gr_eff_vs_false_chi2rz_cuts->SetTitle("False Rate vs Efficiency of chi2rz Cut; Efficiency; False Rate");
    gr_eff_vs_false_chi2rz_cuts->SetName("gr_eff_vs_false_chi2rz_cuts");
    gr_eff_vs_false_chi2rz_cuts->Draw("AC*");
    c.SaveAs(DIR + "/"+ gr_eff_vs_false_chi2rz_cuts->GetName() + ".pdf");
    delete gr_eff_vs_false_chi2rz_cuts;
    
    std::vector<double> false_minD0_cuts;
    for(uint i=0;i<minD0_cuts.size();i++){
      if(all_vert_minD0_cut[i]!=0){
	false_minD0_cuts.push_back(false_vert_minD0_cut[i] / all_vert_minD0_cut[i]);
      }
      else{
	false_minD0_cuts.push_back(0);
      }
    }
    TGraph* gr_false_minD0_cuts = new TGraph(int(minD0_cuts.size()),minD0_cuts.data(),false_minD0_cuts.data());
    gr_false_minD0_cuts->SetTitle("False Rate vs minD0 Cut; minD0 Cut (cm); False Rate");
    gr_false_minD0_cuts->SetName("gr_false_minD0_cuts");
    gr_false_minD0_cuts->Draw("AC*");
    c.SaveAs(DIR + "/"+ gr_false_minD0_cuts->GetName() + ".pdf");
    delete gr_false_minD0_cuts;
    
    TGraph* gr_eff_vs_false_minD0_cuts = new TGraph(int(minD0_cuts.size()),eff_minD0_cuts.data(),false_minD0_cuts.data());
    gr_eff_vs_false_minD0_cuts->SetTitle("False Rate vs Efficiency of minD0 Cut; Efficiency; False Rate");
    gr_eff_vs_false_minD0_cuts->SetName("gr_eff_vs_false_minD0_cuts");
    gr_eff_vs_false_minD0_cuts->Draw("AC*");
    c.SaveAs(DIR + "/"+ gr_eff_vs_false_minD0_cuts->GetName() + ".pdf");
    delete gr_eff_vs_false_minD0_cuts;
    
    std::vector<double> false_stubSum_cuts;
    for(uint i=0;i<stubSum_cuts.size();i++){
      if(all_vert_stubSum_cut[i]!=0){
	false_stubSum_cuts.push_back(false_vert_stubSum_cut[i] / all_vert_stubSum_cut[i]);
      }
      else{
	false_stubSum_cuts.push_back(0);
      }
    }
    TGraph* gr_false_stubSum_cuts = new TGraph(int(stubSum_cuts.size()),stubSum_cuts.data(),false_stubSum_cuts.data());
    gr_false_stubSum_cuts->SetTitle("False Rate vs stubSum Cut; stubSum Cut; False Rate");
    gr_false_stubSum_cuts->SetName("gr_false_stubSum_cuts");
    gr_false_stubSum_cuts->Draw("AC*");
    c.SaveAs(DIR + "/"+ gr_false_stubSum_cuts->GetName() + ".pdf");
    delete gr_false_stubSum_cuts;
    
    TGraph* gr_eff_vs_false_stubSum_cuts = new TGraph(int(stubSum_cuts.size()),eff_stubSum_cuts.data(),false_stubSum_cuts.data());
    gr_eff_vs_false_stubSum_cuts->SetTitle("False Rate vs Efficiency of stubSum Cut; Efficiency; False Rate");
    gr_eff_vs_false_stubSum_cuts->SetName("gr_eff_vs_false_stubSum_cuts");
    gr_eff_vs_false_stubSum_cuts->Draw("AC*");
    c.SaveAs(DIR + "/"+ gr_eff_vs_false_stubSum_cuts->GetName() + ".pdf");
    delete gr_eff_vs_false_stubSum_cuts;
    
    TH2F *array_cut_num = new TH2F("array_cut_num","array_cut_num",20,0,1,20,0,1);
    TH2F *array_cut_i = new TH2F("array_cut_i","array_cut_i",20,0,1,20,0,1);
    TH2F *array_cut_j = new TH2F("array_cut_j","array_cut_j",20,0,1,20,0,1);
    for(uint k_a=0;k_a<dz_cuts.size();k_a++){
      for(uint k_b=0;k_b<cos_T_cuts.size();k_b++){
	for(uint k_c=0;k_c<d_T_cuts.size();k_c++){
	  for(uint k_d=0;k_d<R_T_cuts.size();k_d++){
	    for(uint k_e=0;k_e<chi2rz_cuts.size();k_e++){
	      for(uint k_f=0;k_f<minD0_cuts.size();k_f++){
		for(uint k_g=0;k_g<stubSum_cuts.size();k_g++){
		  double eff = correct_vert_array_cut[k_a][k_b][k_c][k_d][k_e][k_f][k_g] / true_vertices;
		  double fake = (all_vert_array_cut[k_a][k_b][k_c][k_d][k_e][k_f][k_g]!=0) ? (false_vert_array_cut[k_a][k_b][k_c][k_d][k_e][k_f][k_g]/all_vert_array_cut[k_a][k_b][k_c][k_d][k_e][k_f][k_g]) : 0;
		  array_cut_num->Fill(fake,eff,1);
		  if(array_cut_num->GetBinContent(int(fake/0.05)+1,int(eff/0.05)+1)==1){
		    int indices_i = k_a * 1000 + k_b * 100 + k_c * 10 + k_d; 
		    int indices_j = k_e * 100 + k_f * 10 + k_g;
		    array_cut_i->Fill(fake,eff,indices_i);
		    array_cut_j->Fill(fake,eff,indices_j);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    array_cut_i->SetMarkerSize(1.0);
    array_cut_j->SetMarkerSize(1.0);
    array_cut_i->SetBarOffset(0.2);
    array_cut_j->SetBarOffset(-0.2);
    array_cut_num->GetXaxis()->SetTitle("Fake Rate");
    array_cut_num->GetYaxis()->SetTitle("Efficiency");
    array_cut_num->Draw("COLZ");
    array_cut_i->Draw("TEXT SAME");
    array_cut_j->Draw("TEXT SAME");
    c.SaveAs(DIR + "/"+ "effVsFakeArray.pdf");
  }
  
  h_all_trackVertex_pt->Sumw2();
  h_false_trackVertex_pt->Sumw2();
  TH1F* h_fake_trackVertex_pt = (TH1F*)h_false_trackVertex_pt->Clone();
  h_fake_trackVertex_pt->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_fake_trackVertex_pt);
  h_fake_trackVertex_pt->SetName("fake_trackVertex_pt");
  h_fake_trackVertex_pt->GetXaxis()->SetTitle("Track p_{T} (GeV)");
  h_fake_trackVertex_pt->GetYaxis()->SetTitle("Fake Rate");
  h_fake_trackVertex_pt->Divide(h_false_trackVertex_pt,h_all_trackVertex_pt, 1.0, 1.0, "B");
  h_fake_trackVertex_pt->SetStats(0);
  h_fake_trackVertex_pt->Draw();
  h_fake_trackVertex_pt->SetAxisRange(0, 1.1, "Y");
  mySmallText(0.4, 0.82, 1, ctxt);
  h_fake_trackVertex_pt->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_fake_trackVertex_pt->GetName() + ".pdf");
  delete h_fake_trackVertex_pt;
  delete h_all_trackVertex_pt;
  delete h_false_trackVertex_pt;

  removeFlows(h_all_trueVertex_eta);
  h_all_trueVertex_eta->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_all_trueVertex_eta->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_all_trueVertex_eta->GetName() + ".pdf");

  removeFlows(h_correct_trueVertex_eta);
  h_correct_trueVertex_eta->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_correct_trueVertex_eta->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_correct_trueVertex_eta->GetName() + ".pdf");

  h_all_trueVertex_eta->Sumw2();
  h_correct_trueVertex_eta->Sumw2();
  TH1F* h_eff_trueVertex_eta = (TH1F*)h_correct_trueVertex_eta->Clone();
  h_eff_trueVertex_eta->GetYaxis()->SetNoExponent(kTRUE);
  h_eff_trueVertex_eta->SetName("eff_trueVertex_eta");
  h_eff_trueVertex_eta->GetXaxis()->SetTitle("Tracking Particle #eta");
  h_eff_trueVertex_eta->GetYaxis()->SetTitle("Efficiency");
  h_eff_trueVertex_eta->Divide(h_correct_trueVertex_eta,h_all_trueVertex_eta, 1.0, 1.0, "B");
  h_eff_trueVertex_eta->SetStats(0);
  h_eff_trueVertex_eta->Draw();
  h_eff_trueVertex_eta->SetAxisRange(0, 1.1, "Y");
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_trueVertex_eta->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_trueVertex_eta->GetName() + ".pdf");
  delete h_eff_trueVertex_eta;
  delete h_all_trueVertex_eta;
  delete h_correct_trueVertex_eta;

  h_all_oneMatch_trueVertex_eta->Sumw2();
  h_correct_oneMatch_trueVertex_eta->Sumw2();
  TH1F* h_eff_oneMatch_trueVertex_eta = (TH1F*)h_correct_oneMatch_trueVertex_eta->Clone();
  h_eff_oneMatch_trueVertex_eta->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_oneMatch_trueVertex_eta);
  h_eff_oneMatch_trueVertex_eta->SetName("eff_oneMatch_trueVertex_eta");
  h_eff_oneMatch_trueVertex_eta->GetXaxis()->SetTitle("Tracking Particle #eta");
  h_eff_oneMatch_trueVertex_eta->GetYaxis()->SetTitle("Efficiency");
  h_eff_oneMatch_trueVertex_eta->Divide(h_correct_oneMatch_trueVertex_eta,h_all_oneMatch_trueVertex_eta, 1.0, 1.0, "B");
  h_eff_oneMatch_trueVertex_eta->SetStats(0);
  h_eff_oneMatch_trueVertex_eta->Draw();
  h_eff_oneMatch_trueVertex_eta->SetAxisRange(0, 1.1, "Y");
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_oneMatch_trueVertex_eta->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_oneMatch_trueVertex_eta->GetName() + ".pdf");
  delete h_eff_oneMatch_trueVertex_eta;
  delete h_all_oneMatch_trueVertex_eta;
  delete h_correct_oneMatch_trueVertex_eta;

  h_all_oneMatchAlt_trueVertex_eta->Sumw2();
  h_correct_oneMatchAlt_trueVertex_eta->Sumw2();
  TH1F* h_eff_oneMatchAlt_trueVertex_eta = (TH1F*)h_correct_oneMatchAlt_trueVertex_eta->Clone();
  h_eff_oneMatchAlt_trueVertex_eta->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_oneMatchAlt_trueVertex_eta);
  h_eff_oneMatchAlt_trueVertex_eta->SetName("eff_oneMatchAlt_trueVertex_eta");
  h_eff_oneMatchAlt_trueVertex_eta->GetXaxis()->SetTitle("Tracking Particle #eta");
  h_eff_oneMatchAlt_trueVertex_eta->GetYaxis()->SetTitle("Efficiency");
  h_eff_oneMatchAlt_trueVertex_eta->Divide(h_correct_oneMatchAlt_trueVertex_eta,h_all_oneMatchAlt_trueVertex_eta, 1.0, 1.0, "B");
  h_eff_oneMatchAlt_trueVertex_eta->Draw();
  h_eff_oneMatchAlt_trueVertex_eta->SetAxisRange(0, 1.1, "Y");
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_oneMatchAlt_trueVertex_eta->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_oneMatchAlt_trueVertex_eta->GetName() + ".pdf");
  delete h_eff_oneMatchAlt_trueVertex_eta;
  delete h_all_oneMatchAlt_trueVertex_eta;
  delete h_correct_oneMatchAlt_trueVertex_eta;

  h_all_trackVertex_eta->Sumw2();
  h_false_trackVertex_eta->Sumw2();
  TH1F* h_fake_trackVertex_eta = (TH1F*)h_false_trackVertex_eta->Clone();
  h_fake_trackVertex_eta->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_fake_trackVertex_eta);
  h_fake_trackVertex_eta->SetName("fake_trackVertex_eta");
  h_fake_trackVertex_eta->GetXaxis()->SetTitle("Track #eta");
  h_fake_trackVertex_eta->GetYaxis()->SetTitle("Fake Rate");
  h_fake_trackVertex_eta->Divide(h_false_trackVertex_eta,h_all_trackVertex_eta, 1.0, 1.0, "B");
  h_fake_trackVertex_eta->SetStats(0);
  h_fake_trackVertex_eta->Draw();
  h_fake_trackVertex_eta->SetAxisRange(0, 1.1, "Y");
  mySmallText(0.4, 0.82, 1, ctxt);
  h_fake_trackVertex_eta->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_fake_trackVertex_eta->GetName() + ".pdf");
  delete h_fake_trackVertex_eta;
  delete h_all_trackVertex_eta;
  delete h_false_trackVertex_eta;

  removeFlows(h_all_trueVertex_dxy);
  h_all_trueVertex_dxy->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_all_trueVertex_dxy->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_all_trueVertex_dxy->GetName() + ".pdf");

  removeFlows(h_correct_trueVertex_dxy);
  h_correct_trueVertex_dxy->Draw();
  mySmallText(0.4, 0.82, 1, ctxt);
  h_correct_trueVertex_dxy->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_correct_trueVertex_dxy->GetName() + ".pdf");
   
  h_all_trueVertex_dxy->Sumw2();
  h_correct_trueVertex_dxy->Sumw2();
  TH1F* h_eff_trueVertex_dxy = (TH1F*)h_correct_trueVertex_dxy->Clone();
  h_eff_trueVertex_dxy->GetYaxis()->SetNoExponent(kTRUE);
  h_eff_trueVertex_dxy->SetName("eff_trueVertex_dxy");
  h_eff_trueVertex_dxy->GetXaxis()->SetTitle("Vertex d_{xy} (cm)");
  h_eff_trueVertex_dxy->GetYaxis()->SetTitle("Efficiency");
  h_eff_trueVertex_dxy->Divide(h_correct_trueVertex_dxy,h_all_trueVertex_dxy, 1.0, 1.0, "B");
  h_eff_trueVertex_dxy->SetStats(0);
  h_eff_trueVertex_dxy->Draw();
  h_eff_trueVertex_dxy->SetAxisRange(0, 1.1, "Y");
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_trueVertex_dxy->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_trueVertex_dxy->GetName() + ".pdf");
  delete h_eff_trueVertex_dxy;
  delete h_all_trueVertex_dxy;
  delete h_correct_trueVertex_dxy;

  h_all_oneMatch_trueVertex_dxy->Sumw2();
  h_correct_oneMatch_trueVertex_dxy->Sumw2();
  TH1F* h_eff_oneMatch_trueVertex_dxy = (TH1F*)h_correct_oneMatch_trueVertex_dxy->Clone();
  h_eff_oneMatch_trueVertex_dxy->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_oneMatch_trueVertex_dxy);
  h_eff_oneMatch_trueVertex_dxy->SetName("eff_oneMatch_trueVertex_dxy");
  h_eff_oneMatch_trueVertex_dxy->GetXaxis()->SetTitle("Vertex d_{xy} (cm)");
  h_eff_oneMatch_trueVertex_dxy->GetYaxis()->SetTitle("Efficiency");
  h_eff_oneMatch_trueVertex_dxy->Divide(h_correct_oneMatch_trueVertex_dxy,h_all_oneMatch_trueVertex_dxy, 1.0, 1.0, "B");
  h_eff_oneMatch_trueVertex_dxy->SetStats(0);
  h_eff_oneMatch_trueVertex_dxy->Draw();
  h_eff_oneMatch_trueVertex_dxy->SetAxisRange(0, 1.1, "Y");
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_oneMatch_trueVertex_dxy->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_oneMatch_trueVertex_dxy->GetName() + ".pdf");
  delete h_eff_oneMatch_trueVertex_dxy;
  delete h_all_oneMatch_trueVertex_dxy;
  delete h_correct_oneMatch_trueVertex_dxy;

  h_all_oneMatchAlt_trueVertex_dxy->Sumw2();
  h_correct_oneMatchAlt_trueVertex_dxy->Sumw2();
  TH1F* h_eff_oneMatchAlt_trueVertex_dxy = (TH1F*)h_correct_oneMatchAlt_trueVertex_dxy->Clone();
  h_eff_oneMatchAlt_trueVertex_dxy->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_eff_oneMatchAlt_trueVertex_dxy);
  h_eff_oneMatchAlt_trueVertex_dxy->SetName("eff_oneMatchAlt_trueVertex_dxy");
  h_eff_oneMatchAlt_trueVertex_dxy->GetXaxis()->SetTitle("Vertex d_{xy} (cm)");
  h_eff_oneMatchAlt_trueVertex_dxy->GetYaxis()->SetTitle("Efficiency");
  h_eff_oneMatchAlt_trueVertex_dxy->Divide(h_correct_oneMatchAlt_trueVertex_dxy,h_all_oneMatchAlt_trueVertex_dxy, 1.0, 1.0, "B");
  h_eff_oneMatchAlt_trueVertex_dxy->Draw();
  h_eff_oneMatchAlt_trueVertex_dxy->SetAxisRange(0, 1.1, "Y");
  mySmallText(0.4, 0.82, 1, ctxt);
  h_eff_oneMatchAlt_trueVertex_dxy->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_eff_oneMatchAlt_trueVertex_dxy->GetName() + ".pdf");
  delete h_eff_oneMatchAlt_trueVertex_dxy;
  delete h_all_oneMatchAlt_trueVertex_dxy;
  delete h_correct_oneMatchAlt_trueVertex_dxy;
   
  h_all_trackVertex_dxy->Sumw2();
  h_false_trackVertex_dxy->Sumw2();
  TH1F* h_fake_trackVertex_dxy = (TH1F*)h_false_trackVertex_dxy->Clone();
  h_fake_trackVertex_dxy->GetYaxis()->SetNoExponent(kTRUE);
  removeFlows(h_fake_trackVertex_dxy);
  h_fake_trackVertex_dxy->SetName("fake_trackVertex_dxy");
  h_fake_trackVertex_dxy->GetXaxis()->SetTitle("Vertex d_{xy} (cm)");
  h_fake_trackVertex_dxy->GetYaxis()->SetTitle("Fake Rate");
  h_fake_trackVertex_dxy->Divide(h_false_trackVertex_dxy,h_all_trackVertex_dxy, 1.0, 1.0, "B");
  h_fake_trackVertex_dxy->SetStats(0);
  h_fake_trackVertex_dxy->Draw();
  h_fake_trackVertex_dxy->SetAxisRange(0, 1.1, "Y");
  mySmallText(0.4, 0.82, 1, ctxt);
  h_fake_trackVertex_dxy->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_fake_trackVertex_dxy->GetName() + ".pdf");
  delete h_fake_trackVertex_dxy;
  delete h_all_trackVertex_dxy;
  delete h_false_trackVertex_dxy;
   
  std::stringstream txt;
  float rate;
  if(inputFile.Contains("NeutrinoGun")){
    rate = 40000.0;
    txt << " kHz ";
  }
  else{
    rate = 100.0;
    txt << " % ";
  }

  for(int k=0;k<5;k++){
    for (int l=0; l < 5; l++){
      h_Count_trk_pt_d0->SetBinContent(l+1,k+1,(h_Count_trk_pt_d0->GetBinContent(l+1,k+1) * rate/nevt));
      h_Count_trk_pt_d0_dv->SetBinContent(l+1,k+1,(h_Count_trk_pt_d0_dv->GetBinContent(l+1,k+1) * rate/nevt));
    }
  }
   
  h_Count_trk_pt_d0->Draw("colz");
  h_Count_trk_pt_d0->SetMarkerSize(2);
  h_Count_trk_pt_d0->Draw("textsame");
  h_Count_trk_pt_d0->GetXaxis()->SetNdivisions(5);
  h_Count_trk_pt_d0->GetYaxis()->SetNdivisions(5);
  h_Count_trk_pt_d0->GetXaxis()->CenterLabels();
  h_Count_trk_pt_d0->GetYaxis()->CenterLabels();
  h_Count_trk_pt_d0->GetYaxis()->SetTitle("Transverse Momentum p_T (GeV)");
  h_Count_trk_pt_d0->GetXaxis()->SetTitle("Transverse Impact Parameter d_0 (cm)");
  h_Count_trk_pt_d0->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_Count_trk_pt_d0->GetName() + ".pdf");
  delete h_Count_trk_pt_d0;

  h_Count_trk_pt_d0_dv->Draw("colz");
  h_Count_trk_pt_d0_dv->SetMarkerSize(2);
  h_Count_trk_pt_d0_dv->Draw("textsame");
  h_Count_trk_pt_d0_dv->GetXaxis()->SetNdivisions(5);
  h_Count_trk_pt_d0_dv->GetYaxis()->SetNdivisions(5);
  h_Count_trk_pt_d0_dv->GetXaxis()->CenterLabels();
  h_Count_trk_pt_d0_dv->GetYaxis()->CenterLabels();
  h_Count_trk_pt_d0_dv->GetYaxis()->SetTitle("Transverse Momentum p_T (GeV)");
  h_Count_trk_pt_d0_dv->GetXaxis()->SetTitle("(DV selection) Transverse Impact Parameter d_0 (cm)");
  h_Count_trk_pt_d0_dv->Write("", TObject::kOverwrite);
  c.SaveAs(DIR + "/"+ h_Count_trk_pt_d0_dv->GetName() + ".pdf");
  delete h_Count_trk_pt_d0_dv;

  //Geometric plot of circle projections and vertex locations
  Double_t x_min=0;
  Double_t x_max=0;
  Double_t y_min=0;
  Double_t y_max=0;
  Double_t x_values[6] = {geomTrackVertex.a.x0-geomTrackVertex.a.rho,geomTrackVertex.a.x0+geomTrackVertex.a.rho,geomTrackVertex.b.x0-geomTrackVertex.b.rho,geomTrackVertex.b.x0+geomTrackVertex.b.rho,geomTrackVertex.x_dv,geomTrueVertex.x_dv};
  Double_t y_values[6] = {geomTrackVertex.a.y0-geomTrackVertex.a.rho,geomTrackVertex.a.y0+geomTrackVertex.a.rho,geomTrackVertex.b.y0-geomTrackVertex.b.rho,geomTrackVertex.b.y0+geomTrackVertex.b.rho,geomTrackVertex.y_dv,geomTrueVertex.y_dv};
  for(uint i=0;i<6;i++){
    if(x_values[i]<x_min) x_min = x_values[i];
    if(x_values[i]>x_max) x_max = x_values[i];
    if(y_values[i]<y_min) y_min = y_values[i];
    if(y_values[i]>y_max) y_max = y_values[i];
  }
  x_min*=1.1;
  x_max*=1.1;
  y_min*=1.1;
  y_max*=1.1;
  c.DrawFrame(x_min,y_min,x_max,y_max);
  float trk1_POCA_x = geomTrackVertex.a.d0*sin(geomTrackVertex.a.phi);
  float trk1_POCA_y = -1*geomTrackVertex.a.d0*cos(geomTrackVertex.a.phi);
  float trk2_POCA_x = geomTrackVertex.b.d0*sin(geomTrackVertex.b.phi);
  float trk2_POCA_y = -1*geomTrackVertex.b.d0*cos(geomTrackVertex.b.phi);
  TEllipse *circleTrk1 = new TEllipse(geomTrackVertex.a.x0,geomTrackVertex.a.y0,geomTrackVertex.a.rho,geomTrackVertex.a.rho);
  circleTrk1->SetLineColor(kGreen);
  circleTrk1->SetFillStyle(0);
  TEllipse *circleTrk2 = new TEllipse(geomTrackVertex.b.x0,geomTrackVertex.b.y0,geomTrackVertex.b.rho,geomTrackVertex.b.rho);
  circleTrk2->SetLineColor(kBlack);
  circleTrk2->SetFillStyle(0);
  auto trackTraj1 = new TF1("trackTraj1","((sin([0])/cos([0]))*(x-[1]))+[2]",x_min,x_max);
  trackTraj1->SetParameters(geomTrackVertex.a.phi,trk1_POCA_x,trk1_POCA_y);
  trackTraj1->SetLineColor(kGreen);
  auto trackTraj2 = new TF1("trackTraj2","((sin([0])/cos([0]))*(x-[1]))+[2]",x_min,x_max);
  trackTraj2->SetParameters(geomTrackVertex.b.phi,trk2_POCA_x,trk2_POCA_y);
  trackTraj2->SetLineColor(kBlack);
  TMarker m1(geomTrackVertex.x_dv,geomTrackVertex.y_dv,8);
  TMarker m2(geomTrueVertex.x_dv,geomTrueVertex.y_dv,8);
  TMarker m3(trk1_POCA_x,trk1_POCA_y,5);
  TMarker m4(trk2_POCA_x,trk2_POCA_y,5);
  //std::cout<<"trk1 POCA: "<<trk1_POCA_x<<" "<<trk1_POCA_y<<" trk2 POCA: "<<trk2_POCA_x<<" "<<trk2_POCA_y<<std::endl;
  m1.SetMarkerColor(kRed);
  m2.SetMarkerColor(kBlue);
  m3.SetMarkerColor(kRed);
  m4.SetMarkerColor(kBlue);
  circleTrk1->Draw("SAME");
  circleTrk2->Draw("SAME");
  trackTraj1->Draw("SAME");
  trackTraj2->Draw("SAME");
  m1.Draw("SAME");
  m2.Draw("SAME");
  m3.Draw("SAME");
  m4.Draw("SAME");
  c.SaveAs(DIR + "/h_circleFitGeom.pdf");
  c.SaveAs(DIR + "/h_circleFitGeom.pdf");
  x_min = geomTrackVertex.x_dv;
  x_max = geomTrueVertex.x_dv;
  y_min = geomTrackVertex.y_dv;
  y_max = geomTrueVertex.y_dv;
  if(geomTrueVertex.x_dv<x_min){
    x_max = x_min;
    x_min = geomTrueVertex.x_dv; 
  }
  if(geomTrueVertex.y_dv<y_min){
    y_max = y_min;
    y_min = geomTrueVertex.y_dv; 
  }
  //std::cout<<"geom track vertex: "<<geomTrackVertex.x_dv<<" "<<geomTrackVertex.y_dv<<" geom true vertex: "<<geomTrueVertex.x_dv<<" "<<geomTrueVertex.y_dv<<std::endl;
  //std::cout<<"x_min: "<<x_min<<" x_max: "<<x_max<<" y_min: "<<y_min<<" y_max: "<<y_max<<std::endl;
  x_min-=2;
  x_max+=2;
  y_min-=2;
  y_max+=2;
  c.DrawFrame(x_min,y_min,x_max,y_max);
  circleTrk1->Draw("SAME");
  circleTrk2->Draw("SAME");
  //trackTraj1->Draw("SAME");
  //trackTraj2->Draw("SAME");
  m1.Draw("SAME");
  m2.Draw("SAME");
  m3.Draw("SAME");
  m4.Draw("SAME");
  c.SaveAs(DIR + "/h_circleFitGeomZoom.pdf");
  c.SaveAs(DIR + "/h_circleFitGeomZoom.pdf");
  x_min=geomTrackVertex.x_dv;
  x_max=geomTrackVertex.x_dv;
  y_min=geomTrackVertex.y_dv;
  y_max=geomTrackVertex.y_dv;
  Double_t x_values_POCA[4] = {geomTrackVertex.x_dv,geomTrueVertex.x_dv,trk1_POCA_x,trk2_POCA_x};
  Double_t y_values_POCA[4] = {geomTrackVertex.y_dv,geomTrueVertex.y_dv,trk1_POCA_y,trk2_POCA_y};
  for(uint i=0;i<4;i++){
    if(x_values_POCA[i]<x_min) x_min = x_values_POCA[i];
    if(x_values_POCA[i]>x_max) x_max = x_values_POCA[i];
    if(y_values_POCA[i]<y_min) y_min = y_values_POCA[i];
    if(y_values_POCA[i]>y_max) y_max = y_values_POCA[i];
  }
  x_min-=1;
  x_max+=1;
  y_min-=1;
  y_max+=1;
  c.DrawFrame(x_min,y_min,x_max,y_max);
  circleTrk1->Draw("SAME");
  circleTrk2->Draw("SAME");
  trackTraj1->Draw("SAME");
  trackTraj2->Draw("SAME");
  m1.Draw("SAME");
  m2.Draw("SAME");
  m3.Draw("SAME");
  m4.Draw("SAME");
  c.SaveAs(DIR + "/h_circleFitGeomPOCA.pdf");
  c.SaveAs(DIR + "/h_circleFitGeomPOCA.pdf");
  fout->Close();
}


Double_t dist_TPs(Track_Parameters* a, Track_Parameters* b){  
  float x1 = a->x0; //   Centers of the circles
  float y1 = a->y0; // 
  float x2 = b->x0; // 
  float y2 = b->y0; // 
  float R1 = a->rho;   // Radii of the circles
  float R2 = b->rho;
  float R = dist(x1,y1,x2,y2); // Distance between centers
  if((R>=(R1-R2)) && (R<=(R1+R2))){
    return (0);
  }
  else if(R==0){
    return (-99999.0);
  }
  else{

    return(R-R1-R2);
  }
}

Double_t dist_TPs(Track_Parameters a, Track_Parameters b){  
  float x1 = a.x0; //   Centers of the circles
  float y1 = a.y0; // 
  float x2 = b.x0; // 
  float y2 = b.y0; // 
  float R1 = a.rho;   // Radii of the circles
  float R2 = b.rho;
  float R = dist(x1,y1,x2,y2); // Distance between centers
  if((R>=(R1-R2)) && (R<=(R1+R2))){
    return (0);
  }
  else if(R==0){
    return (-99999.0);
  }
  else{

    return(R-R1-R2);
  }
}

Int_t calcVertex(Track_Parameters a, Track_Parameters b, Double_t &x_vtx, Double_t &y_vtx, Double_t &z_vtx){
  float x1 = a.x0; //   Centers of the circles
  float y1 = a.y0; // 
  float x2 = b.x0; // 
  float y2 = b.y0; // 
  float R1 = a.rho;   // Radii of the circles
  float R2 = b.rho;
  float R = dist(x1,y1,x2,y2); // Distance between centers
  if(R==0) return -1;
  float co1 = (pow(R1,2)-pow(R2,2))/(2*pow(R,2));
  float radicand = (2/pow(R,2))*(pow(R1,2)+pow(R2,2))-(pow(pow(R1,2)-pow(R2,2),2)/pow(R,4))-1;
  float co2 = 0;
  if(radicand>0) co2 = 0.5*TMath::Sqrt(radicand);
  float ix1_x = 0.5*(x1+x2)+co1*(x2-x1)+co2*(y2-y1);
  float ix2_x = 0.5*(x1+x2)+co1*(x2-x1)-co2*(y2-y1);
  float ix1_y = 0.5*(y1+y2)+co1*(y2-y1)+co2*(x1-x2);
  float ix2_y = 0.5*(y1+y2)+co1*(y2-y1)-co2*(x1-x2);
  float ix1_z1 = a.z(ix1_x,ix1_y);
  float ix1_z2 = b.z(ix1_x,ix1_y);
  float ix1_delz = fabs(ix1_z1-ix1_z2); 
  float ix2_z1 = a.z(ix2_x,ix2_y);
  float ix2_z2 = b.z(ix2_x,ix2_y);
  float ix2_delz = fabs(ix2_z1-ix2_z2); 
  //std::cout<<"R: "<<R<<" co1: "<<co1<<" co2: "<<co2<<" ix1 delz: "<<ix1_delz<<" ix2 delz: "<<ix2_delz<<std::endl;
  //std::cout<<"ix1_x: "<<ix1_x<<" ix1_y: "<<ix1_y<<" ix2_x: "<<ix2_x<<" ix2_y: "<<ix2_y<<std::endl;
  //std::cout<<"ix1_z1: "<<ix1_z1<<" ix1_z2: "<<ix1_z2<<" ix2_z1: "<<ix2_z1<<" ix2_z2: "<<ix2_z2<<" trk 1 z0: "<<a.z0<<" trk 2 z0: "<<b.z0<<std::endl;
  float trk1_POCA[2] = {a.d0*sin(a.phi),-1*a.d0*cos(a.phi)};
  float trk2_POCA[2] = {b.d0*sin(b.phi),-1*b.d0*cos(b.phi)};
  float trk1_ix1_delxy[2] = {ix1_x-trk1_POCA[0],ix1_y-trk1_POCA[1]};
  float trk1_ix2_delxy[2] = {ix2_x-trk1_POCA[0],ix2_y-trk1_POCA[1]};
  float trk2_ix1_delxy[2] = {ix1_x-trk2_POCA[0],ix1_y-trk2_POCA[1]};
  float trk2_ix2_delxy[2] = {ix2_x-trk2_POCA[0],ix2_y-trk2_POCA[1]};
  float trk1_traj[2] = {cos(a.phi),sin(a.phi)};
  float trk2_traj[2] = {cos(b.phi),sin(b.phi)};
  bool trk1_ix1_inTraj = ((trk1_ix1_delxy[0]*trk1_traj[0]+trk1_ix1_delxy[1]*trk1_traj[1])>0) ? true : false;
  bool trk1_ix2_inTraj = ((trk1_ix2_delxy[0]*trk1_traj[0]+trk1_ix2_delxy[1]*trk1_traj[1])>0) ? true : false;
  bool trk2_ix1_inTraj = ((trk2_ix1_delxy[0]*trk2_traj[0]+trk2_ix1_delxy[1]*trk2_traj[1])>0) ? true : false;
  bool trk2_ix2_inTraj = ((trk2_ix2_delxy[0]*trk2_traj[0]+trk2_ix2_delxy[1]*trk2_traj[1])>0) ? true : false;
  //std::cout<<"ix1 inTraj: "<<trk1_ix1_inTraj<<" "<<trk2_ix1_inTraj<<" ix2 inTraj: "<<trk1_ix2_inTraj<<" "<<trk2_ix2_inTraj<<std::endl;
  if(trk1_ix1_inTraj&&trk2_ix1_inTraj&&trk1_ix2_inTraj&&trk2_ix2_inTraj){
    if(ix1_delz<ix2_delz){
      x_vtx = ix1_x;
      y_vtx = ix1_y;
      z_vtx = (ix1_z1+ix1_z2)/2;
      return 0;
    }
    else{
      x_vtx = ix2_x;
      y_vtx = ix2_y;
      z_vtx = (ix2_z1+ix2_z2)/2;
      return 0;
    }
  }
  if(trk1_ix1_inTraj&&trk2_ix1_inTraj){
    x_vtx = ix1_x;
    y_vtx = ix1_y;
    z_vtx = (ix1_z1+ix1_z2)/2;
    return 1;
  }
  if(trk1_ix2_inTraj&&trk2_ix2_inTraj){
    x_vtx = ix2_x;
    y_vtx = ix2_y;
    z_vtx = (ix2_z1+ix2_z2)/2;
    return 2;
  }
  else{
    if(ix1_delz<ix2_delz){
      x_vtx = ix1_x;
      y_vtx = ix1_y;
      z_vtx = (ix1_z1+ix1_z2)/2;
      return 3;
    }
    else{
      x_vtx = ix2_x;
      y_vtx = ix2_y;
      z_vtx = (ix2_z1+ix2_z2)/2;
      return 3;
    }
  }
  return 4;
}

Int_t Vertex(Track_Parameters a, Track_Parameters b, Double_t &x_vtx, Double_t &y_vtx, Double_t &z_vtx){
  float x1 = a.x0; //   Centers of the circles
  float y1 = a.y0; // 
  float x2 = b.x0; // 
  float y2 = b.y0; // 
  float R1 = a.rho;   // Radii of the circles
  float R2 = b.rho;
  float del1, del2, x_11, x_12, x_21, x_22, y_11, y_12, y_21, y_22;
  float R = dist(x1,y1,x2,y2); // Distance between centers
  float centerdx = x1 - x2;
  float centerdy = y1 - y2;
  int retval = -1;

  if((R>=(R1-R2)) && (R<=(R1+R2))){
    // Circles Intersect
    float R4 = R*R*R*R;
    float A = (R1*R1 - R2*R2) / (2 * R*R);
    float r2r2 = (R1*R1 - R2*R2);
    float C = TMath::Sqrt(2 * (R1*R1 + R2*R2) / (R*R) - (r2r2 * r2r2) / R4 - 1);
      
    float fx = (x1+x2) / 2 + A * (x2 - x1);
    float gx = C * (y2 - y1) / 2;
    float ix1 = fx + gx;
    float ix2 = fx - gx;

    float fy = (y1+y2) / 2 + A * (y2 - y1);
    float gy = C * (x1 - x2) / 2;
    float iy1 = fy + gy;
    float iy2 = fy - gy;
      
      
    float ap1 = a.phi_T(ix1,iy1);
    float bp1 = b.phi_T(ix1,iy1);
    float ap2 = a.phi_T(ix2,iy2);
    float bp2 = b.phi_T(ix2,iy2);/*
				   float z11 = a->z(ap1);
				   float z21 = b->z(bp1);
				   float z12 = a->z(ap2);
				   float z22 = b->z(bp2);
				 */
    float z11 = a.z(ix1,iy1);
    float z21 = b.z(ix1,iy1);
    float z12 = a.z(ix2,iy2);
    float z22 = b.z(ix2,iy2);
      
    float delz1 = fabs(z11 - z21);//fabs(fabs(ap1 - bp1)-TMath::Pi());// fabs(dist3(ix1,iy1,z11)-dist3(ix1,iy1,z21));
    float delz2 = fabs(z12 - z22);//fabs(fabs(ap2 - bp2)-TMath::Pi());// fabs(dist3(ix2,iy2,z12)-dist3(ix2,iy2,z22));
      
    if(VERBOSE[3]){// &&(fabs(z11-z_vtx)>0.1 && fabs(z12-z_vtx)>0.1 && fabs(z21-z_vtx)>0.1 && fabs(z22-z_vtx)>0.1)){
      std::cout<<Form("ix1 = %5.2f    |   iy1 = %5.2f    |  ix2 = %5.2f    |   iy2 = %5.2f",ix1,iy1,ix2,iy2)<<endl;
      // std::cout<<Form("ap1 = %5.2f    |   bp1 = %5.2f    |  ap2 = %5.2f    |   bp2 = %5.2f",ap1,bp1,ap2,bp2)<<endl;
      std::cout<<Form("z11 = %5.2f    |   z21 = %5.2f    |  z12 = %5.2f    |   z22 = %5.2f    |   delz1 = %5.2f    |   delz2 = %5.2f",z11,z21,z12,z22,delz1,delz2)<<endl;//<<" dxy = "<<tp_dxy->at(it)<<"\t \t dist = "<<TMath::Sqrt((*selectedTPs)[j]->x*(*selectedTPs)[j]->x + (*selectedTPs)[j]->y*(*selectedTPs)[j]->y)<<" \t eta = "<<tp_eta->at(it)<<" \t phi = "<<tp_phi->at(it)<<" \t pt = "<<tp_pt->at(it)<<endl;
    }
    retval = 0;
    if(gx==0 && gy==0){ // Only 1 intersection
      x_vtx = ix1;
      y_vtx = iy1;
      if(delz1<CUTOFF){//(fabs(z11-tp_z)<1.0 || fabs(z21-tp_z)<1.0){//
	retval = 1;
      }
      z_vtx = z11;//fabs(z11)<fabs(z21)?z11:z21;//fabs(a->z0-z11)<fabs(a->z0-z21)?z11:z21;
      // cout<<"----1 intersection ----";
      return retval;//true;
    }
      
    if(dist(ix1,iy1)>20){
      x_vtx = ix2;
      y_vtx = iy2;
      if(delz2<CUTOFF){//(fabs(z12-tp_z)<1.0 || fabs(z22-tp_z)<1.0){//
	retval = 2;
      }
      z_vtx = z12;//fabs(z12)<fabs(z22)?z12:z22;//fabs(a->z0-z12)<fabs(a->z0-z22)?z12:z22;// dist3(ix2,iy2,z12)<dist3(ix2,iy2,z22)?z12:z22;//fabs(a->z0-z12)<fabs(b->z0-z22)?z12:z22;//z12;
      return retval;//true;
    }
    if(dist(ix2,iy2)>20){
      x_vtx = ix1;
      y_vtx = iy1;
      if(delz1<CUTOFF){//(fabs(z11-tp_z)<1.0 || fabs(z21-tp_z)<1.0){//
	retval = 2;
      }
      z_vtx = z11;//fabs(z11)<fabs(z21)?z11:z21;//fabs(a->z0-z11)<fabs(a->z0-z21)?z11:z21;//dist3(ix1,iy1,z11)<dist3(ix1,iy1,z21)?z11:z21;//fabs(a->z0-z11)<fabs(b->z0-z21)?z11:z21;//z11;
      return retval;//true;
    }
      
    // Finding minimum z distance to decide between the 2 intersection vertex

    // if((fabs(z11-tp_z)<1.0 || fabs(z12-tp_z)<1.0 || fabs(z21-tp_z)<1.0 || fabs(z22-tp_z)>1.0)){
    retval = 3;
    // }
 
    if (delz1<=delz2){
      x_vtx = ix1;
      y_vtx = iy1;
      z_vtx = z11;//fabs(z11)<fabs(z21)?z11:z21;// fabs(a->z0-z11)<fabs(a->z0-z21)?z11:z21;
      // cout<<"----2 intersection ----";
      return retval;//true;         
    }
    else{
      x_vtx = ix2;
      y_vtx = iy2;
      z_vtx = z12;//fabs(z12)<fabs(z22)?z12:z22;// fabs(a->z0-z12)<fabs(a->z0-z22)?z12:z22;
      return retval;//true;
    }
      

  }
  else if(R==0){
    // Circles with Same center. 
    x_vtx = -999.0;//a->x();  // Considering First track coordinates as vertex
    y_vtx = -999.0;//a->y();
    return -2;
  }/*
     x_vtx = -99999.0;
     y_vtx = -99999.0;
     z_vtx = -99999.0;
     return false;
     else{
     x_vtx = -9999.0;
     y_vtx = -9999.0;
     return;
     }
   */
   //* Circles don't intersect. 

  if(x1==x2){
    if(y1==y2){
      // Circles with Same center. 
      //x_vtx = a->x();  // Considering First track coordinates as vertex
      // y_vtx = a->y();
      return -2;
    }
    x_11 = x1;
    x_12 = x1;

    x_21 = x2;
    x_22 = x2;

    y_11 = y1 + R1;
    y_12 = y1 - R1;

    y_21 = y2 + R2;
    y_22 = y2 - R2;

  }
  else{

    del1 = R1 / (TMath::Sqrt(1+((y1-y2)*(y1-y2)/((x1-x2)*(x1-x2)))));
    del2 = R2 / (TMath::Sqrt(1+((y1-y2)*(y1-y2)/((x1-x2)*(x1-x2)))));

    x_11 = x1 + del1;
    x_12 = x1 - del1;

    x_21 = x2 + del2;
    x_22 = x2 - del2;

    y_11 = y1 + (x_11-x1) * ((y1-y2)/(x1-x2));
    y_12 = y1 + (x_12-x1) * ((y1-y2)/(x1-x2));

    y_21 = y2 + (x_21-x2) * ((y1-y2)/(x1-x2));
    y_22 = y2 + (x_22-x2) * ((y1-y2)/(x1-x2));
  }

  if(dist(x_11,y_11,x2,y2) <= dist(x_12,y_12,x2,y2)){
    x_vtx = x_11;
    y_vtx = y_11;
  }
  else{
    x_vtx = x_12;
    y_vtx = y_12;
  }

  if(dist(x_21,y_21,x1,y1) <= dist(x_22,y_22,x1,y1)){
    x_vtx = (x_vtx + x_21)/2;
    y_vtx = (y_vtx + y_21)/2;
  }
  else{
    x_vtx = (x_vtx + x_22)/2;
    y_vtx = (y_vtx + y_22)/2;
  }
  //! Does this make sense for no intersection case???
  float az = a.z(x_vtx,y_vtx);
  float bz = b.z(x_vtx,y_vtx);
  z_vtx = (az+bz)/2.0;
  if(VERBOSE[3]){
    std::cout<<Form("z_vtx = %5.1f  |  az = %5.1f  |  bz = %5.1f",z_vtx,az,bz)<<endl;
  }
  return -1;
}

void SetPlotStyle()
{
  // from ATLAS plot style macro

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
  gStyle->SetPaperSize(20, 26);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  gStyle->SetTitleXOffset(1.4);
  gStyle->SetTitleYOffset(1.4);

  // use large fonts
  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.05);
  gStyle->SetLabelFont(42, "x");
  gStyle->SetTitleFont(42, "x");
  gStyle->SetLabelFont(42, "y");
  gStyle->SetTitleFont(42, "y");
  gStyle->SetLabelFont(42, "z");
  gStyle->SetTitleFont(42, "z");
  gStyle->SetLabelSize(0.05, "x");
  gStyle->SetTitleSize(0.05, "x");
  gStyle->SetLabelSize(0.05, "y");
  gStyle->SetTitleSize(0.05, "y");
  gStyle->SetLabelSize(0.05, "z");
  gStyle->SetTitleSize(0.05, "z");

  // use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetHistLineWidth(2.);
  gStyle->SetLineStyleString(2, "[12 12]");

  // get rid of error bar caps
  gStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(0);
  //gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
}

void mySmallText(Double_t x, Double_t y, Color_t color, char *text)
{
  Double_t tsize = 0.044;
  TLatex l;
  l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x, y, text);
}

void removeFlows(TH1F* h)
{
  int nbins = h->GetNbinsX();
  double underflow = h->GetBinContent(0);
  double overflow = h->GetBinContent(nbins+1);
  h->AddBinContent(1,underflow);
  h->AddBinContent(nbins,overflow);
  h->SetBinContent(0,0);
  h->SetBinContent(nbins+1,0);
}

void removeFlows(TH2F* h)
{
}
