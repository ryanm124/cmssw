#ifndef __L1Trigger_TrackFindingTracklet_DisplacedVertexProducer_h__
#define __L1Trigger_TrackFindingTracklet_DisplacedVertexProducer_h__

#include "DataFormats/L1Trigger/interface/DisplacedVertex.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"
#include "TMath.h" 
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <valarray>
#include <deque>

using namespace std;

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
  float tp_pt;
  float x0;
  float y0;
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
    return (z0+t*r); // can do higher order terms if necessary from displaced math
  }
  Track_Parameters(float pt_in, float d0_in, float z0_in, float eta_in, float phi_in, float charge_in, int index_in = -1, int pdgid_in = -99999, float tp_pt_in = 0){
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
    tp_pt = tp_pt_in;
    rho = 100*pt / (0.3 * 3.8);
    x0 =  (d0 + charge * rho)*TMath::Sin(phi);
    y0 = -(d0 + charge * rho)*TMath::Cos(phi);
  }
  ~Track_Parameters(){};
};

class Vertex_Parameters
{
 public:
  Double_t x_dv;
  Double_t y_dv;
  Double_t z_dv;
  Track_Parameters a;
  Track_Parameters b;
  bool matched = false;
  int errorCode = 0;
 Vertex_Parameters(Double_t x_dv_in, Double_t y_dv_in, Double_t z_dv_in, Track_Parameters a_in, Track_Parameters b_in):
  a(a_in),
  b(b_in)
  {
    x_dv = x_dv_in;
    y_dv = y_dv_in;
    z_dv = z_dv_in;
  }
  ~Vertex_Parameters(){};
};

class DisplacedVertexProducer : public edm::global::EDProducer<> {
 public:
  explicit DisplacedVertexProducer(const edm::ParameterSet&);
  ~DisplacedVertexProducer() override {}

 private:
  typedef edm::View<TTTrack<Ref_Phase2TrackerDigi_>> TTTrackCollectionView;

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

 private:
  const edm::EDGetTokenT<TTTrackAssociationMap<Ref_Phase2TrackerDigi_>> ttTrackMCTruthToken_;
  const edm::EDGetTokenT<TTTrackCollectionView> l1TracksToken_;
  const std::string outputCollectionName_;
  const std::string qualityAlgorithm_;
  const std::string ONNXmodel_;
  const std::string ONNXInputName_;
  const std::vector<std::string> featureNames_;
};

#endif
