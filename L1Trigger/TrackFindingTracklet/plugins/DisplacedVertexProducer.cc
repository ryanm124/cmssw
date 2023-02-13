#include "L1Trigger/TrackFindingTracklet/interface/DisplacedVertexProducer.h"

bool ComparePtTrack(edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_>> a, edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_>> b) { return a->momentum().perp() > b->momentum().perp(); }

Double_t dist(Double_t x1, Double_t y1 , Double_t x2=0, Double_t y2=0){ // Distance between 2 points
  return (TMath::Sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)));
}

std::valarray<float> calcPVec(Track_Parameters a, double_t v_x, double_t v_y){
  std::valarray<float> r_vec = {float(v_x-a.x0),float(v_y-a.y0)};
  std::valarray<float> p_vec = {-r_vec[1],r_vec[0]};
  if(a.charge>0){
    p_vec *= -1;
  }
  p_vec /= TMath::Sqrt(pow(p_vec[0],2)+pow(p_vec[1],2));
  p_vec *= a.pt;
  return p_vec;
}

Double_t dist_TPs(Track_Parameters a, Track_Parameters b){
  float x1 = a.x0; //   Centers of the circles
  float y1 = a.y0; //
  float x2 = b.x0; //
  float y2 = b.y0; //
  float R1 = a.rho;// Radii of the circles
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

Int_t Vertex(Track_Parameters a, Track_Parameters b, float &x_vtx, float &y_vtx, float &z_vtx){
  float CUTOFF = 1.0;
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
    float bp2 = b.phi_T(ix2,iy2);
    float z11 = a.z(ix1,iy1);
    float z21 = b.z(ix1,iy1);
    float z12 = a.z(ix2,iy2);
    float z22 = b.z(ix2,iy2);
    float delz1 = fabs(z11 - z21);
    float delz2 = fabs(z12 - z22);
    retval = 0;
    if(gx==0 && gy==0){ // Only 1 intersection
      x_vtx = ix1;
      y_vtx = iy1;
      if(delz1<CUTOFF){
	retval = 1;
      }
      z_vtx = z11;
      return retval;
    }
    if(dist(ix1,iy1)>20){
      x_vtx = ix2;
      y_vtx = iy2;
      if(delz2<CUTOFF){
	retval = 2;
      }
      z_vtx = z12;
      return retval;
    }
    if(dist(ix2,iy2)>20){
      x_vtx = ix1;
      y_vtx = iy1;
      if(delz1<CUTOFF){
	retval = 2;
      }
      z_vtx = z11;
      return retval;
    }
    else{
      x_vtx = ix2;
      y_vtx = iy2;
      z_vtx = z12;
      return retval;
    }
  }
  else if(R==0){
    // Circles with Same center.
    x_vtx = -999.0;
    y_vtx = -999.0;
    return -2;
  }
  if(x1==x2){
    if(y1==y2){
      // Circles with Same center.
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
  float az = a.z(x_vtx,y_vtx);
  float bz = b.z(x_vtx,y_vtx);
  z_vtx = (az+bz)/2.0;
  return -1;
}

DisplacedVertexProducer::DisplacedVertexProducer(const edm::ParameterSet& iConfig)
  : ttTrackMCTruthToken_(consumes<TTTrackAssociationMap<Ref_Phase2TrackerDigi_> >(iConfig.getParameter<edm::InputTag>("mcTruthTrackInputTag"))),
    l1TracksToken_(consumes<TTTrackCollectionView>(iConfig.getParameter<edm::InputTag>("l1TracksInputTag"))),
    outputCollectionName_(iConfig.getParameter<std::string>("l1VertexCollectionName")),
    qualityAlgorithm_(iConfig.getParameter<std::string>("qualityAlgorithm")),
    ONNXmodel_(iConfig.getParameter<std::string>("ONNXmodel")),
    ONNXInputName_(iConfig.getParameter<std::string>("ONNXInputName")),
    featureNames_(iConfig.getParameter<std::vector<std::string>>("featureNames"))
{
  //--- Define EDM output to be written to file (if required)
  produces<l1t::DisplacedVertexCollection>(outputCollectionName_);
}

void DisplacedVertexProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  edm::Handle<TTTrackAssociationMap<Ref_Phase2TrackerDigi_> > MCTruthTTTrackHandle;
  iEvent.getByToken(ttTrackMCTruthToken_, MCTruthTTTrackHandle);
  edm::Handle<TTTrackCollectionView> l1TracksHandle;
  iEvent.getByToken(l1TracksToken_, l1TracksHandle);
  std::deque<edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_>>> selectedTracks;
  float TP_maxEta = 2.4;
  float TP_minPt = 2.0;
  float TP_maxPt = 10000.0;
  float TP_maxD0 = 10.0;
  float TP_minD0 = 0.0004196;
  float TP_maxZ0 = 20.0;
  float maxChi2rzdof = 5.0;
  float maxBendChi2 = 9.0;

  for (const auto& track : l1TracksHandle->ptrs()) {
    float pt = track->momentum().perp();
    float eta = track->eta();
    float phi = track->phi();
    float z0 = track->z0();
    float d0 = track->d0();
    std::vector<edm::Ref<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_> >, TTStub<Ref_Phase2TrackerDigi_> > > stubRefs = track->getStubRefs();
    float chi2rzdof = track->chi2Z() / (2*stubRefs.size()-5);
    float bendchi2 = track->stubPtConsistency();
    if(fabs(eta)<TP_maxEta && pt>TP_minPt && pt<TP_maxPt && fabs(d0)<TP_maxD0 && fabs(d0)>TP_minD0 && z0<TP_maxZ0 && chi2rzdof<maxChi2rzdof && bendchi2<maxBendChi2){
      selectedTracks.push_back(track);  
    }
  }
  sort(selectedTracks.begin(), selectedTracks.end(), ComparePtTrack);
  if(selectedTracks.size()>20){
    selectedTracks.erase(selectedTracks.begin()+20,selectedTracks.end());
  } 
  float x_dv_trk = -9999.0;
  float y_dv_trk = -9999.0;
  float z_dv_trk = -9999.0;
  std::unique_ptr<l1t::DisplacedVertexCollection> product(new std::vector<l1t::DisplacedVertex>());
  while(selectedTracks.size()>1){
    for(uint j=1;j<selectedTracks.size();j++){
      Track_Parameters trackParams_0 = Track_Parameters(selectedTracks[0]->momentum().perp(), selectedTracks[0]->d0(), selectedTracks[0]->z0(), selectedTracks[0]->eta(), selectedTracks[0]->phi(), -selectedTracks[0]->rInv());
      Track_Parameters trackParams_j = Track_Parameters(selectedTracks[j]->momentum().perp(), selectedTracks[j]->d0(), selectedTracks[j]->z0(), selectedTracks[j]->eta(), selectedTracks[j]->phi(), -selectedTracks[j]->rInv());
      if( dist_TPs( trackParams_0, trackParams_j ) == 0 ){
	edm::Ptr<TrackingParticle> tp_0 = MCTruthTTTrackHandle->findTrackingParticlePtr(selectedTracks[0]);
	edm::Ptr<TrackingParticle> tp_j = MCTruthTTTrackHandle->findTrackingParticlePtr(selectedTracks[j]);
	bool isPV = false;
	if(!tp_0.isNull() && !tp_j.isNull()){
	  if(tp_0->eventId().event()==0 && tp_j->eventId().event()==0 && tp_0->vx()==tp_j->vx() && tp_0->vy()==tp_j->vy() && tp_0->vz()==tp_j->vz()){
	    isPV = true;
	  }
	}
	int Vertex_check = Vertex(trackParams_0,trackParams_j,x_dv_trk,y_dv_trk,z_dv_trk);
	float delta_z = fabs(trackParams_0.z(x_dv_trk,y_dv_trk)-trackParams_j.z(x_dv_trk,y_dv_trk));
	std::valarray<float> p_trk_1 = calcPVec(trackParams_0,x_dv_trk,y_dv_trk);
	std::valarray<float> p_trk_2 = calcPVec(trackParams_j,x_dv_trk,y_dv_trk);
	std::valarray<float> p_tot = p_trk_1+p_trk_2;
	float R_T = TMath::Sqrt(pow(x_dv_trk,2)+pow(y_dv_trk,2));
	float cos_T = (p_tot[0]*x_dv_trk+p_tot[1]*y_dv_trk)/(R_T*TMath::Sqrt(pow(p_tot[0],2)+pow(p_tot[1],2)));
	float theta = atan2(p_tot[1],p_tot[0]);
	float d_T = std::abs(cos(theta)*y_dv_trk-sin(theta)*x_dv_trk);
	std::vector vertTracks = {selectedTracks[0],selectedTracks[j]};
	l1t::DisplacedVertex vert(vertTracks, delta_z, R_T, cos_T, d_T, isPV);
	vert.setPos({x_dv_trk,y_dv_trk,z_dv_trk});
	std::vector<std::string> ortinput_names;
	std::vector<std::string> ortoutput_names;
	cms::Ort::FloatArrays ortinput;
	cms::Ort::FloatArrays ortoutputs;
	std::vector<float> Transformed_features = {delta_z, R_T, cos_T, d_T, vert.chi2rzdofSum(), float(vert.numStubsSum()), vert.chi2rphidofSum(), vert.minD0(), vert.sumPt()};
	cms::Ort::ONNXRuntime Runtime(this->ONNXmodel_);  //Setup ONNX runtime
	ortinput_names.push_back(this->ONNXInputName_);
	ortoutput_names = Runtime.getOutputNames();
	ortinput.push_back(Transformed_features);
	int batch_size = 1;
	ortoutputs = Runtime.run(ortinput_names, ortinput, {}, ortoutput_names, batch_size);
	vert.setScore(ortoutputs[0][0]);
	product->emplace_back(vert);
      }
    }
    selectedTracks.pop_front();
  }

  // //=== Store output
  iEvent.put(std::move(product), outputCollectionName_);
}

DEFINE_FWK_MODULE(DisplacedVertexProducer);
