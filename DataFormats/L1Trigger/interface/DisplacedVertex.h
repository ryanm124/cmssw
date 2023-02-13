#ifndef DataFormats_L1TVertex_DisplacedVertex_h
#define DataFormats_L1TVertex_DisplacedVertex_h

#include <vector>

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"

namespace l1t {

  class DisplacedVertex {
  public:
    typedef TTTrack<Ref_Phase2TrackerDigi_> Track_t;
  DisplacedVertex(const std::vector<edm::Ptr<Track_t>>& tracks, float delta_z, float R_T, float cos_T, float d_T, bool isPV) : tracks_(tracks), delta_z_(delta_z), R_T_(R_T), cos_T_(cos_T), d_T_(d_T), isPV_(isPV) {
      float maxPt = 0;
      float maxPt_i = 0;
      float chi2rzdofSum = 0;
      int numStubsSum = 0;
      float chi2rphidofSum = 0;
      float minD0 = 9999;
      float sumPt = 0;
      for(uint i=0; i<tracks.size(); i++){
	float pt = tracks[i]->momentum().perp();
	std::vector<edm::Ref<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_> >, TTStub<Ref_Phase2TrackerDigi_> > > 
	  stubRefs = tracks[i]->getStubRefs();
	int nstubs = stubRefs.size();
	chi2rzdofSum += tracks[i]->chi2Z() / (2*nstubs-5);
	numStubsSum += nstubs;
	chi2rphidofSum += tracks[i]->chi2XY() / (2*nstubs-5);
	float phi = tracks[i]->momentum().phi();
	float d0 = tracks[i]->POCA().x() * sin(phi) - tracks[i]->POCA().y() * cos(phi);
	sumPt += pt;
	if(pt>maxPt){
	  maxPt = pt;
	  maxPt_i = i;
	}
	if(d0<minD0){
	  minD0 = d0;
	}
      }
      pt_ = maxPt;
      z0_ = tracks[maxPt_i]->z0();
      eta_ = tracks[maxPt_i]->momentum().eta();
      phi_ = tracks[maxPt_i]->momentum().phi();
      d0_ = tracks[maxPt_i]->POCA().x() * sin(phi_) - tracks[maxPt_i]->POCA().y() * cos(phi_);
      chi2rzdofSum_ = chi2rzdofSum;
      numStubsSum_ = numStubsSum;
      chi2rphidofSum_ = chi2rphidofSum;
      minD0_ = minD0;
      sumPt_ = sumPt;
    }
    DisplacedVertex() {}
    ~DisplacedVertex() {}

    float pt() const { return pt_; }
    float z0() const { return z0_; }
    float d0() const { return d0_; }
    float eta() const { return eta_; }
    float phi() const { return phi_; }
    float delta_z() const { return delta_z_; }
    float R_T() const { return R_T_; }
    float cos_T() const { return cos_T_; }
    float d_T() const { return d_T_; }
    float chi2rzdofSum() const { return chi2rzdofSum_; }
    int numStubsSum() const { return numStubsSum_; }
    float chi2rphidofSum() const { return chi2rphidofSum_; }
    float minD0() const { return minD0_; }
    float sumPt() const { return sumPt_; }
    void setScore(float score) { score_ = score; }
    float score() const { return score_; }
    void setPos(std::vector<float> pos) { pos_ = pos; }
    std::vector<float> pos() const { return pos_; }
    const std::vector<edm::Ptr<Track_t>>& tracks() const { return tracks_; }
    bool isPV() const { return isPV_; }

  private:
    std::vector<edm::Ptr<Track_t>> tracks_;
    float pt_;
    float z0_;
    float d0_;
    float eta_;
    float phi_;
    float delta_z_;
    float R_T_;
    float cos_T_;
    float d_T_;
    float chi2rzdofSum_;
    int numStubsSum_;
    float chi2rphidofSum_;
    float minD0_;
    float sumPt_;
    float score_ = -1.0;
    bool isPV_;
    std::vector<float> pos_ = {0,0,0};
    
  };

  typedef std::vector<DisplacedVertex> DisplacedVertexCollection;

}  // namespace l1t

#endif
