import FWCore.ParameterSet.Config as cms
import os

DisplacedVertexProducer = cms.EDProducer('DisplacedVertexProducer',
  l1TracksInputTag = cms.InputTag("TTTracksFromExtendedTrackletEmulation", "Level1TTTracks"),
  l1VertexCollectionName = cms.string("dispVertices"),
  mcTruthTrackInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigisExtended", "Level1TTTracks"),
  tpInputTag = cms.InputTag("mix", "MergedTrackTruth"),
  qualityAlgorithm = cms.string("NN"),
  ONNXmodel = cms.string("/afs/cern.ch/user/r/rmccarth/private/dispVert/DisplacedMuon/src/model.onnx"),
  ONNXInputName = cms.string("dense_input"),
  featureNames = cms.vstring(["delta_z","R_T","cos_T","d_T","chi2rzdofSum","numStubsSum","chi2rphidofSum","minD0","sumPt"])
)
