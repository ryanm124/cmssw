import FWCore.ParameterSet.Config as cms

DisplacedVertexProducer = cms.EDProducer('DisplacedVertexProducer',
  l1TracksInputTag = cms.InputTag("l1tTTTracksFromExtendedTrackletEmulation", "Level1TTTracks"),
  l1TrackVertexCollectionName = cms.string("dispVertices"),
  mcTruthTrackInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigisExtended", "Level1TTTracks"),
  ONNXmodel = cms.string("/afs/cern.ch/user/r/rmccarth/private/dispVert/CMSSW_14_0_0_pre2/src/L1Trigger/L1TTrackMatch/test/model.onnx"),
  ONNXInputName = cms.string("dense_input"),
  featureNames = cms.vstring(["delta_z","R_T","cos_T","d_T","chi2rzdofSum","numStubsSum","chi2rphidofSum","minD0","sumPt"])
)
