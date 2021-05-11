import FWCore.ParameterSet.Config as cms

def customisePFForPixelTracks(process):

    process.hltPFMuonMerging.TrackProducers = cms.VInputTag("hltIterL3MuonTracks", "hltPixelTracksClean")
    process.hltPFMuonMerging.selectedTrackQuals = cms.VInputTag("hltIterL3MuonTracks", "hltPixelTracksClean")

    return process
