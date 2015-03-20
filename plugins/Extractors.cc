#include <Extractors/PatExtractor/interface/ElectronExtractor.h>
#include <Extractors/PatExtractor/interface/EventExtractor.h>
#include <Extractors/PatExtractor/interface/HLTExtractor.h>
#include <Extractors/PatExtractor/interface/JetMETExtractor.h>
#include <Extractors/PatExtractor/interface/MCExtractor.h>
#include <Extractors/PatExtractor/interface/MuonExtractor.h>
#include <Extractors/PatExtractor/interface/PFpartExtractor.h>
#include <Extractors/PatExtractor/interface/PhotonExtractor.h>
#include <Extractors/PatExtractor/interface/TrackExtractor.h>
#include <Extractors/PatExtractor/interface/VertexExtractor.h>

DEFINE_EDM_PLUGIN(PatExtractorExtractorFactory, ElectronExtractor, "electron_extractor");
DEFINE_EDM_PLUGIN(PatExtractorExtractorReadOnlyFactory, ElectronExtractor, "electron_extractor");

DEFINE_EDM_PLUGIN(PatExtractorExtractorFactory, EventExtractor, "event_extractor");
DEFINE_EDM_PLUGIN(PatExtractorExtractorReadOnlyFactory, EventExtractor, "event_extractor");

DEFINE_EDM_PLUGIN(PatExtractorExtractorFactory, HLTExtractor, "hlt_extractor");
DEFINE_EDM_PLUGIN(PatExtractorExtractorReadOnlyFactory, HLTExtractor, "hlt_extractor");

DEFINE_EDM_PLUGIN(PatExtractorExtractorFactory, JetMETExtractor, "jet_met_extractor");
DEFINE_EDM_PLUGIN(PatExtractorExtractorReadOnlyFactory, JetMETExtractor, "jet_met_extractor");

DEFINE_EDM_PLUGIN(PatExtractorExtractorFactory, MCExtractor, "mc_extractor");
DEFINE_EDM_PLUGIN(PatExtractorExtractorReadOnlyFactory, MCExtractor, "mc_extractor");

DEFINE_EDM_PLUGIN(PatExtractorExtractorFactory, MuonExtractor, "muon_extractor");
DEFINE_EDM_PLUGIN(PatExtractorExtractorReadOnlyFactory, MuonExtractor, "muon_extractor");

DEFINE_EDM_PLUGIN(PatExtractorExtractorFactory, PFpartExtractor, "pfparticle_extractor");
DEFINE_EDM_PLUGIN(PatExtractorExtractorReadOnlyFactory, PFpartExtractor, "pfparticle_extractor");

DEFINE_EDM_PLUGIN(PatExtractorExtractorFactory, PhotonExtractor, "photon_extractor");
DEFINE_EDM_PLUGIN(PatExtractorExtractorReadOnlyFactory, PhotonExtractor, "photon_extractor");

DEFINE_EDM_PLUGIN(PatExtractorExtractorFactory, TrackExtractor, "track_extractor");
DEFINE_EDM_PLUGIN(PatExtractorExtractorReadOnlyFactory, TrackExtractor, "track_extractor");

DEFINE_EDM_PLUGIN(PatExtractorExtractorFactory, VertexExtractor, "vertex_extractor");
DEFINE_EDM_PLUGIN(PatExtractorExtractorReadOnlyFactory, VertexExtractor, "vertex_extractor");
