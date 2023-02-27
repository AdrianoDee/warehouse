#include <cassert>

#include "MyGPUAnalyzer/MyGPUData/interface/DataSoAHost.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "TTree.h"

class DiMuonSoAAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  DiMuonSoAAnalyzer(edm::ParameterSet const& config)
      : pVToken_(consumes<reco::VertexCollection>(config.getParameter < edm::InputTag > ("primaryVertices"))),
        triggerToken_(consumes<edm::TriggerResults>(config.getParameter < edm::InputTag > ("triggerResults"))),
        hltPaths_(config.getParameter<std::vector<std::string>>("HLTs")),
        diCandToken_{consumes(config.getParameter<edm::InputTag>("diMuonSOA"))},
        candsToken_{consumes(config.getParameter<edm::InputTag>("candsSOA"))},
        diMuonMassCuts_(config.getParameter<std::vector<double>>("diMuonMassCut")),
        keepSameSign_(config.getParameter<bool>("keepSameSign")) 
        {

  
        }

  void beginJob() override {
    
          edm::Service<TFileService> fs;
          candTree_ = fs->make<TTree>("CandidateTree","CandidateTree");

          candTree_->Branch("run",                &run,                "run/I");
          candTree_->Branch("event",              &ev,              "ev/I");
          candTree_->Branch("nCandPerEvent", &nCandPerEvent,  "nCandPerEvent/I");
          candTree_->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/I");
          candTree_->Branch("trigger",            &trigger,            "trigger/I");

          candTree_->Branch("muon1_p4",   "math::XYZTLorentzVector", &muon1_p4);
          candTree_->Branch("muon2_p4",   "math::XYZTLorentzVector", &muon2_p4);
          candTree_->Branch("dimuon_p4",  "math::XYZTLorentzVector", &dimuon_p4);
          
          candTree_->Branch("diMuonMass",      &diMuonMass,          "diMuonMass/F");
          candTree_->Branch("diMuonCharge",      &diMuonCharge,          "diMuonMass/F");
          candTree_->Branch("diMuonDR",      &diMuonDR,          "diMuonDR/F");
  }


  void analyze(const edm::Event& event, const edm::EventSetup& es) override {
    
    TwoDaughterCandHost const& product = event.get(diCandToken_);
    CandidateWithMassHost const& candProduct = event.get(candsToken_);

    auto const& view = product.const_view();
    auto const& candView = candProduct.const_view();

    edm::Handle < reco::VertexCollection  > primaryVertices;
    event.getByToken(pVToken_, primaryVertices);

    edm::Handle < edm::TriggerResults > triggerResults_handle;
    event.getByToken(triggerToken_, triggerResults_handle);

    numPrimaryVertices = primaryVertices->size();
    run = event.id().run();
    ev = event.id().event();

    trigger = 0;

    if (diMuonMassCuts_.size()!=2) {
      throw cms::Exception("Assert") << "diMuonMassCuts input shoud be (min,max),\n got a vector of size"
                                     << diMuonMassCuts_.size() << ".\n";
    }
    
    auto minMass = diMuonMassCuts_[0];
    auto maxMass = diMuonMassCuts_[1];

    if (triggerResults_handle.isValid()) {
        const edm::TriggerNames & TheTriggerNames = event.triggerNames(*triggerResults_handle);
        for (unsigned int i = 0; i < hltPaths_.size(); i++) {
          for (int version = 1; version < 20; version++) {
              std::stringstream ss;
              ss << hltPaths_[i] << "_v" << version;
              unsigned int bit = TheTriggerNames.triggerIndex(edm::InputTag(ss.str()).label());
              if (bit < triggerResults_handle->size() && triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
                trigger += (1<<i);
                break;
              }
          }
        }
      } else std::cout << "*** NO triggerResults found " << event.id().run() << "," << event.id().event() << std::endl;
    
    nCandPerEvent = 0;

    for (int32_t i = 0; i < view.metadata().size(); ++i) {

      auto el = view[i];
      auto mass = el.mass();
      auto charge = el.charge();

      if(mass < minMass or mass > maxMass)
        continue;
      if(not keepSameSign_ and charge!=0)
        continue;

      ++nCandPerEvent;
      dimuon_p4.SetXYZT(el.p().x,el.p().y,el.p().z,mass);
      diMuonMass = mass;
      diMuonCharge = charge;
      diMuonDR = el.dr();

      auto m1 = candView[el.index1()];
      muon1_p4.SetXYZT(m1.p().x,m1.p().y,m1.p().z,m1.mass());

      auto m2 = candView[el.index2()];
      muon2_p4.SetXYZT(m2.p().x,m2.p().y,m2.p().z,m2.mass());

    }

    nCands_ += nCandPerEvent;

    if(nCandPerEvent>0)
      nEvents_++;
    
    candTree_->Fill();

  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("diMuonSOA",edm::InputTag("diMuonGPUProducer"));
    desc.add<edm::InputTag>("candsSOA",edm::InputTag("diMuonGPUProducer"));
    desc.add<std::vector<std::string>>("HLTs",{""});
    desc.add<edm::InputTag>("primaryVertices",edm::InputTag("offlineSlimmedPrimaryVertices"));
    desc.add<edm::InputTag>("triggerResults",edm::InputTag("TriggerResults", "", "HLT"));
    desc.add<bool>("keepSameSign",false);
    desc.add<std::vector<double>>("diMuonMassCut",{2.8,3.3});
    descriptions.addWithDefaultLabel(desc);
  }

  void endJob(){
    std::cout << "###########################" << std::endl;
    std::cout << "DiMuonSoAAnalyzer report:" << std::endl;
    std::cout << "###########################" << std::endl;
    std::cout << "Found " << nEvents_ << " Events with a DiMuon cand." << std::endl;
    std::cout << "###########################" << std::endl;
    std::cout << "Found " << nCands_ << " DiMuon Cands." << std::endl;
    std::cout << "###########################" << std::endl;
 
}


private:
  const edm::EDGetTokenT<reco::VertexCollection>            pVToken_;
  const edm::EDGetTokenT<edm::TriggerResults>               triggerToken_;
  const std::vector<std::string>  hltPaths_;
  const edm::EDGetTokenT<TwoDaughterCandHost> diCandToken_;
  const edm::EDGetTokenT<CandidateWithMassHost> candsToken_;
  const std::vector<double> diMuonMassCuts_;
  const bool keepSameSign_;

  TTree* candTree_;
  unsigned int nCands_ = 0, nEvents_ = 0;

  unsigned int run, ev, nCandPerEvent, numPrimaryVertices, trigger;
  float diMuonMass, diMuonDR, diMuonCharge;

  math::XYZTLorentzVector muon1_p4, muon2_p4, dimuon_p4;

};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DiMuonSoAAnalyzer);