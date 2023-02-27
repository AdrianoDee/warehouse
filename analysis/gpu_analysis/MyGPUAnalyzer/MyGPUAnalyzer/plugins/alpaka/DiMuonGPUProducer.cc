#include "DataFormats/Portable/interface/Product.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/global/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

#include "MyGPUAnalyzer/MyGPUData/interface/DataSoAHost.h"
#include "MyGPUAnalyzer/MyGPUData/interface/alpaka/DataSoADevice.h"

#include "TwoCandsAlgo.h"

#include <alpaka/alpaka.hpp>

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class DiMuonGPUProducer : public global::EDProducer<> {

  public:
    DiMuonGPUProducer(edm::ParameterSet const& config)
        : deviceDiMuons_{produces()}, 
          deviceCands_{produces()},
          // ptCut_(config.getParameter<double>("ptCut_")) 
          muonCollection_(consumes(config.getParameter<edm::InputTag>("muons")))
    {}
    
    void produce(edm::StreamID sid, device::Event& event, device::EventSetup const&) const override {

      using namespace myGPUData;

      const auto cands = event.get(muonCollection_);
      
      int size = cands.size();
      auto queue = event.queue();
      CandidateWithMassHost hostCands(size, queue);

      auto candViewHost = hostCands.view();
      int twoSize = 0, realSize = 0;

      //Let's fill the portable SoA on host (CPU) ... 
      for (int i = 0; i < size; i++)
      {
          auto cand = cands.at(i);

          // if (cand.pt() < ptCut_)
          //   continue
          
          twoSize += i;
          
          candViewHost[i].mass() = muonmass;

          candViewHost[i].charge() = cand.charge();

          candViewHost[i].p().x = cand.px();
          candViewHost[i].p().y = cand.py();
          candViewHost[i].p().z = cand.pz();

          candViewHost[i].e() = cand.energy();

          // std::cout << i << " - " << candViewHost[i].mass() << " - "
          //           << candViewHost[i].p().x << " - "
          //           << candViewHost[i].p().y << " - "
          //           << candViewHost[i].p().z << " - "
          //           << candViewHost[i].e() << " - "
          //           << candViewHost[i].charge() << std::endl;
      }
      
      CandidateWithMassDevice deviceCands(size, queue);
      //... and copy it on device
      alpaka::memcpy(queue, deviceCands.buffer(), hostCands.buffer());

      TwoDaughterCandDevice deviceDiMuons(twoSize, queue);
      // run the dimuon builder on device (GPU)
      algo_.run(event.queue(), deviceCands, deviceDiMuons, muonmass, muonmass);

      // put the device product into the event (will be copied on host transparentely)
      event.emplace(deviceDiMuons_, std::move(deviceDiMuons));
      event.emplace(deviceCands_, std::move(deviceCands));
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("muons",edm::InputTag("slimmedMuons"));
      descriptions.addWithDefaultLabel(desc);
    }

  private:

    constexpr static float muonmass = 0.1056583715;

    const device::EDPutToken<TwoDaughterCandDevice> deviceDiMuons_;
    const device::EDPutToken<CandidateWithMassDevice> deviceCands_;

    // const double ptCut_;
    const edm::EDGetTokenT<pat::MuonCollection> muonCollection_;

    // implementation of the dimuon builder algorithm
    TwoCandsAlgo algo_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(DiMuonGPUProducer);