#ifndef MyGPUAnalyzer_MyGPUAnalyzer_plugins_alpaka_TwoCandsAlgo_h
#define MyGPUAnalyzer_MyGPUAnalyzer_plugins_alpaka_TwoCandsAlgo_h

#include "MyGPUAnalyzer/MyGPUData/interface/alpaka/DataSoADevice.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class TwoCandsAlgo {
  public:
    void run(Queue& queue, 
            CandidateWithMassDevice& input, 
            TwoDaughterCandDevice& output, 
            double candMass1,
            double candMass2) const;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // MyGPUAnalyzer_MyGPUAnalyzer_plugins_alpaka_TwoCandsAlgo_h