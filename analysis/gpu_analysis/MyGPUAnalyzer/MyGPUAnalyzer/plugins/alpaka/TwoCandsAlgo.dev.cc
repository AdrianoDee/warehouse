// Check that ALPAKA_HOST_ONLY is not defined during device compilation:
// #ifdef ALPAKA_HOST_ONLY
// #error ALPAKA_HOST_ONLY defined in device compilation
// #endif

#include <alpaka/alpaka.hpp>

#include "MyGPUAnalyzer/MyGPUData/interface/alpaka/DataSoADevice.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/traits.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "TwoCandsAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace cms::alpakatools;

  namespace utilities
  {
    constexpr float m_pi = 3.14159265358979f;
    
    using ConstView = CandidateWithMass::ConstView;  

    ALPAKA_FN_ACC float deltaR(const ConstView &view, int32_t i, int32_t j)
    {
      float dp = std::abs(view[i].phi()-view[j].phi());
      if (dp>float(m_pi)) dp -= 2.f*m_pi;
      float de = view[i].eta() - view[j].eta();
      
      return sqrt(de*de + dp*dp);
    }
  }

  class TwoCandKernel {

  public:
    
    using ConstView = CandidateWithMass::ConstView;  
    
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  CandidateWithMassDevice::ConstView input_view,
                                  TwoDaughterCand::View output_view,
                                  double candMass1,
                                  double candMass2,
                                  int32_t size) const {
      // global index of the thread within the grid
      const int32_t thread = alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0u];

      // set this only once in the whole kernel grid
      if (thread == 0) {
        output_view.mass1() = candMass1;
        output_view.mass2() = candMass2;
      }
      
      auto sqM1 = candMass1 * candMass1;
      auto sqM2 = candMass2 * candMass2;
      // make a strided loop over the kernel grid, covering up to "size" elements
      for (int32_t i : elements_with_stride(acc, size)) {

        auto cand1 = input_view[i];

        float px1 = cand1.p().x;
        float py1 = cand1.p().y;
        float pz1 = cand1.p().z;

        float e1 = sqrt(sqM1 + px1*px1 + py1*py1 + pz1*pz1);

        float q1 = cand1.charge();

        auto start = candMass1==candMass2 ? i + 1 : 0;

        for (int j = start; j < size; j++)
        {
          if(i==j)
            continue;

          int I = (candMass1==candMass2) ? (size*(size-1)/2) - (size-i)*((size-i)-1)/2 + j - i - 1 : (size*i) + j;

          auto cand2 = input_view[j];

          float px2 = cand2.p().x;
          float py2 = cand2.p().y;
          float pz2 = cand2.p().z;
          
          float e2 = sqrt(sqM2 + cand2.p().x*cand2.p().x + cand2.p().y*cand2.p().y + cand2.p().z*cand2.p().z);
          
          auto twoCand = output_view[I];

          float e = e1 + e2;
          
          twoCand.e() = e;
          twoCand.dr() = utilities::deltaR(input_view,i,j);
          
          twoCand.p().x = px1 + px2;
          twoCand.p().y = py1 + py2;
          twoCand.p().z = pz1 + pz2;

          twoCand.charge() = q1 + cand2.charge();

          twoCand.mass() = e*e - twoCand.p().x*twoCand.p().x - twoCand.p().y*twoCand.p().y - twoCand.p().z*twoCand.p().z;
          
          twoCand.index1() = i;
          twoCand.index2() = j;

        }
      }
    }

  };

  void TwoCandsAlgo::run(Queue& queue, 
            CandidateWithMassDevice& input, 
            TwoDaughterCandDevice& output, 
            double candMass1,
            double candMass2) const {
    // use 64 items per group (this value is arbitrary, but it's a reasonable starting point)
    uint32_t items = 64;

    // use as many groups as needed to cover the whole problem
    uint32_t groups = divide_up_by(input->metadata().size(), items);

    // map items to
    //   - threads with a single element per thread on a GPU backend
    //   - elements within a single thread on a CPU backend
    auto workDiv = make_workdiv<Acc1D>(groups, items);

    alpaka::exec<Acc1D>(queue, workDiv, TwoCandKernel{}, input.const_view(), output.view(), 
                        candMass1, candMass2, input->metadata().size());
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE