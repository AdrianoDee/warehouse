#ifndef MyGPUAnalyzer_MyGPUData_interface_alpaka_DataSoADevice_h
#define MyGPUAnalyzer_MyGPUData_interface_alpaka_DataSoADevice_h

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "MyGPUAnalyzer/MyGPUData/interface/DataSoA.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
    
    using CandidateWithMassDevice = PortableCollection<CandidateWithMass>;
    using TwoDaughterCandDevice = PortableCollection<TwoDaughterCand>;

}
#endif