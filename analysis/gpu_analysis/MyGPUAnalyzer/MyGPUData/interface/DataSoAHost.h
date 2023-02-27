#ifndef MyGPUAnalyzer_MyGPUData_interface_DataSoAHost_h
#define MyGPUAnalyzer_MyGPUData_interface_DataSoAHost_h

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataSoA.h"

using CandidateWithMassHost = PortableHostCollection<CandidateWithMass>;
using TwoDaughterCandHost = PortableHostCollection<TwoDaughterCand>;

#endif