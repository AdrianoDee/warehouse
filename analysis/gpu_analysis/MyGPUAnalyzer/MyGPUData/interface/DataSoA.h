#ifndef MyGPUAnalyzer_MyGPUData_interface_DataSoAs_h
#define MyGPUAnalyzer_MyGPUData_interface_DataSoAs_h

#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/SoATemplate/interface/SoAView.h"

namespace myGPUData {

  struct P
  {
    float x;
    float y;
    float z;
  };

}

// A simple candidate
GENERATE_SOA_LAYOUT(CandidateWithMassLayout, 
          SOA_COLUMN(float, mass),
          SOA_COLUMN(float, e),
          SOA_COLUMN(myGPUData::P,p),
          SOA_COLUMN(float, eta),
          SOA_COLUMN(float, phi),
          SOA_COLUMN(float, charge));

// A simple two daughter candidate
GENERATE_SOA_LAYOUT(TwoDaughterCandLayout, 
          SOA_SCALAR(float, mass1), //mass 1
          SOA_SCALAR(float, mass2), //mass 2
          
          SOA_COLUMN(float, mass),
          SOA_COLUMN(float, e),
          SOA_COLUMN(float, charge),

          SOA_COLUMN(float, dr),

          SOA_COLUMN(myGPUData::P,p),
          SOA_COLUMN(myGPUData::P,p1), //daughter 1 p
          SOA_COLUMN(myGPUData::P,p2), //daughter 2 p

          SOA_COLUMN(uint32_t,index1),
          SOA_COLUMN(uint32_t,index2));

using CandidateWithMass = CandidateWithMassLayout<>;
using TwoDaughterCand = TwoDaughterCandLayout<>;
  

#endif