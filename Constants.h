/// \file Constants
/// \brief General constants
/// \author ruben.shahoyan@cern.ch

#ifndef ALICEO2_BASE_CONSTANTS
#define ALICEO2_BASE_CONSTANTS

namespace AliceO2 {
  namespace Base {
    namespace Constants {    
      
      const float kAlmost0 = 1.17549e-38;
      const float kAlmost1 = 1.f-kAlmost0;
      const float kVeryBig = 1.f/kAlmost0;

      const float kPI     = 3.14159274101257324e+00f;
      const float k2PI    = 2.f*kPI;
      const float kPIhalf = 0.5f*kPI;      

      // conversion from B(kGaus) to curvature for 1GeV pt
      const float kB2C     = -0.299792458e-3;

    }
  }
}
#endif