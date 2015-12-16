/// \file Utils
/// \brief General auxilliary methods
/// \author ruben.shahoyan@cern.ch

#ifndef ALICEO2_BASE_CONSTANTS
#define ALICEO2_BASE_CONSTANTS

#include "Constants.h"

namespace AliceO2 {
  namespace Base {
    namespace Utils {    

      using AliceO2::Base::Constants;

      inline void  BringTo02Pi(float &phi) { 
	// ensure angle in [0:2pi] for the input in [-pi:pi] or [0:pi]
	if (phi < 0) phi += k2PI; 
      }

      inline void  BringTo02PiGen(float &phi) { 
	// ensure angle in [0:2pi] for the any input angle
	while(phi<0)    {phi += k2PI;}
	while(phi>k2PI) {phi -= k2PI;}	
      }

      inline void  BringToPMPi(float &phi) { 
	// ensure angle in [-pi:pi] for the input in [-pi:pi] or [0:pi]
	if (phi > kPI) phi -= k2PI; 
      }

      inline void  BringToPMPiGen(float &phi) { 
	// ensure angle in [-pi:pi] for any input angle
	while(phi<-kPI)   {phi += k2PI;}
	while(phi> kPI)   {phi -= k2PI;}
      }

      inline vois sincosf(float ang, float& s, float &c) {
	// consider speedup for simultaneus calculation
	s = sinf(ang);
	c = cosf(ang);
      }

    }
  }
}
