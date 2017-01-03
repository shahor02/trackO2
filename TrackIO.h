/// \file TrackIO
/// \brief Test persistent track model
/// \author ruben.shahoyan@cern.ch
#ifndef ALICEO2_BASE_TRACKIO
#define ALICEO2_BASE_TRACKIO

#include "TObject.h"
#include "Track.h"

namespace AliceO2 {
  namespace Base {
    namespace TrackIO {

      using namespace AliceO2::Base::Track;
      
      class Track : public TObject {
      public:
      Track() : mTrack(), mChi2(0) {} 
      Track(const TrackParCov& trc, float chi2=0.f) : mTrack(trc), mChi2(chi2) {} 
	const TrackParCov& GetTrack() const {return mTrack;}
	float GetChi2()               const {return mChi2;}
	void Print(Option_t* opt="")  const;
      protected:
	TrackParCov mTrack;
	Float_t     mChi2;
	//
	ClassDef(Track,1)
      }; 

    }
  }
}

#endif
