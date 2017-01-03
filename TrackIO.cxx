#include "TrackIO.h"

using namespace AliceO2::Base;

ClassImp(TrackIO::Track);

void TrackIO::Track::Print(Option_t*) const
{
  mTrack.Print();
  printf("Chi2: %f\n",mChi2);
}
  
