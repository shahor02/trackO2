#include "Track.h"
#include "TrackIO.h"
#include "TFile.h"
#include <math.h>
#include <memory>

using namespace AliceO2::Base::Track;
using namespace AliceO2::Base::TrackIO;
using namespace std;

void CompPar(const TrackParBase* t0,const TrackParBase* t1);

// Test to run under root6, w/o aliroot
// Load 1st compiled library by root [0] gSystem->Load("libTrackO2.so");
void testR6()
{
  const float  xyzF[3] = {0,0,0}, pxyzF[3] = {0.5,0.5,0.5};
  const float covF[21]  = {1e-6,
			   0,1e-6,
			   0,0,1e-6,
			   0,0,0,1e-4,
			   0,0,0,0,1e-4,
			   0,0,0,0,0,1e-4};

  unique_ptr<TrackParCov>  trcO2(new TrackParCov(xyzF,pxyzF,covF,1,false));
  unique_ptr<TrackPar> trpO2(new TrackPar(xyzF,pxyzF,1,false));
  //
  printf("\nAt creation:\n");
  printf("[0] TrackParCov:\n");
  trcO2->Print();
  printf("[1] TrackPar:\n");
  trpO2->PrintParam();
  //
  const float xmax = 400.,xstep=4.;;
  float bF[3] = {0.01,0.01,5.0};
  const float x2x0 = 0.01;
  const float xrho = xstep*2.0;
  const float mass = 0.14;
  //
  float x = trcO2->GetX();
  int cnt = 0;
  float chi0=0,chi3=0;
  float measErrF[3]={1e-6,0.25e-5,1e-4};
  int nupd = 0;
  while (x<xmax) {
    x += xstep;
    if (!trcO2->PropagateTo(x,bF[2]))           {printf("[0] Fail on prop.%f\n",x); exit(1);}
    if (!trpO2->PropagateParamTo(x,bF[2]))      {printf("[1] Fail on prop.%f\n",x); exit(1);}
    cnt++;
    if (cnt%2==0) {
      if (!trcO2->CorrectForMaterial(x2x0,xrho,mass,false))  {printf("[0] Fail on matcor. %f\n",x); exit(1);}
    }

    if (cnt%3==0) {
      float  measF[2] = {trcO2->GetY()+sqrtf(measErrF[0]),trcO2->GetZ()+sqrtf(measErrF[2])};
      chi0 += trcO2->GetPredictedChi2(measF,measErrF);
      if (!trcO2->Update(measF,measErrF))  {printf("[0] Fail on update. %f\n",x); exit(1);}
      nupd++;
    }
    //

    if (cnt%4 == 0 || fabs(x-xmax)<1e-4f) {
      float alp = trcO2->GetPhi();
      if (!trcO2->Rotate(alp))       {printf("[0] Fail on rot.%f\n",x); exit(1);}
      if (!trpO2->RotateParam(alp))  {printf("[1] Fail on rot.%f\n",x); exit(1);}
    }

  }

  printf("\nAt Max X:\n");
  printf("[0] TrackParCov:\n");
  trcO2->Print();
  printf("[1] TrackPar:\n");
  trpO2->PrintParam();

  printf("\nDifference between full and param-only propagation\n");
  CompPar(trpO2.get(),trcO2.get());

  // propagate back with 3d field
  while (x>0) {
    x -= xstep;
    if (!trcO2->PropagateTo(x,bF))                 {printf("[0] Fail on prop3d.%f\n",x); exit(1);}
    if (!trpO2->PropagateParamTo(x,bF))            {printf("[1] Fail on prop3d.%f\n",x); exit(1);}
  }

  // bring to the same frame
  float alp = trcO2->GetPhi();
  if (!trcO2->Rotate(alp))       {printf("[0] Fail on rot.%f\n",x); exit(1);}
  if (!trpO2->RotateParam(alp))  {printf("[1] Fail on rot.%f\n",x); exit(1);}
  //
  x = 0;
  if (!trcO2->PropagateTo(x,bF))                 {printf("[0] Fail on prop3d.%f\n",x); exit(1);}
  if (!trpO2->PropagateParamTo(x,bF))            {printf("[1] Fail on prop3d.%f\n",x); exit(1);}
  //
  printf("\n\nBackToVertex with 3D field, w/o updates:\n");
  printf("[0] TrackParCov:\n");
  trcO2->Print();
  printf("[1] TrackPar:\n");
  trpO2->PrintParam();
  //
  printf("\nDifference between full and param-only back-propagation\n");
  CompPar(trpO2.get(),trcO2.get());
  //
  printf("Sizes of TrackPar: %lu TrackParCov: %lu bytes\n",sizeof(TrackPar),sizeof(TrackParCov));
  //
  printf("doing output test>>\n");
  unique_ptr<Track>  trcIO(new Track(*trcO2.get())); 
  TFile* fout = TFile::Open("trcOut.root","recreate");
  trcIO.get()->Write("trcIO");
  fout->Close();
  delete fout;
  printf("Output test done>>\n");  
  //
  TFile* finp = TFile::Open("trcOut.root");
  Track* trcInp = (Track*)finp->Get("trcIO");
  finp->Close();
  delete finp;
  trcInp->Print();
  printf("Input test done>>\n");  
  //
}

void CompPar(const TrackParBase* t0,const TrackParBase* t1)
{
  printf("Diff: ");
  for (int i=0;i<5;i++) printf("%+.3e ",t0->GetParam()[i]-t1->GetParam()[i]);
  printf("\n");
}
