#include "AliExternalTrackParam.h"
#include "Track.h"

using namespace AliceO2::Base::Track;

void CompPar(const AliExternalTrackParam* etpar,const TrackPar* tp);
void CompCov(const AliExternalTrackParam* etpar,const TrackParCov* tp);

void test()
{
  double xyzD[3] = {0,0,0}, pxyzD[3] = {0.5,0.5,0.5};
  const float  xyzF[3] = {0,0,0}, pxyzF[3] = {0.5,0.5,0.5};
  double covD[21] = {1e-6,
		     0,1e-6,
		     0,0,1e-6,
		     0,0,0,1e-4,
		     0,0,0,0,1e-4,
		     0,0,0,0,0,1e-4};
  const float covF[21]  = {1e-6,
			   0,1e-6,
			   0,0,1e-6,
			   0,0,0,1e-4,
			   0,0,0,0,1e-4,
			   0,0,0,0,0,1e-4};

  TrackParCov *trcO2 = new TrackParCov(xyzF,pxyzF,covF,1,false);
  TrackPar *trpO2 = new TrackPar(xyzF,pxyzF,1,false);
  AliExternalTrackParam* etp = new AliExternalTrackParam(xyzD,pxyzD,covD,1);
  AliExternalTrackParam* etpp = new AliExternalTrackParam(xyzD,pxyzD,covD,1);
  //
  printf("\nAt creation:\n");
  printf("[0] TrackParCov:\n");
  trcO2->Print();
  printf("[1] TrackPar:\n");
  trpO2->Print();
  printf("[2] AliExternalTrackParam:\n");
  etp->Print();
  printf("[3] AliExternalTrackParam: param only\n");
  etpp->Print();
  //
  const float xmax = 400.,xstep=4.;;
  float bF[3] = {0.01,0.01,5.0};
  double bD[3] = {0.01,0.01,5.0};
  const float x2x0 = 0.01;
  const float xrho = xstep*2.0;
  const float mass = 0.14;
  //
  float x = trcO2->GetX();
  int cnt = 0;
  float chi0=0,chi3=0;
  float measErrF[3]={1e-6,0.25e-5,1e-4};
  double measErrD[3]={1e-6,0.25e-5,1e-4};
  int nupd = 0;
  while (x<xmax) {
    x += xstep;
    if (!trcO2->PropagateTo(x,bF[2]))           {printf("[0] Fail on prop.%f\n",x); exit(1);}
    if (!trpO2->PropagateParamTo(x,bF[2]))      {printf("[1] Fail on prop.%f\n",x); exit(1);}
    if (!etp->PropagateTo(x,bD[2]))             {printf("[2] Fail on prop.%f\n",x); exit(1);}
    if (!etpp->PropagateParamOnlyTo(x,bD[2]))   {printf("[3] Fail on prop.%f\n",x); exit(1);}

    cnt++;

    if (cnt%2==0) {
      if (!trcO2->CorrectForMaterial(x2x0,xrho,mass,false))  {printf("[0] Fail on matcor. %f\n",x); exit(1);}
      if (!etp->CorrectForMeanMaterial(x2x0,xrho,mass,false))    {printf("[2] Fail on matcor. %f\n",x); exit(1);}
    }

    if (cnt%3==0) {
      double measD[2] = {etp->GetY()+TMath::Sqrt(measErrD[0]),etp->GetZ()+TMath::Sqrt(measErrD[2])};
      float  measF[2] = {etp->GetY()+TMath::Sqrt(measErrF[0]),etp->GetZ()+TMath::Sqrt(measErrF[2])};
      chi0 += trcO2->GetPredictedChi2(measF,measErrF);
      chi3 += etp->GetPredictedChi2(measD,measErrD);
      if (!trcO2->Update(measF,measErrF))  {printf("[0] Fail on update. %f\n",x); exit(1);}
      if (!etp->Update(measD, measErrD))  {printf("[2] Fail on update. %f\n",x); exit(1);}
      nupd++;
    }
    //

    if (cnt%4 == 0 || TMath::Abs(x-xmax)<1e-4) {
      float alp = etp->Phi();
      if (!trcO2->Rotate(alp))       {printf("[0] Fail on rot.%f\n",x); exit(1);}
      if (!trpO2->RotateParam(alp))  {printf("[1] Fail on rot.%f\n",x); exit(1);}
      if (!etp->Rotate(alp))              {printf("[2] Fail on rot.%f\n",x); exit(1);}
      if (!etpp->RotateParamOnly(alp))    {printf("[3] Fail on rot.%f\n",x); exit(1);}
    }

  }

  printf("\nAt Max X:\n");
  printf("[0] TrackParCov:\n");
  trcO2->Print();
  printf("[1] TrackPar:\n");
  trpO2->Print();
  printf("[2] AliExternalTrackParam:\n");
  etp->Print();
  printf("[3] AliExternalTrackParam: param only\n");
  etpp->Print();

  printf("\nDifference for full propagation\n");
  CompPar(etp,(TrackPar*)trcO2);
  //CompPar(etp,trcO2);
  CompCov(etp,trcO2);

  printf("\nDifference for ParamOnly propagation\n");
  CompPar(etpp,trpO2);


  // propagate back with 3d field
  while (x>0) {
    x -= xstep;
    if (!trcO2->PropagateTo(x,bF))                 {printf("[0] Fail on prop3d.%f\n",x); exit(1);}
    if (!trpO2->PropagateParamTo(x,bF))            {printf("[1] Fail on prop3d.%f\n",x); exit(1);}
    if (!etp->PropagateToBxByBz(x,bD))             {printf("[2] Fail on prop3d.%f\n",x); exit(1);}
    if (!etpp->PropagateParamOnlyBxByBzTo(x,bD))   {printf("[3] Fail on prop3d.%f\n",x); exit(1);}
  }

  // bring to the same frame
  float alp = etp->Phi();
  if (!trcO2->Rotate(alp))       {printf("[0] Fail on rot.%f\n",x); exit(1);}
  if (!trpO2->RotateParam(alp))  {printf("[1] Fail on rot.%f\n",x); exit(1);}
  if (!etp->Rotate(alp))              {printf("[2] Fail on rot.%f\n",x); exit(1);}
  if (!etpp->RotateParamOnly(alp))    {printf("[3] Fail on rot.%f\n",x); exit(1);}
  //
  x = 0;
  if (!trcO2->PropagateTo(x,bF))                 {printf("[0] Fail on prop3d.%f\n",x); exit(1);}
  if (!trpO2->PropagateParamTo(x,bF))            {printf("[1] Fail on prop3d.%f\n",x); exit(1);}
  if (!etp->PropagateToBxByBz(x,bD))             {printf("[2] Fail on prop3d.%f\n",x); exit(1);}
  if (!etpp->PropagateParamOnlyBxByBzTo(x,bD))   {printf("[3] Fail on prop3d.%f\n",x); exit(1);}
  //
  printf("\n\nBackToVertex with 3D field, w/o updates:\n");
  printf("[0] TrackParCov:\n");
  trcO2->Print();
  printf("[1] TrackPar:\n");
  trpO2->Print();
  printf("[2] AliExternalTrackParam:\n");
  etp->Print();
  printf("[3] AliExternalTrackParam: param only\n");
  etpp->Print();
  if (nupd) {
    printf("[0] chi2/nupdate = %f\n",chi0/nupd);
    printf("[3] chi2/nupdate = %f\n",chi3/nupd);
  }
  //
  printf("\nDifference for full propagation\n");
  CompPar(etp,(TrackPar*)trcO2);
  //CompPar(etp,trcO2);
  CompCov(etp,trcO2);

  printf("\nDifference for ParamOnly propagation\n");
  CompPar(etpp,trpO2);

}

void CompPar(const AliExternalTrackParam* etpar,const TrackPar* tp)
{
  printf("DiffParam: dX: %+.3e dAlp: %+.3e | DPar: ",etpar->GetX()-tp->GetX(),etpar->GetAlpha()-tp->GetAlpha());
  for (int i=0;i<5;i++) printf("%+.3e ",etpar->GetParameter()[i]-(*tp)[i+2]);
  printf("\n");
  float xyzF[3];
  double xyzD[3];
  tp->GetXYZ(xyzF);
  etpar->GetXYZ(xyzD);
  printf("Diff XYZ: ");
  for (int i=0;i<3;i++) printf("%+e ",xyzD[i]-xyzF[i]); printf("\n");
}

void CompCov(const AliExternalTrackParam* etpar,const TrackParCov* tp)
{
  printf("DiffCov: \n");
  int cnt=0;
  for (int i=0;i<5;i++) {
    for (int j=0;j<=i;j++) {
      printf("%+.3e ",etpar->GetCovariance()[cnt]-(*tp)[7+cnt]);
      cnt++;
    }
    printf("\n");
  }
}
