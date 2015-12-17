#include "Track.h"

using namespace AliceO2::Base;


//______________________________________________________________
Track::TrackParam::TrackParam(const float xyz[3],const float pxpypz[3],int sign, bool sectorAlpha)
{
  // construct track from kinematics

  // Alpha of the frame is defined as:
  // sectorAlpha == false : -> angle of pt direction
  // sectorAlpha == true : ->  angle of the sector from X,Y coordinate for r>1 
  //                           angle of pt direction for r==0
  //
  //
  const float kSafe = 1e-5f;
  float radPos2 = xyz[0]*xyz[0]+xyz[1]*xyz[1];  
  float alp = 0;
  if (sectorAlpha || radPos2<1) alp = atan2f(pxpypz[1],pxpypz[0]);
  else                          alp = atan2f(xyz[1],xyz[0]);
  if (sectorAlpha) alp = Angle2Alpha(alp);
  //
  float sn,cs; 
  sincosf(alp,sn,cs);
  // protection:  avoid alpha being too close to 0 or +-pi/2
  if (fabs(sn)<2*kSafe) {
    if (alp>0) alp += alp< kPIHalf ?  2*kSafe : -2*kSafe;
    else       alp += alp>-kPIHalf ? -2*kSafe :  2*kSafe;
    sincosf(alp,sn,cs);
  }
  else if (fabs(cs)<2*kSafe) {
    if (alp>0) alp += alp> kPIHalf ? 2*kSafe : -2*kSafe;
    else       alp += alp>-kPIHalf ? 2*kSafe : -2*kSafe;
    sincosf(alp,sn,cs);
  }
  // Get the vertex of origin and the momentum
  TVector3 ver(xyz[0],xyz[1],xyz[2]);
  TVector3 mom(pxpypz[0],pxpypz[1],pxpypz[2]);
  //
  // Rotate to the local coordinate system
  ver.RotateZ(-fAlpha);
  mom.RotateZ(-fAlpha);

  //
  // x of the reference plane
  fX = ver.X();

  Double_t charge = (Double_t)sign;

  fP[0] = ver.Y();
  fP[1] = ver.Z();
  fP[2] = TMath::Sin(mom.Phi());
  fP[3] = mom.Pz()/mom.Pt();
  fP[4] = TMath::Sign(1/mom.Pt(),charge);
  //
  if      (TMath::Abs( 1-fP[2]) < kSafe) fP[2] = 1.- kSafe; //Protection
  else if (TMath::Abs(-1-fP[2]) < kSafe) fP[2] =-1.+ kSafe; //Protection
  //
  // Covariance matrix (formulas to be simplified)
  Double_t pt=1./TMath::Abs(fP[4]);
  Double_t r=TMath::Sqrt((1.-fP[2])*(1.+fP[2]));
  //
  Double_t cv34 = TMath::Sqrt(cv[3 ]*cv[3 ]+cv[4 ]*cv[4 ]);
  //
  Int_t special = 0;
  double sgcheck = r*sn + fP[2]*cs;
  if (TMath::Abs(sgcheck)>=1-kSafe) { // special case: lab phi is +-pi/2
    special = 1;
    sgcheck = TMath::Sign(1.0,sgcheck);
  }
  else if (TMath::Abs(sgcheck)<kSafe) {
    sgcheck = TMath::Sign(1.0,cs);
    special = 2;   // special case: lab phi is 0
  }
  //
}

//______________________________________________________________
bool Track::RotateParam(TrackPar& track, float alpha)
{
  // rotate to alpha frame
  if (fabs(track[kSnp]) > kAlmost1) {
    //FairLogger::GetLogger()->Error(MESSAGE_ORIGIN, 
    printf("Precondition is not satisfied: |sin(phi)|>1 ! %f\n",track[kSnp]); 
    return false;
  }
  //
  BringToPMPi(alpha);
  //
  float ca=0,sa=0;
  sincosf(alpha-track[kAlpha],sa,ca);
  float snp = track[kSnp], csp = sqrtf((1.f-snp)*(1.f+snp)); // Improve precision
  // RS: check if rotation does no invalidate track model (cos(local_phi)>=0, i.e. particle
  // direction in local frame is along the X axis
  if ((csp*ca+snp*sa)<0) {
    //FairLogger::GetLogger()->Warning(MESSAGE_ORIGIN,
    printf("Rotation failed: local cos(phi) would become %.2f\n",csp*ca+snp*sa);
    return false;
  }
  //
  float tmp = snp*ca - csp*sa;
  if (fabs(tmp) > kAlmost1) {
    //FairLogger::GetLogger()->Warning(MESSAGE_ORIGIN,
    printf("Rotation failed: new snp %.2f\n",tmp);
    return false;
  }
  float xold = track[kX], yold = track[kY];
  track[kAlpha] = alpha;
  track[kX]  =  xold*ca + yold*sa;
  track[kY]  = -xold*sa + yold*ca;
  track[kSnp]=  tmp;
  return true;
}

//______________________________________________________________
bool Track::Rotate(TrackParCov& track, float alpha)
{
  // rotate to alpha frame
  if (fabs(track[kSnp]) > kAlmost1) {
    //FairLogger::GetLogger()->Error(MESSAGE_ORIGIN, 
    printf("Precondition is not satisfied: |sin(phi)|>1 ! %f\n",track[kSnp]); 
    return false;
  }
  //
  BringToPMPi(alpha);
  //
  float ca=0,sa=0;
  sincosf(alpha-track[kAlpha],sa,ca);
  float snp = track[kSnp], csp = sqrtf((1.f-snp)*(1.f+snp)); // Improve precision
  // RS: check if rotation does no invalidate track model (cos(local_phi)>=0, i.e. particle
  // direction in local frame is along the X axis
  if ((csp*ca+snp*sa)<0) {
    //FairLogger::GetLogger()->Warning(MESSAGE_ORIGIN,
    printf("Rotation failed: local cos(phi) would become %.2f\n",csp*ca+snp*sa);
    return false;
  }
  //
  float tmp = snp*ca - csp*sa;
  if (fabs(tmp) > kAlmost1) {
    //FairLogger::GetLogger()->Warning(MESSAGE_ORIGIN,
    printf("Rotation failed: new snp %.2f\n",tmp);
    return false;
  }
  float xold = track[kX], yold = track[kY];
  track[kAlpha] = alpha;
  track[kX]  =  xold*ca + yold*sa;
  track[kY]  = -xold*sa + yold*ca;
  track[kSnp]=  tmp;

  if (fabs(csp)<kAlmost0) {
    printf("Too small cosine value %f\n",csp);
    csp = kAlmost0;
  } 

  float rr=(ca+snp/csp*sa);  

  track[kSigY2]      *= (ca*ca);
  track[kSigZY]      *= ca;
  track[kSigSnpY]    *= ca*rr;
  track[kSigSnpZ]    *= rr;
  track[kSigSnp2]    *= rr*rr;
  track[kSigTglY]    *= ca;
  track[kSigTglSnp]  *= rr;
  track[kSigQ2PtY]   *= ca;
  track[kSigQ2PtSnp] *= rr;

  CheckCovariance(track);
  return true;
}

//______________________________________________________________
void Track::InvertParam(TrackPar& track) 
{
  // Transform this track to the local coord. system rotated by 180 deg. 
  track[kX] = -track[kX];
  track[kAlpha] += kPI;
  BringToPMPi(track[kAlpha]);
  //
  track[0] = -track[0];
  track[3] = -track[3];
  track[4] = -track[4];
  //
}

//______________________________________________________________
void Track::Invert(TrackParCov& track) 
{
  // Transform this track to the local coord. system rotated by 180 deg. 
  InvertParam(track);
  // since the fP1 and fP2 are not inverted, their covariances with others change sign
  track[kSigZY]      = -track[kSigZY];
  track[kSigSnpY]    = -track[kSigSnpY];
  track[kSigTglZ]    = -track[kSigTglZ];
  track[kSigTglSnp]  = -track[kSigTglSnp];
  track[kSigQ2PtZ]   = -track[kSigQ2PtZ];
  track[kSigQ2PtSnp] = -track[kSigQ2PtSnp];
}

//____________________________________________________________
bool Track::PropagateParamTo(TrackPar &track,float xk, float b) 
{
  //----------------------------------------------------------------
  // Propagate this track to the plane X=xk (cm) in the field "b" (kG)
  // Only parameters are propagated, not the matrix. To be used for small 
  // distances only (<mm, i.e. misalignment)
  //----------------------------------------------------------------
  float dx=xk-track[kX];
  if (fabs(dx)<kAlmost0)  return true;
  float crv = fabs(b) ? track.GetCurvature(b) : 0.f;
  if (fabs(b) < kAlmost0) crv=0.f;
  else crv = track.GetCurvature(b);
  float x2r = crv*dx;
  float f1 = track[kSnp], f2=f1 + x2r;
  if (fabs(f1) > kAlmost1) return false;
  if (fabs(f2) > kAlmost1) return false;
  if (fabs(track[kQ2Pt])< kAlmost0) return false;
  float r1=sqrtf((1.f-f1)*(1.f+f1)), r2=sqrtf((1.f-f2)*(1.f+f2));
  if (fabs(r1)<kAlmost0)  return false;
  if (fabs(r2)<kAlmost0)  return false;
  track[kX] = xk;
  double dy2dx = (f1+f2)/(r1+r2);
  track[kY] += dx*dy2dx;
  track[kSnp] += x2r;
  if (fabs(x2r)<0.05f) track[kZ] += dx*(r2 + f2*dy2dx)*track[kTgl];
  else { 
    // for small dx/R the linear apporximation of the arc by the segment is OK,
    // but at large dx/R the error is very large and leads to incorrect Z propagation
    // angle traversed delta = 2*asin(dist_start_end / R / 2), hence the arc is: R*deltaPhi
    // The dist_start_end is obtained from sqrt(dx^2+dy^2) = x/(r1+r2)*sqrt(2+f1*f2+r1*r2)
    //    double chord = dx*TMath::Sqrt(1+dy2dx*dy2dx);   // distance from old position to new one
    //    double rot = 2*TMath::ASin(0.5*chord*crv); // angular difference seen from the circle center
    //    track1 += rot/crv*track3;
    // 
    float rot = asinf(r1*f2 - r2*f1); // more economic version from Yura.
    if (f1*f1+f2*f2>1.f && f1*f2<0.f) {          // special cases of large rotations or large abs angles
      if (f2>0.f) rot = kPI - rot;    //
      else       rot = -kPI - rot;
    }
    track[kZ] += track[kTgl]/crv*rot; 
  }
  return true;
}

//______________________________________________________________
bool Track::PropagateTo(TrackParCov &track, float xk, float b) 
{
  //----------------------------------------------------------------
  // Propagate this track to the plane X=xk (cm) in the field "b" (kG)
  //----------------------------------------------------------------
  float crv,dx=xk-track[kX];
  if (fabs(dx)< kAlmost0)  return true;      
  if (fabs(b) < kAlmost0) crv=0.f;
  else crv = track.GetCurvature(b);
  float x2r = crv*dx;
  float f1 = track[kSnp], f2=f1 + x2r;
  if (fabs(f1) > kAlmost1) return false;
  if (fabs(f2) > kAlmost1) return false;
  if (fabs(track[kQ2Pt])< kAlmost0) return false;
  float r1=sqrtf((1.f-f1)*(1.f+f1)), r2=sqrtf((1.f-f2)*(1.f+f2));
  if (fabs(r1)<kAlmost0)  return false;
  if (fabs(r2)<kAlmost0)  return false;
  track[kX] = xk;
  double dy2dx = (f1+f2)/(r1+r2);
  track[kY] += dx*dy2dx;
  track[kSnp] += x2r;
  if (fabs(x2r)<0.05f) track[kZ] += dx*(r2 + f2*dy2dx)*track[kTgl];
  else { 
    // for small dx/R the linear apporximation of the arc by the segment is OK,
    // but at large dx/R the error is very large and leads to incorrect Z propagation
    // angle traversed delta = 2*asin(dist_start_end / R / 2), hence the arc is: R*deltaPhi
    // The dist_start_end is obtained from sqrt(dx^2+dy^2) = x/(r1+r2)*sqrt(2+f1*f2+r1*r2)
    //    double chord = dx*TMath::Sqrt(1+dy2dx*dy2dx);   // distance from old position to new one
    //    double rot = 2*TMath::ASin(0.5*chord*crv); // angular difference seen from the circle center
    //    track1 += rot/crv*track3;
    // 
    float rot = asinf(r1*f2 - r2*f1); // more economic version from Yura.
    if (f1*f1+f2*f2>1.f && f1*f2<0.f) {          // special cases of large rotations or large abs angles
      if (f2>0.f) rot = kPI - rot;    //
      else       rot = -kPI - rot;
    }
    track[kZ] += track[kTgl]/crv*rot; 
  }
    
  float
    &c00=track[kSigY2],
    &c10=track[kSigZY],    &c11=track[kSigZ2],
    &c20=track[kSigSnpY],  &c21=track[kSigSnpZ],  &c22=track[kSigSnp2],
    &c30=track[kSigTglY],  &c31=track[kSigTglZ],  &c32=track[kSigTglSnp],  &c33=track[kSigTgl2],  
    &c40=track[kSigQ2PtY], &c41=track[kSigQ2PtZ], &c42=track[kSigQ2PtSnp], &c43=track[kSigQ2PtTgl], &c44=track[kSigQ2Pt2];
    
  // evaluate matrix in double prec.
  double rinv  = 1./r1;
  double r3inv = rinv*rinv*rinv;
  double f24   = dx*b*kB2C; // x2r/track[kQ2Pt];
  double f02   = dx*r3inv;
  double f04   = 0.5*f24*f02;
  double f12   = f02*track[kTgl]*f1;
  double f14   = 0.5*f24*f02*track[kTgl]*f1;
  double f13   = dx*rinv;
    
  //b = C*ft
  double b00=f02*c20 + f04*c40, b01=f12*c20 + f14*c40 + f13*c30;
  double b02=f24*c40;
  double b10=f02*c21 + f04*c41, b11=f12*c21 + f14*c41 + f13*c31;
  double b12=f24*c41;
  double b20=f02*c22 + f04*c42, b21=f12*c22 + f14*c42 + f13*c32;
  double b22=f24*c42;
  double b40=f02*c42 + f04*c44, b41=f12*c42 + f14*c44 + f13*c43;
  double b42=f24*c44;
  double b30=f02*c32 + f04*c43, b31=f12*c32 + f14*c43 + f13*c33;
  double b32=f24*c43;
  
  //a = f*b = f*C*ft
  double a00=f02*b20+f04*b40,a01=f02*b21+f04*b41,a02=f02*b22+f04*b42;
  double a11=f12*b21+f14*b41+f13*b31,a12=f12*b22+f14*b42+f13*b32;
  double a22=f24*b42;
    
  //F*C*Ft = C + (b + bt + a)
  c00 += b00 + b00 + a00;
  c10 += b10 + b01 + a01; 
  c20 += b20 + b02 + a02;
  c30 += b30;
  c40 += b40;
  c11 += b11 + b11 + a11;
  c21 += b21 + b12 + a12;
  c31 += b31; 
  c41 += b41;
  c22 += b22 + b22 + a22;
  c32 += b32;
  c42 += b42;
    
  CheckCovariance(track);
    
  return true;
}

//____________________________________________________________
bool Track::PropagateParamBxByBzTo(TrackPar& track, float xk, const float b[3])
{
  //----------------------------------------------------------------
  // Extrapolate this track params (w/o cov matrix) to the plane X=xk in the field b[].
  //
  // X [cm] is in the "tracking coordinate system" of this track.
  // b[]={Bx,By,Bz} [kG] is in the Global coordidate system.
  //----------------------------------------------------------------

  float dx=xk-track.GetX();
  if (fabs(dx)<kAlmost0)  return true;
  if (fabs(track.GetQ2Pt())<kAlmost0) return true;
  // Do not propagate tracks outside the ALICE detector
  if (fabs(dx)>1e5 || fabs(track.GetY())>1e5 || fabs(track.GetZ())>1e5) {
    printf("Anomalous track, target X:%f\n",xk);
    //    Print();
    return false;
  }
  float crv = fabs(b[2]) ? track.GetCurvature(b[2]) : 0.f;
  float x2r = crv*dx;
  float f1 = track[kSnp], f2 = f1 + x2r;
  if (fabs(f1)>kAlmost1 || fabs(f2)>kAlmost1) return false;
  if (fabs(track[kQ2Pt])<kAlmost0) return false;
  float r1=sqrtf((1.f-f1)*(1.f+f1)), r2=sqrtf((1.f-f2)*(1.f+f2));
  if (fabs(r1)<kAlmost0 || fabs(r2)<kAlmost0)  return false;
  float dy2dx = (f1+f2)/(r1+r2);
  float step = (fabs(x2r)<0.05f) ? dx*fabs(r2 + f2*dy2dx)      // chord
    : 2.f*asinf(0.5f*dx*sqrtf(1.f+dy2dx*dy2dx)*crv)/crv;       // arc
  step *= sqrtf(1.f+ track.GetTgl()*track.GetTgl());
  //
  // Get the track x,y,z,px/p,py/p,pz/p,p,sinAlpha,cosAlpha in the Global System
  float vecLab[9];
  if (!GetPosDir(track, vecLab)) return false;

  // Rotate to the system where Bx=By=0.
  float bxy2 = b[0]*b[0] + b[1]*b[1];
  float bt = sqrtf(bxy2);
  float cosphi=1.f, sinphi=0.f;
  if (bt > kAlmost0) {
    cosphi=b[0]/bt; 
    sinphi=b[1]/bt;
  }
  float bb = sqrtf(bxy2 + b[2]*b[2]);
  float costet=1., sintet=0.;
  if (bb > kAlmost0) {
    costet=b[2]/bb; 
    sintet=bt/bb;
  }
  float vect[7] = {
    costet*cosphi*vecLab[0] + costet*sinphi*vecLab[1] - sintet*vecLab[2],
    -sinphi*vecLab[0] + cosphi*vecLab[1],
    sintet*cosphi*vecLab[0] + sintet*sinphi*vecLab[1] + costet*vecLab[2],
    costet*cosphi*vecLab[3] + costet*sinphi*vecLab[4] - sintet*vecLab[5],
    -sinphi*vecLab[3] + cosphi*vecLab[4],
    sintet*cosphi*vecLab[3] + sintet*sinphi*vecLab[4] + costet*vecLab[5],
    vecLab[6]};

  // Do the helix step
  float sgn = track.GetSign();
  g3helx3(sgn*bb,step,vect);

  // Rotate back to the Global System
  vecLab[0] = cosphi*costet*vect[0] - sinphi*vect[1] + cosphi*sintet*vect[2];
  vecLab[1] = sinphi*costet*vect[0] + cosphi*vect[1] + sinphi*sintet*vect[2];
  vecLab[2] = -sintet*vect[0] + costet*vect[2];

  vecLab[3] = cosphi*costet*vect[3] - sinphi*vect[4] + cosphi*sintet*vect[5];
  vecLab[4] = sinphi*costet*vect[3] + cosphi*vect[4] + sinphi*sintet*vect[5];
  vecLab[5] = -sintet*vect[3] + costet*vect[5];

  // Rotate back to the Tracking System
  float sinalp=vecLab[7],cosalp=-vecLab[8];
  float t = cosalp*vecLab[0] - sinalp*vecLab[1];
  vecLab[1] = sinalp*vecLab[0] + cosalp*vecLab[1];  
  vecLab[0] = t;
  t    = cosalp*vecLab[3] - sinalp*vecLab[4]; 
  vecLab[4] = sinalp*vecLab[3] + cosalp*vecLab[4];
  vecLab[3] = t; 

  // Do the final correcting step to the target plane (linear approximation)
  float x=vecLab[0], y=vecLab[1], z=vecLab[2];
  if (fabs(dx) > kAlmost0) {
    if (fabs(vecLab[3]) < kAlmost0) return false;
    dx = xk - vecLab[0];
    x += dx;
    y += vecLab[4]/vecLab[3]*dx;
    z += vecLab[5]/vecLab[3]*dx;  
  }
  
  // Calculate the track parameters
  t = 1.f/sqrtf(vecLab[3]*vecLab[3] + vecLab[4]*vecLab[4]);
  track[kX]    = x;
  track[kY]    = y;
  track[kZ]    = z;
  track[kSnp]  = vecLab[4]*t;
  track[kTgl]  = vecLab[5]*t; 
  track[kQ2Pt] = sgn*t/vecLab[6];

  return true;
}

//____________________________________________________________
bool Track::PropagateBxByBzTo(TrackParCov& track, float xk, const float b[3])
{
  //----------------------------------------------------------------
  // Extrapolate this track to the plane X=xk in the field b[].
  //
  // X [cm] is in the "tracking coordinate system" of this track.
  // b[]={Bx,By,Bz} [kG] is in the Global coordidate system.
  //----------------------------------------------------------------

  float dx=xk-track.GetX();
  if (fabs(dx)<kAlmost0)  return true;
  if (fabs(track.GetQ2Pt())<kAlmost0) return true;
  // Do not propagate tracks outside the ALICE detector
  if (fabs(dx)>1e5 || fabs(track.GetY())>1e5 || fabs(track.GetZ())>1e5) {
    printf("Anomalous track, target X:%f\n",xk);
    //    Print();
    return false;
  }
  float crv = fabs(b[2]) ? track.GetCurvature(b[2]) : 0.f;
  float x2r = crv*dx;
  float f1 = track[kSnp], f2 = f1 + x2r;
  if (fabs(f1)>kAlmost1 || fabs(f2)>kAlmost1) return false;
  if (fabs(track[kQ2Pt])<kAlmost0) return false;
  float r1=sqrtf((1.f-f1)*(1.f+f1)), r2=sqrtf((1.f-f2)*(1.f+f2));
  if (fabs(r1)<kAlmost0 || fabs(r2)<kAlmost0)  return false;
  float dy2dx = (f1+f2)/(r1+r2);
  float step = (fabs(x2r)<0.05f) ? dx*fabs(r2 + f2*dy2dx)      // chord
    : 2.f*asinf(0.5f*dx*sqrtf(1.f+dy2dx*dy2dx)*crv)/crv;       // arc
  step *= sqrtf(1.f+ track.GetTgl()*track.GetTgl());
  //
  // Get the track x,y,z,px/p,py/p,pz/p,p,sinAlpha,cosAlpha in the Global System
  float vecLab[9];
  if (!GetPosDir(track, vecLab)) return false;
  //
  // matrix transformed with Bz component only
  float
    &c00=track[kSigY2],
    &c10=track[kSigZY],    &c11=track[kSigZ2],
    &c20=track[kSigSnpY],  &c21=track[kSigSnpZ],  &c22=track[kSigSnp2],
    &c30=track[kSigTglY],  &c31=track[kSigTglZ],  &c32=track[kSigTglSnp],  &c33=track[kSigTgl2],  
    &c40=track[kSigQ2PtY], &c41=track[kSigQ2PtZ], &c42=track[kSigQ2PtSnp], &c43=track[kSigQ2PtTgl], &c44=track[kSigQ2Pt2];
  // evaluate matrix in double prec.
  double rinv  = 1./r1;
  double r3inv = rinv*rinv*rinv;
  double f24   = dx*b[2]*kB2C; // x2r/track[kQ2Pt];
  double f02   = dx*r3inv;
  double f04   = 0.5*f24*f02;
  double f12   = f02*track[kTgl]*f1;
  double f14   = 0.5*f24*f02*track[kTgl]*f1;
  double f13   = dx*rinv;
    
  //b = C*ft
  double b00=f02*c20 + f04*c40, b01=f12*c20 + f14*c40 + f13*c30;
  double b02=f24*c40;
  double b10=f02*c21 + f04*c41, b11=f12*c21 + f14*c41 + f13*c31;
  double b12=f24*c41;
  double b20=f02*c22 + f04*c42, b21=f12*c22 + f14*c42 + f13*c32;
  double b22=f24*c42;
  double b40=f02*c42 + f04*c44, b41=f12*c42 + f14*c44 + f13*c43;
  double b42=f24*c44;
  double b30=f02*c32 + f04*c43, b31=f12*c32 + f14*c43 + f13*c33;
  double b32=f24*c43;
  
  //a = f*b = f*C*ft
  double a00=f02*b20+f04*b40,a01=f02*b21+f04*b41,a02=f02*b22+f04*b42;
  double a11=f12*b21+f14*b41+f13*b31,a12=f12*b22+f14*b42+f13*b32;
  double a22=f24*b42;
    
  //F*C*Ft = C + (b + bt + a)
  c00 += b00 + b00 + a00;
  c10 += b10 + b01 + a01; 
  c20 += b20 + b02 + a02;
  c30 += b30;
  c40 += b40;
  c11 += b11 + b11 + a11;
  c21 += b21 + b12 + a12;
  c31 += b31; 
  c41 += b41;
  c22 += b22 + b22 + a22;
  c32 += b32;
  c42 += b42;
    
  CheckCovariance(track);
  //
  // Rotate to the system where Bx=By=0.
  float bxy2 = b[0]*b[0] + b[1]*b[1];
  float bt = sqrtf(bxy2);
  float cosphi=1.f, sinphi=0.f;
  if (bt > kAlmost0) {
    cosphi=b[0]/bt; 
    sinphi=b[1]/bt;
  }
  float bb = sqrtf(bxy2 + b[2]*b[2]);
  float costet=1., sintet=0.;
  if (bb > kAlmost0) {
    costet=b[2]/bb; 
    sintet=bt/bb;
  }
  float vect[7] = {
    costet*cosphi*vecLab[0] + costet*sinphi*vecLab[1] - sintet*vecLab[2],
    -sinphi*vecLab[0] + cosphi*vecLab[1],
    sintet*cosphi*vecLab[0] + sintet*sinphi*vecLab[1] + costet*vecLab[2],
    costet*cosphi*vecLab[3] + costet*sinphi*vecLab[4] - sintet*vecLab[5],
    -sinphi*vecLab[3] + cosphi*vecLab[4],
    sintet*cosphi*vecLab[3] + sintet*sinphi*vecLab[4] + costet*vecLab[5],
    vecLab[6]};

  // Do the helix step
  float sgn = track.GetSign();
  g3helx3(sgn*bb,step,vect);

  // Rotate back to the Global System
  vecLab[0] = cosphi*costet*vect[0] - sinphi*vect[1] + cosphi*sintet*vect[2];
  vecLab[1] = sinphi*costet*vect[0] + cosphi*vect[1] + sinphi*sintet*vect[2];
  vecLab[2] = -sintet*vect[0] + costet*vect[2];

  vecLab[3] = cosphi*costet*vect[3] - sinphi*vect[4] + cosphi*sintet*vect[5];
  vecLab[4] = sinphi*costet*vect[3] + cosphi*vect[4] + sinphi*sintet*vect[5];
  vecLab[5] = -sintet*vect[3] + costet*vect[5];

  // Rotate back to the Tracking System
  float sinalp=vecLab[7],cosalp=-vecLab[8];
  float t = cosalp*vecLab[0] - sinalp*vecLab[1];
  vecLab[1] = sinalp*vecLab[0] + cosalp*vecLab[1];  
  vecLab[0] = t;
  t    = cosalp*vecLab[3] - sinalp*vecLab[4]; 
  vecLab[4] = sinalp*vecLab[3] + cosalp*vecLab[4];
  vecLab[3] = t; 

  // Do the final correcting step to the target plane (linear approximation)
  float x=vecLab[0], y=vecLab[1], z=vecLab[2];
  if (fabs(dx) > kAlmost0) {
    if (fabs(vecLab[3]) < kAlmost0) return false;
    dx = xk - vecLab[0];
    x += dx;
    y += vecLab[4]/vecLab[3]*dx;
    z += vecLab[5]/vecLab[3]*dx;  
  }
  
  // Calculate the track parameters
  t = 1.f/sqrtf(vecLab[3]*vecLab[3] + vecLab[4]*vecLab[4]);
  track[kX]    = x;
  track[kY]    = y;
  track[kZ]    = z;
  track[kSnp]  = vecLab[4]*t;
  track[kTgl]  = vecLab[5]*t; 
  track[kQ2Pt] = sgn*t/vecLab[6];

  return true;

}


//______________________________________________
void Track::CheckCovariance(TrackParCov& track) 
{
  // This function forces the diagonal elements of the covariance matrix to be positive.
  // In case the diagonal element is bigger than the maximal allowed value, it is set to
  // the limit and the off-diagonal elements that correspond to it are set to zero.
  
  track[kSigY2] = fabs(track[kSigY2]);
  if (track[kSigY2]>kCY2max) {
    float scl = sqrtf(kCY2max/track[kSigY2]);
    track[kSigY2]     = kCY2max;
    track[kSigZY]    *= scl;
    track[kSigSnpY]  *= scl;
    track[kSigTglY]  *= scl;
    track[kSigQ2PtY] *= scl;
  }
  track[kSigZ2] = fabs(track[kSigZ2]);
  if (track[kSigZ2]>kCZ2max) {
    float scl = sqrtf(kCZ2max/track[kSigZ2]);
    track[kSigZ2]     = kCZ2max;
    track[kSigZY]    *= scl;
    track[kSigSnpZ]  *= scl;
    track[kSigTglZ]  *= scl;
    track[kSigQ2PtZ] *= scl;
  }
  track[kSigSnp2] = fabs(track[kSigSnp2]);
  if (track[kSigSnp2]>kCSnp2max) {
    float scl = sqrtf(kCSnp2max/track[kSigSnp2]);
    track[kSigSnp2] = kCSnp2max;
    track[kSigSnpY] *= scl;
    track[kSigSnpZ] *= scl;
    track[kSigTglSnp] *= scl;
    track[kSigQ2PtSnp] *= scl;
  }
  track[kSigTgl2] = fabs(track[kSigTgl2]);
  if (track[kSigTgl2]>kCTgl2max) {
    float scl = sqrtf(kCTgl2max/track[kSigTgl2]);
    track[kSigTgl2] = kCTgl2max;
    track[kSigTglY] *= scl;
    track[kSigTglZ] *= scl;
    track[kSigTglSnp] *= scl;
    track[kSigQ2PtTgl] *= scl;
  }
  track[kSigQ2Pt2] = fabs(track[kSigQ2Pt2]);
  if (track[kSigQ2Pt2]>kCQ2Pt2max) {
    float scl = sqrtf(kCQ2Pt2max/track[kSigQ2Pt2]);
    track[kSigQ2Pt2] = kCQ2Pt2max;
    track[kSigQ2PtY] *= scl;
    track[kSigQ2PtZ] *= scl;
    track[kSigQ2PtSnp] *= scl;
    track[kSigQ2PtTgl] *= scl;
  }
}


//=================================================
//
// Aux. methods for tracks manipulation
//
//=================================================


//____________________________________________________
bool Track::TrackPar2Momentum(float p[3], float alpha) 
{
  //----------------------------------------------------------------
  // This function performs local->global transformation of the
  // track momentum.
  // When called, the arguments are:
  //    p[0] = 1/pt * charge of the track;
  //    p[1] = sine of local azim. angle of the track momentum;
  //    p[2] = tangent of the track momentum dip angle;
  //   alpha - rotation angle. 
  // The result is returned as:
  //    p[0] = px
  //    p[1] = py
  //    p[2] = pz
  // Results for (nearly) straight tracks are meaningless !
  //----------------------------------------------------------------
  if (fabs(p[0])<kAlmost0 || fabs(p[1])>kAlmost1) return false;
  float cs,sn, pt=1.f/fabs(p[0]);
  float r = sqrtf((1.f - p[1])*(1.f + p[1]));
  sincosf(alpha,sn,cs);
  p[0] = pt*(r*cs - p[1]*sn); 
  p[1] = pt*(p[1]*cs + r*sn); 
  p[2] = pt*p[2];
  return true;
}

//____________________________________________________
bool Track::GetPosDir(const TrackPar& track, float posdirp[9])
{
  // fill vector with lab x,y,z,px/p,py/p,pz/p,p,sinAlpha,cosAlpha
  float pti = fabs(track.GetQ2Pt());
  float snp = track.GetSnp();
  if (pti<kAlmost0 || fabs(snp)>kAlmost1) return false;
  float sn=posdirp[7],cs=posdirp[8]; 
  float csp = sqrtf((1.f - snp)*(1.f + snp));
  float cstht = sqrtf(1.f+ track.GetTgl()*track.GetTgl());
  float csthti = 1.f/cstht;
  float x = track.GetX();
  sincosf(track.GetAlpha(),sn,cs);
  posdirp[0] = track.GetX()*cs - track.GetY()*sn;
  posdirp[1] = track.GetX()*sn + track.GetY()*cs;
  posdirp[2] = track.GetZ();
  posdirp[3] = (csp*cs - snp*sn)*csthti;  // px/p
  posdirp[4] = (snp*cs + csp*sn)*csthti;  // py/p
  posdirp[5] = track.GetTgl()*csthti;     // pz/p
  posdirp[6] = cstht/pti;                 // p
  return true;
}

