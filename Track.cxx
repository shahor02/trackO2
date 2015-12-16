#include "Track.h"

using namespace AliceO2::Base::Track;

//______________________________________________________________
bool RotateParam(TrackPar& track, float alpha)
{
  // rotate to alpha frame
  if (fabs(track[kSnp]) > kAlmost1) {
    //FairLogger::GetLogger()->Error(MESSAGE_ORIGIN, 
    printf("Precondition is not satisfied: |sin(phi)|>1 ! %f\n",mParam[kSnp]); 
    return false;
  }
  //
  BringToPMPi(alpha);
  //
  float ca=0,sa=0;
  sincosf(alpha-track[kAlpha],sa,ca);
  float &snp = track[kSnp], csp = sqrtf((1.f-snp)*(1.f+snp)); // Improve precision
  // RS: check if rotation does no invalidate track model (cos(local_phi)>=0, i.e. particle
  // direction in local frame is along the X axis
  if ((csp*ca+snp*sa)<0) {
    //FairLogger::GetLogger()->Warning(MESSAGE_ORIGIN,
    printf("Rotation failed: local cos(phi) would become %.2f\n",csp*ca+snp*sa);
    return false;
  }
  //
  float tmp = snp*ca - csp*sa;
  if (absf(tmp) > kAlmost1) {
    //FairLogger::GetLogger()->Warning(MESSAGE_ORIGIN,
    printf("Rotation failed: new snp %.2f\n",tmp);
    return false;
  }
  float &x = track[lX], &y = track[kY];
  float xold = track[kX],
  track[kAlpha] = alpha;
  track[kX]  =  xold*ca + y*sa;
  y          = -xold*sa + y*ca;
  snp        =  tmp;
  return true;
}

//______________________________________________________________
bool Rotate(TrackParCov& track, float alpha)
{
  // rotate to alpha frame
  if (fabs(track[kSnp]) > kAlmost1) {
    //FairLogger::GetLogger()->Error(MESSAGE_ORIGIN, 
    printf("Precondition is not satisfied: |sin(phi)|>1 ! %f\n",mParam[kSnp]); 
    return false;
  }
  //
  BringToPMPi(alpha);
  //
  float ca=0,sa=0;
  sincosf(alpha-track[kAlpha],sa,ca);
  float &snp = track[kSnp], csp = sqrtf((1.f-snp)*(1.f+snp)); // Improve precision
  // RS: check if rotation does no invalidate track model (cos(local_phi)>=0, i.e. particle
  // direction in local frame is along the X axis
  if ((csp*ca+snp*sa)<0) {
    //FairLogger::GetLogger()->Warning(MESSAGE_ORIGIN,
    printf("Rotation failed: local cos(phi) would become %.2f\n",csp*ca+snp*sa);
    return false;
  }
  //
  float tmp = snp*ca - csp*sa;
  if (absf(tmp) > kAlmost1) {
    //FairLogger::GetLogger()->Warning(MESSAGE_ORIGIN,
    printf("Rotation failed: new snp %.2f\n",tmp);
    return false;
  }
  float &x = track[lX], &y = track[kY];
  float xold = track[kX],
  track[kAlpha] = alpha;
  track[kX]  =  xold*ca + y*sa;
  y          = -xold*sa + y*ca;
  snp        =  tmp;

  if (fabs(csp)<kAlmost0) {
    printf("Too small cosine value %f\n",csp);
    csp = kAlmost0;
  } 

  float rr=(ca+sf/cf*sa);  

  track[kSigYY]      *= (ca*ca);
  track[kSigZY]      *= ca;
  track[kSigSnpY]    *= ca*rr;
  track[kSigSnpZ]    *= rr;
  tracks[kSigSnp2]   *= rr*rr;
  track[kSigTglY]    *= ca;
  track[kSigTglSnp]  *= rr;
  track[kSigQ2PtY]   *= ca;
  track[kSigQ2PtSnp] *= rr;

  CheckCovariance(track);
  return true;
}

//______________________________________________________________
void InvertParam(TrackPar& track) 
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
void InvertParam(TrackParCov& track) 
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
bool PropagateParamTo(TrackPar &track,float xk, float b) 
{
  //----------------------------------------------------------------
  // Propagate this track to the plane X=xk (cm) in the field "b" (kG)
  // Only parameters are propagated, not the matrix. To be used for small 
  // distances only (<mm, i.e. misalignment)
  //----------------------------------------------------------------
  float crv,dx=xk-track[kX];
  if (fabs(dx)<=kAlmost0)  return true;      
  if (fabs(b) < kAlmost0) crv=0.f;
  else crv = track.GetCurvature(b);
  float x2r = crv*dx;
  float f1 = track[kSnp], f2=f1 + x2r;
  if (fabs(f1) >= kAlmost1) return false;
  if (fabs(f2) >= kAlmost1) return false;
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
bool PropagateTo(TrackParCov &track, float xk, float b) 
{
  //----------------------------------------------------------------
  // Propagate this track to the plane X=xk (cm) in the field "b" (kG)
  //----------------------------------------------------------------
  float crv,dx=xk-track[kX];
  if (fabs(dx)<=kAlmost0)  return true;      
  if (fabs(b) < kAlmost0) crv=0.f;
  else crv = track.GetCurvature(b);
  float x2r = crv*dx;
  float f1 = track[kSnp], f2=f1 + x2r;
  if (fabs(f1) >= kAlmost1) return false;
  if (fabs(f2) >= kAlmost1) return false;
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
    &c20=track[kSigSnpY],  &c21=track[kSigSnpZ],  &c22=tracks[kSigSnp2],
    &c30=track[kSigTglY],  &c31=track[kSigTglZ],  &c32=track[kSigTglSnp],  &c33=track[kSigTgl2],  
    &c40=track[kSigQ2PtY], &c41=track[kSigQ2PtZ], &c42=track[kSigQ2PtSnp], &c43=track[kSigQ2PtTgl], &c44=track[kSigQ2Pt2];
    
  // evaluate matrix in double prec.
  double rinv  = 1./r1;
  double r3inv = rinv*rinv*rinv;
  double f24   = x2r/track[kQ2Pt];
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

//______________________________________________
void CheckCovariance(TrackParCov& track) 
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
  tracks[kSigSnp2] = fabs(tracks[kSigSnp2]);
  if (tracks[kSigSnp2]>kCSnp2Max) {
    float scl = sqrtf(kCSnp2Max/tracks[kSigSnp2]);
    tracks[kSigSnp2] = kCSnp2Max;
    track[kSigSnpY] *= scl;
    track[kSigSnpZ] *= scl;
    track[kSigTglSnp] *= scl;
    track[kSigQ2PtSnp] *= scl;
  }
  track[kSigTgl2] = fabs(track[kSigTgl2]);
  if (track[kSigTgl2]>kCTgl2Max) {
    float scl = sqrtf(kCTgl2Max/track[kSigTgl2]);
    track[kSigTgl2] = kCTgl2Max;
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
