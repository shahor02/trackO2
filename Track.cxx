#include "Track.h"

using namespace AliceO2::Base;


//______________________________________________________________
Track::TrackPar::TrackPar(const float xyz[3],const float pxpypz[3], int charge, bool sectorAlpha)
{
  // construct track param from kinematics

  // Alpha of the frame is defined as:
  // sectorAlpha == false : -> angle of pt direction
  // sectorAlpha == true  : -> angle of the sector from X,Y coordinate for r>1 
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
  float ver[3] = {xyz[0],xyz[1],xyz[2]};
  float mom[3] = {pxpypz[0],pxpypz[1],pxpypz[2]};
  //
  // Rotate to the local coordinate system
  RotateZ(ver,-alp);
  RotateZ(mom,-alp);
  //
  float ptI      = 1.f/sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
  mParam[kX]     = ver[0];
  mParam[kAlpha] = alp;
  mParam[kY]     = ver[1];
  mParam[kZ]     = ver[2];
  mParam[kSnp]   = mom[1]*ptI;
  mParam[kTgl]   = mom[2]*ptI;
  mParam[kQ2Pt]  = ptI*charge;
  //
  if      (fabs( 1-mParam[kSnp]) < kSafe) mParam[kSnp] = 1.- kSafe; //Protection
  else if (fabs(-1-mParam[kSnp]) < kSafe) mParam[kSnp] =-1.+ kSafe; //Protection
  //
}

//______________________________________________________________
Track::TrackParCov::TrackParCov(const float xyz[3],const float pxpypz[3], 
				const float cv[kLabCovMatSize], int charge, bool sectorAlpha)
{
  // construct track param and covariance from kinematics and lab errors

  // Alpha of the frame is defined as:
  // sectorAlpha == false : -> angle of pt direction
  // sectorAlpha == true  : -> angle of the sector from X,Y coordinate for r>1 
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
  float ver[3] = {xyz[0],xyz[1],xyz[2]};
  float mom[3] = {pxpypz[0],pxpypz[1],pxpypz[2]};
  //
  // Rotate to the local coordinate system
  RotateZ(ver,-alp);
  RotateZ(mom,-alp);
  //
  float pt       = sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
  float ptI      = 1.f/pt;
  mParCov[kX]     = ver[0];
  mParCov[kAlpha] = alp;
  mParCov[kY]     = ver[1];
  mParCov[kZ]     = ver[2];
  mParCov[kSnp]   = mom[1]*ptI; // cos(phi)
  mParCov[kTgl]   = mom[2]*ptI; // tg(lambda)
  mParCov[kQ2Pt]  = ptI*charge;
  //
  if      (fabs( 1-mParCov[kSnp]) < kSafe) mParCov[kSnp] = 1.- kSafe; //Protection
  else if (fabs(-1-mParCov[kSnp]) < kSafe) mParCov[kSnp] =-1.+ kSafe; //Protection
  //
  // Covariance matrix (formulas to be simplified)
  float r=mom[0]*ptI;  // cos(phi)
  float cv34 = sqrtf(cv[3]*cv[3]+cv[4]*cv[4]);
  //
  int special = 0;
  float sgcheck = r*sn + mParCov[kSnp]*cs;
  if (fabs(sgcheck)>1-kSafe) { // special case: lab phi is +-pi/2
    special = 1;
    sgcheck = sgcheck<0 ? -1.f:1.f;
  }
  else if (fabs(sgcheck)<kSafe) {
    sgcheck = cs<0 ? -1.0f:1.0f;
    special = 2;   // special case: lab phi is 0
  }
  //
  mParCov[kSigY2] = cv[0]+cv[2];  
  mParCov[kSigZY] = (-cv[3 ]*sn)<0 ? -cv34 : cv34;
  mParCov[kSigZ2] = cv[5]; 
  //
  float ptI2 = ptI*ptI;
  float tgl2 = mParCov[kTgl]*mParCov[kTgl];
  if (special==1) {
    mParCov[kSigSnpY   ] = cv[6]*ptI;
    mParCov[kSigSnpZ   ] = -sgcheck*cv[8]*r*ptI;
    mParCov[kSigSnp2   ] = fabs(cv[9]*r*r*ptI2);
    mParCov[kSigTglY   ] = (cv[10]*mParCov[kTgl]-sgcheck*cv[15])*ptI/r;
    mParCov[kSigTglZ   ] = (cv[17]-sgcheck*cv[12]*mParCov[kTgl])*ptI;
    mParCov[kSigTglSnp ] = (-sgcheck*cv[18]+cv[13]*mParCov[kTgl])*r*ptI2;
    mParCov[kSigTgl2   ] = fabs( cv[20]-2*sgcheck*cv[19]*mParCov[4]+cv[14]*tgl2)*ptI2;
    mParCov[kSigQ2PtY  ] = cv[10]*ptI2/r*charge;
    mParCov[kSigQ2PtZ  ] = -sgcheck*cv[12]*ptI2*charge;
    mParCov[kSigQ2PtSnp] = cv[13]*r*ptI*ptI2*charge;
    mParCov[kSigQ2PtTgl] = (-sgcheck*cv[19]+cv[14]*mParCov[kTgl])*r*ptI2*ptI;
    mParCov[kSigQ2Pt2  ] = fabs(cv[14]*ptI2*ptI2);
  } else if (special==2) {
    mParCov[kSigSnpY   ] = -cv[10]*ptI*cs/sn;
    mParCov[kSigSnpZ   ] = cv[12]*cs*ptI;
    mParCov[kSigSnp2   ] = fabs(cv[14]*cs*cs*ptI2);
    mParCov[kSigTglY   ] = (sgcheck*cv[6]*mParCov[kTgl]-cv[15])*ptI/sn;
    mParCov[kSigTglZ   ] = (cv[17]-sgcheck*cv[8]*mParCov[kTgl])*ptI;
    mParCov[kSigTglSnp ] = (cv[19]-sgcheck*cv[13]*mParCov[kTgl])*cs*ptI2;
    mParCov[kSigTgl2   ] = fabs( cv[20]-2*sgcheck*cv[18]*mParCov[kTgl]+cv[9]*tgl2)*ptI2;
    mParCov[kSigQ2PtY  ] = sgcheck*cv[6]*ptI2/sn*charge;
    mParCov[kSigQ2PtZ  ] = -sgcheck*cv[8]*ptI2*charge;
    mParCov[kSigQ2PtSnp] = -sgcheck*cv[13]*cs*ptI*ptI2*charge;
    mParCov[kSigQ2PtTgl] = (-sgcheck*cv[18]+cv[9]*mParCov[kTgl])*ptI2*ptI*charge;
    mParCov[kSigQ2Pt2  ] = fabs(cv[9]*ptI2*ptI2);
  }
  else {
    double m00=-sn;// m10=cs;
    double m23=-pt*(sn + mParCov[kSnp]*cs/r), m43=-pt*pt*(r*cs - mParCov[kSnp]*sn);
    double m24= pt*(cs - mParCov[kSnp]*sn/r), m44=-pt*pt*(r*sn + mParCov[kSnp]*cs);
    double m35=pt, m45=-pt*pt*mParCov[kTgl];
    //
    m43 *= charge;
    m44 *= charge;
    m45 *= charge;
    //
    double a1=cv[13]-cv[9]*(m23*m44+m43*m24)/m23/m43;
    double a2=m23*m24-m23*(m23*m44+m43*m24)/m43;
    double a3=m43*m44-m43*(m23*m44+m43*m24)/m23;
    double a4=cv[14]+2.*cv[9];
    double a5=m24*m24-2.*m24*m44*m23/m43;
    double a6=m44*m44-2.*m24*m44*m43/m23;
    //    
    mParCov[kSigSnpY ] = (cv[10]*m43-cv[6]*m44)/(m24*m43-m23*m44)/m00; 
    mParCov[kSigQ2PtY] = (cv[6]/m00-mParCov[kSigSnpY ]*m23)/m43; 
    mParCov[kSigTglY ] = (cv[15]/m00-mParCov[kSigQ2PtY]*m45)/m35; 
    mParCov[kSigSnpZ ] = (cv[12]*m43-cv[8]*m44)/(m24*m43-m23*m44); 
    mParCov[kSigQ2PtZ] = (cv[8]-mParCov[kSigSnpZ]*m23)/m43; 
    mParCov[kSigTglZ ] = cv[17]/m35-mParCov[kSigQ2PtZ]*m45/m35; 
    mParCov[kSigSnp2 ] = fabs((a4*a3-a6*a1)/(a5*a3-a6*a2));
    mParCov[kSigQ2Pt2] = fabs((a1-a2*mParCov[kSigSnp2])/a3);
    mParCov[kSigQ2PtSnp] = (cv[9]-mParCov[kSigSnp2]*m23*m23-mParCov[kSigQ2Pt2]*m43*m43)/m23/m43;
    double b1=cv[18]-mParCov[kSigQ2PtSnp]*m23*m45-mParCov[kSigQ2Pt2]*m43*m45;
    double b2=m23*m35;
    double b3=m43*m35;
    double b4=cv[19]-mParCov[kSigQ2PtSnp]*m24*m45-mParCov[kSigQ2Pt2]*m44*m45;
    double b5=m24*m35;
    double b6=m44*m35;
    mParCov[kSigTglSnp ] = (b4-b6*b1/b3)/(b5-b6*b2/b3);
    mParCov[kSigQ2PtTgl] = b1/b3-b2*mParCov[kSigTglSnp]/b3;
    mParCov[kSigTgl2 ] = fabs((cv[20]-mParCov[kSigQ2Pt2]*(m45*m45)-mParCov[kSigQ2PtTgl]*2.*m35*m45)/(m35*m35));
  }
  CheckCovariance();
}


//______________________________________________________________
bool Track::TrackPar::RotateParam(float alpha)
{
  // rotate to alpha frame
  if (fabs(GetSnp()) > kAlmost1) {
    //FairLogger::GetLogger()->Error(MESSAGE_ORIGIN, 
    printf("Precondition is not satisfied: |sin(phi)|>1 ! %f\n",GetSnp()); 
    return false;
  }
  //
  BringToPMPi(alpha);
  //
  float ca=0,sa=0;
  sincosf(alpha-GetAlpha(),sa,ca);
  float snp = GetSnp(), csp = sqrtf((1.f-snp)*(1.f+snp)); // Improve precision
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
  float xold = GetX(), yold = GetY();
  mParam[kAlpha] = alpha;
  mParam[kX]  =  xold*ca + yold*sa;
  mParam[kY]  = -xold*sa + yold*ca;
  mParam[kSnp]=  tmp;
  return true;
}

//______________________________________________________________
bool Track::TrackParCov::Rotate(float alpha)
{
  // rotate to alpha frame
  if (fabs(mParCov[kSnp]) > kAlmost1) {
    //FairLogger::GetLogger()->Error(MESSAGE_ORIGIN, 
    printf("Precondition is not satisfied: |sin(phi)|>1 ! %f\n",mParCov[kSnp]); 
    return false;
  }
  //
  BringToPMPi(alpha);
  //
  float ca=0,sa=0;
  sincosf(alpha-mParCov[kAlpha],sa,ca);
  float snp = mParCov[kSnp], csp = sqrtf((1.f-snp)*(1.f+snp)); // Improve precision
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
  float xold = mParCov[kX], yold = mParCov[kY];
  mParCov[kAlpha] = alpha;
  mParCov[kX]  =  xold*ca + yold*sa;
  mParCov[kY]  = -xold*sa + yold*ca;
  mParCov[kSnp]=  tmp;

  if (fabs(csp)<kAlmost0) {
    printf("Too small cosine value %f\n",csp);
    csp = kAlmost0;
  } 

  float rr=(ca+snp/csp*sa);  

  mParCov[kSigY2]      *= (ca*ca);
  mParCov[kSigZY]      *= ca;
  mParCov[kSigSnpY]    *= ca*rr;
  mParCov[kSigSnpZ]    *= rr;
  mParCov[kSigSnp2]    *= rr*rr;
  mParCov[kSigTglY]    *= ca;
  mParCov[kSigTglSnp]  *= rr;
  mParCov[kSigQ2PtY]   *= ca;
  mParCov[kSigQ2PtSnp] *= rr;

  CheckCovariance();
  return true;
}

//______________________________________________________________
void Track::TrackPar::InvertParam() 
{
  // Transform this track to the local coord. system rotated by 180 deg. 
  mParam[kX] = -mParam[kX];
  mParam[kAlpha] += kPI;
  BringToPMPi(mParam[kAlpha]);
  //
  mParam[0] = -mParam[0];
  mParam[3] = -mParam[3];
  mParam[4] = -mParam[4];
  //
}

//______________________________________________________________
void Track::TrackParCov::Invert() 
{
  // Transform this track to the local coord. system rotated by 180 deg. 
  InvertParam();
  // since the fP1 and fP2 are not inverted, their covariances with others change sign
  mParCov[kSigZY]      = -mParCov[kSigZY];
  mParCov[kSigSnpY]    = -mParCov[kSigSnpY];
  mParCov[kSigTglZ]    = -mParCov[kSigTglZ];
  mParCov[kSigTglSnp]  = -mParCov[kSigTglSnp];
  mParCov[kSigQ2PtZ]   = -mParCov[kSigQ2PtZ];
  mParCov[kSigQ2PtSnp] = -mParCov[kSigQ2PtSnp];
}

//____________________________________________________________
bool Track::TrackPar::PropagateParamTo(float xk, float b) 
{
  //----------------------------------------------------------------
  // Propagate this track to the plane X=xk (cm) in the field "b" (kG)
  // Only parameters are propagated, not the matrix. To be used for small 
  // distances only (<mm, i.e. misalignment)
  //----------------------------------------------------------------
  float dx=xk-mParam[kX];
  if (fabs(dx)<kAlmost0)  return true;
  float crv = (fabs(b)<kAlmost0) ? 0.f : GetCurvature(b);
  float x2r = crv*dx;
  float f1 = mParam[kSnp], f2=f1 + x2r;
  if (fabs(f1) > kAlmost1) return false;
  if (fabs(f2) > kAlmost1) return false;
  if (fabs(mParam[kQ2Pt])< kAlmost0) return false;
  float r1=sqrtf((1.f-f1)*(1.f+f1)), r2=sqrtf((1.f-f2)*(1.f+f2));
  if (fabs(r1)<kAlmost0)  return false;
  if (fabs(r2)<kAlmost0)  return false;
  mParam[kX] = xk;
  double dy2dx = (f1+f2)/(r1+r2);
  mParam[kY] += dx*dy2dx;
  mParam[kSnp] += x2r;
  if (fabs(x2r)<0.05f) mParam[kZ] += dx*(r2 + f2*dy2dx)*mParam[kTgl];
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
    mParam[kZ] += mParam[kTgl]/crv*rot; 
  }
  return true;
}

//______________________________________________________________
bool Track::TrackParCov::PropagateTo(float xk, float b) 
{
  //----------------------------------------------------------------
  // Propagate this track to the plane X=xk (cm) in the field "b" (kG)
  //----------------------------------------------------------------
  float dx=xk-mParCov[kX];
  if (fabs(dx)< kAlmost0)  return true;      
  float crv = (fabs(b)<kAlmost0) ? 0.f : GetCurvature(b);
  float x2r = crv*dx;
  float f1 = mParCov[kSnp], f2=f1 + x2r;
  if (fabs(f1) > kAlmost1) return false;
  if (fabs(f2) > kAlmost1) return false;
  if (fabs(mParCov[kQ2Pt])< kAlmost0) return false;
  float r1=sqrtf((1.f-f1)*(1.f+f1)), r2=sqrtf((1.f-f2)*(1.f+f2));
  if (fabs(r1)<kAlmost0)  return false;
  if (fabs(r2)<kAlmost0)  return false;
  mParCov[kX] = xk;
  double dy2dx = (f1+f2)/(r1+r2);
  mParCov[kY] += dx*dy2dx;
  mParCov[kSnp] += x2r;
  if (fabs(x2r)<0.05f) mParCov[kZ] += dx*(r2 + f2*dy2dx)*mParCov[kTgl];
  else { 
    // for small dx/R the linear apporximation of the arc by the segment is OK,
    // but at large dx/R the error is very large and leads to incorrect Z propagation
    // angle traversed delta = 2*asin(dist_start_end / R / 2), hence the arc is: R*deltaPhi
    // The dist_start_end is obtained from sqrt(dx^2+dy^2) = x/(r1+r2)*sqrt(2+f1*f2+r1*r2)
    //    double chord = dx*TMath::Sqrt(1+dy2dx*dy2dx);   // distance from old position to new one
    //    double rot = 2*TMath::ASin(0.5*chord*crv); // angular difference seen from the circle center
    //    mParCov1 += rot/crv*mParCov3;
    // 
    float rot = asinf(r1*f2 - r2*f1); // more economic version from Yura.
    if (f1*f1+f2*f2>1.f && f1*f2<0.f) {          // special cases of large rotations or large abs angles
      if (f2>0.f) rot = kPI - rot;    //
      else       rot = -kPI - rot;
    }
    mParCov[kZ] += mParCov[kTgl]/crv*rot; 
  }
    
  float
    &c00=mParCov[kSigY2],
    &c10=mParCov[kSigZY],    &c11=mParCov[kSigZ2],
    &c20=mParCov[kSigSnpY],  &c21=mParCov[kSigSnpZ],  &c22=mParCov[kSigSnp2],
    &c30=mParCov[kSigTglY],  &c31=mParCov[kSigTglZ],  &c32=mParCov[kSigTglSnp],  &c33=mParCov[kSigTgl2],  
    &c40=mParCov[kSigQ2PtY], &c41=mParCov[kSigQ2PtZ], &c42=mParCov[kSigQ2PtSnp], &c43=mParCov[kSigQ2PtTgl], &c44=mParCov[kSigQ2Pt2];
    
  // evaluate matrix in double prec.
  double rinv  = 1./r1;
  double r3inv = rinv*rinv*rinv;
  double f24   = dx*b*kB2C; // x2r/mParCov[kQ2Pt];
  double f02   = dx*r3inv;
  double f04   = 0.5*f24*f02;
  double f12   = f02*mParCov[kTgl]*f1;
  double f14   = 0.5*f24*f02*mParCov[kTgl]*f1;
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
    
  CheckCovariance();
    
  return true;
}

//____________________________________________________________
bool Track::TrackPar::PropagateParamBxByBzTo(float xk, const float b[3])
{
  //----------------------------------------------------------------
  // Extrapolate this track params (w/o cov matrix) to the plane X=xk in the field b[].
  //
  // X [cm] is in the "tracking coordinate system" of this track.
  // b[]={Bx,By,Bz} [kG] is in the Global coordidate system.
  //----------------------------------------------------------------

  float dx=xk-mParam[kX];
  if (fabs(dx)<kAlmost0)  return true;
  // Do not propagate tracks outside the ALICE detector
  if (fabs(dx)>1e5 || fabs(mParam[kY])>1e5 || fabs(mParam[kZ])>1e5) {
    printf("Anomalous track, target X:%f\n",xk);
    //    Print();
    return false;
  }
  float crv = (fabs(b[2])<kAlmost0) ? 0.f : GetCurvature(b[2]);
  float x2r = crv*dx;
  float f1 = mParam[kSnp], f2 = f1 + x2r;
  if (fabs(f1)>kAlmost1 || fabs(f2)>kAlmost1) return false;
  if (fabs(mParam[kQ2Pt])<kAlmost0) return false;
  float r1=sqrtf((1.f-f1)*(1.f+f1)), r2=sqrtf((1.f-f2)*(1.f+f2));
  if (fabs(r1)<kAlmost0 || fabs(r2)<kAlmost0)  return false;
  float dy2dx = (f1+f2)/(r1+r2);
  float step = (fabs(x2r)<0.05f) ? dx*fabs(r2 + f2*dy2dx)      // chord
    : 2.f*asinf(0.5f*dx*sqrtf(1.f+dy2dx*dy2dx)*crv)/crv;       // arc
  step *= sqrtf(1.f+ mParam[kTgl]*mParam[kTgl]);
  //
  // Get the track x,y,z,px/p,py/p,pz/p,p,sinAlpha,cosAlpha in the Global System
  float vecLab[9];
  if (!GetPosDir(vecLab)) return false;

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
  float sgn = GetSign();
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
  mParam[kX]    = x;
  mParam[kY]    = y;
  mParam[kZ]    = z;
  mParam[kSnp]  = vecLab[4]*t;
  mParam[kTgl]  = vecLab[5]*t; 
  mParam[kQ2Pt] = sgn*t/vecLab[6];

  return true;
}

//____________________________________________________________
bool Track::TrackParCov::PropagateBxByBzTo(float xk, const float b[3])
{
  //----------------------------------------------------------------
  // Extrapolate this track to the plane X=xk in the field b[].
  //
  // X [cm] is in the "tracking coordinate system" of this track.
  // b[]={Bx,By,Bz} [kG] is in the Global coordidate system.
  //----------------------------------------------------------------

  float dx=xk-GetX();
  if (fabs(dx)<kAlmost0)  return true;
  // Do not propagate tracks outside the ALICE detector
  if (fabs(dx)>1e5 || fabs(GetY())>1e5 || fabs(GetZ())>1e5) {
    printf("Anomalous track, target X:%f\n",xk);
    //    Print();
    return false;
  }
  float crv = (fabs(b[2])<kAlmost0) ? 0.f : GetCurvature(b[2]);
  float x2r = crv*dx;
  float f1 = mParCov[kSnp], f2 = f1 + x2r;
  if (fabs(f1)>kAlmost1 || fabs(f2)>kAlmost1) return false;
  if (fabs(mParCov[kQ2Pt])<kAlmost0) return false;
  float r1=sqrtf((1.f-f1)*(1.f+f1)), r2=sqrtf((1.f-f2)*(1.f+f2));
  if (fabs(r1)<kAlmost0 || fabs(r2)<kAlmost0)  return false;
  float dy2dx = (f1+f2)/(r1+r2);
  float step = (fabs(x2r)<0.05f) ? dx*fabs(r2 + f2*dy2dx)      // chord
    : 2.f*asinf(0.5f*dx*sqrtf(1.f+dy2dx*dy2dx)*crv)/crv;       // arc
  step *= sqrtf(1.f+ GetTgl()*GetTgl());
  //
  // Get the track x,y,z,px/p,py/p,pz/p,p,sinAlpha,cosAlpha in the Global System
  float vecLab[9];
  if (!Param()->GetPosDir(vecLab)) return false;
  //
  // matrix transformed with Bz component only
  float
    &c00=mParCov[kSigY2],
    &c10=mParCov[kSigZY],   &c11=mParCov[kSigZ2],
    &c20=mParCov[kSigSnpY], &c21=mParCov[kSigSnpZ], &c22=mParCov[kSigSnp2],
    &c30=mParCov[kSigTglY], &c31=mParCov[kSigTglZ], &c32=mParCov[kSigTglSnp], &c33=mParCov[kSigTgl2],  
    &c40=mParCov[kSigQ2PtY],&c41=mParCov[kSigQ2PtZ],&c42=mParCov[kSigQ2PtSnp],&c43=mParCov[kSigQ2PtTgl],&c44=mParCov[kSigQ2Pt2];
  // evaluate matrix in double prec.
  double rinv  = 1./r1;
  double r3inv = rinv*rinv*rinv;
  double f24   = dx*b[2]*kB2C; // x2r/track[kQ2Pt];
  double f02   = dx*r3inv;
  double f04   = 0.5*f24*f02;
  double f12   = f02*GetTgl()*f1;
  double f14   = 0.5*f24*f02*GetTgl()*f1;
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
    
  CheckCovariance();
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
  float sgn = GetSign();
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
  mParCov[kX]    = x;
  mParCov[kY]    = y;
  mParCov[kZ]    = z;
  mParCov[kSnp]  = vecLab[4]*t;
  mParCov[kTgl]  = vecLab[5]*t; 
  mParCov[kQ2Pt] = sgn*t/vecLab[6];

  return true;

}


//______________________________________________
void Track::TrackParCov::CheckCovariance() 
{
  // This function forces the diagonal elements of the covariance matrix to be positive.
  // In case the diagonal element is bigger than the maximal allowed value, it is set to
  // the limit and the off-diagonal elements that correspond to it are set to zero.
  
  mParCov[kSigY2] = fabs(mParCov[kSigY2]);
  if (mParCov[kSigY2]>kCY2max) {
    float scl = sqrtf(kCY2max/mParCov[kSigY2]);
    mParCov[kSigY2]     = kCY2max;
    mParCov[kSigZY]    *= scl;
    mParCov[kSigSnpY]  *= scl;
    mParCov[kSigTglY]  *= scl;
    mParCov[kSigQ2PtY] *= scl;
  }
  mParCov[kSigZ2] = fabs(mParCov[kSigZ2]);
  if (mParCov[kSigZ2]>kCZ2max) {
    float scl = sqrtf(kCZ2max/mParCov[kSigZ2]);
    mParCov[kSigZ2]     = kCZ2max;
    mParCov[kSigZY]    *= scl;
    mParCov[kSigSnpZ]  *= scl;
    mParCov[kSigTglZ]  *= scl;
    mParCov[kSigQ2PtZ] *= scl;
  }
  mParCov[kSigSnp2] = fabs(mParCov[kSigSnp2]);
  if (mParCov[kSigSnp2]>kCSnp2max) {
    float scl = sqrtf(kCSnp2max/mParCov[kSigSnp2]);
    mParCov[kSigSnp2] = kCSnp2max;
    mParCov[kSigSnpY] *= scl;
    mParCov[kSigSnpZ] *= scl;
    mParCov[kSigTglSnp] *= scl;
    mParCov[kSigQ2PtSnp] *= scl;
  }
  mParCov[kSigTgl2] = fabs(mParCov[kSigTgl2]);
  if (mParCov[kSigTgl2]>kCTgl2max) {
    float scl = sqrtf(kCTgl2max/mParCov[kSigTgl2]);
    mParCov[kSigTgl2] = kCTgl2max;
    mParCov[kSigTglY] *= scl;
    mParCov[kSigTglZ] *= scl;
    mParCov[kSigTglSnp] *= scl;
    mParCov[kSigQ2PtTgl] *= scl;
  }
  mParCov[kSigQ2Pt2] = fabs(mParCov[kSigQ2Pt2]);
  if (mParCov[kSigQ2Pt2]>kCQ2Pt2max) {
    float scl = sqrtf(kCQ2Pt2max/mParCov[kSigQ2Pt2]);
    mParCov[kSigQ2Pt2] = kCQ2Pt2max;
    mParCov[kSigQ2PtY] *= scl;
    mParCov[kSigQ2PtZ] *= scl;
    mParCov[kSigQ2PtSnp] *= scl;
    mParCov[kSigQ2PtTgl] *= scl;
  }
}


//=================================================
//
// Aux. methods for tracks manipulation
//
//=================================================


//_______________________________________________________      
bool Track::TrackPar::GetPxPyPz(float pxyz[3]) const 
{
  // track momentum
  if (fabs(GetQ2Pt())<kAlmost0 || fabs(GetSnp())>kAlmost1) return false;
  float cs,sn, pt=fabs(1.f/GetQ2Pt());
  float r = sqrtf((1.f - GetSnp())*(1.f + GetSnp()));
  sincosf(GetAlpha(),sn,cs);
  pxyz[0] = pt*(r*cs - GetSnp()*sn); 
  pxyz[1] = pt*(GetSnp()*cs + r*sn); 
  pxyz[2] = pt*GetTgl();
  return true;
}

//____________________________________________________
bool Track::TrackPar::GetPosDir(float posdirp[9]) const
{
  // fill vector with lab x,y,z,px/p,py/p,pz/p,p,sinAlpha,cosAlpha
  float ptI = fabs(GetQ2Pt());
  float snp = GetSnp();
  if (ptI<kAlmost0 || fabs(snp)>kAlmost1) return false;
  float &sn=posdirp[7],&cs=posdirp[8]; 
  float csp = sqrtf((1.f - snp)*(1.f + snp));
  float cstht = sqrtf(1.f+ GetTgl()*GetTgl());
  float csthti = 1.f/cstht;
  sincosf(GetAlpha(),sn,cs);
  posdirp[0] = GetX()*cs - GetY()*sn;
  posdirp[1] = GetX()*sn + GetY()*cs;
  posdirp[2] = GetZ();
  posdirp[3] = (csp*cs - snp*sn)*csthti;  // px/p
  posdirp[4] = (snp*cs + csp*sn)*csthti;  // py/p
  posdirp[5] = GetTgl()*csthti;           // pz/p
  posdirp[6] = cstht/ptI;                 // p
  return true;
}

void Track::g3helx3(float qfield, 
		    float step,
		    float vect[7]) {
/******************************************************************
 *                                                                *
 *       GEANT3 tracking routine in a constant field oriented     *
 *       along axis 3                                             *
 *       Tracking is performed with a conventional                *
 *       helix step method                                        *
 *                                                                *
 *       Authors    R.Brun, M.Hansroul  *********                 *
 *       Rewritten  V.Perevoztchikov                              *
 *                                                                *
 *       Rewritten in C++ by I.Belikov                            *
 *                                                                *
 *  qfield (kG)       - particle charge times magnetic field      *
 *  step   (cm)       - step length along the helix               *
 *  vect[7](cm,GeV/c) - input/output x, y, z, px/p, py/p ,pz/p, p *
 *                                                                *
 ******************************************************************/
  const int ix=0, iy=1, iz=2, ipx=3, ipy=4, ipz=5, ipp=6;
  const float kOvSqSix=sqrtf(1./6.);

  float cosx=vect[ipx], cosy=vect[ipy], cosz=vect[ipz];

  float rho = qfield*kB2C/vect[ipp]; 
  float tet = rho*step;

  float tsint, sintt, sint, cos1t; 
  if (fabs(tet) > 0.03f) {
     sint  = sinf(tet);
     sintt = sint/tet;
     tsint = (tet - sint)/tet;
     float t=sinf(0.5f*tet);
     cos1t = 2*t*t/tet;
  } else {
     tsint = tet*tet/6.f;
     sintt = (1.f-tet*kOvSqSix)*(1.f+tet*kOvSqSix); // 1.- tsint;
     sint  = tet*sintt;
     cos1t = 0.5f*tet; 
  }

  float f1 = step*sintt;
  float f2 = step*cos1t;
  float f3 = step*tsint*cosz;
  float f4 = -tet*cos1t;
  float f5 = sint;

  vect[ix]  += f1*cosx - f2*cosy;
  vect[iy]  += f1*cosy + f2*cosx;
  vect[iz]  += f1*cosz + f3;

  vect[ipx] += f4*cosx - f5*cosy;
  vect[ipy] += f4*cosy + f5*cosx;  

}
