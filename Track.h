/// \file TrackP
/// \brief Base track model for the Barrel, params only, w/o covariance
/// \author ruben.shahoyan@cern.ch

#ifndef ALICEO2_BASE_TRACK
#define ALICEO2_BASE_TRACK

#include <iostream>
#include <stdio.h>
#include <string.h>
#include "Constants.h"
#include "Utils.h"

namespace AliceO2 {
  namespace Base {
    namespace Track {
      
      using namespace AliceO2::Base::Constants;
      using namespace AliceO2::Base::Utils;
      using namespace std;
      
      // aliases for track elements
      enum {kY,kZ,kSnp,kTgl,kQ2Pt};
      enum {kSigY2,
	    kSigZY,kSigZ2,
	    kSigSnpY,kSigSnpZ,kSigSnp2,
	    kSigTglY,kSigTglZ,kSigTglSnp,kSigTgl2,
	    kSigQ2PtY,kSigQ2PtZ,kSigQ2PtSnp,kSigQ2PtTgl,kSigQ2Pt2};
      enum {kNParams=5,kCovMatSize=15,kLabCovMatSize=21};
      
      const float 
	kCY2max=100*100, // SigmaY<=100cm
	kCZ2max=100*100, // SigmaZ<=100cm
	kCSnp2max=1*1,     // SigmaSin<=1
	kCTgl2max=1*1,     // SigmaTan<=1
	kC1Pt2max=100*100; // Sigma1/Pt<=100 1/GeV
      
      
      class TrackPar { // track parameterization, kinematics only
      public:

	const float* GetParam()              const { return mP; }
	float GetX()                         const { return mX; }
	float GetAlpha()                     const { return mAlpha; }
	float GetY()                         const { return mP[kY]; }
	float GetZ()                         const { return mP[kZ]; }
	float GetSnp()                       const { return mP[kSnp]; }
	float GetTgl()                       const { return mP[kTgl]; }
	float GetQ2Pt()                      const { return mP[kQ2Pt]; }

	// derived getters
	float GetCurvature(float b)          const { return mP[kQ2Pt]*b*kB2C;}
	float GetSign()                      const { return mP[kQ2Pt]>0 ? 1.f:-1.f;}
	float GetPhi()                       const { return asinf(GetSnp()) + GetAlpha();}
	float GetPhiPos()                    const;

	float GetP()                         const;
	float GetPt()                        const;
	void  GetXYZ(float xyz[3])           const;
	bool  GetPxPyPz(float pxyz[3])       const;
	bool  GetPosDir(float posdirp[9])    const;

	// parameters manipulation
	bool  RotateParam(float alpha);
	bool  PropagateParamTo(float xk, float b);
	bool  PropagateParamTo(float xk, const float b[3]);
	void  InvertParam();

	void  PrintParam()                   const;

      protected:
	// to keep this class non-virtual but derivable the c-tors and d-tor are protected
        TrackPar() : mX(0),mAlpha(0) {memset(mP,0,kNParams*sizeof(float));}
	TrackPar(float x,float alpha, const float par[kNParams]);
	TrackPar(const float xyz[3],const float pxpypz[3],int sign, bool sectorAlpha=true);
	TrackPar(const TrackPar& src);
	~TrackPar() {}
	TrackPar& operator=(const TrackPar& src);
	//
	static void g3helx3(float qfield, float step,float vect[7]);
	
      protected:
	float mX;               // X of track evaluation
	float mAlpha;           // track frame angle
	float mP[kNParams];     // 5 parameters: Y,Z,sin(phi),tg(lambda),q/pT
      };

      // rootcint does not swallow final keyword here
      class TrackParCov /*final*/ : public TrackPar { // track+error parameterization
      public:
	TrackParCov() { memset(mC,0,kCovMatSize*sizeof(float)); }
	TrackParCov(float x,float alpha, const float par[kNParams], const float cov[kCovMatSize]);
	TrackParCov(const float xyz[3],const float pxpypz[3],const float[kLabCovMatSize], int sign, bool sectorAlpha=true);
	TrackParCov(const TrackParCov& src);
	~TrackParCov() {}
	TrackParCov& operator=(const TrackParCov& src);

	const float* GetCov()                const { return mC; }	
	float GetSigmaY2()                   const { return mC[kSigY2]; }
	float GetSigmaZY()                   const { return mC[kSigZY]; }
	float GetSigmaZ2()                   const { return mC[kSigZ2]; }
	float GetSigmaSnpY()                 const { return mC[kSigSnpY]; }
	float GetSigmaSnpZ()                 const { return mC[kSigSnpZ]; }
	float GetSigmaSnp2()                 const { return mC[kSigSnp2]; }
	float GetSigmaTglY()                 const { return mC[kSigTglY]; }
	float GetSigmaTglZ()                 const { return mC[kSigTglZ]; }
	float GetSigmaTglSnp()               const { return mC[kSigTglSnp]; }
	float GetSigmaTgl2()                 const { return mC[kSigTgl2]; }
	float GetSigma1PtY()                 const { return mC[kSigQ2PtY]; }
	float GetSigma1PtZ()                 const { return mC[kSigQ2PtZ]; }
	float GetSigma1PtSnp()               const { return mC[kSigQ2PtSnp]; }
	float GetSigma1PtTgl()               const { return mC[kSigQ2PtTgl]; }
	float GetSigma1Pt2()                 const { return mC[kSigQ2Pt2]; }

	void  Print()                        const;

	// parameters + covmat manipulation
	bool  Rotate(float alpha);
	bool  PropagateTo(float xk, float b);
	bool  PropagateTo(float xk, const float b[3]);
	void  Invert();

	float GetPredictedChi2(const float p[2], const float cov[3]) const;
	bool  Update(const float p[2], const float cov[3]);

	bool  CorrectForMaterial(float x2x0,float xrho,float mass,bool anglecorr=false,float dedx=kCalcdEdxAuto);

	void  ResetCovariance(float s2=0);
	void  CheckCovariance();

      protected:
	float mC[kCovMatSize];  // x, alpha + 5 parameters + 15 errors

	static const float kCalcdEdxAuto; // value indicating request for dedx calculation
      };

      class TrackParOnly final : public TrackPar { // track parameterization only
      public:
	TrackParOnly() {}
        TrackParOnly(float x,float alpha, const float par[kNParams]) : TrackPar(x,alpha,par) {}
        TrackParOnly(const float xyz[3],const float pxpypz[3],int sign, bool sectorAlpha=true) : TrackPar(xyz,pxpypz,sign,sectorAlpha) {}
        TrackParOnly(const TrackParOnly& src) : TrackPar(src) {}
        TrackParOnly(const TrackParCov& src) : TrackPar((const TrackPar&)src) {}
	~TrackParOnly() {}
	TrackParOnly& operator=(const TrackParOnly& src) {TrackPar::operator=((const TrackPar&)src); return *this;}
	//
	void  Print() const {PrintParam();}
      };

      //____________________________________________________________
      inline TrackPar::TrackPar(float x, float alpha, const float par[kNParams]) : mX(x), mAlpha(alpha) {
	// explicit constructor
	memcpy(mP,par,kNParams*sizeof(float));
      }

      //____________________________________________________________
      inline TrackPar::TrackPar(const TrackPar& src) : mX(src.mX), mAlpha(src.mAlpha) {
	// copy c-tor
	memcpy(mP,src.mP,kNParams*sizeof(float));
      }

      //____________________________________________________________
      inline TrackPar& TrackPar::operator=(const TrackPar& src) {
	// assignment operator
	if (this!=&src) memcpy(this,&src,sizeof(TrackPar));
	return *this;
      }
      
      //_______________________________________________________
      inline void TrackPar::GetXYZ(float xyz[3]) const {
	// track coordinates in lab frame
	xyz[0] = GetX(); 
	xyz[1] = GetY();
	xyz[2] = GetZ();
	RotateZ(xyz,GetAlpha());
      }

      //_______________________________________________________
      inline float TrackPar::GetPhiPos() const {
	// angle of track position
	float xy[2]={GetX(),GetY()};
	return atan2(xy[1],xy[0]);
      }

      //____________________________________________________________
      inline float TrackPar::GetP() const {
	// return the track momentum
	float ptI = fabs(GetQ2Pt());
	return (ptI>kAlmost0) ? sqrtf(1.f+ GetTgl()*GetTgl())/ptI : kVeryBig;
      }

      //____________________________________________________________
      inline float TrackPar::GetPt() const {
	// return the track transverse momentum
	float ptI = fabs(GetQ2Pt());
	return (ptI>kAlmost0) ? 1.f/ptI : kVeryBig;
      }

      //============================================================

      //____________________________________________________________
      inline TrackParCov::TrackParCov(float x, float alpha, const float par[kNParams], const float cov[kCovMatSize]) : TrackPar(x,alpha,par) {
	// explicit constructor
	memcpy(mC,cov,kCovMatSize*sizeof(float));
      }

      //____________________________________________________________
      inline TrackParCov::TrackParCov(const TrackParCov& src) : TrackPar(src) {
	// copy c-tor
	memcpy(mC,src.GetCov(),sizeof(kCovMatSize));
      }

      //____________________________________________________________
      inline TrackParCov& TrackParCov::operator=(const TrackParCov& src) {
	// assignment operator
	if (this!=&src) memcpy(mC,src.mC,kCovMatSize*sizeof(float)); 
	return *this;
      }

    }  
  }
}


#endif
