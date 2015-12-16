/// \file TrackP
/// \brief Base track model for the Barrel, params only, w/o covariance
/// \author ruben.shahoyan@cern.ch

#ifndef ALICEO2_BASE_TRACK
#define ALICEO2_BASE_TRACK

#include "Constants.h"

namespace AliceO2 {
  namespace Base {
    namespace Track {
      
      using AliceO2::Base::Constants;
      
      // aliases for track elements
      enum {kX,kAlpha,
	    kY,kZ,kSnp,kTgl,kQ2Pt,
	    kSigY2,
	    kSigZY,kSigZ2,
	    kSigSnpY,kSigSnpZ,kSigSnp2,
	    kSigTglY,kSigTglZ,kSigTglSnp,kSigTgl2,
	    kSigQ2PtY,kSigQ2PtZ,kSigQ2PtSnp,kSigQ2PtTgl,kSigQ2Pt2};
      enum {kNParams=5,kNCovMatSize=15,kTrackPSize=kNParams+2,kTrackPCSize=kTrackPSize+kNCovMatSize};

      const float 
	kCY2max=100*100, // SigmaY<=100cm
	kCZ2max=100*100, // SigmaZ<=100cm
	kCSnp2max=1*1,     // SigmaSin<=1
	kCTgl2max=1*1,     // SigmaTan<=1
	kCQ2Pt2max=100*100; // Sigma1/Pt<=100 1/GeV

  
      class TrackPar { // track parameterization, kinematics only
      public:
	TrackPar(float x,float alpha, const float *par);
	TrackPar(const TrackPar& src);
	~TrackPar() {}
	TrackPar& operator=(const TrackPar& src);

	float& operator[](int i)                   { return mParam[i]; }
	float  operator[](int i)             const { return mParam[i]; }
	operator float*()                    const { return (float*)mParam; }

	const float* GetParam()              const { return &mParam[kY]; }
	float GetX()                         const { return mParam[kX]; }
	float GetAlpha()                     const { return mParam[kAlpha]; }
	float GetY()                         const { return mParam[kY]; }
	float GetZ()                         const { return mParam[kZ]; }
	float GetSnp()                       const { return mParam[kSnp]; }
	float GetTgl()                       const { return mParam[kTgl]; }
	float GetQ2Pt()                      const { return mParam[kQ2Pt]; }
	float GetCurvature(float b)          const { return mParam[kQ2Pt]*b*kB2C;}
  
      protected:
	float mParam[kTrackPSize];  // x,alpha + 5 parameters
      };

      class TrackParCov { // track+error parameterization
      public:
	TrackParCov(float x,float alpha, const float *par, const float* cov);
	TrackParCov(const TrackParCov& src);
	~TrackParCov() {}
	TrackParCov& operator=(const TrackParCov& src);

	operator TrackPar*() { return reinterpret_cast<TrackPar*>(this); }
	operator TrackPar()  { return *reinterpret_cast<TrackPar*>(this); }
	operator TrackPar&() { return *reinterpret_cast<TrackPar*>(this); }

	float& operator[](int i)                   { return mParCov[i]; }
	float  operator[](int i)             const { return mParCov[i]; }
	operator float*()                    const { return (float*)mParCov; }
	const float* GetParam()              const { return &mParCov[kY]; }
	const float* GetCov()                const { return &mParCov[kSigY2]; }

	float GetX()                         const { return mParCov[kX]; }
	float GetAlpha()                     const { return mParCov[kAlpha]; }
	float GetY()                         const { return mParCov[kY]; }
	float GetZ()                         const { return mParCov[kZ]; }
	float GetSnp()                       const { return mParCov[kSnp]; }
	float GetTgl()                       const { return mParCov[kTgl]; }
	float GetQ2Pt()                      const { return mParCov[kQ2Pt]; }

	float GetSigmaY2()                   const { return mParCov[kSigY2]; }
	float GetSigmaZY()                   const { return mParCov[kSigZY]; }
	float GetSigmaZ2()                   const { return mParCov[kSigZ2]; }
	float GetSigmaSnpY()                 const { return mParCov[kSigSnpY]; }
	float GetSigmaSnpZ()                 const { return mParCov[kSigSnpZ]; }
	float GetSigmaSnp2()                 const { return mParCov[kSigSnp2]; }
	float GetSigmaTglY()                 const { return mParCov[kSigTglY]; }
	float GetSigmaTglZ()                 const { return mParCov[kSigTglZ]; }
	float GetSigmaTglSnp()               const { return mParCov[kSigTglSnp]; }
	float GetSigmaTgl2()                 const { return mParCov[kSigTgl2]; }
	float GetSigma1PtY()                 const { return mParCov[kSigQ2PtY]; }
	float GetSigma1PtZ()                 const { return mParCov[kSigQ2PtZ]; }
	float GetSigma1PtSnp()               const { return mParCov[kSigQ2PtSnp]; }
	float GetSigma1PtTgl()               const { return mParCov[kSigQ2PtTgl]; }
	float GetSigma1Pt2()                 const { return mParCov[kSigQ2Pt2]; }
	float GetCurvature(float b)          const { return mParCov[kQ2Pt]*b*kB2C;}

      protected:
	float mParCov[kTrackPCSize];  // x, alpha + 5 parameters + 15 errors
      };


      //____________________________________________________________
      inline TrackPar::TrackPar(float x, float alpha, const float *par) {
	// explicit constructor
	mParam[kX] = x;
	mParam[kAlpha] = alpha;
	memcpy(&mParam[kY],par,kNParams*sizeof(float));
      }

      //____________________________________________________________
      inline TrackPar::TrackPar(const TrackPar& src) {
	// copy c-tor
	memcpy(mParam,src.mParam,kTrackPSize*sizeof(float));
      }

      //____________________________________________________________
      inline TrackPar& TrackPar::operator=(const TrackPar& src) {
	// assignment operator
	if (this!=&src) memcpy(mParam,src.mParam,kTrackPSize*sizeof(float)); 
	return this;
      }

      //____________________________________________________________
      inline TrackParCov::TrackParCov(float x, float alpha, const float *par, const float *cov) {
	// explicit constructor
	mParCov[kX] = x;
	mParCov[kAlpha] = alpha;
	memcpy(&mParCov[kY],par,kNParams*sizeof(float));
	memcpy(&mParCov[kSigY2],cov,kNCovMatSize*sizeof(float));
      }

      //____________________________________________________________
      inline TrackParCov::TrackParCov(const TrackParCov& src) {
	// copy c-tor
	memcpy(mParCov,src.mParCov,kTrackPCSize*sizeof(float));
      }

      //____________________________________________________________
      inline TrackParCov& TrackParCov::operator=(const TrackParCov& src) {
	// assignment operator
	if (this!=&src) memcpy(mParCov,src.mParCov,kTrackPSize*sizeof(float)); 
	return this;
      }

      //===========================================================
      //
      //           Track manipulation methods
      //
      //===========================================================
      bool RotateParam(TrackPar& track, float alpha);
      bool Rotate(TrackParCov& track, float alpha);

      bool PropagateParamTo(TrackPar &track,float xk, float b);
      bool PropagateTo(TrackParCov &track, float xk, float b) 

      void InvertParam(TrackPar& track);
      void Invert(TrackParCov& track);
      
      
    }  
  }
}


#endif
