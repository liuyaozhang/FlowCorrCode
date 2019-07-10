#ifndef __MYANACONST_h__
#define __MYANACONST_h__

#include <iostream>
#include <string>
#include "TMath.h"

namespace ana{

   const double PI = TMath::Pi();

   const int nEtaBin = 33;
   const double etaBegin = -4.95;
   const double etaEnd = 4.95;

   const int nPhiBin = 32;
//   const double phiBegin = -(0.5-1.0/32)*PI;
//   const double phiEnd = (1.5-1.0/32)*PI;
   const double phiBegin = -0.5*PI;
   const double phiEnd = 1.5*PI;

   const double etaMin_ass_ = -2.4;
   const double etaMax_ass_ = 2.4;
   const double ptMin_ass_ = 0.3;
   const double ptMax_ass_ = 3.0;
   const double multMin_ = 100;//185
   const double multMax_ = 250;
   const int    bkgFactor_ = 20;
   const double mis_d0_range_ = 0.020;
   const bool   rejectDaughter_ = true;
   const double d0_eta_ = 9999.0;
   const double d0_dau_eta_ = 1.5;
   const int    d0_dau_nhit_ = 11;
   const double d0_dau_pterr_ = 0.1;
   const double d0_dau_pt_ = 0.7;
   const double d0_DCA_ = 999999.;
   const double d0_rapidity_min_ = -1.;
   const double d0_rapidity_max_ = 1;

   const int nMass = 14;
   const int nPt = 4;
   const int nY = 1;
   const double massbin[nMass+1] = {1.70, 1.74, 1.78, 1.8, 1.82, 1.84, 1.85, 1.86, 1.865, 
      1.87, 1.88, 1.9, 1.92, 1.96, 2.0};
   //const double massbin[nMass+1] = {1.86, 1.865};
   const double ptbin[nPt+1] = {1.5, 3., 5., 7., 10};
   const double ybin[nY+1] = {0., 2.};
   const double mvaCut[nPt] = {0.52, 0.52, 0.52, 0.52};
   const double mvaCut_PD0[nPt] = {0.52, 0.52, 0.52, 0.52}; 
   const double mvaCut_NPD0[nPt] = {0.52, 0.52, 0.52, 0.52}; 
   const double dcaSep = 0.008; 
   const unsigned int nMixedEvts = 10;

   std::string findDCA(const double&);
   int findMassBin(const double&);
   int findPtBin(const double&);
   int findYBin(const double&);
};

#endif
