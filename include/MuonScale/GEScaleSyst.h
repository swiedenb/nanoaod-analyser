#ifndef HighPtMuonAnalysis_GEScale_H
#define HighPtMuonAnalysis_GEScale_H
// #pragma once
#include <iostream>
#include <cmath>
#include <map>

#include "TLorentzVector.h"

class GEScaleSyst {
  public:
    GEScaleSyst();
    virtual ~GEScaleSyst();

    float GEScaleCorrPt(int icopy, float pt, float eta, float phi, int charge, bool doOpp = false);
    TLorentzVector GEScaleCorrLvec(int icopy, float pt, float eta, float phi, int charge, bool doOpp = false);

    void SetVerbose(int _v) {
      verbose = _v;
    };

    void AddNewMatrix( int new_copy, std::map<int, std::map<int, float> > new_matrix );

  private:
    const double mu_mass = 0.105658;

    float GetGEScaleKappa(int icopy, float eta, float phi);

    const int neta = 6;
    const double etaBinEdge[7] = { -2.4, -2.1, -1.2, 0, 1.2, 2.1, 2.4 };
    const int nphi = 3;
    const double phiBinEdge[4] = { -M_PI, -M_PI/3.0, M_PI/3.0, M_PI };

    // _matrix[copy #][eta bin #][phi bin #]

    // copy #:
    //   yr: 16, 17, 18
    //   yr00:  kb (original)
    //   yr01:  kb + 1 sigma variation
    //   yr02:  kb - 1 sigma variation
    //   yr03:  kb + 2 sigma variation
    //   yr04:  kb - 2 sigma variation
    //   yr0000, yr0001, ...: random copies using gaus(kappa, sigma)
    //   yr1000, yr1001, ...: random copies passing sign contraint, restrict sign of bias if (k + 1 sigma)(k - 1 sigma) > 0

    // eta bin #:
    //   0: [-2.4, -2.1]
    //   1: [-2.1, -1.2]
    //   2: [-1.2, 0.0]
    //   3: [0.0, 1.2]
    //   4: [1.2, 2.1]
    //   5: [2.1, 2.4]

    // phi bin #:
    //   0: [-180, -60]
    //   1: [-60, 60]
    //   2: [60, 180]

    std::map<int, std::map<int, std::map<int, float> > > _matrix;

    int verbose;
};

#endif
