#ifndef KALMANFILTER_H
#define KALMANFILTER_H

#include "TrackRepresentation.h"

class KalmanFilter{

 public:
  KalmanFilter();
  ~KalmanFilter();
  void AddHit(double r, double phi, double z, double re, double phie, double ze);
  void Process();
  void CalcChi2();
  

 public:
  // The track representation object
  TrackRepresentation trackrep;
  // the actual measurement values with uncertainties
  // to be used during the fit
  std::vector<double> m_r, m_phi, m_z, m_re, m_ze, m_phie;
  // the measurement covariance matrix
  TMatrixD V;

  //Kalman gain
  TMatrixD K;

  // chi squared
  double chi2;

};



#endif
