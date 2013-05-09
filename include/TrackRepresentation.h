#ifndef TRACKREPRESENTATION_H
#define TRACKREPRESENTATION_H

// ROOT includes
#include "TMatrixD.h"
#include "TArrayD.h"

// stl includes
#include <vector>

// Track Parameterization is: f = (phi, z, px, py, pz)


class TrackRepresentation{

 public:
  TrackRepresentation();
  ~TrackRepresentation();
  void Extrapolate(int ki, int kf);
  TMatrixD GetMeasurementMatrix(){return m_H;};
  TMatrixD GetCovarianceMatrix(){return m_Cov;};
  void InitTrackRep(double r, double phi, double z, double px, double py, double pz); 
  void UpdateTrackParameters(double phi, double z, double px, double py, double pz); 
  void UpdateCovMatrix(TMatrixD cov){m_Cov = cov;};
  void SetBFieldmap(const char * bmapfile);



  int GetIndex(){return m_k;};
  double GetTime(){return m_t;};
  void GetBField(const double point[3], double *bfield[3]);

 public:
  // the index of the current step
  int m_k; 
  // hypothetical mass
  double m_m;
  // phi, z and r coordinates
  double m_phi, m_z, m_r;
  std::vector<double> m_phi_v;
  std::vector<double> m_z_v;
  std::vector<double> m_r_v;
  // momentum vector
  double m_px, m_py, m_pz;
  std::vector<double> m_px_v, m_py_v, m_pz_v;
  
  // Measurement matrix: 2x5
  TMatrixD m_H;
  // Evolution matrix: 5x5 matrix: \frac{\partial f(phi_k^{k-1}, z_k^{k-1}, px_k^{k-1}, py_k^{k-1}, pz_k^{k-1})}{\partial f(phi_k, z_k, px_k, py_k, pz_k)}
  TMatrixD m_F;
  // Track error/covariance matrix:
  // must be initialized with reasonable large values,
  // and let the Kalman Filter do its evolution
  TMatrixD m_Cov, m_Cov_t;
  
  double m_dt; // The time development step size 
  double m_t; // time
  std::vector<double> m_t_v; // bookkeeping the timestamps
};



#endif
