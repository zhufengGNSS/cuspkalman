#include "KalmanFilter.h"


KalmanFilter::KalmanFilter(){

}

KalmanFilter::~KalmanFilter(){

}

void KalmanFilter::AddHit(double r, double phi, double z, double re, double phie, double ze){
  m_r.push_back(r);m_re.push_back(re);
  m_phi.push_back(phi);m_phie.push_back(phie);
  m_z.push_back(z);m_ze.push_back(ze);

}

void Process();

void CalcChi2();


