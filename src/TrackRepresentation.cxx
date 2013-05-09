#include "TrackRepresentation.h"

#include "fieldmap.h"

TrackRepresentation::TrackRepresentation(){

  // initialize 
  m_px = m_py = m_pz = m_phi = m_z = m_r = m_t = 0;
  m_dt = 0.01;

}

TrackRepresentation::~TrackRepresentation(){
}

// To 
void TrackRepresentation::Extrapolate(int ki, int kf){

  //  m_px = m_px_v[ki] + 

  // const double point[3] = {0.01, 0.01, -0.01};
  // double * bfield = new double[3];
  // GetBField(&point[0], &bfield);



}

void TrackRepresentation::InitTrackRep(double r, double phi, double z, double px, double py, double pz){

  m_px = px;m_py = py; m_pz = pz; m_phi = phi; m_z = z;
  m_px_v.push_back(m_px);
  m_py_v.push_back(m_py);
  m_pz_v.push_back(m_pz);
  m_phi_v.push_back(m_phi);
  

}
void TrackRepresentation::UpdateTrackParameters(double phi, double z, double px, double py, double pz){}


void TrackRepresentation::SetBFieldmap(const char * bmapfile){
  
  cout << "Reading CUSP magnetic field map from ASCII file" << endl;
  ReadMagneticFieldAscii(bmapfile);

}

// input point in meters, output in Tesla
void TrackRepresentation::GetBField(const double point[3], double *bfield[3]){

  std::vector<double> B_v;
  B_v = getMagnAscii(point[0],point[1],point[2]);
  (*bfield)[0] = B_v[0];
  (*bfield)[1] = B_v[1];
  (*bfield)[2] = B_v[2];
  cout << "Point: " << point[0] << ", " << point[1] << ", " << point[2] << endl;
  cout << "\t Bfield is " << (*bfield)[0] << ", " << (*bfield)[1] << ", " << (*bfield)[2] << endl;


}
