//fieldmap.h
/*Still todo. 
        *Build in 2D axial symmetry in Comsol
        *Build in Comsol with units in mm, as it should be, so no conversion should be done here.
        *Make the MatrixSize dynamic, or calculate the stepsize at least!
        *Think about out of bound he, Check if it`s bigger as trapparameters orso?? (via_r_trap best he...)
*/

#ifndef FIELDMAP_H
#define FIELDMAP_H


#include <iostream>
using namespace std;

#include <fstream>
using std::ifstream;
using std::ofstream;

#include <cstdlib> // for exit function
#include <sstream>
#include "math.h"
#include <stdexcept>
#include "matrix.h"
#include <map>
#include <vector>

Matrix<pair<double, double> > * ReadElectricFieldAscii(const char * Erz_filename);
Matrix<pair<double, double> > * ReadElectricFieldAscii(const char *_Er_filename,const char *_Ez_filename);
Matrix<pair<double, double> > * ReadMagneticFieldAscii(const char *_B_filename);


vector<double> getElectrAscii(const double& _x,const double& _y,const double& _z);
vector<double> getMagnAscii(const double& _x,const double& _y,const double& _z);

double get_Ez_min();
double get_Ez_max();
double get_Ez_step();
double get_Er_min();
double get_Er_max();
double ger_Er_step();

double get_Bz_min();
double get_Bz_max();
double get_Bz_step();
double get_Br_min();
double get_Br_max();
double ger_Br_step();


#endif
