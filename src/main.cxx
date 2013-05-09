#include <iostream>

// ROOT includes
#include <TApplication.h>

// local includes
#include "TrackRepresentation.h"
#include "KalmanFilter.h"

using namespace std;

int main(int argc, char * argv[]){


  TrackRepresentation t;
  //  t.SetBFieldmap("Brz_fieldmap.txt");
  KalmanFilter k();

  cout << "Hello World" << std::endl;

  return 0;
}
