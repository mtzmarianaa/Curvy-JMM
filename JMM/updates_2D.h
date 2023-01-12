#pragma once

#include "triMesh_2D.h"


// we need to define a struct with the information needed
// for the input and output, this will make things easier

typedef struct info_update {
  int indexAccepted; // this is the index of x0
  int x1_ind; // this is the index of x1
  int x2_ind; // this is the index of x2 (if same as x1 means that the update is a two point update
  int xHat_ind; // index of xHat, the point we want to update
  double T0; // value of the eikonal at x0 WHICH IS A VALID NODE
  double T1; // value of the eikonal at x1 WHICH IS A VALID NODE
}


