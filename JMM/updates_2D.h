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
  double grad0[2]; // value of the gradient of the eikonal at x0
  double T1; // value of the eikonal at x1 WHICH IS A VALID NODE
  double grad1[2]; // value of the gradient of the eikonal at x1
  double THat; // proposef value of the eikonal at xHat WHICH IS A TRIAL NODE
  double indexRef_01; // index of refraction first section
  double indexRef_02; // index of refraction second section, might be the same as indexRef_01
  double lambda; // optimal lambda found, parametrization from x0 to x1
  double mu; // optimal mu found, parametrization from x0 to x2
} info_updateS;

void info_update_alloc(info_updateS **info_update);

void info_update_dealloc(info_updateS **info_update);

void info_update_init(info_updateS *info_update, int indexAccepted, int x1_ind, int x2_ind, int xHat_ind, double T0, double grad0[2], double T1, double grad1[2], double indexRef_01, double indexRef_02);

void info_update_initCr(info_updateS *info_update, int indexAccepted, int xHat_ind, double T0, double indexMin);

void info_update_initTwo(info_updateS *info_update, int indexAccepted, int x1_ind, int xHat_ind, double T0, double grad0[2], double T1, double grad1[2], double indexRef);

void print_info_update(info_updateS *info_update);

void creepingUpdate(triMesh_2Ds *triM_2D, info_updateS *info_update);

void simple_TwoPointUpdate(triMesh_2Ds *triM_2D, info_updateS *info_update);

void fromBoundary_TwoPointUpdate(triMesh_2Ds *triM_2D, info_updateS *info_update);

void anchorHatonBoundary_freeSpaceUpdate(triMesh_2Ds *triM_2D, info_updateS *info_update, int anchorOnBoundary);


