/* CAUSTICS
   This is the caustics struct, how to refine the mesh and numerically find caustics given that an eikonal solver has been implemented

*/

#include "eik_grid.h"
#include "opti_method.h"
#include "linAlg.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

void caustics_alloc(caustics_s **caus) {
  *caus = malloc(sizeof(caustics_s));
  assert(*caus != NULL);
}

void caustics_dealloc(caustics_s **caus) {
  free(*caus);
  *caus = NULL;
}

void caustics_initFromEik(caustics_s *caus, eik_gridS *eik_g){
  // from running an eikonal solver on a triangle mesh we can initiate the caustics object
}
