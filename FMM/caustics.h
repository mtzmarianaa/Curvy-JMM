#pragma once

typedef struct caustics {
  eik_gridS *eik_g; // where we are going to get the information from, an already solved eikonal in a triangle mesh
} caustics_s;
