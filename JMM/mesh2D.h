#pragma once

#include "coord.h"
#include "faces.h"
#include "neighbors.h"
#include "files_methods.h"


typedef struct {
  // Mesh struct in 2D that also captures information
  // of the boundary (which edges form part of a
  // boundary and the boundary tangent at those points)
  double (*points)[2]; // points on the mesh ALL OF THEM, not only those on the boundaries
  size_t nPoints; // number of points on the mesh
  size_t (*faces)[3]; // indices of the points who form faces SORTED
  size_t nFaces; // number of faces in the mesh
  size_t (*edges)[2]; // pairs of indices of points who form an edge on the mesh (not necesarily on the boundary) SORTED
  size_t nEdges; // number of edges in the mesh
  neighborsRS *neighbors; // for each point i, its neighbors (i.e. there is a face that includes both points)
  neighborsRS *incidentFaces; // for each point i, its incident faces
  double *eta; // indices of refraction on each face
  boundaryCurve *h_i; // boundaryCurve struct for each edge
} mesh2S;


typedef struct {
  // Boundary curve struct in 2D
  size_t i_Edge; // to which edge these tangents belong to
  double B[2][2]; // B1 and B2, tangents at the end point of the edge i
} boundaryCurve;


void mesh2_alloc(mesh2S **mesh2);

void mesh2_dealloc(mesh2S **mesh2);

void boundaryCurve_alloc(boundaryCurve **h_i);

void boundaryCurve_dealloc(mesh2S **h_i);

void boundaryCurve_init(boundaryCurve *h_i, size_t i_Edge, double B[2][2]);

void mesh2_init(mesh2S *mesh2, double (*points)[2], size_t nPoints, size_t (*faces)[3], size_t nFaces, size_t (*edges)[2], size_t nEdges, double *eta, boundaryCurve *h_i);





