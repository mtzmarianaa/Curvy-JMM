/* TRIANGLE 2D MESH STRUCTURE

This is the 2d triangle mesh structure. It assumes that an output from 
meshpy is given

*/

#include "triMesh_2D.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

struct triMesh_2D {
    // ask if everything in here is usefull/necessary or if its too much
    // ask if this should be an opaque type (according to me it shouldn't)
    // inspiration: what we might need + what might be useful to plot the mesh using triplot in Python
  coordS *points; // these are ALL the coordinates + number of points in the mesh
  neighborsRS *neighbors; // for each point i, its neighbors (i.e. there is a face that includes both points)
  coordS *boundaryPoints;  // these are just the coordinates of the boundary points + number of boundary points
  facetsS *facets; // these are the "instructions on how to connect indexed dots"
  facesS *faces; // these are the faces, the triangles in the mesh
  int nPoints;
} ;

void triMesh_2Dalloc(triMesh_2Ds **triM_2D) {
    *triM_2D = malloc(sizeof(triMesh_2Ds));
    assert(*triM_2D != NULL);
}

void triMesh_2Ddalloc(triMesh_2Ds **triM_2D) {
    free(*triM_2D);
    *triM_2D = NULL;
}

void triMesh2_init(triMesh_2Ds *triM_2D, coordS *points, neighborsRS *neighbors, coordS *boundaryPoints, facetsS *facets, facesS *faces, int nPoints){
    triM_2D->points = points;
    triM_2D->neighbors = neighbors;
    triM_2D->boundaryPoints = boundaryPoints;
    triM_2D->facets = facets;
    triM_2D->faces = faces;
    triM_2D->nPoints = nPoints;
}

void triMesh2_init_from_meshpy(triMesh_2Ds *triM_2D, char const *pathPoints, char const *pathNeighbors, char const *pathBoundaryPoints, char const *pathFacets, char const *pathFaces){
    // there are a lot of files needed to be opened

    int nPoints;
    nPoints = numLinesInFile(pathPoints); // set number of points
    triM_2D->nPoints = nPoints; // we need nPoints to use neighborsRSalloc_n

    coordS *points; 
    coord_alloc(&points);

    neighborsRS *neighbors;
    neighborsRSalloc_n(&neighbors, nPoints);

    coordS *boundaryPoints;
    coord_alloc(&boundaryPoints);

    facetsS *facets;
    facets_alloc(&facets);

    facesS *faces;
    faces_alloc(&faces);

    // ...
    // Now that everything is alloc'ed and declared we can open the files and set the parameters of triM_2D

    // For the coordS structs
    coord_initFromFile(points, pathPoints);
    triM_2D->points = points;
    coord_initFromFile(boundaryPoints, pathBoundaryPoints);
    triM_2D->boundaryPoints = boundaryPoints;

    // For the neighborsRS struct
    neighbors_init(neighbors, pathNeighbors, nPoints);
    triM_2D->neighbors = neighbors;

    // For the facetsS struct
    facets_initFromFile(facets, pathFacets);
    triM_2D->facets = facets;

    // For the facesS struct
    faces_initFromFile(faces, pathFaces);
    triM_2D->faces = faces;

}

void printGeneralInfoMesh(triMesh_2Ds *triM_2D) {
    printf("\n GENERAL INFORMATION ABOUT THIS MESH \n\n");
    printf("Number of points in the mesh:  %d.\n", triM_2D->nPoints);
    printf("Number of points that conform the boundary:  %d.\n", triM_2D->boundaryPoints->nPoints);
    printf("Number of faces or triangles in the mesh:  %d.\n", triM_2D->faces->nFaces);
}

void printEverythingInMesh(triMesh_2Ds *triM_2D) {
    printf("\n\n---------------------------------------\n");
    printf("---------------------------------------\n");
    printf("\n EVERYTHING CONTAINED IN THIS MESH \n\n");
    printf("\n\n---------------------------------------\n");
    printf("POINTS\n");
    printf("Number of points in the mesh:  %d.\n", triM_2D->nPoints);
    printf("Such points are the following: \n");
    print_coord(triM_2D->points);
    printf("\n\n---------------------------------------\n");
    printf("NEIGHBORS\n");
    printf("The neighbors for each indexed point in this mesh are the following: \n");
    printAllNeighbors(triM_2D->neighbors, triM_2D->nPoints);
    printf("\n\n---------------------------------------\n");
    printf("BOUNDARY POINTS\n");
    printf("Number of points that conform the boundary:  %d.\n", triM_2D->boundaryPoints->nPoints);
    printf("Such boundary points are the following: \n");
    print_coord(triM_2D->boundaryPoints);
    printf("\n\n---------------------------------------\n");
    printf("FACETS\n");
    printf("The facets conforming the boundary are: \n");
    print_facets(triM_2D->facets);
    printf("\n\n---------------------------------------\n");
    printf("FACES\n");
    printf("Number of faces or triangles in the mesh:  %d.\n", triM_2D->faces->nFaces);
    printf("Such faces or triangles are defined as follows: \n");
    print_faces(triM_2D->faces);
}
