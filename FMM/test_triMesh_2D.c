#include "coord.h"
#include "triMesh_2D.h"
#include "facets.h"
#include "faces.h"
#include "neighbors.h"
#include "files_methods.h"

#include <stdio.h>
#include <stdlib.h>


int main(){

    // NAIVE TESTING
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n NAIVE TESTING, BUILDING EVERYTHING AND SETTING EACH PARAMETER IN THE MESH \n\n\n\n");

    // creating everything
    // points

    coordS *points;
    coord_alloc(&points);
    const char *pathPoints;
    pathPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/MeshPoints.txt";
    coord_initFromFile(points, pathPoints);
    // neighbors
    int nPoints;
    char const *pathNeighbors;
    pathNeighbors = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/Neigh.txt";
    nPoints = numLinesInFile(pathNeighbors);
    neighborsRS *neighbors;
    neighborsRSalloc_n(&neighbors, nPoints);
    neighbors_init(neighbors, pathNeighbors, nPoints);
    //boundaryPoints
    coordS *boundaryPoints;
    coord_alloc(&boundaryPoints);
    const char *pathBoundaryPoints;
    pathBoundaryPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/BoundaryPoints.txt";
    coord_initFromFile(boundaryPoints, pathBoundaryPoints);
    // facets
    facetsS *facets;
    facets_alloc(&facets);
    const char *pathFacets;
    pathFacets = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/Facets.txt";
    facets_initFromFile(facets, pathFacets);
    // faces
    facesS *faces;
    faces_alloc(&faces);
    const char *pathFaces;
    pathFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/Faces.txt";
    faces_initFromFile(faces, pathFaces);
    
    // init the triMesh_2Ds with all these structs we've already built
    triMesh_2Ds *triM_2D_1;
    triMesh_2Dalloc(&triM_2D_1);
    triMesh2_init(triM_2D_1, points, neighbors, boundaryPoints, facets, faces, nPoints);
    // print to see what's in there
    printGeneralInfoMesh(triM_2D_1);
    // print EVERYTHING
    printEverythingInMesh(triM_2D_1);


    // TESTING FROM FILE
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n TESTING FROM FILES \n\n\n\n");
    triMesh_2Ds *triM_2D_2;
    triMesh_2Dalloc(&triM_2D_2);
    triMesh2_init_from_meshpy(triM_2D_2, pathPoints, pathNeighbors, pathBoundaryPoints, pathFacets, pathFaces);
    // print to see what's in there
    printGeneralInfoMesh(triM_2D_2);
    // print EVERYTHING
    printEverythingInMesh(triM_2D_2);

    // TESTING FROM FILES BUT IN THE MESH THAT ALSO INCLUDES A RECTANGLE
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n TESTING FROM FILES WITH RECTANGLE INCLUDED IN MESH \n\n\n\n");
    triMesh_2Ds *triM_2D_3;
    triMesh_2Dalloc(&triM_2D_3);
    const char *pathPoints_sq, *pathNeighbors_sq, *pathBoundaryPoints_sq, *pathFacets_sq, *pathFaces_sq;
    pathPoints_sq = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/MeshPoints_Sq.txt";
    pathNeighbors_sq = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/Neigh_Sq.txt";
    pathBoundaryPoints_sq = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/BoundaryPoints_Sq.txt";
    pathFacets_sq = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/Facets_Sq.txt";
    pathFaces_sq = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/Faces_Sq.txt";
    triMesh2_init_from_meshpy(triM_2D_3, pathPoints_sq, pathNeighbors_sq, pathBoundaryPoints_sq, pathFacets_sq, pathFaces_sq);
    // print to see what's in there
    printGeneralInfoMesh(triM_2D_3);
    // print EVERYTHING
    printEverythingInMesh(triM_2D_3);



}