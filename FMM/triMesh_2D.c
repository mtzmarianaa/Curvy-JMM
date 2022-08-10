
/* TRIANGLE 2D MESH STRUCTURE
This is the 2d triangle mesh structure. It assumes that an output from 
meshpy is given
*/

#include "triMesh_2D.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


void triMesh_2Dalloc(triMesh_2Ds **triM_2D) {
    *triM_2D = malloc(sizeof(triMesh_2Ds));
    assert(*triM_2D != NULL);
}

void triMesh_2Ddalloc(triMesh_2Ds **triM_2D) {
    free(*triM_2D);
    *triM_2D = NULL;
}

void triMesh2_init(triMesh_2Ds *triM_2D, coordS *points, neighborsRS *neighbors, neighborsRS *incidentFaces, coordS *boundaryPoints, facetsS *facets, facesS *faces, int nPoints, int *indexRegions){
    triM_2D->nPoints = nPoints;
    triM_2D->indexRegions = indexRegions;
    triM_2D->points = points;
    triM_2D->neighbors = neighbors;
    triM_2D->incidentFaces = incidentFaces;
    triM_2D->boundaryPoints = boundaryPoints;
    triM_2D->facets = facets;
    triM_2D->faces = faces;
    triM_2D->indexRegions = indexRegions;
}

void triMesh2_init_from_meshpy(triMesh_2Ds *triM_2D, char const *pathPoints, char const *pathNeighbors, char const *pathIncidentFaces, char const *pathBoundaryPoints, char const *pathFacets, char const *pathFaces, char const *pathIndexRegions){
    // there are a lot of files needed to be opened

    int nPoints;
    nPoints = numLinesInFile(pathPoints); // set number of points
    triM_2D->nPoints = nPoints; // we need nPoints to use neighborsRSalloc_n

    coordS *points; 
    coord_alloc(&points);

    neighborsRS *neighbors;
    neighborsRSalloc_n(&neighbors, nPoints);

    neighborsRS *incidentFaces;
    neighborsRSalloc_n(&incidentFaces, nPoints);

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

    // For the neighborsRS structs
    neighbors_init(neighbors, pathNeighbors, nPoints);
    triM_2D->neighbors = neighbors;
    neighbors_init(incidentFaces, pathIncidentFaces, nPoints);
    triM_2D->incidentFaces = incidentFaces;

    // For the facetsS struct
    facets_initFromFile(facets, pathFacets);
    triM_2D->facets = facets;

    // For the facesS struct
    faces_initFromFile(faces, pathFaces);
    triM_2D->faces = faces;

    // for the index Regions
    int *indexRegions;
    indexRegions = malloc(faces->nFaces*sizeof(int));
    readIntColumn(pathIndexRegions, indexRegions);
    triM_2D->indexRegions = indexRegions;

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
    printf("INCIDENT FACES\n");
    printf("The incident faces for each indexed point in this mesh: \n");
    printAllNeighbors(triM_2D->incidentFaces, triM_2D->nPoints);
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
    printf("\n\n---------------------------------------\n");
    printf("FACES INDICES\n");
    for(int i = 0; i<triM_2D->faces->nFaces; i++){
        printf("Face %d, index %d\n", i, triM_2D->indexRegions[i]);
    }
}

int regionBetweenTwoPoints(triMesh_2Ds *triM_2D, int index_from, int index_to){
    // we look in the incidentFaces, both points must share two faces, we are looking for this and for the smaller region (meaning that the index of refraction is smaller there)
    int current_face, i, region, region_test;
    region = triM_2D->indexRegions[ triM_2D->incidentFaces[index_from].neis_i[0] ]; // index corresponding to the first incident face of from
    for (i = 1; i<triM_2D->incidentFaces[index_from].len; i++){
        current_face = triM_2D->incidentFaces[index_from].neis_i[i];
        region_test = triM_2D->indexRegions[ current_face ];
        if ( region != region_test & region_test < region ){
            // if there is an incident face with different index region and that index is smaller (if this happens then we're on the edge and we assume that 
            // we can travel "fast" via the edge, check the idea behind computing the analytic solution)
            if (  triM_2D->faces->points[current_face][0] == index_to |  triM_2D->faces->points[current_face][1] == index_to | triM_2D->faces->points[current_face][2] == index_to ){
                //printf("From %d to %d we consider the face indexed with %d \n", index_from, index_to, current_face);
                region = region_test;
            }
        }
    }
    return region;
}

int faceBetween3Points(triMesh_2Ds *triM_2D, int index1, int index2, int index3){
    // given 3 points that share a face it outputs the index of such face
    int faceIndex;
    faceIndex = triM_2D->faces->nFaces;
    int currentFace;
    for(int i = 0; i<triM_2D->incidentFaces[index1].len; i++){
        currentFace = triM_2D->incidentFaces[index1].neis_i[i];
        if( triM_2D->faces->points[currentFace][0] == index2 | triM_2D->faces->points[currentFace][1] == index2 | triM_2D->faces->points[currentFace][2] == index2  ){
            if(   triM_2D->faces->points[currentFace][0] == index3 | triM_2D->faces->points[currentFace][1] == index3 | triM_2D->faces->points[currentFace][2] == index3   ){
                faceIndex = currentFace;
                //printf("Face between points %d  %d  %d  is  %d", index1, index2, index3, faceIndex);
                break;
            }
        }
    }
    assert(faceIndex != triM_2D->faces->nFaces);
    return faceIndex;
}