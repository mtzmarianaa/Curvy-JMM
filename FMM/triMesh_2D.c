
/* TRIANGLE 2D MESH STRUCTURE
This is the 2d triangle mesh structure. It assumes that an output from 
meshpy is given
*/

#include "triMesh_2D.h"
#include "linAlg.h"
#include "SoSFunction.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>


void triMesh_2Dalloc(triMesh_2Ds **triM_2D) {
    *triM_2D = malloc(sizeof(triMesh_2Ds));
    assert(*triM_2D != NULL);
}

void triMesh_2Ddalloc(triMesh_2Ds **triM_2D) {
    free(*triM_2D);
    *triM_2D = NULL;
}

void infoTwoPartUpdate_alloc(infoTwoPartUpdate **infoOut) {
    *infoOut = malloc(sizeof(infoTwoPartUpdate));
    assert(*infoOut != NULL);
}

void infoTwoPartUpdate_dalloc(infoTwoPartUpdate **infoOut) {
    free(*infoOut);
    *infoOut = NULL;
}

void triMesh2_init(triMesh_2Ds *triM_2D, coordS *points, neighborsRS *neighbors, neighborsRS *incidentFaces, facesS *faces, int nPoints, int *indexRegions){
    triM_2D->nPoints = nPoints;
    triM_2D->indexRegions = indexRegions;
    triM_2D->points = points;
    triM_2D->neighbors = neighbors;
    triM_2D->incidentFaces = incidentFaces;
    triM_2D->faces = faces;
    triM_2D->indexRegions = indexRegions;
}

void triMesh2_init_from_meshpy(triMesh_2Ds *triM_2D, char const *pathPoints, char const *pathNeighbors, char const *pathIncidentFaces, char const *pathFaces, char const *pathIndexRegions){
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

    facesS *faces;
    faces_alloc(&faces);

    // ...
    // Now that everything is alloc'ed and declared we can open the files and set the parameters of triM_2D

    // For the coordS structs
    coord_initFromFile(points, pathPoints);
    triM_2D->points = points;

    // For the neighborsRS structs
    neighbors_init(neighbors, pathNeighbors, nPoints);
    triM_2D->neighbors = neighbors;
    neighbors_init(incidentFaces, pathIncidentFaces, nPoints);
    triM_2D->incidentFaces = incidentFaces;

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
    region = INFINITY;
    for (i = 0; i<triM_2D->incidentFaces[index_from].len; i++){
        // we first iterate through all the points
        current_face = triM_2D->incidentFaces[index_from].neis_i[i];
        if(triM_2D->faces->points[current_face][0] == index_to |  triM_2D->faces->points[current_face][1] == index_to | triM_2D->faces->points[current_face][2] == index_to) {
            // we've found a triangle that contains both index_from and index_to (there are only 2 triangles that contain both of these points
            region_test = triM_2D->indexRegions[current_face]; // the region of one of the two triangles containing both points
            if(region_test<region){
                region = region_test; // would only happen when we find a face and then if the edge between these two points separates two different areas
            }
        }
    }
    return region;
}



int faceBetween3Points(triMesh_2Ds *triM_2D, int index1, int index2, int index3){
    // given 3 points that share a face it outputs the index of such face if found, if not found it outputs -1
    int faceIndex;
    faceIndex = -1;
    int currentFace;
    for(int i = 0; i<triM_2D->incidentFaces[index1].len; i++){
        currentFace = triM_2D->incidentFaces[index1].neis_i[i];
        if( triM_2D->faces->points[currentFace][0] == index2 | triM_2D->faces->points[currentFace][1] == index2 | triM_2D->faces->points[currentFace][2] == index2  ){
            if(   triM_2D->faces->points[currentFace][0] == index3 | triM_2D->faces->points[currentFace][1] == index3 | triM_2D->faces->points[currentFace][2] == index3   ){
                faceIndex = currentFace;
                break;
            }
        }
    }
    return faceIndex;
}

void twoTrianglesFromEdge(triMesh_2Ds *triM_2D, int index0, int index1, int possibleTriangles[2], int possibleThirdVertices[2]){
    int j, currentTriangle;
    j = 0;
    for(int i = 0; i<triM_2D->incidentFaces[index0].len; i++){  // we loop through the incident faces of x0
        currentTriangle = triM_2D->incidentFaces[index0].neis_i[i]; // current incident face to x0
        if( triM_2D->faces->points[currentTriangle][0] == index0 & triM_2D->faces->points[currentTriangle][1] == index1 ){
            possibleTriangles[j] = currentTriangle;
            possibleThirdVertices[j] = triM_2D->faces->points[currentTriangle][2];
            j++;
        }
        else if ( triM_2D->faces->points[currentTriangle][1] == index0 & triM_2D->faces->points[currentTriangle][0] == index1 ){
            possibleTriangles[j] = currentTriangle;
            possibleThirdVertices[j] = triM_2D->faces->points[currentTriangle][2];
            j++;
        }
        else if ( triM_2D->faces->points[currentTriangle][0] == index0 & triM_2D->faces->points[currentTriangle][2] == index1 ){
            possibleTriangles[j] = currentTriangle;
            possibleThirdVertices[j] = triM_2D->faces->points[currentTriangle][1];
            j++;
        }
        else if ( triM_2D->faces->points[currentTriangle][2] == index0 & triM_2D->faces->points[currentTriangle][0] == index1 ){
            possibleTriangles[j] = currentTriangle;
            possibleThirdVertices[j] = triM_2D->faces->points[currentTriangle][1];
            j++;
        }
        else if ( triM_2D->faces->points[currentTriangle][1] == index0 & triM_2D->faces->points[currentTriangle][2] == index1 ){
            possibleTriangles[j] = currentTriangle;
            possibleThirdVertices[j] = triM_2D->faces->points[currentTriangle][0];
            j++;
        }
        else if ( triM_2D->faces->points[currentTriangle][2] == index0 & triM_2D->faces->points[currentTriangle][1] == index1 ){
            possibleTriangles[j] = currentTriangle;
            possibleThirdVertices[j] = triM_2D->faces->points[currentTriangle][0];
            j++;
        }
        // printf("\n\nOne of the possible triangles is: %d    with third edge %d\n", possibleTriangles[0], possibleThirdVertices[0]);
        // printf("One of the possible triangles is: %d    with third edge %d\n", possibleTriangles[1], possibleThirdVertices[1]);
    }
    if( j == 1){
        // this happens when the point is on an edge, we can't march via two directions
        possibleTriangles[1] = possibleTriangles[0];
        possibleThirdVertices[1] = possibleThirdVertices[0];
    }
}

void pointWhereRegionChanges(triMesh_2Ds *triM_2D, int x0_ind, int x1_ind, int xHat, int directionToStart, infoTwoPartUpdate *infoOut) {
    // given the indices of x0, x1, and xHat we "march" along the triangles with x0 as one of their vertices
    // such that we start at a triangle with x0,x1 as part of their vertices to a triangle that has x0, xHat
    // as part of their vertices. If there is a change in region this function outputs the index of a vertex such
    // that the edge x0 xChange defines the edge of the domain, if no change in region is found, this function outputs -1
    // directionToStart is either 0 or 1, the triangle on which to start marching, since there are two triangles
    // in the mesh that have x0, x1 as one of their edges.
    int previousTriangle, currentTriangle, initialTriangles[2], initialx2[2], x2_ind, previous_x2_ind, iterationTriangles[2], iterationx2[2];
    double x0[2], x2[2], x2_prev[2], xhatC[2];

    x0[0] = triM_2D->points->x[x0_ind];
    x0[1] = triM_2D->points->y[x0_ind];
    x2_prev[0] = triM_2D->points->x[x1_ind];
    x2_prev[1] = triM_2D->points->y[x1_ind];
    xhatC[0] = triM_2D->points->x[xHat];
    xhatC[1] = triM_2D->points->y[xHat];
    
    // printf("The coordinates of x0 are: (   %fl   |   %fl   )\n", x0[0], x0[1] );
    // printf("The coordinates of x1 are: (   %fl   |   %fl   )\n", x2_prev[0], x2_prev[1] );
    // printf("The coordinates of xHat are: (   %fl   |   %fl   )\n", xhatC[0], xhatC[1] );

    infoOut->xChange_ind = -1; // so far no change in region
    infoOut->angle_xHat = 0;
    infoOut->angle_xChange = 0;

    // printf("Try to find the first triangle\n");

    // we get the first triangle we're marching along (depends on directionToStart)
    twoTrianglesFromEdge(triM_2D, x0_ind, x1_ind, initialTriangles, initialx2);

    // define the starting triangle depending on the value of directionToStart
    if(directionToStart == 0){
        currentTriangle = initialTriangles[0];
        x2_ind = initialx2[0];
        // printf("The initial triangle is %d with index of refraction %fl\n", currentTriangle, s_function_threeSections(x0, triM_2D->indexRegions[currentTriangle] ));
        // printf("The initial third point is %d with coordinates   (   %fl   |   %fl   )\n", x2_ind, triM_2D->points->x[x2_ind], triM_2D->points->y[x2_ind]);
        previousTriangle = initialTriangles[1];
    }
    else{
        currentTriangle = initialTriangles[1];
        x2_ind = initialx2[1];
        // printf("The initial triangle is %d with index of refraction %fl\n", currentTriangle, s_function_threeSections(x0, triM_2D->indexRegions[currentTriangle] ));
        // printf("The initial third point is %d with coordinates   (   %fl   |   %fl   )\n", x2_ind, triM_2D->points->x[x2_ind], triM_2D->points->y[x2_ind]);
        previousTriangle = initialTriangles[0];
    }
    // get the coordinates of this x0 and add the angle, initialize the indices of refraction
    x2[0] = triM_2D->points->x[x2_ind];
    x2[1] = triM_2D->points->y[x2_ind];
    infoOut->angle_xHat += angleThreePoints(x2, x0, x2_prev);
    infoOut->angle_xChange += angleThreePoints(x2, x0, x2_prev);
    // printf("The initial indices of refraction are:   %fl    %fl\n", s_function_threeSections(x0, triM_2D->indexRegions[currentTriangle] ), s_function_threeSections(x0, triM_2D->indexRegions[currentTriangle] ));
    infoOut->indexRef_01 = s_function_threeSections(x0, triM_2D->indexRegions[currentTriangle]);
    // printf("%fl\n", infoOut->indexRef_01);
    infoOut->indexRef_02 = s_function_threeSections(x0, triM_2D->indexRegions[currentTriangle]);
    // printf("%fl\n", infoOut->indexRef_02);


    // march along the other triangles
    while( x2_ind != xHat  ) { // while we haven't marched along ALL the possible triangles to go from x0x1 to x0xHat:
        // get the next triangle 
        twoTrianglesFromEdge(triM_2D, x0_ind, x2_ind, iterationTriangles, iterationx2);
        if ( iterationTriangles[0] == currentTriangle ){
            // then the next triangle should be iterationTriangles[1] and the previous triangles now changes to the current one
            previousTriangle = currentTriangle;
            previous_x2_ind = x2_ind;
            currentTriangle = iterationTriangles[1];
            x2_ind = iterationx2[1];
        }
        else {
            // then the next triangle should be iterationTriangles[0] and the previous triangle now changes to the current one
            previousTriangle = currentTriangle;
            previous_x2_ind = x2_ind;
            currentTriangle = iterationTriangles[0];
            x2_ind = iterationx2[0];
        }
        if( triM_2D->indexRegions[currentTriangle] != triM_2D->indexRegions[previousTriangle] ){
            // if this happens then there was a change in region refined by the edge x0 previous_x2
            // this is an opportunity if we want to consider more than one change of region
            infoOut->xChange_ind = previous_x2_ind;
            // we also want to get the different index of refraction
            infoOut->indexRef_02 = s_function_threeSections(x0, triM_2D->indexRegions[currentTriangle]);
        }
        // we also need to update the coordinates of x2 and x2_prev
        x2_prev[0] = x2[0];
        x2_prev[1] = x2[1];
        x2[0] = triM_2D->points->x[x2_ind];
        x2[1] = triM_2D->points->y[x2_ind];
        if (infoOut->xChange_ind == -1){
            // meaning that we haven't found an edge where the region changes
            infoOut->angle_xChange += angleThreePoints(x2, x0, x2_prev);
        }
        infoOut->angle_xHat += angleThreePoints(x2, x0, x2_prev); // but we still add this angle because x2 is not yet xHat
        // printf("In this iteration the previous triangle is %d\n", previousTriangle);
        // printf("The current triangle is %d\n", currentTriangle);
        // printf("The previous x2 is %d with coordinates   (   %fl   |   %fl   )\n", previous_x2_ind, triM_2D->points->x[previous_x2_ind], triM_2D->points->y[previous_x2_ind]);
        // printf("The current x2 is %d with coordinates   (   %fl   |   %fl   )\n", x2_ind, triM_2D->points->x[x2_ind], triM_2D->points->y[x2_ind]);
        // printf("The current angle being considered is: %fl\n", infoOut->angle_xHat);
    }


}