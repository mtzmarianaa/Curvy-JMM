#include "mesh2D.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main(){

    mesh2S *mesh2;
    mesh2_alloc(&mesh2);
    const char *pathPoints, *pathFaces, *pathEdges, *pathEdgesInFace, *pathNeighbors;
    const char *pathIncidentFaces, *pathIndices, *pathBoundary;
    
    pathPoints = "/Users/marianamartinez/Documents/Curvy-JMM/TestBaseSnow//H0/H0_MeshPoints.txt";
    pathFaces =  "/Users/marianamartinez/Documents/Curvy-JMM/TestBaseSnow/H0/H0_Faces.txt";
    pathEdges =  "/Users/marianamartinez/Documents/Curvy-JMM/TestBaseSnow/H0/H0_Edges.txt";
    pathEdgesInFace = "/Users/marianamartinez/Documents/Curvy-JMM/TestBaseSnow/H0/H0_EdgesInFace.txt";
    pathNeighbors =  "/Users/marianamartinez/Documents/Curvy-JMM/TestBaseSnow/H0/H0_Neigh.txt";
    pathIncidentFaces =  "/Users/marianamartinez/Documents/Curvy-JMM/TestBaseSnow/H0/H0_IncidentFaces.txt";
    pathIndices =  "/Users/marianamartinez/Documents/Curvy-JMM/TestBaseSnow/H0/H0_Indices.txt";
    pathIndices =  "/Users/marianamartinez/Documents/Curvy-JMM/TestBaseSnow/H0/H0_Indices.txt";
    pathBoundary =  "/Users/marianamartinez/Documents/Curvy-JMM/TestBaseSnow/H0/H0_BoundaryCurve.txt";

    
    mesh2_init_from_meshpy(mesh2, pathPoints, pathFaces, pathEdges, pathEdgesInFace,
			   pathNeighbors, pathIncidentFaces, pathIndices, pathBoundary);

    /* printf("GENERAL INFO \n\n"); */
    /* printGeneralInfoMesh(mesh2); */

    /* printf("\n\n\n\nALL INFO \n\n"); */
    /* printEverythingInMesh(mesh2); */


    // Now we test the triangleFan

    triangleFanS *triFan;
    triangleFan_alloc(&triFan);
    size_t nRegions, index0, index1, indexHat, *listIndicesNodes;

    // going to use:
    //      triangle 247 with points 89, 350, 513
    //      triangle 1017 with points 89, 90, 513

    nRegions = 2;
    index0 = 89; // triangle 247 with points 89, 350, 513
    index1 = 350;
    indexHat = 90;

    listIndicesNodes = malloc(4*sizeof(size_t));
    listIndicesNodes[0] = 89;
    listIndicesNodes[1] = 350;
    listIndicesNodes[2] = 513;
    listIndicesNodes[3] = 90;

    triangleFan_initFromIndices(triFan, mesh2, nRegions, index0, index1, indexHat, listIndicesNodes);

    printEverythingTriFan(triFan);
    

    

    mesh2_dealloc(&mesh2);


}
