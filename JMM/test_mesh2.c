#include "mesh2D.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main(){

    mesh2S *mesh2;
    mesh2_alloc(&mesh2);
    const char *pathPoints, *pathFaces, *pathEdges, *pathNeighbors, *pathIncidentFaces, *pathIndices, *pathBoundary;
    
    pathPoints = "/Users/marianamartinez/Documents/Curvy-JMM/TestBaseSnow//H0/H0_MeshPoints.txt";
    pathFaces =  "/Users/marianamartinez/Documents/Curvy-JMM/TestBaseSnow/H0/H0_Faces.txt";
    pathEdges =  "/Users/marianamartinez/Documents/Curvy-JMM/TestBaseSnow/H0/H0_Edges.txt";
    pathNeighbors =  "/Users/marianamartinez/Documents/Curvy-JMM/TestBaseSnow/H0/H0_Neigh.txt";
    pathIncidentFaces =  "/Users/marianamartinez/Documents/Curvy-JMM/TestBaseSnow/H0/H0_IncidentFaces.txt";
    pathIndices =  "/Users/marianamartinez/Documents/Curvy-JMM/TestBaseSnow/H0/H0_Indices.txt";
    pathIndices =  "/Users/marianamartinez/Documents/Curvy-JMM/TestBaseSnow/H0/H0_Indices.txt";
    pathBoundary =  "/Users/marianamartinez/Documents/Curvy-JMM/TestBaseSnow/H0/H0_BoundaryCurve.txt";

    
    mesh2_init_from_meshpy(mesh2, pathPoints, pathFaces, pathEdges, pathNeighbors, pathIncidentFaces, pathIndices, pathBoundary);

    printf("GENERAL INFO \n\n");
    printGeneralInfoMesh(mesh2);

    printf("\n\n\n\nALL INFO \n\n");
    printEverythingInMesh(mesh2);

    mesh2_dealloc(&mesh2);


}
