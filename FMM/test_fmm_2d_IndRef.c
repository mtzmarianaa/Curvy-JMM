/*
TESTS FOR THE GRID TO TEST THE DIFFERENCE ON INDICES OF REFRACTION
*/
#include "eik_grid.h"
#include "fmm_2d.h"

#include <math.h>
#include <stdio.h>


int main()
{
    // // // TEST GEOMETRY FOR DIFFERENT INDICES OF REFRACTION - SOME PATHS GO AROUND, NOT INTO REG3



    const char *pathPoints_H1, *pathNeighbors_H1, *pathIncidentFaces_H1, *pathBoundaryPoints_H1, *pathFacets_H1, *pathFaces_H1, *pathIndexRegions_H1, *pathToSaveTr_H1_, *pathSaveGradientsTr_H1_;
    const char *pathSavePath;
    pathPoints_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/TestIndex/MeshPoints.txt";
    pathNeighbors_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/TestIndex/Neigh.txt";
    pathIncidentFaces_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/TestIndex/IncidentFaces.txt";
    pathBoundaryPoints_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/TestIndex/BoundaryPoints.txt";
    pathFacets_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/TestIndex/Facets.txt";
    pathFaces_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/TestIndex/Faces.txt";
    pathIndexRegions_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/TestIndex/FacesLabel.txt";

    pathToSaveTr_H1_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/TestIndex/ComputedValues_10.bin";
    pathSaveGradientsTr_H1_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/TestIndex/ComputedGradients_10.bin";
    pathSavePath = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/TestIndex/Paths_10.bin";

    int *start_H1;
    int nstart_H1, s_H1;

    s_H1 = 0;
    start_H1 = &s_H1;
    nstart_H1 = 1;
    // now we test the init with just the path to the files

    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n TESTING FROM FILES OF A TRIANGULAR MESH FOR TEST GEOMERTY difference indices of refraction - THREE SECTIONS \n\n\n\n");
    eik_gridS *eik_g1;
    eik_grid_alloc(&eik_g1);
    eik_grid_initFromFile(eik_g1, start_H1, nstart_H1, pathPoints_H1, pathNeighbors_H1, pathIncidentFaces_H1, pathBoundaryPoints_H1, pathFacets_H1, pathFaces_H1, pathIndexRegions_H1);

    FMM_2D( eik_g1 );
    

    saveComputedValues(eik_g1, pathToSaveTr_H1_);

    saveComputedGradients(eik_g1, pathSaveGradientsTr_H1_);

    savePathsTaken(eik_g1, pathSavePath);

    printGeneralInfo(eik_g1);

    eik_grid_dealloc(&eik_g1);




}