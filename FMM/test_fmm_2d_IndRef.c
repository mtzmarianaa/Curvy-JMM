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
    const char *pathSavePath, *pathSaveLambdas;
    pathPoints_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_MeshPoints.txt";
    pathNeighbors_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_Neigh.txt";
    pathIncidentFaces_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_IncidentFaces.txt";
    pathBoundaryPoints_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_BoundaryPoints.txt";
    pathFacets_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_Facets.txt";
    pathFaces_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_Faces.txt";
    pathIndexRegions_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_FacesLabel.txt";

    pathToSaveTr_H1_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_ComputedValues.bin";
    pathSaveGradientsTr_H1_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_ComputedGradients.bin";
    pathSavePath = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_Parents.bin";
    pathSaveLambdas = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_LambdasOpt.bin";

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

    saveComputedParents(eik_g1, pathSavePath);

    saveComputedLambdas(eik_g1, pathSaveLambdas);

    printGeneralInfo(eik_g1);

    eik_grid_dealloc(&eik_g1);




}