// testing the marcher with hermite interpolation for T(xLambda)

#include "eik_grid.h"
#include "updates_2D.h"

#include <math.h>
#include <stdio.h>
#include <time.h>


int main()
{
    // TEST GEOMETRY FOR DIFFERENT INDICES OF REFRACTION - SOME PATHS GO AROUND, NOT INTO REG3



  const char *pathPoints, *pathNeighbors, *pathIncidentFaces, *pathFaces, *pathIndexRegions, *pathToSaveTr_, *pathSaveGradientsTr_, *path_BoundaryTan, *path_BoundaryChain;
    const char *pathSavePath, *pathSaveLambdas, *pathTimes;
    pathPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H0/H0_MeshPoints.txt";
    pathNeighbors = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H0/H0_Neigh.txt";
    pathIncidentFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H0/H0_IncidentFaces.txt";
    pathFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H0/H0_Faces.txt";
    pathIndexRegions = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H0/H0_FacesLabel.txt";

    pathToSaveTr_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H0/H0_ComputedValuesCubic.bin";
    pathSaveGradientsTr_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H0/H0_ComputedGradientsCubic.bin";
    path_BoundaryTan = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H0/H0_Boundary_tan.txt";
    path_BoundaryChain = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H0/H0_Boundary_chain.txt";
    pathSavePath = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H0/H0_ParentsCubic.bin";
    pathSaveLambdas = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H0/H0_LambdasOptCubic.bin";
    pathTimes = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H0/H0_TimesCubic.bin";

    int *start;
    int nstart, s;

    s = 0;
    start = &s;
    nstart = 1;
    // now we test the init with just the path to the files

    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n TESTING THE UPDATES WITH CURVY TRIANGLES \n\n\n\n");
    eik_gridS *eik_g1;
    eik_grid_alloc(&eik_g1);
    eik_grid_initFromFile(eik_g1, start, nstart, pathPoints, pathNeighbors, pathIncidentFaces, pathFaces, pathIndexRegions, path_BoundaryTan, path_BoundaryChain);

    printAllInfoMesh(eik_g1);


    info_updateS *info_update;
    info_update_alloc(&info_update);

    int indexAccepted, x1_ind, x2_ind, xHat_ind;
    double T0, T1, indexRef_01, indexRef_02;

    indexAccepted = 48;
    x1_ind = 48;
    x2_ind = 48;
    xHat_ind = 49;
    T0 = 14.8753;
    T1 = -1;
    indexRef_01 = regionBetweenTwoPoints(eik_g1->triM_2D, indexAccepted, x1_ind);
    indexRef_02 = regionBetweenTwoPoints(eik_g1->triM_2D, indexAccepted, x2_ind);

    info_update_initCr(info_update, indexAccepted, xHat_ind, T0, indexRef_01);
    
    print_info_update(info_update);

    creepingUpdate(eik_g1->triM_2D, info_update);

    print_info_update(info_update);


    eik_grid_dealloc(&eik_g1);



}
