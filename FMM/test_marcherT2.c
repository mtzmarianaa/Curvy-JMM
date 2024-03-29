// testing the marcher with hermite interpolation for T(xLambda)

#include "eik_grid.h"
#include "marcher_T2.h"

#include <math.h>
#include <stdio.h>
#include <time.h>


int main()
{
    // TEST GEOMETRY FOR DIFFERENT INDICES OF REFRACTION - SOME PATHS GO AROUND, NOT INTO REG3



    const char *pathPoints, *pathNeighbors, *pathIncidentFaces, *pathFaces, *pathIndexRegions, *pathToSaveTr_, *pathSaveGradientsTr_;
    const char *pathSavePath, *pathSaveLambdas, *pathTimes;
    pathPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H25/H25_MeshPoints.txt";
    pathNeighbors = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H25/H25_Neigh.txt";
    pathIncidentFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H25/H25_IncidentFaces.txt";
    pathFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H25/H25_Faces.txt";
    pathIndexRegions = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H25/H25_FacesLabel.txt";

    pathToSaveTr_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H25/H25_ComputedValuesCubic.bin";
    pathSaveGradientsTr_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H25/H25_ComputedGradientsCubic.bin";
    pathSavePath = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H25/H25_ParentsCubic.bin";
    pathSaveLambdas = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H25/H25_LambdasOptCubic.bin";
    pathTimes = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H25/H25_TimesCubic.bin";

    int *start;
    int nstart, s;
    double time_taken[1], time;
    clock_t start_t, end_t;

    s = 0;
    start = &s;
    nstart = 1;
    // now we test the init with just the path to the files

    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n TESTING THE UPDATES WITH CUBIC HERMITE INTERPOLATION \n\n\n\n");
    eik_gridS *eik_g1;
    eik_grid_alloc(&eik_g1);
    eik_grid_initFromFile(eik_g1, start, nstart, pathPoints, pathNeighbors, pathIncidentFaces, pathFaces, pathIndexRegions);

    // printAllInfoMesh(eik_g1);

    start_t = clock();
    marcher_T2( eik_g1, 2.5 );
    end_t = clock();
    time = (double)(end_t - start_t)/ CLOCKS_PER_SEC;
    time_taken[0] = time;
    

    saveComputedValues(eik_g1, pathToSaveTr_);

    saveComputedGradients(eik_g1, pathSaveGradientsTr_);

    saveComputedParents(eik_g1, pathSavePath);

    saveComputedLambdas(eik_g1, pathSaveLambdas);

    saveTimes(time_taken, pathTimes);

    printGeneralInfo(eik_g1);

    eik_grid_dealloc(&eik_g1);



}
