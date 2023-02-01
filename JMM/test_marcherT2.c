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
    const char *pathSavePath, *pathSaveLambdas, *pathTimes, *pathMus, *pathTypes, *pathBoundaryTan, *pathBoundaryChain;
    pathPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H14/H14_MeshPoints.txt";
    pathNeighbors = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H14/H14_Neigh.txt";
    pathIncidentFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H14/H14_IncidentFaces.txt";
    pathFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H14/H14_Faces.txt";
    pathIndexRegions = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H14/H14_FacesLabel.txt";
    pathBoundaryTan = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H14/H14_Boundary_tan.txt";
    pathBoundaryChain = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H14/H14_Boundary_chain.txt";

    pathToSaveTr_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H14/H14_ComputedValues.bin";
    pathSaveGradientsTr_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H14/H14_ComputedGradients.bin";
    pathSavePath = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H14/H14_Parents.bin";
    pathSaveLambdas = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H14/H14_Lambdas.bin";
    pathMus = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H14/H14_Mus.bin";
    pathTimes = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H14/H14_Times.bin";
    pathTypes = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H14/H14_Types.bin";

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
    eik_grid_initFromFile(eik_g1, start, nstart, pathPoints, pathNeighbors, pathIncidentFaces, pathFaces, pathIndexRegions, pathBoundaryTan, pathBoundaryChain);

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

    saveComputedMus(eik_g1, pathMus);

    saveComputedTypesUpdate(eik_g1, pathTypes);

    saveTimes(time_taken, pathTimes);

    printGeneralInfo(eik_g1);

    eik_grid_dealloc(&eik_g1);



}
