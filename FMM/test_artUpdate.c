// testing the artificial triangle update

#include "eik_grid.h"
#include "fmm_2d.h"

#include <math.h>
#include <stdio.h>


int main()
{
    // TEST GEOMETRY FOR DIFFERENT INDICES OF REFRACTION - SOME PATHS GO AROUND, NOT INTO REG3



    const char *pathPoints, *pathNeighbors, *pathIncidentFaces, *pathBoundaryPoints, *pathFacets, *pathFaces, *pathIndexRegions, *pathToSaveTr_, *pathSaveGradientsTr_;
    const char *pathSavePath, *pathSaveLambdas;
    pathPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/TestArtificial/MeshPoints.txt";
    pathNeighbors = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/TestArtificial/Neigh.txt";
    pathIncidentFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/TestArtificial/IncidentFaces.txt";
    pathBoundaryPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/TestArtificial/BoundaryPoints.txt";
    pathFacets = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/TestArtificial/Facets.txt";
    pathFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/TestArtificial/Faces.txt";
    pathIndexRegions = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/TestArtificial/FacesLabel.txt";

    pathToSaveTr_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/TestArtificial/ComputedValues_ARTIFICIAL.bin";
    pathSaveGradientsTr_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/TestArtificial/ComputedGradients_ARTIFICIAL.bin";
    pathSavePath = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/TestArtificial/Parents_ARTIFICIAL.bin";
    pathSaveLambdas = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/TestArtificial/LambdasOpt_ARTIFICIAL.bin";

    int *start;
    int nstart, s;

    s = 0;
    start = &s;
    nstart = 1;
    // now we test the init with just the path to the files

    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n TESTING THE UPDATES WITH ARTIFICIAL TRIANGLES \n\n\n\n");
    eik_gridS *eik_g1;
    eik_grid_alloc(&eik_g1);
    eik_grid_initFromFile(eik_g1, start, nstart, pathPoints, pathNeighbors, pathIncidentFaces, pathBoundaryPoints, pathFacets, pathFaces, pathIndexRegions);

    printAllInfoMesh(eik_g1);

    FMM_2D( eik_g1 , 1);
    

    // saveComputedValues(eik_g1, pathToSaveTr_);

    // saveComputedGradients(eik_g1, pathSaveGradientsTr_);

    // saveComputedParents(eik_g1, pathSavePath);

    // saveComputedLambdas(eik_g1, pathSaveLambdas);

    printGeneralInfo(eik_g1);

    eik_grid_dealloc(&eik_g1);




}