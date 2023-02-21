// testing the marcher with hermite interpolation for T(xLambda)

#include "eik_grid.h"
#include "updates_2D.h"
#include "opti_method.h"

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

    /* int indexAccepted, x1_ind, x2_ind, xHat_ind; */
    /* double T0, T1, indexRef_01, indexRef_02; */

    /* indexAccepted = 48; */
    /* xHat_ind = 49; */
    /* T0 = 14.8753; */
    /* indexRef_01 = regionBetweenTwoPoints(eik_g1->triM_2D, indexAccepted, x1_ind); */

    /* printf("\n\nTESTING CREEPING UPDATE\n\n"); */

    /* info_update_initCr(info_update, indexAccepted, xHat_ind, T0, indexRef_01); */
    
    /* print_info_update(info_update); */

    /* creepingUpdate(eik_g1->triM_2D, info_update); */

    /* print_info_update(info_update); */

    /* printf("\n\nTESTING FREE SPACE TWO POINT UPDATE\n\n"); */

    /* int indexAccepted, x1_ind, xHat_ind; */
    /* double T0, T1, grad0[2], grad1[2], indexRef; */

    /* grad0[0] = -0.228191; */
    /* grad0[1] = -0.973616; */
    /* grad1[0] = -0.409301; */
    /* grad1[1] = -0.9124; */
    /* indexAccepted = 3128; */
    /* x1_ind = 2942; */
    /* xHat_ind = 127; */
    /* T0 = 8.21679; */
    /* T1 = 7.32957; */
    /* indexRef = 1.0; */

    /* printf("T1: %lf\n", T1); */
    /* info_update_initTwo(info_update, indexAccepted, x1_ind, xHat_ind, T0, grad0, T1, grad1, indexRef); */

    /* print_info_update(info_update); */

    /* simple_TwoPointUpdate(eik_g1->triM_2D, info_update); */

    /* print_info_update(info_update); */


    /* printf("\n\nTESTING UPDATE FROM EDGE\n\n"); */

    /* double lambda0, T0, grad0[2], B0[2], T1, grad1[2], B1[2], x0[2], x1[2], xHat[2], indexRef, lam_opti, fObj, der_Obj; */

    /* lambda0 = 0.056389; */
    /* T0 = 0.70710678; */
    /* T1 = 0.70710678; */
    /* x0[0] = -1; */
    /* x0[1] = 1; */
    /* grad0[0] = -0.70710678; */
    /* grad0[1] = 0.70710678; */
    /* x1[0] = 1; */
    /* x1[1] = 1; */
    /* grad1[0] = 0.70710678; */
    /* grad1[1] = 0.70710678; */
    /* B0[0] = 0.89442719; */
    /* B0[1] = 0.4472136; */
    /* B1[0] = 0.89442719; */
    /* B1[1] = -0.4472136; */
    /* xHat[0] = 0; */
    /* xHat[1] = 5; */
    /* indexRef = 1; */
    

    /* fObj = fobjective_fromEdge(0, T0, grad0, B0, T1, grad1, B1, x0, x1, xHat, indexRef); */
    /* der_Obj = der_fromEdge(0, T0, grad0, B0, T1, grad1, B1, x0, x1, xHat, indexRef); */

    /* printf("Starting with objetive value %lf    and derivative  %lf\n", fObj, der_Obj); */

    /* lam_opti = projectedGradient_fromEdge(lambda0, 0, 1,T0, grad0, B0, T1, grad1, B1, x0, x1, xHat, 0.0001, 50, indexRef); */

    /* printf("Optimum lambda %lf\n", lam_opti); */
    /* fObj = fobjective_fromEdge(lam_opti, T0, grad0, B0, T1, grad1, B1, x0, x1, xHat, indexRef); */
    /* der_Obj = der_fromEdge(lam_opti, T0, grad0, B0, T1, grad1, B1, x0, x1, xHat, indexRef); */
    /* printf("Final objective value %lf   and derivative   %lf\n",  fObj, der_Obj); */
    
    /* printf("\n\nTESTING FREE SPACE OPTIMIZATION - PROJECTED GRADIENT DESCENT\n\n"); */
    
    /* double lambda0, lambdaMin, lambdaMax, TA, gradA[2], TB, gradB[2], xA[2], xB[2], xHat[2], tol, maxIter, indexRef, lamOpti; */
    /* double fObj, derObj; */
    
    /* lambda0 = 1.0; */
    /* lambdaMin = 0.6; */
    /* lambdaMax = 1.0; */
    /* TA = 0.70710678; */
    /* gradA[0] = -0.70710678; */
    /* gradA[1] = 0.70710678; */
    /* TB = 0.70710678; */
    /* gradB[0] = 0.70710678; */
    /* gradB[1] = 0.70710678; */
    /* xA[0] = -1; */
    /* xA[1] = 1; */
    /* xB[0] = 1; */
    /* xB[1] = 1; */
    /* xHat[0] = 0; */
    /* xHat[1] = 2; */
    /* tol = 0.0001; */
    /* maxIter = 50; */
    /* indexRef = 1.0; */

    

    /* fObj = fobjective_freeSpace(lambda0, TA, gradA, TB, gradB, xA, xB, xHat, indexRef); */
    /* derObj = der_freeSpace(lambda0, TA, gradA, TB, gradB, xA, xB, xHat, indexRef); */

    /* printf("Starting with objective value  %lf   and derivative  %lf\n", fObj, derObj); */

    /* lamOpti = projectedGradient_freeSpace(lambda0, lambdaMin, lambdaMax, TA, gradA, TB, gradB, xA, xB, xHat, tol, maxIter, indexRef); */

    /* fObj = fobjective_freeSpace(lamOpti, TA, gradA, TB, gradB, xA, xB, xHat, indexRef); */
    /* derObj = der_freeSpace(lamOpti, TA, gradA, TB, gradB, xA, xB, xHat, indexRef); */

    /* printf("\nOptimum lambda found: %lf   with objective value  %lf    and derivative  %lf\n", lamOpti, fObj, derObj);  */

    /* printf("\n\n\nTESTING FREE SPACE UPDATE, ONE ANCHOR POINT AND XHAT ON THE BOUNDARY\n\n\n"); */

    /* int indexAccepted, x1_ind, xHat_ind; */
    /* double T0, T1, grad0[2], grad1[2], indexRef; */

    /* grad0[0] = 0.70710678; */
    /* grad0[1] = 0.70710678; */
    /* grad1[0] = 0.20485516; */
    /* grad1[1] = 0.9787923; */
    /* indexAccepted = 2; */
    /* x1_ind = 468; */
    /* xHat_ind = 1; */
    /* T0 = 1.4142135623730951; */
    /* T1 = 1.9543124768462696; */
    /* indexRef = 1.0; */

    /* printf("T1: %lf\n", T1); */
    /* info_update_initTwo(info_update, indexAccepted, x1_ind, xHat_ind, T0, grad0, T1, grad1, indexRef); */

    /* print_info_update(info_update); */

    /* anchorHatonBoundary_freeSpaceUpdate(eik_g1->triM_2D, info_update, 0); */

    /* print_info_update(info_update); */


    /* printf("\n\n\nTESTING FREE SPACE UPDATE, JUST X0 ON THE BOUNDARY\n\n\n"); */

    /* int indexAccepted, x1_ind, xHat_ind; */
    /* double T0, T1, grad0[2], grad1[2], indexRef; */

    /* grad0[0] = -0.65702763; */
    /* grad0[1] = -0.7538665; */
    /* grad1[0] = -0.95178913; */
    /* grad1[1] = 0.30675307; */
    /* indexAccepted = 1; */
    /* x1_ind = 168; */
    /* xHat_ind = 468; */
    /* T0 = 1.6301719299656714; */
    /* T1 = 0.6327956253704661; */
    /* indexRef = 1.0; */

    /* info_update_initTwo(info_update, indexAccepted, x1_ind, xHat_ind, T0, grad0, T1, grad1, indexRef); */

    /* print_info_update(info_update); */

    /* justx0Boundary_TwoPointUpdate(eik_g1->triM_2D, info_update); */

    /* print_info_update(info_update); */


    /* printf("\n\n\nTESTING FREE SPACE UPDATE, JUST XHAT ON THE BOUNDARY\n\n\n"); */

    /* int indexAccepted, x1_ind, xHat_ind; */
    /* double T0, T1, grad0[2], grad1[2], indexRef; */

    /* grad1[0] = 0; */
    /* grad1[1] = 1; */
    /* grad0[0] = -0.84990039; */
    /* grad0[1] = 0.52694338; */
    /* indexAccepted = 464; */
    /* x1_ind = 614; */
    /* xHat_ind = 18; */
    /* T1 = 4.742027478979956; */
    /* T0 = 0; */
    /* indexRef = 1.0; */

    /* info_update_initTwo(info_update, indexAccepted, x1_ind, xHat_ind, T0, grad0, T1, grad1, indexRef); */

    /* print_info_update(info_update); */

    /* justxHatBoundary_TwoPointUpdate(eik_g1->triM_2D, info_update); */

    /* print_info_update(info_update); */



    /* printf("\n\n\nTESTING FREE SPACE UPDATE, JUST XHAT ON THE BOUNDARY\n\n\n"); */

    /* int indexAccepted, x1_ind, xHat_ind; */
    /* double T0, T1, grad0[2], grad1[2], indexRef; */

    /* grad0[0] = 0; */
    /* grad0[1] = 1; */
    /* grad1[0] = -0.84990039; */
    /* grad1[1] = 0.52694338; */
    /* indexAccepted = 614; */
    /* x1_ind = 464; */
    /* xHat_ind = 18; */
    /* T0 = 4.742027478979956; */
    /* T1 = 0; */
    /* indexRef = 1.0; */

    /* info_update_initTwo(info_update, indexAccepted, x1_ind, xHat_ind, T0, grad0, T1, grad1, indexRef); */

    /* print_info_update(info_update); */

    /* justxHatBoundary_TwoPointUpdate(eik_g1->triM_2D, info_update); */

    /* print_info_update(info_update); */

    /* printf("TESTING FUNCTIONS RELATED TO A TWO STEP UPDATE"); */

    /* double gradient[2], lambda, mu, T0, grad0[2], T1, grad1[2], x0[2], x1[2], x2[2], xHat[2], B0[2], B2[2]; */
    /* double indexRef_01, indexRef_02, optimizers[2], grad_f[2]; */

    /* lambda = 0.3; */
    /* mu = 0.2; */
    /* T0 = 1.4; */
    /* grad0[0] = 0.4472136; */
    /* grad0[1] = 0.89442719; */
    /* grad1[0] = -0.89442719; */
    /* grad1[1] = 0.4472136; */
    /* T1 = 1.6; */
    /* x0[0] = 0; */
    /* x0[1] = 0; */
    /* x1[0] = -3; */
    /* x1[1] = 1; */
    /* x2[0] = 0; */
    /* x2[1] = 0.3; */
    /* xHat[0] = 0.3; */
    /* xHat[1] = 2.5; */
    /* B0[0] = 0.01; */
    /* B0[1] = 1; */
    /* B2[0] = 0.1; */
    /* B2[1] = 0.8; */
    /* indexRef_01 = 1; */
    /* indexRef_02 = 1.6; */

    /* projectedGradient_TwoStep(optimizers, 0, 1, 0, 1, T0, grad0, T1, grad1, x0, x1, x2, xHat, B0, B2, indexRef_01, indexRef_02, 0.000001, 50); */

    /* printf("Optimizers found %lf  %lf\n", optimizers[0], optimizers[1]); */

    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////

    printf("TESTING THE SHOOT + CREEPING RAY UPDATE AND FININD MUMIN\n\n");

    double muMin, x0[2], x1[2], xHat[2], xR[2], BHat[2], BR[2], T0, T1, grad0[2], grad1[2];

    muMin = 0.25;

    x0[0] = 3;
    x0[1] = -2;
    
    x1[0] = -1;
    x1[1] = -2;
    
    xR[0] = 0;
    xR[1] = -1;
    
    xHat[0] = 1;
    xHat[1] = 0;
    
    BHat[0] = 0;
    BHat[1] = 1;
    
    BR[0] = 1;
    BR[1] = 0;

    T0 = 4.123105625617661;
    T1 = 1.0;

    grad0[0] = 0.9701425;
    grad0[1] = 0.24253563;
    grad1[0] = 0.0;
    grad1[1] = 1.0;

    double fObj, mu;
    
    printf("Testing the function value at mu = 0.5\n");

    mu = 0.5;

    fObj = fobjective_shootCr(mu, x0, x1, xHat, xR, BHat, BR, T0, T1, grad0, grad1, 1.0);
    printf("With mu = 0.5: %lf:\n", fObj);
    
    muMin = find_minMu(0, x0, x1, xHat, xR, BHat, BR, 0.00001, 50);
    printf("The minimum mu found is: %lf\n", muMin);

    double fObj_min;
    fObj_min = fobjective_shootCr(muMin, x0, x1, xHat, xR, BHat, BR, T0, T1, grad0, grad1, 1.0);
    printf("With muMin: %lf:\n", fObj_min);

    double muOpt, fObj_opt;
    muOpt = projectedGradient_shootCr(muMin, muMin, 1, x0, x1, xHat, xR, BHat, BR, T0, T1, grad0, grad1, 0.00001, 50, 1.0);
    printf("\n\n\n\n\nOptimum mu found for a shoot and creep update: %lf\n", muOpt);
    fObj_opt = fobjective_shootCr(muOpt, x0, x1, xHat, xR, BHat, BR, T0, T1, grad0, grad1, 1.0);
    printf("With optimum mu: %lf\n", fObj_opt);


    eik_grid_dealloc(&eik_g1);



}
