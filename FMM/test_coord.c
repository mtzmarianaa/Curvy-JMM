#include "coord.h"

#include <stdio.h>
#include <stdlib.h>

int main(){

    // NAIVE TESTING
    coordS *coord1;
    coord_alloc(&coord1);
    double xn[] = {1.0,3.0,5.0,7.0,9.0};
    double yn[] = {2.0,4.0,6.0,8.0,10.0};
    double *x = xn;
    double *y = yn;

    coord_init(coord1, x, y, sizeof(xn)/sizeof(double));

    print_coord(coord1);

    printf("Size of xn: %lu", sizeof(xn)/sizeof(double));

    printf("\nNumber of points: %d", coord1->nPoints);


    // TESTING FROM FILE
    coordS *coord2;
    coord_alloc(&coord2);
    const char *pathPoints;
    pathPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/MeshPoints.txt";
    coord_initFromFile(coord2, pathPoints);
    printf("\n\n\n  PRINTING THE COORDINATES AS FOUND IN FILE \n\n");
    print_coord(coord2);


}