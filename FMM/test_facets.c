#include "facets.h"

#include <stdio.h>
#include <stdlib.h>

int main(){

    // NAIVE TESTING
    facetsS *facets1;
    facets_alloc(&facets1);
    int from_n[] = {1,3,5,7,9};
    int to_n[] = {2,4,6,8,10};
    int *from = from_n;
    int *to = to_n;

    print_facets(facets1);


    // TESTING FROM FILE
    facetsS *facets2;
    facets_alloc(&facets2);
    const char *pathFacets;
    pathFacets = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/Facets.txt";
    facets_initFromFile(facets2, pathFacets);
    printf("\n\n\n  PRINTING THE FACETS AS FOUND IN FILE \n\n");
    print_facets(facets2);



}