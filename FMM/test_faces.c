#include "faces.h"

#include <stdio.h>
#include <stdlib.h>

int main(){

    // NAIVE TESTING
    facesS *faces1;
    faces_alloc(&faces1);
    int (*points)[3];
    points = malloc(3*3*sizeof(int));
    points[0][0] = 1;
    points[0][1] = 2;
    points[0][2] = 3;

    printf("%d \n", points[0][0]);

    points[1][0] = 2;
    points[1][1] = 4;
    points[1][2] = 6;

    printf("%d\n", points[1][1]);

    points[2][0] = 3;
    points[2][1] = 6;
    points[2][2] = 9;

    printf("%d\n\n", points[2][2]);

    printf("Trying to initialize faces1 \n");
    
    faces_init(faces1, points, 3);

    print_faces(faces1);

    printf("\n\n\n\n");

    // TESTING FROM FILE
    facesS *faces2;
    faces_alloc(&faces2);
    const char *pathFaces;
    pathFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/Faces.txt";
    faces_initFromFile(faces2, pathFaces);
    printf("\n\n\n  PRINTING THE FACES AS FOUND IN FILE \n\n");
    print_faces(faces2);



}