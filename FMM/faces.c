/*
In order to plot a triangle mesh in Python using triplot we need the triangles
A triangle is defined by 3 points, i.e. by 3 indeces that represent those points
*/

#include "faces.h"
#include "files_methods.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


void faces_alloc(facesS **faces) {
  *faces = malloc(sizeof(facesS));
  assert(*faces != NULL);
}

void faces_dealloc(facesS **faces) {
    free(faces);
    *faces = NULL;
}

void faces_init(facesS *faces, int (*points)[3], int nFaces) {
    faces->points = points;
    faces->nFaces = nFaces;
}

void print_faces(facesS *faces) {
    printf("Number of faces: %d \n", faces->nFaces);
    for(int i = 0; i< faces->nFaces; i ++) {
        printf("Face %d number has:    ", i);
        printf("%d  | %d  | %d  \n", faces->points[i][0], faces->points[i][1], faces->points[i][2]);
    }
}

void faces_initFromFile(facesS *faces, char const *pathFaces) {
  // This method is going to store what is found un the pathFaces text file 
    FILE *fp; // the file were going to read
    char * line = NULL; // each line of the file
    int i = 0;
    size_t len = 0; // length of each line of the file
    ssize_t read; // reading each line in the file
    int *row; // temporary row of the file we're reading
    row = malloc(3*sizeof(int));
    row[0] = 0;
    row[1] = 0;
    row[2] = 0;
    int nFaces;

    nFaces = numLinesInFile(pathFaces); // we have the number of faces in the file
    faces->nFaces = nFaces;
    faces->points = malloc(nFaces*3*sizeof(int));

    printf("\ntrying to open to file\n");
    // Check if the file exists under that path
    fp = fopen(pathFaces, "r");
    // if the file doesnt exist or for some reason it can't be opened:
    if (fp == NULL) {
        printf("No such file");
        exit(EXIT_FAILURE);
    }
    else{
        printf("\nFile successfully found\n");
    }

    // scan such file

        while ((read = getline(&line, &len, fp)) != -1) {
        printf("Iteration %d \n", i);
        printf("%s", line);
        printf("\n");
        separateARowInt(line, 3, row);
        printf("First element found: %d\n", row[0]);
        printf("Second element found: %d\n", row[1]);
        printf("Third element found: %d\n", row[2]);
        faces->points[i][0] = row[0];
        faces->points[i][1] = row[1];
        faces->points[i][2] = row[2];
        printf("\n\n");
        i ++;
    }
    fclose(fp);
    if (line)
        free(line);
}
