#pragma once


typedef struct{
    int len;
    double *individual_path_x;
    double *individual_path_y;
    int maxSize;
} pathS;

void path_alloc_n(pathS **path, int nPoints);

void path_dealloc(pathS **path);

void grow_This_path( pathS *path, int indexHat );

void path_init(pathS *path, int nPoints);

void insertOneToPath(pathS *path, int indexHat, int indexPath);

void insertAfterAcceptedToPath(pathS *path, int indexHat, int lastInPath);

void printPaths(pathS *path, int nPoints);