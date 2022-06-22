#pragma once

#include "files_methods.h"

typedef struct{
    int len;
    int* neis_i;
} neighborsRS;

void neighborsRSalloc(neighborsRS **neighbors);

void neighborsRSalloc_n(neighborsRS **neighbors, int N);

void neighborsRSdealloc(neighborsRS **neighbors);

void neighbors_init(neighborsRS *neighbors, char const *pathNeighbors, int N);

void printThisLinesNeighbors(int *neighborsRow, int SizeRow);

void printAllNeighbors(neighborsRS *neighbors, int N);

