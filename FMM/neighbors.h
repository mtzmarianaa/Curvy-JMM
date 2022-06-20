#pragma once

typedef struct neighborsR neighborsRS;

void neigborsRSalloc(neighborsRS **neighbors);

void neighborsRSdealloc(neighborsRS **neighbors);

void neighbors_init(neighborsRS *neighbors, char const *pathNeighbors);

int numberNeighborsFound(char *line, int nCharInLine);

void separateARow(char *line, int nNei, int *neighborsRow);

void printThisLinesNeighbors(int *neighborsRow, int SizeRow);

void printAllNeighbors(neighborsRS *neighbors);