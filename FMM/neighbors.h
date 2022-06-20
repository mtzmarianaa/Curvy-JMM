#pragma once

typedef struct neighborsR neighborsRS;

int numberNeighborsFound(char *line, int nCharInLine);

void separateARow(char *line, int nCharInLine, int *neighborsRow);

void printThisLinesNeighbors(int *neighborsRow, int SizeRow);