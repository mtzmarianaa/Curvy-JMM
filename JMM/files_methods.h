#pragma once

/*
General methods for dealing with txt files, specially ints or doubles separated by commas
*/

#include <stddef.h>


void separateARowSize_t(char *line, int nElementsRow, size_t *row);

void separateARowInt(char *line, int nElementsRow, int *row);

void separateARowDb(char *line, int nElementsRow, double *row);

int numberElementsInRow(char *line, int nCharInLine);

int numLinesInFile(const char *pathFile);

void readIntColumn(const char *pathFile, int *column);

void readDbColumn(const char *pathFile, double *column);

void saveTimes(double *times, const char *pathFile);

void read_n2File_double(double *flattenMatrix, char const *pathFile);

void read_n2File_size_t(size_t *flattenMatrix, char const *pathFile);

void read_n2File_int(int *flattenMatrix, char const *pathFile);

void read_n3File_int(int *flattenMatrix, char const *pathFile);

void read_n3File_size_t(size_t *flattenMatrix, char const *pathFile);
