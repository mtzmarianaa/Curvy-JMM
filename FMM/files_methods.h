/*
General methods for dealing with txt files, specially ints or doubles separated by commas
*/
#pragma once

void separateARowInt(char *line, int nElementsRow, int *row);

int numberElementsInRow(char *line, int nCharInLine);

int numLinesInFile(const char *pathFile);