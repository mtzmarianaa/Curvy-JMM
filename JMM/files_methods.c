#include "files_methods.h"

#include "feature_test.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

void separateARowSize_t(char *line, int nElementsRow, size_t *row) {
    // In this method given a line of text read from a file
    // indices of neighbors we separate them and put it in row.
    char *char_part;
    int int_part;
    int_part = strtol(line, &char_part, 10);
    row[0] = (size_t)int_part;
    char *end = char_part;
    char_part = char_part + 1;
    for (int i = 1; i<nElementsRow; i++){
        int_part = strtol(char_part++, &end , 10);
        row[i] = (size_t)int_part;
        char_part = end + 1;
    }
}

void separateARowInt(char *line, int nElementsRow, int *row) {
    // In this method given a line of text read from a file
    // indices of neighbors we separate them and put it in row.
    char *char_part;
    int int_part;
    int_part = strtol(line, &char_part, 10);
    row[0] = int_part;
    char *end = char_part;
    char_part = char_part + 1;
    for (int i = 1; i<nElementsRow; i++){
        int_part = strtol(char_part++, &end , 10);
        row[i] = int_part;
        char_part = end + 1;
    }
}

void separateARowDb(char *line, int nElementsRow, double *row) {
    // In this method given a line of text read from a file
    // indices of neighbors we separate them and put it in row.
    const char *sep = ", ";
    char *token;
    int i;
    token = strtok(line, sep); //first double found in the line
    //printf("First token %s\n", token);
    row[0] = atof(token);
    //printf("First token saved %fl\n", row[0]);
    for (i = 1; i<nElementsRow; i++){
        token = strtok(NULL, sep);
        //printf("Token %d: %s\n",i, token);
        row[i] = atof(token);
        //printf("Token saved %d: %fl\n",i, row[i]);
    }
}

int numberElementsInRow(char *line, int nCharInLine) {
    int i, count;
    count = 0;
    //printf("%s  \n", line);
    for (i=0; i<nCharInLine; i++){
        if ((char) line[i] == ','){
            count+= 1;
        }
    }
    count += 1; // plus 1 because in theory we're counting the ,
    return count;
}

int numLinesInFile(const char *pathFile){
    int nLines = 0;
    FILE *fp; // the file were going to read
    char c;

    fp = fopen(pathFile, "r");
    // if the file doesnt exist or for some reason it can't be opened:
    if (fp == NULL) {
        printf("No such file");
        exit(EXIT_FAILURE);
    }

    for (c = getc(fp); c != EOF; c = getc(fp)){
        if (c == '\n'){
            nLines ++;
        }
    }
    fclose(fp);
    return nLines;
}

void readIntColumn(const char *pathFile, int *column){
    // this method is going to store what is found on the pathFile text file to a column (should only be a column stored in that file)
    FILE *fp;
    char * line = NULL;
    int i = 0;
    size_t len = 0;
    ssize_t read;
    int row[1];
    row[0] = 0;
    printf("\ntrying to open to file\n");
    fp = fopen(pathFile, "r");
    // if the file doesnt exists or it cant be opened
    if(fp == NULL) {
        printf("No such file for integer column");
        exit(EXIT_FAILURE);
    }
    else{
        printf("\nFile successfully found \n");
    }
    // scan the file
    fp = fopen(pathFile, "r");
    while( (read = getline(&line, &len, fp)) != -1 ) {
        separateARowInt(line, 1, row);
        column[i] = row[0];
        i++;
    }
    fclose(fp);
    if(line){
        free(line);
    }
}


void readDbColumn(const char *pathFile, double *column) {
    // this method is going to store what is found on the pathFile text file to a column (should only be a column stored in that file)
    FILE *fp;
    char * line = NULL;
    int i = 0;
    size_t len = 0;
    ssize_t read;
    double row[1];
    row[0] = 0;
    printf("\ntrying to open to file\n");
    fp = fopen(pathFile, "r");
    // if the file doesnt exists or it cant be opened
    if(fp == NULL) {
        printf("No such file for double column");
        exit(EXIT_FAILURE);
    }
    else{
        printf("\nFile successfully found \n");
    }
    // scan the file
    fp = fopen(pathFile, "r");
    while( (read = getline(&line, &len, fp)) != -1 ) {
        separateARowDb(line, 1, row);
        column[i] = row[0];
        i++;
    }
    fclose(fp);
    if(line){
        free(line);
    }
}


void saveTimes(double *times, const char *pathFile){
  FILE *fp;
  fp = fopen(pathFile, "wb");
  fwrite(times, sizeof(double), 1, fp);
  fclose(fp);
}


void read_n2File_double(double *flattenMatrix, char const *pathFile) {
  // this method reads an nx2 double matrix stored in a text file
    FILE *fp;
    char *line = NULL;
    int i = 0;
    size_t len = 0; // length of each line of the file
    ssize_t read; // reading each line in the file
    double row[2]; // temporary row
    row[0] = 0;
    row[1] = 0;

    printf("\ntrying to open to file\n");
    // Check if the file exists under that path
    fp = fopen(pathFile, "r");
    // if the file doesnt exist or for some reason it can't be opened:
    if (fp == NULL) {
        printf("No such file, try again");
        exit(EXIT_FAILURE);
    }
    else{
        printf("\nFile successfully found\n");
    }

    // scan such file
    fp = fopen(pathFile, "r");
    while ((read = getline(&line, &len, fp)) != -1) {
      separateARowDb(line, 2, row);
      flattenMatrix[i] = row[0];
      flattenMatrix[i + 1] = row[1];
      i = i + 2;
    }
    fclose(fp);
    if(line){
        free(line);
    }
}

void read_n2File_size_t(size_t *flattenMatrix, char const *pathFile) {
  // this method reads an nx2 int matrix stored in a text file
    FILE *fp;
    char *line = NULL;
    int i = 0;
    size_t len = 0; // length of each line of the file
    ssize_t read; // reading each line in the file
    size_t row[2]; // temporary row
    row[0] = 0;
    row[1] = 0;

    printf("\ntrying to open to file\n");
    // Check if the file exists under that path
    fp = fopen(pathFile, "r");
    // if the file doesnt exist or for some reason it can't be opened:
    if (fp == NULL) {
        printf("No such file, sorry");
        exit(EXIT_FAILURE);
    }
    else{
        printf("\nFile successfully found\n");
    }

    // scan such file
    fp = fopen(pathFile, "r");
    while ((read = getline(&line, &len, fp)) != -1) {
      separateARowSize_t(line, 2, row);
      flattenMatrix[i] = row[0];
      flattenMatrix[i+1] = row[1];
      i = i + 2;
    }
    fclose(fp);
    if(line){
        free(line);
    }
}

void read_n2File_int(int *flattenMatrix, char const *pathFile) {
  // this method reads an nx2 int matrix stored in a text file
    FILE *fp;
    char *line = NULL;
    int i = 0;
    size_t len = 0; // length of each line of the file
    ssize_t read; // reading each line in the file
    int row[2]; // temporary row
    row[0] = 0;
    row[1] = 0;

    printf("\ntrying to open to file\n");
    // Check if the file exists under that path
    fp = fopen(pathFile, "r");
    // if the file doesnt exist or for some reason it can't be opened:
    if (fp == NULL) {
        printf("No such file, sorry");
        exit(EXIT_FAILURE);
    }
    else{
        printf("\nFile successfully found\n");
    }

    // scan such file
    fp = fopen(pathFile, "r");
    while ((read = getline(&line, &len, fp)) != -1) {
      separateARowInt(line, 2, row);
      flattenMatrix[i] = row[0];
      flattenMatrix[i+1] = row[1];
      i = i + 2;
    }
    fclose(fp);
    if(line){
        free(line);
    }
}



void read_n3File_int(int *flattenMatrix, char const *pathFile) {
  // this method reads an nx3 int matrix stored in a text file
    FILE *fp;
    char *line = NULL;
    int i = 0;
    size_t len = 0; // length of each line of the file
    ssize_t read; // reading each line in the file
    int row[3]; // temporary row
    row[0] = 0;
    row[1] = 0;
    row[2] = 0;

    printf("\ntrying to open to file\n");
    // Check if the file exists under that path
    fp = fopen(pathFile, "r");
    // if the file doesnt exist or for some reason it can't be opened:
    if (fp == NULL) {
        printf("No such file, sorry");
        exit(EXIT_FAILURE);
    }
    else{
        printf("\nFile successfully found\n");
    }

    // scan such file
    fp = fopen(pathFile, "r");
    while ((read = getline(&line, &len, fp)) != -1) {
      separateARowInt(line, 3, row);
      flattenMatrix[i] = row[0];
      flattenMatrix[i+1] = row[1];
      flattenMatrix[i+2] = row[2];
      i = i + 3;
    }
    fclose(fp);
    if(line){
        free(line);
    }
}



void read_n3File_size_t(size_t *flattenMatrix, char const *pathFile) {
  // this method reads an nx3 int matrix stored in a text file
    FILE *fp;
    char *line = NULL;
    int i = 0;
    size_t len = 0; // length of each line of the file
    ssize_t read; // reading each line in the file
    size_t row[3]; // temporary row
    row[0] = 0;
    row[1] = 0;
    row[2] = 0;

    printf("\ntrying to open to file\n");
    // Check if the file exists under that path
    fp = fopen(pathFile, "r");
    // if the file doesnt exist or for some reason it can't be opened:
    if (fp == NULL) {
        printf("No such file, sorry");
        exit(EXIT_FAILURE);
    }
    else{
        printf("\nFile successfully found\n");
    }

    // scan such file
    fp = fopen(pathFile, "r");
    while ((read = getline(&line, &len, fp)) != -1) {
      separateARowSize_t(line, 3, row);
      flattenMatrix[i] = row[0];
      flattenMatrix[i+1] = row[1];
      flattenMatrix[i+2] = row[2];
      i = i + 3;
    }
    fclose(fp);
    if(line){
        free(line);
    }
}
