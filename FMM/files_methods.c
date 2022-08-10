#include "files_methods.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

void separateARowDb(char *line, int nElementsRow, double row[nElementsRow]) {
    // In this method given a line of text read from a file 
    // indices of neighbors we separate them and put it in row.
    const char sep[2] = ", ";
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
        printf("No such file");
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

void saveTimes(double times[7], const char *pathFile){
  FILE *fp;
  fp = fopen(pathFile, "wb");
  fwrite(times, sizeof(double), 7, fp);
  fclose(fp);
}

