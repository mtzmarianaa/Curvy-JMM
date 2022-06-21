#include "neighbors.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// int numberNeighborsFound(char *line, int nCharInLine);

// void separateARow(char *line, int nCharInLine, int *neighborsRow);

// void printThisLinesNeighbors(int *neighborsRow, int SizeRow);

int main(){

    int N;
    char const *pathNeighbors;

    pathNeighbors = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/Neigh.txt";
    N = numLinesInFile(pathNeighbors);

    neighborsRS *neighbors;

    neighborsRSalloc_n(&neighbors, N); // allocate those N lines of organized information


    printf("\n Size 1 %p", &neighbors);
    printf("\n Size 2%lu", sizeof(neighbors));
    

    neighbors_init(neighbors, pathNeighbors, N);

    printf("\n Size 1 %p", &neighbors);
    printf("\n Size 2%lu", sizeof(neighbors));
    printf("\n Size 3%lu", sizeof(*neighbors));

    printf("\n Printing all the neighbors found: \n");
    printAllNeighbors(neighbors, N);


//    FILE *fp;
//    char * line = NULL;
//    char *tempLine;
//     size_t len = 0;
//     ssize_t read;
//     int nNei;
//     int i = 0;
  
//    fp = fopen("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/Neigh.txt","r");
//    if (fp == NULL)
//         exit(EXIT_FAILURE);

//     while ((read = getline(&line, &len, fp)) != -1) {
//         int nCharInLine = (int) read -1;
//         tempLine = line;
//         printf("Retrieved line of length %d:\n", nCharInLine);
//         printf("Iteration %d \n", i);
//         printf("%s", line);
//         printf("Memory allocation of line: %p \n", &line);
//         printf("Memory allocation of templine: %p \n", &tempLine);
//         i ++;
//         printf("Temporary line: %s", tempLine);
//         nNei = numberNeighborsFound(tempLine, nCharInLine);
//         printf("Number of indices found in curreny line: %d :\n", nNei);
//         printf("\n");
//         int *neighborsRow = malloc(nNei*sizeof(int));
//         separateARow(tempLine, nNei, neighborsRow);
//         printThisLinesNeighbors(neighborsRow, nNei);
//         printf("\n\n");
//         }
    

//     fclose(fp);
//     if (line)
//         free(line);
//     exit(EXIT_SUCCESS);

}

// int numberNeighborsFound(char *line, int nCharInLine) {
//     int i, count;
//     count = 0;
//     printf("%s  \n", line);
//     for (i=0; i<nCharInLine; i++){
//         if ((char) line[i] == ','){
//             count+= 1;
//         }
//     }
//     count += 1; // plus 1 because in theory we're counting the ,
//     return count;
// }

// void separateARow(char *line, int nNei, int *neighborsRow) {
//     // In this method given a line of text read from a file with nNeighRow 
//     // indices of neighbors we separate them and put it in neighborsRow.
//     char *char_part;
//     int int_part;
//     int_part = strtol(line, &char_part, 10);
//     neighborsRow[0] = int_part;
//     char *end = char_part;
//     char_part = char_part + 1;
//     for (int i = 1; i<nNei; i++){
//         int_part = strtol(char_part++, &end , 10);
//         neighborsRow[i] = int_part;
//         char_part = end + 1;
//     }
// }



// void printThisLinesNeighbors(int *neighborsRow, int SizeRow) {
//     for (int i = 0; i<SizeRow; i++){
//         printf(" %d ", neighborsRow[i]);
//     }
// }