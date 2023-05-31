
#include "files_methods.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(){

    const char *pathFile;
    int *column;

    pathFile = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/FacesLabel.txt";
    column = malloc(60*sizeof(int));

    readIntColumn(pathFile, column);

    for(int i = 0; i<60; i++){
        printf("Item in column %d \n", column[i]);
    }

}