#include <stdio.h>
#include <stdlib.h>

int numberNeighborsFound(char *line);

int main(){
   FILE *fp;
   char * line = NULL;
    size_t len = 0;
    ssize_t read;
  
   fp = fopen("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/Neigh.txt","r");
   if (fp == NULL)
        exit(EXIT_FAILURE);

    while ((read = getline(&line, &len, fp)) != -1) {
        printf("Retrieved line of length %zu:\n", read-1);
        printf("%s", line);
        printf("\n");
        char *ptr;
        int n = strtol(line, &ptr, 10);
        printf("Integer part found: %d \n", n);
        printf("String part found: %s \n", ptr);
        }
    

    fclose(fp);
    if (line)
        free(line);
    exit(EXIT_SUCCESS);

}

int numberNeighborsFound(char *line) {
    int i, count;
    count = 0;
    for (i=0, count=0; line[i]; i++){
        count += (line[i] == ',');
    }
    count += 1;
    return count;
}