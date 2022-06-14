#include<stdio.h>
#include<string.h>
void main()
{
    int *numbers;
    int i = 0;
    FILE *file;
    if (file = fopen("MeshInfo/Faces.txt", "r"))
    {
        while (fscanf(file, "%d", &numbers[i]) != EOF)
        {
            i++;
        }
        fclose(file);
        numbers[i] = '\0';
        for (i = 0; numbers[i] != '\0'; i++)
        printf("%d\n", numbers[i]);
    }

  return 0;
}