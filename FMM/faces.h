#pragma once

typedef struct {
  int (*points)[3]; // as far as I understand, this is a pointer to an array of 3 integers
  int nFaces;
} facesS;

void faces_alloc(facesS **faces);

void faces_dealloc(facesS **faces);

void faces_init(facesS *faces, int (*points)[3], int nFaces);

void print_faces(facesS *faces);

void faces_initFromFile(facesS *faces, char const *pathFaces);