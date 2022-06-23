#pragma once

typedef struct {
  int *from;
  int *to;
  int nFacets;
} facetsS;

void facets_alloc(facetsS **facets );

void facets_dealloc(facetsS **facets);

void facets_init(facetsS *facets, int *from, int *to, int nFacets);

void print_facets(facetsS *facets);

void facets_initFromFile(facetsS *facets, char const *pathFacets);