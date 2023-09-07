#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>

typedef struct fanCache {
  // memory cache for each node,
  // it saves the information about the child (xHat)
  // it updated, all the information necessary
  // to initiate a triangle fan
  size_t index0;
  size_t index1;
  size_t index2;
  size_t indexChild;
  size_t firstTriangle;
} fanCacheS;


typedef struct memCache{
  // the hash table with information of all of a point's children
  size_t *indexChildren; // indices of the children - binary tree
  fanCacheS *fanChildren; // information for the fan update for its children
  int size; // current occupied size
  int maxSize; // max size allowed currently
} memCacheS;

void fanCache_alloc(fanCacheS **fanCache);

void fanCache_dealloc(fanCacheS **fanCache);

void fanCache_init(fanCacheS *fanCache, size_t index0, size_t index1,
		   size_t index2, size_t indexChild, size_t firstTriangle);

void memCache_alloc(memCacheS **memCache);

void memCache_dealloc(memCacheS **memCache);

void memCache_init(memCacheS *memCache);

void grow_memCache(memCacheS *memCache);

void swap_size_t(size_t *a, size_t *b);

void swap_fanCache(fanCacheS *a, fanCacheS *b);

void heapify_memCache(memCacheS *memCache, int i);

void insert_memCache(memCacheS *memCache, size_t index0, size_t index1, size_t index2,
		     size_t indexChild, size_t firstTriangle);

void delete_child(memCacheS *memCache, size_t indexChild);

bool isChild(memCacheS *memCache, size_t indexKid);


