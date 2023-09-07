/* MEMORY CACHE FOR EACH NODE

   For each node we want to save the information of each
   of its children's update (index0, index1, index2,
   indexChild, firstTriangle).
   This is a hash table, we sort the information of the
   children according to the children's indices. 

 */

#include "memCache.h"
#include <math.h>


void fanCache_alloc(fanCacheS **fanCache) {
  *fanCache = malloc(sizeof(fanCacheS));
  assert(*fanCache != NULL);
}

void fanCache_dealloc(fanCacheS **fanCache) {
  free(*fanCache);
  *fanCache = NULL;
}

void fanCache_init(fanCacheS *fanCache, size_t index0, size_t index1,
		   size_t index2, size_t indexChild, size_t firstTriangle) {
  fanCache->index0 = index0;
  fanCache->index1 = index1;
  fanCache->index2 = index2;
  fanCache->indexChild = indexChild;
  fanCache->firstTriangle = firstTriangle;
}


void memCache_alloc(memCacheS **memCache) {
  *memCache = malloc(sizeof(memCacheS));
  assert(*memCache != NULL);
}

void memCache_dealloc(memCacheS **memCache) {
  free(*memCache);
  *memCache = NULL;
}

void memCache_init(memCacheS *memCache) {
  memCache->maxSize = 2;
  memCache->size = 0;
  memCache->indexChildren = malloc( memCache->maxSize*sizeof(size_t));
  memCache->fanChildren = malloc( memCache->maxSize*sizeof(fanCacheS));
}

void grow_memCache(memCacheS *memCache) {
  memCache->maxSize *= 2;
  memCache->indexChildren = realloc( memCache->indexChildren, memCache->maxSize*sizeof(size_t));
  memCache->fanChildren = realloc( memCache->fanChildren, memCache->maxSize*sizeof(fanCacheS));
}

void swap_size_t(size_t *a, size_t *b) {
  size_t temp = *b;
  *b = *a;
  *a = temp;
}

void swap_fanCache(fanCacheS *a, fanCacheS *b) {
  fanCacheS temp = *b;
  *b = *a;
  *a = temp;
}


void heapify_memCache(memCacheS *memCache, int i) {
  int smallest, l, r;
  if (memCache->size > 1)
  {
    smallest = i;
    l = 2 * i + 1; // left child
    r = 2 * i + 2; // right child
    // find if the children are smaller or not, if they are we need to swap them
    if (l < memCache->size && memCache->indexChildren[l] < memCache->indexChildren[smallest])
      smallest = l;
    if (r < memCache->size && memCache->indexChildren[r] < memCache->indexChildren[smallest])
      smallest = r;
    if (smallest != i)
    {
      swap_size_t(&memCache->indexChildren[i], &memCache->indexChildren[smallest]);
      swap_fanCache(&memCache->fanChildren[i], &memCache->fanChildren[smallest]); 
      heapify_memCache(memCache, smallest); // recurrencia
    }
  }
}

void insert_memCache(memCacheS *memCache, size_t index0, size_t index1, size_t index2,
		     size_t indexChild, size_t firstTriangle) {
  int i, l;
  l = memCache->size;
  if(memCache->size == 0){
    // first time we insert something
    memCache->size = 1;
    memCache->indexChildren[l] = indexChild;
    memCache->fanChildren[l].index0 = index0;
    memCache->fanChildren[l].index1 = index1;
    memCache->fanChildren[l].index2 = index2;
    memCache->fanChildren[l].indexChild = indexChild;
    memCache->fanChildren[l].firstTriangle = firstTriangle;
  }
  else // meaning that there were elements in the tree before
  {
    memCache->size += 1;
    if( memCache->size >= memCache->maxSize){
      // got bigger, we need to reallocate
      grow_memCache(memCache);
    }
    // insert information at the end, then heapify
    memCache->indexChildren[l] = indexChild;
    memCache->fanChildren[l].index0 = index0;
    memCache->fanChildren[l].index1 = index1;
    memCache->fanChildren[l].index2 = index2;
    memCache->fanChildren[l].indexChild = indexChild;
    memCache->fanChildren[l].firstTriangle = firstTriangle;
    for( i = memCache->size/2; i>=0; i--){
      heapify_memCache(memCache, i);
    }
  }
}


void delete_child(memCacheS *memCache, size_t indexChild) {
  // deletes the information of a child (maybe its not longer its child)
  // MAKE SURE THIS CHILD INDEX IS ACTUALLY HERE
  int i, sc;
  sc = memCache->size - 1;
  for( i = 0; i<memCache->size; i++){
    if(indexChild == memCache->indexChildren[i]) // we've found the place of the child
      break;
  }
  // swap this with the last element (we are going to "forget" this last element)
  swap_size_t(&memCache->indexChildren[i], &memCache->indexChildren[sc]); 
  swap_fanCache(&memCache->fanChildren[i], &memCache->fanChildren[sc]);
  memCache->size -= 1;
  for( i = memCache->size/2; i >= 0; i--){
    heapify_memCache(memCache, i);
  }
}


bool isChild(memCacheS *memCache, size_t indexKid) {
  // true or false, if this node is in this memCache (i.e. if this
  // node is the child of this point
  int i;
  bool res = false;
  for( i = 0; i<memCache->size; i++){
    if(memCache->indexChildren[i] == indexKid){
      return true;
    }
  }
  return res;
}




