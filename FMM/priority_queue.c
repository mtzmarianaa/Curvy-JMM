/* PRIORITY QUEUE

Here is the implementation of the prioriry queue for the 2D FMM. Explanaition

   - eik_queue: Eikonal values considered, binary tree in an array, this is the thing to heapify
   - index_queue : index (grid coordinates) if the Eikonal values considered in the queue
   - size : current non empty entries of the binary tree

*/

#include "priority_queue.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

struct Priority_queue {
  double *queue_vals; // if we need more we'll add more
  int *queue_index; // same here
  int size; // current occupied size occupied
  int maxSize; // max sized of queue_vals and queue_index currently allowed
};

void priority_queue_alloc(p_queue **Priority_queue ) {
  *Priority_queue = malloc(sizeof(p_queue));
  assert(*Priority_queue != NULL);
}

void priority_queue_dealloc(p_queue **Priority_queue ) {
  free(*Priority_queue);
  *Priority_queue = NULL;
}

void priority_queue_init( p_queue *Priority_queue  ) 
{
  Priority_queue->maxSize = 2;
  Priority_queue->queue_vals = malloc( Priority_queue->maxSize*sizeof(double)  );
  Priority_queue->queue_index = malloc( Priority_queue->maxSize*sizeof(int) );
  Priority_queue->size = 0;
  assert( Priority_queue != NULL  ); // the queue should not be null if initialized
}


void grow_queue( p_queue *Priority_queue ) 
{
  Priority_queue->maxSize *= 2;
  Priority_queue->queue_vals = realloc( Priority_queue->queue_vals, Priority_queue->maxSize*sizeof(double)  );
  Priority_queue->queue_index = realloc( Priority_queue->queue_index, Priority_queue->maxSize*sizeof(int) );
}


void swap_double(double *a, double *b)
{
  double temp = *b;
  *b = *a;
  *a = temp;
}

void swap_int(int *a, int *b)
{
  int temp = *b;
  *b = *a;
  *a = temp;
}

void heapify(p_queue *Priority_queue, int i)
{
  if (Priority_queue->size == 1)
  {
    // printf("Single element in the heap, heapified done");
  }
  else
  {
    int smallest = i;
    int l = 2 * i + 1; // left child
    int r = 2 * i + 2; // right child
    // find if the children are smaller or not, if they are we need to swap them
    // BUT THE COMPARISON IS DONE IN THE eik_queue
    if (l < Priority_queue->size && Priority_queue->queue_vals[l] < Priority_queue->queue_vals[smallest])
      smallest = l;
    if (r < Priority_queue->size && Priority_queue->queue_vals[r] < Priority_queue->queue_vals[smallest])
      smallest = r;
    if (smallest != i)
    {
      swap_double(&Priority_queue->queue_vals[i], &Priority_queue->queue_vals[smallest]); // swap in the eik_queue
      swap_int(&Priority_queue->queue_index[i], &Priority_queue->queue_index[smallest]); // swap in the index_queue
      heapify(Priority_queue, smallest); // recurrencia
    }
  }
}

void insert(p_queue *Priority_queue, double newNum, int newIndex)
{
  int i;
  if (Priority_queue->size == 0) // First time we insert an element in the tree
  {
    Priority_queue->size = 1;
    Priority_queue->queue_vals[0] = newNum;
    Priority_queue->queue_index[0] = newIndex;
  }
  else // if there were elements in the tree before inserting
  {
    Priority_queue->size += 1;
    if ( Priority_queue->size >= Priority_queue->maxSize  ) {
      grow_queue( Priority_queue );
    }
    Priority_queue->queue_vals[Priority_queue->size -1 ] = newNum;
    Priority_queue->queue_index[Priority_queue->size -1 ] = newIndex;
    for (i = Priority_queue->size / 2 ; i >= 0; i--)
    {
      heapify(Priority_queue, i);
    }
  }
}



void insert_end(p_queue *Priority_queue, double newNum, int newIndex)
{
    Priority_queue->size += 1;
    Priority_queue->queue_vals[Priority_queue->size -1 ] = newNum;
    Priority_queue->queue_index[Priority_queue->size -1] = newIndex;
}



void delete_findValue(p_queue *Priority_queue, double num)
{
  int i;
  for (i = 0; i < Priority_queue->size; i++)
  {
    if (num == Priority_queue->queue_vals[i] ) // find the value
      break;
  }

  swap_double(&Priority_queue->queue_vals[i], &Priority_queue->queue_vals[Priority_queue->size - 1]); 
  swap_int(&Priority_queue->queue_index[i], &Priority_queue->queue_index[Priority_queue->size - 1]); 

  Priority_queue->size -= 1; // make it "forget", after size-1 everything is ignored
  for (i = Priority_queue->size / 2 ; i >= 0; i--)
  {
    heapify(Priority_queue, i);
  }
}

void delete_findIndex(p_queue *Priority_queue, int ind)
{
  int i;
  for (i = 0; i < Priority_queue->size; i++)
  {
    if (ind == Priority_queue->queue_index[i])
      break;
  }

  swap_double(&Priority_queue->queue_vals[i], &Priority_queue->queue_vals[Priority_queue->size - 1]); 
  swap_int(&Priority_queue->queue_index[i], &Priority_queue->queue_index[Priority_queue->size - 1]); 

  Priority_queue->size -= 1; // make it "forget", after size-1 everything is ignored
  for (i = Priority_queue->size / 2 ; i >= 0; i--)
  {
    heapify(Priority_queue, i);
  }
}

int indexRoot(p_queue *Priority_queue){
  return Priority_queue->queue_index[0];
}

double valueRoot(p_queue *Priority_queue){
  return Priority_queue->queue_vals[0];
}

void deleteRoot(p_queue *Priority_queue)
{
  int i;
  // we dont need to look for the index, we know its on the 0th position
  swap_double(&Priority_queue->queue_vals[0], &Priority_queue->queue_vals[Priority_queue->size - 1]); 
  swap_int(&Priority_queue->queue_index[0], &Priority_queue->queue_index[Priority_queue->size - 1]); 
  Priority_queue->size -= 1; 
  for (i = Priority_queue->size / 2; i >= 0; i--)
  {
    heapify(Priority_queue, i);
  }
}


void printeik_queue(p_queue *Priority_queue)
{
  int i;
  printf("Eikonal values:");
  for (i = 0; i < Priority_queue->size; ++i)
    printf("%lf ", Priority_queue->queue_vals[i]);
  printf("\n");
  printf("Indices:");
  for (int i = 0; i < Priority_queue->size; ++i)
    printf("%d ", Priority_queue->queue_index[i]);
  printf("\n");
  printf("Current size: %d", Priority_queue->size);
}

void update(p_queue *Priority_queue, double new_valConsidered, int index)
{
    int i, j;
    printf("Before updating the queue looks like this: \n");
    printeik_queue(Priority_queue);
    printf("\n");
    // First find the current value associated with index
    for (i = 0; i < Priority_queue->size; i++)
    {
        if (index == Priority_queue->queue_index[i])
        break;
    }
    // Then, if the new value considered is smaller than the current value
    if ( Priority_queue->queue_vals[i] > new_valConsidered  )
    {
        printf("New value updated in %d to the queue: %fl\n", index, new_valConsidered);
        Priority_queue->queue_vals[i] = new_valConsidered;
        for (j = Priority_queue->size / 2 ; j >= 0; j--)
        {
            heapify(Priority_queue, j);
        }
        printf("The queue so far looks like this: \n");
        printeik_queue(Priority_queue);
    }
    else{
      printf("\nWe were trying to update at %d the value of %fl but the current value in the queue is better.\n", index, new_valConsidered);
    }
}



 double get_valueAtIndex(p_queue *Priority_queue, int index)
{
    int i;
    for (i = 0; i < Priority_queue->size; i ++)
    {
        if (index == Priority_queue->queue_index[i])
        break;
    }
    return Priority_queue->queue_vals[i];
}

int getSize(p_queue *Priority_queue)
{
  return Priority_queue->size;
}

int getIndicesInQueue(p_queue *Priority_queue)
{
  return *Priority_queue->queue_index;
}