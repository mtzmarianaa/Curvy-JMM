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

struct Priority_queue 
{
  double *queue_vals; // if we need more we'll add more
  int *queue_index; // same here
  int size; // current occupied size occupied
  int maxSize; // max sized of queue_vals and queue_index currently allowed
} ;

void priority_queue_init( p_queue *p_queueImp  ) 
{
  p_queueImp->maxSize = 2;
  p_queueImp->queue_vals = malloc( p_queueImp->maxSize*sizeof(double)  );
  p_queueImp->queue_index = malloc( p_queueImp->maxSize*sizeof(int) );
  p_queueImp->size = 0;
  assert( p_queueImp != NULL  ); // the queue should not be null if initialized
}

void priority_queue_deinit( p_queue *p_queueImp ) 
{
  free( p_queueImp->queue_vals );
  free(p_queueImp->queue_index  );
  p_queueImp->queue_vals = NULL;
  p_queueImp->queue_index = NULL;
  p_queueImp->size = 0;
  p_queueImp->maxSize = 2;
}


void grow_queue( p_queue *p_queueImp ) 
{
  p_queueImp->maxSize *= 2;
  p_queueImp->queue_vals = realloc( p_queueImp->queue_vals, p_queueImp->maxSize*sizeof(double)  );
  p_queueImp->queue_index = realloc( p_queueImp->queue_index, p_queueImp->maxSize*sizeof(int) );
}


static void swap_double(double *a, double *b)
{
  double temp = *b;
  *b = *a;
  *a = temp;
}

static void swap_int(int *a, int *b)
{
  int temp = *b;
  *b = *a;
  *a = temp;
}

static void heapify(p_queue *p_queueImp, int i)
{
  if (p_queueImp->size == 1)
  {
    printf("Single element in the heap, heapified done");
  }
  else
  {
    int smallest = i;
    int l = 2 * i + 1; // left child
    int r = 2 * i + 2; // right child
    // find if the children are smaller or not, if they are we need to swap them
    // BUT THE COMPARISON IS DONE IN THE eik_queue
    if (l < p_queueImp->size && p_queueImp->queue_vals[l] < p_queueImp->queue_vals[smallest])
      smallest = l;
    if (r < p_queueImp->size && p_queueImp->queue_vals[r] < p_queueImp->queue_vals[smallest])
      smallest = r;
    if (smallest != i)
    {
      swap_double(&p_queueImp->queue_vals[i], &p_queueImp->queue_vals[smallest]); // swap in the eik_queue
      swap_int(&p_queueImp->queue_index[i], &p_queueImp->queue_index[smallest]); // swap in the index_queue
      heapify(p_queueImp, smallest); // recurrencia
    }
  }
}

static void insert(p_queue *p_queueImp, double newNum, int newIndex)
{
  int i;
  if (p_queueImp->size == 0) // First time we insert an element in the tree
  {
    p_queueImp->size = 1;
    p_queueImp->queue_vals[0] = newNum;
    p_queueImp->queue_index[0] = newIndex;
  }
  else // if there were elements in the tree before inserting
  {
    p_queueImp->size += 1;
    if ( p_queueImp->size >= p_queueImp->maxSize  ) {
      grow_queue( p_queueImp );
    }
    p_queueImp->queue_vals[p_queueImp->size -1 ] = newNum;
    p_queueImp->queue_index[p_queueImp->size -1 ] = newIndex;
    for (i = p_queueImp->size / 2 ; i >= 0; i--)
    {
      heapify(p_queueImp, i);
    }
  }
}



static void insert_end(p_queue *p_queueImp, double newNum, int newIndex)
{
    p_queueImp->size += 1;
    p_queueImp->queue_vals[p_queueImp->size -1 ] = newNum;
    p_queueImp->queue_index[p_queueImp->size -1] = newIndex;
}



static void delete_findValue(p_queue *p_queueImp, double num)
{
  int i;
  for (i = 0; i < p_queueImp->size; i++)
  {
    if (num == p_queueImp->queue_vals[i] ) // find the value
      break;
  }

  swap_double(&p_queueImp->queue_vals[i], &p_queueImp->queue_vals[p_queueImp->size - 1]); 
  swap_int(&p_queueImp->queue_index[i], &p_queueImp->queue_index[p_queueImp->size - 1]); 

  p_queueImp->size -= 1; // make it "forget", after size-1 everything is ignored
  for (i = p_queueImp->size / 2 ; i >= 0; i--)
  {
    heapify(p_queueImp, i);
  }
}

static void delete_findIndex(p_queue *p_queueImp, int ind)
{
  int i;
  for (i = 0; i < p_queueImp->size; i++)
  {
    if (ind == p_queueImp->queue_index[i])
      break;
  }

  swap_double(&p_queueImp->queue_vals[i], &p_queueImp->queue_vals[p_queueImp->size - 1]); 
  swap_int(&p_queueImp->queue_index[i], &p_queueImp->queue_index[p_queueImp->size - 1]); 

  p_queueImp->size -= 1; // make it "forget", after size-1 everything is ignored
  for (i = p_queueImp->size / 2 ; i >= 0; i--)
  {
    heapify(p_queueImp, i);
  }
}


static void deleteRoot(p_queue *p_queueImp)
{
  int i;
  // we dont need to look for the index, we know its on the 0th position
  swap_double(&p_queueImp->queue_vals[0], &p_queueImp->queue_vals[p_queueImp->size - 1]); 
  swap_int(&p_queueImp->queue_index[0], &p_queueImp->queue_index[p_queueImp->size - 1]); 
  p_queueImp->size -= 1; 
  for (i = p_queueImp->size / 2; i >= 0; i--)
  {
    heapify(p_queueImp, i);
  }
}


static void printeik_queue(p_queue *p_queueImp)
{
  int i;
  printf("Eikonal values:");
  for (i = 0; i < p_queueImp->size; ++i)
    printf("%lf ", p_queueImp->queue_vals[i]);
  printf("\n");
  printf("Indices:");
  for (int i = 0; i < p_queueImp->size; ++i)
    printf("%d ", p_queueImp->queue_index[i]);
  printf("\n");
  printf("Current size: %d", p_queueImp->size);
}

static void update(p_queue *p_queueImp, double new_valConsidered, int index)
{
    int i, j;
    // First find the current value associated with index
    for (i = 0; i < p_queueImp->size; i++)
    {
        if (index == p_queueImp->queue_index[i])
        break;
    }
    // Then, if the new value considered is smaller than the current value
    if ( p_queueImp->queue_vals[i] > new_valConsidered  )
    {
        p_queueImp->queue_vals[i] = new_valConsidered;
        for (j = p_queueImp->size / 2 ; j >= 0; j--)
        {
            heapify(p_queueImp, j);
        }
    }
}



static double get_valueAtIndex(p_queue *p_queueImp, int index)
{
    int i;
    for (i = 0; i < p_queueImp->size; i ++)
    {
        if (index == p_queueImp->queue_index[i])
        break;
    }
    return p_queueImp->queue_vals[i];
}

