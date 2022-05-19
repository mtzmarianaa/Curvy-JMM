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

int size = 0; // Initial size is zero because we need to insert something to initialize the queue


typedef struct Priority_queue {
  double *queue_vals; // if we need more we'll add more
  int *queue_index; // same here
  int size; // current occupied size occupied
} p_queue;

void priority_queue_init( p_queue *p_queueImp  ) {
  p_queueImp->queue_vals = malloc( 16*sizeof(double)  );
  p_queueImp->queue_index = malloc( 16*sizeof(int) );
  p_queueImp->size = 0;
  assert( p_queueImp != NULL  ); // the queue should not be null if initialized
}

void priority_queue_deinit( p_queue *p_queueImp ) {
  free( p_queueImp->queue_vals );
  free(p_queueImp->queue_index  );
  p_queueImp->queue_vals = NULL;
  p_queueImp->queue_index = NULL;
}



void grow_queue( p_queue *p_queueImp ) {
  p_queueImp->queue_vals = realloc( p_queueImp->queue_vals, 2*sizeof(p_queueImp->queue_vals)  );
  p_queueImp->queue_index = realloc(p_queueImp->queue_index, 2*sizeof(p_queueImp->queue_index) );
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
  if (p_queueImp->size == 0) // First time we insert an element in the tree
  {
    p_queueImp->queue_vals[0] = newNum;
    p_queueImp->queue_index[0] = newIndex;
    p_queueImp->size += 1;
  }
  else // if there were elements in the tree before inserting
  {
    if ( p_queueImp->size >=16  ) {
      grow_queue( p_queueImp );
    }
    p_queueImp->queue_vals[p_queueImp->size] = newNum;
    p_queueImp->queue_index[p_queueImp->size] = newIndex;
    p_queueImp->size += 1;
    for (int i = p_queueImp->size / 2 - 1; i >= 0; i--)
    {
      heapify(p_queueImp, i);
    }
  }
}

/*

static void insert_end(p_queue *p_queueImp, double newNum, int newIndex)
{
    p_queueImp->queue_vals[   ]
    eik_queue[size] = newNum;
    index_queue[size] = newIndex;
    size += 1;
}



static void delete_findValue(double eik_queue[], int index_queue[], double num)
{
  int i;
  for (i = 0; i < size; i++)
  {
    if (num == eik_queue[i]) // find the value
      break;
  }

  swap_double(&eik_queue[i], &eik_queue[size - 1]);
  swap_int(&index_queue[i], &index_queue[size - 1]);
  size -= 1;
  for (int i = size / 2 - 1; i >= 0; i--)
  {
    heapify(eik_queue, index_queue, size, i);
  }
}

static void delete_findIndex(double eik_queue[], int index_queue[], int ind)
{
  int i;
  for (i = 0; i < size; i++)
  {
    if (ind == index_queue[i])
      break;
  }

  swap_double(&eik_queue[i], &eik_queue[size - 1]);
  swap_int(&index_queue[i], &index_queue[size - 1]);
  size -= 1;
  for (int i = size / 2 - 1; i >= 0; i--)
  {
    heapify(eik_queue, index_queue, size, i);
  }
}

static void deleteRoot(double eik_queue[], int index_queue[])
{

  swap_double(&eik_queue[0], &eik_queue[size - 1]);
  swap_int(&index_queue[0], &index_queue[size - 1]);
  size -= 1;
  for (int i = size / 2 - 1; i >= 0; i--)
  {
    heapify(eik_queue, index_queue, size, i);
  }
}

static void printeik_queue(double eik_queue[], int index_queue[], int size)
{
  printf("Eikonal values:");
  for (int i = 0; i < size; ++i)
    printf("%lf ", eik_queue[i]);
  printf("\n");
  printf("Indices:");
  for (int i = 0; i < size; ++i)
    printf("%d ", index_queue[i]);
  printf("\n");
}

static void update(double eik_queue[], int index_queue[], double new_valConsidered, int index)
{
    // First find the current value associated with index
    int i;
    for (i = 0; i < size; i++)
    {
        if (index == index_queue[i])
        break;
    }
    // Then, if the new value considered is smaller than the current value
    if ( eik_queue[i] > new_valConsidered  )
    {
        eik_queue[i] = new_valConsidered;
        for (int j = size / 2 - 1; j >= 0; j--)
        {
            heapify(eik_queue, index_queue, size, j);
        }
    }
}

static double get_valueAtIndex(double eik_queue[], int index_queue[], int index, int size)
{
    int i;
    for (i = 0; i < size; i ++)
    {
        if (index == index_queue[i])
        break;
    }
    return eik_queue[i];
}

*/