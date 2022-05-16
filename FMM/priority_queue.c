/* PRIORITY QUEUE

Here is the implementation of the prioriry queue for the 2D FMM. Explanaition

   - eik_queue: Eikonal values considered, binary tree in an array, this is the thing to heapify
   - index_queue : index (grid coordinates) if the Eikonal values considered in the queue
   - size : current non empty entries of the binary tree

*/

#include "priority_queue.h"
#include <stdio.h>

int size = 0; // Initial size is zero because we need to insert something to initialize the queue

/*
typedef struct Priority_queue {
  int queue_vals[10];
  int queue_index[10];
} p_queue;
*/

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

static void heapify(double eik_queue[], int index_queue[], int size, int i)
{
  if (size == 1)
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
    if (l < size && eik_queue[l] < eik_queue[smallest])
      smallest = l;
    if (r < size && eik_queue[r] < eik_queue[smallest])
      smallest = r;
    if (smallest != i)
    {
      swap_double(&eik_queue[i], &eik_queue[smallest]); // swap in the eik_queue
      swap_int(&index_queue[i], &index_queue[smallest]); // swap in the index_queue
      heapify(eik_queue, index_queue, size, smallest); // recurrencia
    }
  }
}

static void insert(double eik_queue[], int index_queue[], double newNum, int newIndex)
{
  if (size == 0) // First time we insert an element in the tree
  {
    eik_queue[0] = newNum;
    index_queue[0] = newIndex;
    size += 1;
  }
  else // if there were elements in the tree before inserting
  {
    eik_queue[size] = newNum;
    index_queue[size] = newIndex;
    size += 1;
    for (int i = size / 2 - 1; i >= 0; i--)
    {
      heapify(eik_queue, index_queue, size, i);
    }
  }
}

static void insert_end(double eik_queue[], int index_queue[], double newNum, int newIndex)
{
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

