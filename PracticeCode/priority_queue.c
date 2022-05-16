/* PRIORITY QUEUE

Here is the implementation of the prioriry queue for the 2D FMM. Explanaition

   - eik_queue: Eikonal values considered, binary tree in an array, this is the thing to heapify
   - index_queue : index (grid coordinates) if the Eikonal values considered in the queue
   - size : current non empty entries of the binary tree

*/


#include <stdio.h>

int size = 0; // Initial size is zero because we need to insert something to initialize the queue

/*
typedef struct Priority_queue {
  int queue_vals[10];
  int queue_index[10];
} p_queue;
*/

void swap(int *a, int *b)
{
  int temp = *b;
  *b = *a;
  *a = temp;
}

void heapify(int eik_queue[], int index_queue[], int size, int i)
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
      swap(&eik_queue[i], &eik_queue[smallest]); // swap in the eik_queue
      swap(&index_queue[i], &index_queue[smallest]); // swap in the index_queue
      heapify(eik_queue, index_queue, size, smallest); // recurrencia
    }
  }
}

void insert(int eik_queue[], int index_queue[], int newNum, int newIndex)
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

void insert_end(int eik_queue[], int index_queue[], int newNum, int newIndex)
{
    eik_queue[size] = newNum;
    index_queue[size] = newIndex;
    size += 1;
}

void delete_findValue(int eik_queue[], int index_queue[], int num)
{
  int i;
  for (i = 0; i < size; i++)
  {
    if (num == eik_queue[i]) // find the value
      break;
  }

  swap(&eik_queue[i], &eik_queue[size - 1]);
  swap(&index_queue[i], &index_queue[size - 1]);
  size -= 1;
  for (int i = size / 2 - 1; i >= 0; i--)
  {
    heapify(eik_queue, index_queue, size, i);
  }
}

void delete_findIndex(int eik_queue[], int index_queue[], int ind)
{
  int i;
  for (i = 0; i < size; i++)
  {
    if (ind == index_queue[i])
      break;
  }

  swap(&eik_queue[i], &eik_queue[size - 1]);
  swap(&index_queue[i], &index_queue[size - 1]);
  size -= 1;
  for (int i = size / 2 - 1; i >= 0; i--)
  {
    heapify(eik_queue, index_queue, size, i);
  }
}

void deleteRoot(int eik_queue[], int index_queue[])
{

  swap(&eik_queue[0], &eik_queue[size - 1]);
  swap(&index_queue[0], &index_queue[size - 1]);
  size -= 1;
  for (int i = size / 2 - 1; i >= 0; i--)
  {
    heapify(eik_queue, index_queue, size, i);
  }
}

void printeik_queue(int eik_queue[], int index_queue[], int size)
{
  printf("Eikonal values:");
  for (int i = 0; i < size; ++i)
    printf("%d ", eik_queue[i]);
  printf("\n");
  printf("Indices:");
  for (int i = 0; i < size; ++i)
    printf("%d ", index_queue[i]);
  printf("\n");
}

int main()
{
  int eik_queue[10], index_queue[10];

  insert(eik_queue, index_queue, 3, 0);
  printeik_queue(eik_queue, index_queue, size);
  insert(eik_queue, index_queue, 1, 1);
  printeik_queue(eik_queue, index_queue, size);
  insert(eik_queue, index_queue, 9, 2);
  printeik_queue(eik_queue, index_queue, size);
  insert(eik_queue, index_queue, 5, 3);
  printeik_queue(eik_queue, index_queue, size);
  insert(eik_queue, index_queue, 2, 4);
  printeik_queue(eik_queue, index_queue, size);
  insert(eik_queue, index_queue, 6, 5);
  printeik_queue(eik_queue, index_queue, size);
  insert(eik_queue,index_queue, 90, 6);
  printeik_queue(eik_queue,index_queue, size);
  insert(eik_queue,index_queue, 80, 7);
  printeik_queue(eik_queue,index_queue, size);
  insert(eik_queue, index_queue, 7, 8);
  printeik_queue(eik_queue, index_queue, size);
  insert(eik_queue, index_queue, 10, 9);

  printf("Min-Heap eik_queue: ");
  printeik_queue(eik_queue, index_queue, size);

  delete_findValue(eik_queue, index_queue, 3);

  printf("After deleting an element: ");

  printeik_queue(eik_queue, index_queue, size);

  delete_findValue(eik_queue, index_queue, 99);

  printf("After deleting the last one: ");

  printeik_queue(eik_queue, index_queue, size);

  deleteRoot(eik_queue,index_queue);

  printf("After deleting the root: ");

  printeik_queue(eik_queue, index_queue, size);

  insert_end(eik_queue, index_queue, 1000, 10);

  printf("After inserting a number at the end of the queue: ");

  printeik_queue(eik_queue, index_queue, size);

  delete_findIndex(eik_queue, index_queue, 10);

  printf("After deleating the value associated with the 10th index: ");

  printeik_queue(eik_queue, index_queue, size);

}