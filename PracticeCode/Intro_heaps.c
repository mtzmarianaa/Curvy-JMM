/* Learning C

             HEAPS

Heap data structure is a complete binary tree that satisfies the heap property, where any given node is
     - always greater than its child node/s and the key of the root 
       node is the largest among all other nodes. This property is also called max heap property 
     - always smaller than the child node/s and the key of the root 
       node is the smallest among all other nodes. This property is also called min heap property.

Here we code different important operations performed on a heap.

     - Heapify is the process of creating a heap data structure from a binary tree. 
       It is used to create a Min-Heap or a Max-Heap.

*/


#include <stdio.h>
int size = 0;
void swap(int *a, int *b)
{
  int temp = *b;
  *b = *a;
  *a = temp;
}

void heapify(int array[], int size, int i)
{
  if (size == 1)
  {
    printf("Single element in the heap");
  }
  else
  {
    int smallest = i;
    int l = 2 * i + 1; // left child
    int r = 2 * i + 2; // right child
    if (l < size && array[l] < array[smallest])
      smallest = l;
    if (r < size && array[r] < array[smallest])
      smallest = r;
    if (smallest != i)
    {
      swap(&array[i], &array[smallest]);
      heapify(array, size, smallest); // recurencia
    }
  }
}

void insert(int array[], int newNum)
{
  if (size == 0)
  {
    array[0] = newNum;
    size += 1;
  }
  else
  {
    array[size] = newNum;
    size += 1;
    for (int i = size / 2 - 1; i >= 0; i--)
    {
      heapify(array, size, i);
    }
  }
}

void delete_findNumber(int array[], int num)
{
  int i;
  for (i = 0; i < size; i++)
  {
    if (num == array[i])
      break;
  }

  swap(&array[i], &array[size - 1]);
  size -= 1;
  for (int i = size / 2 - 1; i >= 0; i--)
  {
    heapify(array, size, i);
  }
}

void delete_index(int array[], int num)
{
  int i;
  for (i = 0; i < size; i++)
  {
    if (num == array[i])
      break;
  }

  swap(&array[i], &array[size - 1]);
  size -= 1;
  for (int i = size / 2 - 1; i >= 0; i--)
  {
    heapify(array, size, i);
  }
}

void deleteRoot(int array[])
{

  swap(&array[0], &array[size - 1]);
  size -= 1;
  for (int i = size / 2 - 1; i >= 0; i--)
  {
    heapify(array, size, i);
  }
}

void printArray(int array[], int size)
{
  for (int i = 0; i < size; ++i)
    printf("%d ", array[i]);
  printf("\n");
}

int main()
{
  int array[10];

  insert(array, 3);
  printArray(array, size);
  insert(array, 1);
  printArray(array, size);
  insert(array, 9);
  printArray(array, size);
  insert(array, 5);
  printArray(array, size);
  insert(array, 2);
  printArray(array, size);
  insert(array, 6);
  printArray(array, size);
  insert(array, 90);
  printArray(array, size);
  insert(array, 80);
  printArray(array, size);
  insert(array, 7);
  printArray(array, size);
  insert(array, 10);

  printf("Min-Heap array: ");
  printArray(array, size);

  delete(array, 3);

  printf("After deleting an element: ");

  printArray(array, size);

  delete(array, 99);

  printf("After deleting the last one: ");

  printArray(array, size);

  deleteRoot(array);

  printf("After deleting the root: ");

  printArray(array, size);

}