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
#include <math.h>

void swap(int *n1, int *n2);
void heapify(int array[], int size, int i);
void insert(int array[], int newNum, int size);
void deleteRoot(int array[], int num, int size);
void printArray(int array[], int size);

int main()
{
    int initial_i, size;
    size = 6;
    int array[size];
    initial_i = 0;
    array[0] = 1;
    array[1] = 99;
    array[2] = 3;
    array[3] = 40;
    array[4] = 5;
    array[5] = 10;


    printf("Original array: \n");
    printArray(array, size);

    heapify(array, size, initial_i);

    printf("Heapified array: \n");
    printArray(array, size);

  
  return 0;

}

void swap(int* n1, int* n2)
{
    /*
    Exchanges two pointers. 
    */
    int temp;
    temp = *n1;
    *n1 = *n2;
    *n2 = temp;
}

void heapify(int array[], int size, int i)
{
    /*
    Creates a heap like structure from data in an array form. RECURSION
    */
  if (size == 1)
  {
    printf("Single element in the heap");
  }
  else
  {
    int largest = i;
    int l = 2 * i + 1;
    int r = 2 * i + 2;
    if (l < size && array[l] > array[largest])
      largest = l;
    if (r < size && array[r] > array[largest])
      largest = r;
    if (largest != i)
    {
      swap(&array[i], &array[largest]);
      heapify(array, size, largest);
    }
  }
}


void insert(int array[], int newNum, int size)
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



void deleteRoot(int array[], int num, int size)
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



void printArray(int array[], int size)
{
  for (int i = 0; i < size; ++i)
    printf("Index, %d, contains: %d  \n", i, array[i]);
  printf("\n");
}