/*
TESTS FOR THE PRIORITY QUEUE
*/

#include "priority_queue.c"

int main()
{
  double eik_queue[10];
  int index_queue[10], ind_found;

  insert(eik_queue, index_queue, 3.0, 0);
  printeik_queue(eik_queue, index_queue, size);
  insert(eik_queue, index_queue, 1.0, 1);
  printeik_queue(eik_queue, index_queue, size);
  insert(eik_queue, index_queue, 9.0, 2);
  printeik_queue(eik_queue, index_queue, size);
  insert(eik_queue, index_queue, 5.0, 3);
  printeik_queue(eik_queue, index_queue, size);
  insert(eik_queue, index_queue, 2.0, 4);
  printeik_queue(eik_queue, index_queue, size);
  insert(eik_queue, index_queue, 6.0, 5);
  printeik_queue(eik_queue, index_queue, size);
  insert(eik_queue,index_queue, 90.0, 6);
  printeik_queue(eik_queue,index_queue, size);
  insert(eik_queue,index_queue, 80.0, 7);
  printeik_queue(eik_queue,index_queue, size);
  insert(eik_queue, index_queue, 7.0, 8);
  printeik_queue(eik_queue, index_queue, size);
  insert(eik_queue, index_queue, 10.0, 9);

  printf("Min-Heap eik_queue: ");
  printeik_queue(eik_queue, index_queue, size);

  delete_findValue(eik_queue, index_queue, 3);

  printf("After deleting an element: ");

  printeik_queue(eik_queue, index_queue, size);

  delete_findValue(eik_queue, index_queue, 99.0);

  printf("After deleting the last one: ");

  printeik_queue(eik_queue, index_queue, size);

  deleteRoot(eik_queue,index_queue);

  printf("After deleting the root: ");

  printeik_queue(eik_queue, index_queue, size);

  insert_end(eik_queue, index_queue, 1000.0, 10);

  printf("After inserting a number at the end of the queue: ");

  printeik_queue(eik_queue, index_queue, size);

  delete_findIndex(eik_queue, index_queue, 10);

  printf("After deleating the value associated with the 10th index: ");

  printeik_queue(eik_queue, index_queue, size);

  update(eik_queue, index_queue, 70.0, 9);

  printf("After comparing the 9th index with a greater value: ");

  printeik_queue(eik_queue, index_queue, size);

  update(eik_queue, index_queue, 4.0, 9);

  printf("After comparing the 9th index with a smaller value: ");

  printeik_queue(eik_queue, index_queue, size);

  printf("Value of index 3 ");

  ind_found = get_valueAtIndex(eik_queue, index_queue, 3, size);

  printf("%d ", ind_found);

}