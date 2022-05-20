/*
TESTS FOR THE PRIORITY QUEUE
*/

#include "priority_queue.c"

int main()
{
  
  p_queue *p_queueImp = malloc(sizeof(p_queue));

  priority_queue_init( p_queueImp  );

  printeik_queue(p_queueImp);

  printf("%d", p_queueImp->size);

  printf("\n--------");

  printf("Insert value \n");

  insert(p_queueImp, 3.0, 0);

  printeik_queue(p_queueImp);

  printf("\n-----");

  printf("\n Insert value \n");

  insert(p_queueImp, 2.0, 1);

  printeik_queue(p_queueImp);

  printf("\n-----");

  printf("\n Insert value \n");

  insert(p_queueImp, 12.5, 2);

  printeik_queue(p_queueImp);

  printf("\n--------");

  printf("\n Insert value \n");

  insert(p_queueImp, 5.0, 3);

  printeik_queue(p_queueImp);

  printf("\n--------");

  printf("\n Insert value \n");

  insert(p_queueImp, 1.1, 4);

  printeik_queue(p_queueImp);

  printf("\n--------");

  printf("\n Insert value \n");

  insert(p_queueImp, 1.0, 5);

  printeik_queue(p_queueImp);

  printf("\n--------");

  printf("\n Insert value \n");

  insert(p_queueImp, 4.0, 6);

  printeik_queue(p_queueImp);

  printf("\n--------");

  printf("\n Insert value \n");

  insert(p_queueImp, 1.2, 7);

  printeik_queue(p_queueImp);

  printf("\n--------");

  printf("\n Delete fourth element inserted, %fl \n", get_valueAtIndex(p_queueImp, 4));

  delete_findIndex(p_queueImp, 4);

  printeik_queue(p_queueImp);

  printf("\n--------");

  printf("\n Delete root \n");

  deleteRoot(p_queueImp);

  printeik_queue(p_queueImp);

  printf("\n--------");

  printf("\n Memory allocation of the whole struct \n");

  printf("%d", &p_queueImp);

  printf("\n--------");

  printf("\n Memory allocation of the list of values of the struct \n");

  printf("%d", &p_queueImp->queue_vals);

  printf("\n--------");

  printf("\n Memory allocation of the list of indices of the struct \n");

  printf("%d", &p_queueImp->queue_index);

  priority_queue_deinit(p_queueImp);

}