/*
TESTS FOR THE PRIORITY QUEUE
*/

#include "priority_queue.c"

int main()
{
  
  p_queue p_queueImp;

  priority_queue_init( &p_queueImp  );

  printf("%d", p_queueImp.size);

  grow_queue( &p_queueImp );

  insert(&p_queueImp, 3.0, 0);

  printf("\n");

  printf("%d", p_queueImp.size);

  printf("\n");

  printf("%d", p_queueImp.queue_index[0]);
  
  printf("\n");

  printf("%lf", p_queueImp.queue_vals[0]);

}