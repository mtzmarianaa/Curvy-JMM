/* Learning C

             POINTERS

If you have a variable var in your program, &var will give you its address in the memory.

Pointers (pointer variables) are special variables that are used to store addresses rather than values.



*/

#include <stdio.h>

void swap(int *n1, int *n2);

int main()
{
  int var = 5;
  printf("var: %d\n", var);

  // Notice the use of & before var
  printf("address of var: %p \n \n", &var); 

  // But if we declare a pointer


  int* pc, c;
  c = 5;
  pc = &c;
  printf("First part, int %d \n", c);  
  printf("First part, pointer %d \n", *pc);
  c = 1;
  printf("Second part, int %d \n", c);  
  printf("Second part pointer %d \n", *pc);
  printf("Second part pointer %p \n \n", pc);


  int num1 = 5, num2 = 10;
  printf("Pre swap num1 = %d\n", num1);
  printf("Pre swap num2 = %d \n", num2);
  printf("Pre swap address num1 = %p\n", &num1);
  printf("Pre swap address num2 = %p \n", & num2);
  // address of num1 and num2 is passed
  swap( &num1, &num2);
  
  printf("num1 = %d\n", num1);
  printf("num2 = %d \n", num2);
  printf("Address num1 = %p\n", &num1);
  printf("Address num2 = %p \n", & num2);

  return 0;

}

void swap(int* n1, int* n2)
{
    int temp;
    temp = *n1;
    *n1 = *n2;
    *n2 = temp;
}
