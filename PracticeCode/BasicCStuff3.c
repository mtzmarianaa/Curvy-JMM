/* Learning C

             FUNCTIONS

The "template" for user define functions is the following:

#include <stdio.h>
void functionName()
{
    ... .. ...
    ... .. ...
}

int main()
{
    ... .. ...
    ... .. ...

    functionName();
    
    ... .. ...
    ... .. ...
}


The execution of a C program begins from the main() function.
When the compiler encounters functionName();, control of the program jumps to void functionName()
And, the compiler starts executing the codes inside functionName().
The control of the program jumps back to the main() function once code inside the function definition is executed.

A function prototype is simply the declaration of a function that specifies function's name, parameters and return type. It doesn't contain function body.
A function prototype gives information to the compiler that the function may later be used in the program. The syntax of the function prototype is:

     returnType functionName(type1 argument1, type2 argument2, ...);

The function prototype is not needed if the user-defined function is defined before the main() function.

Control of the program is transferred to the user-defined function by calling it. The syntax of the function call is:
     
     functionName(argument1, argument2, ...);

Function definition contains the block of code to perform a specific task. When a function is called, the control of the program is transferred 
to the function definition. And, the compiler starts executing the codes inside the body of a function. The syntax of the function definition is:

     returnType functionName(type1 argument1, type2 argument2, ...)
     {
         //body of the function
     }

The type of arguments passed to a function and the formal parameters must match, otherwise, the compiler will throw an error. The return statement 
terminates the execution of a function and returns a value to the calling function. The program control is transferred to the calling function 
after the return statement.



*/

#include <stdio.h>
int checkPrimeNumber(int n);

int main() {

  int n, flag;

  printf("Enter a positive integer: ");
  scanf("%d",&n);

  // n is passed to the checkPrimeNumber() function
  // the returned value is assigned to the flag variable
  flag = checkPrimeNumber(n);

  if(flag == 1)
    printf("%d is not a prime number",n);
  else
    printf("%d is a prime number",n);

  return 0;
}

// int is returned from the function
int checkPrimeNumber(int n) {

  // 0 and 1 are not prime numbers    
  if (n == 0 || n == 1)
    return 1;

  int i;

  for(i=2; i <= n/2; ++i) {
    if(n%i == 0)
      return 1;
  }

  return 0;
}