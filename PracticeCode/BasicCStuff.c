/* LEARNING C

             INTRODUCTION

Identifier refers to name given to entities such as variables, functions, structures etc.
Identifiers must be unique. They are created to give a unique name to an entity to identify 
it during the execution of the program.

Rules for naming variables:
  - A variable name can only have letters (both uppercase and lowercase letters), digits and underscore.
   - The first letter of a variable should be either a letter or an underscore.
   - There is no rule on how long a variable name (identifier) can be.
   - This means that the variable type cannot be changed once it is declared.

 Literals are data used for representing fixed values. They can be used directly in the code.


 Integers:
 An integer is a numeric literal(associated with numbers) without any fractional or exponential part. 
 There are three types of integer literals in C programming: decimal (base 10), octal (base 8), hexadecimal (base 16).
 In C programming, octal starts with a 0, and hexadecimal starts with a 0x.

 Escape sequences:
 Sometimes, it is necessary to use characters that cannot be typed or has special meaning in C programming. 
 For example: newline(enter), tab, question mark etc. In order to use these characters, escape sequences are used.
   - \b backspace
   - \f form feed
   - \n new line
   - \r return
   - \t horizontal tab
   - \v vertical tab
   - \\ backslash
   - \' single quotation mark
   - \" double quotation mark
   - \? question mark
   - \0 null character

 Constants:
 If you want to define a variable whose value cannot be changed, you can use the const keyword. This will create a constant.

 DATA TYPES
  Type           size(bytes)           format specifies
   int         at least 2, usually 4    %d, % i
   char              4                    %c
  float              8                    %f
  double             8                    %lf
   short int         2                    %hd
 unsigned int        2                    %u
 long int         4 or 8                  %ld, %li
 long long int   at least 8               %lld, %lli
 unsigned long int   8                    %lld, %lli
 signed char         1                    %c
  long double     10/12/16                %Lf
 
 The size of float (single precision float data type) is 4 bytes. And the size of double (double precision float data type) is 8 bytes.

 void is an incomplete type. It means "nothing" or "no type". 
 You can think of void as absent. For example, if a function is not 
 returning anything, its return type should be void. Note that, you 
 cannot create variables of void type.
 Tip: You can always check the size of a variable using the sizeof() operator.

 signed - allows for storage of both positive and negative numbers
 unsigned - allows for storage of only positive numbers
 
 Data types that are derived from fundamental data types are derived types. For example: arrays, pointers, function types, structures, etc.

 HEADER
 stdio.h is a header file which has the necessary information to include the input/output related functions in our program. Example printf, scanf etc.
 If we use #include<stdio.h> in your c program, it will include stdio.h file into our source program which has the information for all input, 
 output related functions. Thus a very simple program should look like:
        #include<stdio.h>       header section
        
        int main()             main section
        {
            printf("Hello World");
            return 0;
        }

MAIN FUNCTION
A C program starts with a main() function, usually kept in a file named main.c.
The main() function is the first function in your program that is executed when it begins executing, but it's not the first function executed. 
The first function is _start(), which is typically provided by the C runtime library, linked in automatically when your program is compiled. 
The details are highly dependent on the operating system and compiler toolchain, so I'm going to pretend I didn't mention it.
The main() function has two arguments that traditionally are called argc and argv and return a signed integer. Most Unix environments expect 
programs to return 0 (zero) on success and -1 (negative one) on failure.

   - argc is the argument count, its the length of the argument vectos
   - argv is the argument vector, its the array of character pointers, is a tokenized representation of the command line that invoked your program


LOGICAL OPERATORS

   - && Logical AND. True only if all operands are true
   - || Logical OR. True only if either one operand is true
   - !  Logical NOT. True only if the operand is 0

We can also have bitwise operators
   - & bitwise and
   - | bitwise or
   - ^ bitwise exclusive or
   - ~ bitwise complement
   - << shift left
   - >> shift right


 Example: */
#include <stdio.h>
int main()
{
    char chr = 'm';    
    printf("character = %c", chr);  
    return 0;
} 