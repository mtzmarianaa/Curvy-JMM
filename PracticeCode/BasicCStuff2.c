/* LEARNING C

             FLOW CONTROL

IF STATEMENT

if (test expression) 
{
   // code
}



IF ELSE STATEMENT

if (test expression) 
{
    // run code if test expression is true
}
else 
{
    // run code if test expression is false
}



ELIF

if (test expression1) 
{
   // statement(s)
}
else if(test expression2) 
{
   // statement(s)
}
else if (test expression3) 
{
   // statement(s)
}
.
.
else 
{
   // statement(s)
}


FOR LOOP

for (initializationStatement; testExpression; updateStatement)
{
    // statements inside the body of loop
}


WHILE LOOP

while (testExpression) {
  // the body of the loop 
}

DO WHILE LOOP

do {
  // the body of the loop
}
while (testExpression);




*/
#include <stdio.h>
int main() {

    // example if
    int a = 4;
    int b = 10;
    printf("If loop example \n");
    if (a > b) {
    printf("Hello\n");
    }
    printf("Hi\n");

    // example for

    int i;

    for (i = 1; i < 11; ++i)
    {
        printf("%d\n", i);
    }

    printf("While loop example\n");

    int j = 1;
    while (j < 10)
    {
        printf("%d\n", j);
        ++j;
    }

    printf("Do while example \n");

    int k = 1;
    do
    {
        ++k;
        printf("%d", k);
    } 
    while (2*k < 2);
    


    return 0;
}