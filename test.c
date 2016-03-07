#include <stdio.h>
#include <math.h>
#include "rangen_double.h"



int main(){

double x,i;

for (i=0;i<20;i++){
x = uniformrandom ();

printf("%f\n", x);
}



return 0;
} 

