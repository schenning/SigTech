#include <stdio.h>
#include <stdlib.h>
#include "smmintrin.h"
#include "emmintrin.h"

void m16_vv(int16_t *x, int16_t *y, int16_t *z, int N){
	__m128i *x128, *y128, *z128;
	x128 = (__m128i *)x;
	y128 = (__m128i *)y;
	z128 = (__m128i *)z;
	int i;
	for(i=0;i<(N>>3); i++){
		z128[i] = _mm_mulhi_epi16(x128[i], y128[i]);
	}
}

int main() {


	int16_t x[16]; int16_t y[16]; int16_t z[16];
	int i,j;


	time_t t; 
	/* Intializes random number generator */
    srand((unsigned) time(&t));
	for( j = 0 ; j < 16 ; j++ ){ 
		x[i] = rand()%2;
		y[i] = rand()%2; 
		z[i] = 0;
   	}


	unsigned long long k;
	for(k=0;k< 4294967295; k++){
		

	}

return 0;
}


