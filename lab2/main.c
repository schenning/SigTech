#include <stdio.h>
#include "stdlib.h"
#include "smmintrin.h"
#include "emmintrin.h"
#include "time.h"
#include "time_meas.h"
#include "math.h"

void m16_vv(int16_t *x, int16_t *y, int16_t *z, int N){
	__m128i *x128, *y128, *z128;
	x128 = (__m128i *)x;
	y128 = (__m128i *)y;
	z128 = (__m128i *)z;
	int i;
	for(i=0;i<(N>>3); i++){
		z128[i] =  _mm_slli_epi16( _mm_mulhi_epi16(x128[i], y128[i]),1);
	}
}





void m16_vv_mulhrs(int16_t *x, int16_t *y, int16_t *z, int N){
	__m128i *x128, *y128, *z128;
	x128 = (__m128i *)x;
	y128 = (__m128i *)y;
	z128 = (__m128i *)z;
	int i;
	for(i=0;i<(N>>3); i++){
		z128[i] =  _mm_slli_epi16( _mm_mulhrs_epi16(x128[i], y128[i]),1);
	}
}





int main() {


	int16_t x[2048]; int16_t y[2048]; int16_t z[2048];
	int i,j;


	time_t t; 
	/* Intializes random number generator */
    srand((unsigned) time(&t));
	for( j = 0 ; j < 2048 ; j++ ){ 
		x[j] = rand()%1000;
		y[j] = rand()%1000; 
		z[j] = 0;
   	}

	printf("x = [ ");
	for( j = 0 ; j < 2048 ; j++ ){ 
		printf("%d, ", x[j]);
	
	}
	printf(" ] \n");
	printf("y = [ ");

	for( j = 0 ; j < 2048 ; j++ ){ 
		printf("%d, ", y[j]);	
	}
	printf("]\n");

	m16_vv(x,y,z,2048);
	printf("z = [ ");
	for (i=0;i<16; i++){
	printf("%d, ", z[i]);
    }

	printf(" ]\n");

	printf("should be ");
	
	printf(" = [ ");
	for (i=0;i<16; i++){
		printf("%d, ", (int)((x[i]*y[i])/pow(2,16)) *2  );	
	}
	printf(" ]\n");
	time_stats_t ts;
	reset_meas (&ts);
	start_meas (&ts);

	unsigned long long k;
	for(k=0;k< 4000000; k++){
		m16_vv(x,y,z,2048);

	}
	stop_meas (&ts);
	printf("difference:  %lld \n", ts.diff);

return 0;
}


