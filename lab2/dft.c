#include "emmintrin.h"
#include "xmmintrin.h"



 __m128i cmac_tmp,cmac_tmp_re32,cmac_tmp_im32;

short reflip[8] = {1,-1,1,-1,1,-1,1,-1};

inline void cmac(__m128i a,__m128i b, __m128i *re32, __m128i *im32) {

  cmac_tmp    = _mm_sign_epi16(b,*(__m128i*)reflip);
  cmac_tmp_re32  = _mm_madd_epi16(a,cmac_tmp);


  cmac_tmp    = _mm_shufflelo_epi16(b,_MM_SHUFFLE(2,3,0,1));
  cmac_tmp    = _mm_shufflehi_epi16(cmac_tmp,_MM_SHUFFLE(2,3,0,1));

  cmac_tmp_im32  = _mm_madd_epi16(cmac_tmp,a);

  *re32 = _mm_add_epi32(*re32,cmac_tmp_re32);
  *im32 = _mm_add_epi32(*im32,cmac_tmp_im32);
}


__m128i mmtmpb,cre,cim;

inline void cmult(__m128i a,__m128i b, __m128i *re32, __m128i *im32) {

  mmtmpb    = _mm_sign_epi16(b,*(__m128i*)reflip);
  *re32     = _mm_madd_epi16(a,mmtmpb);
  mmtmpb    = _mm_shufflelo_epi16(b,_MM_SHUFFLE(2,3,0,1));
  mmtmpb    = _mm_shufflehi_epi16(mmtmpb,_MM_SHUFFLE(2,3,0,1));
  *im32  = _mm_madd_epi16(a,mmtmpb);

}

__m128i cpack_tmp1,cpack_tmp2;
inline __m128i cpack(__m128i xre,__m128i xim) {

  cpack_tmp1 = _mm_unpacklo_epi32(xre,xim);
  cpack_tmp1 = _mm_srai_epi32(cpack_tmp1,15);
  cpack_tmp2 = _mm_unpackhi_epi32(xre,xim);
  cpack_tmp2 = _mm_srai_epi32(cpack_tmp2,15);
  return(_mm_packs_epi32(cpack_tmp1,cpack_tmp2));

}

__m128i cre,cim;
inline void packed_cmult(__m128i a,__m128i b, __m128i *c) {

  cmult(a,b,&cre,&cim);
  *c = cpack(cre,cim);

}


short conjugatedft[8]__attribute__((aligned(16))) = {-1,1,-1,1,-1,1,-1,1} ;

// This macro performs 4 4-point radix-4 butterflies in parallel, input in x0,x1,x2,x3 and output in y0,y1,y2,y3

__m128i x1_flip,x3_flip,x3_2;
#define bfly4_tw1(x0,x1,x2,x3,y0,y1,y2,y3)  *(y0) = _mm_adds_epi16(*(x0),_mm_adds_epi16(*(x1),_mm_adds_epi16(*(x2),*(x3)))); \
  x1_flip = _mm_sign_epi16(*(x1),*(__m128i*)conjugatedft);\
  x1_flip = _mm_shufflelo_epi16(x1_flip,_MM_SHUFFLE(2,3,0,1));\
  x1_flip = _mm_shufflehi_epi16(x1_flip,_MM_SHUFFLE(2,3,0,1));\
  x3_flip = _mm_sign_epi16(*(x3),*(__m128i*)conjugatedft);\
  x3_flip = _mm_shufflelo_epi16(x3_flip,_MM_SHUFFLE(2,3,0,1));\
  x3_flip = _mm_shufflehi_epi16(x3_flip,_MM_SHUFFLE(2,3,0,1));\
  *(y1)   = _mm_adds_epi16(*(x0),_mm_subs_epi16(x1_flip,_mm_adds_epi16(*(x2),x3_flip)));\
  *(y2)   = _mm_subs_epi16(*(x0),_mm_subs_epi16(*(x1),_mm_subs_epi16(*(x2),*(x3))));\
  *(y3)   = _mm_subs_epi16(*(x0),_mm_adds_epi16(x1_flip,_mm_subs_epi16(*(x2),x3_flip)));


// This macro performs 4 4-point DFTs in parallel, including twiddle factors
__m128i x0r_2,x0i_2,x1r_2,x1i_2,dy0r,dy1r,dy0i,dy1i;
__m128i x2r_2,x2i_2,x3r_2,x3i_2,dy2r,dy2i,dy3r,dy3i;
__m128i x1_2,x2_2;

#define  bfly4(x0,x1,x2,x3,y0,y1,y2,y3,tw1,tw2,tw3)   cmult(*(x0),*(W0),&x0r_2,&x0i_2);\
  cmult(*(x1),*(tw1),&x1r_2,&x1i_2);\
  cmult(*(x2),*(tw2),&x2r_2,&x2i_2);\
  cmult(*(x3),*(tw3),&x3r_2,&x3i_2);\
  dy0r = _mm_add_epi32(x0r_2,_mm_add_epi32(x1r_2,_mm_add_epi32(x2r_2,x3r_2)));\
  dy0i = _mm_add_epi32(x0i_2,_mm_add_epi32(x1i_2,_mm_add_epi32(x2i_2,x3i_2)));\
  *(y0)  = cpack(dy0r,dy0i);\
  dy1r = _mm_add_epi32(x0r_2,_mm_sub_epi32(x1i_2,_mm_add_epi32(x2r_2,x3i_2)));\
  dy1i = _mm_sub_epi32(x0i_2,_mm_add_epi32(x1r_2,_mm_sub_epi32(x2i_2,x3r_2)));\
  *(y1)  = cpack(dy1r,dy1i);\
  dy2r = _mm_sub_epi32(x0r_2,_mm_sub_epi32(x1r_2,_mm_sub_epi32(x2r_2,x3r_2)));\
  dy2i = _mm_sub_epi32(x0i_2,_mm_sub_epi32(x1i_2,_mm_sub_epi32(x2i_2,x3i_2)));\
  *(y2)  = cpack(dy2r,dy2i);\
  dy3r = _mm_sub_epi32(x0r_2,_mm_add_epi32(x1i_2,_mm_sub_epi32(x2r_2,x3i_2)));\
  dy3i = _mm_add_epi32(x0i_2,_mm_sub_epi32(x1r_2,_mm_add_epi32(x2i_2,x3r_2)));\
  *(y3) = cpack(dy3r,dy3i);

#define transpose16(x0,x1,x2,x3) x0_tmp = _mm_unpacklo_epi32(x0,x1);\
   x1_tmp = _mm_unpackhi_epi32(x0,x1);\
   x2_tmp = _mm_unpacklo_epi32(x2,x3);\
   x3_tmp = _mm_unpackhi_epi32(x2,x3);\
   x0 = _mm_unpacklo_epi64(x0_tmp,x2_tmp);\
   x1 = _mm_unpackhi_epi64(x0_tmp,x2_tmp);\
   x2 = _mm_unpacklo_epi64(x1_tmp,x3_tmp);\
   x3 = _mm_unpackhi_epi64(x1_tmp,x3_tmp);

static short tw16[88]__attribute__((aligned(16))) = {
  32767,0,30272,-12539,23169,-23169,12539,-30272,
  32767,0,23169,-23169,0,-32767,-23169,-23169,
  32767,0,12539,-30272,-23169,-23169,-30272,12539};

// Fill the rest in here yourself
// You can use the stimulus from the first lab session if you like, but the vector

void dft16(int *x,int *y) {

}

main() {

}
