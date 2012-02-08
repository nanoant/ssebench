#include <xmmintrin.h>
#include <stdio.h>
#include <sys/time.h>
#include <string.h>

// http://fhtr.blogspot.com/2010/02/4x4-float-matrix-multiplication-using.html
#if defined(__cplusplus)
#include <algebra.hpp>
__inline void mmul(const float *a, const float *b, float *r)
{
	for(int i = 0; i < 16; i += 4) {
		vec4 rl = vec4(a) * vec4(b[i]);
		for(int j = 1; j < 4; j++) {
			rl += vec4(&a[j*4]) * vec4(b[i+j]);
		}
		rl >> &r[i];
	}
}
#elif TYPE == 0
__inline void mmul(const float *a, const float *b, float *r)
{
	int i, j;
	for (i=0; i<16; i+=4)
		for (j=0; j<4; j++)
			r[i+j] = b[i]*a[j] + b[i+1]*a[j+4] + b[i+2]*a[j+8] + b[i+3]*a[j+12];
}
#elif TYPE == 1
__inline void mmul(const float *a, const float *b, float *r)
{
	__m128 a_line, b_line, r_line;
	float mc[16] __attribute__((aligned(16)));  // 16-byte aligned temp array
	int i, j;
	for (i=0; i<16; i+=4) {
		b_line = _mm_load_ps(&b[i]);              // b_line = vec4(column(b, i))
                                              // remember that we are column-major
		for (j=0; j<4; j++) {
			a_line = _mm_set_ps(a[j+12], a[j+8], a[j+4], a[j]); 
                                              // a_line = vec4(row(a, j))
                                              // note that SSE is little-endian
			r_line = _mm_mul_ps(a_line, b_line);    // r_line = a_line * b_line
			_mm_store_ps(mc, r_line);               // copy r_line to memory
			r[i+j] = mc[0]+mc[1]+mc[2]+mc[3];       // r[i][j] = sum(r_line)
                                              //         = dot(a_line, b_line)
		}
	}
}
#elif TYPE == 2
__inline void mmul(const float * a, const float * b, float * r)
{
	__m128 a_line, b_line, r_line;
	int i, j;
	for (i=0; i<16; i+=4) {
    // unroll the first step of the loop to avoid having to initialize r_line to zero
		a_line = _mm_load_ps(a);         // a_line = vec4(column(a, 0))
		b_line = _mm_set1_ps(b[i]);      // b_line = vec4(b[i][0])
		r_line = _mm_mul_ps(a_line, b_line); // r_line = a_line * b_line
		for (j=1; j<4; j++) {
			a_line = _mm_load_ps(&a[j*4]); // a_line = vec4(column(a, j))
			b_line = _mm_set1_ps(b[i+j]);  // b_line = vec4(b[i][j])
                                     // r_line += a_line * b_line
			r_line = _mm_add_ps(_mm_mul_ps(a_line, b_line), r_line);
		}
		_mm_store_ps(&r[i], r_line);     // r[i] = r_line
	}
}
#endif

int main(int argc, char const *argv[])
{
	float a[16] = {1, 0, 0, 0,
	               0, 1, 0, 0,
	               0, 0, 1, 0,
	               0, 0, 0, 1};
	float b[16] = {0, 0,-1, 0,
	               0, 1, 0, 0,
	               1, 0, 0, 0,
	               0, 0, 0, 1};
	long i;
	long iter = 20000000;
	struct timeval start, end;
	gettimeofday(&start, NULL);
	if(argc >= 2) iter = atol(argv[1]);
	for(i = 0; i < iter; i++) {
		float r[16];
		mmul(a, b, r);
		memcpy(a, r, sizeof(r));
	}
	gettimeofday(&end, NULL);
	printf("done in %.2f seconds\n", end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec) / 1000000.0);
	printf("%.12g\t%.12g\t%.12g\t%.12g\n", a[0],  a[1],  a[2],  a[3]);
	printf("%.12g\t%.12g\t%.12g\t%.12g\n", a[4],  a[5],  a[6],  a[7]);
	printf("%.12g\t%.12g\t%.12g\t%.12g\n", a[8],  a[9],  a[10], a[11]);
	printf("%.12g\t%.12g\t%.12g\t%.12g\n", a[12], a[13], a[14], a[15]);
	return 0;
}
