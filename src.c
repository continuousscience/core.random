#include <math.h>

/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/er->mt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

/* initializes r->mt[N] with a seed */
void init_genrand(Random *r, unsigned long s) {
    r->mt[0]= s & 0xffffffffUL;
    for (r->mti=1; r->mti<N; r->mti++) {
        r->mt[r->mti] = 
	    (1812433253UL * (r->mt[r->mti-1] ^ (r->mt[r->mti-1] >> 30)) + r->mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array r->mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        r->mt[r->mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(Random *r, unsigned long init_key[], int key_length) {
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        r->mt[i] = (r->mt[i] ^ ((r->mt[i-1] ^ (r->mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        r->mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { r->mt[0] = r->mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        r->mt[i] = (r->mt[i] ^ ((r->mt[i-1] ^ (r->mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        r->mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { r->mt[0] = r->mt[N-1]; i=1; }
    }

    r->mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(Random *r) {
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (r->mti >= N) { /* generate N words at one time */
        int kk;

        //if (r->mti == N+1)   /* if init_genrand() has not been called, */
        //    init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (r->mt[kk]&UPPER_MASK)|(r->mt[kk+1]&LOWER_MASK);
            r->mt[kk] = r->mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (r->mt[kk]&UPPER_MASK)|(r->mt[kk+1]&LOWER_MASK);
            r->mt[kk] = r->mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (r->mt[N-1]&UPPER_MASK)|(r->mt[0]&LOWER_MASK);
        r->mt[N-1] = r->mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        r->mti = 0;
    }
  
    y = r->mt[r->mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(Random *r) {
    return genrand_int32(r)*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(Random *r) {
    return genrand_int32(r)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(Random *r) {
    return (((double)genrand_int32(r)) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(Random *r) { 
    unsigned long a=genrand_int32(r)>>5, b=genrand_int32(r)>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */

// ST(Random, Int)
int randi(sil_State *S) {
    size_t len; // Always assume len is wrong!
    Random *v = (Random *)sil_getST(S, &len);

    sil_pushinteger(S, genrand_int32(v));
    return 0;
}

// ST(Random, Float)
int rand(sil_State *S) {
    size_t len; // Always assume len is wrong!
    Random *v = (Random *)sil_getST(S, &len);

    sil_pushdouble(S, genrand_real1(v));
    return 0;
}

// ST(Random, Float)
int randh(sil_State *S) {
    size_t len; // Always assume len is wrong!
    Random *v = (Random *)sil_getST(S, &len);

    sil_pushdouble(S, genrand_real2(v));
    return 0;
}

// ST(Random, Float)
int randc(sil_State *S) {
    size_t len; // Always assume len is wrong!
    Random *v = (Random *)sil_getST(S, &len);

    sil_pushdouble(S, genrand_real3(v));
    return 0;
}

int normal(sil_State *S) {
    size_t len; // Always assume len is wrong!
    Random *v = (Random *)sil_getST(S, &len);

    double x1, x2, w;
    do {
        x1 = 2.*genrand_real3(v) - 1.0;
        x2 = 2.*genrand_real3(v) - 1.0;
        w = x1*x1 + x2*x2;
    } while ( w >= 1.0 );

    w = sqrt( (-2.*log(w)) / w );

    sil_pushdouble(S, x1*w);
    sil_pushdouble(S, x2*w);
    sil_settuple(S, 2);
    return 0;
}

// Int -> Random
int mkRandom(sil_State *S) {
    int seed = sil_tointeger(S, 1);
    Random *v = (Random *)malloc(sizeof(Random));
    init_genrand(v, seed);
    sil_pushRandom(S, v);
    return 0;
}

