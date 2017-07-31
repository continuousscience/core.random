#include <math.h>

// ST(Random, Int)
int rand_int(sil_State *S) {
    size_t len; // Always assume len is wrong!
    Random *v = (Random *)sil_getST(S, &len);

    sil_pushinteger(S, genrand_int32(v));
    return 0;
}

// ST(Random, Float)
int rand_real(sil_State *S) {
    size_t len; // Always assume len is wrong!
    Random *v = (Random *)sil_getST(S, &len);

    sil_pushdouble(S, genrand_real1(v));
    return 0;
}

// ST(Random, Float)
int rand_half(sil_State *S) {
    size_t len; // Always assume len is wrong!
    Random *v = (Random *)sil_getST(S, &len);

    sil_pushdouble(S, genrand_real2(v));
    return 0;
}

// ST(Random, Float)
int rand_open(sil_State *S) {
    size_t len; // Always assume len is wrong!
    Random *v = (Random *)sil_getST(S, &len);

    sil_pushdouble(S, genrand_real3(v));
    return 0;
}

int normal(sil_State *S) {
    size_t len; // Always assume len is wrong!
    Random *v = (Random *)sil_getST(S, &len);
    double y1, y2;
    normaldist(v, &y1, &y2);
    sil_pushdouble(S, y1);
    sil_pushdouble(S, y2);
    sil_settuple(S, 2);
    return 0;
}

// Int -> Random
int mkRandom(sil_State *S) {
    int seed = sil_tointeger(S, 1);
    Random *v = (Random *)malloc(sizeof(Random));
    init_genrand(v, seed);
    sil_settop(S, 0);
    sil_newuserdata(S, random_hash, v);
    return 0;
}

