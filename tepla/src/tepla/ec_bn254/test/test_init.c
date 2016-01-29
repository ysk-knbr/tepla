#include <assert.h>

#include "rdtsc.h"

#include "../ec_bn254_lcl.h"

#define N 1000

//============================================
//   initialize test
//============================================
void test_init_field()
{
    
    fprintf(stderr, "filed init test\n");
    int i;
    unsigned long long int t1, t2;

    Field fa, fb;

    t1 = clock();
    for (i = 0; i < N; i++) {
        field_init(fa, "bn254_fpa");
    }
    t2 = clock();
    printf("field init fpa: %.5lf [msec]\n", (double)(t2 - t1) / N / CLOCKS_PER_SEC * 1000);
    field_clear(fa);

    t1 = clock();
    for (i = 0; i < N; i++) {
        field_init(fa, "bn254_fp2a");
    }
    t2 = clock();
    printf("field init fp2a: %.5lf [msec]\n", (double)(t2 - t1) / N / CLOCKS_PER_SEC * 1000);
    field_clear(fa);
    
    t1 = clock();
    for (i = 0; i < N; i++) {
        field_init(fa, "bn254_fp6a");
    }
    t2 = clock();
    printf("field init fp6a: %.5lf [msec]\n", (double)(t2 - t1) / N / CLOCKS_PER_SEC * 1000);
    field_clear(fa);
    
    t1 = clock();
    for (i = 0; i < N; i++) {
        field_init(fa, "bn254_fp12a");
    }
    t2 = clock();
    printf("field init fp12a: %.5lf [msec]\n", (double)(t2 - t1) / N / CLOCKS_PER_SEC * 1000);
    field_clear(fa);

    t1 = clock();
    for (i = 0; i < N; i++) {
        field_init(fb, "bn254_fpb");
    }
    t2 = clock();
    printf("field init fpb: %.5lf [msec]\n", (double)(t2 - t1) / N / CLOCKS_PER_SEC * 1000);
    field_clear(fb);
    
    t1 = clock();
    for (i = 0; i < N; i++) {
        field_init(fb, "bn254_fp2b");
    }
    t2 = clock();
    printf("field init fp2b: %.5lf [msec]\n", (double)(t2 - t1) / N / CLOCKS_PER_SEC * 1000);
    field_clear(fb);
    
    t1 = clock();
    for (i = 0; i < N; i++) {
        field_init(fb, "bn254_fp6b");
    }
    t2 = clock();
    printf("field init fp6b: %.5lf [msec]\n", (double)(t2 - t1) / N / CLOCKS_PER_SEC * 1000);
    field_clear(fb);
    
    t1 = clock();
    for (i = 0; i < N; i++) {
        field_init(fb, "bn254_fp12b");
    }
    t2 = clock();
    printf("field init fp12b: %.5lf [msec]\n", (double)(t2 - t1) / N / CLOCKS_PER_SEC * 1000);
    field_clear(fb);
}

void test_init_element()
{
    fprintf(stderr, "element init test\n");
    
    int i;
    unsigned long long int t1, t2;

    Field fa, fb;
    Element a;
    
    field_init(fa, "bn254_fpa");
    t1 = clock();
    for (i = 0; i < N; i++) {
        element_init(a, fa);
    }
    t2 = clock();
    printf("element init fpa: %.5lf [msec]\n", (double)(t2 - t1) / N / CLOCKS_PER_SEC * 1000);
    element_clear(a);
    field_clear(fa);

    field_init(fa, "bn254_fp2a");
    t1 = clock();
    for (i = 0; i < N; i++) {
        element_init(a, fa);
    }
    t2 = clock();
    printf("element init fp2a: %.5lf [msec]\n", (double)(t2 - t1) / N / CLOCKS_PER_SEC * 1000);
    element_clear(a);
    field_clear(fa);
    
    field_init(fa, "bn254_fp6a");
    t1 = clock();
    for (i = 0; i < N; i++) {
        element_init(a, fa);
    }
    t2 = clock();
    printf("element init fp6a: %.5lf [msec]\n", (double)(t2 - t1) / N / CLOCKS_PER_SEC * 1000);
    element_clear(a);
    field_clear(fa);

    field_init(fa, "bn254_fp12a");
    t1 = clock();
    for (i = 0; i < N; i++) {
        element_init(a, fa);
    }
    t2 = clock();
    printf("element init fp12a: %.5lf [msec]\n", (double)(t2 - t1) / N / CLOCKS_PER_SEC * 1000);
    element_clear(a);
    field_clear(fa);

    field_init(fb, "bn254_fpb");
    t1 = clock();
    for (i = 0; i < N; i++) {
        element_init(a, fb);
    }
    t2 = clock();
    printf("element init fpb: %.5lf [msec]\n", (double)(t2 - t1) / N / CLOCKS_PER_SEC * 1000);
    element_clear(a);
    field_clear(fb);
    
    field_init(fb, "bn254_fp2b");
    t1 = clock();
    for (i = 0; i < N; i++) {
        element_init(a, fb);
    }
    t2 = clock();
    printf("element init fp2b: %.5lf [msec]\n", (double)(t2 - t1) / N / CLOCKS_PER_SEC * 1000);
    element_clear(a);
    field_clear(fb);
    
    field_init(fb, "bn254_fp6b");
    t1 = clock();
    for (i = 0; i < N; i++) {
        element_init(a, fb);
    }
    t2 = clock();
    printf("element init fp6b: %.5lf [msec]\n", (double)(t2 - t1) / N / CLOCKS_PER_SEC * 1000);
    element_clear(a);
    field_clear(fb);

    field_init(fb, "bn254_fp12b");
    t1 = clock();
    for (i = 0; i < N; i++) {
        element_init(a, fb);
    }
    t2 = clock();
    printf("element init fp12b: %.5lf [msec]\n", (double)(t2 - t1) / N / CLOCKS_PER_SEC * 1000);
    element_clear(a);
    field_clear(fb);
}

void test_init_curve()
{
    fprintf(stderr, "curve init test\n");
    
    int i;
    unsigned long long int t1, t2;

    EC_GROUP fp, tw;
    
    t1 = clock();
    for (i = 0; i < N; i++) {
        curve_init(fp, "ec_bn254_fpa");
    }
    t2 = clock();
    printf("curve init fpa: %.5lf [msec]\n", (double)(t2 - t1) / N / CLOCKS_PER_SEC * 1000);
    curve_clear(fp);

    t1 = clock();
    for (i = 0; i < N; i++) {
        curve_init(tw, "ec_bn254_twa");
    }
    t2 = clock();
    printf("curve init twa: %.5lf [msec]\n", (double)(t2 - t1) / N / CLOCKS_PER_SEC * 1000);
    curve_clear(tw);
    
    t1 = clock();
    for (i = 0; i < N; i++) {
        curve_init(fp, "ec_bn254_fpb");
    }
    t2 = clock();
    printf("curve init fpb: %.5lf [msec]\n", (double)(t2 - t1) / N / CLOCKS_PER_SEC * 1000);
    curve_clear(fp);

    t1 = clock();
    for (i = 0; i < N; i++) {
        curve_init(tw, "ec_bn254_twb");
    }
    t2 = clock();
    printf("curve init twb: %.5lf [msec]\n", (double)(t2 - t1) / N / CLOCKS_PER_SEC * 1000);
    curve_clear(tw);
}

void test_init_point()
{
    fprintf(stderr, "point init test\n");
    
    int i;
    unsigned long long int t1, t2;

    EC_GROUP fp, tw;
    EC_POINT P, Q;

    curve_init(fp, "ec_bn254_fpa");
    t1 = clock();
    for (i = 0; i < N; i++) {
        point_init(P, fp);
    }
    t2 = clock();
    printf("point init fpa: %.5lf [msec]\n", (double)(t2 - t1) / N / CLOCKS_PER_SEC * 1000);
    point_clear(P);
    curve_clear(fp);

    curve_init(tw, "ec_bn254_twa");
    t1 = clock();
    for (i = 0; i < N; i++) {
        point_init(Q, tw);
    }
    t2 = clock();
    printf("point init twa: %.5lf [msec]\n", (double)(t2 - t1) / N / CLOCKS_PER_SEC * 1000);
    point_clear(Q);
    curve_clear(tw);
    
    curve_init(fp, "ec_bn254_fpb");
    t1 = clock();
    for (i = 0; i < N; i++) {
        point_init(P, fp);
    }
    t2 = clock();
    printf("point init fpb: %.5lf [msec]\n", (double)(t2 - t1) / N / CLOCKS_PER_SEC * 1000);
    point_clear(P);
    curve_clear(fp);

    curve_init(tw, "ec_bn254_twb");
    t1 = clock();
    for (i = 0; i < N; i++) {
        point_init(Q, tw);
    }
    t2 = clock();
    printf("curve init twb: %.5lf [msec]\n", (double)(t2 - t1) / N / CLOCKS_PER_SEC * 1000);
    point_clear(Q);
    curve_clear(tw);
}

void test_init_pairing()
{
    fprintf(stderr, "pairing init test\n");
    
    int i;
    unsigned long long int t1, t2;

    EC_PAIRING p1, p2;
    
    t1 = clock();
    for (i = 0; i < N; i++) {
        pairing_init(p1, "ECBN254a");
    }
    t2 = clock();
    printf("pairing init ECBN254a: %.5lf [msec]\n", (double)(t2 - t1) / N / CLOCKS_PER_SEC * 1000);
    pairing_clear(p1);

    t1 = clock();
    for (i = 0; i < N; i++) {
        pairing_init(p2, "ECBN254b");
    }
    t2 = clock();
    printf("pairing init ECBN254b: %.5lf [msec]\n", (double)(t2 - t1) / N / CLOCKS_PER_SEC * 1000);
    pairing_clear(p2);
}
//============================================
// main program
//============================================
int main(void)
{
    test_init_field();
    
    test_init_element();

    test_init_curve();

    test_init_point();

    test_init_pairing();

    fprintf(stderr, "ok\n");

    return 0;
}
