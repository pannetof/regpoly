/*
 * regpoly.c — Miscellaneous utility functions used throughout REGPOLY.
 */

#include <stdio.h>
#include <stdlib.h>
#include "regpoly.h"
#include "vectorsF2.h"

/* Computes the inverse of x in the ring Z/modulo*Z using the extended
   Euclidean algorithm. Exits with an error if no inverse exists. */
int InverseModN(int x, int modulo)
{
  int remainder, quotient, mod;
  int p2 = 0, p1 = 1, p, q2 = 0, q1 = 0;

  mod = modulo;
  do {
    quotient = mod / x;
    remainder = mod - quotient * x;
    p = (2 * modulo + p2 - p1 * q2) % modulo;
    if (p < 0)
      p += modulo;
    q2 = q1;  q1 = quotient;  p2 = p1;  p1 = p;
    mod = x;  x = remainder;
  } while (remainder != 0);

  if (mod != 1) {
    printf("No inverse exists!\n");
    exit(1);
  }
  p = (2 * modulo + p2 - p1 * q2) % modulo;
  if (p < 0)
    p += modulo;
  return p;
}

/* Returns the smaller of A and B. */
int intmin(int A, int B)
{
  if (A > B)
    return B;
  else
    return A;
}

/* Returns the larger of A and B. */
int intmax(int A, int B)
{
  if (A > B)
    return A;
  else
    return B;
}

/* Returns the smaller of A and B (float version). */
float floatmin(float A, float B)
{
  if (A > B)
    return B;
  else
    return A;
}

/* Returns the larger of A and B (float version). */
float floatmax(float A, float B)
{
  if (A > B)
    return A;
  else
    return B;
}

/* Returns gcd(x, y) using the Euclidean algorithm. */
int gcd(int x, int y)
{
  int x1 = x, y1 = y;

  if (x == 0)
    return y;
  while (x1 != 0) {
    y1 = y1 % x1;
    if (y1 != 0)
      x1 = x1 % y1;
    else
      return x1;
  }
  return y1;
}

/* Reads and discards all characters in file f up to and including
   the next newline or end-of-file. */
void ReadLn(FILE *f)
{
  int i;
  while ((i = fgetc(f)) != EOF && i != '\n')
    ;
}

/* Fills the bit vector v with random bits using MRG32k3a. */
void RandVect(BitVect *v)
{
  int j;
  for (j = 0; j < v->n; j++)
    v->vect[j] = MRG32k3a();
}

/* --- MRG32k3a: combined multiple-recursive generator of period ~2^192 --- */

static double S10, S11, S12, S20, S21, S22;

#define norm1   2.32830654983782883e-10
#define norm2   2.328318825240739e-10
#define m1      4294967087.0
#define m2      4294944443.0
#define a12     1403580.0
#define a13n    810728.0
#define a21     527612.0
#define a23n    1370589.0

/* Returns the next 32-bit output of the MRG32k3a generator. */
uint32_t MRG32k3a(void)
{
  long k;
  double p;

  /* Component 1 */
  p = a12 * S11 - a13n * S10;
  k = p / m1;
  p -= k * m1;
  if (p < 0.0)
    p += m1;
  S10 = S11;  S11 = S12;  S12 = p;

  /* Component 2 */
  p = a21 * S22 - a23n * S20;
  k = p / m2;
  p -= k * m2;
  if (p < 0.0)
    p += m2;
  S20 = S21;  S21 = S22;  S22 = p;

  /* Combination */
  if (S12 <= S22)
    return (uint32_t) (S12 - S22 + m1);
  else
    return (uint32_t) (S12 - S22);
}

/* Returns the next output of MRG32k3a as a double in [0, 1). */
double doubleMRG32k3a(void)
{
  return MRG32k3a() * 2.328306436539e-10;
}

/* Sets the initial state (seed) of the MRG32k3a generator. */
void SetMRG32k3a(double x10, double x11, double x12,
                 double x20, double x21, double x22)
{
  S10 = x10;  S11 = x11;  S12 = x12;
  S20 = x20;  S21 = x21;  S22 = x22;
}

/* Returns TRUE if 2^exp - 1 is a known Mersenne prime (<= 2^110503 - 1). */
int ValidMersennePrime(int exp)
{
  if (exp == 2   || exp == 3   || exp == 5    || exp == 7    ||
      exp == 13  || exp == 17  || exp == 19   || exp == 31   ||
      exp == 61  || exp == 89  || exp == 107  || exp == 127  ||
      exp == 521 || exp == 607 || exp == 1279 || exp == 2203 ||
      exp == 2281 || exp == 3217 || exp == 4253 || exp == 4423 ||
      exp == 9689 || exp == 9941 || exp == 11213 || exp == 19937 ||
      exp == 21701 || exp == 23209 || exp == 44497 || exp == 86243 ||
      exp == 110503)
    return 1;
  else
    return 0;
}
