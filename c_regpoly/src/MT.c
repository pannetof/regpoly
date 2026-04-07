/*
 * MT.c — Mersenne Twister (MT) generator for REGPOLY.
 *
 * Implements the Mersenne Twister PRNG parameterized by word size w,
 * register count r, middle distance m, offset p, and twist matrix
 * coefficient A. The characteristic polynomial is computed using the ZEN
 * library for GF(2) polynomial arithmetic.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <limits.h>
#include <zen.h>
#include "MT.h"
#include "regpoly.h"
#include "polynomials.h"

#define MERSNAME "Mersenne Twister"

#define TW    ((paramMT *)(Gen->ParamGen))->w
#define TR    ((paramMT *)(Gen->ParamGen))->r
#define TM    ((paramMT *)(Gen->ParamGen))->m
#define TP    ((paramMT *)(Gen->ParamGen))->p
#define TA    ((paramMT *)(Gen->ParamGen))->a
#define TI    ((paramMT *)(Gen->ParamGen))->i
#define TO    (Gen->L)
#define TU    ((paramMT *)(Gen->ParamGen))->uu
#define TL    ((paramMT *)(Gen->ParamGen))->ll
#define MASKW ((paramMT *)(Gen->ParamGen))->mw
#define ETAT  (Gen->GenState)

static int k, s, j, t, ii, size;
static uint32_t Y, V_i;

/* Returns the name string for the Mersenne Twister generator type. */
char *MTName(void)
{
  return MERSNAME;
}

/* Prints the parameters (w, r, m, p, a) of the Mersenne Twister. */
void DispMT(Generateur *Gen)
{
  printf(" w= %3u  r= %3u  m= %3u  p= %3u a= %10x  \n",
         ((paramMT *)(Gen->ParamGen))->w,
         ((paramMT *)(Gen->ParamGen))->r,
         ((paramMT *)(Gen->ParamGen))->m,
         ((paramMT *)(Gen->ParamGen))->p,
         ((paramMT *)(Gen->ParamGen))->a);
}

/* Initializes all parameters of a Mersenne Twister generator.
   Exits if w > 32. Wires all function pointers for this type. */
void InitParamMT(Generateur *Gen, int w, int r, int m, int p, uint32_t A, int L)
{
  if (w > 32) {
    printf("Error in InitParamMT(): w must be 2, 4, 8, 16 or 32");
    exit(1);
  }
  ((paramMT *)(Gen->ParamGen))->w  = w;
  ((paramMT *)(Gen->ParamGen))->r  = r;
  ((paramMT *)(Gen->ParamGen))->m  = m;
  ((paramMT *)(Gen->ParamGen))->p  = p;
  ((paramMT *)(Gen->ParamGen))->a  = A;
  ((paramMT *)(Gen->ParamGen))->i  = 0;
  ((paramMT *)(Gen->ParamGen))->ll = 0;
  ((paramMT *)(Gen->ParamGen))->uu = 0;
  ((paramMT *)(Gen->ParamGen))->mw = 0;

  Gen->k       = w * r - p;
  Gen->L       = L;
  Gen->Step    = 1;
  Gen->smax    = INT_MAX;
  Gen->DispGen   = DispMT;
  Gen->DispName  = MTName;
  Gen->InitGen   = InitMT;
  Gen->Iteration = MT;
  Gen->CopyGen   = CopyMT;
  Gen->AllocGen  = AllocMT;
  Gen->PolyChar  = CharMT;
  PutBVToZero(&(Gen->GenState));
}

/* Copies all parameters and state from Mersenne Twister Gen2 into Gen1. */
void CopyMT(Generateur *Gen1, Generateur *Gen2)
{
  ((paramMT *)(Gen1->ParamGen))->w  = ((paramMT *)(Gen2->ParamGen))->w;
  ((paramMT *)(Gen1->ParamGen))->r  = ((paramMT *)(Gen2->ParamGen))->r;
  ((paramMT *)(Gen1->ParamGen))->m  = ((paramMT *)(Gen2->ParamGen))->m;
  ((paramMT *)(Gen1->ParamGen))->p  = ((paramMT *)(Gen2->ParamGen))->p;
  ((paramMT *)(Gen1->ParamGen))->a  = ((paramMT *)(Gen2->ParamGen))->a;
  ((paramMT *)(Gen1->ParamGen))->i  = ((paramMT *)(Gen2->ParamGen))->i;
  ((paramMT *)(Gen1->ParamGen))->ll = ((paramMT *)(Gen2->ParamGen))->ll;
  ((paramMT *)(Gen1->ParamGen))->uu = ((paramMT *)(Gen2->ParamGen))->uu;
  ((paramMT *)(Gen1->ParamGen))->mw = ((paramMT *)(Gen2->ParamGen))->mw;

  Gen1->k         = Gen2->k;
  Gen1->L         = Gen2->L;
  Gen1->smax      = Gen2->smax;
  Gen1->Step      = Gen2->Step;
  Gen1->DispGen   = Gen2->DispGen;
  Gen1->DispName  = Gen2->DispName;
  Gen1->InitGen   = Gen2->InitGen;
  Gen1->Iteration = Gen2->Iteration;
  Gen1->PolyChar  = Gen2->PolyChar;
  Gen1->CopyGen   = Gen2->CopyGen;
  Gen1->AllocGen  = Gen2->AllocGen;
  CopyBV(&(Gen1->GenState), &(Gen2->GenState));
}

/* Allocates and returns a new Mersenne Twister generator with state of k bits. */
Generateur *AllocMT(int k)
{
  Generateur *G;
  G = (Generateur *) malloc(sizeof(Generateur));
  AllocGen(G, k);
  G->ParamGen = (paramMT *) malloc(sizeof(paramMT));
  return G;
}

/* Frees all memory associated with the Mersenne Twister generator G. */
void FreeMT(Generateur *G)
{
  free(G->ParamGen);
  FreeGen(G);
}

/* Initializes the Mersenne Twister state from the bit vector init,
   sets up the upper/lower masks TU and TL, and writes L bits to retour. */
void InitMT(Generateur *Gen, BitVect *init, BitVect *retour)
{
  CopyBV(&ETAT, init);

  TI    = 0;
  MASKW = TU = 0xffffffffU >> (WL - TW);
  TU <<= TP;

  if (TP != 0)
    TL = 0xffffffffU >> (WL - TP);
  else
    TL = 0x0U;

  t = 0;
  k = 0;
  if (TW == WL) {
    for (j = 0; j <= (TO - 1) / WL; j++)
      retour->vect[j] = ETAT.vect[j];
  } else {
    s = t - WL * k;
    for (j = 0; j <= (TO - 1) / WL; j++) {
      if (s != 0)
        retour->vect[j] = (ETAT.vect[j] << s) ^ (ETAT.vect[(j + 1) % (TR * TW / WL)] >> (WL - s));
      else
        retour->vect[j] = ETAT.vect[j];
    }
  }
}

/* Returns the i-th w-bit word from the generator state. */
static uint32_t V(Generateur *Gen, int i)
{
  ii = TR - i - 1;
  if (TW == WL)
    return ETAT.vect[ii];
  t = ii * TW;
  k = t / WL;
  s = t - WL * k;
  return (ETAT.vect[k] >> (WL - TW - s)) & MASKW;
}

/* Sets the i-th w-bit word of the generator state to val. */
static void SetV(Generateur *Gen, int i, uint32_t val)
{
  ii = TR - i - 1;
  if (TW == WL) {
    ETAT.vect[ii] = val;
    return;
  }
  t = ii * TW;
  k = t / WL;
  s = t - WL * k;
  val &= MASKW;
  val <<= (WL - TW - s);
  ETAT.vect[k] &= ~(MASKW << (WL - TW - s));
  ETAT.vect[k] |= val;
}

/* Performs one Mersenne Twister recurrence step and writes L output bits
   to retour. */
void MT(Generateur *Gen, BitVect *retour)
{
  Y = (V(Gen, TI % TR) & TU) | (V(Gen, (TI + 1) % TR) & TL);

  if (Y & 1U)
    V_i = V(Gen, (TI + TM) % TR) ^ (Y >> 1) ^ TA;
  else
    V_i = V(Gen, (TI + TM) % TR) ^ (Y >> 1);

  SetV(Gen, TI, V_i);

  ii   = TR - TI - 1;
  t    = ii * TW;
  k    = t / WL;
  size = TW * TR / WL;
  if (TW == WL) {
    for (j = 0; j <= (TO - 1) / WL; j++)
      retour->vect[j] = ETAT.vect[(k + j) % size];
  } else {
    s = t - WL * k;
    for (j = 0; j <= (TO - 1) / WL; j++) {
      if (s != 0)
        retour->vect[j] = (ETAT.vect[(k + j) % size] << s)
                        ^ (ETAT.vect[(k + j + 1) % size] >> (WL - s));
      else
        retour->vect[j] = ETAT.vect[(k + j) % size];
    }
  }
  TI = (TI + 1) % TR;
}

/* Reads Mersenne Twister parameters from filename and adds each generator
   to component E. Exits if the file cannot be opened, w > 32, or allocation
   fails. */
void ReadDataMT(Component *E, char *filename, boolean same, int L)
{
  int nbgen;
  int w, r, p, m;
  uint32_t a;
  FILE *f;
  Generateur *g;

  f = fopen(filename, "r");
  if (f == NULL) {
    printf("File %s not found\n", filename);
    exit(1);
  }
  ReadLn(f);
  fscanf(f, "%d", &nbgen);
  ReadLn(f);
  AllocGensInComponent(E, nbgen, same);

  for (j = 0; j < nbgen; j++) {
    fscanf(f, "%d %d %d %d %x", &w, &r, &m, &p, &a);
    if (w > 32) {
      printf("w must be <= 32 bits\n");
      exit(1);
    }
    if ((g = AllocMT(w * r - p)) == NULL) {
      printf("Error in ReadDataMT()\n");
      exit(1);
    }
    ReadLn(f);
    InitParamMT(g, w, r, m, p, a, L);
    AddGenInComponent(E, g);
    FreeMT(g);
  }
  fclose(f);
}

/* Generates nb random Mersenne Twister parameters with the given Mersenne prime
   exponent exp, word size w, and middle distance m, and writes them to filename.
   Exits if exp is not a valid Mersenne prime exponent. */
void ProduceParamsMT(char *filename, int nb, int exp, int w, int m)
{
  int r, p, i, k, primitif;
  uint32_t LSB, y, maskw, u, ll, a, *init, *x;
  int count = 0;
  int nn = 0;
  FILE *f;

  if (!ValidMersennePrime(exp)) {
    printf("The degree of the Mersenne Twister is not a Mersenne prime exponent\n");
    exit(1);
  }
  SetMRG32k3a((double)time(NULL), (double)time(NULL), (double)time(NULL),
              (double)time(NULL), (double)time(NULL), (double)time(NULL));
  r    = exp / w + 1;
  p    = w * r - exp;
  init = (uint32_t *) malloc(r * sizeof(uint32_t));
  x    = (uint32_t *) malloc(2 * exp * sizeof(uint32_t));
  f    = fopen(filename, "w");

  fprintf(f, "MT\n%d\n", nb);

  maskw = 0xffffffffU >> (32 - w);
  u     = 0xffffffffU << p;
  ll    = ~u;
  LSB   = 0x00000001U << (32 - w);
  for (j = 0; j < r; j++)
    init[j] = j + 1;

  while (count < nb) {
    a  = MRG32k3a() & maskw;
    a |= 1U << (w - 1);

    for (j = 0; j < r; j++)
      x[j] = init[j];
    if ((x[2] & LSB) == (x[3] & LSB))
      exit(1);

    for (i = 0; i < exp; i++) {
      for (j = r; j < 2 * exp; j++) {
        y = (x[j - r] & u) | (x[j - r + 1] & ll);
        if (y & LSB)
          x[j] = x[j + m - r] ^ (y >> 1) ^ a;
        else
          x[j] = x[j + m - r] ^ (y >> 1);
      }
      for (j = 2; j <= exp; j++)
        x[j] = x[2 * j - 1];
      for (k = exp; k >= r; k--) {
        if (x[k - r + 1] & LSB)
          y = x[k] ^ x[k - r + m] ^ a;
        else
          y = x[k] ^ x[k - r + m];
        y = (y << 1) & maskw;
        y |= x[k - r + 1] & LSB;
        x[k - r + 1] = (u & x[k - r + 1]) | (ll & y);
        x[k - r]     = (u & y) | (ll & x[k - r]);
      }
    }

    primitif = 1;
    for (j = 1; j < r; j++) {
      if (init[j] != x[j]) {
        primitif = 0;
        break;
      }
    }
    if (((x[0] & u) == (init[0] & u)) && primitif) {
      count++;
      fprintf(f, "%d %d %d %d %08x\n", w, r, m, p, a);
    } else {
      nn++;
    }
  }
  printf("nn=%d\n", nn);
  fclose(f);
}

/* Computes the characteristic polynomial of the Mersenne Twister using
   the ZEN library. Stores the result in coeff[] and BVPoly. */
void CharMT(Generateur *Gen, int coeff[], BitVect *BVPoly)
{
  int j;
  ZENPoly tntm, tn1tm1, res, temp, temp2, tntmwr;
  ZENRing F2;
  BigNum two;
  BigNumLength twol;
  ZENElt ONE, Coeff;

  ZBNReadFromString(&two, &twol, "2", 10);
  if (ZENBaseRingAlloc(F2, two, twol)) {
    printf("Error in CharMT(): ZENBaseRingAlloc failed\n");
    exit(1);
  }
  ZENEltAlloc(ONE, F2);
  ZENEltSetToOne(ONE, F2);
  ZENEltAlloc(Coeff, F2);

  ZENPolyAlloc(temp,   TR * TW, F2);
  ZENPolyAlloc(temp2,  TR * TW, F2);
  ZENPolyAlloc(res,    TR * TW, F2);

  /* t^n + t^m */
  ZENPolyAlloc(tntm, TR * TW, F2);
  ZENPolySetToXi(tntm, TR, F2);
  ZENPolySetCoeff(tntm, TM, ONE, F2);

  /* t^(n-1) + t^(m-1) */
  ZENPolyAlloc(tn1tm1, TR * TW, F2);
  ZENPolySetToXi(tn1tm1, TR - 1, F2);
  ZENPolySetCoeff(tn1tm1, TM - 1, ONE, F2);

  /* (t^n + t^m)^(w-p) */
  ZENPolyAlloc(tntmwr, TR * TW, F2);
  PowPoly(&tntmwr, &tntm, &F2, TW - TP);

  ZENPolySetToZero(res, F2);
  for (j = 0; j < TP; j++)
    if (TA & (1U << j)) {
      PowPoly(&temp, &tn1tm1, &F2, TP - j - 1);
      ZENPolyMultiply(temp2, temp, tntmwr, F2);
      ZENPolyAdd(res, temp2, F2);
    }
  for (j = TP; j < TW; j++)
    if (TA & (1U << j)) {
      PowPoly(&temp, &tntm, &F2, TW - j - 1);
      ZENPolyAdd(res, temp, F2);
    }
  PowPoly(&temp, &tn1tm1, &F2, TP);
  ZENPolyMultiply(temp2, temp, tntmwr, F2);
  ZENPolyAdd(res, temp2, F2);

  for (j = 0; j <= TR * TW - TP; j++) {
    ZENPolyGetCoeff(Coeff, res, j, F2);
    coeff[j] = ZENEltAreEqual(Coeff, ONE, F2) ? 1 : 0;
  }
  ZENPolyFree(tntm,   F2);
  ZENPolyFree(tn1tm1, F2);
  ZENPolyFree(res,    F2);
  ZENPolyFree(temp,   F2);
  ZENPolyFree(temp2,  F2);
  ZENEltFree(ONE, F2);
  ZBNF(two);
  if (ZENRingClose(F2)) {
    printf("Error in CharMT(): ZENRingClose failed\n");
    exit(1);
  }
  PutBVToZero(BVPoly);
  for (j = 0; j <= Gen->k; j++)
    PutBitBV(BVPoly, j, coeff[j]);
}
