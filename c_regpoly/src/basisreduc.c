#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include "basisreduc.h"
#include "timer.h"
#include "regpoly.h"
#include "polynomials.h"
#include "myzen.h"

#define vperm(i) (B->permutations[i])
#define vb(i) B->vect[i]
#define NormeVect(P) ((P)->deg)
#define DB(i) {printf("%d\n", i); fflush(stdout);}
#define TIMING 0

static void GetPolys(Combinaison *C, BitVect *M, BitVect BVPolys[], int resolution);
static void GetPolysSpecial(Combinaison *C, BitVect *M, BitVect BVPolys[], int resolution, int rotative_shift);
static void AllocPolVect(Base *B, PolVect *P);
static void FreePolVect(Base *B, PolVect *P);
static inline void CopyPolVect(Base *B, PolVect *a, PolVect *b);
static void VectorToZero(Base *B, PolVect *P);
static inline void VectorToZeroLazy(PolVect *P);
static int NormeElementVector(Base *B, PolVect *P, int j);
static void SetElementVectorToZero(PolVect *P, int j);
static void SetElementVectorToOne(PolVect *P, int j);
static void BVToElementVector(Base *B, BitVect *A, PolVect *P, int j);
static inline void PermuteCoord(Base *B, int m);
static inline void renumber(Base *B, int m, int resolution);
static inline void SolveAxb(Base *B, int m, int *x);
static inline void AddPolVect(Base *B, PolVect *res, PolVect *a, PolVect *b);
static inline void AddSelfPolVectMult(Base *B, PolVect *res, PolVect *a, int s);
static inline int ElementVectorNormeEqual(Base *B, PolVect *P, int j, int norme);
static void DispBase(Base *B, int resolution);
static void DispElementVector(Base *B, PolVect *P, int j);
static void DispVector(Base *B, PolVect *P);
static void PolyToElt(ZENElt *e, ZENRing polring, ZENPoly p, ZENRing F2);
static void EltToPoly(ZENPoly *p, ZENRing F2, ZENElt e, ZENRing polring);
static inline void CopyPolVectLazy(PolVect *a, PolVect *b);
static inline void Multbys(Base *B, PolVect *r, PolVect *a, int s);
static int NormeBitVect(Base *B, BitVect *A);

/* Computes the equidistribution of C using polynomial lattice basis reduction;
   fills params->ecart and related fields. Skips if meverif is FALSE. */
void TestMELat(Combinaison *C, paramMecf *params) {
  Base B;
  BitVect *gVect;
  BitVect P;
  timer_Chrono timer;
  int Deg, l, d, *coeff, lenght, i, RES, maxl;

  if (params->meverif == FALSE)
    return;

  if (params->L < C->L) {
    printf("Error in TestME().\n  Resolution of allocated paramMecf < resolution of generator\n");
    exit(1);
  }
  if (TIMING)
    timer_Init(&timer); /* Initialize the timer */
  Deg = C->k_g;
  RES = intmin(params->L, Deg);
  for (l = 1; l <= params->L; l++) /* Reset the ecart array */
    params->ecart[l] = -1;

  gVect = (BitVect *) malloc(RES * sizeof(BitVect));
  for (i = 0; i < RES; i++)
    AllocBV(&(gVect[i]), Deg + 1);
  coeff = (int *) malloc((Deg + 1) * sizeof(int));
  AllocBV(&P, Deg + 1);
  polycharComb(C, coeff, &P);
  free(coeff);
  GetPolys(C, &P, gVect, RES);

  AllocBase(&B, RES, Deg);
  B.degmax = Deg;
  DualBase(&B, gVect, &P, 1);
  SetPsi12(C, params);
  params->se = 0;
  maxl = params->L;
  if (TIMING)
    printf("time=%6.3f\n", timer_Val(timer, timer_sec)); fflush(stdout);
  for (l = 1; l <= RES; l++) {
    if (TIMING)
      timer_Init(&timer);
    lenght = Lenstra(&B, l);
    if (TIMING) {
      printf("len = %d\n", lenght);
      printf("%d %6.5f\n", l, timer_Val(timer, timer_sec));
    }
    d = intmin(lenght, Deg / l);
    params->ecart[l] = Deg / l - d;
    params->se += params->ecart[l];
    if (params->ecart[l] > params->delta[l] || params->se > params->mse) {
      maxl = l;
      break;
    }
    if (l != RES)
      DualBaseIncrease(&B, gVect);
  }
  FreeBase(&B);
  FreeBV(&P);
  for (i = 0; i < RES; i++)
    FreeBV(&(gVect[i]));
  free(gVect);
  params->se = 0;
  for (l = 1; l <= maxl; l++) {
    if (params->ecart[l] == -1)
      params->ecart[l] = 0;
    params->se += params->ecart[l];
  }
  for (; l <= params->L; l++)
    if (params->ecart[l] == -1)
      params->ecart[l] = INT_MAX;
  params->meverified = TRUE;
}

/* Computes the inverse of g1 in the ring F2k and stores it in Invg1. */
void ComputeInvg1(BitVect *g1, ZENRing F2, ZENElt *Invg1, ZENRing *F2k) {
  ZENElt tempg;
  ZENEltAlloc(tempg, *F2k);
  BitVectToElt(&tempg, g1, *F2k, F2);
  ZENEltAlloc(*Invg1, *F2k);
  ZENEltInverse(*Invg1, tempg, *F2k);
  ZENEltFree(tempg, *F2k);
}

/* Computes the generating polynomials g_j(z) for j from startres to endres. */
void FindPolys(Combinaison *C, BitVect *M, BitVect BVPolys[], int startres, int endres) {
  BitVect OUTPUT, BC, OUT, State, TempBV;
  Matrix A;
  int i, k, j, K;

  K = C->k_g;
  AllocMat(&A, endres, K, 1);
  AllocBV(&OUTPUT, endres);
  AllocBV(&OUT, endres);
  PutBVToZero(&OUTPUT);

  /* Initialize each component */
  for (j = 0; j < C->J; j++) {
    AllocBV(&BC, DEGGEN(GEN(C, j)));
    AllocBV(&State, C->Components[j].Gen[C->Components[j].CurrentGen]->L);
    BVCanonic(&BC, 0);
    INITGEN(GEN(C, j), &BC, &State);
    TRANSFORME(C, j, &State, &OUT);
    XORBV(&OUTPUT, &OUTPUT, &OUT);
    FreeBV(&State);
    FreeBV(&BC);
  }

  /* Build matrix A */
  for (i = startres - 1; i < endres; i++)
    PutBitBV(&(A.lignes[i][0]), 0, ValBitBV(&OUTPUT, i));
  for (k = 1; k < K; k++) {
    PutBVToZero(&OUTPUT);
    for (j = 0; j < C->J; j++) {
      AllocBV(&State, C->Components[j].Gen[C->Components[j].CurrentGen]->L);
      ITERATION(GEN(C, j), &State);
      TRANSFORME(C, j, &State, &OUT);
      XORBV(&OUTPUT, &OUTPUT, &OUT);
      FreeBV(&State);
    }
    for (i = startres - 1; i < endres; i++)
      PutBitBV(&(A.lignes[i][0]), k, ValBitBV(&OUTPUT, i));
  }
  FreeBV(&OUT);
  FreeBV(&OUTPUT);

  /* Determine each polynomial */
  AllocBV(&TempBV, K + 1);
  for (i = startres - 1; i < endres; i++) {
    AllocBV(&BVPolys[i], K + 1);
    PutBVToZero(&BVPolys[i]);
    for (k = 0; k < K; k++) {
      if (ValBitBV(&(A.lignes[i][0]), k) == 1) {
        BVLShift(&TempBV, M, k + 1);
        XORBV(&BVPolys[i], &BVPolys[i], &TempBV);
      }
    }
    ANDBVMask(&BVPolys[i], &BVPolys[i], K);
  }
  FreeBV(&TempBV);
  FreeMat(&A);
}

/* Multiplies each BVPolys[j] by the inverse inv (in F2k) for j from startres to endres. */
void FindPolysTimesInv(BitVect BVPolys[], ZENRing F2, ZENElt *inv, ZENRing *F2k, int startres, int endres) {
  int j;
  ZENElt tempg, g;
  for (j = startres - 1; j < endres; j++) {
    ZENEltAlloc(tempg, *F2k);
    BitVectToElt(&tempg, &BVPolys[j], *F2k, F2);
    ZENEltAlloc(g, *F2k);
    ZENEltMultiply(g, tempg, *inv, *F2k);
    ZENEltFree(tempg, *F2k);
    EltToBitVect(&BVPolys[j], &g, *F2k, F2);
    ZENEltFree(g, *F2k);
  }
}

/* Computes the generating polynomials g_i(z) for each bit, normalised by g_1^{-1}. */
static void GetPolys(Combinaison *C, BitVect *M, BitVect BVPolys[], int resolution) {
  ZENRing F2, polyring;
  ZENElt inv;
  BigNum two;
  BigNumLength twol;

  ZBNReadFromString(&two, &twol, "2", 10);
  if (ZENBaseRingAlloc(F2, two, twol)) {
    printf("ERROR!\n");
    exit(1);
  }
  ZBNF(two);

  /* Find the generating functions g_i(z) for each bit */
  FindPolys(C, M, BVPolys, 1, resolution);
  /* Build the finite field F[z]/M(z) */
  ConstructZENRing(M, C->k_g, &polyring);
  /* Find the inverse of the first polynomial g_1(z) */
  ComputeInvg1(&BVPolys[0], F2, &inv, &polyring);
  /* Multiply each generating function by g_1^{-1}(z) */
  FindPolysTimesInv(BVPolys, F2, &inv, &polyring, 1, resolution);
  ZENEltFree(inv, polyring);
  ZENRingClose(polyring);
  ZENRingClose(F2);
}

/* Allocates a PolVect with B->maxresolution BitVect coefficients of degree B->degmax. */
static void AllocPolVect(Base *B, PolVect *P) {
  int j;
  P->coeffs = (BitVect *) malloc(B->maxresolution * sizeof(BitVect));
  for (j = 0; j < B->maxresolution; j++)
    AllocBV(&(P->coeffs[j]), B->degmax + 1);
  P->deg = INT_MIN;
  P->indicemaxdeg = INT_MIN;
}

/* Frees all BitVect coefficients of a PolVect and resets its degree fields. */
static void FreePolVect(Base *B, PolVect *P) {
  int j;
  for (j = 0; j < B->maxresolution; j++)
    FreeBV(&(P->coeffs[j]));
  free(P->coeffs);
  P->deg = INT_MIN;
  P->indicemaxdeg = INT_MIN;
}

/* Allocates a Base with given resolution and max polynomial degree. */
void AllocBase(Base *B, int resolution, int degmax) {
  int i;
  B->degmax = degmax;
  B->resolution = 0;
  B->maxresolution = resolution;
  B->vect = (PolVect *) malloc(resolution * sizeof(PolVect));
  B->permutations = (int *) malloc(resolution * sizeof(int));
  B->invpermutations = (int *) malloc(resolution * sizeof(int));
  for (i = 0; i < resolution; i++)
    AllocPolVect(B, &(B->vect[i]));
}

/* Copies Base B2 into Base B1 (both must have compatible dimensions). */
void CopyBase(Base *B1, Base *B2) {
  int j;
  B1->degmax = B2->degmax;
  B1->resolution = B2->resolution;
  B1->maxresolution = B2->maxresolution;
  for (j = 0; j < B2->maxresolution; j++)
    CopyPolVect(B2, &(B1->vect[j]), &(B2->vect[j]));
  for (j = 0; j < B2->maxresolution; j++) {
    B1->permutations[j] = B2->permutations[j];
    B1->invpermutations[j] = B2->invpermutations[j];
  }
}

/* Frees all memory allocated for a Base and resets its fields. */
void FreeBase(Base *B) {
  int i;
  for (i = 0; i < B->maxresolution; i++)
    FreePolVect(B, &(B->vect[i]));
  free(B->vect);
  free(B->permutations);
  free(B->invpermutations);
  B->resolution = 0;
  B->degmax = INT_MIN;
}

/* Sets the j-th (original, unpermuted) element of P to zero; assumes not yet used. */
static void SetElementVectorToZero(PolVect *P, int j) {
  int i;
  for (i = 0; i <= P->deg / WL; i++)
    P->coeffs[j].vect[i] = 0U;
}

/* Sets the j-th (original, unpermuted) element of P to the constant 1. */
static void SetElementVectorToOne(PolVect *P, int j) {
  SetElementVectorToZero(P, j);
  PutBitBV(&(P->coeffs[j]), 0, 1);
  if (P->deg == INT_MIN) {
    P->deg = 0;
    P->indicemaxdeg = j;
  }
}

/* Copies BitVect A into the j-th element of P and updates the degree if needed. */
static void BVToElementVector(Base *B, BitVect *A, PolVect *P, int j) {
  int deg;
  CopyBV(&(P->coeffs[j]), A);
  deg = NormeBitVect(B, &(P->coeffs[j]));
  if (P->deg < deg) {
    P->deg = deg;
    P->indicemaxdeg = j;
  }
}

/* Sets all elements of the vector P to zero (dense reset). */
static void VectorToZero(Base *B, PolVect *P) {
  int j;
  for (j = 0; j < B->resolution; j++)
    PutBVToZero(&(P->coeffs[j]));
}

/* Lazily marks the vector P as zero by resetting degree fields only. */
static inline void VectorToZeroLazy(PolVect *P) {
  P->deg = INT_MIN;
  P->indicemaxdeg = INT_MIN;
}

/* Initializes the dual basis of resolution res from the row polynomials R0 and char poly M. */
void DualBase(Base *B, BitVect *R0, BitVect *M, int res) {
  int i;
  if (res > B->maxresolution) {
    printf("Error DualBase()!\n");
    exit(1);
  }
  for (i = 0; i < B->maxresolution; i++) {
    B->permutations[i] = i;
    B->invpermutations[i] = i;
  }
  B->resolution = res;
  VectorToZeroLazy(&(vb(0)));
  BVToElementVector(B, M, &(vb(0)), 0);

  /* Build the remaining dual basis vectors */
  for (i = 1; i < res; i++) {
    VectorToZeroLazy(&(vb(i)));
    BVToElementVector(B, &R0[i], &(vb(i)), 0);
    SetElementVectorToOne(&(vb(i)), i);
  }
}

/* Increases the dual basis resolution by one, adding the next row polynomial R0. */
void DualBaseIncrease(Base *B, BitVect *R0) {
  int j;
  if (B->resolution == B->maxresolution) {
    printf("Error UpdateDualBase()!\n");
    exit(1);
  }
  B->resolution++;
  for (j = 0; j <= B->resolution - 2; j++)
    SetElementVectorToZero(&(vb(j)), B->resolution - 1);
  VectorToZeroLazy(&(vb(B->resolution - 1)));
  BVToElementVector(B, &R0[B->resolution - 1], &(vb(B->resolution - 1)), 0);
  SetElementVectorToOne(&(vb(B->resolution - 1)), B->resolution - 1);
}

/* Swaps the coordinate at position m+1 with the coordinate of the leading-degree element. */
static inline void PermuteCoord(Base *B, int m) {
  int j, temp;

  j = B->invpermutations[(vb(m + 1)).indicemaxdeg];
  if (j != (m + 1)) {
    temp = vperm(m + 1);
    vperm(m + 1) = vperm(j);
    vperm(j) = temp;
    temp = B->invpermutations[vperm(j)];
    B->invpermutations[vperm(j)] = B->invpermutations[vperm(m + 1)];
    B->invpermutations[vperm(m + 1)] = temp;
  }
}

/* Runs Lenstra's algorithm on basis B up to resolution res and returns the first successive minimum. */
int Lenstra(Base *B, int res) {
  int i, s, m, *x, NormeTbmp1, Normebmp1, Normebi, resolution;
  int L, OK, newm, swap;
  PolVect sum, temp;
  resolution = res - 1;
  m = -1;
  if (resolution > 0)
    x = (int*) malloc((resolution) * sizeof(int));
  AllocPolVect(B, &temp);
  AllocPolVect(B, &sum);
  while (m < resolution) {
    renumber(B, m, resolution);
    if (m > -1)
      SolveAxb(B, m, x);
    Normebmp1 = NormeVect(&(vb(m + 1)));
    VectorToZeroLazy(&sum);
    for (i = 0; i <= m; i++)
      if (x[i]) {
        Normebi = NormeVect(&(vb(i)));
        s = Normebmp1 - Normebi;
        AddSelfPolVectMult(B, &sum, &(vb(i)), s);
      }
    if (sum.deg == INT_MIN) {
      CopyPolVectLazy(&temp, &(vb(m + 1)));
      swap = 1;
    } else {
      AddPolVect(B, &temp, &(vb(m + 1)), &sum);
      swap = 0;
    }

    NormeTbmp1 = NormeVect(&temp);
    if (NormeTbmp1 == Normebmp1) {
      CopyPolVectLazy(&(vb(m + 1)), &temp);
      PermuteCoord(B, m);
      m++;
    } else {
      if (NormeTbmp1 < Normebmp1) {
        CopyPolVectLazy(&(vb(m + 1)), &temp);
        OK = 0;
        for (L = 0; L <= m; L++)
          if (NormeVect(&(vb(L))) <= NormeTbmp1) {
            newm = L;
            OK = 1;
          }
        if (!OK)
          m = -1;
        else
          m = newm;
      } else if (swap == 1) {
        CopyPolVectLazy(&(vb(m + 1)), &temp);
      }
    }
  }
  if (resolution > 0)
    free(x);
  FreePolVect(B, &sum);
  FreePolVect(B, &temp);
  return NormeVect(&(vb(0)));
}

/* Reorders basis vectors from m+1 to resolution so that the smallest-norm vector is at m+1. */
static inline void renumber(Base *B, int m, int resolution) {
  int j, normebj, minj, min;
  PolVect T;
  min = INT_MAX;
  for (j = m + 1; j <= resolution; j++) {
    normebj = NormeVect(&(vb(j)));
    if (normebj < min) {
      min = normebj;
      minj = j;
    }
  }
  if (minj != m + 1) {
    T.deg = B->vect[m + 1].deg;
    B->vect[m + 1].deg = B->vect[minj].deg;
    B->vect[minj].deg = T.deg;
    T.coeffs = B->vect[m + 1].coeffs;
    B->vect[m + 1].coeffs = B->vect[minj].coeffs;
    B->vect[minj].coeffs = T.coeffs;
    T.indicemaxdeg = B->vect[m + 1].indicemaxdeg;
    B->vect[m + 1].indicemaxdeg = B->vect[minj].indicemaxdeg;
    B->vect[minj].indicemaxdeg = T.indicemaxdeg;
  }
}

/* Returns 1 if the j-th element of P has a nonzero bit at position norme, 0 otherwise. */
static inline int ElementVectorNormeEqual(Base *B, PolVect *P, int j, int norme) {
  if (ValBitBV(&(P->coeffs[vperm(j)]), norme) == 1)
    return 1;
  else
    return 0;
}

/* Solves the system A*x = b (mod 2) needed by Lenstra's algorithm for basis vectors 0..m. */
static inline void SolveAxb(Base *B, int m, int *x) {
  int i, j, sum;
  Matrix A;

  AllocMat(&A, m + 1, m + 2, 1);
  for (j = 0; j <= m; j++)
    PutBVToZero(&(A.lignes[j][0]));
  for (i = 0; i <= m; i++) {
    for (j = i; j <= m; j++)
      if (ElementVectorNormeEqual(B, &(B->vect[i]), j, NormeVect(&(B->vect[i]))))
        PutBitBV(&(A.lignes[j][0]), i, 1);
  }
  i = m + 1;
  for (j = 0; j <= m; j++) {
    if (ElementVectorNormeEqual(B, &(B->vect[i]), j, NormeVect(&(B->vect[i]))))
      PutBitBV(&(A.lignes[j][0]), i, 1);
  }
  for (j = 0; j <= m; j++)
    x[j] = 0;
  for (j = 0; j <= m; j++) {
    sum = ValBitBV(&(A.lignes[j][0]), m + 1);
    for (i = 0; i < j; i++)
      sum += ValBitBV(&(A.lignes[j][0]), i) * x[i];
    x[j] = sum & 1U;
  }
  FreeMat(&A);
}

/* Updates the degree of a after an in-place modification by scanning from the current degree downward. */
static inline void UpdatePolVectDeg(Base *B, PolVect *a) {
  int k, j, ll, K, LL;
  extern uint32_t MMC[WL];
  uint32_t BLOB;
  K = (a->deg) / WL;
  LL = a->deg - WL * K;

  for (k = (a->deg) / WL; k >= 0; k--) {
    BLOB = 0U;
    for (j = 0; j < B->resolution; j++)
      BLOB |= a->coeffs[j].vect[k];
    for (ll = LL; ll >= 0; ll--)
      if (BLOB & MMC[ll]) {
        a->deg = k * WL + ll;
        for (j = 0; j < B->resolution; j++)
          if (a->coeffs[j].vect[k] & MMC[ll]) {
            a->indicemaxdeg = j;
            return;
          }
      }
    LL = WL - 1;
  }
}

/* Computes res = a XOR b (polynomial vector addition over F_2). */
static inline void AddPolVect(Base *B, PolVect *res, PolVect *a, PolVect *b) {
  int k, j, K, LL;

  if (a->deg == b->deg) {
    res->deg = a->deg;
    for (k = a->deg / WL; k >= 0; k--)
      for (j = 0; j < B->resolution; j++)
        res->coeffs[j].vect[k] = a->coeffs[j].vect[k] ^ b->coeffs[j].vect[k];
    UpdatePolVectDeg(B, res);
    return;
  }
  if (a->deg > b->deg) {
    res->deg = a->deg;
    res->indicemaxdeg = a->indicemaxdeg;
    K = b->deg / WL;
    LL = b->deg - K * WL;
    for (j = 0; j < B->resolution; j++) {
      for (k = a->deg / WL; k > K; k--)
        res->coeffs[j].vect[k] = a->coeffs[j].vect[k];
      b->coeffs[j].vect[k] &= 0xffffffffU << (WL - LL - 1);
      for (; k >= 0; k--)
        res->coeffs[j].vect[k] = a->coeffs[j].vect[k] ^ b->coeffs[j].vect[k];
    }
    return;
  }
}

/* Computes r = a * X^S (left shift of polynomial vector by S bit positions). */
static inline void Multbys(Base *B, PolVect *r, PolVect *a, int S) {
  int i, ii, nbblocks, nbbits;
  int WLmn;
  int N;

  N = (a->deg + S) / WL;
  nbblocks = S / WL;
  nbbits = S - nbblocks * WL;

  if (nbbits != 0) {
    WLmn = WL - nbbits;
    for (ii = 0; ii < B->resolution; ii++) {
      for (i = 0; i < nbblocks; i++)
        r->coeffs[ii].vect[i] = 0U;
      r->coeffs[ii].vect[i++] = a->coeffs[ii].vect[0] >> nbbits;
      for (; i <= N; i++)
        r->coeffs[ii].vect[i] = (a->coeffs[ii].vect[i - nbblocks] >> nbbits) ^ (a->coeffs[ii].vect[i - nbblocks - 1] << (WLmn));
    }
  } else {
    for (ii = 0; ii < B->resolution; ii++) {
      for (i = 0; i < nbblocks; i++)
        r->coeffs[ii].vect[i] = 0U;
      r->coeffs[ii].vect[i++] = a->coeffs[ii].vect[0];
      for (; i <= N; i++)
        r->coeffs[ii].vect[i] = (a->coeffs[ii].vect[i - nbblocks]);
    }
  }
  r->deg = a->deg + S;
  r->indicemaxdeg = a->indicemaxdeg;
}

/* Computes res XOR= a * X^S in-place, handling degree cases relative to res->deg. */
static inline void Multbys2(Base *B, PolVect *r, PolVect *a, int S) {
  int i, ii, nbblocks, nbbits;
  int WLmn;
  int N;

  nbblocks = S / WL;
  nbbits = S - nbblocks * WL;

  if (a->deg + S == r->deg) {
    N = (r->deg) / WL;
    if (nbbits != 0) {
      WLmn = WL - nbbits;
      for (ii = 0; ii < B->resolution; ii++) {
        for (i = 0; i < nbblocks; i++)
          r->coeffs[ii].vect[i] ^= 0U;
        r->coeffs[ii].vect[i++] ^= a->coeffs[ii].vect[0] >> nbbits;
        for (; i <= N; i++)
          r->coeffs[ii].vect[i] ^= (a->coeffs[ii].vect[i - nbblocks] >> nbbits) ^ (a->coeffs[ii].vect[i - nbblocks - 1] << (WLmn));
      }
    } else {
      for (ii = 0; ii < B->resolution; ii++) {
        for (i = 0; i < nbblocks; i++)
          r->coeffs[ii].vect[i] ^= 0U;
        r->coeffs[ii].vect[i++] ^= a->coeffs[ii].vect[0];
        for (; i <= N; i++)
          r->coeffs[ii].vect[i] ^= (a->coeffs[ii].vect[i - nbblocks]);
      }
    }
    r->indicemaxdeg = a->indicemaxdeg;
    UpdatePolVectDeg(B, r);
    return;
  }
  if (a->deg + S < r->deg) {
    N = (a->deg + S) / WL;
    if (nbbits != 0) {
      WLmn = WL - nbbits;
      for (ii = 0; ii < B->resolution; ii++) {
        for (i = 0; i < nbblocks; i++)
          r->coeffs[ii].vect[i] ^= 0U;
        r->coeffs[ii].vect[i++] ^= a->coeffs[ii].vect[0] >> nbbits;
        for (; i < N; i++)
          r->coeffs[ii].vect[i] ^= (a->coeffs[ii].vect[i - nbblocks] >> nbbits) ^ (a->coeffs[ii].vect[i - nbblocks - 1] << (WLmn));
        r->coeffs[ii].vect[i] ^= (a->coeffs[ii].vect[i - nbblocks - 1] << (WLmn));
      }
    } else {
      for (ii = 0; ii < B->resolution; ii++) {
        for (i = 0; i < nbblocks; i++)
          r->coeffs[ii].vect[i] ^= 0U;
        r->coeffs[ii].vect[i++] ^= a->coeffs[ii].vect[0];
        for (; i <= N; i++)
          r->coeffs[ii].vect[i] ^= (a->coeffs[ii].vect[i - nbblocks]);
      }
    }
    r->deg = a->deg + S;
    r->indicemaxdeg = a->indicemaxdeg;
    UpdatePolVectDeg(B, r);
    return;
  }
}

/* Computes res = res + a * X^s, allocating lazily if res is currently zero. */
static inline void AddSelfPolVectMult(Base *B, PolVect *res, PolVect *a, int s) {
  if (res->deg == INT_MIN)
    Multbys(B, res, a, s);
  else
    Multbys2(B, res, a, s);
}

/* Deep-copies polynomial vector b into a (uses maxresolution and degmax). */
static inline void CopyPolVect(Base *B, PolVect *a, PolVect *b) {
  int ii, k;
  for (ii = 0; ii < B->maxresolution; ii++)
    for (k = 0; k <= B->degmax / WL; k++)
      a->coeffs[ii].vect[k] = b->coeffs[ii].vect[k];
  a->deg = b->deg;
  a->indicemaxdeg = b->indicemaxdeg;
}

/* Swaps the contents of polynomial vectors a and b by exchanging their pointers and degrees. */
static inline void CopyPolVectLazy(PolVect *a, PolVect *b) {
  PolVect T;

  T.deg = b->deg;
  b->deg = a->deg;
  a->deg = T.deg;

  T.coeffs = b->coeffs;
  b->coeffs = a->coeffs;
  a->coeffs = T.coeffs;

  T.indicemaxdeg = b->indicemaxdeg;
  b->indicemaxdeg = a->indicemaxdeg;
  a->indicemaxdeg = T.indicemaxdeg;
}

/* Displays the j-th element of polynomial vector P (using current permutation). */
static void DispElementVector(Base *B, PolVect *P, int j) {
  int k;
  if (P->deg == INT_MIN) {
    printf("Vecteur 0");
    return;
  }
  DispBitVect(&(P->coeffs[vperm(j)]), P->deg + 1);
  for (k = P->deg + 1; k <= B->degmax; k++)
    printf("0");
}

/* Displays the entire basis B up to the given resolution. */
static void DispBase(Base *B, int resolution) {
  register int i, j, n;

  printf("\n");
  for (i = 0; i < resolution; i++)
    printf("%d ", B->permutations[i]);
  printf("\n");
  for (i = 0; i < resolution; i++) {
    for (j = 0; j < resolution; j++) {
      printf("[");
      DispElementVector(B, &(vb(i)), j);
      n = NormeElementVector(B, &(vb(i)), j);
      if (n < 0)
        printf("]-I   ");
      else
        printf("]%2d   ", n);
    }
    printf(" |b%d|=%d\n", i + 1, NormeVect(&(vb(i))));
  }
}

/* Displays all elements of polynomial vector P under the current basis resolution. */
static void DispVector(Base *B, PolVect *P) {
  int j;
  for (j = 0; j < B->resolution; j++) {
    DispElementVector(B, P, j);
    printf("\n");
  }
  printf("\n");
}

/* Converts a ZENPoly p (in F2) to a ZENElt e in the polynomial ring polring. */
static void PolyToElt(ZENElt *e, ZENRing polring, ZENPoly p, ZENRing F2) {
  char *polystring;
  polystring = ZENPolyPutToString(p, F2);
  ZENEltGetFromString(*e, polystring, polring);
  free(polystring);
}

/* Converts a ZENElt e (in polring) to a ZENPoly p over F2. */
static void EltToPoly(ZENPoly *p, ZENRing F2, ZENElt e, ZENRing polring) {
  char *polystring;
  polystring = ZENEltPutToString(e, polring);
  ZENPolyGetFromString(*p, polystring, F2);
  free(polystring);
}

/* Converts a BitVect A to a ZENElt e in the given polynomial ring polring over F2. */
void BitVectToElt(ZENElt *e, BitVect *A, ZENRing polring, ZENRing F2) {
  int D, flag, k;
  ZENPoly temp;
  ZENElt ONE;

  ZENEltAlloc(ONE, F2);
  ZENEltSetToOne(ONE, F2);
  ZENPolyAlloc(temp, (D = ZENRingDeg(polring)) - 1, F2);
  flag = 0;
  for (k = D - 1; k >= 0; k--)
    if (ValBitBV(A, k) == 1) {
      if (flag == 0) {
        ZENPolySetToXi(temp, k, F2);
        flag = 1;
      } else {
        ZENPolySetCoeff(temp, k, ONE, F2);
      }
    }

  PolyToElt(e, polring, temp, F2);
  ZENPolyFree(temp, F2);
  ZENEltFree(ONE, F2);
}

/* Converts a ZENElt e (in polring over F2) back to a BitVect A. */
void EltToBitVect(BitVect *A, ZENElt *e, ZENRing polring, ZENRing F2) {
  int k;
  ZENPoly temp;
  ZENElt Coeff;

  ZENPolyAlloc(temp, ZENRingDeg(polring) - 1, F2);
  EltToPoly(&temp, F2, *e, polring);
  PutBVToZero(A);
  for (k = ZENPolyDeg(temp, F2); k >= 0; k--) {
    Coeff = ZENPolyGetCoeffPtr(temp, k, F2);
    if (ZENEltIsZero(Coeff, F2))
      PutBitBV(A, k, 0);
    else
      PutBitBV(A, k, 1);
  }
  ZENPolyFree(temp, F2);
}

/* Returns the degree of the j-th element of P (using current permutation), or INT_MIN if zero. */
static int NormeElementVector(Base *B, PolVect *P, int j) {
  int k;
  k = P->deg;
  while (k >= 0) {
    if (ValBitBV(&(P->coeffs[vperm(j)]), k) == 0)
      k--;
    else
      break;
  }
  if (k < 0)
    return INT_MIN;
  else
    return k;
}

/* Returns the degree of BitVect A within the basis B (highest set bit), or INT_MIN if zero. */
static int NormeBitVect(Base *B, BitVect *A) {
  int k;
  k = B->degmax;
  while (k >= 0) {
    if (ValBitBV(A, k) == 0)
      k--;
    else
      break;
  }
  if (k < 0)
    return INT_MIN;
  else
    return k;
}

/* Same as TestMELat but with a rotated coordinate mapping (rotative_shift). */
void TestMELatSpecial(Combinaison *C, paramMecf *params, int rotative_shift) {
  Base B;
  BitVect *gVect;
  BitVect P;
  timer_Chrono timer;
  int Deg, l, d, *coeff, lenght, i, RES, maxl;

  if (params->meverif == FALSE)
    return;

  if (params->L < C->L) {
    printf("Error in TestME().\n  Resolution of allocated paramMecf < resolution of generator\n");
    exit(1);
  }
  if (TIMING)
    timer_Init(&timer); /* Initialize the timer */
  Deg = C->k_g;
  RES = intmin(params->L, Deg);
  for (l = 1; l <= params->L; l++) /* Reset the ecart array */
    params->ecart[l] = -1;

  gVect = (BitVect *) malloc(RES * sizeof(BitVect));
  for (i = 0; i < RES; i++)
    AllocBV(&(gVect[i]), Deg + 1);
  coeff = (int *) malloc((Deg + 1) * sizeof(int));
  AllocBV(&P, Deg + 1);
  polycharComb(C, coeff, &P);
  free(coeff);

  GetPolysSpecial(C, &P, gVect, RES, rotative_shift);

  AllocBase(&B, RES, Deg);
  B.degmax = Deg;
  DualBase(&B, gVect, &P, 1);
  SetPsi12(C, params);
  params->se = 0;
  maxl = params->L;
  for (l = 1; l <= RES; l++) {
    lenght = Lenstra(&B, l);
    d = intmin(lenght, Deg / l);
    params->ecart[l] = Deg / l - d;
    params->se += params->ecart[l];
    if (l != RES)
      DualBaseIncrease(&B, gVect);
  }
  FreeBase(&B);
  FreeBV(&P);
  for (i = 0; i < RES; i++)
    FreeBV(&(gVect[i]));
  free(gVect);
  params->se = 0;
  for (l = 1; l <= maxl; l++) {
    if (params->ecart[l] == -1)
      params->ecart[l] = 0;
    params->se += params->ecart[l];
  }
  for (; l <= params->L; l++)
    if (params->ecart[l] == -1)
      params->ecart[l] = INT_MAX;
  params->meverified = TRUE;
}

/* Same as GetPolys but applies a rotative_shift to the coordinate mapping before normalization. */
static void GetPolysSpecial(Combinaison *C, BitVect *M, BitVect BVPolys[], int resolution, int rotative_shift) {
  ZENRing F2, polyring;
  ZENElt inv;
  BigNum two;
  BigNumLength twol;
  BitVect *Temp;
  int i;

  Temp = (BitVect *) malloc(resolution * sizeof(BitVect));
  for (i = 0; i < resolution; i++)
    AllocBV(&Temp[i], C->k_g + 1);

  ZBNReadFromString(&two, &twol, "2", 10);
  if (ZENBaseRingAlloc(F2, two, twol)) {
    printf("ERROR!\n");
    exit(1);
  }
  ZBNF(two);

  /* Find the generating functions g_i(z) for each bit */
  FindPolys(C, M, Temp, 1, resolution);

  for (i = 0; i < resolution; i++)
    CopyBV(&BVPolys[i], &Temp[(i + rotative_shift) % resolution]);

  /* Build the finite field F[z]/M(z) */
  ConstructZENRing(M, C->k_g, &polyring);
  /* Find the inverse of the first polynomial g_1(z) */
  ComputeInvg1(&BVPolys[0], F2, &inv, &polyring);
  /* Multiply each generating function by g_1^{-1}(z) */
  FindPolysTimesInv(BVPolys, F2, &inv, &polyring, 1, resolution);
  ZENEltFree(inv, polyring);
  ZENRingClose(polyring);
  ZENRingClose(F2);
}
