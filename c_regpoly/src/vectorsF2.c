/*
 * vectorsF2.c — Bit vector and matrix operations over GF(2).
 *
 * Provides allocation, arithmetic, shift, logical, and display operations
 * on multi-word bit vectors (BitVect) and matrices of bit vectors (Matrix),
 * including Gaussian elimination and matrix inversion over GF(2).
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <vectorsF2.h>

/* High bit of a 32-bit word; used as the most-significant-bit mask for
   row pivot searches in elimination routines. */
#define MC 0x80000000U

uint32_t MMC[WL] = {
  MC,    MC>>1,  MC>>2,  MC>>3,  MC>>4,  MC>>5,  MC>>6,  MC>>7,
  MC>>8,  MC>>9,  MC>>10, MC>>11, MC>>12, MC>>13, MC>>14, MC>>15,
  MC>>16, MC>>17, MC>>18, MC>>19, MC>>20, MC>>21, MC>>22, MC>>23,
  MC>>24, MC>>25, MC>>26, MC>>27, MC>>28, MC>>29, MC>>30, MC>>31
};

/* Computes the inverse of the square matrix M over GF(2) using augmented
   Gaussian elimination, storing the result in InvM. Returns TRUE if M is
   invertible (full rank), FALSE otherwise. */
boolean InverseMatrix(Matrix *InvM, Matrix *M)
{
  Matrix Temp;
  int j, rang;

  if (M->nblignes != M->l) {
    printf("Matrix M is not square!\n");
    exit(1);
  }
  AllocMat(&Temp, M->nblignes, M->l, 2);
  for (j = 0; j < M->l; j++)
    CopyBV(&(Temp.lignes[j][0]), &(M->lignes[j][0]));
  for (j = 0; j < M->l; j++)
    BVCanonic(&(Temp.lignes[j][1]), j);
  rang = CompleteElimination(&Temp, M->nblignes, M->l, 2);
  for (j = 0; j < M->l; j++)
    CopyBV(&(InvM->lignes[j][0]), &(Temp.lignes[j][1]));
  FreeMat(&Temp);
  return (rang == M->l);
}

/* Performs complete (reduced) Gaussian elimination on the first nblignes
   rows of matrix m, treating each row as t groups of l bits. Eliminates
   both above and below each pivot, returning the rank of the matrix. */
int CompleteElimination(Matrix *m, int nblignes, int l, int t)
{
  int i, j, cl, rang;

  rang = 0;
  j = 0;
  while (j < t) {
    cl = 0;
    while (cl < l) {
      /* Search column j from row 'rang' for a pivot at bit position cl. */
      i = rang;
      while ((i < nblignes) && (m->lignes[i][j].vect[cl/WL] < MMC[cl%WL]))
        i++;
      if (i < nblignes) {
        ExchangeVect(m, rang, i);
        for (i = 0; i < nblignes; i++)
          if (i != rang)
            if (m->lignes[i][j].vect[cl/WL] & MMC[cl%WL])
              XorVect(m, i, rang, 0, m->t);
        rang++;
        if (rang == nblignes)
          return rang;
      } else {
        return rang;
      }
      cl++;
    }
    j++;
  }
  return rang;
}

/* Returns the value (0 or 1) of bit noBit in bit vector A.
   Bits are indexed from 0 (most significant). */
int ValBitBV(BitVect *A, int noBit)
{
  int k;
  uint32_t mask;

  k = noBit / WL;
  mask = 0x80000000U >> (noBit - k * WL);
  if (A->vect[k] & mask)
    return 1;
  else
    return 0;
}

/* Sets bit noBit in bit vector A to valBit (0 or 1).
   Bits are indexed from 0 (most significant). */
void PutBitBV(BitVect *A, int noBit, int valBit)
{
  int k;
  uint32_t mask;

  k = noBit / WL;
  if (valBit == 1) {
    mask = 0x80000000U >> (noBit - k * WL);
    A->vect[k] |= mask;
  } else {
    mask = 0xffffffffU ^ (0x80000000U >> (noBit - k * WL));
    A->vect[k] &= mask;
  }
}

/* Sets all words of bit vector A to zero. */
void PutBVToZero(BitVect *A)
{
  int i;

  for (i = 0; i < A->n; i++)
    A->vect[i] = 0U;
}

/* Copies bit vector B into A (A = B). Both must have the same word count. */
void CopyBV(BitVect *A, BitVect *B)
{
  int i;

  if (A->n != B->n) {
    printf("Error in CopyBV(): vectors of different dimensions! (%d and %d bits)\n",
           A->n * WL, B->n * WL);
    exit(1);
  }
  if (B->n == 0) {
    printf("Nothing to copy!\n");
    exit(1);
  }
  for (i = 0; i < B->n; i++)
    A->vect[i] = B->vect[i];
}

/* Copies the first l bits of B into A, masking out any higher bits.
   A must be large enough to hold ceil(l/WL) words. */
void CopyBVPart(BitVect *A, BitVect *B, int l)
{
  int i, n;

  n = (l - 1) / WL + 1;
  if (A->n < n) {
    printf("Error in CopyBVPart() : The vector A is not large enough!\n");
    exit(1);
  }
  if (B->n == 0) {
    printf("Nothing to copy!\n");
    exit(1);
  }
  for (i = 0; i < n; i++)
    A->vect[i] = B->vect[i];
  if (l % WL) {
    BitVect m;
    AllocBV(&m, A->n * WL);
    mask(&m, l);
    ANDBVSelf(A, &m);
    FreeBV(&m);
  }
}

/* Computes A = B XOR C (bitwise). All three must have the same word count. */
void XORBV(BitVect *A, BitVect *B, BitVect *C)
{
  int i;

  if ((A->n != B->n) || (B->n != C->n)) {
    printf("Error in XORBV(): Vectors of different sizes\n");
    exit(1);
  }
  for (i = 0; i < B->n; i++)
    A->vect[i] = B->vect[i] ^ C->vect[i];
}

/* Computes A = B XOR C XOR D (bitwise). All four must have the same word count. */
void XOR2BV(BitVect *A, BitVect *B, BitVect *C, BitVect *D)
{
  int i;

  if ((A->n != B->n) || (B->n != C->n) || (C->n != D->n)) {
    printf("Error in XOR2BV(): Vectors of different sizes\n");
    exit(1);
  }
  for (i = 0; i < B->n; i++)
    A->vect[i] = B->vect[i] ^ C->vect[i] ^ D->vect[i];
}

/* Computes A = B AND C (bitwise). All three must have the same word count. */
void ANDBV(BitVect *A, BitVect *B, BitVect *C)
{
  int i;

  if ((A->n != B->n) || (B->n != C->n)) {
    printf("Error in ANDBV(): Vectors of different sizes\n");
    exit(1);
  }
  for (i = 0; i < B->n; i++)
    A->vect[i] = B->vect[i] & C->vect[i];
}

/* Copies B into A keeping only the first t bits (masks out bits t and above). */
void ANDBVMask(BitVect *A, BitVect *B, int t)
{
  int n, m, j;

  if (A->n != B->n) {
    printf("Error in ANDBVMask(): Vectors of different sizes\n");
    exit(1);
  }
  if (t > B->n * WL) {
    CopyBV(A, B);
  } else if (t == 0) {
    PutBVToZero(A);
  } else {
    n = t / WL;
    m = t - n * WL;
    for (j = 0; j < n; j++)
      A->vect[j] = B->vect[j];
    if (m != 0) {
      A->vect[j] = B->vect[j] & (0xffffffffU << (WL - m));
      j++;
    }
    for (; j < A->n; j++)
      A->vect[j] = 0U;
  }
}

/* Copies B into A keeping only bits from position t onward (masks out bits
   0..t-1), i.e., the complement operation to ANDBVMask. */
void ANDBVInvMask(BitVect *A, BitVect *B, int t)
{
  int n, m, j;

  if (A->n != B->n) {
    printf("Error in ANDBV(): Vectors of different sizes\n");
    exit(1);
  }
  if (t > B->n * WL) {
    PutBVToZero(A);
  } else if (t == 0) {
    CopyBV(A, B);
  } else {
    n = t / WL;
    m = t - n * WL;
    for (j = 0; j < n; j++)
      A->vect[j] = 0U;
    if (m == 0)
      A->vect[j] = B->vect[j];
    else
      A->vect[j] = B->vect[j] & (0xffffffffU >> m);
    j++;
    for (; j < A->n; j++)
      A->vect[j] = B->vect[j];
  }
}

/* Computes A &= B (bitwise AND in place). Both must have the same word count. */
void ANDBVSelf(BitVect *A, BitVect *B)
{
  int i;

  if (A->n != B->n) {
    printf("Error in ANDBVSelf(): Vectors of different sizes\n");
    exit(1);
  }
  for (i = 0; i < B->n; i++)
    A->vect[i] &= B->vect[i];
}

/* Computes A ^= B (bitwise XOR in place). Both must have the same word count. */
void XORBVSelf(BitVect *A, BitVect *B)
{
  int i;

  if (A->n != B->n) {
    printf("Error in XORBVSelf(): Vectors of different sizes\n");
    exit(1);
  }
  for (i = 0; i < B->n; i++)
    A->vect[i] ^= B->vect[i];
}

/* Computes R = A << n (logical left shift by n bit positions). R and A must
   have the same word count. */
void BVLShift(BitVect *R, BitVect *A, int n)
{
  int i;
  int WLmn;
  uint32_t temp;

  if (R->n != A->n) {
    printf("Error in BVLShift(): Vectors of different sizes\n");
    exit(1);
  }
  for (i = 0; i < A->n; i++)
    R->vect[i] = A->vect[i];
  while (n >= 32) {
    for (i = 1; i < A->n; i++)
      R->vect[i - 1] = R->vect[i];
    R->vect[A->n - 1] = 0U;
    n -= 32;
  }
  if (n > 0) {
    WLmn = WL - n;
    R->vect[0] <<= n;
    for (i = 1; i < A->n; i++) {
      temp = R->vect[i] >> WLmn;
      R->vect[i - 1] |= temp;
      R->vect[i] <<= n;
    }
  }
}

/* Rotates the first w bits of A left by n positions, leaving bits above
   position w unchanged, and stores the result in R. */
void BVLRotativeShift(BitVect *R, BitVect *A, int n, int w)
{
  BitVect Temp;
  int W;

  W = A->n * WL;
  AllocBV(&Temp, W);
  BVLShift(&Temp, A, n);       /* Temp = A << n */
  BVRShift(R, A, w - n);       /* R = A >> (w-n) */
  ANDBVMask(&Temp, &Temp, w - n); /* Temp &= mask(w-n): keep first w-n bits */
  XORBVSelf(R, &Temp);         /* R ^= Temp */
  ANDBVMask(R, R, w);          /* R &= mask(w): keep first w bits */
  ANDBVInvMask(&Temp, A, w);   /* Temp = A & ~mask(w): bits above w */
  XORBVSelf(R, &Temp);         /* R ^= Temp: restore bits above w */
  FreeBV(&Temp);
}

/* Computes R = A >> n (logical right shift by n bit positions). R and A must
   have the same word count. */
void BVRShift(BitVect *R, BitVect *A, int n)
{
  int i;
  int WLmn;
  uint32_t temp;

  if (R->n != A->n) {
    printf("Error in BVRShift(): Vectors of different sizes\n");
    exit(1);
  }
  for (i = 0; i < A->n; i++)
    R->vect[i] = A->vect[i];
  while (n >= 32) {
    for (i = A->n; i > 1; i--)
      R->vect[i - 1] = R->vect[i - 2];
    R->vect[0] = 0U;
    n -= 32;
  }
  if (n > 0) {
    WLmn = WL - n;
    R->vect[A->n - 1] >>= n;
    for (i = A->n - 2; i >= 0; i--) {
      temp = R->vect[i] << WLmn;
      R->vect[i + 1] |= temp;
      R->vect[i] >>= n;
    }
  }
}

/* Computes R <<= n (left shift in place by n bit positions). */
void BVLShiftSelf(BitVect *R, int n)
{
  int i;
  int WLmn;
  uint32_t temp;

  while (n >= 32) {
    for (i = 1; i < R->n; i++)
      R->vect[i - 1] = R->vect[i];
    R->vect[R->n - 1] = 0U;
    n -= 32;
  }
  if (n > 0) {
    WLmn = WL - n;
    R->vect[0] <<= n;
    for (i = 1; i < R->n; i++) {
      temp = R->vect[i] >> WLmn;
      R->vect[i - 1] |= temp;
      R->vect[i] <<= n;
    }
  }
}

/* Computes R <<= 1 (left shift in place by 1). Specialised version of
   BVLShiftSelf for single-bit shifts used in inner loops. */
void BVLS1Self(BitVect *R)
{
  int i;

  R->vect[0] <<= 1;
  for (i = 1; i < R->n; i++) {
    if (R->vect[i] & MC)
      R->vect[i - 1] |= 0x1U;
    R->vect[i] <<= 1;
  }
}

/* Computes R >>= n (right shift in place by n bit positions). */
void BVRShiftSelf(BitVect *R, int n)
{
  int i;
  int WLmn;
  uint32_t temp;

  while (n >= 32) {
    for (i = R->n - 1; i > 0; i--)
      R->vect[i] = R->vect[i - 1];
    R->vect[0] = 0U;
    n -= 32;
  }
  if (n > 0) {
    WLmn = WL - n;
    R->vect[R->n - 1] >>= n;
    for (i = R->n - 2; i >= 0; i--) {
      temp = R->vect[i] << WLmn;
      R->vect[i + 1] |= temp;
      R->vect[i] >>= n;
    }
  }
}

/* Bitwise complement of A in place (A = ~A). */
void InverseBV(BitVect *A)
{
  int i;

  for (i = 0; i < A->n; i++)
    A->vect[i] = ~A->vect[i];
}

/* Returns TRUE if bit vectors ds1 and ds2 share at least one common set bit
   (i.e., their bitwise AND is non-zero). */
boolean VerifBitsCommuns(BitVect *ds1, BitVect *ds2)
{
  int i;
  uint32_t temp = 0U;

  if (ds1->n != ds2->n) {
    printf("Error in VerifBitsCommuns(): Vectors of different sizes\n");
    exit(1);
  }
  for (i = 0; i < ds1->n; i++)
    temp |= (ds1->vect[i] & ds2->vect[i]);
  if (temp)
    return TRUE;
  else
    return FALSE;
}

/* --- Matrix manipulation functions --- */

/* Sets A to the l-th canonical basis vector (all zeros except bit l = 1). */
void BVCanonic(BitVect *A, int l)
{
  int n;

  PutBVToZero(A);
  n = l / WL;
  if (n > A->n) {
    printf("Error in BVCanonic(): vector A is not long enough to store BVCanonic[%d].\n", l);
    exit(1);
  }
  A->vect[n] = 0x80000000U >> (l - n * WL);
}

/* Sets A to the mask with the first l bits set to 1 and the rest to 0. */
void mask(BitVect *A, int l)
{
  invmask(A, l);
  InverseBV(A);
}

/* Sets A to the inverse mask with the first l bits set to 0 and the rest to 1. */
void invmask(BitVect *A, int l)
{
  AllOnes(A);
  BVRShiftSelf(A, l);
}

/* Sets all bits of A to 1. */
void AllOnes(BitVect *A)
{
  int i;

  for (i = 0; i < A->n; i++)
    A->vect[i] = 0xffffffffU;
}

/* Allocates a bit vector of l bits (rounded up to whole 32-bit words),
   initialised to zero. */
void AllocBV(BitVect *A, int l)
{
  int n;

  n = (l - 1) / WL + 1;
  A->vect = (uint32_t *) calloc((size_t)n, sizeof(uint32_t));
  A->n = n;
}

/* Frees the storage of bit vector A and resets its word count to zero. */
void FreeBV(BitVect *A)
{
  if (A->vect != NULL)
    free(A->vect);
  A->n = 0;
}

/* Allocates a matrix of nblines rows, each row containing t bit vectors of
   l bits. All entries are initialised to zero. */
void AllocMat(Matrix *m, int nblines, int l, int t)
{
  int i, j;

  m->lignes = (BitVect **) calloc((size_t)nblines, sizeof(BitVect *));
  for (i = 0; i < nblines; i++) {
    if (!(m->lignes[i] = (BitVect *) calloc((size_t)t, sizeof(BitVect)))) {
      printf("\n*** Insufficient memory in AllocMat()! nblines=%d ***\n", nblines);
      exit(1);
    }
    for (j = 0; j < t; j++)
      AllocBV(&(m->lignes[i][j]), l);
  }
  m->nblignes = nblines;
  m->t = t;
  m->l = l;
}

/* Frees all bit vectors and row arrays of matrix m, then resets its fields. */
void FreeMat(Matrix *m)
{
  int i, j;

  for (i = 0; i < m->nblignes; i++) {
    for (j = 0; j < m->t; j++)
      FreeBV(&(m->lignes[i][j]));
    free(m->lignes[i]);
  }
  free(m->lignes);
  m->nblignes = 0;
  m->l = 0;
  m->t = 0;
}

/* Copies the first nl rows and t bit-vector columns of matrix ms into m.
   m must already be allocated and large enough. */
void CopyMat(Matrix *m, Matrix *ms, int nl, int t)
{
  int i, j;

  if (m == NULL) {
    AllocMat(m, ms->nblignes, ms->l, ms->t);
  } else if ((ms->nblignes < nl) || (ms->t < t)) {
    printf("Error in CopyMat(): source matrix too small %d\n", ms->nblignes / ms->t);
    exit(1);
  } else if ((m->nblignes < nl) || (m->t < t)) {
    printf("Error in CopyMat(): destination matrix too small\n");
    exit(1);
  }
  for (i = 0; i < nl; i++)
    for (j = 0; j < t; j++)
      CopyBV(&(m->lignes[i][j]), &(ms->lignes[i][j]));
}

/* Swaps rows i and j of matrix m in place. */
void ExchangeVect(Matrix *m, int i, int j)
{
  BitVect *p;

  if (i != j) {
    p            = m->lignes[i];
    m->lignes[i] = m->lignes[j];
    m->lignes[j] = p;
  }
}

/* Computes the transpose of matrix M into T (T must already be allocated
   with compatible dimensions). Works for single-word bit vectors only. */
void TransposeMatrices(Matrix *T, Matrix *M)
{
  int s, l, m;

  if (!((T->nblignes == M->l) && (T->l == M->nblignes) && (T->t == M->t))) {
    printf("Error in TransposeMatrices!\n");
    exit(1);
  }
  for (s = 0; s < M->t; s++) {
    for (l = 0; l < M->l; l++)
      PutBVToZero(&T->lignes[l][s]);
    for (m = 0; m < M->nblignes; m++)
      for (l = 0; l < M->l; l++)
        if (M->lignes[m][s].vect[0] & (0x80000000U >> l))
          T->lignes[l][s].vect[0] |= (0x80000000U >> m);
  }
}

/* XORs row s into row r of matrix m for bit-vector columns min..max-1
   (i.e., m[r][j] ^= m[s][j] for j in [min, max)). */
void XorVect(Matrix *m, int r, int s, int min, int max)
{
  int j;

  for (j = min; j < max; j++)
    XORBVSelf(&(m->lignes[r][j]), &(m->lignes[s][j]));
}

/* Displays matrix m (kg rows, t groups of l bits each) to stdout. */
void DispMat(Matrix *m, int t, int l, int kg)
{
  int i, j;

  printf("    ");
  for (i = 0; i < t * l; i++)
    printf("%d", (i / 100) % 10);
  printf("\n    ");
  for (i = 0; i < t * l; i++)
    printf("%d", (i / 10) % 10);
  printf("\n    ");
  for (i = 0; i < t * l; i++)
    printf("%d", i % 10);
  printf("\n");
  for (i = 0; i < kg; i++) {
    printf("%3d[", i);
    for (j = 0; j < t; j++)
      DispBitVect(&(m->lignes[i][j]), l);
    printf("]\n");
  }
  printf("\n");
}

/* Displays the first l bits of bit vector A to stdout. */
void DispBitVect(BitVect *A, int l)
{
  int j;
  unsigned Un;

  Un = 1U;
  for (j = 0; j < l; j++)
    printf("%u", (A->vect[j / WL] >> (((WL * A->n) - j - 1) % WL)) & Un);
}
