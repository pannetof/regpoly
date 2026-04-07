#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "GFSRWU.h"
#include "polynomials.h"

#define GFSRNAME "GFSR Pei-Chi Wu"

#define PP     ((paramGFSR *)(Gen->ParamGen))->p
#define QQ     ((paramGFSR *)(Gen->ParamGen))->q
#define RR     ((paramGFSR *)(Gen->ParamGen))->r
#define SS     ((paramGFSR *)(Gen->ParamGen))->s
#define DD     ((paramGFSR *)(Gen->ParamGen))->d
#define IX     ((paramGFSR *)(Gen->ParamGen))->ix
#define BUFFER(i) ((paramGFSR *)(Gen->ParamGen))->Buffer.vect[i]

/* Return the name string of this generator type. */
char *GFSRName(void) {
  return GFSRNAME;
}

/* Initialize all parameters of the GFSR generator pointed by Gen with the given
   recurrence parameters p, q, r, s, d and output resolution L. */
void InitParamGFSR(Generateur *Gen, int p, int q, int r, int s, int d, int L) {
  ((paramGFSR *)(Gen->ParamGen))->p = p;
  ((paramGFSR *)(Gen->ParamGen))->q = q;
  ((paramGFSR *)(Gen->ParamGen))->r = r;
  ((paramGFSR *)(Gen->ParamGen))->s = s;
  ((paramGFSR *)(Gen->ParamGen))->d = d;
  ((paramGFSR *)(Gen->ParamGen))->ix = 0;
  AllocBV(&(((paramGFSR *)(Gen->ParamGen))->Buffer), (p + q + d) * WL);

  Gen->k = p;
  Gen->L = L;
  Gen->Step = 1;
  Gen->smax = INT_MAX;
  Gen->DispGen = DispGFSR;
  Gen->DispName = GFSRName;
  Gen->InitGen = InitGFSR;
  Gen->Iteration = GFSR;
  Gen->CopyGen = CopyGFSR;
  Gen->AllocGen = AllocGFSR;
  Gen->PolyChar = polychar;
  PutBVToZero(&(Gen->GenState));
}

/* Initialize the state of the GFSR generator from the bit vector init, expanding
   it into the internal word buffer using the recurrence, then prime the buffer
   with one step.  Stores the last computed word in retour. */
void InitGFSR(Generateur *Gen, BitVect *init, BitVect *retour) {
  unsigned char *init_buf;
  int i, j, k;
  unsigned long tmp;

  init_buf = (unsigned char *) malloc(PP * sizeof(unsigned char));

  for (j = 0; j < PP; j++)
    init_buf[j] = ValBitBV(init, j);

  IX = 0;
  for (i = 0; i < PP; i++) { /* init P words */
    unsigned long w = 0;
    unsigned char b, t;
    int ixq;
    for (k = 0; k < DD; k++) {
      for (j = 0; j < WL; j++) {
        t = init_buf[IX];
        if ((ixq = IX - QQ) < 0)
          ixq += PP;
        t ^= init_buf[ixq];
        if ((ixq = IX - RR) < 0)
          ixq += PP;
        t ^= init_buf[ixq];
        if ((ixq = IX - SS) < 0)
          ixq += PP;
        t ^= init_buf[ixq];
        b = init_buf[IX] = t;
        if (++IX == PP)
          IX = 0;
        w = (w << 1) | (b & 1);
      }
      if (k == 0)
        BUFFER(i) = w;
    }
  }
  IX = PP;
  for (i = 0; i < QQ + DD; i++)
    BUFFER(PP + i) = BUFFER(i);

  for (i = 0; i < DD; i++)
    tmp = BUFFER(IX + i) = BUFFER(IX + i) ^ BUFFER(IX - QQ + i)
                           ^ BUFFER(IX - RR + i) ^ BUFFER(IX - SS + i);

  IX = IX + DD;
  if (IX >= PP + QQ) {
    for (i = IX; i < PP + QQ + DD; i++)
      BUFFER(i) = BUFFER(i) ^ BUFFER(i - QQ) ^ BUFFER(i - RR) ^ BUFFER(i - SS);
    IX = IX - PP;
    for (i = 0; i < IX; i++)
      BUFFER(i) = BUFFER(i + PP); /* copy the last IX numbers to head of buffer */
  }
  for (i = 0; i < PP; i++)
    printf("%08lx\n", BUFFER(i));
  exit(1);
  retour->vect[0] = tmp;
  printf("--------\n");
  printf("%08lx\n", tmp);
}

/* Perform one iteration of the GFSR recurrence, updating the circular word buffer
   and storing the new output word in retour. */
void GFSR(Generateur *Gen, BitVect *retour) {
  int i;
  unsigned long tmp;

  for (i = 0; i < DD; i++)
    tmp = BUFFER(IX + i) = BUFFER(IX + i) ^ BUFFER(IX - QQ + i)
                           ^ BUFFER(IX - RR + i) ^ BUFFER(IX - SS + i);

  IX = IX + DD;
  if (IX >= PP + QQ) {
    for (i = IX; i < PP + QQ + DD; i++)
      BUFFER(i) = BUFFER(i) ^ BUFFER(i - QQ) ^ BUFFER(i - RR) ^ BUFFER(i - SS);
    IX = IX - PP;
    for (i = 0; i < IX; i++)
      BUFFER(i) = BUFFER(i + PP); /* copy the last IX numbers to head of buffer */
  }
  retour->vect[0] = tmp;
  printf("%08lx\n", tmp);
}

/* Display the parameters of the GFSR generator and its characteristic polynomial. */
void DispGFSR(Generateur *Gen) {
  int *coeff;
  BitVect Poly;
  int sum, j;

  AllocBV(&Poly, Gen->k + 1);
  coeff = (int *) malloc((Gen->k + 1) * sizeof(int));

  printf("GFSR : p= %3d q= %3d r= %3d s= %3d   d=%d\n",
         ((paramGFSR *)(Gen->ParamGen))->p,
         ((paramGFSR *)(Gen->ParamGen))->q,
         ((paramGFSR *)(Gen->ParamGen))->r,
         ((paramGFSR *)(Gen->ParamGen))->s,
         ((paramGFSR *)(Gen->ParamGen))->d);
  polychar(Gen, coeff, &Poly);
  DispPolynomial(coeff, Gen->k);

  sum = 0;
  for (j = 0; j <= Gen->k; j++) {
    sum += coeff[j];
    if (coeff[j])
      printf("%d ", 521 - j);
  }
  printf("\nNumber of terms:%d\n", sum);
  free(coeff);
  FreeBV(&Poly);
}

/* Copy all parameters from Gen2 into Gen1, including the word buffer and all
   function pointers. */
void CopyGFSR(Generateur *Gen1, Generateur *Gen2) {
  int sum;

  sum =  ((paramGFSR *)(Gen1->ParamGen))->p = ((paramGFSR *)(Gen2->ParamGen))->p;
  sum += ((paramGFSR *)(Gen1->ParamGen))->q = ((paramGFSR *)(Gen2->ParamGen))->q;
         ((paramGFSR *)(Gen1->ParamGen))->r = ((paramGFSR *)(Gen2->ParamGen))->r;
         ((paramGFSR *)(Gen1->ParamGen))->s = ((paramGFSR *)(Gen2->ParamGen))->s;
  sum += ((paramGFSR *)(Gen1->ParamGen))->d = ((paramGFSR *)(Gen2->ParamGen))->d;
         ((paramGFSR *)(Gen1->ParamGen))->ix = ((paramGFSR *)(Gen2->ParamGen))->ix;

  AllocBV(&(((paramGFSR *)(Gen1->ParamGen))->Buffer), sum * WL);
  CopyBV(&(((paramGFSR *)(Gen1->ParamGen))->Buffer),
         &(((paramGFSR *)(Gen2->ParamGen))->Buffer));

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

/* Allocate a new Generateur and its GFSR parameter block for a generator of degree k. */
Generateur *AllocGFSR(int k) {
  Generateur *G;
  G = (Generateur *) malloc(sizeof(Generateur));
  AllocGen(G, k);
  G->ParamGen = (paramGFSR *) malloc(sizeof(paramGFSR));
  return G;
}

/* Free the parameter block and state vector of a GFSR generator. */
void FreeGFSR(Generateur *G) {
  free(G->ParamGen);
  // FreeBV(&(((paramGFSR *)(G->ParamGen))->Buffer));
  FreeGen(G);
}

/*
 * File format for ReadDataGFSR:
 *   <gentype string>
 *   <number of generators N>
 *   p q r s d       (one line per generator)
 *   ...
 */

/* Read GFSR generator parameters from file filepoly, fill Component E with
   the resulting generators, and mark whether this component is the same as
   the previous one (same flag) with output resolution L. */
void ReadDataGFSR(Component *E, char *filepoly, boolean same, int L) {
  FILE *f;
  int j, p, q, r, s, d;
  char gentype[25];
  int nbpoly = 0;
  Generateur *g;

  f = fopen(filepoly, "r");
  if (f == NULL) {
    printf("File \"%s\" not found\n", filepoly);
    exit(1);
  }
  fscanf(f, "%s", gentype);
  ReadLn(f);
  fscanf(f, "%d", &nbpoly);
  ReadLn(f);
  AllocGensInComponent(E, nbpoly, same);

  for (j = 0; j < nbpoly; j++) {
    fscanf(f, "%d %d %d %d %d", &p, &q, &r, &s, &d);

    if ((g = AllocGFSR(p)) == NULL) {
      printf("Error in ReadDataGFSR()\n");
      exit(1);
    }
    InitParamGFSR(g, p, q, r, s, d, L);
    AddGenInComponent(E, g);
    FreeGFSR(g);
  }
  fclose(f);
}
