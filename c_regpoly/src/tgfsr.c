/*
 * tgfsr.c — Twisted GFSR (TGFSR) generator for REGPOLY.
 *
 * Implements the Twisted Generalized Feedback Shift Register. The state
 * consists of r words of w bits; each step applies a twist matrix A to
 * two consecutive words then XORs in the word at distance m.
 * The characteristic polynomial is computed using the ZEN library.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include "tgfsr.h"
#include "polynomials.h"

#define TGFSRNAME "TGFSR"

#define TW   ((paramTGFSR *)(Gen->ParamGen))->w
#define TR   ((paramTGFSR *)(Gen->ParamGen))->r
#define TM   ((paramTGFSR *)(Gen->ParamGen))->m
#define TA   ((paramTGFSR *)(Gen->ParamGen))->a
#define ETAT Gen->GenState

/* Returns the name string for the TGFSR generator type. */
char *TGFSRName(void)
{
  return TGFSRNAME;
}

/* Prints the parameters (w, r, m, a) of the TGFSR generator. */
void DispTGFSR(Generateur *Gen)
{
  printf(" w= %3u  r= %3u  m= %3u  a= %08x\n",
         ((paramTGFSR *)(Gen->ParamGen))->w,
         ((paramTGFSR *)(Gen->ParamGen))->r,
         ((paramTGFSR *)(Gen->ParamGen))->m,
         ((paramTGFSR *)(Gen->ParamGen))->a.vect[0]);
}

/* Initializes all parameters of a TGFSR generator with word size w, register
   count r, middle distance m, twist matrix A, and output resolution L.
   Wires all function pointers for this type. */
void InitParamTGFSR(Generateur *Gen, int w, int r, int m, BitVect *A, int L)
{
  ((paramTGFSR *)(Gen->ParamGen))->w = w;
  ((paramTGFSR *)(Gen->ParamGen))->r = r;
  ((paramTGFSR *)(Gen->ParamGen))->m = m;
  CopyBV(&(((paramTGFSR *)(Gen->ParamGen))->a), A);

  Gen->k    = w * r;
  Gen->L    = L;
  Gen->Step = 1;
  Gen->smax = INT_MAX;
  Gen->DispGen   = DispTGFSR;
  Gen->DispName  = TGFSRName;
  Gen->InitGen   = InitTGFSR;
  Gen->Iteration = TGFSR;
  Gen->CopyGen   = CopyTGFSR;
  Gen->AllocGen  = AllocTGFSR;
  Gen->PolyChar  = CharTGFSR;
  PutBVToZero(&(Gen->GenState));
}

/* Copies all fields of TGFSR generator Gen2 into Gen1,
   including a deep copy of the twist matrix. */
void CopyTGFSR(Generateur *Gen1, Generateur *Gen2)
{
  ((paramTGFSR *)(Gen1->ParamGen))->w = ((paramTGFSR *)(Gen2->ParamGen))->w;
  ((paramTGFSR *)(Gen1->ParamGen))->r = ((paramTGFSR *)(Gen2->ParamGen))->r;
  ((paramTGFSR *)(Gen1->ParamGen))->m = ((paramTGFSR *)(Gen2->ParamGen))->m;
  AllocBV(&(((paramTGFSR *)(Gen1->ParamGen))->a), Gen2->k);
  CopyBV(&(((paramTGFSR *)(Gen1->ParamGen))->a),
         &(((paramTGFSR *)(Gen2->ParamGen))->a));

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

/* Allocates and returns a new TGFSR generator with state width k. */
Generateur *AllocTGFSR(int k)
{
  Generateur *G;
  G = (Generateur *) calloc(1, sizeof(Generateur));
  AllocGen(G, k);
  G->ParamGen = (paramTGFSR *) calloc(1, sizeof(paramTGFSR));
  AllocBV(&(((paramTGFSR *)(G->ParamGen))->a), k);
  return G;
}

/* Frees all memory associated with the TGFSR generator G. */
void FreeTGFSR(Generateur *G)
{
  FreeBV(&(((paramTGFSR *)(G->ParamGen))->a));
  free(G->ParamGen);
  FreeGen(G);
  free(G);
}

/* Initializes the TGFSR state from the bit vector init and copies
   the L low-order bits into retour. */
void InitTGFSR(Generateur *Gen, BitVect *init, BitVect *retour)
{
  CopyBV(&ETAT, init);
  CopyBVPart(retour, init, Gen->L);
}

/* Performs one TGFSR recurrence step and writes L output bits to retour. */
void TGFSR(Generateur *Gen, BitVect *RETOUR)
{
  BitVect Temp, Temp2;

  AllocBV(&Temp,  Gen->k);
  AllocBV(&Temp2, Gen->k);

  BVLShift(&Temp2, &ETAT, TW * (TR - TM - 1));
  BVLShift(&Temp,  &ETAT, TW * (TR - 1));
  ANDBVMask(&Temp2, &Temp2, TW);
  ANDBVMask(&Temp,  &Temp,  TW);

  if (ValBitBV(&Temp, TW - 1) == 0) {
    BVRShiftSelf(&Temp, 1);
    XORBV(&Temp, &Temp2, &Temp);
  } else {
    BVRShiftSelf(&Temp, 1);
    XOR2BV(&Temp, &Temp2, &Temp, &TA);
  }

  BVRShiftSelf(&ETAT, TW);
  ANDBVMask(&Temp, &Temp, TW);
  XORBV(&ETAT, &ETAT, &Temp);
  CopyBVPart(RETOUR, &ETAT, Gen->L);
  FreeBV(&Temp);
  FreeBV(&Temp2);
}

/* Reads TGFSR parameters from filename and adds each generator to component E.
   Exits if the file cannot be opened, w > 32, or allocation fails. */
void ReadDataTGFSR(Component *E, char *filename, boolean same, int L)
{
  int i, nbgen, m;
  int w, r;
  BitVect a_BV;
  FILE *f;
  Generateur *p;
  char gentype[25];

  f = fopen(filename, "r");
  if (f == NULL) {
    printf("File %s not found.\n", filename);
    exit(1);
  }
  fscanf(f, "%s", gentype);
  ReadLn(f);
  fscanf(f, "%d", &w);
  fscanf(f, "%d", &r);
  if (w > 32) {
    printf("Error: w must be <= 32 bits.\n");
    exit(1);
  }
  ReadLn(f);
  fscanf(f, "%d", &nbgen);
  AllocGensInComponent(E, nbgen, same);
  AllocBV(&a_BV, w * r);
  if ((p = AllocTGFSR(r * w)) == NULL) {
    printf("Error in ReadDataTGFSR()\n");
    exit(1);
  }
  ReadLn(f);
  for (i = 0; i < nbgen; i++) {
    fscanf(f, "%x", &(a_BV.vect[0]));
    fscanf(f, "%d",  &m);
    InitParamTGFSR(p, w, r, m, &a_BV, intmin(w, L));
    ReadLn(f);
    AddGenInComponent(E, p);
  }
  FreeBV(&a_BV);
  fclose(f);
  FreeTGFSR(p);
}

/* Generates nb random TGFSR parameter sets with word size w and register count r
   that yield primitive characteristic polynomials, and writes them to filename. */
void ProduceParamsTGFSR(char *filename, int nb, int w, int r)
{
  Generateur *p;
  BitVect a, BC, Dummy;
  FILE *f;
  int i, nbcoeff, c, count = 0, *Poly, m;

  SetMRG32k3a((double)time(NULL), (double)time(NULL), (double)time(NULL),
              (double)time(NULL), (double)time(NULL), (double)time(NULL));
  f = fopen(filename, "w");
  Poly = (int *) malloc((w * r + 1) * sizeof(int));
  fprintf(f, "tgfsr\n%d %d\n", w, r);
  fprintf(f, "%d\n", nb);

  while (count < nb) {
    for (i = 0; i <= w; i++)
      Poly[i] = 0;
    Poly[0] = 1;
    Poly[w] = 1;
    nbcoeff = 2;
    for (i = 1; i < w; i++)
      if ((Poly[i] = MRG32k3a() % 2) == 1)
        nbcoeff++;
    if (!(nbcoeff % 2)) {
      c = MRG32k3a() % (w - 1) + 1;
      if (Poly[c])
        Poly[c] = 0;
      else
        Poly[c] = 1;
    }
    if (IrreduciblePolynomial(Poly, w)) {
      AllocBV(&a, w * r);
      PutBVToZero(&a);
      AllocBV(&BC, w * r);
      BVCanonic(&BC, 0);
      for (i = 0; i < w; i++) {
        if (Poly[i])
          XORBV(&a, &a, &BC);
        BVRShiftSelf(&BC, 1);
      }
      FreeBV(&BC);
      if ((p = AllocTGFSR(w * r)) == NULL) {
        printf("Error in ProduceParamsTGFSR()\n");
        exit(1);
      }
      m = (MRG32k3a() % (r - 1)) + 1;
      InitParamTGFSR(p, w, r, m, &a, WL);
      AllocBV(&Dummy, (w * r) + 1);
      CharTGFSR(p, Poly, &Dummy);

      if (PrimitifPolynomial(Poly, w * r)) {
        fprintf(f, "%08x %d\n", a.vect[0], m);
        count++;
      }
      FreeBV(&a);
      FreeBV(&Dummy);
      FreeGen(p);
    }
  }
  fclose(f);
}

/* Computes the characteristic polynomial of the TGFSR using the ZEN library.
   Stores the result in coeff[] and BVPoly. */
void CharTGFSR(Generateur *Gen, int *coeff, BitVect *BVPoly)
{
  int j;
  ZENPoly tntm, res, temp;
  ZENRing F2;
  BigNum two;
  BigNumLength twol;
  ZENElt ONE, Coeff;

  ZBNReadFromString(&two, &twol, "2", 10);
  if (ZENBaseRingAlloc(F2, two, twol)) {
    printf("Error in CharTGFSR(): ZENBaseRingAlloc failed\n");
    exit(1);
  }
  ZENEltAlloc(ONE, F2);
  ZENEltSetToOne(ONE, F2);
  ZENEltAlloc(Coeff, F2);

  ZENPolyAlloc(temp, TR * TW, F2);
  ZENPolySetToZero(temp, F2);
  ZENPolyAlloc(res, TR * TW, F2);

  /* t^n + t^m */
  ZENPolyAlloc(tntm, TR * TW, F2);
  ZENPolySetToXi(tntm, TR, F2);
  ZENPolySetCoeff(tntm, TM, ONE, F2);

  ZENPolySetToZero(res, F2);
  for (j = 0; j < TW; j++)
    if (ValBitBV(&TA, j)) {
      PowPoly(&temp, &tntm, &F2, j);
      ZENPolyAdd(res, temp, F2);
    }
  PowPoly(&temp, &tntm, &F2, TW);
  ZENPolyAdd(res, temp, F2);

  for (j = 0; j <= TR * TW; j++) {
    ZENPolyGetCoeff(Coeff, res, j, F2);
    coeff[j] = ZENEltAreEqual(Coeff, ONE, F2) ? 1 : 0;
  }

  ZENEltFree(ONE, F2);
  ZBNF(two);
  ZENPolyFree(tntm, F2);
  ZENPolyFree(res,  F2);
  ZENPolyFree(temp, F2);
  if (ZENRingClose(F2)) {
    printf("Error in CharTGFSR(): ZENRingClose failed\n");
    exit(1);
  }
  PutBVToZero(BVPoly);
  for (j = 0; j <= Gen->k; j++)
    PutBitBV(BVPoly, j, coeff[j]);
}
