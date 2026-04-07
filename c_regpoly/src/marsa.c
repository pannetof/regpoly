/*
 * marsa.c — Marsaglia three-shift XOR generator for REGPOLY.
 *
 * Implements a three-shift XOR generator as described by Marsaglia.
 * The state is a single w-bit word updated by three XOR-shift operations
 * with parameters a, b, and c. Positive values mean left shifts;
 * negative values mean right shifts.
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "marsa.h"
#include "polynomials.h"
#include "regpoly.h"

#define MARSANAME "Three-Shift Marsaglia Generator"

#define STATE (&(Gen->GenState))
#define WW    ((paramMarsa *)(Gen->ParamGen))->w
#define AA    ((paramMarsa *)(Gen->ParamGen))->a
#define BB    ((paramMarsa *)(Gen->ParamGen))->b
#define CC    ((paramMarsa *)(Gen->ParamGen))->c

/* Returns the name string for the Marsaglia generator type. */
char *MarsaName(void)
{
  return MARSANAME;
}

/* Initializes the parameters of a Marsaglia generator.
   Sets shift parameters a, b, c, output length L, and wires all
   function pointers for this generator type. */
void InitParamMarsa(Generateur *Gen, int w, int a, int b, int c, int L)
{
  ((paramMarsa *)(Gen->ParamGen))->w = w;
  ((paramMarsa *)(Gen->ParamGen))->a = a;
  ((paramMarsa *)(Gen->ParamGen))->b = b;
  ((paramMarsa *)(Gen->ParamGen))->c = c;

  Gen->k = w;
  Gen->L = L;
  Gen->Step = 1;
  Gen->smax = INT_MAX;
  Gen->DispGen = DispMarsa;
  Gen->DispName = MarsaName;
  Gen->InitGen = InitMarsa;
  Gen->Iteration = Marsa;
  Gen->CopyGen = CopyMarsa;
  Gen->AllocGen = AllocMarsa;
  Gen->PolyChar = polychar;
  PutBVToZero(&(Gen->GenState));
}

/* Initializes the generator state from the bit vector init and copies
   the L low-order bits into retour. */
void InitMarsa(Generateur *Gen, BitVect *init, BitVect *retour)
{
  CopyBV(&(Gen->GenState), init);
  CopyBVPart(retour, init, Gen->L);
}

/* Advances the Marsaglia generator by one step using three XOR-shifts
   with parameters a, b, c, and writes the L output bits to retour. */
void Marsa(Generateur *Gen, BitVect *retour)
{
  BitVect Temp;

  AllocBV(&Temp, WW);

  if (AA > 0)
    BVLShift(&Temp, STATE, AA);
  else
    BVRShift(&Temp, STATE, -AA);
  XORBV(STATE, &Temp, STATE);

  if (BB > 0)
    BVRShift(&Temp, STATE, BB);
  else
    BVLShift(&Temp, STATE, -BB);
  XORBV(STATE, &Temp, STATE);

  if (CC > 0)
    BVLShift(&Temp, STATE, CC);
  else
    BVRShift(&Temp, STATE, -CC);
  XORBV(STATE, &Temp, STATE);

  CopyBVPart(retour, STATE, Gen->L);
  FreeBV(&Temp);
}

/* Prints the shift parameters (w, a, b, c) of the generator. */
void DispMarsa(Generateur *Gen)
{
  printf("w=%d a=%d b=%d c=%d\n", WW, AA, BB, CC);
}

/* Copies all fields of Marsaglia generator Gen2 into Gen1. */
void CopyMarsa(Generateur *Gen1, Generateur *Gen2)
{
  ((paramMarsa *)(Gen1->ParamGen))->w = ((paramMarsa *)(Gen2->ParamGen))->w;
  ((paramMarsa *)(Gen1->ParamGen))->a = ((paramMarsa *)(Gen2->ParamGen))->a;
  ((paramMarsa *)(Gen1->ParamGen))->b = ((paramMarsa *)(Gen2->ParamGen))->b;
  ((paramMarsa *)(Gen1->ParamGen))->c = ((paramMarsa *)(Gen2->ParamGen))->c;

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

/* Allocates and returns a new Marsaglia generator with state width k. */
Generateur *AllocMarsa(int k)
{
  Generateur *G;
  G = (Generateur *) malloc(sizeof(Generateur));
  AllocGen(G, k);
  G->ParamGen = (paramMarsa *) malloc(sizeof(paramMarsa));
  return G;
}

/* Frees all memory associated with the Marsaglia generator G. */
void FreeMarsa(Generateur *G)
{
  free(G->ParamGen);
  FreeGen(G);
}

/* Reads Marsaglia generator parameters from file filepoly and adds all
   variants (all sign combinations of a, b, c) to component E.
   Exits if the file cannot be opened or allocation fails. */
void ReadDataMarsa(Component *E, char *filepoly, boolean same, int L)
{
  FILE *f;
  int j, a, b, c, w, nbgen = 0;
  Generateur *p;
  char gentype[25];

  f = fopen(filepoly, "r");
  if (f == NULL) {
    printf("File \"%s\" not found\n", filepoly);
    exit(1);
  }
  fscanf(f, "%s", gentype);
  ReadLn(f);
  fscanf(f, "%d %d", &nbgen, &w);
  ReadLn(f);
  AllocGensInComponent(E, nbgen * 8, same);

  for (j = 0; j < nbgen; j++) {
    fscanf(f, "%d %d %d", &a, &b, &c);
    ReadLn(f);

    if ((p = AllocMarsa(w)) == NULL) {
      printf("Error in ReadDataMarsa()\n");
      exit(1);
    }
    InitParamMarsa(p, w,  a,  b,  c, L); AddGenInComponent(E, p);
    InitParamMarsa(p, w,  c,  b,  a, L); AddGenInComponent(E, p);
    InitParamMarsa(p, w, -a, -b, -c, L); AddGenInComponent(E, p);
    InitParamMarsa(p, w, -c, -b, -a, L); AddGenInComponent(E, p);
    InitParamMarsa(p, w,  a, -c, -b, L); AddGenInComponent(E, p);
    InitParamMarsa(p, w,  c, -a, -b, L); AddGenInComponent(E, p);
    InitParamMarsa(p, w, -a,  c,  b, L); AddGenInComponent(E, p);
    InitParamMarsa(p, w, -c,  a,  b, L); AddGenInComponent(E, p);
    FreeMarsa(p);
  }
  fclose(f);
}
