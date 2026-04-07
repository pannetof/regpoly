/*
 * polylcg.c — Polynomial LCG generator for REGPOLY.
 *
 * Implements a linear feedback shift register based on a user-supplied
 * primitive polynomial over GF(2). At each step the state is shifted left
 * by one and XORed with the polynomial if the high bit was set.
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include "polylcg.h"
#include "regpoly.h"
#include "polynomials.h"

#define POLYLCGNAME  "Polynomial LCG"
#define POLYLCGSTATE (&(Gen->GenState))
#define POLYLCGPOL   (&(((paramPolyLCG *)(Gen->ParamGen))->Polynome))
#define POLYLCGMASK  (&(((paramPolyLCG *)(Gen->ParamGen))->gen_masque))

/* Returns the name string for the polynomial LCG generator type. */
char *PolyLCGName(void)
{
  return POLYLCGNAME;
}

/* Initializes the parameters of a polynomial LCG generator.
   Copies the characteristic polynomial, computes the state mask,
   and wires all function pointers for this type. */
void InitParamPolyLCG(Generateur *Gen, int k, BitVect *Polynome, int L)
{
  CopyBV(&(((paramPolyLCG *)(Gen->ParamGen))->Polynome), Polynome);
  mask(&(((paramPolyLCG *)(Gen->ParamGen))->gen_masque), k);
  Gen->k = k;
  Gen->L = L;
  Gen->Step = 1;
  Gen->smax = INT_MAX;
  Gen->DispGen = DispPolyLCG;
  Gen->DispName = PolyLCGName;
  Gen->InitGen = InitPolyLCG;
  Gen->Iteration = PolyLCG;
  Gen->CopyGen = CopyPolyLCG;
  Gen->AllocGen = AllocPolyLCG;
  Gen->PolyChar = CharPolyLCG;
  PutBVToZero(&(Gen->GenState));
}

/* Initializes the generator state from the bit vector init and copies
   the L low-order bits into retour. */
void InitPolyLCG(Generateur *Gen, BitVect *init, BitVect *retour)
{
  CopyBV(&(Gen->GenState), init);
  CopyBVPart(retour, init, Gen->L);
}

#define S 1

/* Advances the polynomial LCG by one step: shifts state left by 1,
   XORs with polynomial if the high bit was set, and writes L bits to retour. */
void PolyLCG(Generateur *Gen, BitVect *retour)
{
  int XOR = 0, j;

  for (j = 0; j < S; j++) {
    XOR = (ValBitBV(POLYLCGSTATE, 0) == 1);
    BVLShiftSelf(POLYLCGSTATE, 1);
    if (XOR)
      XORBVSelf(POLYLCGSTATE, POLYLCGPOL);
    ANDBVSelf(POLYLCGSTATE, POLYLCGMASK);
  }
  CopyBVPart(retour, POLYLCGSTATE, Gen->L);
}

/* Prints the characteristic polynomial in both exponent and hex notation. */
void DispPolyLCG(Generateur *Gen)
{
  int i;

  printf(" %d ", Gen->k);
  for (i = 0; i < Gen->k; i++)
    if (ValBitBV(POLYLCGPOL, i) == 1)
      printf("%d ", Gen->k - i - 1);
  printf("\nhexadecimal notation:\n ");
  for (i = 0; i <= (Gen->k - 1) / WL; i++)
    printf("%08x ", POLYLCGPOL->vect[i]);
  printf("\n");
}

/* Copies all fields of polynomial LCG Gen2 into Gen1,
   including deep copies of the polynomial and mask bit vectors. */
void CopyPolyLCG(Generateur *Gen1, Generateur *Gen2)
{
  CopyBV(&(((paramPolyLCG *)(Gen1->ParamGen))->Polynome),
         &(((paramPolyLCG *)(Gen2->ParamGen))->Polynome));
  CopyBV(&(((paramPolyLCG *)(Gen1->ParamGen))->gen_masque),
         &(((paramPolyLCG *)(Gen2->ParamGen))->gen_masque));

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

/* Allocates and returns a new polynomial LCG generator with state width k. */
Generateur *AllocPolyLCG(int k)
{
  Generateur *G;
  G = (Generateur *) malloc(sizeof(Generateur));
  AllocGen(G, k);
  G->ParamGen = (paramPolyLCG *) malloc(sizeof(paramPolyLCG));
  AllocBV(&(((paramPolyLCG *)G->ParamGen)->Polynome), k);
  AllocBV(&(((paramPolyLCG *)G->ParamGen)->gen_masque), k);
  return G;
}

/* Frees all memory associated with the polynomial LCG generator G. */
void FreePolyLCG(Generateur *G)
{
  FreeBV(&(((paramPolyLCG *)G->ParamGen)->gen_masque));
  FreeBV(&(((paramPolyLCG *)G->ParamGen)->Polynome));
  free(G->ParamGen);
  FreeGen(G);
}

/* Reads polynomial LCG parameters from file filepoly and adds each
   generator to component E. The polynomial is given as a list of
   exponents in decreasing order, terminated by 0.
   Exits if the file cannot be opened or allocation fails. */
void ReadDataPolyLCG(Component *E, char *filepoly, boolean same, int L)
{
  FILE *f;
  int j, t = 0;
  BitVect Poly;
  int Deg_Poly, nbpoly = 0;
  Generateur *p;

  f = fopen(filepoly, "r");
  if (f == NULL) {
    printf("File \"%s\" not found\n", filepoly);
    exit(1);
  }
  ReadLn(f);
  fscanf(f, "%d", &nbpoly);
  ReadLn(f);
  AllocGensInComponent(E, nbpoly, same);

  for (j = 0; j < nbpoly; j++) {
    fscanf(f, "%d", &t);
    Deg_Poly = t;
    AllocBV(&Poly, t);
    PutBVToZero(&Poly);
    while (t != 0) {
      fscanf(f, "%d", &t);
      PutBitBV(&Poly, Deg_Poly - t - 1, 1);
    }
    if ((p = AllocPolyLCG(Deg_Poly)) == NULL) {
      printf("Error in ReadDataPolyLCG()\n");
      exit(1);
    }
    InitParamPolyLCG(p, Deg_Poly, &Poly, L);
    FreeBV(&Poly);
    AddGenInComponent(E, p);
    FreePolyLCG(p);
  }
  fclose(f);
}

/* Generates nb random primitive polynomials of degree k and writes them
   to filename in the format expected by ReadDataPolyLCG. */
void ProduceParamsPolyLCG(char *filename, int nb, int k)
{
  FILE *f;
  int *Poly;
  int i, count, c, nbcoeff;

  f = fopen(filename, "w");
  count = 0;
  fprintf(f, "poly\n");
  fprintf(f, "%d\n", nb);
  Poly = (int *) malloc((k + 1) * sizeof(int));

  SetMRG32k3a((double)time(NULL), (double)time(NULL), (double)time(NULL),
              (double)time(NULL), (double)time(NULL), (double)time(NULL));
  while (count < nb) {
    for (i = 0; i < k; i++)
      Poly[i] = 0;
    Poly[0] = 1;
    Poly[k] = 1;
    nbcoeff = 2;

    for (i = 1; i < k; i++)
      if ((Poly[i] = MRG32k3a() % 2) == 1)
        nbcoeff++;
    if (!(nbcoeff % 2)) {
      c = MRG32k3a() % (k - 1) + 1;
      if (Poly[c])
        Poly[c] = 0;
      else
        Poly[c] = 1;
    }
    if (PrimitifPolynomial(Poly, k)) {
      for (i = k; i >= 0; i--)
        if (Poly[i])
          fprintf(f, "%d ", i);
      count++;
      fprintf(f, "\n");
    }
  }
  fclose(f);
  free(Poly);
}

/* Extracts the characteristic polynomial of the generator into coeff[]
   and the bit vector BVPoly. */
void CharPolyLCG(Generateur *Gen, int coeff[], BitVect *BVPoly)
{
  int j;

  for (j = 0; j <= Gen->k; j++)
    coeff[j] = 0;
  for (j = 0; j < Gen->k; j++)
    coeff[j] = ValBitBV(POLYLCGPOL, Gen->k - 1 - j);
  coeff[Gen->k] = 1;

  PutBVToZero(BVPoly);
  for (j = 0; j <= Gen->k; j++)
    PutBitBV(BVPoly, j, coeff[j]);
}
