/*
 * tmsnets.c — Digital net generator (base 2) for REGPOLY.
 *
 * Implements a (t,m,s)-net in base 2 using precomputed generating matrices.
 * The generator enumerates points of the net by reading successive rows of
 * the transposed generating matrices. Each component corresponds to one
 * coordinate of the s-dimensional net.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tmsnets.h"
#include "regpoly.h"

#define TMSNAME  "Digital net in base 2"
#define GMAT     (((paramTMS *)(Gen->ParamGen))->M)
#define ITER     (((paramTMS *)(Gen->ParamGen))->i)
#define CURRBC   (((paramTMS *)(Gen->ParamGen))->currentbc)

/* Returns the name string for the digital net generator type. */
char *TMSName(void)
{
  return TMSNAME;
}

/* Initializes the parameters of a digital net generator.
   Allocates or reallocates the generating matrix, copies m rows from M,
   and wires all function pointers. */
void InitParamTMS(Generateur *Gen, Matrix *m, int mmax, int smax,
                  int Lmax, int L, char name[50])
{
  printf("$");
  if (GMAT.nblignes == 0)
    AllocMat(&GMAT, mmax, L, smax);
  else {
    FreeMat(&GMAT);
    AllocMat(&GMAT, mmax, L, smax);
  }
  ((paramTMS *)(Gen->ParamGen))->mmax = mmax;
  ((paramTMS *)(Gen->ParamGen))->smax = smax;
  ((paramTMS *)(Gen->ParamGen))->Lmax = Lmax;
  CopyMat(&GMAT, m, mmax, smax);
  strcpy(((paramTMS *)(Gen->ParamGen))->name, name);

  Gen->smax = smax;
  Gen->k = mmax;
  Gen->L = intmin(L, Lmax);
  Gen->Step = 1;
  Gen->DispGen = DispTMS;
  Gen->DispName = TMSName;
  Gen->InitGen = InitTMS;
  Gen->Iteration = TMS;
  Gen->CopyGen = CopyTMS;
  Gen->AllocGen = AllocTMS;
  Gen->PolyChar = NULL;
  PutBVToZero(&(Gen->GenState));
}

/* Initializes the generator state from a canonical basis vector init.
   Determines which basis vector (currentbc) was given and copies the
   corresponding row of the generating matrix to retour. */
void InitTMS(Generateur *Gen, BitVect *init, BitVect *retour)
{
  int i, j;

  CopyBV(&(Gen->GenState), init);
  ITER = i = j = 0;
  while (!init->vect[i]) i++;
  while (!(init->vect[i] & (0x80000000UL >> j))) j++;
  CURRBC = i * WL + j;
  CopyBVPart(retour, &(GMAT.lignes[CURRBC][ITER]), Gen->L);
}

/* Advances to the next point of the net and writes L output bits to retour. */
void TMS(Generateur *Gen, BitVect *retour)
{
  CopyBVPart(retour, &(GMAT.lignes[CURRBC][++ITER]), Gen->L);
}

/* Prints the name, m, and s parameters of the digital net generator. */
void DispTMS(Generateur *Gen)
{
  printf("%s    m=%3d  s=%3d\n",
         ((paramTMS *)(Gen->ParamGen))->name,
         ((paramTMS *)(Gen->ParamGen))->mmax,
         ((paramTMS *)(Gen->ParamGen))->smax);
}

/* Copies all fields of digital net generator Gen2 into Gen1,
   including a deep copy of the generating matrix. */
void CopyTMS(Generateur *Gen1, Generateur *Gen2)
{
  if (((paramTMS *)(Gen2->ParamGen))->M.nblignes != 0) {
    AllocMat(&(((paramTMS *)(Gen1->ParamGen))->M),
             ((paramTMS *)(Gen2->ParamGen))->M.nblignes,
             ((paramTMS *)(Gen2->ParamGen))->M.l,
             ((paramTMS *)(Gen2->ParamGen))->M.t);
    CopyMat(&(((paramTMS *)(Gen1->ParamGen))->M),
            &(((paramTMS *)(Gen2->ParamGen))->M),
            ((paramTMS *)(Gen2->ParamGen))->M.nblignes,
            ((paramTMS *)(Gen2->ParamGen))->M.t);
  }
  ((paramTMS *)(Gen1->ParamGen))->mmax = ((paramTMS *)(Gen2->ParamGen))->mmax;
  ((paramTMS *)(Gen1->ParamGen))->smax = ((paramTMS *)(Gen2->ParamGen))->smax;
  strcpy(((paramTMS *)(Gen1->ParamGen))->name,
         ((paramTMS *)(Gen2->ParamGen))->name);

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

/* Allocates and returns a new digital net generator with state width k. */
Generateur *AllocTMS(int k)
{
  Generateur *G;
  G = (Generateur *) malloc(sizeof(Generateur));
  AllocGen(G, k);
  G->ParamGen = (paramTMS *) malloc(sizeof(paramTMS));
  ((paramTMS *)G->ParamGen)->M.nblignes = 0;
  return G;
}

/* Frees all memory associated with the digital net generator G,
   including the generating matrix. */
void FreeTMS(Generateur *G)
{
  printf("T");
  FreeMat(&(((paramTMS *)G->ParamGen)->M));
  free(G->ParamGen);
  FreeGen(G);
}

/* Reads digital net generating matrices from file filename and adds each
   net as a generator in component E. The matrices are transposed before
   storage. Exits if the file cannot be opened, L > 32, or allocation fails. */
void ReadDataTMS(Component *E, char *filename, boolean same, int L)
{
  FILE *f;
  int j, l, s, m, Lmax, mmax, smax = 0;
  int nbnets = 0;
  char name[50];
  ulong mthrow;
  Generateur *p;
  Matrix M, T;

  f = fopen(filename, "r");
  if (f == NULL) {
    printf("File \"%s\" not found\n", filename);
    exit(1);
  }
  if (L > 32) {
    printf("Value of L is too large\n");
    exit(1);
  }
  ReadLn(f);
  fscanf(f, "%d", &nbnets);
  ReadLn(f);
  AllocGensInComponent(E, nbnets, same);

  for (j = 0; j < nbnets; j++) {
    fscanf(f, "%d %d %d", &mmax, &smax, &Lmax);
    ReadLn(f);

    AllocMat(&M, Lmax, mmax, smax);
    for (s = 0; s < smax; s++) {
      for (l = 0; l < Lmax; l++) {
        fscanf(f, "%lu", &mthrow);
        M.lignes[l][s].vect[0] = (mthrow << (WL - mmax));
      }
      ReadLn(f);
    }
    AllocMat(&T, mmax, Lmax, smax);  /* transpose of the generating matrices */
    TransposeMatrices(&T, &M);
    FreeMat(&M);

    if ((p = AllocTMS(mmax)) == NULL) {
      printf("Error in ReadDataTMS()\n");
      exit(1);
    }
    sprintf(name, "%s #%d", filename, j);
    InitParamTMS(p, &T, mmax, smax, Lmax, L, name);
    AddGenInComponent(E, p);
    FreeTMS(p);
    FreeMat(&T);
  }
  fclose(f);
}
