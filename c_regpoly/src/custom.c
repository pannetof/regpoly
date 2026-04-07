/*
 * custom.c — Custom (user-defined) generator for REGPOLY.
 *
 * A custom generator is defined by a sequence of elementary linear
 * operations (LeftShift, RightShift, LeftShiftAND, RightShiftAND, XORIF)
 * whose results are XORed together to produce the next state.
 * Operations are specified in an input file and loaded via ReadDataCust.
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include "custom.h"
#include "polynomials.h"

#define CUSTNAME "Custom Generator"

#define MASK  (((paramCust *)(Gen->ParamGen))->Mask)
#define INV   (((paramCust *)(Gen->ParamGen))->inv)
#define TEMP1 (((paramCust *)(Gen->ParamGen))->temp1)
#define TEMP2 (((paramCust *)(Gen->ParamGen))->temp2)

/* Phase 1: extract the sub-vector starting at bit col with length portee. */
#define phase1(dummy, etat, col, portee, dummy2) \
  { invmask(&INV, col); ANDBVMask(&MASK, &INV, col + portee); ANDBV(etat, etat, &MASK); }

/* Phase 2: shift the sub-vector from column col to row, then mask to portee bits. */
#define phase2(dummy, etat, col, row, portee, dummy2) \
  { \
    invmask(&INV, row); \
    if (col < row) \
      BVRShiftSelf(etat, row - col); \
    else \
      BVLShiftSelf(etat, col - row); \
    ANDBVMask(&MASK, &INV, row + portee); \
    ANDBV(etat, etat, &MASK); \
  }

/* Returns the name string for the custom generator type. */
char *CustName(void)
{
  return CUSTNAME;
}

/* Initializes all parameters of a custom generator, copying the elementary
   operation descriptors from p, with degree k and output resolution L.
   Wires all function pointers for this type. */
void InitParamCust(Generateur *Gen, int k, paramCust *p, int L)
{
  int j;

  ((paramCust *)(Gen->ParamGen))->nbTrans = p->nbTrans;
  ((paramCust *)(Gen->ParamGen))->specCustTrans =
    (custTrans **) malloc(p->nbTrans * sizeof(custTrans *));

  for (j = 0; j < p->nbTrans; j++) {
    ((paramCust *)(Gen->ParamGen))->specCustTrans[j] =
      (custTrans *) malloc(sizeof(custTrans));
    strcpy(((paramCust *)(Gen->ParamGen))->specCustTrans[j]->nomtrans,
           p->specCustTrans[j]->nomtrans);
    ((paramCust *)(Gen->ParamGen))->specCustTrans[j]->row    = p->specCustTrans[j]->row;
    ((paramCust *)(Gen->ParamGen))->specCustTrans[j]->col    = p->specCustTrans[j]->col;
    ((paramCust *)(Gen->ParamGen))->specCustTrans[j]->portee = p->specCustTrans[j]->portee;
    ((paramCust *)(Gen->ParamGen))->specCustTrans[j]->nbbit  = p->specCustTrans[j]->nbbit;
    ((paramCust *)(Gen->ParamGen))->specCustTrans[j]->Trans  = p->specCustTrans[j]->Trans;
    if (((paramCust *)(Gen->ParamGen))->specCustTrans[j]->vecteur.n == 0)
      AllocBV(&(((paramCust *)(Gen->ParamGen))->specCustTrans[j]->vecteur), k);
    CopyBV(&(((paramCust *)(Gen->ParamGen))->specCustTrans[j]->vecteur),
           &(p->specCustTrans[j]->vecteur));
  }

  Gen->k = k;
  Gen->L = L;
  Gen->Step = 1;
  Gen->smax = INT_MAX;
  Gen->DispGen = DispCust;
  Gen->DispName = CustName;
  Gen->InitGen = InitCust;
  Gen->Iteration = Cust;
  Gen->CopyGen = CopyCust;
  Gen->AllocGen = AllocCust;
  Gen->PolyChar = polychar;
  PutBVToZero(&(Gen->GenState));
}

/* Initializes the generator state from init and copies the L low-order bits
   into retour. */
void InitCust(Generateur *Gen, BitVect *init, BitVect *retour)
{
  CopyBV(&(Gen->GenState), init);
  CopyBVPart(retour, init, Gen->L);
}

/* Performs one iteration of the custom generator: applies each elementary
   operation to the current state, XORs all results into the new state,
   and copies L bits to retour. */
void Cust(Generateur *Gen, BitVect *retour)
{
  int i;

  PutBVToZero(&TEMP2);
  for (i = 0; i < ((paramCust *)(Gen->ParamGen))->nbTrans; i++) {
    CopyBV(&TEMP1, &(Gen->GenState));
    ((paramCust *)(Gen->ParamGen))->specCustTrans[i]->Trans(
      Gen,
      &TEMP1,
      ((paramCust *)(Gen->ParamGen))->specCustTrans[i]->row,
      ((paramCust *)(Gen->ParamGen))->specCustTrans[i]->col,
      ((paramCust *)(Gen->ParamGen))->specCustTrans[i]->portee,
      ((paramCust *)(Gen->ParamGen))->specCustTrans[i]->nbbit,
      &(((paramCust *)(Gen->ParamGen))->specCustTrans[i]->vecteur));
    XORBVSelf(&TEMP2, &TEMP1);
  }
  CopyBV(&(Gen->GenState), &TEMP2);
  CopyBVPart(retour, &TEMP2, Gen->L);
}

/* Prints each elementary operation's name and parameters. */
void DispCust(Generateur *Gen)
{
  int i, j;
  BitVect G;

  AllocBV(&G, Gen->k);
  for (i = 0; i < ((paramCust *)(Gen->ParamGen))->nbTrans; i++) {
    printf("%s(row=%d, col=%d, portee=%d, nbbit=%d, BV)\n BV=",
           ((paramCust *)(Gen->ParamGen))->specCustTrans[i]->nomtrans,
           ((paramCust *)(Gen->ParamGen))->specCustTrans[i]->row,
           ((paramCust *)(Gen->ParamGen))->specCustTrans[i]->col,
           ((paramCust *)(Gen->ParamGen))->specCustTrans[i]->portee,
           ((paramCust *)(Gen->ParamGen))->specCustTrans[i]->nbbit);
    CopyBV(&G, &(((paramCust *)(Gen->ParamGen))->specCustTrans[i]->vecteur));
    BVLShiftSelf(&G, ((paramCust *)(Gen->ParamGen))->specCustTrans[i]->col);
    for (j = 0; j < G.n; j++)
      printf("%08lx ", G.vect[j]);
    printf("\n");
  }
  FreeBV(&G);
}

/* Copies all fields of custom generator Gen2 into Gen1,
   including deep copies of all elementary operation descriptors. */
void CopyCust(Generateur *Gen1, Generateur *Gen2)
{
  int nbtrans, j;

  nbtrans = ((paramCust *)(Gen1->ParamGen))->nbTrans =
            ((paramCust *)(Gen2->ParamGen))->nbTrans;
  ((paramCust *)(Gen1->ParamGen))->specCustTrans =
    (custTrans **) malloc(nbtrans * sizeof(custTrans *));

  CopyBV(&(((paramCust *)(Gen1->ParamGen))->inv),
         &(((paramCust *)(Gen2->ParamGen))->inv));
  for (j = 0; j < nbtrans; j++) {
    ((paramCust *)(Gen1->ParamGen))->specCustTrans[j] =
      (custTrans *) malloc(sizeof(custTrans));
    strcpy(((paramCust *)(Gen1->ParamGen))->specCustTrans[j]->nomtrans,
           ((paramCust *)(Gen2->ParamGen))->specCustTrans[j]->nomtrans);
    ((paramCust *)(Gen1->ParamGen))->specCustTrans[j]->row    = ((paramCust *)(Gen2->ParamGen))->specCustTrans[j]->row;
    ((paramCust *)(Gen1->ParamGen))->specCustTrans[j]->col    = ((paramCust *)(Gen2->ParamGen))->specCustTrans[j]->col;
    ((paramCust *)(Gen1->ParamGen))->specCustTrans[j]->portee = ((paramCust *)(Gen2->ParamGen))->specCustTrans[j]->portee;
    ((paramCust *)(Gen1->ParamGen))->specCustTrans[j]->nbbit  = ((paramCust *)(Gen2->ParamGen))->specCustTrans[j]->nbbit;
    ((paramCust *)(Gen1->ParamGen))->specCustTrans[j]->Trans  = ((paramCust *)(Gen2->ParamGen))->specCustTrans[j]->Trans;
    if (((paramCust *)(Gen1->ParamGen))->specCustTrans[j]->vecteur.n == 0)
      AllocBV(&(((paramCust *)(Gen1->ParamGen))->specCustTrans[j]->vecteur), Gen2->k);
    CopyBV(&(((paramCust *)(Gen1->ParamGen))->specCustTrans[j]->vecteur),
           &(((paramCust *)(Gen2->ParamGen))->specCustTrans[j]->vecteur));
  }
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

/* Allocates and returns a new custom generator with state width k,
   including all working bit vectors. */
Generateur *AllocCust(int k)
{
  Generateur *G;
  G = (Generateur *) malloc(sizeof(Generateur));
  AllocGen(G, k);
  G->ParamGen = (paramCust *) malloc(sizeof(paramCust));
  AllocBV(&(((paramCust *)(G->ParamGen))->Mask),  k);
  AllocBV(&(((paramCust *)(G->ParamGen))->temp2), k);
  AllocBV(&(((paramCust *)(G->ParamGen))->temp1), k);
  AllocBV(&(((paramCust *)(G->ParamGen))->inv),   k);
  return G;
}

/* Frees all memory associated with the custom generator G. */
void FreeCust(Generateur *G)
{
  FreeBV(&(((paramCust *)(G->ParamGen))->Mask));
  FreeBV(&(((paramCust *)(G->ParamGen))->inv));
  FreeBV(&(((paramCust *)(G->ParamGen))->temp1));
  FreeBV(&(((paramCust *)(G->ParamGen))->temp2));
  free(G->ParamGen);
  FreeGen(G);
}

/* Reads custom generator definitions from file filecust and adds each
   generator to component E. Each generator is described by a degree k,
   a count of operations, and one line per operation (LS, RS, LSA, RSA, XIF).
   Exits if the file cannot be opened or an unknown operation is encountered. */
void ReadDataCust(Component *E, char *filecust, boolean same, int L)
{
  FILE *f;
  int nbgen = 0, nbtrans;
  char trans[10];
  BitVect vecteur, Mask, inv;
  int i, j, h, row, col, portee, nbbit, k;
  Generateur *p;
  custTrans *pointeurcust;
  paramCust *pointeurparam;

  f = fopen(filecust, "r");
  if (f == NULL) {
    printf("File \"%s\" not found\n", filecust);
    exit(1);
  }
  while (!feof(f)) {
    fscanf(f, "%d %d ", &k, &nbtrans);
    ReadLn(f);
    nbgen++;
    for (i = 0; i < nbtrans; i++)
      ReadLn(f);
  }
  rewind(f);
  AllocGensInComponent(E, nbgen, same);
  pointeurparam = (paramCust *) malloc(sizeof(paramCust));

  for (j = 0; j < nbgen; j++) {
    fscanf(f, "%d %d", &k, &nbtrans);
    ReadLn(f);
    fflush(stdout);
    pointeurparam->nbTrans = nbtrans;
    pointeurparam->specCustTrans = (custTrans **) malloc(nbtrans * sizeof(custTrans *));

    AllocBV(&Mask, k);
    AllocBV(&vecteur, k);
    for (i = 0; i < nbtrans; i++) {
      PutBVToZero(&vecteur);
      fscanf(f, "%s", trans);

      if (strcmp(trans, "LS") == 0) {
        pointeurcust = (custTrans *) malloc(sizeof(custTrans));
        strcpy(pointeurcust->nomtrans, "LeftShift");
        fscanf(f, "%d %d %d %d", &row, &col, &portee, &nbbit);
        pointeurcust->Trans = LeftShift;
        ReadLn(f);
      } else if (strcmp(trans, "RS") == 0) {
        pointeurcust = (custTrans *) malloc(sizeof(custTrans));
        strcpy(pointeurcust->nomtrans, "RightShift");
        fscanf(f, "%d %d %d %d", &row, &col, &portee, &nbbit);
        pointeurcust->Trans = RightShift;
        ReadLn(f);
      } else if (strcmp(trans, "LSA") == 0) {
        pointeurcust = (custTrans *) malloc(sizeof(custTrans));
        strcpy(pointeurcust->nomtrans, "LeftShiftAND");
        fscanf(f, "%d %d %d %d", &row, &col, &portee, &nbbit);
        pointeurcust->Trans = LeftShiftAND;
        for (h = 0; h < (portee - 1) / WL + 1; h++)
          fscanf(f, "%lx", &(vecteur.vect[h]));
        ReadLn(f);
        BVRShift(&vecteur, &vecteur, col);
      } else if (strcmp(trans, "RSA") == 0) {
        pointeurcust = (custTrans *) malloc(sizeof(custTrans));
        strcpy(pointeurcust->nomtrans, "RightShiftAND");
        fscanf(f, "%d %d %d %d", &row, &col, &portee, &nbbit);
        pointeurcust->Trans = RightShiftAND;
        for (h = 0; h < (portee - 1) / WL + 1; h++)
          fscanf(f, "%lx", &(vecteur.vect[h]));
        ReadLn(f);
        BVRShift(&vecteur, &vecteur, col);
      } else if (strcmp(trans, "XIF") == 0) {
        pointeurcust = (custTrans *) malloc(sizeof(custTrans));
        strcpy(pointeurcust->nomtrans, "XORIF");
        fscanf(f, "%d %d %d %d", &row, &col, &portee, &nbbit);
        pointeurcust->Trans = XORIF;
        for (h = 0; h < (portee - 1) / WL + 1; h++)
          fscanf(f, "%lx", &(vecteur.vect[h]));
        ReadLn(f);
        BVRShift(&vecteur, &vecteur, col);
      } else {
        printf("Unknown transformation in custom.c: %s\n", trans);
        exit(1);
      }
      pointeurcust->row    = row;
      pointeurcust->col    = col;
      pointeurcust->portee = portee;
      pointeurcust->nbbit  = nbbit;
      AllocBV(&inv, k);
      invmask(&inv, col);
      ANDBVMask(&Mask, &inv, col + portee);
      ANDBV(&vecteur, &vecteur, &Mask);
      AllocBV(&(pointeurcust->vecteur), k);
      CopyBV(&(pointeurcust->vecteur), &vecteur);
      pointeurparam->specCustTrans[i] = pointeurcust;
      FreeBV(&inv);
    }

    FreeBV(&vecteur);
    FreeBV(&Mask);
    if ((p = AllocCust(k)) == NULL) {
      printf("Error in ReadDataCust()\n");
      exit(1);
    }
    InitParamCust(p, k, pointeurparam, L);
    AddGenInComponent(E, p);
    for (i = 0; i < nbtrans; i++)
      free(pointeurparam->specCustTrans[i]);
    free(pointeurparam->specCustTrans);
    FreeCust(p);
  }
  free(pointeurparam);
  fclose(f);
}

/* Extracts sub-vector [col, col+portee) from etat and shifts left by nbbit. */
void LeftShift(Generateur *Gen, BitVect *etat, int row, int col,
               int portee, int nbbit, BitVect *vecteur)
{
  phase1(Gen, etat, col, portee, vecteur->n * WL);
  BVLShift(etat, etat, nbbit);
  phase2(Gen, etat, col, row, portee, vecteur->n * WL);
}

/* Extracts sub-vector [col, col+portee) from etat and shifts right by nbbit. */
void RightShift(Generateur *Gen, BitVect *etat, int row, int col,
                int portee, int nbbit, BitVect *vecteur)
{
  phase1(Gen, etat, col, portee, vecteur->n * WL);
  BVRShift(etat, etat, nbbit);
  phase2(Gen, etat, col, row, portee, vecteur->n * WL);
}

/* Extracts sub-vector, shifts left by nbbit, then ANDs with vecteur mask. */
void LeftShiftAND(Generateur *Gen, BitVect *etat, int row, int col,
                  int portee, int nbbit, BitVect *vecteur)
{
  phase1(Gen, etat, col, portee, vecteur->n * WL);
  BVLShift(etat, etat, nbbit);
  ANDBV(etat, etat, vecteur);
  phase2(Gen, etat, col, row, portee, vecteur->n * WL);
}

/* Extracts sub-vector, shifts right by nbbit, then ANDs with vecteur mask. */
void RightShiftAND(Generateur *Gen, BitVect *etat, int row, int col,
                   int portee, int nbbit, BitVect *vecteur)
{
  phase1(Gen, etat, col, portee, vecteur->n * WL);
  BVRShift(etat, etat, nbbit);
  ANDBV(etat, etat, vecteur);
  phase2(Gen, etat, col, row, portee, vecteur->n * WL);
}

/* If bit col+nbbit of etat is set, replaces etat with vecteur (shifted to row);
   otherwise clears etat. */
void XORIF(Generateur *Gen, BitVect *etat, int row, int col,
           int portee, int nbbit, BitVect *vecteur)
{
  if (ValBitBV(etat, col + nbbit)) {
    phase1(Gen, etat, col, portee, vecteur->n * WL);
    CopyBV(etat, vecteur);
    phase2(Gen, etat, col, row, portee, vecteur->n * WL);
  } else {
    PutBVToZero(etat);
  }
}
