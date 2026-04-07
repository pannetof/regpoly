#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <regpoly.h>
#include <matsumoto.h>
#include <polynomials.h>

#define MatsumotoNAME "New generator"

#define GENTYPE       ((paramMatsumoto *)(Gen->ParamGen))->type
#define GEN_M         ((paramMatsumoto *)(Gen->ParamGen))->m
#define GEN_N         ((paramMatsumoto *)(Gen->ParamGen))->n
#define GEN_NBPARAMS  ((paramMatsumoto *)(Gen->ParamGen))->nbparamsint
#define GEN_PARAMS    ((paramMatsumoto *)(Gen->ParamGen))->paramsint
#define GEN_NBPARAMSU ((paramMatsumoto *)(Gen->ParamGen))->nbparamsunsigned
#define GEN_PARAMSU   ((paramMatsumoto *)(Gen->ParamGen))->paramsunsigned

/* Returns the name string for the Matsumoto generator type. */
char *MatsumotoName(void) {
  return MatsumotoNAME;
}

/* Initializes all parameters of a Matsumoto generator pointed to by Gen.
   Sets the generator type, n, m, integer and unsigned parameters, dimension k,
   output length L, and all function pointers for iteration, display, copy, etc. */
void InitParamMatsumoto(Generateur *Gen, int type, int n, int m,
                        int nbparamsint, int *paramsint,
                        int nbparamsunsigned, unsigned int *paramsunsigned, int L) {
  int j;
  GENTYPE = type;
  GEN_M = m;
  GEN_N = n;
  GEN_NBPARAMS  = nbparamsint;
  GEN_NBPARAMSU = nbparamsunsigned;
  if (GEN_PARAMS == NULL) {
    GEN_PARAMS  = (int *)          malloc(nbparamsint      * sizeof(int));
    GEN_PARAMSU = (unsigned int *) malloc(nbparamsunsigned * sizeof(unsigned int));
  }
  for (j = 0; j < nbparamsint; j++)
    GEN_PARAMS[j] = paramsint[j];
  for (j = 0; j < nbparamsunsigned; j++)
    GEN_PARAMSU[j] = paramsunsigned[j];

  Gen->k    = (n + 1) * 32;
  Gen->L    = L;
  Gen->Step = 1;
  Gen->smax = INT_MAX;
  Gen->DispGen  = DispMatsumoto;
  Gen->DispName = MatsumotoName;
  Gen->InitGen  = InitMatsumoto;
  if (type == 1)
    Gen->Iteration = Matsumoto1;
  else if (type == 2)
    Gen->Iteration = Matsumoto2;
  else if (type == 3) {
    Gen->Iteration = Matsumoto3;
    Gen->k = (n + 3) * 32;
  } else
    exit(1);
  Gen->CopyGen  = CopyMatsumoto;
  Gen->AllocGen = AllocMatsumoto;
  Gen->PolyChar = polychar;
  PutBVToZero(&(Gen->GenState));
}

/* Initializes the state of generator Gen from the bit vector init,
   and copies the first Gen->L bits back into retour. */
void InitMatsumoto(Generateur *Gen, BitVect *init, BitVect *retour) {
  CopyBV(&(Gen->GenState), init);
  CopyBVPart(retour, init, Gen->L);
}

#define VECT(i) ((Gen->GenState)).vect[i]
#define VECTU  VECT(GEN_N)
#define VECTU1 VECT(GEN_N + 1)
#define VECTU2 VECT(GEN_N + 2)
#define SHIFT(v, i) ((i > 0) ? (v >> i) : (v << (-(i))))

/* ***************************** */
/* Implementation  rand1         */
/* ***************************** */

/* Performs one iteration of Matsumoto generator type 1 and stores
   the first Gen->L output bits in retour. */
void Matsumoto1(Generateur *Gen, BitVect *retour) {
  unsigned int z0, z1, z2, z3, v0;
  z0 = VECTU ^ SHIFT(VECT(0), GEN_PARAMS[0]) ^ SHIFT(VECT(0), GEN_PARAMS[1]);
  z1 = z0 ^ VECT(GEN_M);
  z2 = z1 ^ SHIFT(z1, GEN_PARAMS[2]);
  z3 = VECT(0) ^ z2 ^ SHIFT(z2, GEN_PARAMS[3]);
  BVLShiftSelf(&(Gen->GenState), 32);
  VECT(GEN_N - 1) = z3;
  VECTU = z2;
  CopyBVPart(retour, &(Gen->GenState), Gen->L);
}

/* Performs one iteration of Matsumoto generator type 2 and stores
   the first Gen->L output bits in retour. */
void Matsumoto2(Generateur *Gen, BitVect *retour) {
  unsigned int z0, z1, z2, z3, z4;
  z0 = VECT(GEN_M);
  if (VECTU & 1U)
    z1 = (VECTU >> 1) ^ SHIFT(z0, GEN_PARAMS[0]) ^ SHIFT(z0, GEN_PARAMS[1]) ^ GEN_PARAMSU[0];
  else
    z1 = (VECTU >> 1) ^ SHIFT(z0, GEN_PARAMS[0]) ^ SHIFT(z0, GEN_PARAMS[1]);
  z2 = /*VECTU ^ */ VECT(0);
  z3 = z2 ^ SHIFT(z2, GEN_PARAMS[2]);
  z4 = z3 ^ SHIFT(z3, GEN_PARAMS[3]);
  BVLShiftSelf(&(Gen->GenState), 32);
  VECT(GEN_N - 1) = z1;
  VECTU = z4;
  CopyBVPart(retour, &(Gen->GenState), Gen->L);
}

/*
  u0 ^= x[kk] ^ (x[kk+M] >> p0);
  u1 ^= u0 ^ (u0 << p1);   // u1 = u1 + T1*u0
  u2 ^= u1 ^ (u1 >> p2);   // u2 = u2 + T2*u1
  x[kk] ^= u0 ^ u1 ^ u2;
*/

/* Performs one iteration of Matsumoto generator type 3 and stores
   the first Gen->L output bits in retour. */
void Matsumoto3(Generateur *Gen, BitVect *retour) {
  unsigned int z0, z1, z2, z3, v0;
  v0 = VECT(0);
  z0 = VECTU  ^ v0 ^ SHIFT(VECT(GEN_M), GEN_PARAMS[0]);
  z1 = VECTU1 ^ z0 ^ SHIFT(z0, GEN_PARAMS[1]);
  z2 = VECTU2 ^ z1 ^ SHIFT(z1, GEN_PARAMS[2]);
  BVLShiftSelf(&(Gen->GenState), 32);
  VECT(GEN_N - 1) = v0 ^ z0 ^ z1 ^ z2;
  VECTU           = z0;
  VECTU1          = z1;
  VECTU2          = z2;
  CopyBVPart(retour, &(Gen->GenState), Gen->L);
}

/* Displays the parameters (type, n, m, integer params, unsigned params) of
   the generator pointed to by Gen. */
void DispMatsumoto(Generateur *Gen) {
  int j;
  printf("ran%d  -> n = %d, m = %d, ", GENTYPE, GEN_N, GEN_M);
  for (j = 0; j < GEN_NBPARAMS; j++)
    printf("p%d = %d, ", j, GEN_PARAMS[j]);
  if (GENTYPE == 2) {
    for (j = 0; j < GEN_NBPARAMSU; j++)
      printf("q%d = %08x, ", j, GEN_PARAMSU[j]);
  }
  printf("\n");
}

/* Copies all parameters and state from generator Gen2 into generator Gen1,
   allocating parameter arrays in Gen1 if they have not yet been allocated. */
void CopyMatsumoto(Generateur *Gen1, Generateur *Gen2) {
  int j, nbparamsint, nbparamsunsigned;
  ((paramMatsumoto *)(Gen1->ParamGen))->type = ((paramMatsumoto *)(Gen2->ParamGen))->type;
  ((paramMatsumoto *)(Gen1->ParamGen))->m    = ((paramMatsumoto *)(Gen2->ParamGen))->m;
  ((paramMatsumoto *)(Gen1->ParamGen))->n    = ((paramMatsumoto *)(Gen2->ParamGen))->n;
  nbparamsint      = ((paramMatsumoto *)(Gen1->ParamGen))->nbparamsint      = ((paramMatsumoto *)(Gen2->ParamGen))->nbparamsint;
  nbparamsunsigned = ((paramMatsumoto *)(Gen1->ParamGen))->nbparamsunsigned = ((paramMatsumoto *)(Gen2->ParamGen))->nbparamsunsigned;
  if (((paramMatsumoto *)(Gen1->ParamGen))->paramsint == NULL) {
    ((paramMatsumoto *)(Gen1->ParamGen))->paramsint      = (int *)          malloc(nbparamsint      * sizeof(int));
    ((paramMatsumoto *)(Gen1->ParamGen))->paramsunsigned = (unsigned int *) malloc(nbparamsunsigned * sizeof(unsigned int));
  }
  for (j = 0; j < nbparamsint; j++)
    ((paramMatsumoto *)(Gen1->ParamGen))->paramsint[j] = ((paramMatsumoto *)(Gen2->ParamGen))->paramsint[j];
  for (j = 0; j < nbparamsunsigned; j++)
    ((paramMatsumoto *)(Gen1->ParamGen))->paramsunsigned[j] = ((paramMatsumoto *)(Gen2->ParamGen))->paramsunsigned[j];
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

/* Allocates and returns a new Matsumoto generator with a state of k bits.
   The parameter arrays are initialized to NULL and must be filled by InitParamMatsumoto(). */
Generateur *AllocMatsumoto(int k) {
  Generateur *G;
  G = (Generateur *) malloc(sizeof(Generateur));
  AllocGen(G, k);
  G->ParamGen = (paramMatsumoto *) malloc(sizeof(paramMatsumoto));
  ((paramMatsumoto *)(G->ParamGen))->paramsint = NULL;
  return G;
}

/* Frees all memory associated with the Matsumoto generator pointed to by Gen. */
void FreeMatsumoto(Generateur *Gen) {
  if (GEN_PARAMS != NULL) {
    free(GEN_PARAMS);
    free(GEN_PARAMSU);
  }
  free(Gen->ParamGen);
  FreeGen(Gen);
}

/* Reads a Matsumoto generator data file and fills the Component E with the
   generators it describes. The file format begins with "matsumoto" and a count,
   followed by one line per generator. */
void ReadDataMatsumoto(Component *E, char *filename, boolean same, int L) {
  FILE *f;
  int i, j, n, m, nbparamsint, nbparamsunsigned, nbgens = 0, type;
  int *paramsint, *paramsunsigned;
  Generateur *p;

  f = fopen(filename, "r");
  if (f == NULL) {
    printf("File \"%s\" not found\n", filename);
    exit(1);
  }
  ReadLn(f);
  fscanf(f, "%d", &nbgens);
  ReadLn(f);
  AllocGensInComponent(E, nbgens, same);

  for (j = 0; j < nbgens; j++) {
    fscanf(f, "%d %d %d %d", &type, &n, &m, &nbparamsint);
    paramsint = (int *) malloc(nbparamsint * sizeof(int));
    for (i = 0; i < nbparamsint; i++)
      fscanf(f, "%d ", &paramsint[i]);

    fscanf(f, "%d", &nbparamsunsigned);
    paramsunsigned = (unsigned int *) malloc(nbparamsunsigned * sizeof(unsigned int));
    for (i = 0; i < nbparamsunsigned; i++)
      fscanf(f, "%x ", &paramsunsigned[i]);
    if (type == 3) {
      if ((p = AllocMatsumoto(32 * (n + 3))) == NULL) {
        printf("Error in ProduceMatsumoto()\n");
        exit(1);
      }
    } else {
      if ((p = AllocMatsumoto(32 * (n + 1))) == NULL) {
        printf("Error in ProduceMatsumoto()\n");
        exit(1);
      }
    }
    InitParamMatsumoto(p, type, n, m, nbparamsint, paramsint, nbparamsunsigned, paramsunsigned, L);
    free(paramsint);
    free(paramsunsigned);
    AddGenInComponent(E, p);
    FreeMatsumoto(p);
  }
  fclose(f);
}

/* Generates nb sets of Matsumoto parameters of the given type and writes them to
   filename in a format readable by ReadDataMatsumoto(). Parameters are chosen at
   random using MRG32k3a. If Mexp != 0, only parameters whose characteristic
   polynomial contains that Mersenne exponent are accepted; otherwise, only
   parameters with a primitive characteristic polynomial are accepted. */
void ProduceParamsMatsumoto(char *filename, int type, int nb, int n, int Mexp) {
  FILE *f;
  int count = 0, j, m, nbparamsint, nbparamsunsigned, good, k;
  int *paramsint, *paramsunsigned;
  BitVect Dummy;
  int *Poly;
  Generateur *p;

  if (type == 1) {
    nbparamsint      = 4;
    nbparamsunsigned = 0;
    k = (n + 1) * 32;
  } else if (type == 2) {
    nbparamsint      = 4;
    nbparamsunsigned = 1;
    k = (n + 1) * 32;
  } else if (type == 3) {
    nbparamsint      = 3;
    nbparamsunsigned = 0;
    k = (n + 3) * 32;
  } else
    exit(1);

  paramsint      = (int *)          malloc(nbparamsint      * sizeof(int));
  paramsunsigned = (unsigned int *) malloc(nbparamsunsigned * sizeof(unsigned int));

  f = fopen(filename, "w");
  count = 0;
  fprintf(f, "matsumoto\n");
  fprintf(f, "%d\n", nb);

  SetMRG32k3a((double)time(NULL), (double)time(NULL), (double)time(NULL),
              (double)time(NULL), (double)time(NULL), (double)time(NULL));

  if ((p = AllocMatsumoto(k)) == NULL) {
    printf("Error in ProduceMatsumoto()\n");
    exit(1);
  }
  Poly = (int *) calloc(k + 1, sizeof(int));
  AllocBV(&Dummy, k + 1);

  if (1) {
    while (count < nb) {
      int OK;
      m  = MRG32k3a() % n;
      OK = 0;
      while (!OK) {
        for (j = 0; j < nbparamsint; j++)
          paramsint[j] = (MRG32k3a() % 63) - 31;
        if (type == 3) {
          if (gcd(abs(paramsint[0]), abs(paramsint[1])) == 1 &&
              gcd(abs(paramsint[1]), abs(paramsint[2])) == 1)
            OK = 1;
        } else {
          OK = 1;
        }
      }
      for (j = 0; j < nbparamsunsigned; j++)
        paramsunsigned[j] = MRG32k3a() | 0x80000000U;
      InitParamMatsumoto(p, type, n, m, nbparamsint, paramsint,
                         nbparamsunsigned, paramsunsigned, 32);
      polychar(p, Poly, &Dummy);
      good = 0;
      if (Mexp != 0) {
        if (ContainsMersenneExponent(k, Poly, Mexp))
          good = 1;
      } else {
        if (PrimitifPolynomial(Poly, k))
          good = 1;
      }
      if (good) {
        fprintf(f, "%d %d %d %d ", type, n, m, nbparamsint);
        for (j = 0; j < nbparamsint; j++)
          fprintf(f, "%d ", paramsint[j]);
        fprintf(f, "%d ", nbparamsunsigned);
        for (j = 0; j < nbparamsunsigned; j++)
          fprintf(f, "%08x ", paramsunsigned[j]);
        count++;
        fprintf(f, "\n");
        fflush(f);
      }
    }
  }
  fclose(f);
  free(Poly);
  FreeBV(&Dummy);
  free(paramsint);
  free(paramsunsigned);
}
