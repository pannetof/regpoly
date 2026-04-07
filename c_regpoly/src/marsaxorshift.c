#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <marsaxorshift.h>
#include <polynomials.h>
#include <regpoly.h>
#include <time.h>

#define MARSAXORSHIFTNAME "Marsaglia's Xor-shift RNGs"

static inline ulong ShiftL(ulong *b, int s);
static inline ulong ShiftR(ulong *b, int s);

void MarsaXorshiftType1(Generateur *Gen, BitVect *retour);
void MarsaXorshiftType2(Generateur *Gen, BitVect *retour);
void MarsaXorshiftType3(Generateur *Gen, BitVect *retour);

static void init(Generateur *Gen, int type, int w, int r, int m, int a, int b, int c,
                 int *paramstype2p, int *paramstype2q,
                 int paramstype3N, int *paramstype3No, int *paramstype3Shift,
                 int nbxorshift, int nbmi, int *mi, int *mi_nb, int *ai, int L);

/* Returns the name string for this generator type. */
char* MarsaXorshiftName(void) {
  return MARSAXORSHIFTNAME;
}

/* Initializes the parameters for a Type I (single-word xorshift) generator. */
void InitParamMarsaXorshiftTypeI(Generateur *Gen, int w, int a, int b, int c, int L) {
  init(Gen, TYPE1, w, 0, 0, a, b, c, NULL, NULL, 0, NULL, NULL, 0, 0, NULL, NULL, NULL, L);
}

/* Initializes the parameters for a Type II (two-component xorshift) generator. */
void InitParamMarsaXorshiftTypeII(Generateur *Gen, int type, int w, int r, int m, int *paramstype2p, int *paramstype2q, int L) {
  init(Gen, type, w, r, m, 0, 0, 0, paramstype2p, paramstype2q, 0, NULL, NULL, 0, 0, NULL, NULL, NULL, L);
}

/* Initializes the parameters for a Type III (multi-tap xorshift) generator. */
void InitParamMarsaXorshiftTypeIII(Generateur *Gen, int w, int r, int paramstype3N, int *paramstype3No, int *paramstype3Shift, int L) {
  init(Gen, TYPE3, w, r, 0, 0, 0, 0, NULL, NULL, paramstype3N, paramstype3No, paramstype3Shift, 0, 0, NULL, NULL, NULL, L);
}

/* Initializes the parameters for a Type IV (combined two-component, 4-shift) generator. */
void InitParamMarsaXorshiftTypeIV(Generateur *Gen, int a, int b, int c, int d, int w, int r, int m, int L) {
  int paramstype4p[2];
  int paramstype4q[2];
  paramstype4p[0] = a;
  paramstype4p[1] = b;
  paramstype4q[0] = c;
  paramstype4q[1] = d;
  init(Gen, TYPE4, w, r, m, 0, 0, 0, paramstype4p, paramstype4q, 0, NULL, NULL, 0, 0, NULL, NULL, NULL, L);
}

/* Initializes the parameters for a general xorshift generator with arbitrary tap structure. */
void InitParamMarsaXorshiftTypeGeneral(Generateur *Gen, int w, int r, int nbxorshift, int nbmi, int *mi, int *mi_nb, int *ai, int L) {
  init(Gen, TYPEGENERAL, w, r, 0, 0, 0, 0, 0, 0, 0, NULL, NULL, nbxorshift, nbmi, mi, mi_nb, ai, L);
}

/* Common initializer: sets all fields of the generator parameter structure and function pointers. */
static void init(Generateur *Gen, int type, int w, int r, int m, int a, int b, int c,
                 int *paramstype2p, int *paramstype2q,
                 int paramstype3N, int *paramstype3No, int *paramstype3Shift,
                 int nbxorshift, int nbmi, int *mi, int *mi_nb, int *ai, int L) {
  int j;

  if (type != TYPE1 && w != 32) {
    printf("Error in InitParamMarsaXorshift(). For type=%d, only w=32 is supported so far (you put w=%d)\n", type, w);
    exit(1);
  }

  ((paramMarsaXorshift*)(Gen->ParamGen))->type = type;
  ((paramMarsaXorshift*)(Gen->ParamGen))->w = w;
  ((paramMarsaXorshift*)(Gen->ParamGen))->m = m;
  ((paramMarsaXorshift*)(Gen->ParamGen))->a = a;
  ((paramMarsaXorshift*)(Gen->ParamGen))->b = b;
  ((paramMarsaXorshift*)(Gen->ParamGen))->c = c;
  ((paramMarsaXorshift*)(Gen->ParamGen))->nbmi = nbmi;
  ((paramMarsaXorshift*)(Gen->ParamGen))->nbxorshift = nbxorshift;

  if (type == TYPE1)
    ((paramMarsaXorshift*)(Gen->ParamGen))->r = 1;
  else
    ((paramMarsaXorshift*)(Gen->ParamGen))->r = r;

  if (type >= TYPE21 && type <= TYPE25) {
    if (((paramMarsaXorshift*)(Gen->ParamGen))->paramstype2q == NULL) {
      ((paramMarsaXorshift*)(Gen->ParamGen))->paramstype2q = (int *) calloc(3, sizeof(int));
      ((paramMarsaXorshift*)(Gen->ParamGen))->paramstype2p = (int *) calloc(3, sizeof(int));
    }
    for (j = 0; j < 3; j++) {
      ((paramMarsaXorshift*)(Gen->ParamGen))->paramstype2q[j] = paramstype2q[j];
      ((paramMarsaXorshift*)(Gen->ParamGen))->paramstype2p[j] = paramstype2p[j];
    }
  }

  if (type == TYPE3) {
    ((paramMarsaXorshift*)(Gen->ParamGen))->paramstype3N = paramstype3N;
    if (((paramMarsaXorshift*)(Gen->ParamGen))->paramstype3No == NULL) {
      ((paramMarsaXorshift*)(Gen->ParamGen))->paramstype3No = (int *) calloc(paramstype3N, sizeof(int));
      ((paramMarsaXorshift*)(Gen->ParamGen))->paramstype3Shift = (int *) calloc(paramstype3N, sizeof(int));
    }
    for (j = 0; j < paramstype3N; j++) {
      ((paramMarsaXorshift*)(Gen->ParamGen))->paramstype3No[j] = paramstype3No[j];
      ((paramMarsaXorshift*)(Gen->ParamGen))->paramstype3Shift[j] = paramstype3Shift[j];
    }
  }

  if (type == TYPE4) {
    if (((paramMarsaXorshift*)(Gen->ParamGen))->paramstype2q == NULL) {
      ((paramMarsaXorshift*)(Gen->ParamGen))->paramstype2q = (int *) calloc(2, sizeof(int));
      ((paramMarsaXorshift*)(Gen->ParamGen))->paramstype2p = (int *) calloc(2, sizeof(int));
    }
    for (j = 0; j < 2; j++) {
      ((paramMarsaXorshift*)(Gen->ParamGen))->paramstype2q[j] = paramstype2q[j];
      ((paramMarsaXorshift*)(Gen->ParamGen))->paramstype2p[j] = paramstype2p[j];
    }
  }

  if (type == TYPEGENERAL) {
    if (((paramMarsaXorshift*)(Gen->ParamGen))->mi == NULL) {
      ((paramMarsaXorshift*)(Gen->ParamGen))->mi = (int *) calloc(nbmi, sizeof(int));
      ((paramMarsaXorshift*)(Gen->ParamGen))->mi_nb = (int *) calloc(nbmi, sizeof(int));
      ((paramMarsaXorshift*)(Gen->ParamGen))->ai = (int *) calloc(nbxorshift, sizeof(int));
    }
    for (j = 0; j < nbmi; j++) {
      ((paramMarsaXorshift*)(Gen->ParamGen))->mi[j] = mi[j];
      ((paramMarsaXorshift*)(Gen->ParamGen))->mi_nb[j] = mi_nb[j];
    }
    for (j = 0; j < nbxorshift; j++)
      ((paramMarsaXorshift*)(Gen->ParamGen))->ai[j] = ai[j];
  }

  Gen->k = w * r;
  Gen->L = L;
  Gen->Step = 1;
  Gen->smax = INT_MAX;
  Gen->DispGen = DispMarsaXorshift;
  Gen->DispName = MarsaXorshiftName;
  Gen->InitGen = InitMarsaXorshift;
  if (type == TYPE1) {
    Gen->Iteration = MarsaXorshiftType1;
  } else if (type >= TYPE21 && type <= TYPE25) {
    Gen->Iteration = MarsaXorshiftType2;
  } else if (type == TYPE3) {
    Gen->Iteration = MarsaXorshiftType3;
  } else if (type == TYPE4) {
    Gen->Iteration = MarsaXorshiftType4;
  } else if (type == TYPEGENERAL) {
    Gen->Iteration = MarsaXorshiftTypeGeneral;
  } else {
    printf("Error in InitParamMarsaXorshift(). type = %d not recognized.\n", type);
    exit(1);
  }
  Gen->CopyGen = CopyMarsaXorshift;
  Gen->AllocGen = AllocMarsaXorshift;
  Gen->PolyChar = polychar;
  PutBVToZero(&(Gen->GenState));
}

/* Copies the initial state init into the generator and the first L bits into retour. */
void InitMarsaXorshift(Generateur *Gen, BitVect *init, BitVect *retour) {
  CopyBV(&(Gen->GenState), init);
  CopyBVPart(retour, init, Gen->L);
}

#define STATE   (&(Gen->GenState))
#define WW      ((paramMarsaXorshift*)(Gen->ParamGen))->w
#define AA      ((paramMarsaXorshift*)(Gen->ParamGen))->a
#define BB      ((paramMarsaXorshift*)(Gen->ParamGen))->b
#define CC      ((paramMarsaXorshift*)(Gen->ParamGen))->c
#define RR      ((paramMarsaXorshift*)(Gen->ParamGen))->r
#define MM      ((paramMarsaXorshift*)(Gen->ParamGen))->m

#define TYPE    ((paramMarsaXorshift*)(Gen->ParamGen))->type
#define P(i)    ((paramMarsaXorshift*)(Gen->ParamGen))->paramstype2p[i]
#define Q(i)    ((paramMarsaXorshift*)(Gen->ParamGen))->paramstype2q[i]
#define N3      ((paramMarsaXorshift*)(Gen->ParamGen))->paramstype3N
#define NO3     ((paramMarsaXorshift*)(Gen->ParamGen))->paramstype3No
#define SHIFT3  ((paramMarsaXorshift*)(Gen->ParamGen))->paramstype3Shift
/* Type general macros */
#define NBMI        ((paramMarsaXorshift*)(Gen->ParamGen))->nbmi
#define MI(i)       ((paramMarsaXorshift*)(Gen->ParamGen))->mi[i]
#define MI_NB(i)    ((paramMarsaXorshift*)(Gen->ParamGen))->mi_nb[i]
#define NBXORSHIFT  ((paramMarsaXorshift*)(Gen->ParamGen))->nbxorshift
#define AI(i)       ((paramMarsaXorshift*)(Gen->ParamGen))->ai[i]

/* Returns *b >> s if s > 0, or *b << (-s) if s < 0. */
static ulong ShiftR(ulong *b, int s) {
  if (s < 0)
    return *b << (-s);
  else
    return *b >> s;
}

/* Returns *b << s if s > 0, or *b >> (-s) if s < 0. */
static ulong ShiftL(ulong *b, int s) {
  if (s < 0)
    return *b >> (-s);
  else
    return *b << s;
}

/* Performs one iteration of the general xorshift generator with arbitrary tap/shift structure. */
void MarsaXorshiftTypeGeneral(Generateur *Gen, BitVect *retour) {
  ulong t = 0UL;
  int j;
  ulong temp;
  int i, aicounter = 0;
  for (i = 0; i < NBMI; i++) {
    temp = STATE->vect[MI(i) - 1];
    for (j = 0; j < MI_NB(i); j++)
      temp = temp ^ ShiftR(&temp, AI(aicounter++));
    t ^= temp;
  }
  BVRShiftSelf(STATE, WW);
  STATE->vect[0] = t;
  CopyBVPart(retour, STATE, Gen->L);
}

/* Performs one iteration of the Type I (single-word triple xorshift) generator. */
void MarsaXorshiftType1(Generateur *Gen, BitVect *retour) {
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

/* Performs one iteration of the Type II (two-component lagged xorshift) generator. */
void MarsaXorshiftType2(Generateur *Gen, BitVect *retour) {
  ulong x, y, j;
  y = STATE->vect[MM - 1];
  x = STATE->vect[RR - 1];
  BVRShiftSelf(STATE, WW);
  for (j = 0; j < 3; j++) {
    if (P(j) != 0)
      x ^= ShiftR(&x, P(j));
  }
  for (j = 0; j < 3; j++) {
    if (Q(j) != 0)
      y ^= ShiftR(&y, Q(j));
  }
  STATE->vect[0] = x ^ y;
  CopyBVPart(retour, STATE, Gen->L);
}

/* Performs one iteration of the Type IV (two-component, 2-shift each) generator. */
void MarsaXorshiftType4(Generateur *Gen, BitVect *retour) {
  ulong x, y, j;
  y = STATE->vect[MM - 1];
  x = STATE->vect[RR - 1];
  BVRShiftSelf(STATE, WW);
  for (j = 0; j < 2; j++)
    x ^= ShiftR(&x, P(j));
  for (j = 0; j < 2; j++)
    y ^= ShiftR(&y, Q(j));
  STATE->vect[0] = x ^ y;
  CopyBVPart(retour, STATE, Gen->L);
}

/* Performs one iteration of the Type III (multi-tap single xorshift sum) generator. */
void MarsaXorshiftType3(Generateur *Gen, BitVect *retour) {
  ulong t = 0UL;
  int j;
  for (j = 0; j < N3; j++)
    t ^= STATE->vect[NO3[j] - 1] ^ ShiftR(&(STATE->vect[NO3[j] - 1]), SHIFT3[j]);
  BVRShiftSelf(STATE, WW);
  STATE->vect[0] = t;
  CopyBVPart(retour, STATE, Gen->L);
}

/* Displays the parameters of the generator pointed by Gen. */
void DispMarsaXorshift(Generateur *Gen) {
  int nbxors, i, j, aicounter = 0;
  if (TYPE == TYPEGENERAL) {
    printf("Type = general w=%d  r=%d\n", WW, RR);
    nbxors = 0;
    for (i = 0; i < NBMI; i++) {
      printf("m_%d = %d ", i + 1, MI(i));
      for (j = MI_NB(i) - 1; j >= 0; j--) {
        if (AI(aicounter) < 0)
          printf("(I+L^{%d})", -AI(aicounter));
        else
          printf("(I+R^{%d})", AI(aicounter));
        aicounter++;
      }
      printf("  ");
    }
    printf("\n");
  } else if (TYPE == TYPE3) {
    printf("Type=%d  w=%d  r=%d ", TYPE, WW, RR);
    for (j = 0; j < N3; j++) {
      if (SHIFT3[j] < 0)
        printf(", m_%d = %d  (I+L^{%d})", j + 1, NO3[j], -SHIFT3[j]);
      else if (SHIFT3[j] > 0)
        printf(", m_%d = %d  (I+R^{%d})", j + 1, NO3[j], SHIFT3[j]);
    }
    printf("\n");
  } else if ((TYPE >= TYPE21 && TYPE <= TYPE25) || TYPE == TYPE4) {
    if (TYPE == TYPE4)
      nbxors = 2;
    else
      nbxors = 3;
    printf("Type=%d   w=%d  r=%d  m=%d  ", TYPE, WW, RR, MM);
    printf("B = ");
    for (j = nbxors - 1; j >= 0; j--) {
      if (P(j) > 0)
        printf("(I+D^%d)", P(j));
      else if (P(j) < 0)
        printf("(I+G^%d)", -P(j));
    }
    printf("   A = ");
    for (j = nbxors - 1; j >= 0; j--) {
      if (Q(j) > 0)
        printf("(I+D^%d)", Q(j));
      else if (Q(j) < 0)
        printf("(I+G^%d)", -Q(j));
    }
    printf("\n");
  } else {
    printf("Type=%d  w=%d  r=%d  a=%d  b=%d  c=%d\n", TYPE, WW, RR, AA, BB, CC);
  }
}

/* Copies all parameters and state from Gen2 into Gen1 (both must be MarsaXorshift). */
void CopyMarsaXorshift(Generateur *Gen1, Generateur *Gen2) {
  int j, r;
  ((paramMarsaXorshift *)(Gen1->ParamGen))->w = ((paramMarsaXorshift *)(Gen2->ParamGen))->w;
  r = ((paramMarsaXorshift *)(Gen1->ParamGen))->r = ((paramMarsaXorshift *)(Gen2->ParamGen))->r;
  ((paramMarsaXorshift *)(Gen1->ParamGen))->m = ((paramMarsaXorshift *)(Gen2->ParamGen))->m;
  ((paramMarsaXorshift *)(Gen1->ParamGen))->a = ((paramMarsaXorshift *)(Gen2->ParamGen))->a;
  ((paramMarsaXorshift *)(Gen1->ParamGen))->b = ((paramMarsaXorshift *)(Gen2->ParamGen))->b;
  ((paramMarsaXorshift *)(Gen1->ParamGen))->c = ((paramMarsaXorshift *)(Gen2->ParamGen))->c;
  ((paramMarsaXorshift *)(Gen1->ParamGen))->type = ((paramMarsaXorshift *)(Gen2->ParamGen))->type;

  if (((paramMarsaXorshift *)(Gen2->ParamGen))->type == TYPEGENERAL) {
    if (((paramMarsaXorshift *)(Gen1->ParamGen))->mi != NULL) {
      free(((paramMarsaXorshift *)(Gen1->ParamGen))->mi);
      free(((paramMarsaXorshift *)(Gen1->ParamGen))->mi_nb);
      free(((paramMarsaXorshift *)(Gen1->ParamGen))->ai);
    }
    ((paramMarsaXorshift *)(Gen1->ParamGen))->nbmi = ((paramMarsaXorshift *)(Gen2->ParamGen))->nbmi;
    ((paramMarsaXorshift *)(Gen1->ParamGen))->nbxorshift = ((paramMarsaXorshift *)(Gen2->ParamGen))->nbxorshift;
    ((paramMarsaXorshift *)(Gen1->ParamGen))->mi = (int*) malloc(((paramMarsaXorshift *)(Gen2->ParamGen))->nbmi * sizeof(int));
    ((paramMarsaXorshift *)(Gen1->ParamGen))->mi_nb = (int*) malloc(((paramMarsaXorshift *)(Gen2->ParamGen))->nbmi * sizeof(int));
    ((paramMarsaXorshift *)(Gen1->ParamGen))->ai = (int*) malloc(((paramMarsaXorshift *)(Gen2->ParamGen))->nbxorshift * sizeof(int));

    for (j = 0; j < ((paramMarsaXorshift *)(Gen2->ParamGen))->nbmi; j++) {
      ((paramMarsaXorshift *)(Gen1->ParamGen))->mi[j] = ((paramMarsaXorshift *)(Gen2->ParamGen))->mi[j];
      ((paramMarsaXorshift *)(Gen1->ParamGen))->mi_nb[j] = ((paramMarsaXorshift *)(Gen2->ParamGen))->mi_nb[j];
    }
    for (j = 0; j < ((paramMarsaXorshift *)(Gen2->ParamGen))->nbxorshift; j++)
      ((paramMarsaXorshift *)(Gen1->ParamGen))->ai[j] = ((paramMarsaXorshift *)(Gen2->ParamGen))->ai[j];
  }

  if (((paramMarsaXorshift *)(Gen2->ParamGen))->type == TYPE3) {
    if (((paramMarsaXorshift *)(Gen1->ParamGen))->paramstype3N != 0) {
      free(((paramMarsaXorshift *)(Gen1->ParamGen))->paramstype3No);
      free(((paramMarsaXorshift *)(Gen1->ParamGen))->paramstype3Shift);
    }
    ((paramMarsaXorshift *)(Gen1->ParamGen))->paramstype3N = ((paramMarsaXorshift *)(Gen2->ParamGen))->paramstype3N;
    ((paramMarsaXorshift *)(Gen1->ParamGen))->paramstype3No = (int*) calloc(((paramMarsaXorshift *)(Gen1->ParamGen))->paramstype3N, sizeof(int));
    ((paramMarsaXorshift *)(Gen1->ParamGen))->paramstype3Shift = (int*) calloc(((paramMarsaXorshift *)(Gen1->ParamGen))->paramstype3N, sizeof(int));
    for (j = 0; j < ((paramMarsaXorshift *)(Gen2->ParamGen))->paramstype3N; j++) {
      ((paramMarsaXorshift *)(Gen1->ParamGen))->paramstype3No[j] = ((paramMarsaXorshift *)(Gen2->ParamGen))->paramstype3No[j];
      ((paramMarsaXorshift *)(Gen1->ParamGen))->paramstype3Shift[j] = ((paramMarsaXorshift *)(Gen2->ParamGen))->paramstype3Shift[j];
    }
  }

  if (((paramMarsaXorshift *)(Gen2->ParamGen))->type >= TYPE21 && ((paramMarsaXorshift *)(Gen2->ParamGen))->type <= TYPE25) {
    if (((paramMarsaXorshift *)(Gen1->ParamGen))->paramstype2p != NULL) {
      free(((paramMarsaXorshift *)(Gen1->ParamGen))->paramstype2p);
      free(((paramMarsaXorshift *)(Gen1->ParamGen))->paramstype2q);
    }
    ((paramMarsaXorshift *)(Gen1->ParamGen))->paramstype2p = (int*) calloc(3, sizeof(int));
    ((paramMarsaXorshift *)(Gen1->ParamGen))->paramstype2q = (int*) calloc(3, sizeof(int));
    for (j = 0; j < 3; j++) {
      ((paramMarsaXorshift *)(Gen1->ParamGen))->paramstype2p[j] = ((paramMarsaXorshift *)(Gen2->ParamGen))->paramstype2p[j];
      ((paramMarsaXorshift *)(Gen1->ParamGen))->paramstype2q[j] = ((paramMarsaXorshift *)(Gen2->ParamGen))->paramstype2q[j];
    }
  }

  if (((paramMarsaXorshift *)(Gen2->ParamGen))->type == TYPE4) {
    if (((paramMarsaXorshift *)(Gen1->ParamGen))->paramstype2p != NULL) {
      free(((paramMarsaXorshift *)(Gen1->ParamGen))->paramstype2p);
      free(((paramMarsaXorshift *)(Gen1->ParamGen))->paramstype2q);
    }
    ((paramMarsaXorshift *)(Gen1->ParamGen))->paramstype2p = (int*) calloc(2, sizeof(int));
    ((paramMarsaXorshift *)(Gen1->ParamGen))->paramstype2q = (int*) calloc(2, sizeof(int));
    for (j = 0; j < 2; j++) {
      ((paramMarsaXorshift *)(Gen1->ParamGen))->paramstype2p[j] = ((paramMarsaXorshift *)(Gen2->ParamGen))->paramstype2p[j];
      ((paramMarsaXorshift *)(Gen1->ParamGen))->paramstype2q[j] = ((paramMarsaXorshift *)(Gen2->ParamGen))->paramstype2q[j];
    }
  }

  Gen1->k = Gen2->k;
  Gen1->L = Gen2->L;
  Gen1->smax = Gen2->smax;
  Gen1->Step = Gen2->Step;
  Gen1->DispGen = Gen2->DispGen;
  Gen1->DispName = Gen2->DispName;
  Gen1->InitGen = Gen2->InitGen;
  Gen1->Iteration = Gen2->Iteration;
  Gen1->PolyChar = Gen2->PolyChar;
  Gen1->CopyGen = Gen2->CopyGen;
  Gen1->AllocGen = Gen2->AllocGen;
  CopyBV(&(Gen1->GenState), &(Gen2->GenState));
}

/* Allocates and returns a new MarsaXorshift generator with state size k bits. */
Generateur* AllocMarsaXorshift(int k) {
  Generateur *G;
  G = (Generateur *) calloc(1, sizeof(Generateur));
  AllocGen(G, k);
  G->ParamGen = (paramMarsaXorshift *) calloc(1, sizeof(paramMarsaXorshift));
  return G;
}

/* Frees all memory associated with a MarsaXorshift generator. */
void FreeMarsaXorshift(Generateur *Gen) {
  if (TYPE == TYPEGENERAL) {
    if (((paramMarsaXorshift *)(Gen->ParamGen))->mi != NULL) {
      free(((paramMarsaXorshift *)(Gen->ParamGen))->mi);
      free(((paramMarsaXorshift *)(Gen->ParamGen))->mi_nb);
      free(((paramMarsaXorshift *)(Gen->ParamGen))->ai);
    }
  }

  if (TYPE == TYPE3)
    if (((paramMarsaXorshift *)(Gen->ParamGen))->paramstype3N != 0) {
      free(((paramMarsaXorshift *)(Gen->ParamGen))->paramstype3No);
      free(((paramMarsaXorshift *)(Gen->ParamGen))->paramstype3Shift);
    }

  if ((((paramMarsaXorshift *)(Gen->ParamGen))->type >= TYPE21 && ((paramMarsaXorshift *)(Gen->ParamGen))->type <= TYPE25) ||
      (((paramMarsaXorshift *)(Gen->ParamGen))->type == TYPE4)) {
    if (((paramMarsaXorshift *)(Gen->ParamGen))->paramstype2p != NULL) {
      free(((paramMarsaXorshift *)(Gen->ParamGen))->paramstype2p);
      free(((paramMarsaXorshift *)(Gen->ParamGen))->paramstype2q);
    }
  }
  free(Gen->ParamGen);
  FreeGen(Gen);
  free(Gen);
}

/* Reads generator parameters from a file and populates the Component E. */
void ReadDataMarsaXorshift(Component *E, char *filepoly, boolean same, int L) {
  FILE *f;
  int j, a, b, c, d, w, r, m, i, type, nbgen = 0, ii;
  Generateur *p;
  char gentype[25];
  int type2p[3], type2q[3];
  int originaltype2p[3], originaltype2q[3];
  int numberGenToAllocate;
  int nbmi, nbxorshift, aicounter;
  int *mi, *mi_nb, *ai;
  int type3N, *type3No, *type3Shift;

  f = fopen(filepoly, "r");
  if (f == NULL) {
    printf("File \"%s\" not found\n", filepoly);
    exit(1);
  }
  fscanf(f, "%s", gentype);
  ReadLn(f);
  fscanf(f, "%d %d %d", &nbgen, &r, &w);
  ReadLn(f);
  fscanf(f, "%d", &numberGenToAllocate);
  ReadLn(f);
  AllocGensInComponent(E, numberGenToAllocate, same);

  for (j = 0; j < nbgen; j++) {
    fscanf(f, "%d", &type);

    if (type == TYPEGENERAL) {
      fscanf(f, "%d", &nbmi);
      fscanf(f, "%d", &nbxorshift);
      aicounter = 0;
      mi = (int*) malloc(nbmi * sizeof(int));
      mi_nb = (int*) malloc(nbmi * sizeof(int));
      ai = (int*) malloc(nbxorshift * sizeof(int));
      for (i = 0; i < nbmi; i++) {
        fscanf(f, "%d", &mi[i]);
        fscanf(f, "%d", &mi_nb[i]);
        for (ii = 0; ii < mi_nb[i]; ii++)
          fscanf(f, "%d", &ai[aicounter++]);
      }
    } else if (type == TYPE3) {
      fscanf(f, "%d", &type3N);
      type3No = (int*) calloc(type3N, sizeof(int));
      type3Shift = (int*) calloc(type3N, sizeof(int));
      for (i = 0; i < type3N; i++)
        fscanf(f, "%d %d", &type3No[i], &type3Shift[i]);
    } else if (type >= TYPE21 && type <= TYPE25) {
      fscanf(f, "%d %d %d %d %d %d %d", &m,
             &originaltype2p[0], &originaltype2p[1], &originaltype2p[2],
             &originaltype2q[0], &originaltype2q[1], &originaltype2q[2]);
      for (i = 0; i < 3; i++) {
        type2p[i] = originaltype2p[i];
        type2q[i] = originaltype2q[i];
      }
    } else if (type == TYPE4) {
      fscanf(f, "%d %d %d %d %d %d", &r, &m, &a, &b, &c, &d); /* Note: r is re-read */
    } else if (type == TYPE1) {
      fscanf(f, "%d %d %d", &a, &b, &c);
    } else {
      printf("Type=%d not recognized!\n", type);
      exit(1);
    }
    ReadLn(f);

    if ((p = AllocMarsaXorshift(r * w)) == NULL) {
      printf("Error in ReadDataMarsaXorshift()\n");
      exit(1);
    }

    if (type == TYPE1) {
      InitParamMarsaXorshiftTypeI(p, w, a, b, c, L);   // X1
      AddGenInComponent(E, p);
      InitParamMarsaXorshiftTypeI(p, w, c, b, a, L);   // X2
      AddGenInComponent(E, p);
      InitParamMarsaXorshiftTypeI(p, w, -a, -b, -c, L); // X3
      AddGenInComponent(E, p);
      InitParamMarsaXorshiftTypeI(p, w, a, -c, -b, L);  // X5
      AddGenInComponent(E, p);
    } else if (type == TYPE21) { // (A1,B1) --> (A4,B4)
      InitParamMarsaXorshiftTypeII(p, type, w, r, m, type2p, type2q, L);
      AddGenInComponent(E, p);
      type2p[0] = -originaltype2p[0]; type2p[1] = -originaltype2p[1]; type2q[0] = -originaltype2q[0];
      InitParamMarsaXorshiftTypeII(p, type, w, r, m, type2p, type2q, L);
      AddGenInComponent(E, p);
      type2p[0] = -originaltype2p[1]; type2p[1] = -originaltype2p[0]; type2q[0] = -originaltype2q[0];
      InitParamMarsaXorshiftTypeII(p, type, w, r, m, type2p, type2q, L);
      AddGenInComponent(E, p);
    } else if (type == TYPE22) { // (A5,B5) --> (A8,B8)
      InitParamMarsaXorshiftTypeII(p, type, w, r, m, type2p, type2q, L);
      AddGenInComponent(E, p);
      type2q[0] = -originaltype2q[0]; type2q[1] = -originaltype2q[1]; type2p[0] = -originaltype2p[0];
      InitParamMarsaXorshiftTypeII(p, type, w, r, m, type2p, type2q, L);
      AddGenInComponent(E, p);
      type2q[0] = -originaltype2q[1]; type2q[1] = -originaltype2q[0]; type2p[0] = -originaltype2p[0];
      InitParamMarsaXorshiftTypeII(p, type, w, r, m, type2p, type2q, L);
      AddGenInComponent(E, p);
    } else if (type == TYPE23) { // (A9,B9) --> (A10,B10)
      InitParamMarsaXorshiftTypeII(p, type, w, r, m, type2p, type2q, L);
      AddGenInComponent(E, p);
      type2p[0] = -originaltype2p[0]; type2p[1] = -originaltype2p[1]; type2q[0] = -originaltype2q[0];
      InitParamMarsaXorshiftTypeII(p, type, w, r, m, type2p, type2q, L);
      AddGenInComponent(E, p);
    } else if (type == TYPE24) { // (A11,B11) --> (A12,B12)
      InitParamMarsaXorshiftTypeII(p, type, w, r, m, type2p, type2q, L);
      AddGenInComponent(E, p);
      type2q[0] = -originaltype2q[0]; type2q[1] = -originaltype2q[1]; type2p[0] = -originaltype2p[0];
      InitParamMarsaXorshiftTypeII(p, type, w, r, m, type2p, type2q, L);
      AddGenInComponent(E, p);
    } else if (type == TYPE25) { // (A13,B13) --> (A20,B20)
      InitParamMarsaXorshiftTypeII(p, type, w, r, m, type2p, type2q, L);
      AddGenInComponent(E, p);
      type2p[0] = originaltype2p[2]; type2p[1] = originaltype2p[1]; type2p[2] = originaltype2p[0];
      InitParamMarsaXorshiftTypeII(p, type, w, r, m, type2p, type2q, L);
      AddGenInComponent(E, p);
      type2p[0] = -originaltype2p[0]; type2p[1] = -originaltype2p[1]; type2p[2] = -originaltype2p[2];
      InitParamMarsaXorshiftTypeII(p, type, w, r, m, type2p, type2q, L);
      AddGenInComponent(E, p);
      type2p[0] = originaltype2p[0]; type2p[1] = originaltype2p[2]; type2p[2] = originaltype2p[1];
      InitParamMarsaXorshiftTypeII(p, type, w, r, m, type2p, type2q, L);
      AddGenInComponent(E, p);
    } else if (type == TYPE3) {
      InitParamMarsaXorshiftTypeIII(p, w, r, type3N, type3No, type3Shift, L);
      AddGenInComponent(E, p);
      free(type3No);
      free(type3Shift);
    } else if (type == TYPE4) {
      InitParamMarsaXorshiftTypeIV(p, a, b, c, d, w, r, m, L);
      AddGenInComponent(E, p);
    } else if (type == TYPEGENERAL) {
      InitParamMarsaXorshiftTypeGeneral(p, w, r, nbxorshift, nbmi, mi, mi_nb, ai, L);
      AddGenInComponent(E, p);
      free(mi);
      free(mi_nb);
      free(ai);
    } else {
      printf("Type not recognized!\n");
      exit(1);
    }
    FreeMarsaXorshift(p);
  }
  fclose(f);
}

/* Sets the type2p and type2q shift arrays from (a, b, c) according to the given type. */
static void setParameters(int type, int a, int b, int c, int *type2p, int *type2q) {
  if (type == TYPE21) {
    type2p[0] = -a; type2p[1] = b;  type2p[2] = 0;  type2q[0] = c;  type2q[1] = 0; type2q[2] = 0;
  } else if (type == TYPE22) {
    type2p[0] = c;  type2p[1] = 0;  type2p[2] = 0;  type2q[0] = -a; type2q[1] = b; type2q[2] = 0;
  } else if (type == TYPE23) {
    type2p[0] = a;  type2p[1] = b;  type2p[2] = 0;  type2q[0] = -c; type2q[1] = 0; type2q[2] = 0;
  } else if (type == TYPE24) {
    type2p[0] = -c; type2p[1] = 0;  type2p[2] = 0;  type2q[0] = a;  type2q[1] = b; type2q[2] = 0;
  } else if (type == TYPE25) {
    type2p[0] = -a; type2p[1] = b;  type2p[2] = -c; type2q[0] = 0;  type2q[1] = 0; type2q[2] = 0;
  }
}

/* Searches for nb full-period generators of the given type and writes results to filename. */
void ProduceParamsMarsaXorshift(char *filename, int nb, int type, int w, int r, int m, int nbxorshift) {
  FILE *f, *f2;
  int *Poly, *type2q, *type2p;
  BitVect Dummy;
  Generateur *p;
  int j, jj, count, a = 1, b = 1, c = 1, d = 1, aa, bb, cc, dd, mm;
  int alea_m = FALSE;
  int exhaustive_m = FALSE;
  int exhaustive_type = FALSE;
  int exhaustive_type3 = FALSE;
  char filenametemp[100];
  int numberGenToAllocate = 0;
  int *ai, *mi, *mi_nb, nbmi, aicounter, i;
  char Line[255], *operation;
  int type3N, *type3No, *type3Shift;

  sprintf(filenametemp, "%s.temp", filename);
  f = fopen(filenametemp, "w");
  count = 0;

  SetMRG32k3a((double)time(NULL), (double)time(NULL), (double)time(NULL),
               (double)time(NULL), (double)time(NULL), (double)time(NULL));

  if (type == TYPEGENERAL) {
    ai = (int *)malloc(nbxorshift * sizeof(int));
    operation = (char*)malloc(nbxorshift * sizeof(char));
    /* nbmi not yet known */
    mi = (int*) malloc(nbxorshift * sizeof(int));
    mi_nb = (int*) malloc(nbxorshift * sizeof(int));
  } else if (type == TYPE3) {
    if (m < 0) {
      type3N = -m;
      exhaustive_type3 = TRUE;
      if (r < type3N) {
        printf(" r must be greater than N for type 3 generators\n");
        exit(1);
      }
    } else {
      type3N = m;
    }
    type3No = (int*) calloc(type3N, sizeof(int));
    type3Shift = (int*) calloc(type3N, sizeof(int));
    for (j = 0; j < type3N; j++) {
      type3No[j] = type3N - 1 - j;
      type3Shift[j] = -w + 1;
    }
    type3No[type3N - 1] = r;
  } else if (type >= TYPE21 && type <= TYPE25) {
    type2p = (int*) calloc(3, sizeof(int));
    type2q = (int*) calloc(3, sizeof(int));
    for (j = 0; j < 3; j++) {
      type2q[j] = 0;
      type2p[j] = 0;
    }
    if (m == -1) { /* Random m and parameters */
      alea_m = TRUE;
    } else if (m == -2) { /* Exhaustive search over m, fixed type */
      exhaustive_m = TRUE;
      m = 1;
    } else if (m == -3) { /* Exhaustive search over m and type */
      exhaustive_type = TRUE;
      m = 1;
      type = TYPE21;
    }
  } else if (type == TYPE1) {
    r = 1;
  } else if (type == TYPE4) {
    type2p = (int*) calloc(2, sizeof(int));
    type2q = (int*) calloc(2, sizeof(int));
    for (j = 0; j < 2; j++) {
      type2q[j] = 0;
      type2p[j] = 0;
    }
    if (m == -1) {
      alea_m = TRUE;
    } else if (m == -2) {
      a = b = c = d = -w + 1;
      exhaustive_m = TRUE;
      mm = m = 1;
    }
  } else {
    printf("In ProduceParamsMarsaXorshift() : type not recognized\n");
    exit(1);
  }

  Poly = (int*) calloc((w * r + 1), sizeof(int));
  AllocBV(&Dummy, w * r + 1);

  while (TRUE) {
    if ((p = AllocMarsaXorshift(w * r)) == NULL) {
      printf("Error in ProduceParamsMarsaXorshift()\n");
      exit(1);
    }
    if (type == TYPEGENERAL) {
      { /* Multiple occurrences of mi indices allowed */
        int *indices;
        int ind;
        int nbright;

        if (count == nb)
          break;
        indices = (int*) calloc(r, sizeof(int));
        indices[r - 1] = 1;
        nbmi = 1;
#define PROB_OPERATIONPLUS 0.30
#define PROB_LEFT 0.80
#define NB_RIGHT 3
        nbright = 0;
        for (j = 0; j < nbxorshift - 1; j++) {
          ind = MRG32k3a() % r;
          if (indices[ind] == 0)
            nbmi++;
          indices[ind]++;
          ai[j] = MRG32k3a() % (w - 1) + 1;
          if (nbright < NB_RIGHT) {
            ai[j] *= -1;
            nbright++;
          }
          if (doubleMRG32k3a() < PROB_OPERATIONPLUS)
            operation[j] = '+';
          else
            operation[j] = '*';
        }
        /* One more for r-1 */
        ai[j] = MRG32k3a() % (w - 1) + 1;
        if (doubleMRG32k3a() < PROB_LEFT)
          ai[j] *= -1;
        if (doubleMRG32k3a() < PROB_OPERATIONPLUS)
          operation[j] = '+';
        else
          operation[j] = '*';

        nbmi = 0;
        aicounter = 0;
        for (j = 0; j < r; j++) {
          if (indices[j] > 0) {
            mi[nbmi] = j + 1;
            mi_nb[nbmi] = 0;
            for (i = 0; i < indices[j]; i++) {
              if (mi_nb[nbmi] == 0 || operation[aicounter] == '*') {
                mi_nb[nbmi]++;
              } else { /* operation[aicounter] == '+' */
                nbmi++;
                mi[nbmi] = j + 1;
                mi_nb[nbmi] = 1;
              }
              aicounter++;
            }
            nbmi++;
          }
        }
        free(indices);
      }
      InitParamMarsaXorshiftTypeGeneral(p, w, r, nbxorshift, nbmi, mi, mi_nb, ai, WL);
    } else if (type == TYPE3) {
      if (exhaustive_type3) {
        type3Shift[0]++;
        for (j = 0; j < type3N - 1; j++)
          if (type3Shift[j] == w) { type3Shift[j] = -w + 1; type3Shift[j + 1]++; }
        if (type3Shift[type3N - 1] == w) { type3Shift[type3N - 1] = -w + 1; type3No[0]++; }
        for (j = 0; j < type3N - 1; j++)
          if (type3No[j] == r) { type3No[j + 1]++; type3No[j] = type3No[j + 1] + 1; }
        if (type3No[0] == r)
          break;
      } else {
        if (count == nb)
          break;
        for (j = 0; j < m; j++) {
          type3Shift[j] = 0;
          type3No[j] = (MRG32k3a() % r) + 1;
          type3Shift[j] = (MRG32k3a() % (w - 1)) + 1;
          if (doubleMRG32k3a() < 0.90)
            type3Shift[j] *= -1;
        }
        type3No[m - 1] = r;
      }
      InitParamMarsaXorshiftTypeIII(p, w, r, type3N, type3No, type3Shift, WL);
    } else if (type >= TYPE21 && type <= TYPE25) {
      if (alea_m) { /* Random m and parameters */
        if (count == nb)
          break;
        m = (MRG32k3a() % (r - 1)) + 1;
        a = (MRG32k3a() % (w - 1)) + 1;
        b = (MRG32k3a() % (w - 1)) + 1;
        c = (MRG32k3a() % (w - 1)) + 1;
        setParameters(type, a, b, c, type2p, type2q);
      } else if (exhaustive_m) {
        setParameters(type, a, b, c, type2p, type2q);
        c++;
        if (c == w) { c = 1; b++; }
        if (b == w) { b = 1; a++; }
        if (a == w) { a = 1; m++; }
        if (m == r) break;
      } else if (exhaustive_type) {
        setParameters(type, a, b, c, type2p, type2q);
        c++;
        if (c == w) { c = 1; b++; }
        if (b == w) { b = 1; a++; }
        if (a == w) { a = 1; m++; }
        if (m == r) { m = 1; type++; }
        if (type == TYPE25 + 1) break;
      } else { /* Exhaustive search over a, b, c for fixed m, r, w, type */
        setParameters(type, a, b, c, type2p, type2q);
        c++;
        if (c == w) { c = 1; b++; }
        if (b == w) { b = 1; a++; }
        if (a == w) break;
      }
      InitParamMarsaXorshiftTypeII(p, type, w, r, m, type2p, type2q, WL);
    } else if (type == TYPE1) {
      c++;
      if (c == w) { c = 1; b++; }
      if (b == w) { b = 1; a++; }
      if (a == w) break;
      InitParamMarsaXorshiftTypeI(p, w, a, b, c, WL);
    } else if (type == TYPE4) {
      if (alea_m) { /* Random m and parameters */
        if (count == nb)
          break;
        m = (MRG32k3a() % (r - 1)) + 1;
        a = (MRG32k3a() % (2 * w - 1)) - (w - 1);
        b = (MRG32k3a() % (2 * w - 1)) - (w - 1);
        c = (MRG32k3a() % (2 * w - 1)) - (w - 1);
        d = (MRG32k3a() % (2 * w - 1)) - (w - 1);
        InitParamMarsaXorshiftTypeIV(p, aa = a, bb = b, cc = c, dd = d, w, r, mm = m, WL);
      } else if (exhaustive_m) {
        if (!(((a < 0) && (b < 0) && (c < 0) && (d < 0)) || ((a > 0) && (b > 0) && (c > 0) && (d > 0)))) {
          InitParamMarsaXorshiftTypeIV(p, aa = a, bb = b, cc = c, dd = d, w, r, mm = m, WL);
        } else {
          do {
            d++;
            if (d == w) { d = -w + 1; c++; }
            if (c == w) { c = -w + 1; b++; }
            if (b == w) { b = -w + 1; a++; }
            if (a == w) { a = -w + 1; m++; }
            if (m == r) break;
          } while (((a < 0) && (b < 0) && (c < 0) && (d < 0)) || ((a > 0) && (b > 0) && (c > 0) && (d > 0)));
          InitParamMarsaXorshiftTypeIV(p, aa = a, bb = b, cc = c, dd = d, w, r, mm = m, WL);
        }
        d++;
        if (d == w) { d = -w + 1; c++; }
        if (c == w) { c = -w + 1; b++; }
        if (b == w) { b = -w + 1; a++; }
        if (a == w) { a = -w + 1; m++; }
        if (m == r) break;
      } else { /* Exhaustive search over a, b, c, d for fixed m, r, w */
        if (!(((a < 0) && (b < 0) && (c < 0) && (d < 0)) || ((a > 0) && (b > 0) && (c > 0) && (d > 0)))) {
          InitParamMarsaXorshiftTypeIV(p, aa = a, bb = b, cc = c, dd = d, w, r, mm = m, WL);
        } else {
          do {
            d++;
            if (d == w) { d = -w + 1; c++; }
            if (c == w) { c = -w + 1; b++; }
            if (b == w) { b = -w + 1; a++; }
            if (a == w) break;
          } while (((a < 0) && (b < 0) && (c < 0) && (d < 0)) || ((a > 0) && (b > 0) && (c > 0) && (d > 0)));
          InitParamMarsaXorshiftTypeIV(p, aa = a, bb = b, cc = c, dd = d, w, r, mm = m, WL);
        }
        d++;
        if (d == w) { d = -w + 1; c++; }
        if (c == w) { c = -w + 1; b++; }
        if (b == w) { b = -w + 1; a++; }
        if (a == w) break;
      }
    }

    polychar(p, Poly, &Dummy);
    FreeMarsaXorshift(p);
    if (PrimitifPolynomial(Poly, w * r)) {
      if (type == TYPEGENERAL) {
        fprintf(f, "%d %d %d ", type, nbmi, nbxorshift);
        aicounter = 0;
        for (i = 0; i < nbmi; i++) {
          fprintf(f, "%d ", mi[i]);
          fprintf(f, "%d ", mi_nb[i]);
          for (j = 0; j < mi_nb[i]; j++)
            fprintf(f, "%d ", ai[aicounter++]);
        }
        fprintf(f, "\n");
      }

      if (type >= TYPE21 && type <= TYPE25)
        fprintf(f, "%3d %3d %3d %3d %3d %3d %3d %3d\n", type, m, type2p[0], type2p[1], type2p[2], type2q[0], type2q[1], type2q[2]);
      else if (type == TYPE1)
        fprintf(f, "%3d %3d %3d %3d\n", type, a, b, c);
      else if (type == TYPE3) {
        fprintf(f, "%3d %3d ", type, type3N);
        for (jj = 0; jj < type3N; jj++)
          fprintf(f, "%3d %3d ", type3No[jj], type3Shift[jj]);
        fprintf(f, "\n");
      } else if (type == TYPE4) {
        fprintf(f, "%3d %3d %3d %3d %3d %d %3d\n", 4, r, mm, aa, bb, cc, dd);
      }
      if (type == TYPE1 || type == TYPE21 || type == TYPE22 || type == TYPE25)
        numberGenToAllocate += 4;
      else if (type == TYPE23 || type == TYPE24)
        numberGenToAllocate += 2;
      else if (type == TYPE3 || type == TYPE4 || type == TYPEGENERAL)
        numberGenToAllocate += 1;
      count++;
      fflush(f);
    }
  }
  fclose(f);

  if (type == TYPEGENERAL) {
    free(ai);
    free(operation);
    free(mi);
    free(mi_nb);
  } else if ((type >= TYPE21 && type <= TYPE25) || type == TYPE4) {
    free(type2p);
    free(type2q);
  } else if (type == TYPE3) {
    free(type3No);
    free(type3Shift);
  }

  f = fopen(filenametemp, "r");
  f2 = fopen(filename, "w");
  fprintf(f2, "marsaxorshift\n");
  fprintf(f2, "%d %d %d\n", count, r, w);
  fprintf(f2, "%d # number of generators in this file\n", numberGenToAllocate);
  while (NULL != fgets(Line, 255, f))
    fprintf(f2, "%s", Line);
  fclose(f2);
  fclose(f);
  free(Poly);
  FreeBV(&Dummy);
}
