/*
 * genf2w.c — Generator over the extension field F_{2^w}.
 *
 * Implements PRNG components that operate as linear recurrences over GF(2^w),
 * supporting both polynomial-basis and normal-basis arithmetic. Provides LFSR
 * and polynomial-LCG variants with configurable coefficients and modulus.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <time.h>
#include <genf2w.h>
#include <polynomials.h>

#define GENF2WNAME        "Generator that uses F_{2^w}"
#define GENF2WPOLYLCGNAME "Polynomial LCG in F_{2^w}[z]/P(z)"
#define GENF2WLFSRNAME    "LFSR in F_{2^w}"

static void MakeTable(uint32_t *table, int w, uint32_t modM);
static uint32_t multiplypolynomialbasis(uint32_t a, uint32_t b, int w, uint32_t modM, uint32_t *dummy);
static uint32_t multiplynormalbasis(uint32_t a, uint32_t b, int w, uint32_t dummy, uint32_t *table);
static boolean isNormalBasis(int w, uint32_t modM);
static boolean minimalpolynomialprimitif(int w, uint32_t a, uint32_t modp);
static char* GenF2wPolyLCGName(void) { return GENF2WPOLYLCGNAME; }
static char* GenF2wLFSRName(void) { return GENF2WLFSRNAME; }
void GenF2wLFSR32bit(Generateur *Gen, BitVect *retour);
void InitGenF2wLFSR(Generateur *Gen, BitVect *init, BitVect *retour);

/* Returns the name string for this generator type. */
char* GenF2wName(void) {
  return GENF2WNAME;
}

/* Initializes all parameters of an F_{2^w} generator (LFSR or polynomial LCG). */
void InitParamGenF2w(Generateur *Gen, int w, int r, int nbcoeff, int *nocoeff, uint32_t *coeff,
         uint32_t modM, boolean normalbasis, int type, int step, int L) {
  int j;
  ((paramGenF2w*)(Gen->ParamGen))->type = type;
  ((paramGenF2w*)(Gen->ParamGen))->r = r;
  ((paramGenF2w*)(Gen->ParamGen))->w = w;
  ((paramGenF2w*)(Gen->ParamGen))->nbcoeff = nbcoeff;
  if (((paramGenF2w*)(Gen->ParamGen))->nocoeff != NULL) {
    free(((paramGenF2w*)(Gen->ParamGen))->nocoeff);
    free(((paramGenF2w*)(Gen->ParamGen))->coeff);
    free(((paramGenF2w*)(Gen->ParamGen))->table);
  }

  ((paramGenF2w*)(Gen->ParamGen))->nocoeff = (int *) malloc(nbcoeff * sizeof(int));
  ((paramGenF2w*)(Gen->ParamGen))->coeff = (uint32_t *) malloc(nbcoeff * sizeof(uint32_t));

  if (normalbasis) {
    ((paramGenF2w*)(Gen->ParamGen))->table = (uint32_t *) malloc(w * w * sizeof(uint32_t));
    MakeTable(((paramGenF2w*)(Gen->ParamGen))->table, w, modM);
    ((paramGenF2w*)(Gen->ParamGen))->multiply = multiplynormalbasis;
  } else {
    ((paramGenF2w*)(Gen->ParamGen))->multiply = multiplypolynomialbasis;
  }

  for (j = 0; j < nbcoeff; j++) {
    ((paramGenF2w*)(Gen->ParamGen))->nocoeff[j] = nocoeff[j];
    ((paramGenF2w*)(Gen->ParamGen))->coeff[j] = coeff[j];
  }
  ((paramGenF2w*)(Gen->ParamGen))->step = step;
  ((paramGenF2w*)(Gen->ParamGen))->modM = modM;
  ((paramGenF2w*)(Gen->ParamGen))->masque = (~0U) >> (WL - w);

  Gen->k = r * w;
  Gen->L = intmin(L, r * w);
  Gen->Step = 1;
  Gen->smax = INT_MAX;
  Gen->DispGen = DispGenF2w;
  Gen->DispName = GenF2wName;
  Gen->InitGen = InitGenF2w;

  if (type == GENF2WLFSR) {
    if (r * w < WL) {
      Gen->L = ((L - r * w) / w) * w + r * w;
      Gen->InitGen = InitGenF2wLFSR;
    }
    if (((paramGenF2w*)(Gen->ParamGen))->w == WL) {
      Gen->Iteration = GenF2wLFSR32bit;
      if (step != 1) {
  printf("Stepped iterations not implemented for LFSRs\n");
  exit(1);
      }
    } else {
      Gen->Iteration = GenF2wLFSR;
    }
  } else {
    Gen->Iteration = GenF2wPolyLCG;
  }
  Gen->CopyGen = CopyGenF2w;
  Gen->AllocGen = AllocGenF2w;
  Gen->PolyChar = polychar;
  PutBVToZero(&(Gen->GenState));
}

#define GENTYPE ((paramGenF2w*)(Gen->ParamGen))->type
#define RR      ((paramGenF2w*)(Gen->ParamGen))->r
#define WW      ((paramGenF2w*)(Gen->ParamGen))->w
#define NBCOEFF ((paramGenF2w*)(Gen->ParamGen))->nbcoeff
#define MULTIPLY ((paramGenF2w*)(Gen->ParamGen))->multiply
#define TABLE   ((paramGenF2w*)(Gen->ParamGen))->table
#define COEFF(i) ((paramGenF2w*)(Gen->ParamGen))->coeff[i]
#define NOCOEFF(i) ((paramGenF2w*)(Gen->ParamGen))->nocoeff[i]
#define MODP    ((paramGenF2w*)(Gen->ParamGen))->modM
#define MASKDP  ((paramGenF2w*)(Gen->ParamGen))->masque
#define STATE   Gen->GenState
#define PP      ((paramGenF2w*)(Gen->ParamGen))->p
#define STEP    ((paramGenF2w*)(Gen->ParamGen))->step

/* Displays the parameters of the generator (type, recurrence polynomial, step). */
void DispGenF2w(Generateur *Gen) {
  int j;
  if (GENTYPE == GENF2WLFSR)
    printf("%s\n", GenF2wLFSRName());
  else
    printf("%s\n", GenF2wPolyLCGName());
  printf("z^%d + ", RR);
  for (j = 0; j < NBCOEFF - 1; j++)
    printf("(%08x)z^%d +", COEFF(j), NOCOEFF(j));
  if (NOCOEFF(j) == 0)
    printf("(%08x) -- coeffs in F_2/P(z) where P(z)=%08x\n", COEFF(j), MODP);
  else
    printf("(%08x)z^%d -- coeffs in F_2/P(z) where P(z)=%08x\n", COEFF(j), NOCOEFF(j), MODP);
  printf("Step = %d\n", STEP);
}

/* Appends the generator parameters to a file in SSJ-compatible format. */
void PrintUsableFormatSSJ(Generateur *Gen, char *filename) {
  int j;
  FILE *f;
  f = fopen(filename, "a");
  if (f == NULL) {
    printf("Cannot create file \"%s\"\n", filename);
    exit(1);
  }
  fprintf(f, "%d %d %d %d %d ", WW, RR, MODP, STEP, NBCOEFF);
  for (j = 0; j < NBCOEFF; j++)
    fprintf(f, "%d %d ", COEFF(j), NOCOEFF(j));
  fprintf(f, "\n");
  fclose(f);
}

/* Appends the generator parameters to a file in REGPOLY-compatible format. */
void PrintUsableFormatREGPOLY(Generateur *Gen, char *filename) {
  int j;
  FILE *f;
  f = fopen(filename, "a");
  if (f == NULL) {
    printf("Cannot create file \"%s\"\n", filename);
    exit(1);
  }
  fprintf(f, "%d %d %x %d %d ", WW, RR, MODP, STEP, NBCOEFF);
  for (j = 0; j < NBCOEFF; j++)
    fprintf(f, "%x %d ", COEFF(j), NOCOEFF(j));
  fprintf(f, "\n");
  fclose(f);
}

static int t, k, s;

/* Returns the i-th WW-bit block from the left of the generator state. */
static uint32_t V(Generateur *Gen, int i) {
  if (WW == WL)
    return STATE.vect[i];
  else {
    t = i * WW;
    k = (t) / WL;
    s = (t - WL * k);
    return (STATE.vect[k] >> (WL - WW - s)) & MASKDP;
  }
}

/* Places the WW-bit value val into the i-th block from the left of the generator state. */
static void SetV(Generateur *Gen, int i, uint32_t val) {
  if (WW == WL)
    STATE.vect[i] = val;
  else {
    t = i * WW;
    k = (t) / WL;
    s = (t - WL * k);
    val &= MASKDP;
    val <<= (WL - WW - s);
    STATE.vect[k] &= ~((MASKDP) << (WL - WW - s));
    STATE.vect[k] |= val;
  }
}

/* Initializes the generator state from init and copies the first L bits to retour. */
void InitGenF2w(Generateur *Gen, BitVect *init, BitVect *retour) {
  CopyBV(&STATE, init);
  CopyBVPart(retour, init, Gen->L);
  PP = RR;
}

/* Performs one iteration (or STEP iterations) of the polynomial LCG in F_{2^w}. */
void GenF2wPolyLCG(Generateur *Gen, BitVect *retour) {
  uint32_t VP, res;
  int i, j;

  for (i = 0; i < STEP; i++) {
    VP = V(Gen, RR - 1);
    BVRShiftSelf(&STATE, WW);
    if (VP) {
      for (j = 0; j < NBCOEFF; j++) {
  res = MULTIPLY(VP, COEFF(j), WW, MODP, TABLE);
  res ^= V(Gen, NOCOEFF(j));
  SetV(Gen, NOCOEFF(j), res);
      }
    }
  }
  CopyBVPart(retour, &STATE, Gen->L);
}

/*
 * Old version (reverse order) -- kept for reference:
 *
 * void GenF2wPolyLCG(Generateur *Gen, BitVect *retour) {
 *   uint32_t VP, res;
 *   int i, j;
 *   for (i = 0; i < STEP; i++) {
 *     VP = V(Gen, 0);
 *     BVLShiftSelf(&STATE, WW);
 *     if (VP) {
 *       for (j = 0; j < NBCOEFF; j++) {
 *         res = MULTIPLY(VP, COEFF(j), WW, MODP, TABLE);
 *         res ^= V(Gen, RR - 1 - NOCOEFF(j));
 *         SetV(Gen, RR - 1 - NOCOEFF(j), res);
 *       }
 *     }
 *   }
 *   CopyBVPart(retour, &STATE, Gen->L);
 * }
 */

/* Recurrence: X_n = sum_{j=0}^{nbcoeff-1} COEFF(j) X_{n-r+nocoeff(j)} */
/* Char poly:  X^RR + sum_{j=0}^{nbcoeff-1} COEFF(j) X^{nocoeff(j)}    */

/* Performs a full block update and outputs L bits; used when w == WL (32-bit case). */
void GenF2wLFSR32bit(Generateur *Gen, BitVect *retour) {
  int i, j;
  uint32_t res;

  if (PP == RR) {
    PP = 0;
    for (i = 0; i < RR; i++) {
      res = 0U;
      for (j = 0; j < NBCOEFF; j++)
  res ^= MULTIPLY(STATE.vect[RR - ((i + NOCOEFF(j)) % RR) - 1], COEFF(j), WW, MODP, TABLE);
      STATE.vect[RR - i - 1] = res;
    }
  }
  for (j = 0; j < (Gen->L - 1) / WL + 1; j++)
    retour->vect[j] = STATE.vect[(RR - 1 - PP + j) % STATE.n];
  PP++;
}

/* Initializes the LFSR state and precomputes the extended output portion. */
void InitGenF2wLFSR(Generateur *Gen, BitVect *init, BitVect *retour) {
  int j, m, p;
  uint32_t res;
  CopyBV(&STATE, init);
  p = Gen->L / WW - RR;
  for (m = 0; m < p; m++) {
    res = 0U;
    for (j = 0; j < NBCOEFF; j++)
      res ^= MULTIPLY(V(Gen, NOCOEFF(j) + m), COEFF(j), WW, MODP, TABLE);
    SetV(Gen, RR + m, res);
  }
  CopyBVPart(retour, init, Gen->L);
}

/* Performs one STEP iteration of the LFSR in F_{2^w}; state is (m_n,...,m_{n+r-1}). */
void GenF2wLFSR(Generateur *Gen, BitVect *retour) {
  uint32_t res = 0U;
  int m, p, i, j;
  p = Gen->L / WW - RR;
  for (i = 0; i < STEP; i++) {
    res = 0;
    for (j = 0; j < NBCOEFF; j++)
      res ^= MULTIPLY(V(Gen, NOCOEFF(j)), COEFF(j), WW, MODP, TABLE);
    BVLShiftSelf(&STATE, WW);
    SetV(Gen, RR - 1, res);
  }
  for (m = 0; m < p; m++) {
    res = 0U;
    for (j = 0; j < NBCOEFF; j++)
      res ^= MULTIPLY(V(Gen, NOCOEFF(j) + m), COEFF(j), WW, MODP, TABLE);
    SetV(Gen, RR + m, res);
  }
  CopyBVPart(retour, &STATE, Gen->L);
}

/* Copies all parameters and state from Gen2 into Gen1 (both must be GenF2w). */
void CopyGenF2w(Generateur *Gen1, Generateur *Gen2) {
  int j, w;
  boolean normalbasis;
  ((paramGenF2w *)(Gen1->ParamGen))->type = ((paramGenF2w *)(Gen2->ParamGen))->type;
  w = ((paramGenF2w *)(Gen1->ParamGen))->w = ((paramGenF2w *)(Gen2->ParamGen))->w;
  ((paramGenF2w *)(Gen1->ParamGen))->r = ((paramGenF2w *)(Gen2->ParamGen))->r;
  ((paramGenF2w *)(Gen1->ParamGen))->p = ((paramGenF2w *)(Gen2->ParamGen))->p;
  ((paramGenF2w *)(Gen1->ParamGen))->step = ((paramGenF2w *)(Gen2->ParamGen))->step;
  ((paramGenF2w *)(Gen1->ParamGen))->nbcoeff = ((paramGenF2w *)(Gen2->ParamGen))->nbcoeff;
  normalbasis = ((paramGenF2w *)(Gen1->ParamGen))->normalbasis = ((paramGenF2w *)(Gen2->ParamGen))->normalbasis;
  ((paramGenF2w *)(Gen1->ParamGen))->multiply = ((paramGenF2w *)(Gen2->ParamGen))->multiply;

  if (((paramGenF2w *)(Gen1->ParamGen))->nocoeff != NULL) {
    free(((paramGenF2w *)(Gen1->ParamGen))->nocoeff);
    free(((paramGenF2w *)(Gen1->ParamGen))->coeff);
    free(((paramGenF2w *)(Gen1->ParamGen))->table);
  }
  ((paramGenF2w*)(Gen1->ParamGen))->nocoeff = (int *) malloc(((paramGenF2w *)(Gen2->ParamGen))->nbcoeff * sizeof(int));
  ((paramGenF2w*)(Gen1->ParamGen))->coeff = (uint32_t *) malloc(((paramGenF2w *)(Gen2->ParamGen))->nbcoeff * sizeof(uint32_t));

  if (normalbasis) {
    ((paramGenF2w*)(Gen1->ParamGen))->table = (uint32_t *) malloc(w * w * sizeof(uint32_t));
    for (j = 0; j < w * w; j++)
      ((paramGenF2w*)(Gen1->ParamGen))->table[j] = ((paramGenF2w*)(Gen2->ParamGen))->table[j];
  }

  for (j = 0; j < ((paramGenF2w *)(Gen2->ParamGen))->nbcoeff; j++) {
    ((paramGenF2w *)(Gen1->ParamGen))->nocoeff[j] = ((paramGenF2w *)(Gen2->ParamGen))->nocoeff[j];
    ((paramGenF2w *)(Gen1->ParamGen))->coeff[j] = ((paramGenF2w *)(Gen2->ParamGen))->coeff[j];
  }
  ((paramGenF2w *)(Gen1->ParamGen))->modM = ((paramGenF2w *)(Gen2->ParamGen))->modM;
  ((paramGenF2w *)(Gen1->ParamGen))->masque = ((paramGenF2w *)(Gen2->ParamGen))->masque;

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

/* Allocates and returns a new GenF2w generator with state size k bits. */
Generateur* AllocGenF2w(int k) {
  Generateur *G;
  G = (Generateur *) calloc(1, sizeof(Generateur));
  AllocGen(G, k);
  G->ParamGen = (paramGenF2w *) calloc(1, sizeof(paramGenF2w));
  ((paramGenF2w *)G->ParamGen)->nocoeff = NULL;
  ((paramGenF2w *)G->ParamGen)->coeff = NULL;
  ((paramGenF2w *)G->ParamGen)->table = NULL;
  return G;
}

/* Frees all memory associated with a GenF2w generator. */
void FreeGenF2w(Generateur *G) {
  free(((paramGenF2w *)G->ParamGen)->nocoeff);
  free(((paramGenF2w *)G->ParamGen)->coeff);
  free(((paramGenF2w *)G->ParamGen)->table);
  free(G->ParamGen);
  FreeGen(G);
}

/* Reads generator parameters from a file and populates the Component E. */
void ReadDataGenF2w(Component *E, char *file, boolean same, int L) {
  FILE *f;
  int w, r, i, j, nbgen, nbcoeff, *nocoeff, type, step;
  uint32_t modM, *coeff;
  boolean normalbasis;
  Generateur *p;

  f = fopen(file, "r");
  if (f == NULL) {
    printf("File \"%s\" not found\n", file);
    exit(1);
  }
  ReadLn(f);
  fscanf(f, "%d", &type);
  ReadLn(f);
  fscanf(f, "%d", &normalbasis);
  ReadLn(f);
  fscanf(f, "%d", &nbgen);
  ReadLn(f);
  AllocGensInComponent(E, nbgen, same);

  for (i = 0; i < nbgen; i++) {
    fscanf(f, "%d %d %x %d", &w, &r, &modM, &step);
    if ((p = AllocGenF2w(w * r)) == NULL) {
      printf("Error in ReadDataDoublePoly()\n");
      exit(1);
    }
    fscanf(f, "%d", &nbcoeff);
    coeff = (uint32_t *) malloc(nbcoeff * sizeof(uint32_t));
    nocoeff = (int *) malloc(nbcoeff * sizeof(int));
    for (j = 0; j < nbcoeff; j++)
      fscanf(f, "%x %d", &coeff[j], &nocoeff[j]);
    InitParamGenF2w(p, w, r, nbcoeff, nocoeff, coeff, modM, normalbasis, type, step, L);
    AddGenInComponent(E, p);
    ReadLn(f);
    free(coeff);
    free(nocoeff);
    FreeGenF2w(p);
  }
}

/* Multiplies element a by z^k in GF(2^w) represented by modM (polynomial basis). */
static uint32_t multiplyz(uint32_t a, int k, uint32_t modM) {
  int i;
  if (k == 0)
    return a;
  else {
    for (i = 0; i < k; i++) {
      if (1U & a)
  a = (a >> 1) ^ modM;
      else
  a = a >> 1;
    }
    return a;
  }
}

/* Multiplies two elements a and b of GF(2^w) using the polynomial basis. */
static uint32_t multiplypolynomialbasis(uint32_t a, uint32_t b, int w, uint32_t modM, uint32_t *dummy) {
  int i;
  uint32_t res = 0, verif = 1U;

  for (i = 0; i < w; i++) {
    if (b & verif)
      res ^= multiplyz(a, w - i - 1, modM);
    verif <<= 1;
  }
  return res;
}

/* Multiplies two elements a and b of GF(2^w) using the normal basis (precomputed table). */
static uint32_t multiplynormalbasis(uint32_t a, uint32_t b, int w, uint32_t dummy, uint32_t *table) {
  uint32_t res = 0U, verifa, verifb;
  int j, i;
  if (a == 0 || b == 0)
    return 0U;

  verifa = 1U << (w - 1);
  for (j = 0; j < w; j++) {
    verifb = 1U << (w - 1);
    for (i = 0; i < w; i++) {
      if (verifa & a)
  if (verifb & b)
    res ^= table[i * w + j];
      verifb >>= 1;
    }
    verifa >>= 1;
  }
  return res;
}

/* Builds the multiplication table for the normal basis representation of GF(2^w). */
static void MakeTable(uint32_t *table, int w, uint32_t modM) {
  uint32_t *x, *xx, verif, temp;
  Matrix C, InvC, TransInvC; /* C -> Conversion normal to polynomial basis */
  int i, j, jj;

  x = (uint32_t *) malloc(2 * w * sizeof(uint32_t));  /* x[j]  --> x^j */
  xx = (uint32_t *) malloc(w * sizeof(uint32_t));      /* xx[j] --> x^(2^j) */

  x[0] = 1U << (w - 1); /* x^0 */
  for (j = 0; j < 2 * w - 1; j++) {
    if (x[j] & 1U)
      x[j + 1] = (x[j] >> 1) ^ modM;
    else
      x[j + 1] = (x[j] >> 1);
  }
  xx[0] = x[1];
  for (j = 0; j < w - 1; j++) {
    verif = 1U << (w - 1);
    xx[j + 1] = 0U;
    for (i = 0; i < w; i++) {
      if (xx[j] & verif)
  xx[j + 1] ^= x[2 * i];
      verif >>= 1;
    }
  }
  AllocMat(&C, w, w, 1);
  AllocMat(&InvC, w, w, 1);
  AllocMat(&TransInvC, w, w, 1);

  for (j = 0; j < w; j++) {
    verif = 1U << (w - 1);
    for (i = 0; i < w; i++) {
      if (xx[j] & verif)
  PutBitBV(&(C.lignes[i][0]), j, 1);
      verif >>= 1;
    }
  }

  if (!InverseMatrix(&InvC, &C)) {
    printf("No normal basis!\n");
    exit(1);
  }

  TransposeMatrices(&TransInvC, &InvC);
  for (i = 0; i < w; i++) {
    for (j = 0; j < w; j++) {
      temp = multiplypolynomialbasis(xx[i], xx[j], w, modM, NULL);
      table[i * w + j] = 0U;
      for (jj = 0; jj < w; jj++)
  if ((temp >> (w - 1 - jj)) & 1U)
    table[i * w + j] ^= TransInvC.lignes[jj][0].vect[0];
    }
  }
  free(x);
  free(xx);
  FreeMat(&InvC);
  FreeMat(&TransInvC);
  FreeMat(&C);
}

/* Generates nb full-period F_{2^w} generators and writes their parameters to filename. */
void ProduceParamsGenF2w(char *filename, int nb, int w, int r, int nbcoeff, boolean normalbasis, int sizetable) {
  int *nocoeff;
  Generateur *p;
  BitVect Dummy, NoCoeffs;
  uint32_t modM, maskw, *coeff;
  int iteration;
  FILE *f;
  boolean goahead;
  int i, j, nbc, c, count = 0, *Poly, nc;

  SetMRG32k3a((double)time(NULL), (double)time(NULL), (double)time(NULL),
         (double)time(NULL), (double)time(NULL), (double)time(NULL));
  coeff = (uint32_t *) malloc(nbcoeff * sizeof(uint32_t));
  nocoeff = (int *) malloc(nbcoeff * sizeof(int));
  Poly = (int *) malloc((w + 1) * sizeof(int));

  if (!normalbasis)
    goahead = TRUE;
  do {
    do {
      for (i = 0; i <= w; i++)
  Poly[i] = 0;
      Poly[0] = 1;
      Poly[w] = 1;
      nbc = 2;
      for (i = 1; i < w; i++)
  if ((Poly[i] = MRG32k3a() % 2) == 1)
    nbc++;
      if (!(nbc % 2)) {
  c = MRG32k3a() % (w - 1) + 1;
  if (Poly[c])
    Poly[c] = 0;
  else
    Poly[c] = 1;
      }
    } while (!IrreduciblePolynomial(Poly, w));

    modM = 0U;
    for (i = 0; i < w; i++)
      if (Poly[i])
  modM ^= 0x80000000U >> i;
    modM >>= WL - w;

    if (normalbasis)
      goahead = isNormalBasis(w, modM);
  } while (!goahead && normalbasis);

  f = fopen(filename, "w");
  fprintf(f, "genf2w\n");
  fprintf(f, "%d\n", 0);
  fprintf(f, "%d\n", normalbasis);
  fprintf(f, "%d\n", nb);
  fclose(f);

  free(Poly);
  Poly = (int *) malloc((w * r + 1) * sizeof(int));
  AllocBV(&Dummy, w * r + 1);
  maskw = 0xffffffffU >> (WL - w);
  iteration = 0;
  while (count < nb) {
    AllocBV(&NoCoeffs, r);
    PutBVToZero(&NoCoeffs);
    for (j = 0; j < nbcoeff - 1; j++) {
      while (ValBitBV(&NoCoeffs, nc = MRG32k3a() % r) == 1);
      PutBitBV(&NoCoeffs, nc, 1);
    }
    PutBitBV(&NoCoeffs, 0, 1);
    i = 0;
    for (j = r - 1; j >= 0; j--) {
      if (ValBitBV(&NoCoeffs, j) == 1) {
  coeff[i] = MRG32k3a() & maskw & (0xffffffffU << (w - sizetable));
  nocoeff[i++] = j;
      }
    }

    /* Necessary condition for primitivity:
     * the constant term p_0 must be a primitive element of GF(2^w). */
    while (!minimalpolynomialprimitif(w, coeff[nbcoeff - 1], modM))
      coeff[nbcoeff - 1] = MRG32k3a() & (maskw & (0xffffffffU << (w - sizetable)));

    if ((p = AllocGenF2w(w * r)) == NULL) {
      printf("Error in ProduceParamsTGFSR()\n");
      exit(1);
    }

    InitParamGenF2w(p, w, r, nbcoeff, nocoeff, coeff, modM, normalbasis, GENF2WPOLYLCG, 1, WL);
    polychar(p, Poly, &Dummy);
    DispPolynomial(Poly, w * r);
    if (Poly[r * w] == 1)
      if (PrimitifPolynomial(Poly, r * w)) {
  f = fopen(filename, "a");
  fprintf(f, "%d %d %x %d %d ", w, r, modM, 1, nbcoeff);
  for (j = 0; j < nbcoeff; j++)
    fprintf(f, "%x %d ", coeff[j], nocoeff[j]);
  fprintf(f, "\n");
  fclose(f);
  count++;
      }
    FreeGenF2w(p);
  }
  free(Poly);
  FreeBV(&Dummy);
  free(coeff);
  printf("iterations:%d\n", iteration);
}

/* Returns TRUE if the polynomial defined by modM yields a normal basis for GF(2^w). */
static boolean isNormalBasis(int w, uint32_t modM) {
  uint32_t *x, *xx, verif;
  Matrix C, InvC; /* C -> Conversion normal to polynomial basis */
  int i, j;
  boolean res;

  x = (uint32_t *) malloc(2 * w * sizeof(uint32_t));  /* x[j]  --> x^j */
  xx = (uint32_t *) malloc(w * sizeof(uint32_t));      /* xx[j] --> x^(2^j) */

  x[0] = 1U << (w - 1); /* x^0 */
  for (j = 0; j < 2 * w - 1; j++) {
    if (x[j] & 1U)
      x[j + 1] = (x[j] >> 1) ^ modM;
    else
      x[j + 1] = (x[j] >> 1);
  }
  xx[0] = x[1];
  for (j = 0; j < w - 1; j++) {
    verif = 1U << (w - 1);
    xx[j + 1] = 0U;
    for (i = 0; i < w; i++) {
      if (xx[j] & verif)
  xx[j + 1] ^= x[2 * i];
      verif >>= 1;
    }
  }
  AllocMat(&C, w, w, 1);
  AllocMat(&InvC, w, w, 1);

  for (j = 0; j < w; j++) {
    verif = 1U << (w - 1);
    for (i = 0; i < w; i++) {
      if (xx[j] & verif)
  PutBitBV(&(C.lignes[i][0]), j, 1);
      verif >>= 1;
    }
  }
  free(x);
  free(xx);

  res = InverseMatrix(&InvC, &C);
  FreeMat(&InvC);
  FreeMat(&C);
  return res;
}

/* Returns TRUE if the minimal polynomial of a (in GF(2^w) mod modp) is primitive. */
static boolean minimalpolynomialprimitif(int w, uint32_t a, uint32_t modp) {
  char *Sequence;
  int j, *coeff;
  BitVect BVPoly;
  boolean primitif;
  uint32_t g = 1;

  Sequence = (char *) malloc(2 * w * sizeof(char));
  coeff = (int *) malloc((w + 1) * sizeof(int));
  AllocBV(&BVPoly, w + 1);
  for (j = 0; j < 2 * w; j++) {
    g = multiplypolynomialbasis(g, a, w, modp, NULL);
    Sequence[j] = g & 1U;
  }
  SequenceMinimalPolynomial(w, Sequence, coeff, &BVPoly);
  primitif = PrimitifPolynomial(coeff, w);
  return primitif;
}

/* Generates one F_{2^w} generator and adds it to the first component of C. */
void ProduceOneGenF2w(Combinaison *C, int genttype, int w, int r, int nbcoeff, boolean normalbasis, int sizetable, int stepmax, boolean primitive, int Lmax) {
  uint32_t modM, maskw, *coeff;
  int *nocoeff;
  int i, j, nbc, c, *Poly, nc;
  boolean goahead, finished;
  BitVect Dummy, NoCoeffs;
  Generateur *p;

  SetMRG32k3a((double)time(NULL), (double)time(NULL), (double)time(NULL),
         (double)time(NULL), (double)time(NULL), (double)time(NULL));
  coeff = (uint32_t *) malloc(nbcoeff * sizeof(uint32_t));
  nocoeff = (int *) malloc(nbcoeff * sizeof(int));
  AllocGensInComponent(&(C->Components[0]), 1, FALSE);

  Poly = (int *) malloc((w + 1) * sizeof(int));
  if (!normalbasis)
    goahead = TRUE;
  do {
    do {
      for (i = 0; i <= w; i++)
  Poly[i] = 0;
      Poly[0] = 1;
      Poly[w] = 1;
      nbc = 2;
      for (i = 1; i < w; i++)
  if ((Poly[i] = MRG32k3a() % 2) == 1)
    nbc++;
      if (!(nbc % 2)) {
  c = MRG32k3a() % (w - 1) + 1;
  if (Poly[c])
    Poly[c] = 0;
  else
    Poly[c] = 1;
      }
    } while (!IrreduciblePolynomial(Poly, w));
    modM = 0U;
    for (i = 0; i < w; i++)
      if (Poly[i])
  modM ^= 0x80000000U >> i;
    modM >>= WL - w;
    if (normalbasis)
      goahead = isNormalBasis(w, modM);
  } while (!goahead && normalbasis);
  free(Poly);

  maskw = 0xffffffffU >> (WL - w);
  AllocBV(&NoCoeffs, r);

  finished = FALSE;
  do {
    PutBVToZero(&NoCoeffs);
    for (j = 0; j < nbcoeff - 1; j++) {
      while (ValBitBV(&NoCoeffs, nc = (MRG32k3a() % (r - 1)) + 1) == 1);
      PutBitBV(&NoCoeffs, nc, 1);
    }
    PutBitBV(&NoCoeffs, 0, 1);
    i = 0;
    for (j = r - 1; j >= 0; j--) {
      if (ValBitBV(&NoCoeffs, j) == 1) {
  coeff[i] = MRG32k3a() & maskw & (0xffffffffU << (w - sizetable));
  nocoeff[i++] = j;
      }
    }
    if ((p = AllocGenF2w(w * r)) == NULL) {
      printf("Error in ProduceParamsTGFSR()\n");
      exit(1);
    }

    InitParamGenF2w(p, w, r, nbcoeff, nocoeff, coeff, modM, normalbasis, genttype, 1, WL);
    if (primitive) {
      Poly = (int *) malloc((w * r + 1) * sizeof(int));
      AllocBV(&Dummy, w * r + 1);
      polychar(p, Poly, &Dummy);
      if (Poly[r * w] == 1)
  if (PrimitifPolynomial(Poly, r * w))
    finished = TRUE;
    } else {
      finished = TRUE;
    }
    if (!finished)
      FreeGenF2w(p);
  } while (!finished);
  ((paramGenF2w*)((p)->ParamGen))->step = (MRG32k3a() % stepmax) + 1;
  free(coeff);
  free(nocoeff);
  FreeBV(&NoCoeffs);
  AddGenInComponent(&(C->Components[0]), p);
  FreeGenF2w(p);
}

/* Updates the parameters of Gen by random or exhaustive search; returns FALSE when search space exhausted. */
boolean UpdateParamsGenF2w(Generateur *Gen, int normalbasis, int sizetable, int stepmax, boolean primitive, int modM, int exhaustive) {
  uint32_t maskw;
  int i, j, nbc, c, *Poly, nc;
  boolean goahead, finished;
  BitVect Dummy, NoCoeffs;
  static uint32_t counter = 0;

  if (modM == -1) {
    Poly = (int *) malloc((WW + 1) * sizeof(int));
    if (!normalbasis)
      goahead = TRUE;
    do {
      do {
  for (i = 0; i <= WW; i++)
    Poly[i] = 0;
  Poly[0] = 1;
  Poly[WW] = 1;
  nbc = 2;
  for (i = 1; i < WW; i++)
    if ((Poly[i] = MRG32k3a() % 2) == 1)
      nbc++;
  if (!(nbc % 2)) {
    c = MRG32k3a() % (WW - 1) + 1;
    if (Poly[c])
      Poly[c] = 0;
    else
      Poly[c] = 1;
  }
      } while (!IrreduciblePolynomial(Poly, WW));
      MODP = 0U;
      for (i = 0; i < WW; i++)
  if (Poly[i])
    MODP ^= 0x80000000U >> i;
      MODP >>= WL - WW;
      if (normalbasis)
  goahead = isNormalBasis(WW, MODP);
    } while (!goahead && normalbasis);
    free(Poly);
  } else if (modM != -2) { /* -2 means keep the same modM */
    MODP = modM;
  }

  maskw = 0xffffffffU >> (WL - WW);
  AllocBV(&NoCoeffs, RR);

  finished = FALSE;
  do {
    if (exhaustive) {
      if (NBCOEFF != RR) {
  printf("Error. Cannot do exhaustive search if nbcoeff!=r.\n");
  exit(1);
      }
      i = 0;
      for (j = 0; j < RR; j++) {
  COEFF(i) = (counter >> (j * WW)) & (0xffffffffU >> (WL - WW));
  NOCOEFF(i++) = j;
      }
      counter++;
      if (counter > (1U << WW * RR))
  return FALSE;
    } else {
      PutBVToZero(&NoCoeffs);
      for (j = 0; j < NBCOEFF - 1; j++) {
  while (ValBitBV(&NoCoeffs, nc = MRG32k3a() % RR) == 1);
  PutBitBV(&NoCoeffs, nc, 1);
      }
      PutBitBV(&NoCoeffs, 0, 1);
      i = 0;
      for (j = RR - 1; j >= 0; j--) {
  if (ValBitBV(&NoCoeffs, j) == 1) {
    COEFF(i) = MRG32k3a() & maskw & (0xffffffffU << (WW - sizetable));
    NOCOEFF(i++) = j;
  }
      }
    }

    if (primitive) {
      Poly = (int *) malloc((WW * RR + 1) * sizeof(int));
      AllocBV(&Dummy, WW * RR + 1);
      STEP = 1;
      polychar(Gen, Poly, &Dummy);
      FreeBV(&Dummy);
      if (Poly[WW * RR] == 1)
  if (PrimitifPolynomial(Poly, WW * RR))
    finished = TRUE;
      free(Poly);
    } else {
      finished = TRUE;
    }
  } while (!finished);
  STEP = (MRG32k3a() % stepmax) + 1;
  FreeBV(&NoCoeffs);
  return TRUE;
}
