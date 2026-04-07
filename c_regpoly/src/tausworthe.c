#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <time.h>
#include "tausworthe.h"
#include "regpoly.h"
#include "polynomials.h"

#define TAUSNAME "Tausworthe Generator"

typedef BitVect *DivSets;

static void ProducePolyQuick(int *Poly, int degre, int nbcoeff);
static void ProducePoly(int *Poly, int degre, int nbcoeff);

/* IMPORTANT REMARK:
   Two files are essential for this module:
   'seektaus.dt1' contains, at line n, the factors of n.  It has NLIGNES lines.
   'seektaus.dt2' contains, at line n, the factors of 2^n-1 that are less than NLIGNES.
   The constant MAXNBCOEFF in tausworthe.h must also be updated if these files change. */
#define NLIGNES 256

/* Convenience macros to access Tausworthe parameters from a Generateur pointer. */
#define TAUSQUICKTAUS ((paramTaus *)(Gen->ParamGen))->quicktaus
#define TAUSK         Gen->k
#define TAUSQ         ((paramTaus *)(Gen->ParamGen))->Q
#define TAUSS         ((paramTaus *)(Gen->ParamGen))->s
#define TAUSNBC       ((paramTaus *)(Gen->ParamGen))->NbCoeff
#define TAUSKMS       ((paramTaus *)(Gen->ParamGen))->gen_kms
#define TAUSL         ((paramTaus *)(Gen->ParamGen))->StateL
#define TAUSSTATE     Gen->GenState

static void InitDiv(char *fich, DivSets ds);

/* Return the name string of the Tausworthe generator type. */
char *TausName(void) {
  return TAUSNAME;
}

/* Initialize the parameters of the Tausworthe generator pointed by Gen, with
   characteristic polynomial defined by Q[] (NbCoeff non-zero coefficients) and
   step parameter s.  Sets up function pointers for QuickTaus or general algorithm. */
void InitParamTaus(Generateur *Gen, int k, int Q[MAXNBCOEFF], int s, int NbCoeff,
                   boolean quicktaus, int L) {
  int i;

  TAUSS = s;
  TAUSNBC = NbCoeff;
  TAUSQUICKTAUS = quicktaus;
  for (i = 0; i < NbCoeff; i++)
    TAUSQ[i] = Q[i];

  if (k > L) {
    TAUSL = L;
  } else {
    TAUSL = L;
    FreeBV(&TAUSSTATE);
    AllocBV(&TAUSSTATE, TAUSL);
  }
  if (quicktaus) {
    TAUSKMS = k - s;
    Gen->Iteration = QuickTaus;
  } else {
    Gen->Iteration = GeneralTaus;
  }
  Gen->InitGen = InitGeneralTaus;
  Gen->k = k;
  Gen->L = L;
  Gen->Step = s;
  Gen->smax = INT_MAX;
  Gen->DispGen = DispTaus;
  Gen->DispName = TausName;
  Gen->CopyGen = CopyTaus;
  Gen->AllocGen = AllocTaus;
  Gen->PolyChar = polychar;
  PutBVToZero(&TAUSSTATE);
}

/* Copy all parameters from generator Gen2 into generator Gen1; both must be Tausworthe. */
void CopyTaus(Generateur *Gen1, Generateur *Gen2) {
  int i;
  ((paramTaus *)Gen1->ParamGen)->quicktaus = ((paramTaus *)Gen2->ParamGen)->quicktaus;
  for (i = 0; i < ((paramTaus *)Gen2->ParamGen)->NbCoeff; i++)
    ((paramTaus *)Gen1->ParamGen)->Q[i] = ((paramTaus *)Gen2->ParamGen)->Q[i];
  ((paramTaus *)Gen1->ParamGen)->s       = ((paramTaus *)Gen2->ParamGen)->s;
  ((paramTaus *)Gen1->ParamGen)->NbCoeff = ((paramTaus *)Gen2->ParamGen)->NbCoeff;
  ((paramTaus *)Gen1->ParamGen)->StateL  = ((paramTaus *)Gen2->ParamGen)->StateL;

  if (((paramTaus *)Gen2->ParamGen)->quicktaus) {
    ((paramTaus *)Gen1->ParamGen)->gen_kms = ((paramTaus *)Gen2->ParamGen)->gen_kms;
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

/* Allocate a new Generateur structure for a Tausworthe generator of degree k. */
Generateur *AllocTaus(int k) {
  Generateur *G;
  G = (Generateur *) malloc(sizeof(Generateur));
  AllocGen(G, k);
  G->ParamGen = (paramTaus *) malloc(sizeof(paramTaus));
  return G;
}

/* Free the memory associated with a Tausworthe generator previously allocated by AllocTaus(). */
void FreeTaus(Generateur *G) {
  free(G->ParamGen);
  FreeGen(G);
}

/* Initialize the state of the QuickTaus generator pointed by Gen to the bit vector init;
   place the resulting L-bit state in retour. */
void InitQuickTaus(Generateur *Gen, BitVect *init, BitVect *retour) {
  int ik, ir, r, i, j;
  BitVect iB, ib[MAXNBCOEFF], seed;

  AllocBV(&iB, TAUSL);
  AllocBV(&seed, TAUSL);

  for (j = 0; j < TAUSNBC - 1; j++)
    AllocBV(&ib[j], TAUSL);

  CopyBV(&seed, init);

  ANDBVMask(&TAUSSTATE, &seed, TAUSK);
  ik = TAUSK;
  ir = r = TAUSK - TAUSQ[TAUSNBC - 2];
  while (ik < TAUSL) {
    PutBVToZero(&iB);
    for (j = 0; j < TAUSNBC - 1; j++)
      BVLShift(&ib[j], &TAUSSTATE, TAUSQ[j]);

    for (j = 0; j < TAUSNBC - 1; j++)
      XORBV(&iB, &iB, &ib[j]);

    ANDBVMask(&iB, &iB, TAUSK - TAUSL + ir);
    BVRShiftSelf(&iB, TAUSK);
    for (i = 0; i < iB.n; i++)
      TAUSSTATE.vect[i] |= iB.vect[i];
    ik += r;
    ir += r;
  }
  CopyBV(retour, &TAUSSTATE);
  FreeBV(&iB);
  FreeBV(&seed);

  for (j = 0; j < TAUSNBC - 1; j++)
    FreeBV(&ib[j]);
}

/* Perform one QuickTaus iteration on the generator pointed by Gen; return the next
   L-bit output state in State. */
void QuickTaus(Generateur *Gen, BitVect *State) {
  BitVect B[MAXNBCOEFF];
  BitVect gen_B;
  register int j;

  for (j = 1; j < TAUSNBC - 1; j++)
    AllocBV(&B[j], TAUSL);
  AllocBV(&gen_B, TAUSL);

  PutBVToZero(&gen_B);
  for (j = 1; j < TAUSNBC - 1; j++) {
    BVLShift(&B[j], &TAUSSTATE, TAUSQ[j]);
  }
  for (j = 1; j < TAUSNBC - 1; j++)
    XORBVSelf(&gen_B, &B[j]);
  XORBVSelf(&gen_B, &TAUSSTATE);
  BVRShiftSelf(&gen_B, TAUSKMS);
  ANDBVMask(&TAUSSTATE, &TAUSSTATE, TAUSK);

  BVLShiftSelf(&TAUSSTATE, TAUSS);
  XORBVSelf(&TAUSSTATE, &gen_B);
  CopyBV(State, &TAUSSTATE);

  for (j = 1; j < TAUSNBC - 1; j++)
    FreeBV(&B[j]);
  FreeBV(&gen_B);
}

/* Display the characteristic polynomial and the step parameter s of the generator Gen. */
void DispTaus(Generateur *Gen) {
  int i;
  char output[3000] = "\0";
  char temp[3000];

  for (i = TAUSNBC - 1; i > 0; i--)
    if (TAUSQ[i] == 1)
      strcat(output, " x +");
    else {
      sprintf(temp, " x^%d +", TAUSQ[i]);
      strcat(output, temp);
    }
  strcat(output, " 1 ");
  printf("%-40s   s=%d\n", output, TAUSS);
}

/* Read a Tausworthe data file and fill the Component pointed by E with all valid
   (polynomial, s) pairs.  The parameter same indicates that the generators are the same
   as in the previous Component.  L is the output resolution. */
void ReadDataTaus(Component *E, char *nomfich, boolean meme, int L) {
  FILE *f;
  int nbpoly, kk, kmq, j, t, quicktaus, smax;
  int Count = 0;
  Generateur *pTaus;
  int k, s, NbCoeff, R[MAXNBCOEFF], Q[MAXNBCOEFF];
  DivSets s_div, twokm1_div;

  s_div = (BitVect *) malloc(NLIGNES * sizeof(BitVect));
  twokm1_div = (BitVect *) malloc(NLIGNES * sizeof(BitVect));

  for (j = 0; j < NLIGNES; j++) {
    AllocBV(&s_div[j], NLIGNES);
    AllocBV(&twokm1_div[j], NLIGNES);
  }

  f = fopen(nomfich, "r");
  if (f == NULL) {
    printf("file %s not found\n", nomfich);
    exit(1);
  }
  InitDiv("seektaus.dt1", s_div);
  InitDiv("seektaus.dt2", twokm1_div);
  ReadLn(f);
  fscanf(f, "%d %d", &nbpoly, &quicktaus);
  if (!quicktaus) {
    fscanf(f, "%d", &smax);
    if (smax > NLIGNES) {
      printf("The value of smax must be less than %d.\n\n", NLIGNES);
      exit(1);
    }
  }
  ReadLn(f);
  for (j = 0; j < nbpoly; j++) {
    fscanf(f, "%d", &t);
    k = t;
    fscanf(f, "%d", &t);
    kmq = k - t;
    ReadLn(f);
    if (quicktaus) {
      for (s = 1; s <= kmq; s++)
        if (k > NLIGNES)
          Count++;
        else if (!VerifBitsCommuns(&s_div[s], &twokm1_div[k]))
          Count++;
    } else {
      for (s = 1; s <= smax; s++)
        if (k > NLIGNES)
          Count++;
        else if (!VerifBitsCommuns(&s_div[s], &twokm1_div[k]))
          Count++;
    }
  }
  rewind(f);
  AllocGensInComponent(E, Count, meme);
  ReadLn(f);
  ReadLn(f);

  for (j = 0; j < nbpoly; j++) {
    NbCoeff = 0;
    t = 1;
    while (t != 0) {
      fscanf(f, "%d", &t);
      R[NbCoeff++] = t;
    }
    ReadLn(f);
    for (kk = 0; kk < NbCoeff; kk++)
      Q[NbCoeff - kk - 1] = R[kk];
    k = Q[NbCoeff - 1];
    kmq = Q[NbCoeff - 1] - Q[NbCoeff - 2];
    if ((pTaus = AllocTaus(k)) == NULL) {
      printf("Error in ReadDataTaus()\n");
      exit(1);
    }
    if (quicktaus) {
      for (s = 1; s <= kmq; s++)
        if (k > NLIGNES) {
          InitParamTaus(pTaus, k, Q, s, NbCoeff, TRUE, L);
          AddGenInComponent(E, pTaus);
        } else if (!VerifBitsCommuns(&s_div[s], &twokm1_div[k])) {
          InitParamTaus(pTaus, k, Q, s, NbCoeff, TRUE, L);
          AddGenInComponent(E, pTaus);
        }
    } else {
      for (s = 1; s <= smax; s++) {
        if (k > NLIGNES) {
          InitParamTaus(pTaus, k, Q, s, NbCoeff, FALSE, L);
          AddGenInComponent(E, pTaus);
        } else if (!VerifBitsCommuns(&s_div[s], &twokm1_div[k])) {
          InitParamTaus(pTaus, k, Q, s, NbCoeff, FALSE, L);
          AddGenInComponent(E, pTaus);
        }
      }
    }
    FreeTaus(pTaus);
  }

  fclose(f);
  for (j = 0; j < NLIGNES; j++) {
    FreeBV(&s_div[j]);
    FreeBV(&twokm1_div[j]);
  }
  free(s_div);
  free(twokm1_div);
}

/* Read the divisor file fich and populate the DivSets array ds with bit vectors
   encoding the divisors of each index. */
static void InitDiv(char *fich, DivSets ds) {
  BitVect BC;
  register int i, j, k;
  int n, p;
  FILE *f;

  AllocBV(&BC, NLIGNES);

  if (!(f = fopen(fich, "r"))) {
    printf("\n*** File %s not found. ***\n", fich);
    exit(1);
  }
  PutBVToZero(&ds[0]);
  for (i = 1; i < NLIGNES; i++) {
    PutBVToZero(&ds[i]);
    fscanf(f, "%d", &n);
    for (j = 0; j < n; j++) {
      fscanf(f, "%d", &p);
      k = p / WL;
      BVCanonic(&BC, p);
      ds[i].vect[k] |= BC.vect[k];
    }
  }
  fclose(f);
  FreeBV(&BC);
}

/* Generate nb random Tausworthe polynomials of degree k with nbcoeff non-zero coefficients
   and write them in a format readable by ReadDataTaus() to the file filename. */
void ProduceParamsTaus(char *filename, int nb, int k, int nbcoeff,
                       int prim_or_irr, boolean quicktaus) {
  FILE *f;
  int *Poly;
  int i, count;
  boolean OK;

  Poly = (int *) malloc((k + 1) * sizeof(int));
  f = fopen(filename, "w");

  SetMRG32k3a((double) time(NULL), (double) time(NULL), (double) time(NULL),
              (double) time(NULL), (double) time(NULL), (double) time(NULL));
  if ((nbcoeff < 3) && ((nbcoeff % 2) != 1) && (nbcoeff != -1)) {
    printf("Error in ProduceParamsTaus().\n");
    printf("The number of coefficients must be odd and greater than 2.\n");
    exit(1);
  }
  count = 0;
  fprintf(f, "taus2\n%d\n", nb);
  while (count < nb) {
    if (quicktaus)
      ProducePolyQuick(Poly, k, nbcoeff);
    else
      ProducePoly(Poly, k, nbcoeff);

    if (prim_or_irr == PRIMITIVE) {
      if (PrimitifPolynomial(Poly, k))
        OK = TRUE;
      else
        OK = FALSE;
    } else if (prim_or_irr == IRREDUCIBLE) {
      if (IrreduciblePolynomial(Poly, k))
        OK = TRUE;
      else
        OK = FALSE;
    } else {
      exit(1);
    }

    if (OK) {
      for (i = k; i >= 0; i--)
        if (Poly[i])
          fprintf(f, "%d ", i);
      count++;
      fprintf(f, "\n");
    }
  }
  fclose(f);
}

/* Initialize the state of a general Tausworthe generator pointed by Gen to init;
   fill retour with the resulting L-bit state. */
void InitGeneralTaus(Generateur *Gen, BitVect *init, BitVect *retour) {
  int i, j, sum;

  CopyBVPart(&TAUSSTATE, init, TAUSK);

  for (j = TAUSK; j < TAUSL; j++) {
    sum = 0;
    for (i = 0; i < TAUSNBC - 1; i++)
      sum += ValBitBV(&TAUSSTATE, j - (TAUSK - TAUSQ[i]));
    if (sum %= 2)
      PutBitBV(&TAUSSTATE, j, 1);
  }
  CopyBVPart(retour, &TAUSSTATE, TAUSL);
}

/* Perform one iteration of the general Tausworthe generator pointed by Gen;
   place the next L-bit output state in retour. */
void GeneralTaus(Generateur *Gen, BitVect *retour) {
  int i, j, sum, ss, m;
  /* TAUSSTATE is the generator state (extended to L bits).
     TAUSL is the resolution (typically 32 bits).
     TAUSK is the degree of the generator.
     ss is the number of bits to produce per inner loop iteration.
     TAUSS is the step parameter s. */
  ss = intmin(TAUSS, TAUSL - TAUSK);

  m = TAUSS;
  do { /* This loop produces ss bits of the recurrence at a time. */
    BVLShiftSelf(&TAUSSTATE, ss);
    for (j = TAUSL - ss; j < TAUSL; j++) {
      sum = 0;
      for (i = 0; i < TAUSNBC - 1; i++)
        sum += ValBitBV(&TAUSSTATE, j - (TAUSK - TAUSQ[i]));
      if (sum %= 2)
        PutBitBV(&TAUSSTATE, j, 1);
    }
    m -= ss;
    if (m < ss)
      ss = m;
  } while (m > 0);
  /* At this point, TAUSSTATE contains the next generator state (extended to L bits). */
  CopyBVPart(retour, &TAUSSTATE, TAUSL);
}

/* Fill Poly with a random polynomial of given degree and nbcoeff non-zero terms,
   using the general (unrestricted) algorithm. */
static void ProducePoly(int *Poly, int degre, int nbcoeff) {
  int i, c, nbc = 2;

  for (i = 0; i < degre; i++)
    Poly[i] = 0;
  Poly[0] = 1;
  Poly[degre] = 1;
  if (nbcoeff == -1)
    nbcoeff = (MRG32k3a() % ((degre - 1) / 2)) * 2 + 3;
  while (nbcoeff != nbc) {
    c = MRG32k3a() % degre;
    if (!Poly[c]) {
      Poly[c] = 1;
      nbc++;
    }
  }
}

/* Fill Poly with a random polynomial of given degree and nbcoeff non-zero terms,
   using the QuickTaus-compatible (symmetric coefficient placement) algorithm. */
static void ProducePolyQuick(int *Poly, int degree, int nbcoeff) {
  int i, nbc, c;

  for (i = 0; i < degree; i++)
    Poly[i] = 0;
  Poly[0] = 1;
  Poly[degree] = 1;
  nbc = 2;
  while (nbcoeff != nbc) {
    c = MRG32k3a() % ((degree + 1) / 2 - 1) + 1;
    if (!Poly[c]) {
      Poly[c] = 1;
      nbc++;
    }
  }
}
