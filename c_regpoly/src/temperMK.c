#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include "temperMK.h"
#include "basisreduc.h"
#include "polynomials.h"
#include "myzen.h"

static Matrix mat, Mat;
static boolean firsttime = TRUE;

/* Return the integer ID of the MK tempering transformation type. */
int GiveTemperMKID(void) {
  return TEMPMKID;
}

/* Convenience macros for accessing MK tempering parameters through a Transformation pointer. */
#define BB   ((paramTempMK *)(T->ParamTrans))->b
#define CC   ((paramTempMK *)(T->ParamTrans))->c
#define ETA  ((paramTempMK *)(T->ParamTrans))->eta
#define MU   ((paramTempMK *)(T->ParamTrans))->mu
#define UU   ((paramTempMK *)(T->ParamTrans))->u
#define LL   ((paramTempMK *)(T->ParamTrans))->l
#define TYPE ((paramTempMK *)(T->ParamTrans))->type

/* Initialize a Transformation structure with MK tempering parameters; set up all
   function pointers and copy the initial b and c bit vectors. */
void InitTemperMK(Transformation *T, int w, int eta, int mu, int u, int l, int alea,
                  BitVect *b, BitVect *c, boolean opt, int type,
                  boolean DispProgress, int limitv) {
  T->w_original = w;
  T->w = w;
  T->DispLinTrans = DispTemperMK;
  T->Trans = TemperingMK;
  T->ChangeParamTrans = UpdateTemperMK;
  T->AllocTrans = AllocTemperMK;
  T->CopyTrans = CopyTemperMK;
  T->GiveTransID = GiveTemperMKID;
  T->InverseTrans = InverseTemperMK;

  if (((paramTempMK *)T->ParamTrans)->b.n != 0) {
    FreeBV(&(((paramTempMK *)T->ParamTrans)->b));
    FreeBV(&(((paramTempMK *)T->ParamTrans)->b_original));
    FreeBV(&(((paramTempMK *)T->ParamTrans)->BestB));
    FreeBV(&(((paramTempMK *)T->ParamTrans)->c));
    FreeBV(&(((paramTempMK *)T->ParamTrans)->c_original));
    FreeBV(&(((paramTempMK *)T->ParamTrans)->BestC));
  }
  AllocBV(&(((paramTempMK *)T->ParamTrans)->b),          T->w);
  AllocBV(&(((paramTempMK *)T->ParamTrans)->b_original), T->w);
  AllocBV(&(((paramTempMK *)T->ParamTrans)->BestB),      T->w);
  AllocBV(&(((paramTempMK *)T->ParamTrans)->c),          T->w);
  AllocBV(&(((paramTempMK *)T->ParamTrans)->c_original), T->w);
  AllocBV(&(((paramTempMK *)T->ParamTrans)->BestC),      T->w);

  ((paramTempMK *)(T->ParamTrans))->type        = type;
  ((paramTempMK *)(T->ParamTrans))->l           = l;
  ((paramTempMK *)(T->ParamTrans))->u           = u;
  ((paramTempMK *)(T->ParamTrans))->eta         = eta;
  ((paramTempMK *)(T->ParamTrans))->mu          = mu;
  ((paramTempMK *)(T->ParamTrans))->random      = alea;
  ((paramTempMK *)(T->ParamTrans))->Optimize    = opt;
  ((paramTempMK *)(T->ParamTrans))->DispProgress = DispProgress;
  ((paramTempMK *)(T->ParamTrans))->limitv      = limitv;
  CopyBV(&(((paramTempMK *)(T->ParamTrans))->b),          b);
  CopyBV(&(((paramTempMK *)(T->ParamTrans))->b_original), b);
  CopyBV(&(((paramTempMK *)(T->ParamTrans))->BestB),      b);
  CopyBV(&(((paramTempMK *)(T->ParamTrans))->c),          c);
  CopyBV(&(((paramTempMK *)(T->ParamTrans))->BestC),      c);
  CopyBV(&(((paramTempMK *)(T->ParamTrans))->c_original), c);
}

/* Copy all fields of the MK tempering Transformation T2 into T1. */
void CopyTemperMK(Transformation *T1, Transformation *T2) {
  if (((paramTempMK *)T1->ParamTrans)->b.n != 0) {
    FreeBV(&(((paramTempMK *)(T1->ParamTrans))->b));
    FreeBV(&(((paramTempMK *)(T1->ParamTrans))->b_original));
    FreeBV(&(((paramTempMK *)(T1->ParamTrans))->BestB));
    FreeBV(&(((paramTempMK *)(T1->ParamTrans))->c));
    FreeBV(&(((paramTempMK *)(T1->ParamTrans))->c_original));
    FreeBV(&(((paramTempMK *)(T1->ParamTrans))->BestC));
  }
  AllocBV(&(((paramTempMK *)(T1->ParamTrans))->b),          T2->w);
  AllocBV(&(((paramTempMK *)(T1->ParamTrans))->b_original), T2->w);
  AllocBV(&(((paramTempMK *)(T1->ParamTrans))->BestB),      T2->w);
  AllocBV(&(((paramTempMK *)(T1->ParamTrans))->c),          T2->w);
  AllocBV(&(((paramTempMK *)(T1->ParamTrans))->c_original), T2->w);
  AllocBV(&(((paramTempMK *)(T1->ParamTrans))->BestC),      T2->w);

  ((paramTempMK *)T1->ParamTrans)->type         = ((paramTempMK *)T2->ParamTrans)->type;
  ((paramTempMK *)T1->ParamTrans)->l            = ((paramTempMK *)T2->ParamTrans)->l;
  ((paramTempMK *)T1->ParamTrans)->u            = ((paramTempMK *)T2->ParamTrans)->u;
  ((paramTempMK *)T1->ParamTrans)->eta          = ((paramTempMK *)T2->ParamTrans)->eta;
  ((paramTempMK *)T1->ParamTrans)->mu           = ((paramTempMK *)T2->ParamTrans)->mu;
  ((paramTempMK *)T1->ParamTrans)->random       = ((paramTempMK *)T2->ParamTrans)->random;
  ((paramTempMK *)T1->ParamTrans)->Optimize     = ((paramTempMK *)T2->ParamTrans)->Optimize;
  ((paramTempMK *)T1->ParamTrans)->DispProgress = ((paramTempMK *)T2->ParamTrans)->DispProgress;
  ((paramTempMK *)T1->ParamTrans)->limitv       = ((paramTempMK *)T2->ParamTrans)->limitv;

  CopyBV(&(((paramTempMK *)T1->ParamTrans)->b),          &(((paramTempMK *)T2->ParamTrans)->b));
  CopyBV(&(((paramTempMK *)T1->ParamTrans)->b_original), &(((paramTempMK *)T2->ParamTrans)->b_original));
  CopyBV(&(((paramTempMK *)T1->ParamTrans)->BestB),      &(((paramTempMK *)T2->ParamTrans)->BestB));
  CopyBV(&(((paramTempMK *)T1->ParamTrans)->c),          &(((paramTempMK *)T2->ParamTrans)->c));
  CopyBV(&(((paramTempMK *)T1->ParamTrans)->c_original), &(((paramTempMK *)T2->ParamTrans)->c_original));
  CopyBV(&(((paramTempMK *)T1->ParamTrans)->BestC),      &(((paramTempMK *)T2->ParamTrans)->BestC));

  T1->w               = T2->w;
  T1->w_original      = T2->w_original;
  T1->DispLinTrans    = T2->DispLinTrans;
  T1->Trans           = T2->Trans;
  T1->ChangeParamTrans = T2->ChangeParamTrans;
  T1->GiveTransID     = T2->GiveTransID;
  T1->AllocTrans      = T2->AllocTrans;
  T1->CopyTrans       = T2->CopyTrans;
  T1->InverseTrans    = T2->InverseTrans;
}

/* Allocate a new Transformation structure for an MK tempering with all bit vectors
   initialized to zero size. */
Transformation *AllocTemperMK(void) {
  Transformation *T;
  T = (Transformation *) malloc(sizeof(Transformation));
  T->ParamTrans = (paramTempMK *) malloc(sizeof(paramTempMK));
  ((paramTempMK *)T->ParamTrans)->b.n          = 0;
  ((paramTempMK *)T->ParamTrans)->b_original.n = 0;
  ((paramTempMK *)T->ParamTrans)->BestB.n      = 0;
  ((paramTempMK *)T->ParamTrans)->c.n          = 0;
  ((paramTempMK *)T->ParamTrans)->c_original.n = 0;
  ((paramTempMK *)T->ParamTrans)->BestC.n      = 0;
  return T;
}

/* Free all memory associated with an MK tempering Transformation previously allocated
   by AllocTemperMK(). */
void FreeTemperMK(Transformation *T) {
  FreeBV(&(((paramTempMK *)T->ParamTrans)->b));
  FreeBV(&(((paramTempMK *)T->ParamTrans)->b_original));
  FreeBV(&(((paramTempMK *)T->ParamTrans)->BestB));
  FreeBV(&(((paramTempMK *)T->ParamTrans)->c));
  FreeBV(&(((paramTempMK *)T->ParamTrans)->c_original));
  FreeBV(&(((paramTempMK *)T->ParamTrans)->BestC));
  free(T->ParamTrans);
  free(T);
}

/* Compute the inverse of the type-I MK tempering stored in T and apply it in-place
   to the bit vector RES.  Only works for type I. */
void InverseTemperMK(BitVect *RES, Transformation *T) {
  BitVect Inverse1, Inverse2, mask, Temp, RES2;
  int i;

  if (((paramTempMK *)(T->ParamTrans))->type == 2) {
    printf("InverseMK: Only implemented for type I!\n");
    exit(1);
  }
  AllocBV(&RES2,     RES->n * WL);
  AllocBV(&Inverse1, RES->n * WL);
  AllocBV(&Inverse2, RES->n * WL);
  AllocBV(&mask,     RES->n * WL);
  AllocBV(&Temp,     RES->n * WL);

  CopyBV(&RES2, RES);

  ANDBVMask(RES, RES, T->w);

  CopyBV(&mask, &CC);
  CopyBV(&Inverse1, RES);
  i = 1;
  do {
    BVLShift(&Temp, RES, i * MU);
    ANDBV(&Temp, &Temp, &mask);
    XORBV(&Inverse1, &Inverse1, &Temp);
    BVLShift(&Temp, &CC, i * MU);
    ANDBV(&mask, &mask, &Temp);
    i++;
  } while (i * MU <= T->w);

  CopyBV(&mask, &BB);
  CopyBV(&Inverse2, &Inverse1);
  i = 1;
  do {
    BVLShift(&Temp, &Inverse1, i * ETA);
    ANDBV(&Temp, &Temp, &mask);
    XORBV(&Inverse2, &Inverse2, &Temp);
    BVLShift(&Temp, &BB, i * ETA);
    ANDBV(&mask, &mask, &Temp);
    i++;
  } while (i * ETA < T->w);
  ANDBVInvMask(&RES2, &RES2, T->w);
  XORBV(RES, &RES2, &Inverse2);
  FreeBV(&RES2);
  FreeBV(&Inverse1);
  FreeBV(&Inverse2);
  FreeBV(&mask);
  FreeBV(&Temp);
}

/* Apply the MK tempering transformation stored in T to the first w bits of RES,
   implementing either type I or type II of the Matsumoto-Kurita tempering. */
void TemperingMK(BitVect *RES, Transformation *T) {
  BitVect g_state_Tempering, Temp;
  uint32_t x;

  if (T->w != WL) {
    AllocBV(&g_state_Tempering, RES->n * WL);
    AllocBV(&Temp, RES->n * WL);
    CopyBV(&g_state_Tempering, RES);
    ANDBVMask(&g_state_Tempering, &g_state_Tempering, T->w);
    ANDBVInvMask(RES, RES, T->w);

    if (TYPE == 2) {
      BVRShift(&Temp, &g_state_Tempering, UU);
      XORBV(&g_state_Tempering, &g_state_Tempering, &Temp);
      ANDBVMask(&g_state_Tempering, &g_state_Tempering, T->w);
    }
    BVLShift(&Temp, &g_state_Tempering, ETA);
    ANDBV(&Temp, &Temp, &BB);
    XORBV(&g_state_Tempering, &g_state_Tempering, &Temp);

    BVLShift(&Temp, &g_state_Tempering, MU);
    ANDBV(&Temp, &Temp, &CC);
    XORBV(&g_state_Tempering, &g_state_Tempering, &Temp);

    if (TYPE == 2) {
      BVRShift(&Temp, &g_state_Tempering, LL);
      XORBV(&g_state_Tempering, &g_state_Tempering, &Temp);
      ANDBVMask(&g_state_Tempering, &g_state_Tempering, T->w);
    }

    XORBV(RES, &g_state_Tempering, RES);
    FreeBV(&g_state_Tempering);
    FreeBV(&Temp);
  } else {
    x = RES->vect[0];
    if (TYPE == 2)
      x ^= x >> UU;
    x ^= (x << ETA) & BB.vect[0];
    x ^= (x << MU) & CC.vect[0];
    if (TYPE == 2)
      x ^= x >> LL;
    RES->vect[0] = x;
  }
}

/* Update the b and c bit vectors of the MK tempering T: either draw new random vectors
   (if random == -1) or reset to the original values stored in b_original and c_original. */
void UpdateTemperMK(Transformation *T) {
  int NbBVTempering, k;

  NbBVTempering = (T->w - 1) / WL + 1;
  if (((paramTempMK *)(T->ParamTrans))->random == -1) {
    for (k = 0; k < NbBVTempering; k++) {
      ((paramTempMK *)(T->ParamTrans))->BestB.vect[k] =
        ((paramTempMK *)(T->ParamTrans))->b.vect[k] = (uint32_t) MRG32k3a();
      ((paramTempMK *)(T->ParamTrans))->BestC.vect[k] =
        ((paramTempMK *)(T->ParamTrans))->c.vect[k] = (uint32_t) MRG32k3a();
    }
  } else {
    for (k = 0; k < NbBVTempering; k++) {
      ((paramTempMK *)(T->ParamTrans))->BestB.vect[k] =
        ((paramTempMK *)(T->ParamTrans))->b.vect[k] =
        ((paramTempMK *)(T->ParamTrans))->b_original.vect[k];
      ((paramTempMK *)(T->ParamTrans))->BestC.vect[k] =
        ((paramTempMK *)(T->ParamTrans))->c.vect[k] =
        ((paramTempMK *)(T->ParamTrans))->c_original.vect[k];
    }
  }
  ANDBVMask(&(((paramTempMK *)(T->ParamTrans))->BestB), &(((paramTempMK *)(T->ParamTrans))->BestB), T->w);
  ANDBVMask(&(((paramTempMK *)(T->ParamTrans))->BestC), &(((paramTempMK *)(T->ParamTrans))->BestC), T->w);
  ANDBVMask(&(((paramTempMK *)(T->ParamTrans))->b),     &(((paramTempMK *)(T->ParamTrans))->b),     T->w);
  ANDBVMask(&(((paramTempMK *)(T->ParamTrans))->c),     &(((paramTempMK *)(T->ParamTrans))->c),     T->w);
}

/* Display the type and parameters of the MK tempering stored in T, including the
   current best b and c vectors in hexadecimal. */
void DispTemperMK(Transformation *T) {
  int i, NbBVTempering;

  if (((paramTempMK *)(T->ParamTrans))->type == 1) {
    printf("Matsumoto-Kurita Tempering (I)(%d,%d)(w=%d)",
           ((paramTempMK *)(T->ParamTrans))->eta,
           ((paramTempMK *)(T->ParamTrans))->mu,
           T->w);
  } else {
    printf("Matsumoto-Nishimura Tempering (II)(%d,%d,%d,%d)(w=%d)",
           ((paramTempMK *)(T->ParamTrans))->u,
           ((paramTempMK *)(T->ParamTrans))->eta,
           ((paramTempMK *)(T->ParamTrans))->mu,
           ((paramTempMK *)(T->ParamTrans))->l,
           T->w);
  }
  NbBVTempering = (T->w - 1) / WL + 1;
  printf("  b = ");
  for (i = 0; i < NbBVTempering; i++)
    printf("%08x ", ((paramTempMK *)(T->ParamTrans))->BestB.vect[i]);

  printf("  c = ");
  for (i = 0; i < NbBVTempering; i++)
    printf("%08x ", ((paramTempMK *)(T->ParamTrans))->BestC.vect[i]);
  printf("\n");
}

/* Macros used to access tempering parameters for component j in OptimizeTemper(). */
#define Temp_w(j)        C->Components[j].Trans[C->Components[j].NbTrans - 1]->w
#define Temp_type(j)     ((paramTempMK *)(C->Components[j].Trans[C->Components[j].NbTrans - 1]->ParamTrans))->type
#define Temp_u(j)        ((paramTempMK *)(C->Components[j].Trans[C->Components[j].NbTrans - 1]->ParamTrans))->u
#define Temp_mu(j)       ((paramTempMK *)(C->Components[j].Trans[C->Components[j].NbTrans - 1]->ParamTrans))->mu
#define Temp_eta(j)      ((paramTempMK *)(C->Components[j].Trans[C->Components[j].NbTrans - 1]->ParamTrans))->eta
#define Temp_l(j)        ((paramTempMK *)(C->Components[j].Trans[C->Components[j].NbTrans - 1]->ParamTrans))->u
#define Temp_b(j)        ((paramTempMK *)(C->Components[j].Trans[C->Components[j].NbTrans - 1]->ParamTrans))->b
#define Temp_c(j)        ((paramTempMK *)(C->Components[j].Trans[C->Components[j].NbTrans - 1]->ParamTrans))->c
#define Temp_Bb(j)       ((paramTempMK *)(C->Components[j].Trans[C->Components[j].NbTrans - 1]->ParamTrans))->BestB
#define Temp_Bc(j)       ((paramTempMK *)(C->Components[j].Trans[C->Components[j].NbTrans - 1]->ParamTrans))->BestC
#define gen_L(j)         C->Components[j].Gen[C->Components[j].CurrentGen]->L
#define is_temper(j)     (C->Components[j].Trans != NULL) && \
                         (C->Components[j].Trans[C->Components[j].NbTrans - 1]->GiveTransID() == TEMPMKID)
#define is_temperopt(j)  ((paramTempMK *)(C->Components[j].Trans[C->Components[j].NbTrans - 1]->ParamTrans))->Optimize
#define DISPPROGRESS(j)  ((paramTempMK *)(C->Components[j].Trans[C->Components[j].NbTrans - 1]->ParamTrans))->DispProgress
#define LIMITV(j)        ((paramTempMK *)(C->Components[j].Trans[C->Components[j].NbTrans - 1]->ParamTrans))->limitv
#define METHOD           (params)->method

/* Static state variables used by OptimizeTemper(). */
static int *cvmt, *cvmtl, *bv, *bvml, previousv;
static int max_v,  /* bit position reached furthest in the optimization */
    essais = 0;    /* total number of random perturbation attempts across all recursions */
static BitVect *mask_c, /* mask keeping bits of c that maintain equidistribution up to bit v-1 */
               *mask_b; /* same for b */
static int *Bestse;     /* Bestse[v]: best sum of gaps when optimizing the v-th bit */
static BitVect AleaBV, TemporaireBV, M, *BVPolys;
static Base *StackBase, CurrentBase;
static ZENRing F2, polyring;
static ZENElt inv;
static int limite_v, affprog;

/* Recursively optimize the b and c vectors of all MK tempering transformations in C
   to minimize the equidistribution gap at bit v, up to resolution limite_v. */
static void OptimizeTemper(Combinaison *C, paramMecf *params, int v, int t,
                           int limite_v, int affprog) {
  int t1, v1, W, i, j, x, length,
      Lim; /* number of random perturbation attempts with the current masks */

  if (v < limite_v / 2)
    Lim = 5;
  else if (v == max_v)
    Lim = 2;
  else
    Lim = 2;

  for (x = 0; (x < Lim) && (essais < 400); x++) {
    essais++;
    for (j = 0; j < C->J; j++) {
      if (is_temper(j))
        if (is_temperopt(j)) {
          invmask(&mask_c[j], v - 1);
          W = intmin(gen_L(j), Temp_w(j)); /* number of bits used for tempering */
          CopyBV(&mask_b[j], &mask_c[j]);
          for (i = Temp_mu(j) + 1; (i < v + Temp_mu(j)) && (i < W); i++)
            PutBitBV(&mask_b[j], i - 1, 0);
          if (v - Temp_mu(j) > 0)
            cvmt[j] = ValBitBV(&Temp_c(j), v - Temp_mu(j) - 1);
          bv[j] = ValBitBV(&Temp_b(j), v - 1);
          if (Temp_type(j) == 2) {
            if (v - Temp_l(j) - Temp_mu(j) > 0)
              cvmtl[j] = ValBitBV(&Temp_c(j), v - Temp_mu(j) - Temp_l(j) - 1);
            bvml[j] = ValBitBV(&Temp_b(j), v - Temp_l(j) - 1);
          }
          RandVect(&AleaBV);
          ANDBV(&TemporaireBV, &AleaBV, &mask_b[j]);
          XORBVSelf(&Temp_b(j), &TemporaireBV);
          RandVect(&AleaBV);
          ANDBV(&TemporaireBV, &AleaBV, &mask_c[j]);
          XORBVSelf(&Temp_c(j), &TemporaireBV);
          ANDBVMask(&Temp_b(j), &Temp_b(j), W);
          ANDBVMask(&Temp_c(j), &Temp_c(j), W);
          if (v > W - Temp_mu(j))
            PutBitBV(&Temp_c(j), v - 1, 0);
          if (v > W - Temp_eta(j))
            PutBitBV(&Temp_b(j), v - 1, 0);
          if ((cvmt[j] == 1) && (v - Temp_mu(j) > 0))
            PutBitBV(&Temp_b(j), v - 1, bv[j]);
          if ((v + Temp_mu(j) > W - Temp_eta(j)) || (ValBitBV(&Temp_c(j), v - 1) == 0))
            if (v + Temp_mu(j) <= W)
              PutBitBV(&Temp_b(j), v + Temp_mu(j) - 1, 0);
        }
    }
    if (METHOD == TestME) {
      PrepareMat(C, &Mat, C->k_g);
      CopyMat(&mat, &Mat, C->k_g, t);
      params->ecart[v] = t - DimensionEquid(&mat, C->k_g, v, C->smax);
    } else if (METHOD == TestMELat) {
      FindPolys(C, &M, BVPolys, v, v);
      if (v == 1) {
        ZENEltFree(inv, polyring);
        ComputeInvg1(&BVPolys[0], F2, &inv, &polyring);
        DualBase(&StackBase[0], BVPolys, &M, 1);
        CopyBase(&CurrentBase, &StackBase[0]);
      } else if (v > previousv) {
        FindPolysTimesInv(BVPolys, F2, &inv, &polyring, v, v);
        CopyBase(&StackBase[v - 2], &CurrentBase);
        DualBaseIncrease(&CurrentBase, BVPolys);
      } else {
        FindPolysTimesInv(BVPolys, F2, &inv, &polyring, v, v);
        CopyBase(&CurrentBase, &StackBase[v - 2]);
        DualBaseIncrease(&CurrentBase, BVPolys);
      }
      previousv = v;
      length = Lenstra(&CurrentBase, v);
      params->ecart[v] = C->k_g / v - intmin(length, C->k_g / v);
    }

    params->se = 0;
    for (i = 1; i <= v; i++)
      params->se += params->ecart[i];

    if (params->se < Bestse[v] && v >= max_v) {
      essais = 0;
      for (j = 0; j < C->J; j++) {
        if (is_temper(j))
          if (is_temperopt(j)) {
            CopyBV(&Temp_Bb(j), &Temp_b(j));
            CopyBV(&Temp_Bc(j), &Temp_c(j));
          }
      }
      Bestse[v] = params->se;
      max_v = v;
      if (affprog) {
        fprintf(stdout, "v=%d gap sum=%d tries=%d %d\n",
                v, params->se, essais, params->delta[limite_v]);
      }
    } else if (affprog) {
      fprintf(stdout, "v=%d gap sum=%d tries=%d  BestB=%08x  BestC=%08x\n",
              v, params->se, essais, Temp_Bb(0).vect[0], Temp_Bc(0).vect[0]);
    }
    fflush(stdout);
    if ((params->ecart[v] <= params->delta[v]) && (params->se <= params->mse)) {
      if (v >= limite_v) {
        break;
      } else {
        v1 = v + 1;
        t1 = (C->k_g / v1);
        OptimizeTemper(C, params, v1, t1, limite_v, affprog);
      }
    }
    if ((params->ecart[limite_v] <= params->delta[limite_v]) &&
        (params->se <= params->mse) &&
        (params->ecart[limite_v] >= 0)) {
      break;
    }
  }
}

/* Same as TestME() but also optimizes the b and c vectors of all MK tempering
   transformations marked with Optimize=TRUE, using the algorithm described in
   Matsumoto & Kurita (ACM TOMACS 1994). */
void TestMETemperMK(Combinaison *C, paramMecf *params) {
  int *coeff, j, t, l;

  if (METHOD == TestMELat) {
    BVPolys = (BitVect *) malloc(params->L * sizeof(BitVect));
    for (j = 0; j < params->L; j++)
      AllocBV(&(BVPolys[j]), C->k_g + 1);
    AllocBV(&M, C->k_g + 1);
    StackBase = (Base *) malloc(params->L * sizeof(Base));
    for (j = 0; j < params->L; j++)
      AllocBase(&StackBase[j], params->L, C->k_g);
    AllocBase(&CurrentBase, params->L, C->k_g);
    coeff = (int *) malloc((C->k_g + 1) * sizeof(int));
    polycharComb(C, coeff, &M);
    ConstructZENF2(&F2);
    ConstructZENRing(&M, C->k_g, &polyring);
    ZENEltAlloc(inv, polyring);
    free(coeff);
  } else if (METHOD == TestME) {
    if (firsttime) {
      AllocMat(&Mat, C->k_g, C->L, intmin(C->k_g, C->smax));
      AllocMat(&mat, C->k_g, C->L, intmin(C->k_g, C->smax));
      firsttime = FALSE;
    } else {
      if ((Mat.nblignes < C->k_g) || (Mat.t < C->smax)) {
        if (Mat.nblignes != 0) {
          FreeMat(&Mat);
          FreeMat(&mat);
        }
        AllocMat(&Mat, C->k_g, C->L, intmin(C->k_g, C->smax));
        AllocMat(&mat, C->k_g, C->L, intmin(C->k_g, C->smax));
      }
    }
  }
  SetPsi12(C, params);
  Bestse   = (int *) malloc((C->L + 1) * sizeof(int));
  cvmt     = (int *) malloc(C->J * sizeof(int));
  cvmtl    = (int *) malloc(C->J * sizeof(int));
  bv       = (int *) malloc(C->J * sizeof(int));
  bvml     = (int *) malloc(C->J * sizeof(int));
  mask_c   = (BitVect *) malloc(C->J * sizeof(BitVect));
  mask_b   = (BitVect *) malloc(C->J * sizeof(BitVect));
  for (j = 0; j < C->J; j++) {
    AllocBV(&mask_c[j], C->L);
    AllocBV(&mask_b[j], C->L);
  }
  AllocBV(&AleaBV, C->L);
  AllocBV(&TemporaireBV, C->L);

  limite_v = INT_MAX;
  for (j = 0; j < C->J; j++)
    if (limite_v > LIMITV(j))
      limite_v = LIMITV(j);

  affprog = FALSE;
  for (j = 0; j < C->J; j++) {
    if (limite_v > C->Components[j].Gen[C->Components[j].CurrentGen]->L)
      limite_v = C->Components[j].Gen[C->Components[j].CurrentGen]->L;
    affprog |= DISPPROGRESS(j);
  }

  for (j = 0; j < C->L; j++)
    Bestse[j] = INT_MAX;
  max_v = 0;

  for (j = 0; j <= limite_v; j++) /* Reset the ecart array. */
    params->ecart[j] = -1;
  params->ecart[1] = 0;
  params->se = 0;
  l = 1;
  t = C->k_g; /* t associated with l=1 */
  essais = 0;
  OptimizeTemper(C, params, l, t, limite_v, affprog);

  for (j = 0; j < C->J; j++) {
    if (is_temper(j))
      if (is_temperopt(j)) {
        CopyBV(&Temp_b(j), &Temp_Bb(j));
        CopyBV(&Temp_c(j), &Temp_Bc(j));
      }
  }
  free(Bestse);
  free(cvmt);
  free(cvmtl);
  free(bv);
  free(bvml);
  for (j = 0; j < C->J; j++) {
    FreeBV(&mask_c[j]);
    FreeBV(&mask_b[j]);
  }
  free(mask_c);
  free(mask_b);
  FreeBV(&AleaBV);
  FreeBV(&TemporaireBV);
  if (METHOD == TestMELat) {
    for (j = 0; j < params->L; j++)
      FreeBV(&(BVPolys[j]));
    free(BVPolys);
    FreeBV(&M);
    ZENEltFree(inv, polyring);
    for (j = 0; j < params->L; j++)
      FreeBase(&StackBase[j]);
    free(StackBase);
    FreeBase(&CurrentBase);
    ZENRingClose(F2);
    ZENRingClose(polyring);
  }
  TestEquid(C, params);
}
