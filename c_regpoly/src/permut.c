/*
 * permut.c — Bit permutation tempering transformation for REGPOLY generators.
 *
 * Two modes are supported: a linear bit permutation i -> (i*p + q) mod w,
 * and a cyclic rotation (when p == 10) of the w low-order bits by q positions.
 */

#include <stdio.h>
#include <stdlib.h>
#include "permut.h"
#include "regpoly.h"

#define P ((parampermut *)(T->ParamTrans))->p
#define Q ((parampermut *)(T->ParamTrans))->q

/* Returns the identifier for the permutation transformation type. */
int GivePermutID(void)
{
  return PERMID;
}

/* Prints the parameters (p, q) of a permutation transformation T. */
void DispPermut(Transformation *T)
{
  printf("Permutation(%d,%d)\n",
         ((parampermut *)(T->ParamTrans))->p,
         ((parampermut *)(T->ParamTrans))->q);
}

/* Applies the permutation transformation to the w low-order bits of A.
   If p == 10, performs a cyclic rotation by Q bits; otherwise applies
   the linear bit permutation i -> (i*p + q) mod w. */
void Permut(BitVect *A, Transformation *T)
{
  BitVect Temp, Temp2;
  int i, u = 0;

  AllocBV(&Temp, T->w);
  AllocBV(&Temp2, T->w);
  PutBVToZero(&Temp);
  if (P == 10) {
    ANDBVMask(&Temp, A, T->w);
    BVLShift(&Temp, A, Q);
    BVRShift(&Temp2, A, T->w - Q);
    XORBV(&Temp, &Temp, &Temp2);
  } else {
    for (i = 0; i < T->w; i++) {
      u = (i * P + Q) % T->w;
      PutBitBV(&Temp, i, ValBitBV(A, u));
    }
  }
  ANDBVInvMask(A, A, T->w);
  XORBV(A, A, &Temp);
  FreeBV(&Temp);
  FreeBV(&Temp2);
}

/* Computes the inverse of the linear bit permutation transformation.
   Uses the modular inverse of p to reverse the mapping. */
void InversePermut(BitVect *A, Transformation *T)
{
  BitVect Temp;
  int Invp, i, u;

  AllocBV(&Temp, T->w);
  PutBVToZero(&Temp);
  Invp = InverseModN(P, T->w);
  for (i = 0; i < T->w; i++) {
    u = ((i + T->w - Q) * Invp) % T->w;
    PutBitBV(&Temp, i, ValBitBV(A, u));
  }
  ANDBVInvMask(A, A, T->w);
  XORBV(A, A, &Temp);
  FreeBV(&Temp);
}

/* Updates the permutation parameters p and q.
   If p_original (or q_original) is -1, a random valid value is chosen;
   otherwise the original value is used. */
void UpdatePermut(Transformation *T)
{
  int param1, param2;

  if (((parampermut *)(T->ParamTrans))->p_original == -1) {
    param1 = T->w;
    while (gcd(param1, T->w) != 1)
      param1 = MRG32k3a() % T->w;
  } else {
    param1 = ((parampermut *)(T->ParamTrans))->p_original;
  }

  if (((parampermut *)(T->ParamTrans))->q_original == -1)
    param2 = MRG32k3a() % T->w;
  else
    param2 = ((parampermut *)(T->ParamTrans))->q_original;

  ((parampermut *)(T->ParamTrans))->p = param1;
  ((parampermut *)(T->ParamTrans))->q = param2;
}

/* Initializes transformation T as a bit permutation with output width w,
   permutation multiplier p, and offset q. Exits if p and w are not coprime.
   Wires all function pointers for this type. */
void InitPermut(Transformation *T, int w, int p, int q)
{
  ((parampermut *)(T->ParamTrans))->p = p;
  ((parampermut *)(T->ParamTrans))->p_original = p;
  ((parampermut *)(T->ParamTrans))->q = q;
  ((parampermut *)(T->ParamTrans))->q_original = q;

  if (p != -1 && gcd(p, w) != 1) {
    printf("Invalid value for Permutation: gcd(p, w) != 1\n");
    exit(1);
  }
  T->w = w;
  T->w_original = w;
  T->DispLinTrans = DispPermut;
  T->Trans = Permut;
  T->ChangeParamTrans = UpdatePermut;
  T->AllocTrans = AllocPermut;
  T->CopyTrans = CopyPermut;
  T->GiveTransID = GivePermutID;
  T->InverseTrans = InversePermut;
}

/* Allocates and returns a new uninitialized permutation transformation. */
Transformation *AllocPermut(void)
{
  Transformation *T;
  T = (Transformation *) malloc(sizeof(Transformation));
  T->ParamTrans = (parampermut *) malloc(sizeof(parampermut));
  return T;
}

/* Frees all memory associated with the permutation transformation T. */
void FreePermut(Transformation *T)
{
  free(T->ParamTrans);
  free(T);
}

/* Copies all fields of permutation transformation T2 into T1. */
void CopyPermut(Transformation *T1, Transformation *T2)
{
  ((parampermut *)T1->ParamTrans)->p          = ((parampermut *)T2->ParamTrans)->p;
  ((parampermut *)T1->ParamTrans)->p_original = ((parampermut *)T2->ParamTrans)->p_original;
  ((parampermut *)T1->ParamTrans)->q          = ((parampermut *)T2->ParamTrans)->q;
  ((parampermut *)T1->ParamTrans)->q_original = ((parampermut *)T2->ParamTrans)->q_original;
  T1->w = T2->w;
  T1->w_original = T2->w_original;
  T1->DispLinTrans = T2->DispLinTrans;
  T1->Trans = T2->Trans;
  T1->ChangeParamTrans = T2->ChangeParamTrans;
  T1->GiveTransID = T2->GiveTransID;
  T1->AllocTrans = T2->AllocTrans;
  T1->CopyTrans = T2->CopyTrans;
  T1->InverseTrans = T2->InverseTrans;
}
