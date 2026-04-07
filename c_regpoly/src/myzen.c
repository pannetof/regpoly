/*
 * myzen.c — Wrappers for constructing ZEN finite field structures
 *           used by the ZEN library.
 */

#include <stdio.h>
#include <stdlib.h>
#include "myzen.h"

/* Constructs the base field GF(2) as a ZENRing and stores it in *F2. */
void ConstructZENF2(ZENRing *F2)
{
  BigNum two;
  BigNumLength twol;

  ZBNReadFromString(&two, &twol, "2", 10);
  if (ZENBaseRingAlloc(*F2, two, twol)) {
    printf("Error: ZENBaseRingAlloc failed for GF(2)\n");
    exit(1);
  }
  ZBNF(two);
}

/* Constructs the extension field GF(2^degree) defined by the primitive
   polynomial encoded in the bit vector M, and stores it in *F2k. */
void ConstructZENRing(BitVect *M, int degree, ZENRing *F2k)
{
  ZENElt ONE;
  ZENPoly ZenM;
  ZENRing F2;
  int k;

  ConstructZENF2(&F2);
  ZENEltAlloc(ONE, F2);
  ZENEltSetToOne(ONE, F2);
  ZENPolyAlloc(ZenM, degree, F2);
  ZENPolySetToXi(ZenM, degree, F2);
  for (k = 0; k < degree; k++)
    if (ValBitBV(M, k) == 1)
      ZENPolySetCoeff(ZenM, k, ONE, F2);
  if (ZENExtRingAlloc(*F2k, ZenM, F2)) {
    printf("Error: ZENExtRingAlloc failed\n");
    exit(1);
  }
  ZENPolyFree(ZenM, F2);
  ZENEltFree(ONE, F2);
  ZENRingClose(F2);
}
