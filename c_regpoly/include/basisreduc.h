

#include <zen.h>
#include "vectorsF2.h"
#include "combinaisons.h"
#include "mecf.h"

 
#ifndef REDUC_H
#define REDUC_H
 



void TestMELat ( Combinaison *C,
                 paramMecf *params
               );


typedef struct {
  BitVect *coeffs;
  int deg;
  int indicemaxdeg;
} PolVect;



typedef struct{
  PolVect *vect ;
  int *permutations;
  int *invpermutations;
  int resolution;
  int maxresolution;
  int degmax;
} Base;



void AllocBase ( Base *B,
                 int resolution,
                 int degmax
               );



void FreeBase ( Base *B
              );



void CopyBase ( Base *B1,
                Base *B2
              );



void ComputeInvg1 ( BitVect *g1,
                    ZENRing F2,
                    ZENElt *Invg1,
                    ZENRing *F2k
		  );




void FindPolys ( Combinaison *C,
                 BitVect *M,
                 BitVect BVPolys[],
                 int startres,
                 int endres
               );



void FindPolysTimesInv ( BitVect BVPolys[],
                         ZENRing F2,
                         ZENElt *inv,
                         ZENRing *F2k,
                         int startres,
                         int endres
                       );



int Lenstra ( Base *B,
              int res
            );



void DualBase ( Base *B,
                BitVect *R0,
                BitVect *M,
                int res
              );



void DualBaseIncrease ( Base *B,
                        BitVect *R0
                      );




void EltToBitVect(BitVect *A, ZENElt *e, ZENRing polring, ZENRing F2);



void BitVectToElt(ZENElt *e, BitVect *A, ZENRing polring, ZENRing F2);



void TestMELatSpecial ( Combinaison *C, paramMecf *params, int rotative_shift);


#endif

