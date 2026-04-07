
#include "vectorsF2.h"
#include "combinaisons.h"
#include "tempering.h"

 
#ifndef PERMUT_H
#define PERMUT_H
 


#define PERMID 3
#define PERMNAME "Permutation"



typedef struct {
  int p;
  int p_original;
  int q;
  int q_original;
} parampermut;



int GivePermutID ( void
                 );



void InitPermut ( Transformation *T,
                  int w,
                  int p,
                  int q
                );



void CopyPermut ( Transformation *T1,
                  Transformation *T2
                 );



Transformation* AllocPermut ( void
                            );



void FreePermut ( Transformation *T
                );



void UpdatePermut ( Transformation *T
                  );



void Permut ( BitVect *A,
              Transformation *T
            );



void InversePermut ( BitVect *A,
                     Transformation *T
                   );



void DispPermut ( Transformation *T
               );

 
#endif
 

