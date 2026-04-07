
#include "vectorsF2.h"
#include "combinaisons.h"
#include <stdint.h>


#ifndef MT_H
#define MT_H
 



typedef struct {
  int w;
  int r;
  int m;
  int p;
  uint32_t a;   int i;  uint32_t ll;  uint32_t uu;  uint32_t mw; 
} paramMT;



char* MTName ( void
             );



Generateur* AllocMT ( int k
                    );



void FreeMT ( Generateur *Gen
            );



void InitParamMT ( Generateur *Gen,
                   int w,
                   int r,
                   int m,
                   int p,
                   uint32_t a,
                   int L
                 );



void InitMT ( Generateur *Gen,
              BitVect *init,
              BitVect *retour
            );



void CopyMT ( Generateur *Gen1,
              Generateur *Gen2
            );



void MT ( Generateur *Gen,
          BitVect *retour
        );



void DispMT ( Generateur *Gen
            );



void ReadDataMT ( Component *E,
                  char *nomfich,
                  boolean same,
                  int L
                );



void ProduceParamsMT ( char *filename,
                       int nb,
                       int exp,
                       int w,
                       int m
                     );



void CharMT ( Generateur *Gen,
              int coeff[],
              BitVect *BVPoly
            );

 
#endif
 

