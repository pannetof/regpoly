
#include "vectorsF2.h"
#include "combinaisons.h"

 
#ifndef TGFSR_H
#define TGFSR_H
 



typedef struct {
  int w;
  int r;
  int m;
  BitVect a;
} paramTGFSR;



char* TGFSRName ( void
                );



void InitParamTGFSR ( Generateur *Gen,
                      int w,
                      int r,
                      int m,
                      BitVect *A,
                      int L
                    );


void CopyTGFSR ( Generateur *TG1,
                 Generateur *TG2
                );



Generateur* AllocTGFSR ( int k
                       );



void FreeTGFSR( Generateur *G
              );



void InitTGFSR ( Generateur *Gen,
                 BitVect *init,
                 BitVect *retour
               );



void  TGFSR ( Generateur *Gen,
              BitVect *RETOUR
            );



void DispTGFSR ( Generateur *Gen
              );



void ReadDataTGFSR ( Component *E,
                     char *nomfich,
                     boolean meme,
                     int L
                   );



void ProduceParamsTGFSR ( char *filename,
                          int nb,
                          int w,
                          int r
                        );



void CharTGFSR ( Generateur *Gen,
                 int *coeff,
                 BitVect *BVPoly
               );

 
#endif
 
