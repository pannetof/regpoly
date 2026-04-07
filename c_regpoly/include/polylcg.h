
#include "vectorsF2.h"
#include "combinaisons.h"

 
#ifndef POLYLCG_H
#define POLYLCG_H
 



typedef struct {
  BitVect Polynome;
 
  BitVect gen_masque;
 
} paramPolyLCG;



char* PolyLCGName ( void
               );



void InitParamPolyLCG ( Generateur *Gen,
                        int k,
                        BitVect *Polynome,
                        int L
                      );



void InitPolyLCG ( Generateur *Gen,
                   BitVect *init,
                   BitVect *retour
                 );



void CopyPolyLCG ( Generateur *Gen1,
                   Generateur *Gen2
                 );



Generateur* AllocPolyLCG ( int k
                         );



void FreePolyLCG ( Generateur *G
                 );



void PolyLCG ( Generateur *Gen,
               BitVect *retour
             );



void DispPolyLCG ( Generateur *Gen
                 );



void ReadDataPolyLCG ( Component *E,
                       char *nomfich,
                       boolean meme,
                       int L
                     );



void ProduceParamsPolyLCG ( char *filename,
                            int nb,
                            int k
                          );



void CharPolyLCG ( Generateur *Gen,
                   int coeff[],
                   BitVect *BVPoly
                 );

 
#endif
 

