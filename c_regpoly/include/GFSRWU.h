
#include "vectorsF2.h"
#include "combinaisons.h"

 
#ifndef GFSR_H
#define GFSR_H
 



typedef struct {
  int p;
  int q;
  int r;
  int s;
  int d;
  BitVect Buffer;
  int ix;
} paramGFSR;



char* GFSRName ( void
               );



void InitParamGFSR ( Generateur *Gen,
                     int p,
                     int q,
                     int r,
                     int s,
                     int d,
                     int L
                   );



void InitGFSR ( Generateur *Gen,
                BitVect *init,
                BitVect *retour
              );



void CopyGFSR ( Generateur *Gen1,
                Generateur *Gen2
              );



Generateur* AllocGFSR ( int k
                      );



void FreeGFSR ( Generateur *G
              );



void GFSR ( Generateur *Gen,
            BitVect *g_etat
          );



void DispGFSR ( Generateur *Gen
              );



void ReadDataGFSR ( Component *E,
                    char *nomfich,
                    boolean meme,
                    int L
                  );



void ProduceParamsGFSR ( char *filename,
                         int nb,
                         int k
                       );



void CharGFSR ( Generateur *Gen,
                int coeff[],
                BitVect *BVPoly
              );

 
#endif
 

