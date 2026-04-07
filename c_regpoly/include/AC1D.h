
#include "vectorsF2.h"
#include "combinaisons.h"

 
#ifndef AC1D_H
#define AC1D_H
 



typedef struct {
  int N;
  BitVect *MaskShift;
  BitVect *MatrixR;
  int taillevoisinage;
} paramAC1D;



char* AC1DName ( void
               );



void InitParamAC1D ( Generateur *Gen,
                     int N,
                     BitVect *MaskShift,
                     BitVect *MatrixR,
                     int taillevoisinage,
                     int L
                   );



void InitAC1D ( Generateur *Gen,
                BitVect *init,
                BitVect *retour
              );



void CopyAC1D ( Generateur *Gen1,
                Generateur *Gen2
              );



Generateur* AllocAC1D ( int k
                      );



void FreeAC1D ( Generateur *Gen
              );



void AC1D ( Generateur *Gen,
            BitVect *gen_etat
          );



void DispAC1D ( Generateur *Gen
             );



void ReadDataAC1D ( Component *E,
                    char *nomfich,
                    boolean meme,
                    int L
                  );



void ProduceParamsAC1D ( char *filename,
                         int nb,
                         char F,
                         int n,
                         int m,
                         int neighborhood
                       );

 
#endif
 

