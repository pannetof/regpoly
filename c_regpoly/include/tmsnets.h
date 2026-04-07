
#include "vectorsF2.h"
#include "combinaisons.h"


 
#ifndef TMSNETS_H
#define TMSNETS_H
 



typedef struct {
  Matrix M;
  int smax;
  int mmax;
  int Lmax;
  int i;   int currentbc; 
  char name[50];
} paramTMS;



char* TMSName ( void
              );



Generateur* AllocTMS ( int k
                     );



void FreeTMS ( Generateur *G
             );



void InitParamTMS ( Generateur *Gen,
                    Matrix *m,
                    int mmax,
                    int smax,
                    int Lmax,
                    int L,
                    char name[50]
                  );



void InitTMS ( Generateur *Gen,
               BitVect *init,
               BitVect *retour
             );



void CopyTMS ( Generateur *Gen1,
               Generateur *Gen2
              );



void TMS ( Generateur *Gen,
           BitVect *retour
         );



void DispTMS ( Generateur *Gen
            );



void ReadDataTMS ( Component *E,
                   char *nomfich,
                   boolean meme,
                   int L
                 );

 
#endif
 

