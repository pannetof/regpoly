
#include "vectorsF2.h"
#include "combinaisons.h"

 
#ifndef MARSA_H
#define MARSA_H
 



typedef struct {
  int a,b,c;
  int w;
} paramMarsa;



char* MarsaName ( void
                );



void InitParamMarsa ( Generateur *Gen,
                      int w,
                      int a,
                      int b,
                      int c,
                      int L
                    );



void InitMarsa ( Generateur *Gen,
                 BitVect *init,
                 BitVect *retour
               );



void CopyMarsa ( Generateur *Gen1,
                 Generateur *Gen2
               );



Generateur* AllocMarsa ( int k
                       );



void FreeMarsa ( Generateur *G
               );



void Marsa ( Generateur *Gen,
             BitVect *retour
           );



void DispMarsa ( Generateur *Gen
               );



void ReadDataMarsa ( Component *E,
                     char *nomfich,
                     boolean meme,
                     int L
                   );

 
#endif
 

