
#include "vectorsF2.h"
#include "combinaisons.h"

 
#ifndef MATSUMOTO_H
#define MATSUMOTO_H
 



typedef struct {
   int type;
   int n;
   int m;
   int nbparamsint;
   int *paramsint;
   int nbparamsunsigned;
   unsigned int *paramsunsigned;
} paramMatsumoto;



char* MatsumotoName ( void
                    );



void InitParamMatsumoto ( Generateur *Gen,
                        int type,
                        int n,
                        int m,
                        int nbparams,
                        int *params,
                        int nbparamsunsigned,
                        unsigned int *paramsunsigned,
                        int L
                      );



void InitMatsumoto ( Generateur *Gen,
                   BitVect *init,
                   BitVect *retour
                 );



void CopyMatsumoto ( Generateur *Gen1,
                         Generateur *Gen2
                 );



Generateur* AllocMatsumoto ( int k
                         );



void FreeMatsumoto ( Generateur *G
                 );



void Matsumoto1 ( Generateur *Gen,
               BitVect *retour
             );
void Matsumoto2 ( Generateur *Gen,
               BitVect *retour
             );
void Matsumoto3 ( Generateur *Gen,
               BitVect *retour
             );




void DispMatsumoto ( Generateur *Gen
                 );



void ReadDataMatsumoto ( Component *E,
                       char *nomfich,
                       boolean meme,
                       int L
                     );



void ProduceParamsMatsumoto ( char *filename,
                            int type,
                            int nb,
                            int n,
                            int Mexp
                          );

 
#endif
 

