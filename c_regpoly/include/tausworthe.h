
#include "vectorsF2.h"
#include "combinaisons.h"

 
#ifndef TAUSWORTHE_H
#define TAUSWORTHE_H
 



#define PRIMITIVE 1
#define IRREDUCIBLE 0




#define MAXNBCOEFF 100



typedef struct {
  boolean quicktaus;
  int Q[MAXNBCOEFF];
  int s;
  int NbCoeff;   int gen_kms;   
  int StateL;
} paramTaus;



char* TausName ( void
               );



void InitParamTaus ( Generateur *Gen,
                     int k,
                     int Q[MAXNBCOEFF],
                     int s,
                     int NbCoeff,
                     boolean quicktaus,
                     int L
                   );




void CopyTaus ( Generateur *Gen1,
                Generateur *Gen2
              );



Generateur* AllocTaus ( int k
                      );



void FreeTaus ( Generateur *G
              );



void InitQuickTaus ( Generateur *Gen,
                     BitVect *init,
                     BitVect *retour
                   );



void QuickTaus ( Generateur *Gen,
                 BitVect *retour
               );



void InitGeneralTaus ( Generateur *Gen,
                       BitVect *init,
                       BitVect *retour
                     );



void GeneralTaus ( Generateur *Gen,
                   BitVect *retour
                 );



void DispTaus ( Generateur *Gen
              );



void ReadDataTaus ( Component *E,
                    char *nomfich,
                    boolean same,
                    int L
                  );



void ProduceParamsTaus ( char *filename,
                         int nb,
                         int k,
                         int nbcoeff,
                         int prim_or_irr,
                         boolean quicktaus
                       );



void CharTaus ( Generateur *Gen,
                int coeff[],
                BitVect *BVPoly
              );

 
#endif
 

