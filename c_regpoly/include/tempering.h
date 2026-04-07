
#include "vectorsF2.h"
#include "combinaisons.h"


 
#ifndef TEMPERING_H
#define TEMPERING_H
 


#define TRANS(Comb,no)  (Comb)->Components[no].Trans 
#define TRANSNB(Comb, no)  (Comb)->Components[no].NbTrans 
#define TRANSFORME(Comb, no, stateBV, outputBV)   Transform(TRANS(Comb,no),TRANSNB(Comb,no),stateBV,outputBV) 



typedef struct TransType{
  void* ParamTrans;
  void (*DispLinTrans)(struct TransType *T);
  void (*Trans)(BitVect *A, struct TransType *T);
  void (*InverseTrans)(BitVect *A, struct TransType *T);
  void (*ChangeParamTrans)(struct TransType *T);
  int (*GiveTransID)(void);
  struct TransType* (*AllocTrans)(void);
  void (*CopyTrans)(struct TransType *S1, struct TransType *S2);
  int w;
  int w_original;
} Transformation;



void AllocTransInComponent ( Component *E,
                             int nbtrans
                           );



void FreeTransInComponent ( Component *E
                          );



void AddTransInComponent ( Component *E,
                           Transformation *T
                         );



void DispTrans ( Component *E
               );



void Transform ( Transformation **T,
                 int NbTrans,
                 BitVect *A,
                 BitVect *Out
               );



void TransformInverse ( Transformation **T,
                        int NbTrans,
                        BitVect *A,
                        BitVect *Out
                      );



void UpdateAllTrans ( Combinaison *C
                    );



void UpdateTrans ( Component *E
                 );

 
#endif
 

