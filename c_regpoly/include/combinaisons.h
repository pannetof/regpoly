

#include <stdio.h>
#include <stdlib.h>

#include "vectorsF2.h"
 
#ifndef INCLUDETEMPERING_H
#define INCLUDETEMPERING_H 
#include "tempering.h"
#endif

 
#ifndef ETAGES_H
#define ETAGES_H
 
#ifndef USE_BOOLEAN
#define USE_BOOLEAN 
typedef int    boolean; 
#endif 

#define GEN(Comb,noComponent)  (((Comb)->Components[noComponent].Gen[(Comb)->Components[noComponent].CurrentGen])) 
#define ITERATION(generateur, outputBV)  (generateur)->Iteration(generateur, outputBV) 
#define INITGEN(generateur, initBV, outputBV)  (generateur)->InitGen(generateur,initBV,outputBV) 
#define DEGGEN(generateur)  (generateur)->k 



typedef struct GenType{
  void* ParamGen;
  void (*DispGen)(struct GenType *Gen);
  char* (*DispName)(void);
  void (*InitGen)(struct GenType *Gen, BitVect *init, BitVect *retour);
  void (*Iteration)(struct GenType *Gen, BitVect *retour);
  struct GenType* (*AllocGen)(int k);
  void (*CopyGen)(struct GenType *G1, struct GenType *G2);
  void (*PolyChar)(struct GenType *Gen, int coeff[], BitVect *BVPoly);
  BitVect GenState;
  int smax;
  int k;
  int Step;
  int L;
} Generateur;



typedef struct {
  Generateur **Gen;
  struct TransType **Trans;
  int CurrentGen;
  int NbTrans;
  int NbTransAlloc;
  int NbGen;
  int NbGenAlloc;
  boolean same;
} Component;



typedef struct{
  Component* Components;
  int J;
  int k_g;
  int smax;
  int L;
  int Lmax;
} Combinaison;



void AllocGen ( Generateur *Gen,
                int k
              );



void FreeGen ( Generateur *Gen
             );



void AddGenInComponent ( Component *E,
                         Generateur *Gen
                       );



void CopyGen ( Generateur *Gen1,
               Generateur *Gen2
             );



void AllocComponentsInCombinaison ( Combinaison *C,
                                    int J,
                                    int Lmax
                                  );



void FreeComponentsInCombinaison ( Combinaison *C
                                 );



void AllocGensInComponent ( Component *E,
                            int nbgen,
                            boolean same
                          );



void FreeGensInComponent  ( Component *E
                          );



boolean FirstGen ( Combinaison *C
                 );



boolean NextGen ( Combinaison *C
                );



void DispCurrentComb ( Combinaison *C
                     );

 
#endif
 

