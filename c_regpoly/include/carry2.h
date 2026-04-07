
#include "vectorsF2.h"
#include "combinaisons.h"

 
#ifndef CARRY_H
#define CARRY_H
#define NBMAT 8
#define NBTYPE 8
#define NBPARAMSINT 3
#define NBPARAMSULONG 3
 



typedef struct {
  int w;
  int r;
  int p;
  int m1;
  int m2;
  int m3;
  uint32_t paramsulong[NBMAT][NBPARAMSULONG];
  int   paramsint[NBMAT][NBPARAMSINT];
  int   type[NBMAT];
  uint32_t (*Matrices[NBTYPE])(uint32_t, int paramsint[NBPARAMSINT],
                            uint32_t paramsulong[NBPARAMSULONG]);
  void (*DisplayParams[NBTYPE])(int paramsint[NBPARAMSINT],
                                uint32_t paramsulong[NBPARAMSULONG]);
  int (*Cost[NBTYPE])(void);
  int wordno;
  uint32_t maskp;
  int i;
} paramCARRY;



char* CarryName ( void
                );



void InitParamCarry ( Generateur *Gen,
                      int w,
                      int r,
                      int p,
                      int m1,
                      int m2,
                      int m3,
                      int paramsint[NBMAT][NBPARAMSINT],
                      uint32_t paramsulong[NBMAT][NBPARAMSULONG],
                      int *type,
                      int wordno,
                      int L
                   );



void InitCarry ( Generateur *Gen,
                 BitVect *init,
                 BitVect *retour
               );



void CopyCarry ( Generateur *Gen1,
                 Generateur *Gen2
               );



Generateur* AllocCarry ( int k
                       );



void FreeCarry ( Generateur *G
               );



void Carry ( Generateur *Gen,
             BitVect *retour
           );



void DispCarry ( Generateur *Gen
               );



void ReadDataCarry ( Component *E,
                     char *nomfich,
                     boolean meme,
                     int L
                   );



void ProduceParamsCarry ( char *filename,
                          int nb,
                          int w,
                          int r,
                          int p,
                          int maxcost,
                          double seedmodifier
                       );



void CharCarry ( Generateur *Gen,
                 int coeff[],
                 BitVect *BVPoly
               );

 
#endif
 

