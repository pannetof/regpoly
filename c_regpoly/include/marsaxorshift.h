
#include "vectorsF2.h"
#include "combinaisons.h"

 
#ifndef MARSAXORSHIFT_H
#define MARSAXORSHIFT_H
 


#define TYPE1 1
#define TYPE21 21
#define TYPE22 22
#define TYPE23 23
#define TYPE24 24
#define TYPE25 25
#define TYPE3 3
#define TYPE4 4
#define TYPEGENERAL 5


typedef struct {
  int type;
  int w;
  int r;
  int nbxorshift;
  int nbmi;
  int *mi;
  int *mi_nb;
  int *ai;
  int m;
  int a,b,c;
  int *paramstype2p;
  int *paramstype2q;
  int paramstype3N;
  int *paramstype3No;
  int *paramstype3Shift;
} paramMarsaXorshift;



char* MarsaXorshiftName ( void
                        );



void InitParamMarsaXorshiftTypeI ( Generateur *Gen,
                                   int w,
                                   int a,
                                   int b,
                                   int c,
                                   int L
                                 );



void InitParamMarsaXorshiftTypeII ( Generateur *Gen,
                                    int type,
                                    int w,
                                    int r,
                                    int m,
                                    int *paramstype2p,
                                    int *paramstype2q,
                                    int L
                                  );



void InitParamMarsaXorshiftTypeIII ( Generateur *Gen,
                                     int w,
                                     int r,
                                     int paramstype3N,
                                     int *paramstype3No,
                                     int *paramstype3Shift,
                                     int L
                                   );



void InitParamMarsaXorshiftTypeIV ( Generateur *Gen,
                                    int a,
                                    int b,
                                    int c,
                                    int d,
                                    int w,
                                    int r,
                                    int m,
                                    int L
                                   );

 void InitParamMarsaXorshiftTypeGeneral ( Generateur *Gen,
                                               int w,
                                               int r,
                                               int nbxorshift,
                                               int nbmi,
                                               int *mi,
                                               int *mi_nb,
                                               int *ai,
                                               int L
                                             );



void InitMarsaXorshift ( Generateur *Gen,
                         BitVect *init,
                         BitVect *retour
                       );



void CopyMarsaXorshift ( Generateur *Gen1,
                         Generateur *Gen2
                       );



Generateur* AllocMarsaXorshift ( int k
                               );



void FreeMarsaXorshift ( Generateur *G
                       );



void MarsaXorshiftType1 ( Generateur *Gen,
                          BitVect *retour
                        );



void MarsaXorshiftType2 ( Generateur *Gen,
                          BitVect *retour
                        );



void MarsaXorshiftType3 ( Generateur *Gen,
                          BitVect *retour
                        );



void MarsaXorshiftType4 ( Generateur *Gen,
                          BitVect *retour
                        );



void MarsaXorshiftTypeGeneral ( Generateur *Gen,
                                BitVect *retour
                              );



void DispMarsaXorshift ( Generateur *Gen
                       );



void ReadDataMarsaXorshift ( Component *E,
                             char *nomfich,
                             boolean meme,
                             int L
                           );



void ProduceParamsMarsaXorshift ( char *filename,
                                  int type,
                                  int nb,
                                  int w,
                                  int r,
                                  int m,
                                  int nbxorshift
                                );

 
#endif
 

