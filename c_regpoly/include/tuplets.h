

#include "vectorsF2.h"
#include "combinaisons.h"

 
#ifndef TUPLETS_H
#define TUPLETS_H
 


#define EQUIDISTRIBUTION 0
#define TVALUE           1
#define NEIGHBOR         2
#define MINDIST          3

#define SUM              0
#define MAX              1



typedef struct{
  int test;
  int testtype;
  boolean tupletsverif;
  boolean tupletsverified;
  int tupd;
  int *tuph;
  float *DELTA;
  float *gap;
  float *pourcentage;
  float treshold;
  int indice_max;
  float firstpart_sum;
  float secondpart_sum;
  float firstpart_max;
  float secondpart_max;
} paramTuplets;



void AllocTuplets ( paramTuplets *param,
                    int d
                  );



void FreeTuplets ( paramTuplets *param
                 );



void InitParamTuplets ( paramTuplets *params,
                        int d,
                        int *t,
                        float treshold,
                        int tupletsverif,
                        int test,
                        int testtype
                      );



void ChangeTest( paramTuplets *params,
                 int test,
                 int testtype
               );



void ChangeTreshold ( paramTuplets *params,
                      float treshold
                    );



void TestTuplets ( Combinaison *C,
                   paramTuplets *params
                 );



boolean isTuplets ( paramTuplets *params
                  );



void DispTuplets ( paramTuplets *params
                 );



void DispTupletsTest(paramTuplets *params);

 
#endif
 

