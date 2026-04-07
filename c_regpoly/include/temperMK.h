
#include "vectorsF2.h"
#include "combinaisons.h"
#include "tempering.h"
#include "mecf.h"

 
#ifndef TEMPERMK_H
#define TEMPERMK_H
 


#define TEMPMKID 1
#define TEMPNAME "Tempering MK"


typedef struct {
  int type;
  int u;
  int eta;
  int mu;
  int l;
  int random;
  BitVect b;
  BitVect b_original;
  BitVect c;
  BitVect c_original;
  BitVect BestB;
  BitVect BestC;
  boolean Optimize;
  boolean DispProgress;
  int limitv;
} paramTempMK;


int GiveTemperMKID (void);


void InitTemperMK (Transformation *T,
                   int w, int eta, int mu, int u, int l,
                   int random, BitVect *b, BitVect *c,
                   boolean Optimize, int type,
                   boolean DispProgress, int limitv);



void UpdateTemperMK ( Transformation *T
                    );



void TemperingMK ( BitVect *RES,
                   Transformation *T
                 );



void InverseTemperMK ( BitVect *RES,
                       Transformation *T
                     );



void CopyTemperMK ( Transformation *T1,
                    Transformation *T2
                  );



Transformation* AllocTemperMK ( void
                              );



void FreeTemperMK ( Transformation *T
                  );



void DispTemperMK ( Transformation *T
                  );



void TestMETemperMK ( Combinaison *C,
                      paramMecf *params
                    );

 
#endif
 

