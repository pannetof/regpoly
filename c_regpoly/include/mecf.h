
#include "combinaisons.h"

 
#ifndef MECF_H
#define MECF_H
 


typedef boolean *Set;


typedef struct MECF {
  boolean meverif;
  boolean cfverif;
  boolean meverified;
  boolean cfverified;
  Set phi12;
  Set phi4;
  Set psi12;
  int L;
  void (*method)(Combinaison *C, struct MECF *params);
  int *ecart;
  int *delta;
  int *Lambda;
  int *ecartCF;
  int se;
  int secf;
  int mse;
  int msecf;
} paramMecf;


void AllocMecf (paramMecf *params, int Lmax);


void FreeMecf (paramMecf *params);


void InitParamMecf (paramMecf *params,
                    int *delta,
                    int mse, int msecf,
                    int meverif, int cfverif);


#define METHOD_MATRICIAL   0
#define METHOD_DUALLATTICE 1
#define METHOD_NOTHING     2

void SetMethod (paramMecf *params, int method);


#define TestEquid(C, params)   (params)->method(C, params) 


void TestME (Combinaison *C, paramMecf *params);


void TestNothing (Combinaison *C, paramMecf *params);


void TestCF (Combinaison *C, paramMecf *params);


boolean isPresqueME (paramMecf *params);


boolean isQuasiME (paramMecf *params);


boolean isME (paramMecf *params);


boolean isQuasiCF (paramMecf *params);


boolean isCF (paramMecf *params);


int DispTable (Combinaison *C, paramMecf *params, char type);


void ConvGaps (Combinaison *C, paramMecf *params);


void DispME (paramMecf *params);


void DispCF (paramMecf *params);


void PrepareMat (Combinaison *C, Matrix *Mat, int maxdim);


void PrepareMat2 (Combinaison *C, Matrix *Mat, int maxdim, int rotative_shift);


int DimensionEquid (Matrix *m, int kg, int l, int smax);


int ResolutionEquid (Matrix *m, int kg, int t, int *indices);


void ComputeSets (Combinaison *C, paramMecf *params, int indmax);


void SetPsi12 (Combinaison *C, paramMecf *params);


void TestMESpecial (Combinaison *C, paramMecf *params, int rotative_shift);


int GapSum (paramMecf *params);

 
#endif
 

