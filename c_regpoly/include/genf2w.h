
#include "vectorsF2.h"
#include "combinaisons.h"
#include <stdint.h>


#ifndef GENF2W_H
#define GENF2W_H
 


#define GENF2WLFSR     1
#define GENF2WPOLYLCG  0

typedef struct {
  int type;
  int r;
  int w;
  int nbcoeff;
  uint32_t *coeff;
  int *nocoeff;
  boolean normalbasis;
  int step;
  uint32_t modM;   int p; uint32_t masque;  uint32_t *table;  uint32_t (*multiply)(uint32_t a, uint32_t b, int w, uint32_t modM, uint32_t *table);  
} paramGenF2w;



char* GenF2wName ( void
                 );



Generateur* AllocGenF2w ( int k
                        );



void FreeGenF2w ( Generateur *G
                );



void InitParamGenF2w ( Generateur *Gen,
                       int w,
                       int r,
                       int nbcoeff,
                       int *nocoeff,
                       uint32_t *coeff,
                       uint32_t modM,
                       boolean normalbasis,
                       int type,
                       int step,
                       int L
                     );



void InitGenF2w ( Generateur *Gen,
                  BitVect *init,
                  BitVect *retour
                );



void CopyGenF2w ( Generateur *Gen1,
                  Generateur *Gen2
                );



void GenF2wLFSR32bit ( Generateur *Gen,
                       BitVect *retour
                     );
void GenF2wLFSR ( Generateur *Gen,
                  BitVect *retour
                );



void GenF2wPolyLCG ( Generateur *Gen,
                     BitVect *retour
                   );



void DispGenF2w ( Generateur *Gen
                );



void ReadDataGenF2w ( Component *E,
                      char *nomfich,
                      boolean same,
                      int L
                    );



void ProduceParamsGenF2w ( char *filename,
                           int nb,
                           int w,
                           int r,
                           int nbcoeff,
                           boolean normalbasis,
                           int sizetable
                         );



void ProduceOneGenF2w ( Combinaison *C,
                        int type,
                        int w,
                        int r,
                        int nbcoeff,
                        boolean normalbasis,
                        int sizetable,
                        int stepmax,
                        boolean primitive,
                        int Lmax
                      );



boolean UpdateParamsGenF2w( Generateur *Gen, int normalbasis, int sizetable, int stepmax, boolean primitive, int modM, int exhaustive);



void PrintUsableFormatREGPOLY ( Generateur *Gen, char *filename );
void PrintUsableFormatSSJ ( Generateur *Gen, char *filename  );

 
#endif
 

