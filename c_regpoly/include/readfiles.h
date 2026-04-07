

#include "combinaisons.h"
#include "mecf.h"
#include "tuplets.h"



 
#ifndef UTIL_H
#define UTIL_H
 



void ReadSearch ( Combinaison *C,
                  int NbComp,
                  char *testfile,
                  paramMecf *pmecf,
                  paramTuplets *ptuplets,
                  int *nbessais,
                  boolean *MKopt,
                  int *Lmax
                );



void ReadGenDataFiles( Combinaison *C,
                       int J,
                       char **gendatafile,
                       int Lmax
                     );



void ReadTempering ( Component *E,
                     char *filename,
                     boolean *MKopt
                   );

 
#endif
 

