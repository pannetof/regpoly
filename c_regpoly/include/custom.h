
#include "vectorsF2.h"
#include "combinaisons.h"

 
#ifndef CUST_H
#define CUST_H
 




typedef struct {
  char nomtrans[15];
  int row;
  int col;
  int portee;
  int nbbit;
  BitVect vecteur;
  void (*Trans) ( Generateur *Gen,
                  BitVect *etat,
                  int row,
                  int col,
                  int portee,
                  int nbbit,
                  BitVect *vecteur
                );
} custTrans;



typedef struct {
  int nbTrans;
  custTrans **specCustTrans;    BitVect Mask;   BitVect inv; BitVect temp1; BitVect temp2; 
} paramCust;



char* CustName ( void
               );



void InitParamCust ( Generateur *Gen,
                     int k,
                     paramCust *p,
                     int L
		   );



void InitCust ( Generateur *Gen,
                BitVect *init,
                BitVect *retour
              );



void CopyCust ( Generateur *G1,
                Generateur *G2
              );



Generateur* AllocCust ( int k
                      );



void FreeCust( Generateur *G
             );



void Cust ( Generateur *Gen,
            BitVect *g_etat
          );



void DispCust ( Generateur *Gen
              );



void ReadDataCust ( Component *E,
                    char *nomfich,
                    boolean meme,
                    int L
                  );



void LeftShift ( Generateur *Gen,
                 BitVect *etat,
                        int row,
                        int col,
                        int portee,
                        int nbbit,
                        BitVect *vecteur
                      );



void RightShift ( Generateur *Gen,
                  BitVect *etat,
                         int row,
                         int col,
                         int portee,
                         int nbbit,
                         BitVect *vecteur
                       );



void LeftShiftAND ( Generateur *Gen,
                    BitVect *etat,
                           int row,
                           int col,
                           int portee,
                           int nbbit,
                           BitVect *vecteur
                         );



void RightShiftAND ( Generateur *Gen,
                     BitVect *etat,
                            int row,
                            int col,
                            int portee,
                            int nbbit,
                            BitVect *vecteur
                          );



void XORIF ( Generateur *Gen,
             BitVect *etat,
                    int row,
                    int col,
                    int portee,
                    int nbbit,
                    BitVect *vecteur
                  );

 
#endif
 

