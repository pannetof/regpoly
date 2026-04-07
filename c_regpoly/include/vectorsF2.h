
#include <stdlib.h>
#include <stdint.h>


#ifndef VECTORSF2_H
#define VECTORSF2_H
 



#define WL      32

 
#ifndef DEFTRUEFALSE
#define DEFTRUEFALSE 

#define TRUE     1
#define FALSE    0 
#endif 


 
#ifndef USE_BOOLEAN
#define USE_BOOLEAN 
typedef int    boolean; 
#endif 



typedef struct{
   int n;
   uint32_t *vect;
} BitVect;



void AllocBV ( BitVect *A,
               int l
             );



void FreeBV ( BitVect *A
            );



void BVCanonic ( BitVect *B,
                 int t
               );



void AllOnes ( BitVect *B
             );



void mask ( BitVect *B,
            int t
          );
void invmask ( BitVect *B,
               int t
             );



int ValBitBV ( BitVect *A,
               int noBit
             );



void PutBitBV ( BitVect *A,
                int noBit,
                int valBit
              );



void PutBVToZero ( BitVect *A
                 );



void CopyBV ( BitVect *A,
              BitVect *B
            );
  


void CopyBVPart ( BitVect *A,
                  BitVect *B,
                  int l
                );



void XORBV ( BitVect *A,
             BitVect *B,
             BitVect *C
           );



void XOR2BV ( BitVect *A,
              BitVect *B,
              BitVect *C,
              BitVect *D
            );



void ANDBV ( BitVect *A,
             BitVect *B,
             BitVect *C
           );



void ANDBVSelf ( BitVect *A,
                 BitVect *B
               );



void ANDBVMask ( BitVect *A,
                 BitVect *B,
                 int t
               );



void ANDBVInvMask ( BitVect *A,
                    BitVect *B,
                    int t
                  );



void XORBVSelf ( BitVect *A,
                 BitVect *B
               );



void BVLShift ( BitVect *R,
                BitVect *A,
                int n
              );



void BVRShift ( BitVect *R,
                BitVect *A,
                int n
              );



void BVLShiftSelf ( BitVect *R,
                    int n
                  );



void BVLS1Self ( BitVect *R
               );



void BVRShiftSelf ( BitVect *R,
                    int n
                  );


void BVLRotativeShift ( BitVect *R,
                        BitVect *A,
                        int n,
                        int w
                      );



void InverseBV ( BitVect *A
               );



void DispBitVect ( BitVect *A,
                   int l
                 );



void RandVect ( BitVect *v
              );



boolean VerifBitsCommuns ( BitVect *ds1,
                           BitVect *ds2
                         );



typedef struct{
  BitVect **lignes;
  int nblignes;
  int t;
  int l;
} Matrix;



void AllocMat ( Matrix* m,
                int nblines,
                int l,
                int t
              );



void FreeMat ( Matrix *m
             );



void CopyMat ( Matrix *m,
               Matrix *ms,
               int nl,
               int t
             );


int CompleteElimination ( Matrix *m,
                          int nblignes,
                          int l,
                          int t
                        );







void TransposeMatrices ( Matrix *T,
                         Matrix *M                       );



void ExchangeVect ( Matrix *m,
                    int i,
                    int j
                  );



void XorVect ( Matrix *m,
               int r,
               int s,
               int min,
               int max
             );



void DispMat ( Matrix *m,
               int t,
               int l,
               int nblines
             );



boolean InverseMatrix ( Matrix *MInv,
                        Matrix *M
                      );

 
#endif
 

