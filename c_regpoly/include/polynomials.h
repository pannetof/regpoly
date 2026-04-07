

#include <zen.h>
#include "vectorsF2.h"
#include "combinaisons.h"
#include "regpoly.h"

 
#ifndef POLYCHAR_H_
#define POLYCHAR_H_
 

void polychar ( Generateur *Gen,
                int *coeff,
                BitVect *BVPoly
              );



void polycharComb ( Combinaison *C,
                    int *coeff,
                    BitVect *BVPoly
                  );



void DispPolynomial ( int *coeff,
                      int degree
                    );



void TransitionMatrix ( Generateur *G
                      );



boolean IrreduciblePolynomial (  int *coeff,
                                 int degree
                              );



int PrimitifPolynomial ( int *coeff,
                         int degree
                       );



int PrimitifPolynomialSieving ( int *coeff,
                                int degree
                              );



void SequenceMinimalPolynomial ( int K,
                                 char * Sequence,
                                 int *coeff,
                                 BitVect *BVPoly
                               );



boolean ContainsMersenneExponent(int degre, int *poly, int Mexp);

  void PowPoly(ZENPoly *a, ZENPoly *b, ZENRing *F2, int pow);
void PowPuis2Poly(ZENPoly *a, ZENPoly *b, ZENRing *F2, int i);
 
 
#endif
 
