#include <stdint.h>

#ifndef REGPOLY_H
#define REGPOLY_H
 

 
#ifndef DEFTRUEFALSE
#define DEFTRUEFALSE 
#define FALSE    0
#define TRUE     1
 
#endif 


 
#ifndef USE_BOOLEAN
#define USE_BOOLEAN 
typedef int    boolean;
 
#endif 



int InverseModN ( int x,
                  int modulo
                );



uint32_t MRG32k3a ( void
               );



double doubleMRG32k3a ( void
               );



void SetMRG32k3a ( double S1,
                   double S2,
                   double S3,
                   double S4,
                   double S5,
                   double S6
                 );



float floatmin ( float a,
                 float b
               );



float floatmax ( float a,
                 float b
               );



int intmin ( int a,
             int b
           );



int intmax ( int A,
             int B
           );



int gcd ( int x,
          int y
        );



void ReadLn ( FILE *f
            );



boolean ValidMersennePrime ( int exp
                           );

 
#endif
 

