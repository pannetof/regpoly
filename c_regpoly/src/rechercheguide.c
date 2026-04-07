#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vectorsF2.h"
#include "combinaisons.h"
#include "tempering.h"
#include "timer.h"
#include "mecf.h"
#include "tuplets.h"
#include "readfiles.h"
#include "temperMK.h"
#define SEPARATION printf("\n\n++++++++++++++++++++++++++++++++\
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")

void Seek(int NbComp,char *testfile, char **gendatafile);

int main (int argc, char *argv[] ){
  FILE *input;                       
  char **gendatafile,      // Array of names of data files for each component.        
       testfile[200];      // Name of the test file.
  int NbComp,              // Number of components for the combined generators to
                           // be looked at.
      j;
  if( argc < 4 ) {
    printf("\nUSAGE:  POL NbComp <test_file> <gen_data_file1> ... <gen_data_file(NbComp)> \n\n") ;
    exit(1) ;
  }
  NbComp = atoi(argv[1]);
  if( !(input = fopen(argv[2],"r")) ) {
    printf("\n Test file \" %s \" not found.\n\n", argv[2]) ;
    exit(1) ;
  }
  fclose(input);
  gendatafile = (char **) malloc(sizeof(char *)*NbComp);
  
  for(j=0;j<NbComp;j++){ 
    gendatafile[j] = (char *) malloc(sizeof(char)*200);
    strcpy(gendatafile[j],argv[j+3]);
  }
  strcpy(testfile,argv[2]);
  Seek (NbComp, testfile, gendatafile); // Launch of the search.
  return 0;
}

void Seek (int NbComp, char *testfile, char **gendatafile)
{
  Combinaison Combin;             // Contains all the data on the generators to 
                                  // look at.
  int nbgen,                      // Number of generators verified in total
    nbsel,                        // Number of generators that meets the criteria
    nbME,                         // Number of ME generators
    nbCF;                         // Number of CF generators
  boolean Continue;               // This variable is true until all generators have 
                                  // have been verified.
  timer_Chrono timer;             // Time taken to do the search
  int no_try;                     // Number of tempering already tried so far
  int nbtries;                    // Number of tempering to be tried per combined generator
  paramTuplets ptuplets;          // Structure that stores the criterion on non-successive
                                  // output values
  paramMecf pMecf;                // Structure that stores the criterion on successive output
                                  // values
  boolean MKopt;                  // Boolean that tells if optimized MK tempering is to be done
  int Lmax;                       // Maximum output size (bits)

  nbgen = nbsel = nbME = nbCF = 0;
  /* Read data in the test file (see module readfiles) */
  ReadSearch(&Combin, NbComp, testfile, &pMecf, &ptuplets, &nbtries, &MKopt, &Lmax);
  ReadGenDataFiles(&Combin, NbComp, gendatafile, Lmax);
  /* Combin is set to the first combined generator */
  if(!FirstGen(&Combin)){
    printf("First combined generator not found or invalid Combinaison\n");
    exit(1);
  }
  no_try=1;
  timer_Init(&timer);             // Initialization of the timer

  do {
    UpdateAllTrans(&Combin);      // Update all tempering transformation on all
                                  // components of Combin
    nbgen++;     
    if(MKopt)                     // If MKopt==TRUE, tries to optimize the MK tempering.
      TestMETemperMK(&Combin, &pMecf);      
    else                          // else computes the equidistribution of the combined gen.
      TestEquid(&Combin, &pMecf);

    if(isPresqueME(&pMecf)){     // If it meets the criterion defined by pMecf
      TestTuplets(&Combin,&ptuplets); // Test for the equidistribution of non-successive
                                  // output values as defined by ptuplets.
      if(isTuplets(&ptuplets)){   // If the current combined generator meets the
                                  // criterion on non-successive output values
                                  // or it has not been tested.
	DispCurrentComb(&Combin); // Display the current combined generator
	if(!isME(&pMecf)){       // If the generator is not ME
                                  // then the equidistribution is displayed on the screen
	  printf("\n  Dimension gaps for every resolution");
	  DispTable(&Combin,&pMecf,'l');// a mettre a 'l' 
	} 
	else{                     // If the generator is ME, prints a message saying so.
	  DispME(&pMecf); 
	  nbME++;
	} 
	DispTuplets(&ptuplets);  // Display the results for non-successive output values
                                 // if tested
	nbsel++; 
	SEPARATION; 
      }   
    }  
    fflush(stdout);              // To make it print the results right away
                                 // useful for rare good generators.
    if(nbtries==no_try){      // If nbtries tempering transformation have
      Continue=NextGen(&Combin); // been tried on one combined generator,
      no_try=1 ;               // we go to the next combined generator in Combin
    }else{
      Continue=TRUE;
      no_try++;
    }    
  } while (Continue);            // We keep looking until all combined generators
                                 // have been tested

  /* Displays a summary of the search */
  printf("\n===========================\n") ;
  printf("   Total   =  %10u  \n", nbgen) ;
  printf("\n") ;
  printf("     ME    =  %10u  \n", nbME) ;
  printf("   CF-ME   =  %10u  \n", nbCF) ;
  printf("  retained =  %10u  \n", nbsel) ;
  printf("---------------------------\n") ;
  printf(" CPU (sec) =   %5.2f       \n",timer_Val(timer,timer_sec) ) ;
  printf("===========================\n") ;
}      










