/*
 * AC1D.c — One-dimensional additive cellular automaton generator.
 *
 * Implements PRNG components based on 1D and 2D linear cellular automata
 * (rules 90/150 and generalizations). Supports construction from rule
 * vectors, transition-matrix diagonalization, and parameter generation.
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <time.h>
#include "AC1D.h"
#include "polynomials.h"
#include "regpoly.h"

#define AC1DNAME "AC-1D"

static void ConstructD(BitVect *MaskShift, BitVect *MatrixR, int Dimension) {
  
  int i,j;
  for (i=0;i<(2*Dimension)-1;i++)
    PutBVToZero(&MaskShift[i]);
  for (i=0;i<Dimension;i++)
    for (j=0;j<Dimension;j++)
      PutBitBV(&MaskShift[(j-i+2*Dimension-1)%(2*Dimension-1)],i,ValBitBV(&MatrixR[i],j));
}

static void ConstructMatrix(BitVect *MaskShift, BitVect *MatrixR, int Dimension) {
  int i,j;
  for (j=0;j<Dimension;j++)
    PutBVToZero(&MatrixR[j]);
  for (i=0;i<Dimension;i++)
    for (j=0;j<Dimension;j++)
      PutBitBV(&MatrixR[i],j,ValBitBV(&MaskShift[(j-i+2*Dimension-1)%(2*Dimension-1)],i));
}
/* Returns the name string for this generator type. */
char* AC1DName(void) {
  return AC1DNAME;
}

/* Displays the transition matrix of the AC1D generator in Mathematica list format. */
void DispAC1D(Generateur *Gen)
{
  int N,j;
  N=((paramAC1D *)(((Generateur*)Gen)->ParamGen))->N;
  for (j=0;j<N;j++) {
    DispBitVect(&(((paramAC1D *)Gen->ParamGen)->MatrixR[j]),N);
    printf(",\n");
  }
  /*
  printf("\n");
  for (j=0;j<2*N-1;j++) {
    DispBitVect(&(((paramAC1D *)Gen->ParamGen)->MaskShift[j]),N);
    printf("\n");
  }
  */
}
/* Initializes all parameters of an AC1D generator component. */
void InitParamAC1D(Generateur *Gen, int N, BitVect *MaskShift, BitVect *MatrixR, int taillevoisinage, int L) {
  int i;

  ((paramAC1D *)Gen->ParamGen)->N=N;
  for (i=0;i<2*N-1;i++)
    CopyBV(&(((paramAC1D *)Gen->ParamGen)->MaskShift[i]),&MaskShift[i]);
  for (i=0;i<N;i++)
    CopyBV(&(((paramAC1D *)Gen->ParamGen)->MatrixR[i]),  &MatrixR[i]);


  ((paramAC1D *)(Gen->ParamGen))->taillevoisinage=taillevoisinage;

  Gen->k=N;
  Gen->L=L;
  Gen->Step=1;    
  Gen->smax=INT_MAX;
  Gen->DispGen=DispAC1D;
  Gen->DispName=AC1DName; 
  Gen->InitGen=InitAC1D;
  Gen->Iteration=AC1D;
  Gen->CopyGen=CopyAC1D;
  Gen->AllocGen = AllocAC1D;
  Gen->PolyChar = polychar; //fonction generale
  PutBVToZero(&(Gen->GenState));
    
}
void CopyAC1D(Generateur *Gen1, Generateur *Gen2) {
  int i,N;

  N=((paramAC1D *)(Gen1->ParamGen))->N             =  ((paramAC1D *)(Gen2->ParamGen))->N;
  ((paramAC1D *)(Gen1->ParamGen))->taillevoisinage =  ((paramAC1D *)(Gen2->ParamGen))->taillevoisinage;
  for (i=0;i<2*N-1;i++)
    CopyBV(&(((paramAC1D *)(Gen1->ParamGen))->MaskShift[i]),&(((paramAC1D *)(Gen2->ParamGen))->MaskShift[i]));
  for (i=0;i<N;i++)
    CopyBV(&(((paramAC1D *)(Gen1->ParamGen))->MatrixR[i]),&(((paramAC1D *)(Gen2->ParamGen))->MatrixR[i]));

  Gen1->k             =Gen2->k;    
  Gen1->L             =Gen2->L;
  Gen1->smax          =Gen2->smax;
  Gen1->Step          =Gen2->Step;    
  Gen1->DispGen       =Gen2->DispGen;
  Gen1->DispName      =Gen2->DispName;
  Gen1->InitGen       =Gen2->InitGen;
  Gen1->Iteration     =Gen2->Iteration;
  Gen1->PolyChar      =Gen2->PolyChar;
  Gen1->CopyGen       =Gen2->CopyGen;
  Gen1->AllocGen      =Gen2->AllocGen;
  CopyBV(&(Gen1->GenState),&(Gen2->GenState));
}

Generateur* AllocAC1D(int k) {
  Generateur *G;
  int i;
  G = (Generateur *) malloc (sizeof(Generateur));
  AllocGen(G,k);
  G->ParamGen= (paramAC1D *) malloc (sizeof(paramAC1D));
  ((paramAC1D *)(G->ParamGen))->MaskShift       = (BitVect *) malloc(sizeof(BitVect)*(2*k-1));
  ((paramAC1D *)(G->ParamGen))->MatrixR         = (BitVect *) malloc(sizeof(BitVect)*k);

  for (i=0;i<2*k-1;i++)
    AllocBV(&(  ((paramAC1D *)G->ParamGen)->MaskShift[i]),k);
  for (i=0;i<k;i++)
    AllocBV(&(  ((paramAC1D *)G->ParamGen)->MatrixR[i]),k);
  return G;
}

void FreeAC1D(Generateur *G) {
  int i;
  for (i=0;i<2*G->k-1;i++)
    FreeBV(&(  ((paramAC1D *)G->ParamGen)->MaskShift[i]));
  free( ((paramAC1D *)G->ParamGen)->MaskShift);

  for (i=0;i<G->k;i++)
    FreeBV(&(  ((paramAC1D *)G->ParamGen)->MatrixR[i]));
  free( ((paramAC1D *)G->ParamGen)->MatrixR);
  
  free(G->ParamGen);
  FreeGen(G);
  free(G);
}

/* Initializes the generator state from init and copies the first L bits to retour. */


#define PPPPP 16

void InitAC1D(Generateur *Gen, BitVect *init, BitVect *retour)
{
  CopyBV (&(Gen->GenState), init);
  CopyBVPart (retour,init, Gen->L);
  if (0) {
    BitVect temp;
    AllocBV(&temp,Gen->L);
    BVRShift(&temp, retour, PPPPP);
    BVLShift(&temp,retour,Gen->L-PPPPP);
    XORBVSelf(retour,&temp);
    FreeBV(&temp);
  }
}

#define STATE Gen->GenState
#define DIM ((paramAC1D *)(Gen->ParamGen))->N
#define DD(q) ((paramAC1D *)(Gen->ParamGen))->MaskShift[q]
#define VOIS ((paramAC1D *)(Gen->ParamGen))->taillevoisinage
void  AC1D (Generateur *Gen, BitVect *retour)
{
  int i;
  BitVect Temp;
  BitVect NewState;
  AllocBV(&Temp,Gen->k);
  AllocBV(&NewState,Gen->k);

  PutBVToZero(&NewState);
  for (i=0;i<=VOIS;i++) {
    BVLShift(&Temp,&(STATE),i);
    ANDBV(&Temp,&Temp,&(DD(i)));
    XORBV(&NewState,&NewState,&Temp);
  }
  for (i=2*DIM-2;i>=2*DIM-1-VOIS;i--) {
    BVRShift(&Temp,&(STATE),(2*DIM-1)-i);
    ANDBV(&Temp,&Temp,&(DD(i)));
    XORBV(&NewState,&NewState,&Temp);
  }

  CopyBV(&(STATE),&NewState);
  CopyBVPart(retour,&(STATE),Gen->L);
  if (0) {
    BitVect temp;
    AllocBV(&temp,Gen->L);
    BVRShift(&temp, retour, PPPPP);
    BVLShift(&temp,retour,Gen->L-PPPPP);
    XORBVSelf(retour,&temp);
    FreeBV(&temp);
  }
  FreeBV(&Temp);
  FreeBV(&NewState);
  //  printf("------------->");DispBitVect(&(STATE),20,0);printf("\n");
}
/* Reads the generator data file *f containing the search criteria.
   The file format must be strictly followed. */
#define SELF 16UL
#define TOP 8UL
#define LEFT 4UL
#define BOT 2UL
#define RIGHT 1UL

void ReadDataAC1D(Component *E, char *filename, boolean same, int L)
{
  int M, taillevoisinage,i,j,rule,n,m,S,k;
  BitVect *Temp, *MatrixR;
  char *line, *word;
  int R;
  char F;
  FILE *f;
  Generateur *p;

  f=fopen(filename,"r");
  if (f==NULL) {
    printf("File not found: %s\n",filename);
    exit(1);
  }

  fscanf(f, "%d", &M);          /* read M */
  ReadLn(f) ;

  F = fgetc(f) ;

  AllocGensInComponent(E, M, same);
  
  //  printf("Ce module est bugg�.  Il faut le corriger.  La matrice de transition TransitionMatrix()\n");
  //printf("et celle de MatrixR ne correspondent pas\n");
  
  // exit(1);

  if (F == 'W') {
    fscanf(f,"%d",&n);
    ReadLn(f);



    taillevoisinage=1;
    for (i=0;i<M;i++) {
      if ( (p = AllocAC1D(n))==NULL) {
  printf("Error in ReadDataAC1D()\n");
  exit(1);
      }
      Temp = (BitVect *) malloc(sizeof(BitVect)*(2*n+1));
      MatrixR = (BitVect *) malloc(sizeof(BitVect)*n);

      for (j=0;j<2*n-1;j++)
  AllocBV(&Temp[j],n);
      for (j=0;j<n;j++)
  AllocBV(&MatrixR[j],n);
      
      for (j=0;j<n;j++)
  PutBVToZero(&MatrixR[j]);
      for (j=0;j<n;j++) {
  fscanf(f,"%d",&rule);
  if (rule==90 || rule==150 || rule==60 || rule==240)
    if (j!=0) 
      PutBitBV(&MatrixR[j],j-1,1);
  if (rule==90 || rule==150 || rule==102 || rule==170)
    if (j!=n-1) 
      PutBitBV(&MatrixR[j],j+1,1);
  if (rule==204 || rule==150 || rule==60 || rule==102)
    PutBitBV(&MatrixR[j],j,1);
  if (rule!=90 && rule!=150 && rule!=204 && rule!=60 &&rule!=102 &&rule!=170 &&rule!=240 &&rule!=0) {
    printf("Rule %d is non-linear\n",rule);
    exit(1);
  }  
      }
      ReadLn(f);
      ConstructD(Temp,MatrixR,n);
      InitParamAC1D(p,n,Temp,MatrixR,taillevoisinage,L);
      AddGenInComponent(E,p);
      FreeAC1D(p);
      
      for (j=0;j<2*n-1;j++)
  FreeBV(&Temp[j]);
      for (j=0;j<n;j++)
  FreeBV(&MatrixR[j]);
      free(Temp);
      free(MatrixR);
    }
    
  } 
  else if (F == 'C') {
    fscanf(f,"%d",&n);
    ReadLn(f);


    line =(char*) malloc (sizeof(char)*n*12);
    for (i=0;i<M;i++) {
      if ( (p = AllocAC1D(n))==NULL) {
  printf("Error in ReadDataAC1D()\n");
  exit(1);
      }
      Temp = (BitVect *) malloc(sizeof(BitVect)*(2*n+1));
      MatrixR = (BitVect *) malloc(sizeof(BitVect)*n);

      for (j=0;j<2*n-1;j++)
  AllocBV(&Temp[j],n);
      for (j=0;j<n;j++)
  AllocBV(&MatrixR[j],n);
      
      for (j=0;j<n;j++)
  PutBVToZero(&MatrixR[j]);
      taillevoisinage = 0;
      for (j=0;j<n;j++) {
  fgets(line,n*12,f);
  
  //  printf("%s",line);fflush(stdout);

  word = strtok(line," ");
  while (word!=NULL) {
    sscanf(word,"%d",&R);
    if ((j+R>= 0) && j+R<n) {
      PutBitBV(&MatrixR[j],j+R,1);
      if (R<0)
        R=-R;
      if (R> taillevoisinage)
        taillevoisinage =R;
    }
    word = strtok(NULL," ");
  }

      }
      
      ConstructD(Temp,MatrixR,n);
      InitParamAC1D(p,n,Temp,MatrixR,taillevoisinage,L);
      AddGenInComponent(E,p);
      FreeAC1D(p);

      for (j=0;j<2*n-1;j++)
  FreeBV(&Temp[j]);
      for (j=0;j<n;j++)
  FreeBV(&MatrixR[j]);
      free(Temp);
      free(MatrixR);
    }
  }
  else if (F == 'A') {
    fscanf(f,"%d",&n);
    ReadLn(f);
    fscanf(f, "%d", &taillevoisinage);
    for (i=0;i<M;i++) {
      if ( (p = AllocAC1D(n))==NULL) {
  printf("Error in ReadDataAC1D()\n");
  exit(1);
      }

      Temp = (BitVect *) malloc(sizeof(BitVect)*(2*n+1));
      MatrixR = (BitVect *) malloc(sizeof(BitVect)*n);

      for (j=0;j<2*n-1;j++)
  AllocBV(&Temp[j],n);
      for (j=0;j<n;j++)
  AllocBV(&MatrixR[j],n);
      

      for (j=0;j<n;j++) {
  while (BVisZero(&MatrixR[j])) {
    RandVect(&MatrixR[j]);
    ANDBVMask(&MatrixR[j],&MatrixR[j],taillevoisinage*2+1);
    if (j<taillevoisinage)
      BVLShiftSelf(&MatrixR[j],taillevoisinage-j);
    else
      BVRShiftSelf(&MatrixR[j],j-taillevoisinage);
    ANDBVMask(&MatrixR[j],&MatrixR[j],n);
  }
      } 
       /*
       for (j=0;j<2*n-1;j++)
   PutBVToZero(&Temp[j]);
       for (j=(2*n-1-taillevoisinage)%(2*n-1) ; j!= taillevoisinage ; j=(j+1)%(2*n-1)) {
    RandVect(&Temp[j]);
    if (j>n)
      ANDBVInvMask(&Temp[j],&Temp[j],(2*n-1)-j);
    else
      ANDBVMask(&Temp[j],&Temp[j],n-j);
      }
       RandVect(&Temp[taillevoisinage]);
       ANDBVMask(&Temp[taillevoisinage],&Temp[taillevoisinage],n-taillevoisinage);
      ConstructMatrix(Temp,MatrixR,n);
      */
      ConstructD(Temp,MatrixR,n);
      InitParamAC1D(p,n,Temp,MatrixR,taillevoisinage,L);
      AddGenInComponent(E,p);
      FreeAC1D(p);

    }
    for (j=0;j<2*n-1;j++)
      FreeBV(&Temp[j]);
    for (j=0;j<n;j++)
      FreeBV(&MatrixR[j]);
    free(Temp);
    free(MatrixR);
  }
  else if (F == '2' || F =='R') {
    fscanf(f,"%d %d",&n, &m);
    S=n*m;  /* matrix R will be S x S */
    taillevoisinage = m;

    Temp = (BitVect *) malloc(sizeof(BitVect)*(2*S+1));
    MatrixR = (BitVect *) malloc(sizeof(BitVect)*S);
    for (j=0;j<2*S-1;j++)
      AllocBV(&Temp[j],S);
    for (j=0;j<S;j++)
      AllocBV(&MatrixR[j],S);

    for (k=0;k<M;k++) { 

      if ( (p = AllocAC1D(S))==NULL) {
  printf("Error in ReadDataAC1D()\n");
  exit(1);
      }     
      for (j=0;j<S;j++)
  PutBVToZero(&MatrixR[j]);
      for (i=0;i<n;i++)
  for (j=0;j<m;j++) {
    if (F=='2')
      fscanf(f,"%d",&rule);
    else
      rule = MRG32k3a() & 0x1fUL;
    if (rule & SELF) {
      PutBitBV(&MatrixR[i*m+j],i*m+j,1);
    }
    if (rule & LEFT) {
      if (j!=0)
        PutBitBV(&MatrixR[i*m+j],i*m+j-1,1);
    }
    if (rule & RIGHT) {
      if (j!=m-1)
        PutBitBV(&MatrixR[i*m+j],i*m+j+1,1);
    }
    if (rule & TOP) {
      if (i!=0)
        PutBitBV(&MatrixR[i*m+j],(i-1)*m+j,1);
    }
    if (rule & BOT) {
      if (i!=n-1)
        PutBitBV(&MatrixR[i*m+j],(i+1)*m+j,1);
    }
    //printf("\n"); 
    
  }
      ReadLn(f); 
      ConstructD(Temp,MatrixR,S); 
      InitParamAC1D(p,S,Temp,MatrixR,taillevoisinage,L); 
      AddGenInComponent(E,p);       
      FreeAC1D(p); 
    }
    for (j=0;j<2*S-1;j++) 
      FreeBV(&Temp[j]); 
    for (j=0;j<S;j++) 
      FreeBV(&MatrixR[j]);
    free(Temp); 
    free(MatrixR);
  }
  else printf("Unknown character '%c' in file %s\n", F, filename);
  fclose(f);


}

void ProduceParamsAC1D ( char *filename, int nb, char F, int n, int m, int neighborhood ) {

  FILE *f;
  int R,*Poly,i,j,count,*rule2D,*rule,rules[2] = {90,150}; /* see Palash Sarkar */
  BitVect *Temp, *MatrixR, BVPoly, LITest;
  Generateur *p;
  boolean ok;


  f = fopen (filename, "w");
  count=0;
  fprintf(f,"%d\n",nb);
  SetMRG32k3a((double)time(NULL),(double)time(NULL),(double)time(NULL),
        (double)time(NULL),(double)time(NULL),(double)time(NULL));

  if (F!='2')
    m=1;

  Poly = (int*) malloc((n*m+1)*sizeof(int));
  AllocBV(&BVPoly,n*m+1);
  Temp = (BitVect *) malloc(sizeof(BitVect)*(2*n*m+1));
  MatrixR = (BitVect *) malloc(sizeof(BitVect)*n*m);
  for (j=0;j<2*n*m-1;j++)
    AllocBV(&Temp[j],n*m);
  for (j=0;j<n*m;j++)
    AllocBV(&MatrixR[j],n*m);
  AllocBV(&LITest,n*m);

  if (F == 'C') {
    fprintf(f,"C %d\n",n);

    while (count<nb) {
      if ( (p = AllocAC1D(n))==NULL) {
  printf("Error in ProduceParamsAC1D()\n");
  exit(1);
      }
      for (j=0;j<n;j++)
  PutBVToZero(&MatrixR[j]);
      for (j=0;j<2*n-1;j++)
  PutBVToZero(&Temp[j]);
      
      for (j=0;j<n;j++) {
  ok=FALSE;
  while (!ok) {
    for (R = intmax(0,j-neighborhood); R< intmin(n,j+neighborhood+1);R++) {
      if (MRG32k3a()%2) {
        PutBitBV(&MatrixR[j],R,1);
        ok = TRUE;
      }
    }
    if (j>0) {
      XORBV(&LITest,&MatrixR[j],&MatrixR[j-1]);
      if (BVisZero(&LITest)) {
        PutBVToZero(&MatrixR[j]);
        ok = FALSE;
      }
    }
  }
      }
      ConstructD(Temp,MatrixR,n);
      InitParamAC1D(p,n,Temp,MatrixR,neighborhood,32);
      polychar ( p, Poly, &BVPoly);
      FreeAC1D(p);            
      if (PrimitifPolynomial(Poly,n)) {
  for (j=0;j<n;j++) {
    for (R = intmax(0,j-neighborhood); R< intmin(n,j+neighborhood+1);R++)
      if ( ValBitBV(&MatrixR[j],R) )
        fprintf(f,"%d ",R-j);
    fprintf(f,"\n");

  }
  count++;
      }
    }
  } else if (F=='W') {
    fprintf(f,"W %d\n",n);
    rule = (int*) malloc(n * sizeof(int));
    while (count<nb) {
      if ( (p = AllocAC1D(n))==NULL) {
  printf("Error in ProduceParamsAC1D()\n");
  exit(1);
      }
      for (j=0;j<n;j++)
  PutBVToZero(&MatrixR[j]);
      for (j=0;j<2*n-1;j++)
  PutBVToZero(&Temp[j]);

      for (j=0;j<n;j++) {
  do {
    ok = TRUE;
    rule[j] = rules[MRG32k3a() % 2];
    if (rule[j]==90 || rule[j]==150)// || rule[j]==60 || rule[j]==240)
      if (j!=0) 
        PutBitBV(&MatrixR[j],j-1,1);
    if (rule[j]==90 || rule[j]==150)// || rule[j]==102 || rule[j]==170)
      if (j!=n-1) 
        PutBitBV(&MatrixR[j],j+1,1);
    if (rule[j]==150)// || rule[j]==204 || rule[j]==60 || rule[j]==102)
      PutBitBV(&MatrixR[j],j,1);
    if (j==0 || j==(n-1))
      if (BVisZero(&MatrixR[j]))
        ok = FALSE;
      
    if (j>0) {
      XORBV(&LITest,&MatrixR[j],&MatrixR[j-1]);
      if (BVisZero(&LITest)) {
        ok = FALSE;
        PutBVToZero(&MatrixR[j]);
      }
    }
  }while (!ok);
      }
      ConstructD(Temp,MatrixR,n);
      InitParamAC1D(p,n,Temp,MatrixR,1,32);
      //      DispAC1D(p);printf("\n");
      polychar ( p, Poly, &BVPoly);

      FreeAC1D(p);            
      if (PrimitifPolynomial(Poly,n)) {
  for (j=0;j<n;j++)
    fprintf(f,"%d ",rule[j]);
  fprintf(f,"\n");
  count++;
      }
    }
    free(rule);


  } else if (F=='2') {
    fprintf(f,"2 %d %d\n",n,m);
    rule2D = (int*) malloc(sizeof(int)*n*m);
    while (count<nb) {
      if ( (p = AllocAC1D(n*m))==NULL) {
  printf("Error in ProduceParamsAC1D()\n");
  exit(1);
      }
      for (j=0;j<n*m;j++)
  PutBVToZero(&MatrixR[j]);
      for (j=0;j<2*n*m-1;j++)
  PutBVToZero(&Temp[j]);

      for (i=0;i<n;i++) {
  for (j=0;j<m;j++) {
    while ( (*(rule2D+i*m+j) = MRG32k3a() & 0x1fUL) == SELF);
    *(rule2D+i*m+j) |= (LEFT ^ RIGHT); /* LEFT and RIGHT are mandatory */
    *(rule2D+i*m+j) |= (TOP ^BOT);

    //printf("%d\n",*(rule2D+i*m+j));
    if (*(rule2D+i*m+j) & SELF) {
      PutBitBV(&MatrixR[i*m+j],i*m+j,1);  
    }
    if (*(rule2D+i*m+j) & LEFT) {
      if (j!=0)
        PutBitBV(&MatrixR[i*m+j],i*m+j-1,1);
    }
    if (*(rule2D+i*m+j) & RIGHT) {
      if (j!=m-1)
        PutBitBV(&MatrixR[i*m+j],i*m+j+1,1);
    }
    if (*(rule2D+i*m+j) & TOP) {
      if (i!=0)
        PutBitBV(&MatrixR[i*m+j],(i-1)*m+j,1);
    }
    if (*(rule2D+i*m+j) & BOT) {
      if (i!=n-1)
        PutBitBV(&MatrixR[i*m+j],(i+1)*m+j,1);
    }
    if (BVisZero(&MatrixR[i*m+j]))
      j--;
  }
      }
      /*      PutBitBV(&MatrixR[0],1,1);
      PutBitBV(&MatrixR[1],0,1);
      *(rule2D) |= RIGHT;
      *(rule2D+1) |= LEFT;
      */
      ConstructD(Temp,MatrixR,n*m);
      InitParamAC1D(p,n*m,Temp,MatrixR,m,32);

       DispAC1D(p);printf("\n");
       //TransitionMatrix(p);
      polychar ( p, Poly, &BVPoly);

      FreeAC1D(p);            
      if (PrimitifPolynomial(Poly,n*m)) {
  for (i=0;i<n;i++)
    for (j=0;j<m;j++)
      fprintf(f,"%d ",*(rule2D+i*m+j));
  fprintf(f,"\n");
  count++;
      }
    }
    free(rule2D);





  }
  for (j=0;j<2*n-1;j++)
    FreeBV(&Temp[j]);
  for (j=0;j<n;j++)
    FreeBV(&MatrixR[j]);
  FreeBV(&LITest);
  free(Temp);
  free(MatrixR);

  free(Poly);
  FreeBV(&BVPoly);
  fclose(f);
}  

