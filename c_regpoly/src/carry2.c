/*
 * carry2.c — Carry generator (add-with-carry / multiply-with-carry type).
 *
 * Implements a combined generator whose components apply linear matrix
 * transformations of configurable types (shift-based, companion matrix, etc.)
 * to 32-bit words in GF(2). Supports parameter search and file-based I/O.
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include "carry2.h"
#include "polynomials.h"
#include "regpoly.h"

#define CARRYNAME "Carry Generator"
void InitCarry32(Generateur *Gen, BitVect *, BitVect *retour);
void Carry32(Generateur *Gen, BitVect *retour);
static void TransitionMatrixCarry(Generateur *Gen);
static void DispMatSpecial (Matrix *m , int l,   int kg);
static void PrintUsableParameters(Generateur *Gen);
static uint32_t type0(uint32_t, int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]);
static uint32_t type1(uint32_t, int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]);
static uint32_t type2(uint32_t, int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]);
static uint32_t type3(uint32_t, int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]);
static uint32_t type4(uint32_t, int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]);
static uint32_t type5(uint32_t, int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]);
static uint32_t type6(uint32_t, int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]);
static uint32_t type7(uint32_t, int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]);


static void selecttype0(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]);
static void selecttype1(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]);
static void selecttype2(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]);
static void selecttype3(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]);
static void selecttype4(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]);
static void selecttype5(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]);
static void selecttype6(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]);
static void selecttype7(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]);


static void displaytype0(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]);
static void displaytype1(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]);
static void displaytype2(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]);
static void displaytype3(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]);
static void displaytype4(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]);
static void displaytype5(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]);
static void displaytype6(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]);
static void displaytype7(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]);

static int cost0(void);
static int cost1(void);
static int cost2(void);
static int cost3(void);
static int cost4(void);
static int cost5(void);
static int cost6(void);
static int cost7(void);


/* Returns the name string for this generator type. */
char* CarryName(void) {
  return CARRYNAME;
}
#define TM1 ((paramCARRY *)(Gen->ParamGen))->m1
#define TM2 ((paramCARRY *)(Gen->ParamGen))->m2
#define TM3 ((paramCARRY *)(Gen->ParamGen))->m3
#define TTYPE(j) ((paramCARRY *)(Gen->ParamGen))->type[j]
#define TPARINT(j) ((paramCARRY *)(Gen->ParamGen))->paramsint[j]
#define TPARULONG(j) ((paramCARRY *)(Gen->ParamGen))->paramsulong[j]
#define PARINT(i,j) ((paramCARRY *)(Gen->ParamGen))->paramsint[j][i]
#define PARULONG(i,j) ((paramCARRY *)(Gen->ParamGen))->paramsulong[j][i]



void DispCarry(Generateur *Gen)
{
  int j,type,cost;
  int *Poly;
  int count = 0;
  BitVect Dummy;
  //  Gen->k = 96;
  Poly = (int*) malloc((DEGGEN(Gen)+1)*sizeof(int));
  AllocBV(&Dummy, DEGGEN(Gen)+1);
  printf("%d\n",DEGGEN(Gen));fflush(stdout);

  polychar(Gen, Poly, &Dummy);

  for (j=0;j<=DEGGEN(Gen);j++) {
  
    if (Poly[j]==1)
      count++;
  }
  //  DispPolynomial(Poly,DEGGEN(Gen));
  /*
  if (PrimitifPolynomial(Poly,DEGGEN(Gen)))
    printf("Full Period.\n");
  else{
    printf("Not Full Period\n");
    //    exit(1);
  } 
  */
  free(Poly);
  FreeBV(&Dummy);

  printf(" w= %3u  r=%3d  p= %3d  m1=%3d  m2=%3d  m3=%3d  wordno= %3d  hamingweight poly = %d\n",
   ((paramCARRY *)(((Generateur*)Gen)->ParamGen))->w,
   ((paramCARRY *)(((Generateur*)Gen)->ParamGen))->r,
   ((paramCARRY *)(((Generateur*)Gen)->ParamGen))->p,
   ((paramCARRY *)(((Generateur*)Gen)->ParamGen))->m1,
   ((paramCARRY *)(((Generateur*)Gen)->ParamGen))->m2,
   ((paramCARRY *)(((Generateur*)Gen)->ParamGen))->m3,
   ((paramCARRY *)(((Generateur*)Gen)->ParamGen))->wordno, count);
  cost=0;
  for (j=0;j<NBMAT;j++) {
    printf("A_%d = ",j);
    type = ((paramCARRY *)(((Generateur*)Gen)->ParamGen))->type[j];
    cost += ((paramCARRY *)(((Generateur*)Gen)->ParamGen))->Cost[type]();
    ((paramCARRY *)(((Generateur*)Gen)->ParamGen))->DisplayParams[type]
      (((paramCARRY *)(((Generateur*)Gen)->ParamGen))->paramsint[j],((paramCARRY *)(((Generateur*)Gen)->ParamGen))->paramsulong[j]);
    printf("\n");
  }
  printf("Cost = %d\n",cost);
  PrintUsableParameters(Gen);

}

static void PrintUsableParameters(Generateur *Gen) {
  int i,j;
  printf("%3d %3d %3d ",TM1,TM2,TM3);
  for (j=0;j<NBMAT;j++) {
    printf("%3d ",TTYPE(j));
    for (i=0;i<NBPARAMSINT;i++)
      printf("%3d ",PARINT(i,j));
    for (i=0;i<NBPARAMSULONG;i++)
      printf("%08x ",PARULONG(i,j));
  }
  printf("\n");
}
/* Initializes all parameters of a Carry generator component. */
void InitParamCarry ( Generateur *Gen, int w, int r, int p, int m1, int m2, int m3,  
          int paramsint[NBMAT][NBPARAMSINT],
          uint32_t paramsulong[NBMAT][NBPARAMSULONG],
          int *type,
          int wordno,int L) {
  int i,j;
  
  if (w!=WL) {
    printf("w!=32 is not supported yet!\n");
    exit(1);
  }
  
  ((paramCARRY *)(Gen->ParamGen))->w=w;
  ((paramCARRY *)(Gen->ParamGen))->r=r;
  ((paramCARRY *)(Gen->ParamGen))->p=p;
  if (p==0)
    ((paramCARRY *)(Gen->ParamGen))->maskp = 0x0U;
  else
    ((paramCARRY *)(Gen->ParamGen))->maskp = 0xFFFFFFFFU>>(WL-p);
  ((paramCARRY *)(Gen->ParamGen))->m1=m1;
  ((paramCARRY *)(Gen->ParamGen))->m2=m2;
  ((paramCARRY *)(Gen->ParamGen))->m3=m3;
  for (j=0;j<NBMAT;j++) {
    ((paramCARRY *)(Gen->ParamGen))->type[j]   = type[j];
    for (i=0;i<NBPARAMSINT;i++)      
      ((paramCARRY *)(Gen->ParamGen))->paramsint[j][i] = paramsint[j][i];
    for (i=0;i<NBPARAMSULONG;i++)
      ((paramCARRY *)(Gen->ParamGen))->paramsulong[j][i] = paramsulong[j][i];
  }
  
  ((paramCARRY *)(Gen->ParamGen))->wordno=wordno;
  ((paramCARRY *)(Gen->ParamGen))->Matrices[0] = type0;
  ((paramCARRY *)(Gen->ParamGen))->Matrices[1] = type1;
  ((paramCARRY *)(Gen->ParamGen))->Matrices[2] = type2;
  ((paramCARRY *)(Gen->ParamGen))->Matrices[3] = type3;
  ((paramCARRY *)(Gen->ParamGen))->Matrices[4] = type4;
  ((paramCARRY *)(Gen->ParamGen))->Matrices[5] = type5;
  ((paramCARRY *)(Gen->ParamGen))->Matrices[6] = type6;
  ((paramCARRY *)(Gen->ParamGen))->Matrices[7] = type7;

  ((paramCARRY *)(Gen->ParamGen))->DisplayParams[0] = displaytype0;
  ((paramCARRY *)(Gen->ParamGen))->DisplayParams[1] = displaytype1;
  ((paramCARRY *)(Gen->ParamGen))->DisplayParams[2] = displaytype2;
  ((paramCARRY *)(Gen->ParamGen))->DisplayParams[3] = displaytype3;
  ((paramCARRY *)(Gen->ParamGen))->DisplayParams[4] = displaytype4;
  ((paramCARRY *)(Gen->ParamGen))->DisplayParams[5] = displaytype5;
  ((paramCARRY *)(Gen->ParamGen))->DisplayParams[6] = displaytype6;
  ((paramCARRY *)(Gen->ParamGen))->DisplayParams[7] = displaytype7;

  ((paramCARRY *)(Gen->ParamGen))->Cost[0] = cost0;
  ((paramCARRY *)(Gen->ParamGen))->Cost[1] = cost1;
  ((paramCARRY *)(Gen->ParamGen))->Cost[2] = cost2;
  ((paramCARRY *)(Gen->ParamGen))->Cost[3] = cost3;
  ((paramCARRY *)(Gen->ParamGen))->Cost[4] = cost4;
  ((paramCARRY *)(Gen->ParamGen))->Cost[5] = cost5;
  ((paramCARRY *)(Gen->ParamGen))->Cost[6] = cost6;
  ((paramCARRY *)(Gen->ParamGen))->Cost[7] = cost7;

  Gen->k=w*r-p;
  Gen->L=L;
  Gen->Step=1;
  Gen->smax = INT_MAX;
  Gen->DispGen=DispCarry;
  Gen->DispName=CarryName; 

  Gen->InitGen=InitCarry32;
  Gen->Iteration=Carry32;
  Gen->CopyGen=CopyCarry;
  Gen->AllocGen = AllocCarry;
  Gen->PolyChar = polychar;
  PutBVToZero(&(Gen->GenState)); 

}

void CopyCarry(Generateur *Gen1, Generateur *Gen2) {
  int i,j;
  ((paramCARRY*)(Gen1->ParamGen))->w       = ((paramCARRY*)(Gen2->ParamGen))->w ;
  ((paramCARRY*)(Gen1->ParamGen))->r       = ((paramCARRY*)(Gen2->ParamGen))->r ; 
  ((paramCARRY*)(Gen1->ParamGen))->p       = ((paramCARRY*)(Gen2->ParamGen))->p ; 
  ((paramCARRY*)(Gen1->ParamGen))->m1      = ((paramCARRY*)(Gen2->ParamGen))->m1 ; 
  ((paramCARRY*)(Gen1->ParamGen))->m2      = ((paramCARRY*)(Gen2->ParamGen))->m2; 
  ((paramCARRY*)(Gen1->ParamGen))->m3      = ((paramCARRY*)(Gen2->ParamGen))->m3; 
  ((paramCARRY*)(Gen1->ParamGen))->maskp   = ((paramCARRY*)(Gen2->ParamGen))->maskp ; 
  for (j=0;j<NBMAT;j++) {
    ((paramCARRY *)(Gen1->ParamGen))->type[j]   = ((paramCARRY*)(Gen2->ParamGen))->type[j];
    for (i=0;i<NBPARAMSINT;i++)
      ((paramCARRY *)(Gen1->ParamGen))->paramsint[j][i] = ((paramCARRY*)(Gen2->ParamGen))->paramsint[j][i];
    for (i=0;i<NBPARAMSULONG;i++)
      ((paramCARRY *)(Gen1->ParamGen))->paramsulong[j][i] = ((paramCARRY*)(Gen2->ParamGen))->paramsulong[j][i];

  }
  ((paramCARRY*)(Gen1->ParamGen))->wordno = ((paramCARRY*)(Gen2->ParamGen))->wordno;

  for (j=0;j<NBTYPE;j++) {
    ((paramCARRY *)(Gen1->ParamGen))->Matrices[j]      =((paramCARRY *)(Gen2->ParamGen))->Matrices[j];
    ((paramCARRY *)(Gen1->ParamGen))->DisplayParams[j] =((paramCARRY *)(Gen2->ParamGen))->DisplayParams[j];
    ((paramCARRY *)(Gen1->ParamGen))->Cost[j]          =((paramCARRY *)(Gen2->ParamGen))->Cost[j];
  }
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

Generateur* AllocCarry ( int k ) {
  Generateur *G;
  G = (Generateur *) malloc (sizeof(Generateur));
  AllocGen(G,k);
  G->ParamGen= (paramCARRY *) malloc (sizeof(paramCARRY));
  return G;                      
}

void FreeCarry ( Generateur *G) {
  free(G->ParamGen);
  FreeGen(G);
  free(G);
}

#define TW ((paramCARRY *)(Gen->ParamGen))->w
#define TR ((paramCARRY *)(Gen->ParamGen))->r
#define TP ((paramCARRY *)(Gen->ParamGen))->p
#define WORDNO ((paramCARRY *)(Gen->ParamGen))->wordno
#define TI ((paramCARRY *)(Gen->ParamGen))->i
#define UPPER ((paramCARRY *)(Gen->ParamGen))->maskp
#define LOWER ~(((paramCARRY *)(Gen->ParamGen))->maskp)
#define TMAT(j,v) ((paramCARRY *)(Gen->ParamGen))->Matrices[TTYPE(j)](v,TPARINT(j),TPARULONG(j))
#define ETAT (Gen->GenState)
#define S -1 /* shift value in the transition matrix */
#define V0      ETAT.vect[TI]
#define V1      ETAT.vect[(TI+1)%TR]
#define V2      ETAT.vect[(TI+2)%TR]
#define VM1     ETAT.vect[(TI+TM1)%TR]
#define VM2     ETAT.vect[(TI+TM2)%TR]
#define VM3     ETAT.vect[(TI+TM3)%TR]
#define Vrm1    ETAT.vect[(TI+TR-1)%TR]
#define Vrm2    ETAT.vect[(TI+TR-2)%TR]
#define newV0   ETAT.vect[(TI+TR+S)%TR]
#define newV1   ETAT.vect[(TI+TR+S+1)%TR]
#define newVrm3 ETAT.vect[(TI+TR+S-3)%TR]
#define newVrm2 ETAT.vect[(TI+TR+S-2)%TR]
#define newVrm1 ETAT.vect[(TI+TR+S-1)%TR]
void InitCarry32(Generateur *Gen, BitVect *init, BitVect *retour)
{
  int j;
  CopyBV (&ETAT, init);
  TI = 0;
  for (j=0;j<Gen->L/32;j++)
    retour->vect[j] = ETAT.vect[(WORDNO+j)%(TR)];
}
/**********************************
La forme de la matrice de transition est

 T   T1  T2  T3  T4  T5  M0x 0x
 T   T   T7  T8  T9  T10 T   M2x
 I   0   0   0   0   0   0   0x          
 0   I   0   0   0   0   0   0x
 0   0   I   0   0   0   0   0x
 0   0   0   I   0   0   0   0x
 0   0   0   0   I   0   0   0x
 0   0   0   0       I   0   0x

La matrice M4: TMAT(4, )
La matrice M0: TMAT(0, )
Les matrices T1-T13 sont des compositions de matrices TMAT( ).
**************************************/
void Carry32(Generateur *Gen, BitVect *retour) {
    uint32_t z0,z1,z2,z3,z4;
  int j;

  z0 = (Vrm1 & LOWER)|(Vrm2 & UPPER) ;
  z1 = TMAT(0,V0)  ^ TMAT(1,VM1);
  z2 = TMAT(2,VM2) ^ TMAT(3,VM3);
  z3 = z1 ^ z2 ;
  z4 = TMAT(4,z0) ^ TMAT(5,z1)^  TMAT(6,z2)^  TMAT(7,z3);
  newV0 = z4;
  //newV0=z0;
  newV1 = z3;
  //  newVrm1&=0xff800000U; pas vraiment necessaire

  TI=(TI+TR+S) % (TR);
  for (j=0;j<=(Gen->L-1)/WL;j++)
    retour->vect[j] = ETAT.vect[(TI+j+WORDNO)%(TR)];
}

void ProduceParamsCarry( char *filename, int nb, int w, int r, int p, int maxcost, double seedmodifier) {

  int  i,j,count=0,*Poly,tried2,tried,m1,m2,m3;
  FILE *f;
  int paramsint[NBMAT][NBPARAMSINT];
  uint32_t paramsulong[NBMAT][NBPARAMSULONG];
  int type[NBMAT];
  Generateur *gen;
  BitVect BVPoly;
  int rand,cost,degre;
  void (*SelectParams[NBMAT])(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]);

  SelectParams[0] = selecttype0;
  SelectParams[1] = selecttype1;
  SelectParams[2] = selecttype2;
  SelectParams[3] = selecttype3;
  SelectParams[4] = selecttype4;
  SelectParams[5] = selecttype5;
  SelectParams[6] = selecttype6;
  SelectParams[7] = selecttype7;

  f=fopen(filename,"w"); 
  fprintf(f,"carry\n%d %d %d\n",w,r,p); 
  fprintf(f, "%d\n",nb); 
  if (WL - (WL/w)*w) {
    printf("The parameter w must be a divisor of 32.\n");
    exit(1);
  }
  degre = w*r-p;

  //   printf("Seed is NOT set by the clock\n");
  //    SetMRG32k3a(12345.0, 123123.0, 42232.0, 12345.0, 123123.0, 42232.0 );

  SetMRG32k3a((double)time(NULL)+seedmodifier,(double)time(NULL)+seedmodifier,(double)time(NULL)+seedmodifier,
        (double)time(NULL)+seedmodifier,(double)time(NULL)+seedmodifier,(double)time(NULL)+seedmodifier);

  Poly = (int*) malloc((degre+1)*sizeof(int));
  AllocBV(&BVPoly,degre+1);
  tried = tried2= 0;
  fclose(f);
  if ( (gen = AllocCarry(degre))==NULL) {
    printf("Error in ProduceParamsCarry()\n");
    exit(1);
  }
  while (count<nb) {
    cost=0;
    for (j=0;j<NBMAT;j++) {
      if (j==0 || j==4) {
  if (1) {
    rand =  (int)MRG32k3a() % 100;
    if (rand<16) {
      type[j] = 0;
      cost+=cost0();
    }
    else if (rand<32) {
      type[j] = 1;
      cost+=cost1();
    }
    else if (rand<48) {
      type[j] = 2;
      cost+=cost2();
    }
    else if (rand<64) {
      type[j] = 4;
      cost+=cost4();
    }
    else if (rand<80) {
      type[j] = 5; /* type 6 is excluded */
      cost+=cost5();
    }
    else {
      type[j]=5;
      cost+=cost5();
    }
  } else {
    type[j]=1;
    cost+=cost1();
  }
      }
      else{
  rand =  (int)MRG32k3a() % 100;
  if (rand<12) {
    type[j] = 0;
    cost+=cost0();
  }
  else if (rand<25) {
    type[j] = 1;
    cost+=cost1();
  }
  else if (rand<37) {
    type[j] =2;
    cost+=cost2();
  }
  else if (rand<50) {
    type[j] =3;
    cost+=cost3();
  }
  else if (rand<62) {
    type[j] =4;
    cost+=cost4();
  }
  else if (rand<75) {
    type[j] = 5;
    cost+=cost5();
  }
  else if (rand<87) {
    type[j] = 5;  //On exclut le type 6 !!!!!!!!!!!!!!!
    cost+=cost5();
  }
  else {
    type[j] = 7;
    cost+=cost7();
  }
      }
      if (cost>maxcost)
  break;
      SelectParams[type[j]](paramsint[j], paramsulong[j]);      
    }
    if (cost<maxcost) {
      m1=m2=m3=0;
      while ((m1 == 0) || (m1==1) || (m1==(r-1)))      m1 = MRG32k3a()%r;
      while ((m2 == 0) || (m2==1) || (m2==(r-1)))     m2 = MRG32k3a()%r;
      while ((m3 == 0) || (m3==1) || (m3==(r-1)))      m3 = MRG32k3a()%r;
      
      InitParamCarry(gen,w,r,p,m1,m2,m3,paramsint,paramsulong,type,0,WL);
      polychar(gen, Poly, &BVPoly);

      //      DispPolynomial(Poly,DEGGEN(gen));
      
      if (PrimitifPolynomialSieving(Poly,degre)) {
  f = fopen(filename,"a");
  fprintf(f,"%3d %3d %3d ",m1,m2,m3);
  for (j=0;j<NBMAT;j++) {
    fprintf(f,"%3d ",type[j]);
    for (i=0;i<NBPARAMSINT;i++)
      fprintf(f,"%3d ",paramsint[j][i]);
    for (i=0;i<NBPARAMSULONG;i++)
      fprintf(f,"%08x ",paramsulong[j][i]);
  }
  fprintf(f,"\n");
  count++;
  fclose(f);
      }
      tried2++;
    }
    //        else printf("%d\n",cost);
    tried++;

  }
  FreeCarry(gen);
  printf("tried=%d (%d)\n",tried,tried2);

  FreeBV(&BVPoly);
  free(Poly);
}

void ReadDataCarry ( Component *E,  char *filename, boolean same, int L ) {

  int n,i,j, nb;
  int w,r,p,m1,m2,m3,degre;
  int paramsint[NBMAT][NBPARAMSINT];
  uint32_t paramsulong[NBMAT][NBPARAMSULONG];
  int type[NBMAT];
  FILE *f;
  Generateur *g;
  
  f=fopen(filename,"r");
  if (f==NULL) {
    printf("File not found: %s\n",filename);
    exit(1);
  }
  ReadLn(f);
  fscanf(f, "%d", &w);          /*  Lecture de w  */
  fscanf(f, "%d", &r);          /*  Lecture de r  */    
  fscanf(f, "%d", &p);          /*  Lecture de r  */    
  degre = w*r-p;
  if (w>32) {
    printf("w doit etre plus petit que 32-bits");
    exit(1);
  }
  ReadLn(f) ;

  fscanf(f, "%d", &nb);
  ReadLn(f) ;
  if ( (g = AllocCarry(degre))==NULL) {
     printf("Error in ReadDataCarry()\n");
     exit(1);
  }
  AllocGensInComponent(E, nb, same);
  for (n = 0; n < nb; n++ ) {
    fscanf(f, "%d %d %d", &m1,&m2,&m3);
    for (j=0;j<NBMAT;j++) {
      fscanf(f,"%d",&type[j]);
      for (i=0;i<NBPARAMSINT;i++)
  fscanf(f,"%d",&paramsint[j][i]);
  for (i=0;i<NBPARAMSULONG;i++) {
    fscanf(f,"%x",&paramsulong[j][i]);
  }
    }
    ReadLn(f);
    for (j=0;j<1;j++) {
      InitParamCarry(g,w,r,p,m1,m2,m3,paramsint,paramsulong,type,j,L);
      AddGenInComponent(E, g);
    }
  }  
  fclose(f);
  FreeCarry(g);

  return;
  
}


/* type 0: Av = v ^ (v<<t) */
static uint32_t type0(uint32_t v, int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]) {
  if (paramsint[0]<0)
    return v^(v<<(-paramsint[0]));
  else
    return v^(v>>paramsint[0]);
}
/* type 1: Av = v */
static uint32_t type1(uint32_t v,  int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]) {
  return v;
}

/* type 2: A = companion(a) */
static uint32_t type2(uint32_t v,  int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]) {
  if (v & 1UL)
    return (v>>1) ^ paramsulong[0];
  else
    return (v>>1);
}
/* type 3: Av = (v<<t) */
static uint32_t type3(uint32_t v,  int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]) {
  if (paramsint[0]<0)
    return (v<<(-paramsint[0]));
  else
    return (v>>paramsint[0]);
}

// type 4: Av = v^ ((v<<t) & b)
static uint32_t type4(uint32_t v,  int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]) {
  if (paramsint[0]<0)
    return v^((v<<(-paramsint[0])) & paramsulong[0]);
  else
    return v^((v>>paramsint[0]) & paramsulong[0]);
}

//type 5: Av = (((v << t) ^ (v >> (32-t)))& b) (^c) selon une condition
static uint32_t type5(uint32_t v, int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]) {

  if (v & paramsulong[2]) {
    return (((v<<paramsint[0]) ^ (v>>(WL-paramsint[0]))) & paramsulong[1]) ^ paramsulong[0];
  }
  else{
    return (((v<<paramsint[0]) ^ (v>>(WL-paramsint[0]))) & paramsulong[1]);

  }

}
//type 6: three-shift marsaglia
static uint32_t type6 ( uint32_t v, int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]) {
  int s;
  uint32_t y;

  if ((s=paramsint[0])>0)
    y = v ^ (v>>s);
  else
    y = v ^ (v <<(-s));

  if ((s=paramsint[1])>0)
    y ^= (y>>s);
  else
    y ^= (y <<(-s));

  if ((s=paramsint[2])>0)
    y ^= (y>>s);
  else
    y ^= (y <<(-s));
  //  exit(1);
  return y;
}

// The ZERO matrix 
static uint32_t type7 (uint32_t v, int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]) {
  return 0U;
}



static void ResetParams(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]) {
  int j;
  for (j=0;j<NBPARAMSINT;j++)
    paramsint[j]=0;
  for (j=0;j<NBPARAMSULONG;j++)
    paramsulong[j]=0U;
}

static void selecttype0(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]) {
  ResetParams(paramsint,paramsulong);paramsint[0] = (MRG32k3a()%63)-31;
}

static void selecttype1(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]) {
   ResetParams(paramsint,paramsulong);
   return;
}

static void selecttype2(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]) {
   ResetParams(paramsint,paramsulong);
   paramsulong[0] = MRG32k3a() | 0x80000000U;
}

static void selecttype3(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]) {
   ResetParams(paramsint,paramsulong);
   paramsint[0] = (MRG32k3a()%63)-31;
}

static void selecttype4(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]) {
   ResetParams(paramsint,paramsulong);
   paramsint[0] = (MRG32k3a()%63)-31;
  paramsulong[0] = MRG32k3a();
}


static void selecttype5(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]) {
  int j,p,q,r,s,t;
  uint32_t ds,a,dt;
  ResetParams(paramsint,paramsulong);
  p = WL;
  while (gcd(p,WL)!=1)    p = MRG32k3a() % WL;
  q = MRG32k3a() % WL;

  r=1;
  while ( ((r*p)% WL) != 1) r++; /* r = p^-1 mod 32 */
  t = ((32-q)*r ) % WL;
  s = ((64-1-q)*r ) % WL;


  ds = ~(0x80000000U>>s);
  a = MRG32k3a() | (0x80000000U>>s);
  dt = 0x80000000U>>t;
  paramsint[0] = r;
  paramsulong[0] = a;
  paramsulong[1] = ds;
  paramsulong[2] = dt;
}
static void selecttype6(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]) {
  ResetParams(paramsint,paramsulong);
  paramsint[0] = (MRG32k3a()%63)-31;
  paramsint[1] = (MRG32k3a()%63)-31;
  paramsint[2] = (MRG32k3a()%63)-31;
}

static void selecttype7(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]) {
   ResetParams(paramsint,paramsulong);
  return;
}

static void displaytype0(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]) {
  printf("T0(%d)",paramsint[0]);
}
static void displaytype1(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]) {
  printf("Identity");
}
static void displaytype2(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]) {
  printf("T2(%08x)",paramsulong[0]);
}
static void displaytype3(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]) {
  printf("T3(%d)",paramsint[0]);
}
static void displaytype4(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]) {
  printf("T4(%d,%08x)",paramsint[0],paramsulong[0]);
}
static void displaytype5(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]) {
  int s,t,j,p,q;
  printf("T5(%d,%08x,%08x,%08x)",paramsint[0], paramsulong[0],paramsulong[1],paramsulong[2]);
  /*
  Not necessarily the correct values of p and q due to a historical bug;
  files created after July 8 do not have this problem.
  for (j=0;j<WL;j++)
    if ((0x80000000U >> j) & paramsulong[2])
      break;
  s=j;
  for (j=0;j<WL;j++)
    if ((0x80000000U >> j) & (~paramsulong[1]))
      break;
  t=j;
  p  = InverseModN(t-s,WL);
  q  = (WL*WL - t*p) % WL;
  printf(" p = %d  q = %d", p, q );
  */

}
static void displaytype6(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]) {
  printf("T6(%d,%d,%d)",paramsint[0],paramsint[1], paramsint[2]);
}

static void displaytype7(int paramsint[NBPARAMSINT], uint32_t paramsulong[NBPARAMSULONG]) {
  printf("ZERO");
}

static int cost0(void) {
  return 3;
}
static int cost1(void) {
  return 1;
}
static int cost2(void) {
  return 5;
}
static int cost3(void) {
  return 2;
}
static int cost4(void) {
  return 4;
}
static int cost5(void) {
  return 8;
}
static int cost6(void) {
  return 7;
}
static int cost7(void) {
  return 0;
}

void OCarry32(Generateur *Gen, BitVect *retour) {
  
  uint32_t z1,z2,z3,z4; 
  int j; 
 
  z1 = TMAT(0,V1)  ^ TMAT(1,VM1); 
  z2 = TMAT(2,VM2) ^ TMAT(3,VM3);  
  z3 = z1 ^ z2 ; 
  z4 = TMAT(4,V0) ^ TMAT(5,z1) ^ TMAT(6,z2) ^ TMAT(7,z3); 
  newV0 =  Vrm1; 
  newVrm2 = z3 ;
  newVrm1 = z4 ;

  TI=(TI+1)%(TR); 
  for (j=0;j<=(Gen->L-1)/WL;j++) 
    retour->vect[j] = ETAT.vect[(TI+j+WORDNO)%(TR)];
}
static void TransitionMatrixCarry(Generateur *Gen) {
  Matrix A;
  BitVect State,BC;
  int i,k,K;
  if (TW*TR < Gen->k) {
    printf("Cannot display the generating matrix properly.\n");
    printf("Resolution of the generator smaller than degree (L<k)\n");
    exit(1);
  }
  K= TW*TR-TP;
  AllocBV(&BC,K);
  AllocBV(&State,K);
  BVCanonic(&BC,0);
  AllocMat(&A,K,K,1);
  for (i=0;i<K;i++) {
    INITGEN(Gen,&BC,&State);
    BVRShift(&BC,&BC,1);
    ITERATION(Gen,&State);

    for (k=0;k<K;k++)
      PutBitBV(&(A.lignes[k][0]),i,ValBitBV(&State,k));
  }
  DispMatSpecial(&A,K,K);
  FreeBV(&BC);
  FreeBV(&State);
}

static void DispMatSpecial (Matrix *m , int l,   int kg) { 

   int ii, k;    
  
   for (ii=0;ii<kg;ii++) {
     printf("[");     
     for (k=0;k<=(l-1)/WL;k++)
       printf("%08x ",(m->lignes[ii][0]).vect[k]);
     printf ("]\n");
   }
   printf ("l=%d kg=%d\n",l,kg);
}
