#include <stdlib.h>
#include <stdio.h>
#include "polynomials.h"
#include <zenfact.h>
#include "timer.h"
#include "myzen.h"
typedef struct {
  BigNum fact[100];
  BigNumLength fl[100];
  int multi[100];
  int nbfact;
} Facteurs;

//static void pow2minus1(BigNum *n, BigNumLength *nl, int degre);
static void MersenneNumFromFile(int degre, Facteurs *Fact);//,ZENRing *POLYRING,ZENRing *F2);
static int isPrimitive(ZENRing *POLYRING,/* ZENRing *F2,*/Facteurs *Fact, BigNumLength periodl);
static void FunctionF(ZENPoly *G, ZENPoly *F, ZENPoly *P2nm1, ZENRing *F2);

boolean ContainsMersenneExponent(int degre, int *poly, int Mexp){
  ZENRing F2,POLYRING;
  ZENElt ONE, COEFF;
  ZENPoly T, Tprime, F, P2n;
  ZENPrc prc;
  int j,n;
  boolean good;
  ConstructZENF2(&F2);
  ZENEltAlloc(ONE,F2);
  ZENEltSetToOne(ONE,F2);
  ZENPolyAlloc(F,degre,F2);
  ZENPolySetToXi(F,degre,F2);
  for(j=0;j<degre;j++){
    if(poly[j]){
      ZENPolySetCoeff(F,j,ONE,F2);
    }
  }
  if(ZENExtRingAlloc(POLYRING,F,F2)){
    printf("Error\n");
    exit(1);
  }      
  ZENPrcSetAll(prc);
  ZENRingAddPrc(POLYRING,prc);
  
  ZENPolyAlloc(P2n,degre,F2);
  ZENPolyAlloc(Tprime,degre,F2);
  ZENPolyAlloc(T,degre,F2);
  ZENPolySetToXi(P2n,2,F2);
  ZENPolyAssign(T,F,F2);
  n=1;
//  DispPolynomial(poly,degre);
  while(1){
//     printf("n=%d deg=%d\n",n,ZENPolyDeg(T,F2));
     if(ZENPolyDeg(T,F2) == Mexp){
        good =1;
        break;
     }
     if(ZENPolyDeg(T,F2)<Mexp){
        good = 0;
        break;
     }
     if(ZENPolyDeg(T,F2)-Mexp<n){
        good = 0;
        break;
     }
     FunctionF(&Tprime, &T, &P2n, &F2);
     ZENPolyAssign(T,Tprime,F2);
     n++;
  }
  if(good){
    printf("***************************************************\n");fflush(stdout);
    for(j=0;j<=Mexp;j++){
      ZENPolyGetCoeff(COEFF, T, j, F2);
      if(ZENEltIsOne(COEFF,F2))
         poly[j] = 1;
      else
         poly[j] = 0;
    }
    good = IrreduciblePolynomial(poly,Mexp);
  }
//  printf("n=%d deg=%d\n",n,ZENPolyDeg(T,F2));
  ZENEltFree(COEFF,F2);   
  ZENEltFree(ONE,F2);
  ZENPolyFree(P2n,F2);
  ZENPolyFree(F,F2);
  ZENPolyFree(T,F2);
  ZENPolyFree(Tprime,F2);

  if(ZENRingClose(POLYRING)){
    printf("Error\n");
    exit(1);
  }
  if(ZENRingClose(F2)){
    printf("Error\n");
    exit(1);
  }
//  printf("==>%d\n",good);
  return good;
}


static void FunctionF(ZENPoly *G, ZENPoly *F, ZENPoly *P2nm1, ZENRing *F2){
  ZENPoly P2n, P2ntemp;
  ZENPoly X,Q,PP,GCD,Dummy, Gtemp;
  int degre;

  degre = ZENPolyDeg(*F,*F2);
  ZENPolyAlloc(PP,degre,*F2);
  ZENPolyAlloc(GCD,degre,*F2);
  ZENPolyAlloc(Dummy,degre,*F2);
  ZENPolyAlloc(X,1,*F2);
  ZENPolyAlloc(P2n,2*degre,*F2);
  ZENPolyAlloc(P2ntemp,2*degre,*F2);
  ZENPolyAlloc(Q,degre,*F2);

  ZENPolyAssign(*G,*F,*F2);
  ZENPolySquare(P2ntemp,*P2nm1,*F2);  // X^{2^{n-1}}^2 = X^{2^n} 
  ZENPolyDivide(Dummy, P2n, P2ntemp, *F, *F2); // P2n = X^{2^n} mod F
  
  ZENPolyAssign(PP,P2n,*F2);
  ZENPolySetToXi(X,1,*F2);
  ZENPolyAdd(PP,X,*F2);               //  PP = X^{2^n} mod F + X
  ZENPolyGcd(Q, PP, *F,*F2);  // Q =gcd(F,X^{2^n} mod F + X)
  
   while(1){
    ZENPolyAlloc(Gtemp,degre,*F2);
    ZENPolyAssign(Gtemp,*G,*F2);
    ZENPolyGcd(GCD,Gtemp,Q,*F2);
    if(ZENPolyIsXi(GCD,0,*F2)){
      ZENPolyFree(Gtemp,*F2);
      break;
    }
    ZENPolyDivide(*G, Dummy, Gtemp, GCD, *F2); 
    ZENPolyFree(Gtemp,*F2);
  } 
  ZENPolyFree(PP,*F2);
  ZENPolyFree(GCD,*F2);
  ZENPolyFree(Dummy,*F2);
  ZENPolyFree(X,*F2);
  ZENPolyFree(P2n,*F2);
  ZENPolyFree(P2ntemp,*F2);
  ZENPolyFree(Q,*F2);
}



static void MersenneNumFromFile(int degre, Facteurs *Fact){//, ZENRing *POLYRING,ZENRing *F2){
        
  FILE *f;
  char filename[100];
  char r[25] = "/u/panneton/fact2^m-1/";
  
  char nbre_str[5000];
  int i=0;
  BigNum n;
  BigNumLength nl;

  sprintf(filename,"%s%d.txt",r,degre); 
  if(  (f=fopen(filename,"r"))==NULL){
    printf("Prime factor file not found: %s\n", filename);
    exit(1);
  }
  fscanf(f,"%s",nbre_str);
  ReadLn(f);
  while (!feof(f)){
    ZBNReadFromString(&n,&nl,nbre_str,10);
    if(i==0){ 
      Fact->fact[++i] = ZBNC(nl);
      ZBNAssign(Fact->fact[i],n,nl);
      ZBNF(n);
      Fact->fl[i]=nl;
      Fact->multi[i]=1;
    }
    else if(ZBNCompare(Fact->fact[i],Fact->fl[i],n,nl)!=0){
      Fact->fact[++i] = ZBNC(nl);
      ZBNAssign(Fact->fact[i],n,nl);
      ZBNF(n);
      Fact->fl[i]=nl;
      Fact->multi[i]=1;
    } 
    else{
      Fact->multi[i]++;
      ZBNF(n);
    }
  
    fscanf(f,"%s",nbre_str);
    ReadLn(f);
  }
  fclose(f);
  Fact->nbfact=i;
}

void DispPolynomial(int *coeff, int degre){

  int j;
  if(coeff[0])  
    printf("1+");  
  for(j=1;j<degre;j++)   
    if(coeff[j]) 
      printf("x^%d +",j);  
  if(coeff[degre]) 
    printf("x^%d\n\n",degre);
  fflush(stdout);
}

static int isPrimitive(ZENRing *POLYRING ,/*ZENRing *F2,*/ Facteurs *Fact, BigNumLength periodl){
  // Implemente l'algorithme d'Andreas Rieke
  BigNum n,mult,temp;
  BigNumLength nl,multl,templ;
  int j=1,primitive=1,i;
  ZENElt L[100],res,PX,ONE;

  if(Fact->nbfact>=100)
    printf("ERROR in isPrimitive()\n");
  multl=nl=templ=periodl;
  n=ZBNC(nl);
  mult=ZBNC(multl);
  ZBNSetToOne(n,nl);
  ZENEltAlloc(PX,*POLYRING);
  ZENEltSetToGenerator(PX,*POLYRING); 
  ZENEltAlloc(ONE,*POLYRING);
  ZENEltSetToOne(ONE,*POLYRING);
  //  printf("--------\n");
  //On multiplie tous les facteurs de multiplicit� >1 
  // i.e : {\prod_{k=1}^nbfact p_k^{r_k-1}
  for(j=1;j<=Fact->nbfact;j++){
    if(Fact->multi[j]>1){
      for(i=1;i<Fact->multi[j];i++){
	templ= nl+Fact->fl[j];
	temp = ZBNC(templ);
	ZBNMult(temp, n,nl, Fact->fact[j],Fact->fl[j]);
	ZBNF(n);
	nl=templ;
	n=ZBNC(nl);
	ZBNAssign(n,temp,templ);
	nl=ZBNNumDigits(n,nl);
	ZBNF(temp);
      }
    }
  }
  ZENEltAlloc(L[Fact->nbfact], *POLYRING);
  // x^{\prod_{k=1}^nbfact p_k^{r_k-1}
  ZENEltExp(L[Fact->nbfact],n,nl,PX,*POLYRING);
  for(j=Fact->nbfact-1;j>=1;j--){
    ZENEltAlloc(L[j],*POLYRING);
    ZENEltExp(L[j],Fact->fact[j+1],Fact->fl[j+1],L[j+1],*POLYRING);
  }
  ZENEltAlloc(res,*POLYRING);
  ZENEltExp(res,Fact->fact[1],Fact->fl[1],L[1],*POLYRING);
  j=Fact->nbfact;
  if(ZENEltAreEqual(res,ONE,*POLYRING))
    do{
      ZBNSetToOne(mult,multl);      
      for(i=j-1;i>0;i--){
	templ= multl+Fact->fl[j];
	temp = ZBNC(templ);
	ZBNMult(temp, mult,multl, Fact->fact[i],Fact->fl[i]);
	ZBNF(mult);
	multl=templ;
	mult=ZBNC(multl);
	ZBNAssign(mult,temp,templ);
	multl=ZBNNumDigits(mult,multl);
	ZBNF(temp);
      }
      ZENEltExp(res,mult,multl,L[j],*POLYRING);
      if(ZENEltAreEqual(res,ONE,*POLYRING)){
	primitive=0;
      }
      j--;   
    }while((j>=1) && primitive);
  else{
    primitive=0;
  }
  ZBNF(n);
  ZENEltFree(PX,*POLYRING);
  ZENEltFree(res,*POLYRING);
  ZENEltFree(ONE,*POLYRING);
  ZBNF(mult);
  for(j=1;j<=Fact->nbfact;j++)
    ZENEltFree(L[j],*POLYRING);
  return primitive;
}

#define VERBOSE 0
//Sieving propose dans
// @misc{ brent-fast,
//  author = "Richard P. Brent and Samuli Larvala and Paul Zimmermann",
//  title = "A Fast Algorithm For Testing Reducibility Of Trinomials Mod 2 And Some
//    New Primitive Trinomials Of Degree 3021377",
//  url = "citeseer.nj.nec.com/501076.html" }
static boolean Sieving(int degre, ZENRing *POLYRING, ZENRing *F2, double **ts, double *tf, int *r){
  ZENPoly P2n, P2ntemp;
  ZENPoly X,PP,PGCD,Dummy;
  int n;

  boolean reducible;
  timer_Chrono timer;
  double time;

  //  ZBNReadFromString(&four,&fourl,"4",10);
  if(degre!=*r && *r!=0){
    free(*ts);
    *ts = (double *) calloc (degre/2+1,sizeof(double));
    *r = degre;
  }
  else if(*r==0){
    *ts = (double *) calloc (degre/2+1,sizeof(double));
    *r = degre;
  }

  ZENPolyAlloc(P2n,degre,*F2);
  ZENPolyAlloc(X,1,*F2);
  ZENPolyAlloc(P2ntemp,2*degre,*F2);
  ZENPolyAlloc(PGCD,degre,*F2);
  ZENPolyAlloc(Dummy,degre,*F2);
  ZENPolySetToXi(P2n,4,*F2);
  ZENPolySetToXi(X,1,*F2);
  n=2;
  while(1){
    timer_Init(&timer); // Initialise le chronom�tre               
#if VERBOSE==1
    printf("n=%d  ts(n)=%f\n",n,(*ts)[n]);fflush(stdout);
#endif
    ZENPolyAlloc(PP,degre,*F2);
    ZENPolyAssign(PP,P2n,*F2);
    ZENPolyAdd(PP,X,*F2);
    ZENPolyGcd(PGCD,PP,ZENRingPol(*POLYRING),*F2);
    if(!ZENPolyIsXi(PGCD,0,*F2)){
      reducible=1;
      break;
    }
    ZENPolyFree(PP,*F2);
    ZENPolySquare(P2ntemp,P2n,*F2);
    ZENPolyDivide(Dummy,P2n,P2ntemp,ZENRingPol(*POLYRING),*F2);
    time = timer_Val(timer,timer_sec);
#if VERBOSE==1
    printf("time = %5.2f  %f\n",time,*tf);
#endif
    if(time > (*ts)[n] )  
      (*ts)[n]=time;
    if(n*((*ts)[n])>=*tf){
      reducible=0;
      break;
    }
    n++;      
    if(n>=degre/2){
      reducible=0;
      break;
    }
  } 
  ZENPolyFree(X,*F2);
  ZENPolyFree(Dummy,*F2);
  ZENPolyFree(PGCD,*F2);
  ZENPolyFree(P2n,*F2);
  ZENPolyFree(P2ntemp,*F2);
  return reducible;
}

int PrimitifPolynomial(int *tab, int degre){
  int j,ret,sum,reducible;
  ZENRing F2,POLYRING;
  ZENElt ONE;
  ZENPoly P;
  Facteurs Fact;

  if(tab[0]==0) return FALSE;
  for(sum=0,j=0;j<=degre;j++) sum += tab[j];
  if((sum%2)==0) return FALSE;
  ConstructZENF2(&F2);
  ZENEltAlloc(ONE,F2);
  ZENEltSetToOne(ONE,F2);
  ZENPolyAlloc(P,degre,F2);
  ZENPolySetToXi(P,degre,F2);
  for(j=0;j<degre;j++){
    if(tab[j]){
      ZENPolySetCoeff(P,j,ONE,F2);
    }
  }
  if(ZENExtRingAlloc(POLYRING,P,F2)){
    printf("Error\n");
    exit(1);
  }
  reducible = !IrreduciblePolynomial(tab,degre);

  ret = false;
  if(!reducible){
    MersenneNumFromFile(degre,&Fact);//,&POLYRING,&F2);
    ret = isPrimitive(&POLYRING,/*&F2,*/&Fact, degre/32+1);
    for(j=1;j<=Fact.nbfact;j++)
      ZBNF(Fact.fact[j]);
  }
  ZENEltFree(ONE,F2);
  ZENPolyFree(P,F2);
  if(ZENRingClose(POLYRING)){
    printf("Error\n");
    exit(1);
  }
  if(ZENRingClose(F2)){
    printf("Error\n");
    exit(1);
  }
  return ret;
}
int PrimitifPolynomialSieving(int *tab, int degre){
  static double *ts, tf=0.0; 
  static int r=0;
  int j,ret,sum,reducible;
  timer_Chrono timer;
  ZENRing F2,POLYRING;
  ZENElt ONE;
  ZENPoly P;
  Facteurs Fact;
  double time;
  ZENPrc prc;

  if(tab[0]==0) return FALSE; 
  for(sum=0,j=0;j<=degre;j++) sum += tab[j];
  if((sum%2)==0)  return FALSE;
  ConstructZENF2(&F2);
  ZENEltAlloc(ONE,F2);
  ZENEltSetToOne(ONE,F2);
  ZENPolyAlloc(P,degre,F2);
  ZENPolySetToXi(P,degre,F2);
  for(j=0;j<degre;j++){
    if(tab[j]){
      ZENPolySetCoeff(P,j,ONE,F2);
    }
  }
  if(ZENExtRingAlloc(POLYRING,P,F2)){
    printf("Error\n");
    exit(1);
  }
  ZENPrcSetAll(prc);
  ZENRingAddPrc(POLYRING,prc);
  timer_Init(&timer); // Initialise le chronom�tre 
  reducible =  Sieving(degre,&POLYRING,&F2, &ts, &tf, &r);

  if(!reducible){
    timer_Init(&timer); // Initialise le chronom�tre 
    reducible = !IrreduciblePolynomial(tab,degre);
    time = timer_Val(timer,timer_sec);
    if(time > tf )
      tf=time;
  }
  ret=false;
  if(!reducible){
    if(!ValidMersennePrime(degre)){
      MersenneNumFromFile(degre,&Fact);//,&POLYRING,&F2);
      ret = isPrimitive(&POLYRING,/*&F2,*/&Fact, degre/32+1);
      for(j=1;j<=Fact.nbfact;j++){
	ZBNF(Fact.fact[j]);
      }
      
    }else 
      ret=true;
  }
  ZENEltFree(ONE,F2);
  ZENPolyFree(P,F2);

  if(ZENRingClose(POLYRING)){
    printf("Error\n");
    exit(1);
  }
  if(ZENRingClose(F2)){
    printf("Error\n");
    exit(1);
  }
  return ret;

}
int IrreduciblePolynomial(int *tab, int degre){
  int j,ret;
  ZENRing F2;
  ZENElt ONE;
  ZENPoly P;
  ConstructZENF2(&F2);
  ZENEltAlloc(ONE,F2);
  ZENEltSetToOne(ONE,F2);
  ZENPolyAlloc(P,degre,F2);
  ZENPolySetToXi(P,degre,F2);
  for(j=0;j<degre;j++){
    if(tab[j])
      ZENPolySetCoeff(P,j,ONE,F2);
  }
  ret=(ZENPOLY_IS_IRREDUCIBLE == ZENPolyIsNotIrreducible(P,F2));
  ZENEltFree(ONE,F2);
  ZENPolyFree(P,F2); 
  if(ZENRingClose(F2)){
    printf("Error\n");
    exit(1);
  }
  return ret;

}
/* Computes the characteristic polynomial using the Massey-Berlekamp algorithm. */
void polychar ( Generateur *G, int *coeff, BitVect *BVPoly){
  int K,k;
  char *Sequence;
  BitVect BC,State;

  K=DEGGEN(G);
  AllocBV(&BC,K);
  AllocBV(&State,G->L);
  BVCanonic(&BC,0);
  Sequence = (char *) malloc (sizeof(char)*(2*K));
  INITGEN(G,&BC,&State);
  for(k=0;k<(2*K);k++){
    ITERATION(G,&State);
    Sequence[k] = ValBitBV(&State,0);
  }   
  SequenceMinimalPolynomial ( K, Sequence, coeff, BVPoly);
  FreeBV(&State);
  FreeBV(&BC);
  free(Sequence);
}

void SequenceMinimalPolynomial ( int K, char * Sequence, int *coeff, BitVect *BVPoly){
  int L=0,N=0,x=1,k,j,sum;
  BitVect B,C,Temp1,Temp2;
  AllocBV(&B,K+1);
  AllocBV(&C,K+1);
  AllocBV(&Temp1,K+1);
  AllocBV(&Temp2,K+1);
  BVCanonic(&B,0);
  BVCanonic(&C,0);
  while(N!=2*K){
    sum=Sequence[N];
    for(j=1;j<=L;j++)
      sum+=ValBitBV(&C,j)*Sequence[N-j];
    sum%=2;
    if(sum==0){
      x++;
    } else{
      //      printf("N=%d\n",N);
      if(2*L>N){
	BVRShift(&Temp1,&B,x);
	XORBV(&C,&C,&Temp1);
	x++;
      } else {
	CopyBV(&Temp2,&C);
	BVRShift(&Temp1,&B,x);
	XORBV(&C,&C,&Temp1);
	L=N+1-L;
	CopyBV(&B,&Temp2);
	x=1;
      }
    }
    N++;
  }

  for(k=0;k<=K;k++)
    PutBitBV(BVPoly,k,coeff[k]=ValBitBV(&C,K-k));
  FreeBV(&B);
  FreeBV(&C);
  FreeBV(&Temp1);
  FreeBV(&Temp2);
}

void polycharComb ( Combinaison *C, int *coeff,  BitVect *BVPoly ){
  ZENPoly Pj,P,temp;
  ZENRing F2;
  BigNum two;
  BigNumLength twol;
  ZENElt ONE,Coeff;
  BitVect Poly;
  int DegP=0,Degj,j,i;

  for(j=0;j<C->J;j++)
    DegP+=DEGGEN(GEN(C,j));
  AllocBV(&Poly,DegP+1);
  ZBNReadFromString(&two,&twol,"2",10);
  if(ZENBaseRingAlloc(F2,two,twol)){
    printf("Error\n");
    exit(1);
  }
  ZBNF(two);
  ZENEltAlloc(ONE,F2);
  ZENEltSetToOne(ONE,F2);
  ZENPolyAlloc(P,0,F2);
  ZENPolySetToXi(P,0,F2);
  DegP=0;
  for(j=0;j<C->J;j++){
    GEN(C,j)->PolyChar(GEN(C,j), coeff, &Poly);
    Degj=DEGGEN(GEN(C,j));
    ZENPolyAlloc(Pj,Degj,F2);
    ZENPolySetToXi(Pj,Degj,F2);
    for(i=0;i<Degj;i++)
      if(ValBitBV(&Poly,i))
	ZENPolySetCoeff(Pj,i,ONE,F2);
    ZENPolyAlloc(temp,DegP+=Degj,F2);
    ZENPolyMultiply(temp,Pj,P,F2);
    ZENPolyFree(P,F2);
    ZENPolyFree(Pj,F2);
    ZENPolyAlloc(P,DegP,F2);
    ZENPolyAssign(P,temp,F2);
    ZENPolyFree(temp,F2);
  }
  PutBVToZero(BVPoly);
  for(i=0;i<=DegP;i++){
    Coeff=ZENPolyGetCoeffPtr(P,i,F2);
    if(ZENEltIsOne(Coeff,F2)){
      PutBitBV(BVPoly,i,1);
      coeff[i]=1;
    }
    else{
      PutBitBV(BVPoly,i,0);
      coeff[i]=0;
    }
  } 
  ZENPolyFree(P,F2);
  ZENEltFree(ONE,F2);
  ZENRingClose(F2);
  FreeBV(&Poly);
}
void TransitionMatrix(Generateur *G){
  Matrix A;
  BitVect State,BC;
  int i,k,K;
  if(G->L < G->k){
    printf("Cannot display the generating matrix properly.\n");
    printf("Resolution of the generator smaller than degree (L<k)\n");
    exit(1);
  }
  K= G->k;
  AllocBV(&BC,K);
  AllocBV(&State,G->L);
  BVCanonic(&BC,0);
  AllocMat(&A,K,K,1);
  for(i=0;i<K;i++) {
    INITGEN(G,&BC,&State);
    BVRShift(&BC,&BC,1);
    ITERATION(G,&State);

    for(k=0;k<K;k++)
      PutBitBV(&(A.lignes[k][0]),i,ValBitBV(&State,k));
  }
  DispMat(&A,1,K,K);
  FreeMat(&A);
  FreeBV(&State);
  FreeBV(&BC);
}

void PowPuis2Poly(ZENPoly *a, ZENPoly *b, ZENRing *F2, int i){
  int j,puis=1, db;
  ZENElt Coeff,ONE;

  ZENEltAlloc(Coeff,*F2);
  ZENEltAlloc(ONE,*F2);
  ZENEltSetToOne(ONE,*F2);

  puis= 1UL << i;
  db = ZENPolyDeg(*b,*F2);
  ZENPolySetToXi(*a, db*puis, *F2);
  for(j=0;j<db;j++){
    ZENPolyGetCoeff(Coeff,*b,j,*F2);
    if(ZENEltAreEqual(Coeff,ONE,*F2))
      ZENPolySetCoeff(*a,puis*j,ONE,*F2);
  }
  ZENEltFree(Coeff,*F2);
  ZENEltFree(ONE,*F2);
}

void PowPoly(ZENPoly *a, ZENPoly *b, ZENRing *F2, int pow){
  int j,d;
  ZENPoly temp,temp2;
  if(pow==0){
    ZENPolySetToXi(*a,0,*F2);
    return;
  }
  d = ZENPolyDeg(*b,*F2);
  ZENPolySetToXi(*a, 0,*F2);
  ZENPolyAlloc(temp,d*pow,*F2);
  ZENPolyAlloc(temp2,d*pow,*F2);
  for(j=0;j<32;j++)
    if(pow &(1UL<<j)){
      PowPuis2Poly(&temp, b, F2,j);
      ZENPolyMultiply(temp2, *a, temp, *F2);
      *a = ZENPolyCopy(temp2,*F2);
    }
  ZENPolyFree(temp,*F2);
  ZENPolyFree(temp2,*F2);
}


