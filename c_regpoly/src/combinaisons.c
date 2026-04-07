#include "combinaisons.h"
#include <polynomials.h>
#include <limits.h>

static boolean ValidCurrent(Combinaison *C);
static boolean Update(Combinaison *C, int j);

/* Allocate storage for the state vector of a generator with k bits. */
void AllocGen(Generateur *Gen, int k) {
  AllocBV(&(Gen->GenState), k);
}

/* Free the state vector of a generator. */
void FreeGen(Generateur *Gen) {
  FreeBV(&(Gen->GenState));
}

/* Copy all parameters from Gen1 to Gen2 using the generator's own copy function. */
void CopyGen(Generateur *Gen1, Generateur *Gen2) {
  Gen2->CopyGen(Gen1, Gen2);
}

/* Add a copy of Gen into the Component E, allocating space for it. */
void AddGenInComponent(Component *E, Generateur *Gen) {
  if (E->NbGen < E->NbGenAlloc) {
    if ((E->Gen[E->NbGen] = Gen->AllocGen(Gen->k)) == NULL) {
      printf("Error in AddGenInComponent(). Memory Allocation did not work.\n");
      exit(1);
    }
    CopyGen((E->Gen[E->NbGen]), Gen);
    E->NbGen++;
  } else {
    printf("Error in AddGenInComponent(). Not enough space in the component for new generator\n");
    exit(1);
  }
}

/* Allocate space for nbJ components in the Combinaison C, setting the maximum resolution Lmax. */
void AllocComponentsInCombinaison(Combinaison *C, int nbJ, int Lmax) {
  C->J = nbJ;
  C->Lmax = Lmax;

  C->Components = (Component *) calloc((size_t)nbJ, sizeof(Component));
  if (C->Components == NULL) {
    printf("Error: out of memory in AllocComponentsInCombinaison()\n");
    exit(1);
  }
}

/* Free the components array of a Combinaison, including all generators and transformations. */
void FreeComponentsInCombinaison(Combinaison *C) {
  int j;

  if (C != NULL) {
    for (j = 0; j < C->J; j++) {
      FreeGensInComponent(&(C->Components[j]));
      FreeTransInComponent(&(C->Components[j]));
    }
    free(C->Components);
  }
}

/* Allocate space for nbgen generators in Component E, marking whether it is identical to the previous one. */
void AllocGensInComponent(Component *E, int nbgen, boolean same) {
  if (E != NULL) {
    E->Gen = (Generateur **) calloc((size_t)nbgen, sizeof(Generateur *));
    if (E->Gen == NULL) {
      printf("Error in AllocGensInComponent()\n");
      exit(1);
    }
    E->NbGenAlloc = nbgen;
    E->NbGen = 0;
    E->same = same;
  } else {
    printf("Error in AllocGensInComponent(): NULL pointer\n");
    exit(1);
  }
}

/* Free all generators stored in Component E and reset its counts. */
void FreeGensInComponent(Component *E) {
  int i;
  if (E != NULL) {
    if (E->Gen != NULL) {
      for (i = 0; i < E->NbGen; i++)
        free(E->Gen[i]->ParamGen); // TODO: handle BitVect members in params
      free(E->Gen);
    }
    E->NbGenAlloc = 0;
    E->NbGen = 0;
    E->same = 0;
  } else {
    printf("Error in FreeGensInComponent(): NULL pointer\n");
    exit(1);
  }
}

/* Set the current generator index in each component to the first valid combined generator.
   Returns TRUE on success, FALSE if no valid combination exists. */
boolean FirstGen(Combinaison *C) {
  int j, minL = INT_MAX;
  boolean OK = 0;

  C->k_g = 0;
  C->smax = INT_MAX;
  for (j = 0; j < C->J; j++)
    OK |= (C->Components[j].NbGen > 0);
  if (!OK)
    return FALSE;
  for (j = 0; j < C->J; j++) {
    if (C->Components[j].same)
      C->Components[j].CurrentGen = C->Components[j - 1].CurrentGen + 1;
    else
      C->Components[j].CurrentGen = 0;
  }
  OK = ValidCurrent(C);

  if (!OK)
    OK = NextGen(C);
  else {
    for (j = 0; j < C->J; j++) {
      C->k_g += C->Components[j].Gen[0]->k;

      if (C->Components[j].Gen[C->Components[j].CurrentGen]->smax < C->smax)
        C->smax = C->Components[j].Gen[C->Components[j].CurrentGen]->smax;

      if (C->Components[j].Gen[C->Components[j].CurrentGen]->L < minL)
        minL = C->Components[j].Gen[C->Components[j].CurrentGen]->L;
    }

    if (minL > C->Lmax)
      C->L = C->Lmax;
    else
      C->L = minL;
  }

  return OK;
}

/* Advance the current combined generator to the next valid combination.
   Returns TRUE on success, FALSE when all combinations have been exhausted. */
boolean NextGen(Combinaison *C) {
  int j, minL = INT_MAX;
  boolean UpdateOK, OK;

  do {
    OK = TRUE;
    UpdateOK = Update(C, C->J - 1);
    if (UpdateOK)
      OK = ValidCurrent(C);
  } while (!OK && UpdateOK);

  if (UpdateOK) {
    C->k_g = 0;
    C->smax = INT_MAX;
    for (j = 0; j < C->J; j++) {
      C->k_g += C->Components[j].Gen[C->Components[j].CurrentGen]->k;
      if (C->Components[j].Gen[C->Components[j].CurrentGen]->L < minL)
        minL = C->Components[j].Gen[C->Components[j].CurrentGen]->L;
      if (C->Components[j].Gen[C->Components[j].CurrentGen]->smax < C->smax)
        C->smax = C->Components[j].Gen[C->Components[j].CurrentGen]->smax;
    }
    if (minL > C->Lmax)
      C->L = C->Lmax;
    else
      C->L = minL;

    return TRUE;
  }
  return FALSE;
}

/* Check that the current combination is valid, i.e., no two components share the same period (k). */
static boolean ValidCurrent(Combinaison *C) {
  int j, jj;
  boolean OK = TRUE;

  for (j = 0; j < C->J; j++)
    for (jj = 0; jj < C->J; jj++) {
      if ((jj != j) &&
          (C->Components[j].Gen[C->Components[j].CurrentGen]->k ==
           C->Components[jj].Gen[C->Components[jj].CurrentGen]->k))
        OK = FALSE;
    }
  return OK;
}

/* Recursively increment the current generator index at position j, wrapping and
   propagating carries to earlier positions.  Returns FALSE when the last combination
   has been passed. */
static boolean Update(Combinaison *C, int j) {
  boolean OK;

  C->Components[j].CurrentGen++;

  /* Case 1: reached the upper bound of this component */
  if (C->Components[j].CurrentGen == C->Components[j].NbGen) {
    if (j != 0) {
      OK = Update(C, j - 1);
      if (C->Components[j].same)
        C->Components[j].CurrentGen = C->Components[j - 1].CurrentGen + 1;
      else
        C->Components[j].CurrentGen = 0;
      return (OK && TRUE);
    } else
      return FALSE;
  }

  /* Case 2: reached the effective upper bound because the next component is 'same' */
  if ((j < C->J - 1) && C->Components[j + 1].same &&
      (C->Components[j].CurrentGen == (C->Components[j].NbGen - C->J + j + 1))) {
    if (j != 0) {
      OK = Update(C, j - 1);
      if (C->Components[j].same)
        C->Components[j].CurrentGen = C->Components[j - 1].CurrentGen + 1;
      else
        C->Components[j].CurrentGen = 0;
      return (OK && TRUE);
    } else
      return FALSE;
  }

  return TRUE;
}

/* Display the parameters of the current combined generator (generators and tempering)
   to standard output, including Hamming weight of the characteristic polynomial. */
void DispCurrentComb(Combinaison *C) {
  int j;

  printf("==================================================================\n");
  printf("Number of points   : 2^(%3d)\n", C->k_g);
  if (C->smax != INT_MAX)
    printf("Maximum dimension : %3d\n", C->smax);
  printf("\n");

  {
    int *Poly;
    int count = 0;
    BitVect Dummy;
    Poly = (int *) malloc((DEGGEN(GEN(C, 0)) + 1) * sizeof(int));
    AllocBV(&Dummy, DEGGEN(GEN(C, 0)) + 1);

    polychar(GEN(C, 0), Poly, &Dummy);

    for (j = 0; j <= DEGGEN(GEN(C, 0)); j++) {
      if (Poly[j] == 1)
        count++;
    }
    printf("hammingweigth = %d\n", count);
    fflush(stdout);
    free(Poly);
    FreeBV(&Dummy);
  }

  for (j = 0; j < C->J; j++) {
    printf("%s:\n", C->Components[j].Gen[C->Components[j].CurrentGen]->DispName());
    C->Components[j].Gen[C->Components[j].CurrentGen]->DispGen(
        (C->Components[j].Gen[C->Components[j].CurrentGen]));
    DispTrans(&(C->Components[j]));
    if (j != C->J - 1)
      printf("------------------------------------------------------------------\n");
    else
      printf("==================================================================\n");
  }
}
