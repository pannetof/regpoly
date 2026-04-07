#include <stdio.h>
#include <stdlib.h>
#include "tempering.h"
#include "regpoly.h"

/* Display all tempering transformations associated with Component E. */
void DispTrans(Component *E) {
  int i;

  for (i = 0; i < E->NbTrans; i++)
    E->Trans[i]->DispLinTrans((E->Trans[i]));
}

/* Apply all NbTrans tempering transformations in T to g_state in order, then copy the
   result into Out. */
void Transform(Transformation **T, int NbTrans, BitVect *g_state, BitVect *Out) {
  int i;
  for (i = 0; i < NbTrans; i++)
    T[i]->Trans(g_state, T[i]);
  for (i = 0; i < Out->n; i++)
    Out->vect[i] = g_state->vect[i];
}

/* Apply all NbTrans tempering transformations in T to g_state in reverse order (i.e.,
   the inverse tempering), then copy the result into Out. */
void TransformInverse(Transformation **T, int NbTrans, BitVect *g_state, BitVect *Out) {
  int i;
  for (i = NbTrans - 1; i >= 0; i--)
    T[i]->InverseTrans(g_state, T[i]);
  for (i = 0; i < Out->n; i++)
    Out->vect[i] = g_state->vect[i];
}

/* Update the tempering parameters for all components of the Combinaison C. */
void UpdateAllTrans(Combinaison *C) {
  int j;
  for (j = 0; j < C->J; j++)
    UpdateTrans(&(C->Components[j]));
}

/* Update the tempering parameters for Component E, adjusting the width w
   according to the current generator's resolution L. */
void UpdateTrans(Component *E) {
  int i;
  for (i = 0; i < E->NbTrans; i++) {
    if (E->Trans[i]->w_original == -1)
      E->Trans[i]->w = E->Gen[E->CurrentGen]->L;
    else
      E->Trans[i]->w = intmin(E->Trans[i]->w_original, E->Gen[E->CurrentGen]->L);
    E->Trans[i]->ChangeParamTrans((E->Trans[i]));
  }
}

/* Allocate space for nbtrans tempering transformations in Component E. */
void AllocTransInComponent(Component *E, int nbtrans) {
  if (E != NULL) {
    E->Trans = (Transformation **) calloc((size_t)nbtrans, sizeof(Transformation *));
    if (E->Trans == NULL) {
      printf("Error in AllocTransInComponent()\n");
      exit(1);
    }
    E->NbTransAlloc = nbtrans;
    E->NbTrans = 0;
  } else {
    printf("Error in AllocTransInComponent(): NULL pointer\n");
    exit(1);
  }
}

/* Free all tempering transformations stored in Component E and reset its counts. */
void FreeTransInComponent(Component *E) {
  int i;
  if (E != NULL) {
    if (E->Trans != NULL) {
      for (i = 0; i < E->NbTrans; i++)
        free(E->Trans[i]->ParamTrans);
      free(E->Trans);
      E->NbTransAlloc = 0;
      E->NbTrans = 0;
    } else {
      printf("Error in FreeTransInComponent(): Trans is NULL\n");
      exit(1);
    }
  } else {
    printf("Error in FreeTransInComponent(): NULL pointer\n");
    exit(1);
  }
}

/* Add a copy of the tempering transformation T to Component E. */
void AddTransInComponent(Component *E, Transformation *T) {
  if (E->NbTrans < E->NbTransAlloc) {
    if ((E->Trans[E->NbTrans] = T->AllocTrans()) == NULL) {
      printf("Error in AddTransInComponent()\n");
      exit(1);
    }
    T->CopyTrans(E->Trans[E->NbTrans], T);
    E->NbTrans++;
  } else {
    printf("Error in AddTransInComponent(): no space left\n");
    exit(1);
  }
}


