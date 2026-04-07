#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <tuplets.h>
#include <mecf.h>
#include <timer.h>
#include <regpoly.h>

static int nextCombination(int *i, int dim, int range);
static void ResetCombination(void);

/* Allocate the arrays inside a paramTuplets structure for a criterion of depth d. */
void AllocTuplets(paramTuplets *param, int d) {
  param->tuph        = (int *)   malloc((d + 1) * sizeof(int));
  param->gap         = NULL;
  param->DELTA       = (float *) malloc((d + 1) * sizeof(float));
  param->pourcentage = (float *) malloc((d + 1) * sizeof(float));
}

/* Free the memory allocated inside a paramTuplets structure. */
void FreeTuplets(paramTuplets *param) {
  free(param->tuph);
}

/* Change the test type. */
void ChangeTest(paramTuplets *params, int test, int testtype) {
  params->test = test;
  params->testtype = testtype;
}

/* Update the acceptance threshold. */
void ChangeTreshold(paramTuplets *params, float treshold) {
  params->treshold = treshold;
}

/* Initialize the paramTuplets structure with dimension parameters, threshold, and test type. */
void InitParamTuplets(paramTuplets *params, int d, int *t, float treshold,
                      int tupletsverif, int test, int testtype) {
  int i, dim;

  if (tupletsverif) {
    params->test = test;
    params->testtype = testtype;
    if (params->gap == NULL)
      params->gap = (float *) malloc((t[1] + 1) * sizeof(float));

    params->tupd = d;
    for (i = 1; i <= d; i++) {
      params->tuph[i] = t[i];
      params->DELTA[i] = INT_MIN;
    }
    params->indice_max = params->tuph[1];
    for (dim = 2; dim <= params->tupd; dim++) {
      if (params->indice_max < params->tuph[dim])
        params->indice_max = params->tuph[dim];
    }
    if (testtype == MAX)
      params->treshold = treshold;
    else
      params->treshold = -treshold;
  }
  params->tupletsverified = FALSE;
  params->tupletsverif = tupletsverif;
}

/* Return TRUE if the generator passes the tuplets criterion;
   returns TRUE without checking if tupletsverif is FALSE. */
boolean isTuplets(paramTuplets *params) {
  if (params->tupletsverif == FALSE)
    return TRUE;
  else {
    if ((params->testtype == MAX) &&
        (params->treshold >= floatmax(params->firstpart_max, params->secondpart_max)))
      return TRUE;
    if ((params->testtype == SUM) &&
        (params->treshold >= (params->firstpart_sum + params->secondpart_sum)))
      return TRUE;
    return FALSE;
  }
}

/* Compute the Delta(h_1,...,h_d, G) EQUIDISTRIBUTION criterion for the current generator. */
void TestTuplets(Combinaison *C, paramTuplets *params) {
  int *indices;
  int nbposs, nbeczero, dim, j;
  double gap;
  double bound;
  double firstpartmax, secondpartmax;
  double firstpartsum, secondpartsum;
  int maxindsucc;
  static int firsttime = TRUE;
  static Matrix tMat;

  if (params->tupletsverif == FALSE)
    return;
  if (params->test != EQUIDISTRIBUTION) {
    fprintf(stderr, "TestTuplets: only EQUIDISTRIBUTION test is supported\n");
    exit(1);
  }

  if (!firsttime) {
    if ((tMat.nblignes != C->k_g) || (tMat.t != params->indice_max) || (tMat.l != C->L)) {
      FreeMat(&tMat);
      AllocMat(&tMat, C->k_g, C->L, params->indice_max);
    }
  } else {
    AllocMat(&tMat, C->k_g, C->L, params->indice_max);
    firsttime = FALSE;
  }
  PrepareMat(C, &tMat, params->indice_max);

  /* Computation of the first part of the criterion (successive dimensions) */
  firstpartmax = INT_MIN;
  firstpartsum = 0.0;
  indices = (int *) malloc(params->tuph[1] * sizeof(int));
  if (params->tuph[1] > C->k_g)
    maxindsucc = C->k_g;
  else
    maxindsucc = params->tuph[1];

  for (dim = 1; dim <= maxindsucc; dim++) {
    indices[dim - 1] = dim - 1;
    params->gap[dim] = (float) intmin(C->L, C->k_g / dim) -
                       ResolutionEquid(&tMat, C->k_g, dim, indices);

    firstpartsum += params->gap[dim];
    if (params->gap[dim] > firstpartmax)
      firstpartmax = params->gap[dim];

    if ((params->testtype == MAX) && (firstpartmax > params->treshold)) {
      firstpartmax = INT_MAX;
      break;
    }
    if ((params->testtype == SUM) && (firstpartsum > params->treshold)) {
      firstpartsum = INT_MAX;
      break;
    }
  }

  params->firstpart_max = firstpartmax;
  params->firstpart_sum = firstpartsum;
  free(indices);

  if (((params->testtype == MAX) && (firstpartmax <= params->treshold)) ||
      ((params->testtype == SUM) && (firstpartsum <= params->treshold))) {
    secondpartmax = INT_MIN;
    secondpartsum = 0.0;
    for (dim = 2; dim <= params->tupd; dim++) {
      params->DELTA[dim] = INT_MIN;
      indices = (int *) malloc(dim * sizeof(int));
      bound = (float) intmin(C->L, C->k_g / dim);
      nbposs = nbeczero = 0;
      ResetCombination();

      while (nextCombination(indices, dim, params->tuph[dim] - 1)) {
        gap = bound - ResolutionEquid(&tMat, C->k_g, dim, indices);
        if (gap > params->DELTA[dim])
          params->DELTA[dim] = gap;
        if (gap == 0.0)
          nbeczero++;
        else
          secondpartsum += gap;
        nbposs++;
        if ((params->testtype == MAX) && (params->DELTA[dim] > params->treshold)) {
          secondpartmax = params->DELTA[dim] = INT_MAX;
          break;
        }
        if ((params->testtype == SUM) && (secondpartsum > params->treshold)) {
          secondpartsum = params->DELTA[dim] = INT_MAX;
          break;
        }
      }
      params->pourcentage[dim] = (float) nbeczero / (float) nbposs;

      if (secondpartmax < params->DELTA[dim])
        secondpartmax = params->DELTA[dim];

      if (params->DELTA[dim] == INT_MAX)
        break;

      free(indices);
    }
  } else {
    secondpartmax = INT_MAX;
    secondpartsum = INT_MAX;
  }
  params->secondpart_max = secondpartmax;
  params->firstpart_max  = firstpartmax;
  params->secondpart_sum = secondpartsum;
  params->firstpart_sum  = firstpartsum;

  params->tupletsverified = TRUE;
}

/* Print the type of test done (which Delta criterion was computed). */
void DispTupletsTest(paramTuplets *params) {
  int j;
  printf("Test done: DELTA(");
  for (j = 1; j <= params->tupd - 1; j++)
    printf(" %d ,", params->tuph[j]);
  printf(" %d )\n", params->tuph[j]);
}

/* Print tables of the computed gaps. */
void DispTuplets(paramTuplets *params) {
  int i, j, T, maxD;

  if (params->tupletsverified == FALSE)
    return;

  printf("\n  Tables of gaps obtained (successive dimensions)");
  T = 1;
  maxD = params->tuph[1];
  printf("\n");
  while (maxD > 0) {
    printf("======");
    for (j = 1; j <= intmin(maxD, 16); j++)
      printf("+=====");
    printf("+\nDIM   ");
    for (j = 1; j <= intmin(maxD, 16); j++)
      printf("|%5d", T + j - 1);
    printf("|\n------");
    for (j = 1; j <= intmin(maxD, 16); j++)
      printf("+-----");
    printf("|\nGAP   ");
    for (j = 1; j <= intmin(maxD, 16); j++)
      if (params->gap[T + j - 1] == 0.0)
        printf("|     ");
      else
        printf("|%5d", (int) params->gap[T + j - 1]);
    printf("|\n======");
    for (j = 1; j <= intmin(maxD, 16); j++)
      printf("+=====");
    printf("+\n");
    maxD = maxD - 16;
    T = T + 16;
  }
  printf("\n");

  printf("\n\n   Tables of maximal gaps obtained (non-successive dimensions)\n");
  printf("===========+");
  for (i = 2; i <= params->tupd; i++)
    printf("======+");
  printf("\n");
  printf("Dimension  |");
  for (i = 2; i <= params->tupd; i++)
    printf("%6d|", i);
  printf("\n-----------+");
  for (i = 2; i <= params->tupd; i++)
    printf("------+");
  printf("\n");
  printf("t_i        |");
  for (i = 2; i <= params->tupd; i++)
    printf("%6d|", params->tuph[i]);
  printf("\n-----------+");
  for (i = 2; i <= params->tupd; i++)
    printf("------+");
  printf("\n");
  printf("GAP        |");
  for (i = 2; i <= params->tupd; i++)
    printf(" %5d|", (int) params->DELTA[i]);
  printf("\n-----------+");
  for (i = 2; i <= params->tupd; i++)
    printf("------+");
  printf("\n");
  printf("percentage |");
  for (i = 2; i <= params->tupd; i++)
    printf("%6.2f|", params->pourcentage[i] * 100);
  printf("\n===========+");
  for (i = 2; i <= params->tupd; i++)
    printf("======+");
  printf("\n");
  printf("Value of DELTA( ");
  for (i = 1; i < params->tupd; i++)
    printf("%d, ", params->tuph[i]);
  printf("%d ) = %5.3f\n", params->tuph[params->tupd],
         floatmax(params->firstpart_max, params->secondpart_max));
  printf("Sum of all gaps observed = %5.3f\n",
         params->firstpart_sum + params->secondpart_sum);
  fflush(stdout);
}

#define IMAX(j) (range - dim + j)

static int resetcombination = TRUE;

/* Advance the combination index array to the next combination in lexicographic order;
   return 0 when all combinations are exhausted. */
static int nextCombination(int *i, int dim, int range) {
  int j, k;

  if (resetcombination == TRUE) {
    resetcombination = FALSE;
    for (j = 0; j < dim; j++)
      i[j] = j;
  } else {
    i[dim - 1]++;
    for (j = dim - 1; j > 0; j--)
      if (i[j] > IMAX(j) + 1) {
        i[j - 1]++;
        for (k = j; k < dim; k++)
          i[k] = i[k - 1] + 1;
      }
    if (i[0] > 0)
      resetcombination = TRUE;
  }
  return !resetcombination;
}

/* Reset the combination iterator to its initial state. */
static void ResetCombination(void) {
  resetcombination = TRUE;
}
