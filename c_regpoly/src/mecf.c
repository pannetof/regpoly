
#include <stdint.h>
#include <mecf.h>
#include <basisreduc.h>
#include <timer.h>
#include <regpoly.h>
#include <limits.h>

static Matrix Mat, mat;
static boolean firsttime = TRUE;

/* Mask for the most-significant bit, used to diagonalize the matrix in RangCF(). */
#define MC 0x80000000U

static int RangCF(Matrix *m, int kg, int t, int l);
static void SetPhi12(Combinaison *C, paramMecf *params);
void ConvEcarts(Combinaison *C, paramMecf *params);
static void SetPhi4(Combinaison *C, paramMecf *params);

extern uint32_t MMC[WL];

/* Set the equidistribution test method in the params structure. */
void SetMethod(paramMecf *params, int method) {
  if (method == METHOD_MATRICIAL)
    params->method = TestME;
  else if (method == METHOD_DUALLATTICE)
    params->method = TestMELat;
  else if (method == METHOD_NOTHING) {
    params->method = TestNothing;
    params->meverif = FALSE;
  } else {
    printf("Unknown method parameter!\n");
    exit(1);
  }
}

/* Allocate arrays in a paramMecf structure for testing up to resolution Lmax. */
void AllocMecf(paramMecf *params, int Lmax) {
  params->L = Lmax;
  params->delta = (int *) malloc((Lmax + 1) * sizeof(int));
  params->ecart = (int *) malloc((Lmax + 1) * sizeof(int));
  params->psi12 = NULL;
  params->phi12 = NULL;
  params->phi4  = NULL;
}

/* Free all dynamically allocated arrays inside a paramMecf structure. */
void FreeMecf(paramMecf *params) {
  free(params->delta);
  free(params->ecart);
  free(params->psi12);
  if (params->phi12 != NULL)
    free(params->phi12);
  if (params->phi4 != NULL)
    free(params->phi4);
}

/* Initialize a paramMecf structure with the given delta array and threshold values. */
void InitParamMecf(paramMecf *params,
                   int *delta,
                   int mse, int msecf, int meverif, int cfverif) {
  int i;
  params->meverif = meverif;
  params->meverified = FALSE;
  params->cfverif = cfverif;
  params->cfverified = FALSE;
  params->se = INT_MAX;
  params->secf = INT_MAX;
  params->msecf = msecf;
  params->mse = mse;
  for (i = 0; i <= params->L; i++)
    params->delta[i] = delta[i];
}

/* Return TRUE if the generator is presque-ME, i.e. quasi-ME and all per-resolution gaps
   are within the allowed delta bounds. */
boolean isPresqueME(paramMecf *params) {
  int l;
  boolean near = TRUE;
  if (params->meverif == FALSE)
    return TRUE;
  if (params->meverified == FALSE)
    return FALSE;
  for (l = 1; l <= params->L; l++)
    if (params->psi12[l]) {
      near = near && (params->ecart[l] <= params->delta[l]);
    }
  return near && (params->se <= params->mse);
}

/* Return TRUE if the generator is quasi-ME, i.e. the sum of gaps is within the mse bound. */
boolean isQuasiME(paramMecf *params) {
  if (params->meverified == FALSE)
    return FALSE;
  return (params->se <= params->mse);
}

/* Return TRUE if the generator is maximally equidistributed (all gaps are zero). */
boolean isME(paramMecf *params) {
  if (params->meverified == FALSE)
    return FALSE;
  return (params->se == 0);
}

/* Return TRUE if the generator is quasi-CF, i.e. ME and the sum of CF gaps is within msecf. */
boolean isQuasiCF(paramMecf *params) {
  if (params->cfverified == FALSE)
    return FALSE;
  return isME(params) && (params->secf <= params->msecf);
}

/* Return TRUE if the generator is collision-free (ME and all CF gaps are zero). */
boolean isCF(paramMecf *params) {
  if (params->cfverified == FALSE)
    return FALSE;
  return isME(params) && (params->secf == 0);
}

/* Do nothing; used as a placeholder when equidistribution testing is disabled. */
void TestNothing(Combinaison *C, paramMecf *params) {
  return;
}

/* Verify the equidistribution of the current generator in C; update params with the
   dimension gaps Delta_l for all l in Psi_12. */
void TestME(Combinaison *C, paramMecf *params) {
  int l, maxl;
  boolean verif;

  if (params->meverif == FALSE)
    return;

  if (params->L < C->L) {
    printf("Error in TestME().\n  Resolution of allocated paramMecf < resolution of generator\n");
    exit(1);
  }

  if (firsttime) {
    AllocMat(&Mat, C->k_g, params->L, intmin(C->k_g, C->smax));
    AllocMat(&mat, C->k_g, params->L, intmin(C->k_g, C->smax));
    firsttime = FALSE;
  } else {
    if ((Mat.nblignes < C->k_g) || (Mat.t < C->smax)) {
      if (Mat.nblignes != 0) {
        FreeMat(&Mat);
        FreeMat(&mat);
      }
      AllocMat(&Mat, C->k_g, params->L, intmin(C->k_g, C->smax));
      AllocMat(&mat, C->k_g, params->L, intmin(C->k_g, C->smax));
    }
  }

  for (l = 1; l <= params->L; l++) /* Reset the ecart array. */
    params->ecart[l] = -1;

  PrepareMat(C, &Mat, -1);

  params->se = 0;

  SetPsi12(C, params);

  verif = FALSE;
  maxl = params->L;
  for (l = 1; l <= params->L; l++) {
    if (params->ecart[l] == -1) { // Check if not already computed
      if (params->psi12[l] || verif) { /* IF ( l IN psi12 ) */
        timer_Chrono timer;

        CopyMat(&mat, &Mat, C->k_g, intmin(C->k_g / l, C->smax));
#define TIMING 0
        if (TIMING)
          timer_Init(&timer);
        params->ecart[l] = intmin(C->k_g / l, C->smax) - DimensionEquid(&mat, C->k_g, l, C->smax);
        if (TIMING) {
          printf("%d %f\n", l, timer_Val(timer, timer_sec));
          fflush(stdout);
        }

        params->se += params->ecart[l];
        if (params->ecart[l] > params->delta[l] || params->se > params->mse) {
          maxl = l;
          break;
        }

        if (params->ecart[l] != 0) { // If gap is non-zero, verify l-1 if not already done
          verif = TRUE;
          if (l != 1)
            l = l - 2;
        } else {
          verif = FALSE;
        }
      }
    }
  }

  params->se = 0;
  for (l = 1; l <= maxl; l++) {
    if (params->ecart[l] == -1)
      params->ecart[l] = 0;
    params->se += params->ecart[l];
  }
  for (; l <= params->L; l++)
    if (params->ecart[l] == -1)
      params->ecart[l] = INT_MAX;

  params->meverified = TRUE;
}

/* Compute t_l (the maximum equidistribution dimension at resolution l) by Gaussian
   elimination on the generator matrix; return the rank found. */
int DimensionEquid(Matrix *m, int kg, int l, int smax) {
  int i, j, cl, rang;
  int cldiv, clmod;
  int t;

  t = intmin(kg / l, smax);
  rang = 0;

  /* Diagonalize the matrix on entries (i,j) over l bits,
     with 0 <= i < kg and 0 <= j < t. */
  for (j = 0; j < t; j++) {
    for (cl = 0; cl < l; cl++) {
      cldiv = cl / WL;
      clmod = cl % WL;
      /* Search column j starting at row rang for the first entry
         whose most-significant bit is non-zero (pivot search). */
      i = rang;
      while ((i < kg) && (m->lignes[i][j].vect[cldiv] < MMC[clmod]))
        i++;
      if (i < kg) { /* pivot found */
        ExchangeVect(m, rang, i);
        for (i = rang + 1; i < kg; i++) {
          if (m->lignes[i][j].vect[cldiv] & MMC[clmod]) {
            XorVect(m, i, rang, j, t); // faster when t is small (large l)
          }
        }
        rang++;
      } else { /* no pivot found => matrix not full rank */
        return j;
      }
    }
  }
  return j; /* all pivots found => full rank */
}

/* Compute l_t (the maximum equidistribution resolution at dimension t) by Gaussian
   elimination; return the rank found. */
int ResolutionEquid(Matrix *m, int kg, int t, int *indices) {
  int i, j, cl, rang;
  int cldiv, clmod;
  int l;

  l = intmin(kg / t, m->l);
  rang = 0;

  /* Diagonalize the matrix on entries (i,j) over l bits,
     with 0 <= i < kg and 0 <= j < t. */
  for (cl = 0; cl < l; cl++) {
    for (j = 0; j < t; j++) {
      cldiv = cl / WL;
      clmod = cl % WL;
      /* Search column j starting at row rang for the first entry
         whose most-significant bit is non-zero (pivot search). */
      i = rang;
      while ((i < kg) && (m->lignes[i][indices[j]].vect[cldiv] < MMC[clmod]))
        i++;
      if (i < kg) { /* pivot found */
        ExchangeVect(m, rang, i);
        for (i = rang + 1; i < kg; i++) {
          if (m->lignes[i][indices[j]].vect[cldiv] & MMC[clmod])
            XorVect(m, i, rang, 0, m->t);
        }
        rang++;
      } else { /* no pivot found => matrix not full rank */
        return cl;
      }
    }
  }
  return cl; /* all pivots found => full rank */
}

/* Verify the collision-free property of the current generator in C; compute the
   CF gaps epsilon_t for all t in Phi_4 and update params accordingly. */
void TestCF(Combinaison *C, paramMecf *params) {
  int t, l;

  if (params->cfverif == FALSE)
    return;

  if (firsttime) {
    AllocMat(&Mat, C->k_g, params->L, C->k_g);
    AllocMat(&mat, C->k_g, params->L, C->k_g);
    firsttime = FALSE;
  } else {
    if (mat.nblignes < C->k_g) {
      if (mat.nblignes != 0) {
        FreeMat(&Mat);
        FreeMat(&mat);
      }
      AllocMat(&Mat, C->k_g, params->L, C->k_g);
      AllocMat(&mat, C->k_g, params->L, C->k_g);
    }
  }

  if (params->ecartCF != NULL)
    free(params->ecartCF);
  params->ecartCF = (int *) malloc((C->k_g + 1) * sizeof(int));
  for (l = 1; l <= C->k_g; l++) /* Reset the ecartCF array. */
    params->ecartCF[l] = -1;

  params->secf = 0;
  t = C->k_g;

  while (t > 1) {
    if (params->phi4[t]) { /* IF ( t IN phi4 ) */
      l = C->k_g / t; /* l associated with the current t */
      CopyMat(&mat, &Mat, C->k_g, t);
      params->secf += (params->ecartCF[l] = C->k_g - RangCF(&mat, C->k_g, t, l + 1));
    }
    t--;
  }
  params->cfverified = TRUE;
}

/* Build the generator matrix B~_{t,l} for the current generator in C; used by
   TestME() and related functions. */
void PrepareMat(Combinaison *C, Matrix *PMat, int maxdim) {
  BitVect State, BC;
  int gl, ml = 0, mc, j, indice_max;

  if (maxdim == -1)
    indice_max = intmin(C->smax, C->k_g);
  else
    indice_max = maxdim;

  for (j = 0; j < C->J; j++) {
    AllocBV(&BC, DEGGEN(GEN(C, j)));
    AllocBV(&State, C->Components[j].Gen[C->Components[j].CurrentGen]->L);
    BVCanonic(&BC, 0);
    for (gl = 0; gl < C->Components[j].Gen[C->Components[j].CurrentGen]->k; gl++, ml++) {
      INITGEN(GEN(C, j), &BC, &State);
      BVRShiftSelf(&BC, 1);
      TRANSFORME(C, j, &State, &((PMat->lignes[ml])[0]));
      for (mc = 1; mc < indice_max; mc++) {
        ITERATION(GEN(C, j), &State);
        TRANSFORME(C, j, &State, &((PMat->lignes[ml])[mc]));
      }
    }
    FreeBV(&BC);
    FreeBV(&State);
  }
}

/* Build the generator matrix with a rotative left-shift of rotative_shift bits applied
   to each output before storing it in PMat. */
void PrepareMat2(Combinaison *C, Matrix *PMat, int maxdim, int rotative_shift) {
  BitVect State, BC, Temp;
  int gl, ml = 0, mc, j, indice_max;

  if (maxdim == -1)
    indice_max = intmin(C->smax, C->k_g);
  else
    indice_max = maxdim;

  for (j = 0; j < C->J; j++) {
    AllocBV(&BC, DEGGEN(GEN(C, j)));
    AllocBV(&State, C->Components[j].Gen[C->Components[j].CurrentGen]->L);
    AllocBV(&Temp, C->Components[j].Gen[C->Components[j].CurrentGen]->L);
    BVCanonic(&BC, 0);
    for (gl = 0; gl < C->Components[j].Gen[C->Components[j].CurrentGen]->k; gl++, ml++) {
      INITGEN(GEN(C, j), &BC, &State);
      BVRShiftSelf(&BC, 1);
      TRANSFORME(C, j, &State, &Temp);
      BVLRotativeShift(&((PMat->lignes[ml])[0]), &Temp, rotative_shift,
                       C->Components[j].Gen[C->Components[j].CurrentGen]->L);
      for (mc = 1; mc < indice_max; mc++) {
        ITERATION(GEN(C, j), &State);
        TRANSFORME(C, j, &State, &Temp);
        BVLRotativeShift(&((PMat->lignes[ml])[mc]), &Temp, rotative_shift,
                         C->Components[j].Gen[C->Components[j].CurrentGen]->L);
      }
    }
    FreeBV(&BC);
    FreeBV(&Temp);
    FreeBV(&State);
  }
}

/* Display a table of dimension or resolution gaps for the current generator; return the
   sum of gaps over Psi_12 (type='l') or Phi_12 (type='t'). */
int DispTable(Combinaison *C, paramMecf *params, char type) {
  int i, j, max, somme, smax;
  int *table;

  smax = intmin(C->k_g, C->smax);
  if (params->meverified == FALSE)
    return -1;

  if (type == 'l') {
    max = intmin(C->L, params->L);
    table = (int *) malloc((params->L + 1) * sizeof(int));
    for (i = 0; i <= params->L; i++)
      table[i] = params->ecart[i];
  } else if (type == 't') {
    SetPhi12(C, params);
    ConvEcarts(C, params);
    max = smax;
    table = (int *) malloc((smax + 1) * sizeof(int));
    for (i = 0; i <= smax; i++)
      table[i] = params->Lambda[i];
  } else {
    printf("Error in DispTable().  Display type is not valid\n");
    exit(1);
  }

  for (j = 0; j <= (max - 1) / 16; j++) {
    printf("\n=======");
    for (i = j * 16 + 1; i <= (j + 1) * 16 && i <= max; i++)
      printf("+=====");
    if (type == 'l')
      printf("+\nRESOL  ");
    else
      printf("+\nDIM    ");
    for (i = j * 16 + 1; i <= (j + 1) * 16 && i <= max; i++)
      printf("|%5d", i);
    printf("|\n-------");
    for (i = j * 16 + 1; i <= (j + 1) * 16 && i <= max; i++)
      printf("+-----");
    printf("|\nECART  ");
    for (i = j * 16 + 1; i <= (j + 1) * 16 && i <= max; i++)
      if (table[i] == 0) {
        printf("|     ");
      } else {
        printf("|%5d", table[i]);
      }

    printf("|\n-------");
    for (i = j * 16 + 1; i <= (j + 1) * 16 && i <= max; i++)
      printf("+-----");
    if (type == 'l')
      printf("|\nDIM    ");
    else
      printf("|\nRESOL  ");
    for (i = j * 16 + 1; i <= (j + 1) * 16 && i <= max; i++) {
      if (type == 'l') {
        printf("|%5d", intmin(smax, C->k_g / i) - table[i]);
      } else {
        printf("|%5d", intmin(params->L, C->k_g / i) - table[i]);
      }
    }
    printf("|\n=======");
    for (i = j * 16 + 1; i <= (j + 1) * 16 && i <= max; i++)
      printf("+=====");
  }
  printf("+\n");

  somme = 0;
  if (type == 't') {
    for (i = 1; i <= smax; i++)
      if (params->phi12[i])
        somme += table[i];
  } else {
    somme = params->se;
  }

  if (type == 'l')
    printf("--------------------------->DIMENSION GAPS SUM (Psi_12) = %d\n", somme);
  else
    printf("--------------------------->RESOLUTION GAPS SUM (Phi_12) = %d\n", somme);

  return somme;
}

/* Print a message if the generator is maximally equidistributed. */
void DispME(paramMecf *params) {
  if (isME(params))
    printf("\n ===> ME GENERATOR\n");
}

/* Print a message if the generator is collision-free. */
void DispCF(paramMecf *params) {
  if (isCF(params))
    printf("\n ===> CF GENERATOR\n");
}

/* Convert the dimension gaps Delta_l into resolution gaps Lambda_t and store them in
   params->Lambda. */
void ConvEcarts(Combinaison *C, paramMecf *params) {
  int t, l, i, t_i, smax;

  smax = intmin(C->k_g, C->smax);
  if (params->Lambda != NULL)
    free(params->Lambda);
  params->Lambda = (int *) malloc((smax + 1) * sizeof(int));

  for (t = 1; t <= smax; t++) {
    params->Lambda[t] = -1;
    l = C->k_g / t;
    if (l > params->L)
      l = params->L;

    for (i = 1; i <= l; i++) {
      t_i = (intmin((C->k_g) / i, smax) - params->ecart[i]);
      if (t <= t_i)
        params->Lambda[t] = l - i;
    }
  }
}

/* Compute the set Phi_12 (values of t for which the resolution gap must be checked
   for the ME property) and store it in params->phi12. */
static void SetPhi12(Combinaison *C, paramMecf *params) {
  register int t, r, l, m, smax;

  smax = intmin(C->k_g, C->smax);

  if (params->phi12 != NULL)
    free(params->phi12);
  params->phi12 = (boolean *) malloc((smax + 1) * sizeof(boolean));

  for (t = 0; t <= smax; t++)
    params->phi12[t] = FALSE;

  /* PHI 1 part */
  r = (int) floor(sqrt((double) C->k_g));

  m = C->k_g / params->L;
  if (m < 2)
    m = 2;
  if (m > smax)
    m = smax;

  if (r > smax)
    r = smax;

  for (t = m; t <= r; t++)
    params->phi12[t] = TRUE;

  /* PHI 2 part */
  r = (int) floor(sqrt((double)(C->k_g - 1)));

  for (l = 1; l <= r; l++)
    params->phi12[intmin(smax, C->k_g / l)] = TRUE;
}

/* Compute the set Phi_4 (values of t for which the CF property must be checked)
   and store it in params->phi4. */
static void SetPhi4(Combinaison *C, paramMecf *params) {
  register int t, lt;

  if (params->phi4 != NULL)
    free(params->phi4);
  params->phi4 = (boolean *) malloc((C->k_g + 1) * sizeof(boolean));

  for (t = 0; t <= C->k_g; t++)
    params->phi4[t] = FALSE;
  for (t = 2; t <= C->k_g; t++)
    if (C->k_g % t) {
      lt = C->k_g / t;
      if (lt < params->L && C->k_g % (lt + 1) && C->k_g / (t - 1) > lt) {
        params->phi4[t] = TRUE;
      }
    }
}

/* Compute the set Psi_12 (values of l for which the dimension gap must be checked
   for the ME property) and store it in params->psi12. */
void SetPsi12(Combinaison *C, paramMecf *params) {
  int l, r, m, t;

  if (params->psi12 != NULL)
    free(params->psi12);
  params->psi12 = (boolean *) malloc((params->L + 1) * sizeof(boolean));
  for (l = 0; l <= params->L; l++)
    params->psi12[l] = FALSE;

  r = (int) floor(sqrt((float) C->k_g));
  for (l = 1; l <= intmin(r, params->L); l++)
    params->psi12[l] = TRUE;

  r = (int) floor(sqrt((float) C->k_g - 1));
  m = C->k_g / params->L;

  if (m < 2)
    m = 2;

  for (t = m; t <= r; t++)
    params->psi12[intmin(C->k_g / t, params->L)] = TRUE;
}

/* Compute the rank of matrix m (kg rows, t columns over l bits each) by Gaussian
   elimination; used for the collision-free property check. */
static int RangCF(Matrix *m, int kg, int t, int l) {
  int i, j, cl, rang;

  rang = 0;

  /* Diagonalize the matrix on entries (i,j) over l bits,
     with 0 <= i < kg and 0 <= j < t. */
  for (j = 0; j < t; j++) {
    cl = 1;
    while (cl <= l) {
      if (cl > 1) {
        /* Left-shift by one bit for the whole current column. */
        if (l > 32)
          for (i = 0; i < kg; i++)
            BVLS1Self(&((m->lignes[i])[j]));
        else
          for (i = 0; i < kg; i++)
            ((m->lignes[i])[j]).vect[0] <<= 1;
      }
      /* Search column j starting at row rang for the first entry
         whose most-significant bit is non-zero (pivot search). */
      i = rang;
      while ((i < kg) && (((m->lignes[i])[j]).vect[0] < MC))
        i++;
      if (i < kg) { /* pivot found */
        ExchangeVect(m, rang, i);
        for (i = rang + 1; i < kg; i++)
          if (((m->lignes[i])[j]).vect[0] & MC)
            XorVect(m, i, rang, j, t);
        rang++;
        if (rang == kg)
          return rang; /* full rank reached */
      }
      cl++;
    }
  }
  return rang;
}

/* Verify equidistribution of the current generator after applying a rotative left-shift
   of rotative_shift bits to the output; update params accordingly. */
void TestMESpecial(Combinaison *C, paramMecf *params, int rotative_shift) {
  int l, maxl;
  boolean verif;

  if (params->meverif == FALSE)
    return;

  if (params->L < C->L) {
    printf("Error in TestME().\n  Resolution of allocated paramMecf < resolution of generator\n");
    exit(1);
  }

  if (firsttime) {
    AllocMat(&Mat, C->k_g, params->L, intmin(C->k_g, C->smax));
    AllocMat(&mat, C->k_g, params->L, intmin(C->k_g, C->smax));
    firsttime = FALSE;
  } else {
    if ((Mat.nblignes < C->k_g) || (Mat.t < C->smax)) {
      if (Mat.nblignes != 0) {
        FreeMat(&Mat);
        FreeMat(&mat);
      }
      AllocMat(&Mat, C->k_g, params->L, intmin(C->k_g, C->smax));
      AllocMat(&mat, C->k_g, params->L, intmin(C->k_g, C->smax));
    }
  }

  for (l = 1; l <= params->L; l++) /* Reset the ecart array. */
    params->ecart[l] = -1;

  PrepareMat2(C, &Mat, -1, rotative_shift);

  params->se = 0;

  SetPsi12(C, params);

  verif = FALSE;
  maxl = params->L;
  for (l = 1; l <= params->L; l++) {
    if (params->ecart[l] == -1) { // Check if not already computed
      if (params->psi12[l] || verif) { /* IF ( l IN psi12 ) */
        CopyMat(&mat, &Mat, C->k_g, intmin(C->k_g / l, C->smax));
        params->ecart[l] = intmin(C->k_g / l, C->smax) - DimensionEquid(&mat, C->k_g, l, C->smax);
        params->se += params->ecart[l];

        if (params->ecart[l] != 0) { // If gap is non-zero, verify l-1 if not already done
          verif = TRUE;
          if (l != 1)
            l = l - 2;
        } else {
          verif = FALSE;
        }
      }
    }
  }

  params->se = 0;
  for (l = 1; l <= maxl; l++) {
    if (params->ecart[l] == -1)
      params->ecart[l] = 0;
    params->se += params->ecart[l];
  }
  for (; l <= params->L; l++)
    if (params->ecart[l] == -1)
      params->ecart[l] = INT_MAX;

  params->meverified = TRUE;
}
