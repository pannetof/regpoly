#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <readfiles.h>
#include <regpoly.h>
#include <tausworthe.h>
#include <polylcg.h>
#include <tgfsr.h>
#include <genf2w.h>
#include <MT.h>
#include <marsaxorshift.h>
#include <matsumoto.h>
#include <AC1D.h>
#include <carry2.h>
#include <marsa.h>
#include <custom.h>
#include <tmsnets.h>
#include <GFSRWU.h>
#include <permut.h>
#include <temperMK.h>

/* Read the main search parameter file testfile, initialize the Combinaison C with NbComp
   components, and fill pmecf, ptuplets, nbessais, MKopt, and Lmax from the file contents. */
void ReadSearch(Combinaison *C, int NbComp, char *testfile, paramMecf *pmecf,
                paramTuplets *ptuplets, int *nbessais, int *MKopt,
                int *Lmax) {
  boolean tupletsverif;
  int mse, d;
  float mDD;
  int *delta, *s;
  int i, j, jmin, jmax, nb_lignes_ecarts, ec;

  char method[15];
  char transfile[50];
  int deftrans;
  long seed1, seed2;   /* Seeds for the random number generator used during the parameter search */
  double Seed1, Seed2, Seed3, Seed4, Seed5, Seed6;
  boolean Readnbessais = FALSE, MK;
  boolean meverif;
  FILE *f;

  *MKopt = FALSE;
  *nbessais = 0;
  /* Open the main data file */
  f = fopen(testfile, "r");
  if (f == NULL) {
    printf("File \"%s\" not found\n", testfile);
    exit(1);
  }

  /* Read the seed for the random number generator used for tempering parameter selection */
  fscanf(f, "%ld", &seed1);
  if (seed1 == -1) {
    Seed1 = (double) time(NULL);
    Seed2 = Seed1 + 111119903.0;
  } else {
    fscanf(f, "%ld", &seed2);
    Seed1 = (double) seed1;
    Seed2 = (double) seed2;
  }
  Seed3 = Seed2 + 42221220.0;
  Seed4 = Seed1 + 733223331.0;
  Seed5 = 232555333446.0;
  Seed6 = Seed1 + 48392903.0;

  SetMRG32k3a(Seed1, Seed2, Seed3, Seed4, Seed5, Seed6);
  for (j = 0; j < 100; j++)
    MRG32k3a(); /* discard initial values to break initialization structure */
  ReadLn(f);
  fscanf(f, "%d", Lmax);
  ReadLn(f);
  delta = (int *) malloc((*Lmax + 1) * sizeof(int));
  /* Allocate NbComp Components in the Combinaison */
  AllocComponentsInCombinaison(C, NbComp, *Lmax);
  printf("====================================================================\n");
  printf("SUMMARY OF THE SEARCH PARAMETERS\n\n");
  printf("Computer : ");
  fflush(stdout);
  system("hostname");
  printf("\n");
  printf("Seed of RNG for tempering parameters = ( %12.0f, %12.0f )\n\n", Seed1, Seed2);
  if (NbComp == 1)
    printf("1 component:\n");
  else
    printf("%d components:\n", NbComp);
  for (j = 0; j < NbComp; j++) {
    fscanf(f, "%d", &deftrans);
    /* if deftrans==1, apply tempering to this component */
    if (deftrans == TRUE) {
      fscanf(f, "%s", transfile);
      Readnbessais = TRUE;
    }
    ReadLn(f);
    if (deftrans == TRUE) {
      ReadTempering(&(C->Components[j]), transfile, &MK);
      if (MK == TRUE)
        *MKopt = TRUE;
    }
  }
  fscanf(f, "%d", nbessais);
  ReadLn(f);
  if (Readnbessais) {
    printf("Number of tries per combined generator : %d\n", *nbessais);
  } else {
    *nbessais = 1;
  }

  /* Maximum allowed sum of gaps for the equidistribution criterion */
  fscanf(f, "%d", &mse);

  meverif = TRUE;
  if (mse == -1) {
    meverif = FALSE;
    mse = INT_MAX;
    SetMethod(pmecf, METHOD_NOTHING);
  }
  if (meverif) {
    fscanf(f, "%s", method);
    if (strcmp("matrix", method) == 0)
      SetMethod(pmecf, METHOD_MATRICIAL);
    else if (strcmp("lattice", method) == 0)
      SetMethod(pmecf, METHOD_DUALLATTICE);
    else {
      printf("In file %s, Method %s is unrecognized\n.  Search stopped.", testfile, method);
      exit(1);
    }
    printf("Upperbound for the sum of dimension gaps for resolutions in psi_12 : %d\n", mse);
  }
  ReadLn(f);
  /* Number of lines specifying the maximum allowed gaps per resolution range */
  fscanf(f, "%d", &nb_lignes_ecarts);
  ReadLn(f);
  if (nb_lignes_ecarts > 0) {
    printf("delta for particular bits:\n");
    for (i = 0; i < nb_lignes_ecarts; i++) {
      fscanf(f, "%d %d %d", &jmin, &jmax, &ec);
      printf("   >>Bits from %d to %d delta_i=%d\n", jmin, jmax, ec);
      ReadLn(f);
      for (j = jmin; j <= jmax; j++)
        delta[j] = ec;
    }
  } else {
    for (i = 0; i <= *Lmax; i++) {
      delta[i] = 10000000;
    }
  }
  AllocMecf(pmecf, *Lmax);
  InitParamMecf(pmecf, delta, mse, 0, meverif, FALSE);

  mDD = 0.0;
  /* tupletsverif == 1 means compute the Delta(t_1,..,t_d) criterion */
  fscanf(f, "%d", &tupletsverif);
  if (tupletsverif == 1) {
    fscanf(f, "%d", &d); /* value of d in Delta(t_1,..,t_d) */
    s = (int *) malloc((d + 1) * sizeof(int));
    for (j = 1; j <= d; j++)
      fscanf(f, "%d", &s[j]); /* values t[j] in Delta(t_1,..,t_d) */
    fscanf(f, "%f", &mDD);    /* maximum allowed value of DD */
    printf("Verification of DELTA( ");
    for (j = 1; j < d; j++)
      printf("%d, ", s[j]);
    printf("%d)\n", s[d]);
    AllocTuplets(ptuplets, d);
    InitParamTuplets(ptuplets, d, s, mDD, tupletsverif, EQUIDISTRIBUTION, MAX);
  } else {
    InitParamTuplets(ptuplets, 0, NULL, mDD, FALSE, EQUIDISTRIBUTION, MAX);
  }
  ReadLn(f);
  printf("====================================================================\n");
  fflush(stdout);
  fclose(f);
}

/* Fill the Combinaison C with NbComp components by reading each generator data file
   listed in gendatafile and dispatching to the appropriate ReadData function. */
void ReadGenDataFiles(Combinaison *C, int NbComp, char **gendatafile, int Lmax) {
  int j;
  FILE *f2;
  boolean same;
  char gen_type[25], old_gen[50];

  if (C->J != NbComp)
    AllocComponentsInCombinaison(C, NbComp, Lmax);

  for (j = 0; j < NbComp; j++) {
    /* If gendatafile[j] is "same", the component uses the same generators as component j-1. */
    if (strncmp(gendatafile[j], "same", 4) == 0) {
      same = TRUE;
      strcpy(gendatafile[j], old_gen);
    } else {
      same = FALSE;
      strcpy(old_gen, gendatafile[j]);
    }
    /* Open the file to read the generator type tag */
    f2 = fopen(gendatafile[j], "r");
    if (f2 == NULL) {
      printf("File \"%s\" not found\n", gendatafile[j]);
      exit(1);
    }
    fscanf(f2, "%s", gen_type);
    fclose(f2);

    /* Read data for the generators to be tested */
    if (strcmp("taus", gen_type) == 0) {
      ReadDataTaus(&(C->Components[j]), gendatafile[j], same, Lmax);
      printf("- Component %d: Tausworthe \n", j + 1);
    } else if (strcmp("polylcg", gen_type) == 0) {
      ReadDataPolyLCG(&(C->Components[j]), gendatafile[j], same, Lmax);
      printf("- Component %d: Polynomial LCG \n", j + 1);
    } else if (strcmp("tms", gen_type) == 0) {
      ReadDataTMS(&(C->Components[j]), gendatafile[j], same, Lmax);
      printf("- Component %d: (t,m,s)-net\n", j + 1);
    } else if (strcmp("carry", gen_type) == 0) {
      ReadDataCarry(&(C->Components[j]), gendatafile[j], same, Lmax);
      printf("- Component %d: carry\n", j + 1);
    } else if (strcmp("GFSR", gen_type) == 0) {
      ReadDataGFSR(&(C->Components[j]), gendatafile[j], same, Lmax);
      printf("- Component %d: GFSR\n", j + 1);
    } else if (strcmp("marsa", gen_type) == 0) {
      ReadDataMarsa(&(C->Components[j]), gendatafile[j], same, Lmax);
      printf("- Component %d: Marsaglia\n", j + 1);
    } else if (strcmp("marsaxorshift", gen_type) == 0) {
      ReadDataMarsaXorshift(&(C->Components[j]), gendatafile[j], same, Lmax);
      printf("- Component %d: Marsaglia Xor Shift RNG\n", j + 1);
    } else if (strcmp("tgfsr", gen_type) == 0) {
      ReadDataTGFSR(&(C->Components[j]), gendatafile[j], same, Lmax);
      printf("- Component %d: TGFSR \n", j + 1);
    } else if (strcmp("genf2w", gen_type) == 0) {
      ReadDataGenF2w(&(C->Components[j]), gendatafile[j], same, Lmax);
      printf("- Component %d: Generator in F_{2^w} \n", j + 1);
    } else if (strcmp("cust", gen_type) == 0) {
      ReadDataCust(&(C->Components[j]), gendatafile[j], same, Lmax);
      printf("- Component %d: Custom \n", j + 1);
    } else if (strcmp("ac1d", gen_type) == 0) {
      ReadDataAC1D(&(C->Components[j]), gendatafile[j], same, Lmax);
      printf("- Component %d: AC-1D \n", j + 1);
    } else if (strcmp("MT", gen_type) == 0) {
      ReadDataMT(&(C->Components[j]), gendatafile[j], same, Lmax);
      printf("- Component %d: Mersenne Twister \n", j + 1);
    } else if (strcmp("matsumoto", gen_type) == 0) {
      ReadDataMatsumoto(&(C->Components[j]), gendatafile[j], same, Lmax);
      printf("- Component %d: New Generator with M. Matsumoto \n", j + 1);
    } else {
      printf("Unsupported generator= (%s)\n", gen_type);
      exit(1);
    }
  }
}

/* Read the tempering transformation file filename and fill the Component E with the
   transformations described; update MKopt if any optimized MK tempering is found. */
void ReadTempering(Component *E, char *filename, int *MKopt) {
  int i, k, param1, param2, precision;
  FILE *f;
  int NbTrans, aleatoire, M;
  char lintrans[25];
  BitVect init_b, init_c;
  boolean Temp_opt, D;
  int limitv = INT_MAX;      /* Maximum resolution up to which MK tempering is optimized */
  boolean affprogress = FALSE; /* Whether to display optimization progress */
  Transformation *T;

  *MKopt = FALSE;
  /* Open the tempering transformation file */
  f = fopen(filename, "r");
  if (f == NULL) {
    printf("File \"%s\" not found\n", filename);
    exit(1);
  }
  fscanf(f, "%d", &NbTrans);
  /* Allocate NbTrans linear transformations for Component E */
  AllocTransInComponent(E, NbTrans);
  ReadLn(f);
  if (NbTrans > 0) {
    printf("  Tempering transformations:\n");
    for (i = 0; i < NbTrans; i++) {
      /* Read the i-th linear transformation descriptor */
      fscanf(f, "%s %d %d", lintrans, &precision, &param1);
      if (strcmp("permut", lintrans) == 0) {
        if ((T = AllocPermut()) == NULL) {
          printf("Error in ReadTempering()\n");
          exit(1);
        }
        printf("   * Permutation\n");
        fscanf(f, "%d", &param2);
        ReadLn(f);
        InitPermut(T, precision, param1, param2);
        /* Add permutation to Component E */
        AddTransInComponent(E, T);
        FreePermut(T);
      } else if (strcmp("tempMK",    lintrans) == 0 ||
                 strcmp("tempMKopt", lintrans) == 0 ||
                 strcmp("tempMK2",   lintrans) == 0 ||
                 strcmp("tempMK2opt",lintrans) == 0) {
        int param_l, param_u, type, param_tt;

        AllocBV(&init_b, precision);
        AllocBV(&init_c, precision);
        if ((T = AllocTemperMK()) == NULL) {
          printf("Error in ReadTempering()\n");
          exit(1);
        }
        fscanf(f, "%d", &param_tt);
        if (param_tt == -1 || param1 == -1) {
          printf("s and t must be defined for the MK tempering\n");
          exit(1);
        }
        if (strcmp("tempMK2", lintrans) == 0 || strcmp("tempMK2opt", lintrans) == 0) {
          type = 2;
          fscanf(f, "%d %d", &param_u, &param_l);
          printf("   * Matsumoto-Kurita Tempering(II)\n");
        } else {
          printf("   * Matsumoto-Kurita(I) Tempering\n");
          type = 1;
          param_u = param_l = 0;
        }
        fscanf(f, "%d", &aleatoire);

        Temp_opt = FALSE;
        if (strcmp("tempMKopt", lintrans) == 0 || strcmp("tempMK2opt", lintrans) == 0) {
          *MKopt = TRUE;
          Temp_opt = TRUE;
          fscanf(f, "%d %d", &D, &M);
          affprogress |= D;
          if (M < limitv) {
            limitv = M;
          }
        }
        ReadLn(f);
        if (aleatoire != -1) {
          for (k = 0; k < aleatoire; k++)
            fscanf(f, "%x", &init_b.vect[k]);
          ReadLn(f);
          for (k = 0; k < aleatoire; k++)
            fscanf(f, "%x", &init_c.vect[k]);
          ReadLn(f);
        } else {
          PutBVToZero(&init_b);
          PutBVToZero(&init_c);
        }

        InitTemperMK(T, precision, param1, param_tt, param_u, param_l, aleatoire,
                     &init_b, &init_c, Temp_opt, type, affprogress, limitv);
        /* Add MK tempering to Component E */
        AddTransInComponent(E, T);
        FreeTemperMK(T);
        FreeBV(&init_b);
        FreeBV(&init_c);
      } else {
        printf("Unrecognized Tempering Transformation : %s ", lintrans);
        exit(1);
      }
    }
  }
  fclose(f);
}
