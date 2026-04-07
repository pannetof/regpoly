/*
 * test_transition_matrix_generator.c — Display the transition matrix of every
 * generator in a data file.
 *
 * Usage:
 *     ./test_transition_matrix_generator <generator_file>
 *
 * The file type is detected from its first token (e.g. "polylcg", "taus",
 * "taus2").  Each generator is loaded into a Component and its transition
 * matrix is computed and printed in turn.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "combinaisons.h"
#include "polynomials.h"
#include "polylcg.h"
#include "tausworthe.h"
#include "tgfsr.h"
#include "MT.h"
#include "genf2w.h"
#include "carry2.h"
#include "marsaxorshift.h"
#include "matsumoto.h"

int main(int argc, char *argv[])
{
    FILE      *f;
    char       gen_type[25];
    char       class_name[25];
    Combinaison C;
    Component  *comp;
    int        i, Lmax = WL;   /* WL = 32: prevents GenState reallocation mismatch in CopyTaus */

    if (argc != 2) {
        fprintf(stderr, "Usage: %s <generator_file>\n", argv[0]);
        return 1;
    }

    /* Read the type tag from the first token of the file */
    f = fopen(argv[1], "r");
    if (f == NULL) {
        fprintf(stderr, "File \"%s\" not found\n", argv[1]);
        return 1;
    }
    fscanf(f, "%s", gen_type);
    fclose(f);

    /* Map type tag to class name */
    if (strcmp(gen_type, "polylcg") == 0) {
        strcpy(class_name, "PolyLCG");
    } else if (strcmp(gen_type, "taus") == 0 || strcmp(gen_type, "taus2") == 0) {
        strcpy(class_name, "Tausworthe");
    } else if (strcmp(gen_type, "tgfsr") == 0) {
        strcpy(class_name, "TGFSR");
    } else if (strcmp(gen_type, "MT") == 0) {
        strcpy(class_name, "MersenneTwister");
    } else if (strcmp(gen_type, "genf2w") == 0) {
        strcpy(class_name, "GenF2w");
    } else if (strcmp(gen_type, "carry") == 0) {
        strcpy(class_name, "Carry");
    } else if (strcmp(gen_type, "marsaxorshift") == 0) {
        strcpy(class_name, "MarsaXorshift");
    } else if (strcmp(gen_type, "matsumoto") == 0) {
        strcpy(class_name, "Matsumoto");
    } else {
        fprintf(stderr, "Unknown generator type '%s'\n", gen_type);
        return 1;
    }

    /* Allocate a Combinaison with one Component */
    AllocComponentsInCombinaison(&C, 1, Lmax);
    comp = &C.Components[0];

    /* Load generators into the Component */
    if (strcmp(gen_type, "polylcg") == 0) {
        ReadDataPolyLCG(comp, argv[1], FALSE, Lmax);
    } else if (strcmp(gen_type, "tgfsr") == 0) {
        ReadDataTGFSR(comp, argv[1], FALSE, Lmax);
        for (i = 0; i < comp->NbGen; i++) {
            if (comp->Gen[i]->L < comp->Gen[i]->k) {
                comp->Gen[i]->L = comp->Gen[i]->k;
                FreeBV(&(comp->Gen[i]->GenState));
                AllocBV(&(comp->Gen[i]->GenState), comp->Gen[i]->k);
            }
        }
    } else if (strcmp(gen_type, "MT") == 0) {
        ReadDataMT(comp, argv[1], FALSE, Lmax);
        for (i = 0; i < comp->NbGen; i++) {
            if (comp->Gen[i]->L < comp->Gen[i]->k) {
                comp->Gen[i]->L = comp->Gen[i]->k;
                FreeBV(&(comp->Gen[i]->GenState));
                AllocBV(&(comp->Gen[i]->GenState), comp->Gen[i]->k);
            }
        }
    } else if (strcmp(gen_type, "genf2w") == 0) {
        ReadDataGenF2w(comp, argv[1], FALSE, Lmax);
        for (i = 0; i < comp->NbGen; i++) {
            if (comp->Gen[i]->L < comp->Gen[i]->k) {
                comp->Gen[i]->L = comp->Gen[i]->k;
                FreeBV(&(comp->Gen[i]->GenState));
                AllocBV(&(comp->Gen[i]->GenState), comp->Gen[i]->k);
            }
        }
    } else if (strcmp(gen_type, "carry") == 0) {
        ReadDataCarry(comp, argv[1], FALSE, Lmax);
        for (i = 0; i < comp->NbGen; i++) {
            if (comp->Gen[i]->L < comp->Gen[i]->k) {
                comp->Gen[i]->L = comp->Gen[i]->k;
                FreeBV(&(comp->Gen[i]->GenState));
                AllocBV(&(comp->Gen[i]->GenState), comp->Gen[i]->k);
            }
        }
    } else if (strcmp(gen_type, "marsaxorshift") == 0) {
        ReadDataMarsaXorshift(comp, argv[1], FALSE, Lmax);
        for (i = 0; i < comp->NbGen; i++) {
            if (comp->Gen[i]->L < comp->Gen[i]->k) {
                comp->Gen[i]->L = comp->Gen[i]->k;
                FreeBV(&(comp->Gen[i]->GenState));
                AllocBV(&(comp->Gen[i]->GenState), comp->Gen[i]->k);
            }
        }
    } else if (strcmp(gen_type, "matsumoto") == 0) {
        ReadDataMatsumoto(comp, argv[1], FALSE, Lmax);
        for (i = 0; i < comp->NbGen; i++) {
            if (comp->Gen[i]->L < comp->Gen[i]->k) {
                comp->Gen[i]->L = comp->Gen[i]->k;
                FreeBV(&(comp->Gen[i]->GenState));
                AllocBV(&(comp->Gen[i]->GenState), comp->Gen[i]->k);
            }
        }
    } else {
        ReadDataTaus(comp, argv[1], FALSE, Lmax);
    }

    /* Print header */
    printf("File   : %s\n", argv[1]);
    printf("Type   : %s  (%s)\n", gen_type, class_name);
    printf("Generators: %d\n", comp->NbGen);
    printf("\n");

    /* Loop through all generators */
    for (i = 0; i < comp->NbGen; i++) {
        comp->Gen[i]->DispGen(comp->Gen[i]);
        TransitionMatrix(comp->Gen[i]);
        printf("\n");
    }

    return 0;
}
