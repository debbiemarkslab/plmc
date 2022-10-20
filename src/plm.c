/*
 *  plmc
 *  Copyright (c) 2016, John Ingraham
 *  john.ingraham@gmail.com
 */

#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include "include/twister.h"
#include "include/plm.h"
#include "include/inference.h"
#include "include/weights.h"

/* Internal functions to MSARead */
void MSAReadSeq(char *seq, FILE *fpAli);
letter_t MSAReadCode(char c, char *alphabet, int nCodes);

/* Global verbosity & profiling options */
int verbose = 2;

/* Reference amino acid indexing */
const char *codesAA = "-ACDEFGHIKLMNPQRSTVWY";

/* Default parameters */
const numeric_t REGULARIZATION_LAMBDA_H = 0.01;
const numeric_t REGULARIZATION_LAMBDA_E = 100.0;
const numeric_t REGULARIZATION_LAMBDA_GROUP = 0.0;
const numeric_t REWEIGHTING_THETA = 0.20;
const numeric_t REWEIGHTING_SCALE = 1.0;
const int ZERO_APC_PRIORS = 0;
const int SGD_BATCH_SIZE = 2048;
const int REWEIGHTING_SAMPLES = 5000;

options_t* default_options() {
    /* Default options */
    options_t *options = (options_t *) malloc(sizeof(options_t));

    options->lambdaH = REGULARIZATION_LAMBDA_H;
    options->lambdaE = REGULARIZATION_LAMBDA_E;
    options->lambdaGroup = REGULARIZATION_LAMBDA_GROUP;
    options->scale = REWEIGHTING_SCALE;
    options->zeroAPC = 0;
    options->maxIter = 0;
    options->usePairs = 1;
    options->estimator = INFER_MAP;
    options->estimatorMAP = INFER_MAP_PLM;
    options->target = NULL;
    options->alphabet = (char *) codesAA;

    /* Sequence weight options */
    options->fastWeights = 0;
    options->theta = REWEIGHTING_THETA;

    /* SGD options */
    options->sgd = 0;
    options->sgdBatchSize = SGD_BATCH_SIZE;

    return options;
}

void run_plmc(char *alignFile, char* outputFile, char *couplingsFile,
    char *weightsFile, char *weightsOutputFile, options_t *options) {

    /* Initialize PRNG */
    init_genrand(42);

    /* Read multiple sequence alignment */
    alignment_t *ali = MSARead(alignFile, options);

    if (weightsFile != NULL) {
        ReadCustomWeightsFile(weightsFile, ali);
    } else {
        /* Reweight sequences by inverse neighborhood density */
        MSAReweightSequences(ali, options);
    }
    if (weightsOutputFile != NULL) {
        WriteWeightsFile(weightsOutputFile, ali);
    }

    /* Compute sitewise and pairwise marginal distributions */
    MSACountMarginals(ali, options);

    /* Infer model parameters */
    numeric_t *x = InferPairModel(ali, options);

    /* (Optionally) Output estimated parameters and coupling scores */
    if (outputFile != NULL)
        OutputParametersFull(outputFile, x, ali, options);
    if (couplingsFile != NULL)
        OutputCouplingScores(couplingsFile, x, ali, options);

    /* Free alignment and options */
    MSAFree(ali, options);
}

alignment_t *MSARead(char *alignFile, options_t *options) {
    /* Read FASTA-formatted alignment */
    FILE *fpAli = NULL;
    if (alignFile != NULL) {
        fpAli = fopen(alignFile, "r");
    } else {
        fprintf(stderr, "Must specify alignment file: -a ALIGN_FILE\n");
        exit(1);
    }
    if (fpAli == NULL) {
        fprintf(stderr, "Error opening alignment file\n");
        exit(1);
    }

    /* Allocate alignment */
    alignment_t *ali = (alignment_t *) malloc(sizeof(alignment_t));
    ali->nSeqs = ali->nSites = ali->nCodes = 0;
    ali->alphabet = options->alphabet;
    ali->names = NULL;
    ali->sequences = NULL;
    ali->target = -1;
    ali->offsets = NULL;
    ali->nEff = 0;
    ali->weights = ali->fi = ali->fij = NULL;
    ali->nParams = 0;

    /* Verify alignment dimensions and structure (first pass through file) */
    char name[BUFFER_SIZE];
    char seq[BUFFER_SIZE];
    /* Read first line as name */
    fgetstr(name, fpAli);
    if (*name == '>') {
        MSAReadSeq(seq, fpAli);
    } else {
        fprintf(stderr, "Error reading alignment:"
                        " First line should start with >\n");
        exit(1);
    }
    ali->nCodes = strlen(ali->alphabet);
    ali->nSites = strlen(seq);
    ali->nSeqs = 1;
    while (!feof(fpAli)) {
        char c = fgetc(fpAli);
        if (c == '>') {
            /* Read name and sequence pair */
            fgetstr(name, fpAli);
            MSAReadSeq(seq, fpAli);
        } else {
            fprintf(stderr, "Error reading alignment:"
                            " sequence records should start with >\n");
            exit(1);
        }

        /* Validate sequence length */
        if (strlen(seq) != ali->nSites) {
            fprintf(stderr,
                "Incompatible sequence length (%lu should be %d) for %s:\n%s\n",
                strlen(seq), ali->nSites, name, seq);
            exit(1);
        }
        ali->nSeqs++;
    }

    /* Encode full alignment block (second pass through file) */
    ali->sequences = (letter_t *)
        malloc(ali->nSites * ali->nSeqs * sizeof(letter_t));
    ali->names = (char **) malloc(ali->nSeqs * sizeof(char *));
    for (int s = 0; s < ali->nSeqs; s++)
        for (int i = 0; i < ali->nSites; i++) seq(s, i) = 0;
    for (int s = 0; s < ali->nSeqs; s++) ali->names[s] = NULL;
    rewind(fpAli);
    for (int s = 0; s < ali->nSeqs; s++) {
        /* >Name */
        getc(fpAli);
        fgetstr(name, fpAli);
        ali->names[s] = (char *) malloc((strlen(name) + 1) * sizeof(char));
        strcpy(ali->names[s], name);

        /* Sequence */
        MSAReadSeq(seq, fpAli);
        for (int i = 0; i < ali->nSites; i++)
            seq(s, i) = MSAReadCode(seq[i], ali->alphabet, ali->nCodes);
    }

    /* --------------------------------_DEBUG_--------------------------------*/
    /* Alignment to stderr */
    // for (int s = 0; s < 10; s++) {
    // for (int s = 0; s < ali->nSeqs; s++) {
    //     for (int i = 0; i < ali->nSites; i++)
    //         if (seq(s, i) >= 0 && seq(s, i) < ali->nCodes) {
    //             fprintf(stderr, "%c", ali->alphabet[seq(s, i)]);
    //         } else if (seq(s, i) >= -ali->nCodes && seq(s, i) < 0) {
    //             fprintf(stderr, "%c",
    //                 tolower(ali->alphabet[seq(s, i) + ali->nCodes]));
    //         } else {
    //             fprintf(stderr, "*%d*", seq(s, i));
    //         }
    //     fprintf(stderr, "\n");
    // }
    // exit(0);
    /* --------------------------------^DEBUG^--------------------------------*/

    /* Focus mode: If a focus sequence (target) is provided, locate it */
    if (options->target != NULL) {
        for (int s = 0; s < ali->nSeqs; s++)
            if (strncmp(options->target, ali->names[s],
                strlen(options->target)) == 0) {
                if (ali->target >= 0) {
                    fprintf(stderr,
                        "Multiple sequences start with %s, picking sequence %d\n",
                        options->target, s + 1);
                } else {
                    ali->target = s;
                }
            }
        if (ali->target >= 0) {
            fprintf(stderr, "Found focus %s as sequence %d\n", options->target,
                ali->target + 1);
        } else {
            fprintf(stderr,
                "Could not find %s, proceeding without focus sequence\n",
                options->target);
        }
    }

    /* Always discard any sequences (rows) with out-of-alphabet characters */
    int* seqValid = (int *) malloc(ali->nSeqs * sizeof(int));
    for (int s = 0; s < ali->nSeqs; s++) seqValid[s] = 0;
    for (int s = 0; s < ali->nSeqs; s++)
        for (int i = 0; i < ali->nSites; i++)
            if ((seq(s, i) >= -ali->nCodes) && (seq(s, i) < ali->nCodes))
                seqValid[s]++;
    int nValidSeqs = 0;
    for (int s = 0; s < ali->nSeqs; s++)
        if (seqValid[s] == ali->nSites) nValidSeqs++;
    fprintf(stderr, "%d valid sequences out of %d \n", nValidSeqs, ali->nSeqs);
    
    /* Record indices of skipped sequences */
    ali->nSkippedSeqs = ali->nSeqs - nValidSeqs;
    ali->skippedSeqs = (int *) malloc(ali->nSkippedSeqs * sizeof(int));
    for (int s = 0, skipIndex = 0; s < ali->nSeqs; s++)
        if (seqValid[s] != ali->nSites) ali->skippedSeqs[skipIndex++] = s;

    /* Focus mode: select only focus columns (criteria below) */
    int nValidSites = ali->nSites;
    int* siteValid = (int *) malloc(ali->nSites * sizeof(int));
    for (int i = 0; i < ali->nSites; i++) siteValid[i] = 1;
    if (ali->target >= 0) {
        for (int i = 0; i < ali->nSites; i++) {
            /* For proteins, remove lower case and gap columns */
            if ((ali->alphabet == codesAA) 
                && (seq(ali->target, i) < 0))
                siteValid[i] = 0;
            /* Discard gaps */
            if ((ali->alphabet == codesAA)
                || (options->estimatorMAP == INFER_MAP_PLM_GAPREDUCE))
                if (seq(ali->target, i) == 0) siteValid[i] = 0;
        }
        nValidSites = 0;
        for (int i = 0; i < ali->nSites; i++)
            if (siteValid[i] == 1) nValidSites++;
        fprintf(stderr,
            "%d sites out of %d\n", nValidSites, ali->nSites);
    } else {
        fprintf(stderr,
            "%d sites\n", ali->nSites);
    }

    /* Focus mode: parse region (NAME/START_IX-END_IX) and map offsets */
    int leftOffset = 0;
    if (ali->target >= 0) {
        char *focusName = ali->names[ali->target];
        /* Name should be immediately followed by '/' */
        if (strlen(focusName) > strlen(options->target) + 1
            && focusName[strlen(options->target)] == '/') {
            /* Attempt to read integer region start */
            int regLeft = strlen(options->target) + 1;
            int ix = 0;
            if (isdigit(focusName[regLeft])) {
                while (regLeft + ix < strlen(focusName)
                       && isdigit(focusName[regLeft + ix + 1])) ix++;
                int tens = 1;
                leftOffset = -1;
                for (int i = ix; i >= 0; i--) {
                    leftOffset += tens * (focusName[regLeft + i] - '0');
                    tens *= 10;
                }
                fprintf(stderr, "Region starts at %d\n", leftOffset + 1);
            } else {
                fprintf(stderr, "Error parsing region, assuming start at 1");
            }
        }

        /* Map the offsets */
        ali->offsets = (int *) malloc(nValidSites * sizeof(int));
        for (int i = 0; i < nValidSites; i++) ali->offsets[i] = i + 1;
        int ix = 0;
        for (int i = 0; i < ali->nSites; i++)
            if (siteValid[i] == 1) {
                ali->offsets[ix] = i + 1 + leftOffset;
                ix++;
            }

        /* Reposition the target for reduced alignment */
        int targetShift = -1;
        for (int i = 0; i <= ali->target; i++)
            if (seqValid[i] == ali->nSites) targetShift++;
        ali->target = targetShift;
    }

    /* Copy only selected rows and columns */
    if (nValidSeqs < ali->nSeqs || nValidSites < ali->nSites) {
        letter_t *seqsReduced = (letter_t *)
            malloc(nValidSites * nValidSeqs * sizeof(letter_t));
        for (int i = 0; i < nValidSites * nValidSeqs; i++) seqsReduced[i] = 0;
        int sx = 0;
        for (int s = 0; s < ali->nSeqs; s++)
            if (seqValid[s] == ali->nSites) {
                int ix = 0;
                for (int i = 0; i < ali->nSites; i++) {
                    if (siteValid[i] == 1) {
                        seqsReduced[ix + sx * nValidSites] = seq(s, i);
                        ix++;
                    }
                }
                sx++;
            }

        /* Reallocate alignment with reduced dimensions */
        free(ali->sequences);
        ali->nSeqs = nValidSeqs;
        ali->nSites = nValidSites;
        ali->sequences = (letter_t *)
            malloc(nValidSites * nValidSeqs * sizeof(letter_t));
        for (int i = 0; i < nValidSites * nValidSeqs; i++)
            ali->sequences[i] = 0;
        for (int s = 0; s < nValidSeqs; s++)
            for (int i = 0; i < nValidSites; i++)
                seq(s, i) = seqsReduced[i + s * nValidSites];
        free(seqsReduced);
    }

    /* Shift any lowercase codes back to uppercase */
    for (int s = 0; s < ali->nSeqs; s++)
        for (int i = 0; i < ali->nSites; i++)
            if (seq(s, i) < 0) seq(s, i) += ali->nCodes;

    /* Intialize weights to 1.0 */
    ali->weights = (numeric_t *) malloc(ali->nSeqs * sizeof(numeric_t));
    for (int s = 0; s < ali->nSeqs; s++) ali->weights[s] = 1.0;
    ali->nEff = (numeric_t) ali->nSeqs;

    /* --------------------------------_DEBUG_--------------------------------*/
    /* Display offset map */
    // for (int i = 0; i < ali->nSites; i++) {
    //     fprintf(stderr, "%d : %d : %c\n", i + 1, ali->offsets[i],
    //             ali->alphabet[seq(ali->target, i)]);
    // }
    // exit(0);
    /* Display target */
    // for (int i = 0; i < ali->nSites; i++) {
    //     fprintf(stderr, "%c", ali->alphabet[seq(ali->target, i)]);
    // }
    // fprintf(stderr, "\n");
    // exit(0);
    /* --------------------------------^DEBUG^--------------------------------*/

    /* --------------------------------_DEBUG_--------------------------------*/
    // for (int s = 0; s < ali->nSeqs; s++) {
    //     fprintf(stderr, ">%s\n", ali->names[s]);
    //     for (int i = 0; i < ali->nSites; i++)
    //         fprintf(stderr, "%c", ali->alphabet[seq(s, i)]);
    //     fprintf(stderr, "\n");
    // }
    /* --------------------------------^DEBUG^--------------------------------*/
    return ali;
}

void MSAReadSeq(char *seq, FILE *fpAli) {
    /* Read sequence from the current line(s) */
    char buf[BUFFER_SIZE];
    /* Look ahead one character */
    char c = fgetc(fpAli); 
    ungetc(c, fpAli);
    seq[0] = '\0';
    while (c != '>' && !feof(fpAli)) {
        fgetstr(buf, fpAli);
        strcat(seq, buf);
        /* Look ahead one character */
        c = fgetc(fpAli);
        ungetc(c, fpAli);
    }
}

letter_t MSAReadCode(char c, char *alphabet, int nCodes) {
    /* Encode a character as an integer between -nCodes and +nCodes
          In alphabet:                     store index           [0, nCodes - 1]
          Lowercase version of alphabet:   downshift by nCodes   [-nCodes, -1]
          Out of alphabet:                 store nCodes          [nCodes]
     */
    letter_t i = 0;

    /* Protein-specific treatment of '.' */
    if (alphabet == codesAA) if (c == '.') c = '-';

    /* Store lowercase characters as down-shifted by nCodes */
    while ((i < nCodes - 1) && toupper(c) != alphabet[i]) i++;
    if (c != alphabet[i] && toupper(c) == alphabet[i]) i -= nCodes;

    /* Encode out-of-alphabet characters by [nCodes] */
    if (i > 0 && toupper(c) != alphabet[i]) i = nCodes;
    return i;
}

void MSACountMarginals(alignment_t *ali, options_t *options) {
    /* Compute first and second order marginal distributions, according to the 
       sequence weights
     */
    if (options->estimatorMAP == INFER_MAP_PLM_GAPREDUCE) {
        /* Condition the marginals on ungapped */
        ali->nCodes = strlen(ali->alphabet) - 1;

        /* First-order marginals P_i(Ai) */
        int nFi = ali->nSites * ali->nCodes;
        ali->fi = (numeric_t *) malloc(nFi * sizeof(numeric_t));
        for (int i = 0; i < nFi; i++) ali->fi[i] = 0.0;

        for (int s = 0; s < ali->nSeqs; s++)
            for (int i = 0; i < ali->nSites; i++)
                if (seq(s, i) > 0)
                    fi(i, seq(s, i) - 1) += ali->weights[s];

        /* Second-order marginals P_ij(Ai, Aj) */
        int nFij = ali->nSites * (ali->nSites - 1) / 2 * ali->nCodes * ali->nCodes;
        ali->fij = (numeric_t *) malloc(nFij * sizeof(numeric_t));
        for (int i = 0; i < nFij; i++) ali->fij[i] = 0.0;

        for (int s = 0; s < ali->nSeqs; s++)
            for (int i = 0; i < ali->nSites - 1; i++)
                for (int j = i + 1; j < ali->nSites; j++)
                    if (seq(s, i) > 0) if(seq(s, j) > 0)
                        fij(i, j, seq(s, i) - 1, seq(s, j) - 1)
                            += ali->weights[s];

        /* Normalize conditional distributions */
        for (int i = 0; i < ali->nSites; i++) {
            double fsum = 0.0;
            for (int ai = 0; ai < ali->nCodes; ai++)
                fsum += fi(i, ai);
            if (fsum != 0) {
                double fsumInv = 1.0 / fsum;
                for (int ai = 0; ai < ali->nCodes; ai++)
                    fi(i, ai) *= fsumInv;
            } else {
                /* Handle empty columns */
                numeric_t flatF = 1.0 / ((numeric_t) ali->nCodes);
                for (int ai = 0; ai < ali->nCodes; ai++)
                    fi(i, ai) = flatF;
            }
        }
        for (int i = 0; i < ali->nSites - 1; i++)
            for (int j = i + 1; j < ali->nSites; j++) {
                double fsum = 0.0;
                for (int ai = 0; ai < ali->nCodes; ai++)
                    for (int aj = 0; aj < ali->nCodes; aj++)
                        fsum += fij(i, j, ai, aj);
                if (fsum != 0) {
                    double fsumInv = 1.0 / fsum;
                    for (int ai = 0; ai < ali->nCodes; ai++)
                        for (int aj = 0; aj < ali->nCodes; aj++)
                            fij(i, j, ai, aj) *= fsumInv;
                } else {
                    /* Handle pairs of empty columns */
                    numeric_t flatF = 1.0 / ((numeric_t) (ali->nCodes * ali->nCodes));
                    for (int ai = 0; ai < ali->nCodes; ai++)
                        for (int aj = 0; aj < ali->nCodes; aj++)
                            fij(i, j, ai, aj) = flatF;
                }
            }

    } else {
        /* Compute regular marginals */
        numeric_t Zinv = 1.0 / ali->nEff;

        /* First-order marginals P_i(Ai) */
        int nFi = ali->nSites * ali->nCodes;
        ali->fi = (numeric_t *) malloc(nFi * sizeof(numeric_t));
        for (int i = 0; i < nFi; i++) ali->fi[i] = 0.0;

        for (int s = 0; s < ali->nSeqs; s++)
            for (int i = 0; i < ali->nSites; i++)
                fi(i, seq(s, i)) += ali->weights[s] * Zinv;

        /* Second-order marginals P_ij(Ai, Aj) */
        int nFij = ali->nSites * (ali->nSites - 1) / 2 * ali->nCodes * ali->nCodes;
        ali->fij = (numeric_t *) malloc(nFij * sizeof(numeric_t));
        for (int i = 0; i < nFij; i++) ali->fij[i] = 0.0;

        for (int s = 0; s < ali->nSeqs; s++)
            for (int i = 0; i < ali->nSites - 1; i++)
                for (int j = i + 1; j < ali->nSites; j++)
                    fij(i, j, seq(s, i), seq(s, j)) += ali->weights[s] * Zinv;
    }
}

void MSAFree(alignment_t *ali, options_t *options) {
    /* Free alignment and options */
    if (ali->names && ali->names[0])
        for (int i = 0; i < ali->nSeqs; i++) free(ali->names[i]);
    free(ali->names);
    free(ali->sequences);
    free(ali->weights);
    free(ali->fi);
    free(ali->fij);

    /* Note: options->target and options->alphabet are never allocated */
    free(options);
}

#define OUTPUT_PRECISION float
void OutputParametersSite(char *outputFile, const numeric_t *x,
    alignment_t *ali) {
    FILE *fpOutput = NULL;
    fpOutput = fopen(outputFile, "w");
    if (fpOutput != NULL) {
        /* 1: nSites */
        fwrite(&(ali->nSites), sizeof(ali->nSites), 1, fpOutput);

        /* 2: (Focus mode only) target sequence */
        if (ali->target >= 0) {
            for (int i = 0; i < ali->nSites; i++) {
                char c = (char) ali->alphabet[seq(ali->target, i)];
                fwrite(&c, sizeof(char), 1, fpOutput);
            }
        } else {
            char c = ali->alphabet[0];
            for (int i = 0; i < ali->nSites; i++)
                fwrite(&c, sizeof(c), 1, fpOutput);
        }

        /* 3: (Focus mode only) offset map */
        if (ali->target >= 0) {
            for (int i = 0; i < ali->nSites; i++) {
                int ix = ali->offsets[i];
                fwrite(&ix, sizeof(ix), 1, fpOutput);
            }
        } else {
            for (int i = 0; i < ali->nSites; i++) {
                int ix = i + 1;
                fwrite(&ix, sizeof(ix), 1, fpOutput);
            }
        }

        /* 4,5: sitewise marginals fi, twice */
        for (int x = 0; x < 2; x++)
            for (int i = 0; i < ali->nSites; i++)
                for (int ai = 0; ai < ali->nCodes; ai++) {
                    OUTPUT_PRECISION f = (OUTPUT_PRECISION) fi(i, ai);
                    fwrite(&f, sizeof(f), 1, fpOutput);
                }

        /* 6: sitewise parameters hi */
        for (int i = 0; i < ali->nSites; i++)
            for (int ai = 0; ai < ali->nCodes; ai++) {
                OUTPUT_PRECISION h = (OUTPUT_PRECISION) xHi(i, ai);
                fwrite(&h, sizeof(h), 1, fpOutput);
            }

        fclose(fpOutput);
    } else {
        fprintf(stderr, "Error writing parameters\n");
        exit(1);
    }
}

void OutputParametersFull(char *outputFile, const numeric_t *x,
    alignment_t *ali, options_t *options) {
/* File format */
    FILE *fpOutput = NULL;
    fpOutput = fopen(outputFile, "w");
    if (fpOutput != NULL) {
        /* 1: nSites */
        int32_t nSites = (int32_t) ali->nSites;
        fwrite(&nSites, sizeof(nSites), 1, fpOutput);

        /* 2: nCodes */
        int32_t nCodes = (int32_t) ali->nCodes;
        fwrite(&nCodes, sizeof(nCodes), 1, fpOutput);

        /* 3: nSeqs */
        int32_t nSeqs = (int32_t) ali->nSeqs;
        fwrite(&nSeqs, sizeof(nSeqs), 1, fpOutput);

        /* 4: nSkippedSeqs */
        int32_t nSkippedSeqs = (int32_t) ali->nSkippedSeqs;
        fwrite(&nSkippedSeqs, sizeof(nSkippedSeqs), 1, fpOutput);

        /* 5: number of iterations */
        int32_t maxIter = (int32_t) options->maxIter;
        fwrite(&maxIter, sizeof(maxIter), 1, fpOutput);

        /* 6: theta */
        OUTPUT_PRECISION theta = (OUTPUT_PRECISION) options->theta;
        fwrite(&theta, sizeof(theta), 1, fpOutput);

        /* 7: lambda for fields (lh) */
        OUTPUT_PRECISION lh = (OUTPUT_PRECISION) options->lambdaH;
        fwrite(&lh, sizeof(lh), 1, fpOutput);

        /* 8: lambda for couplings (le) */
        OUTPUT_PRECISION le = (OUTPUT_PRECISION) options->lambdaE;
        fwrite(&le, sizeof(le), 1, fpOutput);

        /* 9: group lambda for couplings (lg) */
        OUTPUT_PRECISION lg = (OUTPUT_PRECISION) options->lambdaGroup;
        fwrite(&lg, sizeof(lg), 1, fpOutput);

        /* 10: effective sample size (nEff) */
        OUTPUT_PRECISION nEff = (OUTPUT_PRECISION) ali->nEff;
        fwrite(&nEff, sizeof(nEff), 1, fpOutput);

        /* 11: alphabet */
        int isGapped = (options->estimatorMAP == INFER_MAP_PLM_GAPREDUCE);
        for (int i = 0; i < ali->nCodes; i++) {
            int8_t letter = (int8_t) ali->alphabet[i + isGapped];
            fwrite(&letter, sizeof(letter), 1, fpOutput);
        }

        /* 12: sequence number of neighbors (self included) */
        int skipix = 0, reducedix = 0;
        for (int s = 0; s < ali->nSeqs + ali->nSkippedSeqs; s++) {
            if (skipix < ali->nSkippedSeqs && s == ali->skippedSeqs[skipix]) {
                /* Skip skipped sequences */
                OUTPUT_PRECISION w = (OUTPUT_PRECISION) 0;
                fwrite(&w, sizeof(w), 1, fpOutput);
                skipix++;
            } else {
                numeric_t nNeighbors = ali->weights[reducedix];
                nNeighbors = 1.0 / (nNeighbors * options->scale);
                OUTPUT_PRECISION w = (OUTPUT_PRECISION) nNeighbors;
                fwrite(&w, sizeof(w), 1, fpOutput);
                reducedix++;
            }
        }

        /* 13: (Focus mode) target sequence */
        if (ali->target >= 0) {
            for (int i = 0; i < ali->nSites; i++) {
                int8_t c = (int8_t) ali->alphabet[seq(ali->target, i)];
                fwrite(&c, sizeof(c), 1, fpOutput);
            }
        } else {
            int8_t c = (int8_t) ali->alphabet[0];
            for (int i = 0; i < ali->nSites; i++)
                fwrite(&c, sizeof(c), 1, fpOutput);
        }

        /* 14: (Focus mode) offset map */
        if (ali->target >= 0) {
            for (int i = 0; i < ali->nSites; i++) {
                int32_t ix = (int32_t) ali->offsets[i];
                fwrite(&ix, sizeof(ix), 1, fpOutput);
            }
        } else {
            for (int i = 0; i < ali->nSites; i++) {
                int32_t ix = (int32_t) i + 1;
                fwrite(&ix, sizeof(ix), 1, fpOutput);
            }
        }

        /* 15: sitewise marginals fi */
        for (int i = 0; i < ali->nSites; i++)
            for (int ai = 0; ai < ali->nCodes; ai++) {
                OUTPUT_PRECISION f = (OUTPUT_PRECISION) fi(i, ai);
                fwrite(&f, sizeof(f), 1, fpOutput);
            }

        /* 16: sitewise parameters hi */
        for (int i = 0; i < ali->nSites; i++)
            for (int ai = 0; ai < ali->nCodes; ai++) {
                OUTPUT_PRECISION h = (OUTPUT_PRECISION) xHi(i, ai);
                fwrite(&h, sizeof(h), 1, fpOutput);
            }

        /* 17: pairwise marginals fij */
        for (int i = 0; i < ali->nSites - 1; i++)
            for (int j = i + 1; j < ali->nSites; j++)
                for (int ai = 0; ai < ali->nCodes; ai++)
                    for (int aj = 0; aj < ali->nCodes; aj++) {
                        OUTPUT_PRECISION f =
                            (OUTPUT_PRECISION) fij(i, j, ai, aj);
                        fwrite(&f, sizeof(f), 1, fpOutput);
                    }

        /* 18: couplings eij */
        for (int i = 0; i < ali->nSites - 1; i++)
            for (int j = i + 1; j < ali->nSites; j++)
                for (int ai = 0; ai < ali->nCodes; ai++)
                    for (int aj = 0; aj < ali->nCodes; aj++) {
                        OUTPUT_PRECISION e =
                            (OUTPUT_PRECISION) xEij(i, j, ai, aj);
                        fwrite(&e, sizeof(e), 1, fpOutput);
                    }
        fclose(fpOutput);
    } else {
        fprintf(stderr, "Error writing parameters\n");
        exit(1);
    }
}
#undef OUTPUT_PRECISION

void OutputCouplingScores(char *couplingsFile, const numeric_t *x,
    alignment_t *ali, options_t *options) {
    FILE *fpOutput = NULL;
    fpOutput = fopen(couplingsFile, "w");
    if (fpOutput != NULL) {
        /* Compute the norm of the coupling parameters between each pair */
        numeric_t *couplings =
        (numeric_t *) malloc((ali->nSites * (ali->nSites - 1) / 2)
                * sizeof(numeric_t));

        for (int i = 0; i < ali->nSites * (ali->nSites - 1) / 2;
            i++) couplings[i] = 0;
        for (int i = 0; i < ali->nSites - 1; i++)
            for (int j = i + 1; j < ali->nSites; j++) {
                /* Norm(eij) over ai, aj */
                numeric_t norm = 0.0;
                for (int ai = 0; ai < ali->nCodes; ai++)
                    for (int aj = 0; aj < ali->nCodes; aj++)
                        norm += xEij(i, j, ai, aj) * xEij(i, j, ai, aj);
                norm = sqrt(norm);
                coupling(i, j) = norm;
            }
        numeric_t nPairs =
            ((numeric_t) ((ali->nSites) * (ali->nSites - 1))) / 2.0;

        /* Remove first component of the norms (Average Product Correction) */
        if (!options->zeroAPC) {
            /* Determine the site-wise statistics of the norms */
            numeric_t C_avg = 0.0;
            numeric_t *C_pos_avg =
                (numeric_t *) malloc(ali->nSites * sizeof(numeric_t));
            for (int i = 0; i < ali->nSites; i++) {
                C_pos_avg[i] = 0.0;
            }
            for (int i = 0; i < ali->nSites - 1; i++) {
                for (int j = i + 1; j < ali->nSites; j++) {
                    C_pos_avg[i] +=
                        coupling(i, j) / (numeric_t) (ali->nSites - 1);
                    C_pos_avg[j] +=
                        coupling(i, j) / (numeric_t) (ali->nSites - 1);
                    C_avg += coupling(i, j) / nPairs;
                }
            }

            /* Remove the first component */
            for (int i = 0; i < ali->nSites - 1; i++)
                for (int j = i + 1; j < ali->nSites; j++)
                    coupling(i, j) =
                        coupling(i, j) - C_pos_avg[i] * C_pos_avg[j] / C_avg;    
        }

        /* Output scores */
        if (ali->target >= 0) {
            /* Focus mode */
            for (int i = 0; i < ali->nSites - 1; i++)
                for (int j = i + 1; j < ali->nSites; j++) {
                    char ai = (char) ali->alphabet[seq(ali->target, i)];
                    char aj = (char) ali->alphabet[seq(ali->target, j)];
                    fprintf(fpOutput, "%d %c %d %c 0 %f\n",
                        ali->offsets[i], ai, ali->offsets[j], aj,
                        coupling(i, j));
                }
        } else {
            for (int i = 0; i < ali->nSites - 1; i++)
                for (int j = i + 1; j < ali->nSites; j++)
                    fprintf(fpOutput, "%d - %d - 0 %f\n", i + 1, j + 1,
                        coupling(i, j));
        }

        fclose(fpOutput);
    } else {
        fprintf(stderr, "Error writing coupling scores\n");
        exit(1);
    }
}