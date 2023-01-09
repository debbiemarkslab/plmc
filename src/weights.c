#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include "include/weights.h"
#include "include/twister.h"

#ifdef USE_FLOAT
char* weights_fmt = "%0.8e\n";
#else
char* weights_fmt = "%0.16e\n";
#endif

void MSAReweightSequences(alignment_t *ali, options_t *options) {
    /* Reweight seqeuences by their inverse neighborhood size. Each sequence's
       weight is the inverse of the number of neighboring sequences with less
       than THETA percent divergence
    */
    for (int i = 0; i < ali->nSeqs; i++) ali->weights[i] = 1.0;

    /* Only apply reweighting if theta is on [0,1] */
    if (options->theta >= 0 && options->theta <= 1) {
        /* The neighborhood size of each sequence is the number of sequences
           in the alignment within theta percent divergence */

        if (options->fastWeights > 0 && options->fastWeights < ali->nSeqs) {
            /* Cluster the sequences with k-consensus */
            int nClusters = options->fastWeights;
            int nIterations = 10;
            int nSeqs = ali->nSeqs;
            int nCodes = ali->nCodes;
            int nSites = ali->nSites;

#define COUNTS(i,j,a) counts[i * nSites * nCodes + j * nCodes + a]
#define CONSENSUS(i,j) consensus[i * nSites + j]
#define ALI(i,j) aliPermute[i * nSites + j]

            /* Pick initial clusters with Reservoir sampling */
            int *clusters = (int *) malloc(nClusters * sizeof(int));
            letter_t *consensus = (letter_t *) malloc(nClusters * nSites * sizeof(letter_t));
            for (int i = 0; i < nClusters; i++) clusters[i] = i;
            for (int i = nClusters; i < nSeqs; i++) {
                int ix = genrand_int32() % (i);
                if (ix < nClusters) clusters[ix] = i;
            }
            for (int i = 0; i < nClusters; i++)
                for (int j = 0; j < nSites; j++)
                    CONSENSUS(i,j) = seq(clusters[i], j);
            free(clusters);

            /* EM steps */
            int *assignment = (int *) malloc(nSeqs * sizeof(int));
            int *counts = (int *) malloc(nClusters * nSites * nCodes * sizeof(int));
            int *radii = (int *) malloc(nClusters * sizeof(int));
            for (int i = 0; i < nSeqs; i++) assignment[i] = 0;
            fprintf(stderr, "Clustering");
            for (int t = 0; t < nIterations; t++) {
                fprintf(stderr, ".");
                /* Step 1. Update the assignments */
                for (int i = 0; i < nClusters; i++) radii[i] = 0;
#pragma omp parallel for
                for (int s = 0; s < nSeqs; s++) {
                    int ixOld = assignment[s];
                    /* Current distance to current assignment */
                    numeric_t distance = 0;
                    for (int j = 0; j < nSites; j++)
                        distance += (CONSENSUS(ixOld, j) != seq(s, j));
                    /* Find closest */
                    int ixNew = ixOld;
                    for (int i = 0; i < nClusters; i++) {
                        numeric_t distanceI = 0;
                        for (int j = 0; j < nSites; j++)
                            distanceI += (CONSENSUS(i, j) != seq(s, j));
                        if (distanceI < distance) {
                            ixNew = i;
                            distance = distanceI;
                        }
                    }
                    if (ixNew != ixOld) assignment[s] = ixNew;
                    if (radii[ixNew] < distance) radii[ixNew] = distance;
                }
                /* --------------------------_DEBUG_--------------------------*/
                // for (int s = 0; s < nClusters; s++) {
                //     int size = 0;
                //     for (int i = 0; i < nSeqs; i++) size += (assignment[i] == s);
                //     fprintf(stderr, ">Cluster %d, %d members, radius %d\n", s, size, radii[s]);
                //     for (int i = 0; i < ali->nSites; i++)
                //         if (CONSENSUS(s,i) >= 0) {
                //             fprintf(stderr, "%c", ali->alphabet[CONSENSUS(s,i)]);
                //         } else {
                //             fprintf(stderr, " ");
                //         }
                //     fprintf(stderr, "\n");
                // }
                /* --------------------------^DEBUG^--------------------------*/

                /* Step 2. Update the consensus sequences */
                /* Update the counts */
                if (t < nIterations - 1) {
                    for (int i = 0; i < nClusters * nSites * nCodes; i++)
                        counts[i] = 0;
                    for (int s = 0; s < nSeqs; s++)
                        for (int j = 0; j < nSites; j++)
                            COUNTS(assignment[s], j, seq(s, j)) += 1;
#pragma omp parallel for
                    for (int i = 0; i < nClusters; i++)
                        for (int j = 0; j < nSites; j++) {
                            int topCode = 0;
                            int topCounts = COUNTS(i, j, 0);
                            for (int b = 1; b < nCodes; b++)
                                if (COUNTS(i, j, b) > topCounts) {
                                    topCode = b;
                                    topCounts = COUNTS(i, j, b);
                                }
                            CONSENSUS(i ,j) = topCode;
                        }
                }
            }
            fprintf(stderr, "\n");

            /* Profile-profile distances */
            numeric_t *clusterID = (numeric_t *) malloc(nClusters * nClusters * sizeof(numeric_t));
            for (int i = 0; i < nClusters * nClusters; i++) clusterID[i] = 0;
#pragma omp parallel for
            for (int pi = 0; pi < nClusters; pi++)
                for (int pj = 0; pj < nClusters; pj++)
                    for (int j = 0; j < nSites; j++)
                        clusterID[pi + pj * nClusters] +=
                                (CONSENSUS(pi,j) == CONSENSUS(pj,j));

            free(consensus);
            free(counts);

            /* Permute alignment */
            int *clusterSizes = (int *) malloc(nClusters * sizeof(int));
            int *clusterStart = (int *) malloc(nClusters * sizeof(int));
            int *clusterEnd = (int *) malloc(nClusters * sizeof(int));
            int *permuteMap = (int *) malloc(nSeqs * sizeof(int));
            numeric_t *weightsP = (numeric_t *) malloc(nSeqs * sizeof(numeric_t));
            letter_t *aliPermute = (letter_t *) malloc(nSeqs * nSites * sizeof(letter_t));
            for (int i = 0; i < nClusters; i++) clusterSizes[i] = 0;
            for (int s = 0; s < ali->nSeqs; s++)
                clusterSizes[assignment[s]] += 1;
            int ix = 0;
            for (int i = 0; i < nClusters; i++) {
                clusterStart[i] = ix;
                ix += clusterSizes[i];
                clusterEnd[i] = ix;
            }
            ix = 0;
            for (int i = 0; i < nClusters; i++)
                for (int s = 0; s < ali->nSeqs; s++)
                    if (assignment[s] == i) {
                        for (int j = 0; j < nSites; j++)
                            ALI(ix,j) = seq(s,j);
                        permuteMap[ix] = s;
                        ix++;
                    }

            /* ----------------------------_DEBUG_----------------------------*/
            // for (int s = 0; s < nSeqs; s++) {
            //     fprintf(stdout, ">Seq %d\n", s);
            //     for (int i = 0; i < ali->nSites; i++)
            //             fprintf(stdout, "%c", ali->alphabet[ALI(s,i)]);
            //     fprintf(stdout, "\n");
            // }
            /* ----------------------------^DEBUG^----------------------------*/

            /* Sequence weights */
            numeric_t cutoff = (numeric_t) ((1 - options->theta) * ali->nSites);
            for (int s = 0; s < nSeqs; s++) weightsP[s] = 1;
#pragma omp parallel for
            for (int ci = 0; ci < nClusters; ci++)
                for (int cj = 0; cj < nClusters; cj++)
                    if (clusterID[ci * nClusters + cj] >= 0.9 * cutoff)
                        for (int s = clusterStart[ci]; s < clusterEnd[ci]; s++)
                            for (int t = clusterStart[cj]; t < clusterEnd[cj]; t++)
                                if (s != t) {
                                    int id = 0;
                                    for (int n = 0; n < ali->nSites; n++)
                                        id += (ALI(s, n) == ALI(t, n));
                                    if (id >= cutoff) weightsP[s] += 1.0;
                                }
            for (int s = 0; s < nSeqs; s++)
                ali->weights[permuteMap[s]] = weightsP[s];

#undef COUNTS
#undef CONSENSUS
#undef ALI

            free(clusterSizes);
            free(clusterStart);
            free(clusterEnd);
            free(permuteMap);
            free(weightsP);
            free(radii);
            free(aliPermute);
        } else {
            /* Deterministic sequence weights */
#if defined(_OPENMP)
            /* Naive parallelization is faster ignoring symmetry */
            #pragma omp parallel for
            for (int s = 0; s < ali->nSeqs; s++)
                for (int t = 0; t < ali->nSeqs; t++)
                    if (s != t) {
                        int id = 0;
                        for (int n = 0; n < ali->nSites; n++)
                            id += (seq(s, n) == seq(t, n));
                        if (id >= ((1 - options->theta) * ali->nSites))
                            ali->weights[s] += 1.0;
                    }
#else
            /* For a single core, take advantage of symmetry */
            for (int s = 0; s < ali->nSeqs - 1; s++)
                for (int t = s + 1; t < ali->nSeqs; t++) {
                    int id = 0;
                    for (int n = 0; n < ali->nSites; n++)
                        id += (seq(s, n) == seq(t, n));
                    if (id >= ((1 - options->theta) * ali->nSites)) {
                        ali->weights[s] += 1.0;
                        ali->weights[t] += 1.0;
                    }
                }
#endif
        }

        /* Reweight sequences by the inverse of the neighborhood size */
        for (int i = 0; i < ali->nSeqs; i++)
            ali->weights[i] = 1.0 / ali->weights[i];
    }

    /* Scale sets the effective number of samples per neighborhood */
    for (int i = 0; i < ali->nSeqs; i++)
        ali->weights[i] *= options->scale;

    /* The effective number of sequences is then the sum of the weights */
    ali->nEff = 0;
    for (int i = 0; i < ali->nSeqs; i++) ali->nEff += ali->weights[i];

    if (options->theta >= 0 && options->theta <= 1) {
        fprintf(stderr,
                "Effective number of samples: %.1f\t(%.0f%% identical neighborhood = %.3f samples)\n",
                ali->nEff, 100 * (1 - options->theta), options->scale);
    } else {
        fprintf(stderr,
                "Theta not between 0 and 1, no sequence reweighting applied (N = %.2f)\n",
                ali->nEff);
    }
}

int ValidateCustomWeightsFile(char *weightsFile, alignment_t *ali) {
    /* Check that the weights file exists */
    /* Remember to close file pointer before returning. */
    FILE *fp = fopen(weightsFile, "r");
    if (fp == NULL) {
        fprintf(stderr, "Error: weights file %s does not exist\n", weightsFile);
        fclose(fp);
        return 1;
    }

    /* Count number of lines in file */
    int nLines = 0;
    int MAX_LINE_LENGTH = 1024;
    char line[MAX_LINE_LENGTH];
    while (fgets(line, MAX_LINE_LENGTH, fp) != NULL) nLines++;

    int nSeqsRaw = ali->nSeqs + ali->nSkippedSeqs;
    if (nLines != nSeqsRaw) {
        fprintf(stderr, "Error: weights file %s has %d lines, but alignment has %d sequences\n",
                weightsFile, nLines, nSeqsRaw);
        fclose(fp);
        return 1;
    }

    fclose(fp);
    return 0;
}

void ReadCustomWeightsFile(char *weightsFile, alignment_t *ali) {
    /* Note: Not using options->scale (or options->theta) for now (assuming this is done in original weights calc).*/
    /* Most of this is copied from MSAReweightSequences() */

    int validCode = ValidateCustomWeightsFile(weightsFile, ali);
    if (validCode != 0) {
        fprintf(stderr, "Error: weights file %s is invalid\n", weightsFile);
        exit(1);
    }

    /* Load weights (float array) into ali->weights and set ali->nEff */
    FILE *fp = fopen(weightsFile, "r");
    if (fp == NULL) {
        fprintf(stderr, "Error: could not open weights file %s\n", weightsFile);
        exit(1);
    }
    /* Reinitialize array just in case */
    for (int i = 0; i < ali->nSeqs; i++) ali->weights[i] = 1.0;

    /* Read weights, one float per line */
    int skippedIdx = 0, reducedIdx = 0, nWarnings = 0, maxWarnings = 64;
    int nSeqsTotal = ali->nSeqs + ali->nSkippedSeqs;
    for (int i = 0; i < nSeqsTotal; i++) {
        long double w;
        if (fscanf(fp, "%Lf", &w) != 1) {
            fprintf(stderr, "Error reading weights file %s at position %d\n", weightsFile, i);
            exit(1);
        }
        // Skip invalid sequence weights
        if ((skippedIdx < ali->nSkippedSeqs) && (i == ali->skippedSeqs[skippedIdx])) {
            if (w > 0 && nWarnings < maxWarnings) {
                fprintf(stderr, "Warning: Skipped nonzero weight in file %s at position %d\n", weightsFile, i);
                nWarnings++;
            }
            skippedIdx++;
            continue;
        } else {
            ali->weights[reducedIdx] = (numeric_t)w;
            reducedIdx++;
        }
    }
    fclose(fp);

    /* The effective number of sequences is then the sum of the weights */
    ali->nEff = 0;
    for (int i = 0; i < ali->nSeqs; i++) ali->nEff += ali->weights[i];

    fprintf(stderr,
            "Weights loaded successfully. Effective number of samples (to 1 decimal place): %.1f.\n",
            ali->nEff);
}

void WriteWeightsFile(char *weightsFile, alignment_t *ali) {
    // Note: Ignoring options->scale and options->theta here, writing out raw weights
    /* Write weights to file */
    FILE *fpOutput = fopen(weightsFile, "w");
    if (fpOutput == NULL) {
        fprintf(stderr, "Error: could not open weights file %s\n", weightsFile);
        exit(1);
    }
    // Write out weights, one float per line, and include invalid seqs as weight 0
    // Copied from plm.OutputParametersFull()
    int skipix = 0, reducedix = 0;
    for (int i = 0; i < (ali->nSeqs + ali->nSkippedSeqs); i++) {
        if (skipix < ali->nSkippedSeqs && i == ali->skippedSeqs[skipix]) {
            /* Skip skipped sequences */
            fprintf(fpOutput, "0.0\n");
            skipix++;
        } else {
            numeric_t w = ali->weights[reducedix];
            fprintf(fpOutput, weights_fmt, w);
            reducedix++;
        }
    }
    fclose(fpOutput);
}
