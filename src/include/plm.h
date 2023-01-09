#ifndef PLM_H
#define PLM_H

#include "lbfgs.h"
#include <sys/time.h>

#ifdef USE_FLOAT
typedef float numeric_t;
#else
typedef double numeric_t;
#endif
typedef int letter_t;

/** 
 * Modes of inference
 */
enum {
    /* Maximum a posteriori (MAP) */
    INFER_MAP
};

/* Methods for MAP estimates  */
enum {
    /* Maximum Pseudolikelihood (PLM), site-parallelized */
    INFER_MAP_PLM,
    /* Maximum Pseudolikelihood (PLM), site-parallelized, no gaps */
    INFER_MAP_PLM_GAPREDUCE,
    /* Maximum Pseudolikelihood (PLM), sequence-parallelized parallelization */
    INFER_MAP_PLM_BLOCK,
    /* Maximum Pseudolikelihood (PLM), dropout-regularized */
    INFER_MAP_PLM_DROPOUT,
    /* NOT FULLY IMPLEMENTED: Minimum Probability Flow */
    INFER_MPF
};

/** 
 * Methods for regularization
 */
enum {
    /* L2 (Gaussian priors, lambda = 1/2 inverse variance) */
    REGULARIZE_L2
};

/**
 * User options for data processing & inference
 */
typedef struct {
    /* Alignment processing */
    char *target;
    char *alphabet;

    /* Method for inference */
    int usePairs;
    int estimator;
    int estimatorMAP;
    int maxIter;

    /* SGD options */
    int sgd;
    int sgdBatchSize;

    /* Sequence weights */
    int fastWeights;
    numeric_t theta;
    numeric_t scale;

    /* Regularization */
    numeric_t lambdaH;
    numeric_t lambdaE;
    numeric_t lambdaGroup;

    /* Iterative APC removal */
    int zeroAPC;
} options_t;

options_t *default_options();

/**
 * Multiple sequence alignment
 */
typedef struct {
    /* Alignment dimensions and sequence content */
    int nSeqs;
    int nSites;
    int nCodes;
    char *alphabet;
    char **names;
    letter_t *sequences;

    /* Sequence skipping */
    int nSkippedSeqs;
    int *skippedSeqs;

    /* Focus mode */
    int target;
    int *offsets;

    /* Sequence weights and statistics */
    numeric_t nEff;
    numeric_t *weights;
    numeric_t *fi;
    numeric_t *fij;

    /* Inference */
    int nParams;
    numeric_t negLogLk;
    struct timeval start;  // sys/time.h
} alignment_t;

/* Command-line entrypoint.
 * Could move this whole thing out to main.c instead, but keeping here for minimal changes.
 */
void run_plmc(char *alignFile, char* outputFile, char *couplingsFile,
    char *weightsFile, char *weightsOutputFile, options_t *options);

/* Loads a multiple sequence alignment and encodes it into a specified alphabet.
   Any sequences containing characters outside of the alphabet are discarded.
   By default, these routines are case-insenstive but columns containing
   lower-case characters in a target sequence can optionally be filtered out. 
   If a protein alphabet is used, '.' characters will be treated as '-'. 
 */
alignment_t *MSARead(char *alignFile, options_t *options);

/* Counts empirical sitewise(fi) and pairwise(fij) marginals of the alignment */
void MSACountMarginals(alignment_t *ali, options_t *options);

/* Frees alignment and options */
void MSAFree(alignment_t *ali, options_t *options);

/* Parameter output */
void OutputParametersSite(char *outputFile, const numeric_t *x,
    alignment_t *ali);
void OutputParametersFull(char *outputFile, const numeric_t *x,
    alignment_t *ali, options_t *options);
void OutputCouplingScores(char *couplingsFile, const numeric_t *x,
    alignment_t *ali, options_t *options);


/* File I/O */
#define BUFFER_SIZE 40960
#define fgetstr(str, fp)   {char *endPos;                                   \
                            if (fgets(str, BUFFER_SIZE, fp) != NULL) {      \
                                if ((endPos = strchr(str, '\n')) != NULL)   \
                                    *endPos = '\0';                         \
                                /* Trim trailing carriage returns */        \
                                while (str[strlen(str) - 1] == '\r')        \
                                       str[strlen(str) - 1] = '\0';         \
                                }                                           \
                            }

/* Memory schemes for model parameters and gradients */
#define xHi(i, Ai)             x[i + ali->nSites * (Ai)]
#define xEij(i, j, Ai, Aj)     x[ali->nSites * ali->nCodes + (i < j ? (((j) * (j - 1)/2 + i) * ali->nCodes * ali->nCodes + (Aj) * ali->nCodes + Ai) : (((i)*(i - 1)/2 + j) * ali->nCodes * ali->nCodes + (Ai) * ali->nCodes + Aj))]
#define dHi(i, Ai)             g[i + ali->nSites * (Ai)]
#define dEij(i, j, Ai, Aj)     g[ali->nSites * ali->nCodes + (i < j ? (((j) * (j - 1)/2 + i) * ali->nCodes * ali->nCodes + (Aj) * ali->nCodes + Ai) : (((i)*(i - 1)/2 + j) * ali->nCodes * ali->nCodes + (Ai) * ali->nCodes + Aj))]

#define wHi(w, i, Ai)           w[i + ali->nSites * (Ai)]
#define wEij(w, i, j, Ai, Aj)   w[ali->nSites * ali->nCodes + (i < j ? (((j) * (j - 1)/2 + i) * ali->nCodes * ali->nCodes + (Aj) * ali->nCodes + Ai) : (((i)*(i - 1)/2 + j) * ali->nCodes * ali->nCodes + (Ai) * ali->nCodes + Aj))]
#define wLambdaHi(w, i)         w[i]
#define wLambdaEij(w, i, j)     w[ali->nSites + (i < j ? ((j)*(j - 1)/2 + i) : ((i)*(i - 1)/2 + j))]

/* Memory schemes for site-parallelized conditional loglk calculations
   Layout performances (loop reordering included) vs xHi/xEij on Intel core i7:
        distal site, distal spin, local spin        ~2.4x  [below]
        local spin, distal site, distal spin        ~1.5x
 */
#define siteH(i, a)             Xi[a + ali->nCodes * (a + ali->nCodes * (i))]
#define siteE(j, ai, aj)        Xi[ai + ali->nCodes * (aj + ali->nCodes * (j))]
#define siteDH(i, a)            Di[a + ali->nCodes * (a + ali->nCodes * (i))]
#define siteDE(j, ai, aj)       Di[ai + ali->nCodes * (aj + ali->nCodes * (j))]

/* Memory schemes for sequence-parallelized conditional loglk calculations */
#define Hp(i, ai)               H[ai + ali->nCodes * (i)]
#define Hi(i, ai)               hi[ai + ali->nCodes * (i)]
#define gHi(i, ai)              gHi[ai + ali->nCodes * (i)]
#define Eij(i, ai, j, aj)       eij[aj + ali->nCodes * (j + ali->nSites * (ai + ali->nCodes * (i)))]
#define gEij(i, ai, j, aj)      gEij[aj + ali->nCodes * (j + ali->nSites * (ai + ali->nCodes * (i)))]

/* Memtory scheme for dropout regularization */
#define bitHi(i, Ai)            drop_mask[i + ali->nSites * (Ai)]
#define bitEij(i, j, Ai, Aj)    drop_mask[ali->nSites * ali->nCodes \
                                            + (i < j ? (((j)*(j - 1)/2 + i) * ali->nCodes * ali->nCodes + (Aj) * ali->nCodes + Ai) \
                                                     : (((i)*(i - 1)/2 + j) * ali->nCodes * ali->nCodes + (Ai) * ali->nCodes + Aj))]

/* Memtory scheme for L2 regularization */
#define lambdaHi(i)             lambdas[i]
#define lambdaEij(i,j)          lambdas[ali->nSites + (i < j ? ((j)*(j - 1)/2 + i) : ((i)*(i - 1)/2 + j))]
#define gLambdaHi(i)            gLambdas[i]
#define gLambdaEij(i,j)         gLambdas[ali->nSites + (i < j ? ((j)*(j - 1)/2 + i) : ((i)*(i - 1)/2 + j))]


/* Coupling scores */
#define coupling(i,j)           couplings[(i < j ? ((j)*(j - 1)/2 + i) : ((i)*(i - 1)/2 + j))]

/* Memory scheme for model parameters and gradients */
// #define seq(s, i)               ali->sequences[s + i * ali->nSeqs]
#define seq(s, i)               ali->sequences[i + (s) * ali->nSites]
#define fi(i, Ai)               ali->fi[i + ali->nSites * (Ai)]
#define fij(i, j, Ai, Aj)       ali->fij[(i < j ? (((j)*(j - 1)/2 + i) * ali->nCodes * ali->nCodes + (Aj) * ali->nCodes + Ai) : (((i)*(i - 1)/2 + j) * ali->nCodes * ali->nCodes + (Ai) * ali->nCodes + Aj))]
#define M(s, i, m)              membership_matrix[s + ali->nSeqs * (i + ali->nSites * (m))]
#define g_ij(s, i, m)           g_ij[s + ali->nSeqs * (i + ali->nSites * (m))]

#endif /* PLM_H */
