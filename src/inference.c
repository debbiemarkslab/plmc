#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <sys/time.h>
#include <assert.h>
#include <string.h>

/* Optionally include OpenMP with the -fopenmp flag */
#if defined(_OPENMP)
    #include <omp.h>
#endif

#include "include/lbfgs.h"
#include "include/plm.h"
#include "include/inference.h"

#define PI 3.14159265358979323846

/* Numerical bounds for ZeroAPCPriors */
#define LAMBDA_J_MIN 1E-2
#define LAMBDA_J_MAX 1E4
#define REGULARIZATION_GROUP_EPS 0.001

/* Internal to InferPairModel: 
    MAP estimation of parameters by L-BFGS */
void EstimatePairModelMAP(numeric_t *x, numeric_t *lambdas, alignment_t *ali,
    options_t *options);
/* Internal to EstimatePairModelMAP: 
   Objective functions for point parameter estimates (MAP) */
static lbfgsfloatval_t PLMNegLogPosterior(void *instance,
    const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n,
    const lbfgsfloatval_t step);
static lbfgsfloatval_t PLMNegLogPosteriorGapReduce(void *instance,
    const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n,
    const lbfgsfloatval_t step);
static lbfgsfloatval_t PLMNegLogPosteriorBlock(void *instance,
    const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n,
    const lbfgsfloatval_t step);
static lbfgsfloatval_t PLMNegLogPosteriorDO(void *instance,
    const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n,
    const lbfgsfloatval_t step);
/* Internal to EstimatePairModelMAP: progress reporting */
static int ReportProgresslBFGS(void *instance, const lbfgsfloatval_t *x,
    const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step, int n, int k, int ls);
/* Internal to EstimatePairModelMAP: parameter processing */
void PreCondition(const lbfgsfloatval_t *x, lbfgsfloatval_t *g,
    alignment_t *ali, options_t *options);
lbfgsfloatval_t PostCondition(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, lbfgsfloatval_t fx,
    alignment_t *ali, options_t *options);
void ZeroAPCPriors(alignment_t *ali, options_t *options, numeric_t *lambdas,
    lbfgsfloatval_t *x);
/* Internal to EstimatePairModelMAP: utility functions to L-BFGS */
const char *LBFGSErrorString(int ret);


numeric_t *InferPairModel(alignment_t *ali, options_t *options) {
    /* Estimate the parameters of a maximum entropy model for a
       multiple sequence alignment */

    /* Initialize the regularization parameters */
    numeric_t *lambdas =
    (numeric_t *) malloc((ali->nSites + ali->nSites * (ali->nSites - 1) / 2)
            * sizeof(numeric_t));
    for (int i = 0; i < ali->nSites; i++) lambdaHi(i) = options->lambdaH;
    for (int i = 0; i < ali->nSites - 1; i++)
        for (int j = i + 1; j < ali->nSites; j++)
            lambdaEij(i, j) = options->lambdaE;

    /* For gap-reduced problems, eliminate the gaps and reduce the alphabet */
    if (options->estimatorMAP == INFER_MAP_PLM_GAPREDUCE) {
        ali->nCodes = strlen(ali->alphabet) - 1;
        for (int i = 0; i < ali->nSites; i++)
            for (int s = 0; s < ali->nSeqs; s++)
                seq(s, i) -= 1;
    }

    /* Initialize parameters */
    ali->nParams = ali->nSites * ali->nCodes
        + ali->nSites * (ali->nSites - 1) / 2 * ali->nCodes * ali->nCodes;
    numeric_t *x = (numeric_t *) malloc(sizeof(numeric_t) * ali->nParams);
    if (x == NULL) {
        fprintf(stderr,
            "ERROR: Failed to allocate a memory block for variables.\n");
        exit(1);
    }
    for (int i = 0; i < ali->nParams; i++) x[i] = 0.0;

    /* Initialize site parameters with the ML estimates 
        hi = log(fi) + C
        A single pseudocount is added for stability 
       (Laplace's rule or Morcos et al. with lambda = nCodes) */
    if (options->zeroAPC != 1) {
        numeric_t pseudoC = (numeric_t) ali->nCodes;
        numeric_t Zinv = 1.0 / (ali->nEff + pseudoC);
        for (int i = 0; i < ali->nSites; i++)
            for (int ai = 0; ai < ali->nSites; ai++)
                xHi(i, ai) = Zinv * pseudoC / (numeric_t) ali->nCodes;
        for (int s = 0; s < ali->nSeqs; s++)
            for (int i = 0; i < ali->nSites; i++)
                xHi(i, seq(s, i)) += ali->weights[s] * Zinv;
        for (int i = 0; i < ali->nSites; i++)
            for (int ai = 0; ai < ali->nCodes; ai++)
                xHi(i, ai) = log(xHi(i, ai));
        /* Zero-sum gauge */
        for (int i = 0; i < ali->nSites; i++) {
            numeric_t hSum = 0.0;
            for (int ai = 0; ai < ali->nCodes; ai++) hSum += xHi(i, ai);
            numeric_t hShift = hSum / (numeric_t) ali->nCodes;
            for (int ai = 0; ai < ali->nCodes; ai++)
                xHi(i, ai) -= hShift;
        }
    }

    switch(options->estimator) {
        /* Point estimates */
        case INFER_MAP:
            /* Maximum a posteriori estimates of model parameters */
            EstimatePairModelMAP(x, lambdas, ali, options);
            break;
        /* For: future alternative estimators */
        default:
            /* Maximum a posteriori estimates of model parameters */
            EstimatePairModelMAP(x, lambdas, ali, options);
    }

    /* Restore the alignment encoding after inference */
    if (options->estimatorMAP == INFER_MAP_PLM_GAPREDUCE) {
        for (int i = 0; i < ali->nSites; i++)
            for (int s = 0; s < ali->nSeqs; s++)
                seq(s, i) += 1;
    }

    return (numeric_t *) x;
}

void EstimatePairModelMAP(numeric_t *x, numeric_t *lambdas, alignment_t *ali,
    options_t *options) {
    /* Computes Maximum a posteriori (MAP) estimates for the parameters of 
       and undirected graphical model by L-BFGS */

    /* Start timer */
    gettimeofday(&ali->start, NULL);

    /* Initialize L-BFGS */
    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    param.epsilon = 1E-3;
    param.max_iterations = options->maxIter; /* 0 is unbounded */

    /* Array of void pointers provides relevant data structures */
    void *d[3] = {(void *)ali, (void *)options, (void *)lambdas};

    /* Estimate parameters by optimization */
    static lbfgs_evaluate_t algo;
    switch(options->estimatorMAP) {
        case INFER_MAP_PLM:
            algo = PLMNegLogPosterior;
            break;
        case INFER_MAP_PLM_GAPREDUCE:
            algo = PLMNegLogPosteriorGapReduce;
            break;
        case INFER_MAP_PLM_BLOCK:
            algo = PLMNegLogPosteriorBlock;
            break;
        case INFER_MAP_PLM_DROPOUT:
            algo = PLMNegLogPosteriorDO;
            break;
        default:
            algo = PLMNegLogPosterior;
    }

    if (options->zeroAPC == 1) fprintf(stderr,
            "Estimating coupling hyperparameters le = 1/2 inverse variance\n");

    int ret = 0;
    lbfgsfloatval_t fx;
    ret = lbfgs(ali->nParams, x, &fx, algo, ReportProgresslBFGS,
        (void*)d, &param);
    fprintf(stderr, "Gradient optimization: %s\n", LBFGSErrorString(ret));

    /* Optionally re-estimate parameters with adjusted hyperparameters */
    if (options->zeroAPC == 1) {
        /* Form new priors on the variances */
        ZeroAPCPriors(ali, options, lambdas, x);

        /* Reinitialize coupling parameters */
        for (int i = 0; i < ali->nSites - 1; i++)
            for (int j = i + 1; j < ali->nSites; j++)
                for (int ai = 0; ai < ali->nCodes; ai++)
                    for (int aj = 0; aj < ali->nCodes; aj++)
                        xEij(i, j, ai, aj) = 0.0;

        /* Iterate estimation with new hyperparameter estimates */
        options->zeroAPC = 2;
        ret = lbfgs(ali->nParams, x, &fx, algo,
            ReportProgresslBFGS, (void*)d, &param);
        fprintf(stderr, "Gradient optimization: %s\n", LBFGSErrorString(ret));
    }
}

static lbfgsfloatval_t PLMNegLogPosterior(void *instance,
    const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n,
    const lbfgsfloatval_t step) {
    /* Compute the the negative log posterior, which is the negative 
       penalized log-(pseudo)likelihood and the objective for MAP inference
    */
    void **d = (void **)instance;
    alignment_t *ali = (alignment_t *) d[0];
    options_t *options = (options_t *) d[1];
    numeric_t *lambdas = (numeric_t *) d[2];

    /* Initialize log-likelihood and gradient */
    lbfgsfloatval_t fx = 0.0;
    for (int i = 0; i < ali->nParams; i++) g[i] = 0;

    /* Negative log-pseudolikelihood */
    #pragma omp parallel for
    for (int i = 0; i < ali->nSites; i++) {
        numeric_t *H = (numeric_t *) malloc(ali->nCodes * sizeof(numeric_t));
        numeric_t *P = (numeric_t *) malloc(ali->nCodes * sizeof(numeric_t));

        numeric_t siteFx = 0.0;
        /* Reshape site parameters and gradient into local blocks */
        numeric_t *Xi = (numeric_t *) malloc(ali->nCodes * ali->nCodes
            * ali->nSites * sizeof(numeric_t));
        for (int j = 0; j < i; j++)
            for (int a = 0; a < ali->nCodes; a++)
                for (int b = 0; b < ali->nCodes; b++)
                    siteE(j, a, b) = xEij(i, j, a, b);
        for (int j = i + 1; j < ali->nSites; j++)
            for (int a = 0; a < ali->nCodes; a++)
                for (int b = 0; b < ali->nCodes; b++)
                    siteE(j, a, b) = xEij(i, j, a, b);
        for (int a = 0; a < ali->nCodes; a++) siteH(i, a) = xHi(i, a);

        numeric_t *Di = (numeric_t *) malloc(ali->nCodes * ali->nCodes
        * ali->nSites * sizeof(numeric_t));
        for (int d = 0; d < ali->nCodes * ali->nCodes * ali->nSites; d++)
            Di[d] = 0.0;

        /* Site negative conditional log likelihoods */
        for (int s = 0; s < ali->nSeqs; s++) {
            /* Compute potentials */
            for (int a = 0; a < ali->nCodes; a++) H[a] = siteH(i, a);
            for (int j = 0; j < i; j++)
                for (int a = 0; a < ali->nCodes; a++)
                    H[a] += siteE(j, a, seq(s, j));
            for (int j = i + 1; j < ali->nSites; j++)
                for (int a = 0; a < ali->nCodes; a++)
                    H[a] += siteE(j, a, seq(s, j));

            /* Conditional distribution given sequence background */
            numeric_t scale = H[0];
            for (int a = 1; a < ali->nCodes; a++)
                scale = (scale >= H[a] ? scale : H[a]);
            for (int a = 0; a < ali->nCodes; a++) P[a] = exp(H[a] - scale);
            numeric_t Z = 0;
            for (int a = 0; a < ali->nCodes; a++) Z += P[a];
            numeric_t Zinv = 1.0 / Z;
            for (int a = 0; a < ali->nCodes; a++) P[a] *= Zinv;


            /* Log-likelihood contributions are scaled by sequence weight */
            numeric_t w = ali->weights[s];	
            siteFx -= w * log(P[seq(s, i)]);

            /* Field gradient */
            siteDH(i, seq(s, i)) -= w;
            for (int a = 0; a < ali->nCodes; a++)
                siteDH(i, a) -= -w * P[a];

            /* Couplings gradient */
            int ix = seq(s, i);
            for (int j = 0; j < i; j++)
                siteDE(j, ix, seq(s, j)) -= w;
            for (int j = i + 1; j < ali->nSites; j++)
                siteDE(j, ix, seq(s, j)) -= w;
            for (int j = 0; j < i; j++)
                for (int a = 0; a < ali->nCodes; a++)
                    siteDE(j, a, seq(s, j)) -= -w * P[a];
            for (int j = i + 1; j < ali->nSites; j++)
                for (int a = 0; a < ali->nCodes; a++)
                    siteDE(j, a, seq(s, j)) -= -w * P[a];
        }

        /* Contribute local loglk and gradient to global */
        #pragma omp critical
        {
        fx += siteFx;
        for (int j = 0; j < i; j++)
            for (int a = 0; a < ali->nCodes; a++)
                for (int b = 0; b < ali->nCodes; b++)
                    dEij(i, j, a, b) += siteDE(j, a, b);
        for (int j = i + 1; j < ali->nSites; j++)
            for (int a = 0; a < ali->nCodes; a++)
                for (int b = 0; b < ali->nCodes; b++)
                    dEij(i, j, a, b) += siteDE(j, a, b);
        for (int a = 0; a < ali->nCodes; a++) dHi(i, a) += siteDH(i, a);
        free(Xi);
        free(Di);
        }

        free(H);
        free(P);
    }

    ali->negLogLk = fx;

    /* Gaussian priors */
    for (int i = 0; i < ali->nSites; i++)
        for (int ai = 0; ai < ali->nCodes; ai++) {
            dHi(i, ai) += lambdaHi(i) * 2.0 * xHi(i, ai);
            fx += lambdaHi(i) * xHi(i, ai) * xHi(i, ai);
        }

    for (int i = 0; i < ali->nSites-1; i++)
        for (int j = i + 1; j < ali->nSites; j++)
            for (int ai = 0; ai < ali->nCodes; ai++)
                for (int aj = 0; aj < ali->nCodes; aj++) {
                    dEij(i, j, ai, aj) += lambdaEij(i, j)
                        * 2.0 * xEij(i, j, ai, aj);
                    fx += lambdaEij(i, j)
                        * xEij(i, j, ai, aj) * xEij(i, j, ai, aj);
                }

    fx = PostCondition(x, g, fx, ali, options);
    return fx;
}

static lbfgsfloatval_t PLMNegLogPosteriorGapReduce(void *instance,
    const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n,
    const lbfgsfloatval_t step) {
    /* Compute the the negative log posterior, which is the negative 
       penalized log-(pseudo)likelihood and the objective for MAP inference
    */
    void **d = (void **)instance;
    alignment_t *ali = (alignment_t *) d[0];
    options_t *options = (options_t *) d[1];
    numeric_t *lambdas = (numeric_t *) d[2];

    /* Initialize log-likelihood and gradient */
    lbfgsfloatval_t fx = 0.0;
    for (int i = 0; i < ali->nParams; i++) g[i] = 0;

    /* Negative log-pseudolikelihood */
    #pragma omp parallel for
    for (int i = 0; i < ali->nSites; i++) {
        numeric_t *H = (numeric_t *) malloc(ali->nCodes * sizeof(numeric_t));
        numeric_t *P = (numeric_t *) malloc(ali->nCodes * sizeof(numeric_t));

        numeric_t siteFx = 0.0;
        /* Reshape site parameters and gradient into local blocks */
        numeric_t *Xi = (numeric_t *) malloc(ali->nCodes * ali->nCodes
            * ali->nSites * sizeof(numeric_t));
        for (int j = 0; j < i; j++)
            for (int a = 0; a < ali->nCodes; a++)
                for (int b = 0; b < ali->nCodes; b++)
                    siteE(j, a, b) = xEij(i, j, a, b);
        for (int j = i + 1; j < ali->nSites; j++)
            for (int a = 0; a < ali->nCodes; a++)
                for (int b = 0; b < ali->nCodes; b++)
                    siteE(j, a, b) = xEij(i, j, a, b);
        for (int a = 0; a < ali->nCodes; a++) siteH(i, a) = xHi(i, a);

        numeric_t *Di = (numeric_t *) malloc(ali->nCodes * ali->nCodes
        * ali->nSites * sizeof(numeric_t));
        for (int d = 0; d < ali->nCodes * ali->nCodes * ali->nSites; d++)
            Di[d] = 0.0;

        /* Site negative conditional log likelihoods */
        for (int s = 0; s < ali->nSeqs; s++) {
            /* Only ungapped sites are considered in the model */
            if (seq(s, i) >= 0) {
                /* Compute potentials */
                for (int a = 0; a < ali->nCodes; a++) H[a] = siteH(i, a);
                for (int j = 0; j < i; j++)
                    for (int a = 0; a < ali->nCodes; a++)
                        if (seq(s, j) >= 0)
                            H[a] += siteE(j, a, seq(s, j));
                for (int j = i + 1; j < ali->nSites; j++)
                    for (int a = 0; a < ali->nCodes; a++)
                        if (seq(s, j) >= 0)
                            H[a] += siteE(j, a, seq(s, j));

                /* Conditional distribution given sequence background */
                numeric_t scale = H[0];
                for (int a = 1; a < ali->nCodes; a++)
                    scale = (scale >= H[a] ? scale : H[a]);
                for (int a = 0; a < ali->nCodes; a++) P[a] = exp(H[a] - scale);
                numeric_t Z = 0;
                for (int a = 0; a < ali->nCodes; a++) Z += P[a];
                numeric_t Zinv = 1.0 / Z;
                for (int a = 0; a < ali->nCodes; a++) P[a] *= Zinv;


                /* Log-likelihood contributions are scaled by sequence weight */
                numeric_t w = ali->weights[s];  
                siteFx -= w * log(P[seq(s, i)]);

                /* Field gradient */
                siteDH(i, seq(s, i)) -= w;
                for (int a = 0; a < ali->nCodes; a++)
                    siteDH(i, a) -= -w * P[a];

                /* Couplings gradient */
                int ix = seq(s, i);
                for (int j = 0; j < i; j++)
                    if (seq(s, j) >= 0)
                        siteDE(j, ix, seq(s, j)) -= w;
                for (int j = i + 1; j < ali->nSites; j++)
                    if (seq(s, j) >= 0)
                        siteDE(j, ix, seq(s, j)) -= w;
                for (int j = 0; j < i; j++)
                    if (seq(s, j) >= 0)
                        for (int a = 0; a < ali->nCodes; a++)
                            siteDE(j, a, seq(s, j)) -= -w * P[a];
                for (int j = i + 1; j < ali->nSites; j++)
                    if (seq(s, j) >= 0)
                        for (int a = 0; a < ali->nCodes; a++)
                            siteDE(j, a, seq(s, j)) -= -w * P[a];
            }
        }

        /* Contribute local loglk and gradient to global */
        #pragma omp critical
        {
        fx += siteFx;
        for (int j = 0; j < i; j++)
            for (int a = 0; a < ali->nCodes; a++)
                for (int b = 0; b < ali->nCodes; b++)
                    dEij(i, j, a, b) += siteDE(j, a, b);
        for (int j = i + 1; j < ali->nSites; j++)
            for (int a = 0; a < ali->nCodes; a++)
                for (int b = 0; b < ali->nCodes; b++)
                    dEij(i, j, a, b) += siteDE(j, a, b);
        for (int a = 0; a < ali->nCodes; a++) dHi(i, a) += siteDH(i, a);
        free(Xi);
        free(Di);
        }

        free(H);
        free(P);
    }

    ali->negLogLk = fx;

    /* Gaussian priors */
    for (int i = 0; i < ali->nSites; i++)
        for (int ai = 0; ai < ali->nCodes; ai++) {
            dHi(i, ai) += lambdaHi(i) * 2.0 * xHi(i, ai);
            fx += lambdaHi(i) * xHi(i, ai) * xHi(i, ai);
        }

    for (int i = 0; i < ali->nSites-1; i++)
        for (int j = i + 1; j < ali->nSites; j++)
            for (int ai = 0; ai < ali->nCodes; ai++)
                for (int aj = 0; aj < ali->nCodes; aj++) {
                    dEij(i, j, ai, aj) += lambdaEij(i, j)
                        * 2.0 * xEij(i, j, ai, aj);
                    fx += lambdaEij(i, j)
                        * xEij(i, j, ai, aj) * xEij(i, j, ai, aj);
                }

    fx = PostCondition(x, g, fx, ali, options);
    return fx;
}

static lbfgsfloatval_t PLMNegLogPosteriorBlock(void *instance,
    const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n,
    const lbfgsfloatval_t step) {
    /* Compute the the negative log posterior, which is the negative 
       penalized log-(pseudo)likelihood and the objective for MAP inference
    */
    void **d = (void **)instance;
    alignment_t *ali = (alignment_t *) d[0];
    options_t *options = (options_t *) d[1];
    numeric_t *lambdas = (numeric_t *) d[2];

    /* Initialize log-likelihood and gradient */
    lbfgsfloatval_t fx = 0.0;
    for (int i = 0; i < ali->nParams; i++) g[i] = 0;

    /* Block fields hi */
    numeric_t *hi = (numeric_t *)
        malloc(ali->nSites * ali->nCodes * sizeof(numeric_t));
    numeric_t *gHi = (numeric_t *)
        malloc(ali->nSites * ali->nCodes * sizeof(numeric_t));
    for (int i = 0; i < ali->nSites; i++)
        for (int ai = 0; ai < ali->nCodes; ai++) Hi(i, ai) = xHi(i, ai);
    for (int i = 0; i < ali->nSites * ali->nCodes; i++) gHi[i] = 0;

    /* Block couplings eij */
    numeric_t *eij = (numeric_t *) malloc(ali->nSites * ali->nSites
        * ali->nCodes * ali->nCodes * sizeof(numeric_t));
    numeric_t *gEij = (numeric_t *) malloc(ali->nSites * ali->nSites
        * ali->nCodes * ali->nCodes * sizeof(numeric_t));
    for (int i = 0; i < ali->nSites * ali->nSites * ali->nCodes * ali->nCodes;
        i++) eij[i] = 0.0;
    for (int i = 0; i < ali->nSites * ali->nSites * ali->nCodes * ali->nCodes;
        i++) gEij[i] = 0.0;
    for (int i = 0; i < ali->nSites - 1; i++)
        for (int j = i + 1; j < ali->nSites; j++)
            for (int ai = 0; ai < ali->nCodes; ai++)
                for (int aj = 0; aj < ali->nCodes; aj++)
                    Eij(j, aj, i, ai) = Eij(i, ai, j, aj) = xEij(i, j, ai, aj);


    /* Negative log-pseudolikelihood */
    for (int s = 0; s < ali->nSeqs; s++) {
        /* Form potential for conditional log likelihoods at every site */
        numeric_t *H = (numeric_t *)
            malloc(ali->nCodes * ali->nSites * sizeof(numeric_t));
        numeric_t *Z = (numeric_t *) malloc(ali->nSites * sizeof(numeric_t));

        /* Initialize potentials with fields */
        // memcpy(H, hi, ali->nSites * ali->nCodes * sizeof(numeric_t));
        for(int jx = 0; jx < ali->nSites * ali->nCodes; jx++) H[jx] = hi[jx];

        /* Contribute coupling block due to i, ai */
        for (int i = 0; i < ali->nSites; i++) {
            const letter_t ai = seq(s, i);
            const numeric_t *jB = &(Eij(i, ai, 0, 0));
            for(int jx = 0; jx < ali->nSites * ali->nCodes; jx++)
                H[jx] += jB[jx];
        }

        /* Conditional log likelihoods */
        for (int i = 0; i < ali->nSites * ali->nCodes; i++) H[i] = exp(H[i]);
        for (int i = 0; i < ali->nSites; i++) Z[i] = 0;
        for (int i = 0; i < ali->nSites; i++)
            for (int ai = 0; ai < ali->nSites; ai++) Z[i] += Hp(i, ai);
        for (int i = 0; i < ali->nSites; i++)
            for (int ai = 0; ai < ali->nSites; ai++) Hp(i, ai) /= Z[i];

        numeric_t seqFx = 0;
        for (int i = 0; i < ali->nSites; i++)
            seqFx -= ali->weights[s] * log(Hp(i, seq(s, i)));

        for(int jx = 0; jx < ali->nSites * ali->nCodes; jx++)
            H[jx] *= -ali->weights[s];

        for (int i = 0; i < ali->nSites; i++)
            gHi(i, seq(s, i)) -= ali->weights[s];
        for(int jx = 0; jx < ali->nSites * ali->nCodes; jx++) gHi[jx] -= H[jx];

        for (int i = 0; i < ali->nSites - 1; i++)
            for (int j = i; j < ali->nSites; j++)
                gEij(i, seq(s, i), j, seq(s, j)) -= ali->weights[s];

        for (int i = 0; i < ali->nSites; i++) {
            const letter_t ai = seq(s, i);
            numeric_t *jgBlock = &(gEij(i, ai, 0, 0));
            for (int jx = 0; jx < ali->nSites * ali->nCodes; jx++)
                jgBlock[jx] -= H[jx];
        }

        free(H);
        free(Z);
        fx += seqFx;
    }

    for (int i = 0; i < ali->nSites; i++)
        for (int ai = 0; ai < ali->nCodes; ai++)
            dHi(i, ai) += gHi(i, ai);

    for (int i = 0; i < ali->nSites - 1; i++)
        for (int j = i + 1; j < ali->nSites; j++)
            for (int ai = 0; ai < ali->nCodes; ai++)
                for (int aj = 0; aj < ali->nCodes; aj++)
                    dEij(i, j, ai, aj) += gEij(j, aj, i, ai) + gEij(i, ai, j, aj);
    free(hi);
    free(gHi);
    free(eij);
    free(gEij);

    ali->negLogLk = fx;

    /* Gaussian priors */
    for (int i = 0; i < ali->nSites; i++)
        for (int ai = 0; ai < ali->nCodes; ai++) {
            dHi(i, ai) += lambdaHi(i) * 2.0 * xHi(i, ai);
            fx += lambdaHi(i) * xHi(i, ai) * xHi(i, ai);
        }

    for (int i = 0; i < ali->nSites-1; i++)
        for (int j = i + 1; j < ali->nSites; j++)
            for (int ai = 0; ai < ali->nCodes; ai++)
                for (int aj = 0; aj < ali->nCodes; aj++) {
                    dEij(i, j, ai, aj) += lambdaEij(i, j)
                        * 2.0 * xEij(i, j, ai, aj);
                    fx += lambdaEij(i, j)
                        * xEij(i, j, ai, aj) * xEij(i, j, ai, aj);
                }

    fx = PostCondition(x, g, fx, ali, options);
    return fx;
}

static lbfgsfloatval_t PLMNegLogPosteriorDO(void *instance,
    const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n,
    const lbfgsfloatval_t step) {
    /* Compute the the negative log posterior, which is the negative 
       penalized log-(pseudo)likelihood and the objective for MAP inference
    */
    void **d = (void **)instance;
    alignment_t *ali = (alignment_t *) d[0];
    options_t *options = (options_t *) d[1];
    numeric_t *lambdas = (numeric_t *) d[2];

    /* Initialize log-likelihood and gradient */
    lbfgsfloatval_t fx = 0.0;
    for (int i = 0; i < ali->nParams; i++) g[i] = 0;

    numeric_t *H = (numeric_t *) malloc(ali->nCodes * sizeof(numeric_t));
    numeric_t *P = (numeric_t *) malloc(ali->nCodes * sizeof(numeric_t));
    int *drop_mask = (int *) malloc(ali->nParams * sizeof(int));
    for (int s = 0; s < ali->nSeqs; s++) {
        /* Generate random bit mask over parameters */
        for (int p = 0; p < ali->nParams; p ++)
            drop_mask[p] = (int) rand() % 2;

        /* Pseudolikelihood objective */
        for (int i = 0; i < ali->nSites; i++) {
            for (int a = 0; a < ali->nCodes; a++) H[a] = bitHi(i, a)
                                               * xHi(i, a);
            for (int a = 0; a < ali->nCodes; a++)
                for (int j = 0; j < i; j++)
                    H[a] += bitEij(i, j, a, seq(s, j))
                            * xEij(i, j, a, seq(s, j));
            for (int a = 0; a < ali->nCodes; a++)
                for (int j = i + 1; j < ali->nSites; j++)
                    H[a] += bitEij(i, j, a, seq(s, j))
                            * xEij(i, j, a, seq(s, j));

            /* Compute distribution from potential */
            for (int a = 0; a < ali->nCodes; a++) P[a] = exp(H[a]);
            numeric_t Z = 0;
            for (int a = 0; a < ali->nCodes; a++) Z += P[a];
            numeric_t Zinv = 1.0 / Z;
            for (int a = 0; a < ali->nCodes; a++) P[a] *= Zinv;

            /* Log-likelihood contributions */
            fx -= ali->weights[s] * log(P[seq(s, i)]);

            /* Field gradient */
            dHi(i, seq(s, i)) -= bitHi(i, seq(s, i)) * ali->weights[s];
            for (int a = 0; a < ali->nCodes; a++)
                dHi(i, a) -= -bitHi(i, a) * ali->weights[s] * P[a];

            /* Couplings gradient */
            for (int j = 0; j < i; j++)
                dEij(i, j, seq(s, i), seq(s, j)) -=
                    bitEij(i, j, seq(s, i), seq(s, j)) * ali->weights[s];
            for (int j = i + 1; j < ali->nSites; j++)
                dEij(i, j, seq(s, i), seq(s, j)) -=
                    bitEij(i, j, seq(s, i), seq(s, j)) * ali->weights[s];

            for (int j = 0; j < i; j++)
                for (int a = 0; a < ali->nCodes; a++)
                    dEij(i, j, a, seq(s, j)) -=
                        -bitEij(i, j, a, seq(s, j)) * ali->weights[s] * P[a];
            for (int j = i + 1; j < ali->nSites; j++)
                for (int a = 0; a < ali->nCodes; a++)
                    dEij(i, j, a, seq(s, j)) -=
                        -bitEij(i, j, a, seq(s, j)) * ali->weights[s] * P[a];
        }
    }
    free(H);
    free(P);
    free(drop_mask);

    ali->negLogLk = fx;

    /* Gaussian priors */
    for (int i = 0; i < ali->nSites; i++)
        for (int ai = 0; ai < ali->nCodes; ai++) {
            dHi(i, ai) += lambdaHi(i) * 2.0 * xHi(i, ai);
            fx += lambdaHi(i) * xHi(i, ai) * xHi(i, ai);
        }

    for (int i = 0; i < ali->nSites-1; i++)
        for (int j = i + 1; j < ali->nSites; j++)
            for (int ai = 0; ai < ali->nCodes; ai++)
                for (int aj = 0; aj < ali->nCodes; aj++) {
                    dEij(i, j, ai, aj) += lambdaEij(i, j)
                        * 2.0 * xEij(i, j, ai, aj);
                    fx += lambdaEij(i, j)
                        * xEij(i, j, ai, aj) * xEij(i, j, ai, aj);
                }

    fx = PostCondition(x, g, fx, ali, options);
    return fx;
}

static int ReportProgresslBFGS(void *instance, const lbfgsfloatval_t *x,
    const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step, int n, int k, int ls) {
    void **d = (void **)instance;
    alignment_t *ali = (alignment_t *)d[0];

    /* Compute norms of relevant parameters */
    lbfgsfloatval_t hNorm = 0.0, eNorm = 0.0, hGNorm = 0.0, eGNorm = 0.0;
    for (int i = 0; i < ali->nSites * ali->nCodes; i++)
        hNorm += x[i]*x[i];
    for (int i = 0; i < ali->nSites * ali->nCodes; i++)
        hGNorm += g[i]*g[i];
    for (int i = ali->nSites * ali->nCodes; i < ali->nParams; i++)
        eNorm += x[i]*x[i];
    for (int i = ali->nSites * ali->nCodes; i < ali->nParams; i++)
        eGNorm += g[i]*g[i];
    hNorm = sqrt(hNorm);
    hGNorm = sqrt(hGNorm);
    eNorm = sqrt(eNorm);
    eGNorm = sqrt(eGNorm);

    /* Retrieve elapsed time */
    static struct timeval now;
    gettimeofday(&now, NULL);
    if (now.tv_usec < ali->start.tv_usec) {
        int nsec = (ali->start.tv_usec - now.tv_usec) / 1000000 + 1;
        ali->start.tv_usec -= 1000000 * nsec;
        ali->start.tv_sec += nsec;
    }
    if (now.tv_usec - ali->start.tv_usec > 1000000) {
        int nsec = (now.tv_usec - ali->start.tv_usec) / 1000000;
        ali->start.tv_usec += 1000000 * nsec;
        ali->start.tv_sec -= nsec;
    }
    numeric_t elapsed = (numeric_t) (now.tv_sec - ali->start.tv_sec)
                      + ((numeric_t) (now.tv_usec - ali->start.tv_usec)) / 1E6;

    if (k == 1) fprintf(stderr,
        "iter\ttime\tcond\tfx\t-loglk"
        "\t||h||\t||e||\n");
    fprintf(stderr, "%d\t%.1f\t%.2f\t%.1f\t%.1f\t%.1f\t%.1f\n",
        k, elapsed, gnorm / xnorm, fx, ali->negLogLk, hNorm, eNorm);
    return 0;
}

void PreCondition(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, alignment_t *ali, options_t *options) {
    /* Currently empty */
}

lbfgsfloatval_t PostCondition(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, lbfgsfloatval_t fx, alignment_t *ali, options_t *options) {
    if (options->zeroAPC == 1)
        for (int i = 0; i < ali->nSites; i++)
            for (int ai = 0; ai < ali->nCodes; ai++)
                dHi(i, ai) = 0.0;

    /* Group (L1/L2) regularization  */
    if (options->lambdaGroup > 0)
        for (int i = 0; i < ali->nSites - 1; i++)
            for (int j = i + 1; j < ali->nSites; j++) {
                double l2 = REGULARIZATION_GROUP_EPS;
                for (int ai = 0; ai < ali->nCodes; ai++)
                    for (int aj = 0; aj < ali->nCodes; aj++)
                        l2 += xEij(i, j, ai, aj) * xEij(i, j, ai, aj);
                double l1 = sqrt(l2);
                fx += options->lambdaGroup * l1;
                for (int ai = 0; ai < ali->nCodes; ai++)
                    for (int aj = 0; aj < ali->nCodes; aj++)
                        dEij(i, j, ai, aj) += options->lambdaGroup * xEij(i, j, ai, aj) / l1;
            }

    return fx;
}

void ZeroAPCPriors(alignment_t *ali, options_t *options, numeric_t *lambdas,
    lbfgsfloatval_t *x) {
    /* Compute the variances of the couplings for each pair */
    for (int i = 0; i < ali->nSites - 1; i++)
        for (int j = i + 1; j < ali->nSites; j++) {
            /* Mean(eij) over ai, aj */
            numeric_t mean = 0.0;
            for (int ai = 0; ai < ali->nCodes; ai++)
                for (int aj = 0; aj < ali->nCodes; aj++)
                    mean += xEij(i, j, ai, aj);
            mean *= 1.0 / ((numeric_t) ali->nCodes * ali->nCodes);

            /* Var(eij) over ai, aj */
            numeric_t ssq = 0.0;
            for (int ai = 0; ai < ali->nCodes; ai++)
                for (int aj = 0; aj < ali->nCodes; aj++)
                    ssq += (xEij(i, j, ai, aj) - mean)
                         * (xEij(i, j, ai, aj) - mean);
            /* Use N rather than N-1 since N has better MSE */
            numeric_t var = ssq / ((numeric_t) (ali->nCodes * ali->nCodes));
            lambdaEij(i, j) = var;
        }

    /* Determine the site-wise statistics of the variances */
    numeric_t nPairs =  ((numeric_t) ((ali->nSites) * (ali->nSites - 1))) / 2.0;
    numeric_t V_avg = 0.0;
    numeric_t *V_pos_avg = (numeric_t *) malloc(ali->nSites * sizeof(numeric_t));
    for (int i = 0; i < ali->nSites; i++) {
        V_pos_avg[i] = 0.0;
    }
    for (int i = 0; i < ali->nSites - 1; i++) {
        for (int j = i + 1; j < ali->nSites; j++) {
            V_pos_avg[i] += lambdaEij(i, j) / (numeric_t) (ali->nSites - 1);
            V_pos_avg[j] += lambdaEij(i, j) / (numeric_t) (ali->nSites - 1);
            V_avg += lambdaEij(i, j) / nPairs;
        }
    }

    /* Remove the first component of the variances */
    for (int i = 0; i < ali->nSites - 1; i++)
        for (int j = i + 1; j < ali->nSites; j++)
            lambdaEij(i, j) =
                lambdaEij(i, j) - V_pos_avg[i] * V_pos_avg[j] / V_avg;

    /* Transform and truncate variances into lambda hyperparameters */
    numeric_t pcount = 0.0;
    numeric_t psum = 0.0;
    numeric_t inbounds = 0;
    numeric_t min = LAMBDA_J_MAX;
    numeric_t max = LAMBDA_J_MIN;
    for (int i = 0; i < ali->nSites - 1; i++) {
        for (int j = i + 1; j < ali->nSites; j++) {
            /* Lambda coefficients are 1/2 the inverse variance */
            if (lambdaEij(i, j) > 0) {
                lambdaEij(i, j) = 1.0 / (2.0 * lambdaEij(i, j));
                psum += lambdaEij(i, j);
                pcount += 1.0;
            } else {
                lambdaEij(i, j) = LAMBDA_J_MAX + 1.0;
            }

            /* Truncate lambda for numerical stability */
            if (lambdaEij(i, j) >= LAMBDA_J_MIN && lambdaEij(i, j) <= LAMBDA_J_MAX)
                inbounds += 1.0 / (numeric_t) ((ali->nSites)*(ali->nSites - 1) / 2.0);
            if (lambdaEij(i, j) < 0 || !isfinite(lambdaEij(i, j)))
                lambdaEij(i, j) = LAMBDA_J_MAX;
            if (lambdaEij(i, j) < LAMBDA_J_MIN) lambdaEij(i, j) = LAMBDA_J_MIN;
            if (lambdaEij(i, j) > LAMBDA_J_MAX) lambdaEij(i, j) = LAMBDA_J_MAX;

            /* Track extremes */
            if (lambdaEij(i, j) > max) max = lambdaEij(i, j);
            if (lambdaEij(i, j) < min) min = lambdaEij(i, j);
        }
    }
    fprintf(stderr, "Raw coupling hyperparameter statistics:\n"
                    "\tMean positive lambda: %f\n"
                    "\tPercent of ij's positive: %f\n"
                    "\tPercent in bounds (%f < L < %f): %f\n",
                    psum / pcount,
                    pcount / nPairs,
                    min, max, inbounds);
}

const char *LBFGSErrorString(int ret) {
    const char *p;
    switch(ret) {
        case LBFGSERR_UNKNOWNERROR:
            p = "UNKNOWNERROR";
            break;
        /** Logic error. */
        case LBFGSERR_LOGICERROR:
            p = "LOGICERROR";
            break;
        /** Insufficient memory. */
        case LBFGSERR_OUTOFMEMORY:
            p = "OUTOFMEMORY";
            break;
        /** The minimization process has been canceled. */
        case LBFGSERR_CANCELED:
            p = "CANCELED";
            break;
        /** Invalid number of variables specified. */
        case LBFGSERR_INVALID_N:
            p = "INVALID_N";
            break;
        /** Invalid number of variables (for SSE) specified. */
        case LBFGSERR_INVALID_N_SSE:
            p = "INVALID_N_SSE";
            break;
        /** The array x must be aligned to 16 (for SSE). */
        case LBFGSERR_INVALID_X_SSE:
            p = "INVALID_X_SSE";
            break;
        /** Invalid parameter lbfgs_parameter_t::epsilon specified. */
        case LBFGSERR_INVALID_EPSILON:
            p = "INVALID_EPSILON";
            break;
        /** Invalid parameter lbfgs_parameter_t::past specified. */
        case LBFGSERR_INVALID_TESTPERIOD:
            p = "INVALID_TESTPERIOD";
            break;
        /** Invalid parameter lbfgs_parameter_t::delta specified. */
        case LBFGSERR_INVALID_DELTA:
            p = "INVALID_DELTA";
            break;
        /** Invalid parameter lbfgs_parameter_t::linesearch specified. */
        case LBFGSERR_INVALID_LINESEARCH:
            p = "INVALID_LINESEARCH";
            break;
        /** Invalid parameter lbfgs_parameter_t::max_step specified. */
        case LBFGSERR_INVALID_MINSTEP:
            p = "INVALID_MINSTEP";
            break;
        /** Invalid parameter lbfgs_parameter_t::max_step specified. */
        case LBFGSERR_INVALID_MAXSTEP:
            p = "INVALID_MAXSTEP";
            break;
        /** Invalid parameter lbfgs_parameter_t::ftol specified. */
        case LBFGSERR_INVALID_FTOL:
            p = "INVALID_FTOL";
            break;
        /** Invalid parameter lbfgs_parameter_t::wolfe specified. */
        case LBFGSERR_INVALID_WOLFE:
            p = "INVALID_WOLFE";
            break;
        /** Invalid parameter lbfgs_parameter_t::gtol specified. */
        case LBFGSERR_INVALID_GTOL:
            p = "INVALID_GTOL";
            break;
        /** Invalid parameter lbfgs_parameter_t::xtol specified. */
        case LBFGSERR_INVALID_XTOL:
            p = "INVALID_XTOL";
            break;
        /** Invalid parameter lbfgs_parameter_t::max_linesearch specified. */
        case LBFGSERR_INVALID_MAXLINESEARCH:
            p = "INVALID_MAXLINESEARCH";
            break;
        /** Invalid parameter lbfgs_parameter_t::orthantwise_c specified. */
        case LBFGSERR_INVALID_ORTHANTWISE:
            p = "INVALID_ORTHANTWISE";
            break;
        /** Invalid parameter lbfgs_parameter_t::orthantwise_start specified. */
        case LBFGSERR_INVALID_ORTHANTWISE_START:
            p = "INVALID_ORTHANTWISE_START";
            break;
        /** Invalid parameter lbfgs_parameter_t::orthantwise_end specified. */
        case LBFGSERR_INVALID_ORTHANTWISE_END:
            p = "ORTHANTWISE_END";
            break;
        /** The line-search step went out of the interval of uncertainty. */
        case LBFGSERR_OUTOFINTERVAL:
            p = "OUTOFINTERVAL";
            break;
        /** A logic error occurred; alternatively: the interval of uncertainty
            became too small. */
        case LBFGSERR_INCORRECT_TMINMAX:
            p = "INCORRECT_TMINMAX";
            break;
        /** A rounding error occurred; alternatively: no line-search step
            satisfies the sufficient decrease and curvature conditions. */
        case LBFGSERR_ROUNDING_ERROR:
            p = "ROUNDING_ERROR";
            break;
        /** The line-search step became smaller than lbfgs_parameter_t::min_step. */
        case LBFGSERR_MINIMUMSTEP:
            p = "MINIMUMSTEP";
            break;
        /** The line-search step became larger than lbfgs_parameter_t::max_step. */
        case LBFGSERR_MAXIMUMSTEP:
            p = "MAXILBFGSERR_MUMSTEP";
            break;
        /** The line-search routine reaches the maximum number of evaluations. */
        case LBFGSERR_MAXIMUMLINESEARCH:
            p = "MAXIMUMLINESEARCH";
            break;
        /** The algorithm routine reaches the maximum number of iterations. */
        case LBFGSERR_MAXIMUMITERATION:
            p = "MAXIMUMITERATION";
            break;
        /** Relative width of the interval of uncertainty is at most
            lbfgs_parameter_t::xtol. */
        case LBFGSERR_WIDTHTOOSMALL:
            p = "WIDTHTOOSMALL";
            break;
        /** A logic error (negative line-search step) occurred. */
        case LBFGSERR_INVALIDPARAMETERS:
            p = "INVALIDPARAMETERS";
            break;
        /** The current search direction increases the objective function value. */
        case LBFGSERR_INCREASEGRADIENT:
            p = "INCREASEGRADIENT";
            break;
        case 0:
            p = "Minimization success";
            break;
        default:
            p = "No detected error";
            break;
    }
    return p;
}