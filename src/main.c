#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include "include/plm.h"

/* Optionally include OpenMP with the -fopenmp flag */
#if defined(_OPENMP)
    #include <omp.h>
#endif

/* Usage pattern */
const char *usage =
"plmc\n"
"\n"
"Usage:\n"
"      plm [options] alignmentfile\n"
"      plm -c couplingsfile alignmentfile\n"
"      plm -o paramfile -c couplingsfile alignmentfile\n"
"      plm [-h | --help]\n"
"      \n"
"    Required input:\n"
"      alignmentfile                    Multiple sequence alignment in FASTA format\n"
"\n"
"    Options, input:\n"
"      -w  --weights    weightsfile     Load sequence weights from file (one weight per line)\n"
"\n"
"    Options, output:\n"
"      -c  --couplings  couplingsfile   Save coupling scores to file (text)\n"
"      -o  --output     paramfile       Save estimated parameters to file (binary)\n"
"      --save-weights   weightsfile     Save sequence weights to file (text)\n"
"\n"
"    Options, alignment processing:\n"
"      -s  --scale      <value>         Sequence weights: neighborhood weight [s > 0]\n"
"      -t  --theta      <value>         Sequence weights: neighborhood divergence [0 < t < 1]\n"
"\n"
"    Options, Maximum a posteriori estimation (L-BFGS, default):\n"
"      -lh --lambdah    <value>         Set L2 lambda for fields (h_i)\n"
"      -le --lambdae    <value>         Set L2 lambda for couplings (e_ij)\n"
"      -lg --lambdag    <value>         Set group L1 lambda for couplings (e_ij)\n"
"\n"
"    Options, general:\n"
"          --fast                       Fast weights and stochastic gradient descent\n"
"      -a  --alphabet   alphabet        Alternative character set to use for analysis\n"
"      -f  --focus      identifier      Select only uppercase, non-gapped sites from a focus sequence\n"
"      -g  --gapignore                  Model sequence likelihoods only by coding, non-gapped portions\n"
"      -m  --maxiter                    Maximum number of iterations\n"
"      -n  --ncores    [<number>|max]   Maximum number of threads to use in OpenMP\n"
"      -h  --help                       Usage\n\n";

int main(int argc, char **argv) {
    char *alignFile = NULL;
    char *outputFile = NULL;
    char *couplingsFile = NULL;
    char *weightsFile = NULL;
    char *weightsOutputFile = NULL;
    options_t *options = default_options();

    /* Print usage if no arguments */
    if (argc == 1) {
        fprintf(stderr, "%s", usage);
        exit(1);
    }

    /* Parse command line arguments */
    for (int arg = 1; arg < argc; arg++) {
        if ((arg < argc-1) && (strcmp(argv[arg], "--output") == 0
                    || strcmp(argv[arg], "-o") == 0)) {
            outputFile = argv[++arg];
        } else if ((arg < argc-1) && (strcmp(argv[arg], "--alphabet") == 0
                    || strcmp(argv[arg], "-a") == 0)) {
            options->alphabet = argv[++arg];
        } else if ((arg < argc-1) && (strcmp(argv[arg], "--couplings") == 0
                    || strcmp(argv[arg], "-c") == 0)) {
            couplingsFile = argv[++arg];
        } else if ((arg < argc-1) && (strcmp(argv[arg], "--lambdah") == 0
                    || strcmp(argv[arg], "-lh") == 0)) {
            options->lambdaH = atof(argv[++arg]);
        } else if ((arg < argc-1) && (strcmp(argv[arg], "--lambdae") == 0
                    || strcmp(argv[arg], "-le") == 0)) {
            options->lambdaE = atof(argv[++arg]);
        } else if ((arg < argc-1) && (strcmp(argv[arg], "--lambdag") == 0
                    || strcmp(argv[arg], "-lg") == 0)) {
            options->lambdaGroup = atof(argv[++arg]);
        } else if ((arg < argc-1) && (strcmp(argv[arg], "--theta") == 0
                    || strcmp(argv[arg], "-t") == 0)) {
            options->theta = atof(argv[++arg]);
        } else if ((arg < argc-1) && (strcmp(argv[arg], "--scale") == 0
                    || strcmp(argv[arg], "-s") == 0)) {
            options->scale = atof(argv[++arg]);
        } else if ((arg < argc-1)  && (strcmp(argv[arg], "--maxiter") == 0
                    || strcmp(argv[arg], "-m") == 0)) {
            options->maxIter = atoi(argv[++arg]);
        } else if ((arg < argc-1) && (strcmp(argv[arg], "--independent") == 0
                    || strcmp(argv[arg], "-i") == 0)) {
            options->usePairs = 0;
            fprintf(stderr, "Independent model not yet implemented\n");
            exit(0);
        } else if ((arg < argc-1) && (strcmp(argv[arg], "--gapreduce") == 0
                    || strcmp(argv[arg], "-g") == 0)) {
            options->estimatorMAP = INFER_MAP_PLM_GAPREDUCE;
        } else if ((arg < argc-1) && (strcmp(argv[arg], "--estimatele") == 0
                    || strcmp(argv[arg], "-ee") == 0)) {
            options->zeroAPC = 1;
        } else if ((arg < argc-1) && (strcmp(argv[arg], "--focus") == 0
                    || strcmp(argv[arg], "-f") == 0)) {
            options->target = argv[++arg];
        } else if ((arg < argc-1) && strcmp(argv[arg], "--fast") == 0) {
            options->sgd = 1;
            options->fastWeights = 100;
        } else if ((arg < argc-1) && (strcmp(argv[arg], "--weights") == 0
                    || strcmp(argv[arg], "-w") == 0)) {
            weightsFile = argv[++arg];
        } else if ((arg < argc-1) && (strcmp(argv[arg], "--save-weights") == 0)) {
            weightsOutputFile = argv[++arg];
        } else if ((arg < argc-1) && (strcmp(argv[arg], "--ncores") == 0
                    || strcmp(argv[arg], "-n") == 0)) {
            #if defined(_OPENMP)
                if (strcmp(argv[arg + 1], "max") == 0) {
                    int maxThreads = omp_get_max_threads();
                    /* Redundant, but serves as sanity check */
                    omp_set_num_threads(maxThreads);
                    fprintf(stderr, "OpenMP: Using %d of %d threads\n",
                        maxThreads, maxThreads);
                } else {
                    int numThreads = atoi(argv[arg + 1]);
                    int maxThreads = omp_get_max_threads();
                    if (numThreads >= 1 && numThreads <= maxThreads) {
                        omp_set_num_threads(numThreads);
                        fprintf(stderr, "OpenMP: Using %d of %d threads\n",
                            numThreads, maxThreads);
                    } else if (numThreads > maxThreads) {
                        omp_set_num_threads(maxThreads);
                        fprintf(stderr, "OpenMP: More threads requested than "
                            "available. Using %d of %d threads instead.\n",
                            maxThreads, maxThreads);
                    } else {
                        omp_set_num_threads(1);
                        fprintf(stderr, "OpenMP: Using 1 of %d threads\n",
                            maxThreads);
                    }
                }
                arg++;
            #else
                fprintf(stderr, "Error (-n/--ncores) only available when "
                    "compiled with OpenMP\n");
                exit(1);
            #endif
        } else if (strcmp(argv[arg], "--help") == 0
                    || strcmp(argv[arg], "-h") == 0) {
            fprintf(stderr, "%s", usage);
            exit(1);
        }
    }
    alignFile = argv[argc - 1];

    run_plmc(alignFile, outputFile, couplingsFile, weightsFile, weightsOutputFile, options);
}