# plmc
plmc infers pairwise undirected graphical models for families of biological sequences. It implements a penalized, maximum-pseudolikelihood approach and supports custom alphabets as well as various forms of regularization. With a multiple sequence alignment as an input, plmc can quantify inferred coupling strengths between all pairs of positions (couplingsfile output) or infer a generative model of the sequences for predicting the effects of mutations (paramfile output).

## Usage
      plmc [options] alignmentfile
      plmc -c couplingsfile alignmentfile
      plmc -o paramfile -c couplingsfile alignmentfile
      plmc [-h | --help]
      
    Required input:
      alignmentfile                    Multiple sequence alignment in FASTA format

    Options, output:
      -c  --couplings  couplingsfile   Save coupling scores to file (text)
      -o  --output     paramfile       Save estimated parameters to file (binary)

    Options, alignment processing:
      -s  --scale      <value>         Sequence weights: neighborhood weight [s > 0]
      -t  --theta      <value>         Sequence weights: neighborhood divergence [0 < t < 1]

    Options, Maximum a posteriori estimation (L-BFGS, default):
      -lh --lambdah    <value>         Set L2 lambda for fields (h_i)
      -le --lambdae    <value>         Set L2 lambda for couplings (e_ij)
      -lg --lambdag    <value>         Set group L1 lambda for couplings (e_ij)

    Options, general:
      -a  --alphabet   alphabet        Alternative alphabet (default 21: -ACDEFGHIKLMNPQRSTVWY)
      -f  --focus      identifier      Select only uppercase, non-gapped sites from a focus sequence
      -g  --gapignore                  Exclude first alphabet character from potential calculations
      -m  --maxiter                    Maximum number of iterations
      -n  --ncores    [<number>|max]   Maximum number of threads to use in OpenMP
      -h  --help                       Usage

## Compilation
plmc can be compiled to single-core without external libraries, but requires OpenMP for multi-core compilation. On Mac OS X, `clang` does not yet support OpenMP, but versions of `GCC` that do can be found as precompiled binaries [here](http://hpc.sourceforge.net/) or downloaded through package managers like homebrew or macports.

**Multicore**. To compile with `gcc` and OpenMP: 

    make all-openmp

**Single core, Linux**. To compile with `gcc`: 

    make all

**Single core, Mac OS X**. To compile with `clang`:

    make all-mac

**Single precision**. All of the above targets compile to double precision (64 bit), but reducing the precision to single (32 bit) increases speed and decreases memory requirements by approximately a factor of two. The fastest compile settings are then:

    make all-openmp32

## Examples
**Protein alignments**. The example directory includes an alignment of the protein [dihdyrofolate reductase](https://en.wikipedia.org/wiki/Dihydrofolate_reductase) (DHFR). To infer a model for this family, we can type the following in the base directory:

    bin/plmc -o example/DHFR/DHFR.eij -le 16.0 -lh 0.01 -m 100 -g -f DYR_ECOLI example/DHFR/DHFR.a2m
The numeric options set a strong L2 regularization for the couplings, λ<sub>e</sub> = 16.0, a weak L2 regularization for the sites, λ<sub>h</sub> = 0.01, and the maximum number of iterations at 100. The focus `-f` option tells plmc to only model columns that are present in the E. coli sequence DYR_ECOLI, and the  `-g` gap-ignoring option ignores gaps by modeling only the coding portions of each sequence. To read the binary paramfile `DHFR.eij` and visualize the couplings, we can type the following in MATLAB from the `scripts` directory:

    plot_corrected_couplings('../example/DHFR/DHFR.eij')

This computes and plots the same (APC-corrected) coupling strengths that would be output to couplingfile, revealing many strongly coupled pairs of positions in the long-term evolution of DHFR:
<p align="center"><img src="example/DHFR/DHFR.png" width="400"></p>

**RNA alignments**. To override the default alphabet ("-ACDEFGHIKLMNPQRSTVWY" for proteins) use the option `-a ALPHABET`. An example RNA alignment is included for the [SAM riboswitch](https://en.wikipedia.org/wiki/SAM_riboswitch_(S_box_leader)). To infer the couplings with an RNA alphabet (".ACGU") type the following in the base directory:

    bin/plmc -c example/RNA/RF00162.EC -o example/RNA/RF00162.eij -a .ACGU -le 20.0 -lh 0.01 -m 50 -g example/RNA/RF00162.fasta

To plot we type the following in MATLAB from the `scripts` directory:

    plot_corrected_couplings('../example/RNA/RF00162.eij')
<p align="center"><img src="example/RNA/RF00162.png" width="500"></p>

**Physical systems**. As an example, simulated draws from a 3-state, 1-dimensional Potts model are provided in the examples folder and encoded by the characters _, *, and ^. The following command would estimate the parameters by running to convergence with λ<sub>e</sub> = 1.0, λ<sub>h</sub> = 1.0 and sequence reweighting disabled:

    bin/plmc -c example/potts/potts3.txt -a _*^ -t -1 -le 1.0 -lh 1.0 example/potts/potts3.a2m
A 1D Potts model will only have interactions between i -> i + 1, which should be evident in the coupling summary scores output to example/potts/potts3.txt 

## Author
plmc was written by [John Ingraham](mailto:john.ingraham@gmail.com) in [Debora Marks' lab](https://marks.hms.harvard.edu/) at Harvard Medical School

## Credits
The MAP-based inference uses a [C implementation of L-BFGS by Naoaki Okazaki](https://github.com/chokkan/liblbfgs "libLBFGS"), which is included in this repo


