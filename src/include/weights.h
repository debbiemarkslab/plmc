#ifndef PLMC_WEIGHTS_H
#define PLMC_WEIGHTS_H

#include "plm.h"

/* Reweights sequences by their inverse neighborhood size */
void MSAReweightSequences(alignment_t *ali, options_t *options);

/* Reads custom weights from a file
   The file should contain a list of floats, one per sequence, separated by lines.
 */
void ReadCustomWeightsFile(char *weightsFile, alignment_t *ali);

int ValidateCustomWeightsFile(char *weightsFile, alignment_t *ali);

void WriteWeightsFile(char *weightsOutputFile, alignment_t *ali);

#endif //PLMC_WEIGHTS_H
