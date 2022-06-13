//
// Created by Lood van Niekerk on 2022/06/13.
//

#ifndef PLMC_WEIGHTS_H
#define PLMC_WEIGHTS_H

#include "plm.h"

/* Reweights sequences by their inverse neighborhood size */
void MSAReweightSequences(alignment_t *ali, options_t *options);

/* Reads custom weights from a file
   The file should contain a list of floats, one per sequence, separated by lines.
 */
void ReadCustomWeightsFile(alignment_t *ali, options_t *options, char *weights_file);

int ValidateCustomWeightsFile(alignment_t *ali, options_t *options, char *weights_file);

#endif //PLMC_WEIGHTS_H
