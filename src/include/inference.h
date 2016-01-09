#ifndef INFERENCE_H
#define INFERENCE_H

/* Estimates parameters of maximum entropy model */
lbfgsfloatval_t *InferPairModel(alignment_t *ali, options_t *options);

#endif /* INFERENCE_H */
