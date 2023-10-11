//
// Created by 83793 on 2022/11/12.
//

#ifndef WT_CONVOLUTION_H
#define WT_CONVOLUTION_H

#include "common.h"

/* This file contains several functions for computing the convolution of a
 * signal with a filter. The general scheme is:
 *   output[o] = sum(filter[j] * input[i-j] for j = [0..F) and i = [0..N))
 * where 'o', 'i' and 'j' may progress at different rates.
 *
 * Most of the code deals with different edge extension modes. Values are
 * computed on-demand, in four steps:
 * 1. Filter extends past signal on the left.
 * 2. Filter completely contained within signal (no extension).
 * 3. Filter extends past signal on both sides (only if F > N).
 * 4. Filter extends past signal on the right.
 *
 * MODE_PERIODIZATION produces different output lengths to other modes, so is
 * implemented as a separate function for each case.
 *
 * See 'common.h' for descriptions of the extension modes.
 */
int downsamplingConvolutionPeriodization(const double *const input, const size_t N,
                                                const double *const filter, const size_t F,
                                                double *const output, const size_t step,
                                                const size_t fstep);


int downsamplingConvolution(const double *const input, const size_t N,
                            const double *const filter, const size_t F,
                            double *const output,
                            const size_t step, MODE mode);

int upsamplingConvolutionFull(const double *input, const size_t N,
                              const double *filter, const size_t F,
                              double *const output, const size_t O);

int upsamplingConvolutionValidSfPeriodization(const double *const input, const size_t N,
                                                     const double *const filter, const size_t F,
                                                     double *const output, const size_t O);


/*
 * performs IDWT for all modes
 *
 * The upsampling is performed by splitting filters to even and odd elements
 * and performing 2 convolutions.  After refactoring the PERIODIZATION mode
 * case to separate function this looks much clearer now.
 */

int upsamplingConvolutionValidSf(const double *const input, const size_t N,
                                 const double *const filter, const size_t F,
                                 double *const output, const size_t O,
                                 MODE mode) ;

/* -> swt - todo */
int upsampledFilterConvolution(const double *const input, const size_t N,
                               const double *const filter, const size_t F,
                               double *const output,
                               const size_t step, MODE mode);

#endif //WT_CONVOLUTION_H
