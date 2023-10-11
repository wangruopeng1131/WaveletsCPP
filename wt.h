#ifndef WT_WT_H
#define WT_WT_H

#include "convolution.h"
#include "wavelets.h"
#include "templating.h"

/* Decomposition of input with lowpass filter */

int downcoefAxis(const double *input, ArrayInfo input_info,
                 const double *output, ArrayInfo output_info,
                 const DiscreteWavelet *wavelet, size_t axis,
                 Coefficient coef, MODE dwt_mode,
                 size_t swt_level,
                 DiscreteTransformType transform);


int idwtAxis(const double *coefs_a, const ArrayInfo *a_info,
             const double *coefs_d, const ArrayInfo *d_info,
             double *output, ArrayInfo output_info,
             const DiscreteWavelet *wavelet,
             size_t axis, MODE mode);


int dec_a(const double *input, size_t input_len,
          const DiscreteWavelet *wavelet,
          double *output, size_t output_len,
          MODE mode);


/* Decomposition of input with highpass filter */

int dec_d(const double *input, size_t input_len,
          const DiscreteWavelet *wavelet,
          double *output, size_t output_len,
          MODE mode);


/* Direct reconstruction with lowpass reconstruction filter */

int rec_a(const double * coeffs_a, size_t coeffs_len,
          const DiscreteWavelet * wavelet,
          double *output, size_t output_len);


/* Direct reconstruction with highpass reconstruction filter */

int rec_d(const double *coeffs_d, size_t coeffs_len,
          const DiscreteWavelet *wavelet,
          double *output, size_t output_len);


/*
 * IDWT reconstruction from approximation and detail coeffs, either of which may
 * be NULL.
 *
 * Requires zero-filled output buffer.
 */
int idwt(const double *coeffs_a, size_t coeffs_a_len,
         const double *coeffs_d, size_t coeffs_d_len,
         double *output, size_t output_len,
         const DiscreteWavelet *wavelet, MODE mode);

/* basic SWT step (TODO: optimize) */
int swt(const double *input, pywt_index_t input_len,
        const double *filter, pywt_index_t filter_len,
        double *output, size_t output_len,
        unsigned int level);

/*
 * Approximation at specified level
 * input - approximation coeffs from upper level or signal if level == 1
 */
int swt_a(const double *input, pywt_index_t input_len,
          const DiscreteWavelet *wavelet,
          double *output, pywt_index_t output_len,
          unsigned int level);

/* Details at specified level
 * input - approximation coeffs from upper level or signal if level == 1
 */
int swt_d(const double *input, pywt_index_t input_len,
          const DiscreteWavelet *wavelet,
          double *output, pywt_index_t output_len,
          unsigned int level);

#endif //WT_WT_H

