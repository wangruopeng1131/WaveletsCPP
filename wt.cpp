//
// Created by wangr on 2023/10/8.
//


#include "wt.h"

/* Decomposition of input with lowpass filter */
int downcoefAxis(const double *input, ArrayInfo input_info,
                 const double *output, ArrayInfo output_info,
                 const DiscreteWavelet *wavelet, size_t axis,
                 Coefficient coef, MODE dwt_mode,
                 size_t swt_level,
                 DiscreteTransformType transform)
{
    size_t i;
    size_t num_loops = 1;
    double *temp_input = nullptr, *temp_output = nullptr;

    // These are boolean values, but MSVC does not have <stdbool.h>
    int make_temp_input, make_temp_output;

    if (input_info.ndim != output_info.ndim)
        return 1;
    if (axis >= input_info.ndim)
        return 2;

    for (i = 0; i < input_info.ndim; ++i)
    {
        if (i == axis)
        {
            switch (transform)
            {
                case DWT_TRANSFORM:
                    if (dwt_buffer_length(input_info.shape[i], wavelet->dec_len,
                                          dwt_mode) != output_info.shape[i])
                        return 3;
                    break;
                case SWT_TRANSFORM:
                    if (swt_buffer_length(input_info.shape[i])
                        != output_info.shape[i])
                        return 4;
                    break;
            }
        }
        else
        {
            if (input_info.shape[i] != output_info.shape[i])
                return 5;
        }
    }

    make_temp_input = input_info.strides[axis] != sizeof(double);
    make_temp_output = output_info.strides[axis] != sizeof(double);
    if (make_temp_input)
        if ((malloc(input_info.shape[axis] * sizeof(double))) == NULL)
            goto cleanup;
    if (make_temp_output)
        if ((malloc(output_info.shape[axis] * sizeof(double))) == NULL)
            goto cleanup;

    for (i = 0; i < output_info.ndim; ++i)
    {
        if (i != axis)
            num_loops *= output_info.shape[i];
    }

    for (i = 0; i < num_loops; ++i)
    {
        size_t j;
        size_t input_offset = 0, output_offset = 0;
        const double *input_row;
        double *output_row;

        // Calculate offset into linear buffer
        {
            size_t reduced_idx = i;
            for (j = 0; j < output_info.ndim; ++j)
            {
                size_t j_rev = output_info.ndim - 1 - j;
                if (j_rev != axis)
                {
                    size_t axis_idx = reduced_idx % output_info.shape[j_rev];
                    reduced_idx /= output_info.shape[j_rev];

                    input_offset += (axis_idx * input_info.strides[j_rev]);
                    output_offset += (axis_idx * output_info.strides[j_rev]);
                }
            }
        }

        // Copy to temporary input if necessary
        if (make_temp_input)
            for (j = 0; j < input_info.shape[axis]; ++j)
                // Offsets are byte offsets, to need to cast to char and back
                temp_input[j] = *(double *) (((char *) input) + input_offset
                                             + j * input_info.strides[axis]);

        // Select temporary or direct output and input
        input_row = make_temp_input ? temp_input
                                    : (const double *) ((const char *) input + input_offset);
        output_row = make_temp_output ? temp_output
                                      : (double *) ((char *) output + output_offset);


        switch (transform)
        {
            case DWT_TRANSFORM:
                // Apply along axis
                switch (coef)
                {
                    case COEF_APPROX:
                        dec_a(input_row, input_info.shape[axis],
                              wavelet,
                              output_row, output_info.shape[axis],
                              dwt_mode);
                        break;
                    case COEF_DETAIL:
                        dec_d(input_row, input_info.shape[axis],
                              wavelet,
                              output_row, output_info.shape[axis],
                              dwt_mode);
                        break;
                }
                break;

            case SWT_TRANSFORM:
                // Apply along axis
                switch (coef)
                {
                    case COEF_APPROX:
                        swt_a(input_row, input_info.shape[axis],
                              wavelet,
                              output_row, output_info.shape[axis],
                              swt_level);
                        break;
                    case COEF_DETAIL:
                        swt_d(input_row, input_info.shape[axis],
                              wavelet,
                              output_row, output_info.shape[axis],
                              swt_level);
                        break;
                }
                break;
        }

        // Copy from temporary output if necessary
        if (make_temp_output)
            for (j = 0; j < output_info.shape[axis]; ++j)
                // Offsets are byte offsets, to need to cast to char and back
                *(double *) ((char *) output + output_offset
                             + j * output_info.strides[axis]) = output_row[j];
    }

    free(temp_input);
    free(temp_output);
    return 0;

    cleanup:
    free(temp_input);
    free(temp_output);
    return 6;
}


int idwtAxis(const double *coefs_a, const ArrayInfo *a_info,
             const double *coefs_d, const ArrayInfo *d_info,
             double *output, ArrayInfo output_info,
             const DiscreteWavelet *wavelet,
             size_t axis, MODE mode)
{
    size_t i;
    size_t num_loops = 1;
    double *temp_coefs_a = NULL, *temp_coefs_d = NULL, *temp_output = NULL;

    // These are boolean values, but MSVC does not have <stdbool.h>
    int make_temp_coefs_a, make_temp_coefs_d, make_temp_output;
    int have_a = ((coefs_a != NULL) && (a_info != NULL));
    int have_d = ((coefs_d != NULL) && (d_info != NULL));


    if (!have_a && !have_d)
        return 3;

    if ((have_a && (a_info->ndim != output_info.ndim)) ||
        (have_d && (d_info->ndim != output_info.ndim)))
        return 1;
    if (axis >= output_info.ndim)
        return 1;

    for (i = 0; i < output_info.ndim; ++i)
    {
        if (i == axis)
        {
            size_t input_shape;
            if (have_a && have_d &&
                (d_info->shape[i] != a_info->shape[i]))
                return 1;
            input_shape = have_a ? a_info->shape[i] : d_info->shape[i];

            /* TODO: reconstruction_buffer_length should take a & d shapes
             *       - for odd output_len, d_len == (a_len - 1)
             */
            if (idwt_buffer_length(input_shape, wavelet->rec_len, mode)
                != output_info.shape[i])
                return 1;
        }
        else
        {
            if ((have_a && (a_info->shape[i] != output_info.shape[i])) ||
                (have_d && (d_info->shape[i] != output_info.shape[i])))
                return 1;
        }
    }

    make_temp_coefs_a = have_a && a_info->strides[axis] != sizeof(double);
    make_temp_coefs_d = have_d && d_info->strides[axis] != sizeof(double);
    make_temp_output = output_info.strides[axis] != sizeof(double);
    if (make_temp_coefs_a)
        if ((malloc(a_info->shape[axis] * sizeof(double))) == NULL)
            goto cleanup;
    if (make_temp_coefs_d)
        if ((malloc(d_info->shape[axis] * sizeof(double))) == NULL)
            goto cleanup;
    if (make_temp_output)
        if ((malloc(output_info.shape[axis] * sizeof(double))) == NULL)
            goto cleanup;

    for (i = 0; i < output_info.ndim; ++i)
    {
        if (i != axis)
            num_loops *= output_info.shape[i];
    }

    for (i = 0; i < num_loops; ++i)
    {
        size_t j;
        size_t a_offset = 0, d_offset = 0, output_offset = 0;
        double *output_row;

        // Calculate offset into linear buffer
        {
            size_t reduced_idx = i;
            for (j = 0; j < output_info.ndim; ++j)
            {
                size_t j_rev = output_info.ndim - 1 - j;
                if (j_rev != axis)
                {
                    size_t axis_idx = reduced_idx % output_info.shape[j_rev];
                    reduced_idx /= output_info.shape[j_rev];

                    if (have_a)
                        a_offset += (axis_idx * a_info->strides[j_rev]);
                    if (have_d)
                        d_offset += (axis_idx * d_info->strides[j_rev]);
                    output_offset += (axis_idx * output_info.strides[j_rev]);
                }
            }
        }

        // Copy to temporary input if necessary
        if (make_temp_coefs_a)
            for (j = 0; j < a_info->shape[axis]; ++j)
                // Offsets are byte offsets, to need to cast to char and back
                temp_coefs_a[j] = *(double *) ((char *) coefs_a + a_offset
                                               + j * a_info->strides[axis]);
        if (make_temp_coefs_d)
            for (j = 0; j < d_info->shape[axis]; ++j)
                // Offsets are byte offsets, to need to cast to char and back
                temp_coefs_d[j] = *(double *) ((char *) coefs_d + d_offset
                                               + j * d_info->strides[axis]);

        // Select temporary or direct output
        output_row = make_temp_output ? temp_output
                                      : (double *) ((char *) output + output_offset);

        // upsampling_convolution adds to input, so zero
        memset(output_row, 0, output_info.shape[axis] * sizeof(double));

        if (have_a)
        {
            // Pointer arithmetic on NULL is undefined
            const double *a_row = make_temp_coefs_a ? temp_coefs_a
                                                    : (const double *) ((const char *) coefs_a + a_offset);
            upsamplingConvolutionValidSf(a_row, a_info->shape[axis],
                                         wavelet->rec_lo_double, wavelet->rec_len,
                                         output_row, output_info.shape[axis],
                                         mode);
        }
        if (have_d)
        {
            // Pointer arithmetic on NULL is undefined
            const double *d_row = make_temp_coefs_d ? temp_coefs_d
                                                    : (const double *) ((const char *) coefs_d + d_offset);
            upsamplingConvolutionValidSf(d_row, d_info->shape[axis],
                                         wavelet->rec_hi_double, wavelet->rec_len,
                                         output_row, output_info.shape[axis],
                                         mode);
        }

        // Copy from temporary output if necessary
        if (make_temp_output)
            for (j = 0; j < output_info.shape[axis]; ++j)
                // Offsets are byte offsets, to need to cast to char and back
                *(double *) ((char *) output + output_offset
                             + j * output_info.strides[axis]) = output_row[j];
    }

    free(temp_coefs_a);
    free(temp_coefs_d);
    free(temp_output);
    return 0;

    cleanup:
    free(temp_coefs_a);
    free(temp_coefs_d);
    free(temp_output);
    return 2;
}

int dec_a(const double *input, size_t input_len,
          const DiscreteWavelet *wavelet,
          double *output, size_t output_len,
          MODE mode)
{

    /* check output length */
    int l = dwt_buffer_length(input_len, wavelet->dec_len, mode);
    if (output_len != l)
    {
        return -1;
    }

    return downsamplingConvolution(input, input_len,
                                   wavelet->dec_lo_double,
                                   wavelet->dec_len, output,
                                   2, mode);
}


/* Decomposition of input with highpass filter */

int dec_d(const double *input, size_t input_len,
          const DiscreteWavelet *wavelet,
          double *output, size_t output_len,
          MODE mode)
{

    /* check output length */
    if (output_len != dwt_buffer_length(input_len, wavelet->dec_len, mode))
        return -1;

    return downsamplingConvolution(input, input_len,
                                   wavelet->dec_hi_double,
                                   wavelet->dec_len, output,
                                   2, mode);
}


/* Direct reconstruction with lowpass reconstruction filter */

int rec_a(const double *coeffs_a, size_t coeffs_len,
          const DiscreteWavelet *wavelet,
          double *output, size_t output_len)
{

    /* check output length */
    if (output_len != reconstruction_buffer_length(coeffs_len, wavelet->rec_len))
        return -1;

    return upsamplingConvolutionFull(coeffs_a, coeffs_len,
                                     wavelet->rec_lo_double,
                                     wavelet->rec_len, output,
                                     output_len);
}


/* Direct reconstruction with highpass reconstruction filter */

int rec_d(const double *coeffs_d, size_t coeffs_len,
          const DiscreteWavelet *wavelet,
          double *output, size_t output_len)
{

    /* check for output length */
    if (output_len != reconstruction_buffer_length(coeffs_len, wavelet->rec_len))
        return -1;

    return upsamplingConvolutionFull(coeffs_d, coeffs_len,
                                     wavelet->rec_hi_double,
                                     wavelet->rec_len, output,
                                     output_len);
}


/*
 * IDWT reconstruction from approximation and detail coeffs, either of which may
 * be NULL.
 *
 * Requires zero-filled output buffer.
 */
int idwt(const double *coeffs_a, size_t coeffs_a_len,
         const double *coeffs_d, size_t coeffs_d_len,
         double *output, size_t output_len,
         const DiscreteWavelet *wavelet, MODE mode)
{
    size_t input_len;
    if (coeffs_a != NULL && coeffs_d != NULL)
    {
        if (coeffs_a_len != coeffs_d_len)
            goto error;
        input_len = coeffs_a_len;
    }
    else if (coeffs_a != NULL)
    {
        input_len = coeffs_a_len;
    }
    else if (coeffs_d != NULL)
    {
        input_len = coeffs_d_len;
    }
    else
    {
        goto error;
    }

    /* check output size */
    if (output_len != idwt_buffer_length(input_len, wavelet->rec_len, mode))
        goto error;

    /*
     * Set output to zero (this can be omitted if output array is already
     * cleared) memset(output, 0, output_len * sizeof(T));
     */

    /* reconstruct approximation coeffs with lowpass reconstruction filter */
    if (coeffs_a)
    {
        if (upsamplingConvolutionValidSf(coeffs_a, input_len,
                                         wavelet->rec_lo_double,
                                         wavelet->rec_len, output,
                                         output_len, mode) < 0)
        {
            goto error;
        }
    }
    /*
     * Add reconstruction of details coeffs performed with highpass
     * reconstruction filter.
     */
    if (coeffs_d)
    {
        if (upsamplingConvolutionValidSf(coeffs_d, input_len,
                                         wavelet->rec_hi_double,
                                         wavelet->rec_len, output,
                                         output_len, mode) < 0)
        {
            goto error;
        }
    }

    return 0;

    error:
    return -1;
}

/* basic SWT step (TODO: optimize) */
int swt(const double *input, pywt_index_t input_len,
        const double *filter, pywt_index_t filter_len,
        double *output, size_t output_len,
        unsigned int level)
{

    double *e_filter;
    pywt_index_t i, e_filter_len, fstep;
    int ret;

    if (level < 1)
        return -1;

    if (level > swt_max_level(input_len))
        return -2;

    if (output_len != swt_buffer_length(input_len))
        return -1;

    /* TODO: quick hack, optimize */
    if (level > 1)
    {
        /* allocate filter first */
        e_filter_len = filter_len << (level - 1);
        if ((wtcalloc(e_filter_len, sizeof(double))) == NULL)
            goto cleanup;
        if (e_filter == NULL)
            return -1;
        fstep = 1 << (level - 1);  // spacing between non-zero filter entries

        /* compute upsampled filter values */
        for (i = 0; i < filter_len; ++i)
        {
            e_filter[i << (level - 1)] = filter[i];
        }
        ret = downsamplingConvolutionPeriodization(input, input_len, e_filter,
                                                   e_filter_len, output, 1,
                                                   fstep);
        wtfree(e_filter);
        return ret;

    }
    else
    {
        return downsamplingConvolutionPeriodization(input, input_len, filter,
                                                    filter_len, output, 1,
                                                    1);
    }
    cleanup:
    wtfree(e_filter);
    return -3;
}

/*
 * Approximation at specified level
 * input - approximation coeffs from upper level or signal if level == 1
 */
int swt_a(const double *input, pywt_index_t input_len,
          const DiscreteWavelet *wavelet,
          double *output, pywt_index_t output_len,
          unsigned int level)
{
    return swt(input, input_len, wavelet->dec_lo_double,
               wavelet->dec_len, output, output_len, level);
}

/* Details at specified level
 * input - approximation coeffs from upper level or signal if level == 1
 */
int swt_d(const double *input, pywt_index_t input_len,
          const DiscreteWavelet *wavelet,
          double *output, pywt_index_t output_len,
          unsigned int level)
{
    return swt(input, input_len, wavelet->dec_hi_double,
               wavelet->dec_len, output, output_len, level);
}


