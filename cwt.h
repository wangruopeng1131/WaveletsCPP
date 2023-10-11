/* Copyright (c) 2016 Holger Nahrstaedt */
/* See COPYING for license details. */


#include "convolution.h"
#define _USE_MATH_DEFINES

#include "math.h"
template<typename T>
T pow(const T x, const T y) {
    if (sizeof(T) == sizeof(double))
        return pow(x, y);
    else
        return powf(x, y);
}

template<typename T>
T sqrt(const T x) {
    if (sizeof(T) == sizeof(double))
        return sqrt(x);
    else
        return sqrtf(x);
}

template<typename T>
T exp(const T x) {
    if (sizeof(T) == sizeof(double))
        return exp(x);
    else
        return expf(x);
}

template<typename T>
T cos(const T x) {
    if (sizeof(T) == sizeof(double))
        return cos(x);
    else
        return cosf(x);
}

template<typename T>
T sin(const T x) {
    if (sizeof(T) == sizeof(double))
        return sin(x);
    else
        return sinf(x);
}

template<typename T>
void gaus(const T *const input,
          T *const output, const size_t N,
          const size_t number) {
    size_t i = 0;
    for (i = 0; i < N; i++) {
        switch (number) {
            case 1:
                output[i] = -2 * input[i] * exp((-pow(input[i], 2.0))) /
                            sqrt(sqrt(M_PI / 2));
                break;
            case 2:
                output[i] = -2 * (2 * pow(input[i], 2.0) - 1) * exp(-pow(input[i], 2.0)) /
                            sqrt(3 * sqrt(M_PI / 2));
                break;
            case 3:
                output[i] = -4 * -2 * pow(input[i], 3.0) + 3 * input[i] *
                                                           exp(-pow(input[i], 2.0)) /
                                                           sqrt(15 * sqrt(M_PI / 2));
                break;
            case 4:
                output[i] = 4 * (-12 * pow(input[i], 2.0) + 4 * pow(input[i], 4.0) + 3) *
                            exp(-CATpow(input[i], 2.0)) /
                            sqrt(105 * sqrt(M_PI / 2));
                break;
            case 5:
                output[i] = 8 * (-4 * pow(input[i], 5.0) + 20 * pow(input[i], 3.0) - 15 * input[i]) *
                            exp(-pow(input[i], 2.0)) /
                            sqrt(105 * 9 * sqrt(M_PI / 2));
                break;
            case 6:
                output[i] = -8 * (8 * pow(input[i], 6.0) - 60 * pow(input[i], 4.0) +
                                  90 * pow(input[i], 2.0) - 15) * exp(-pow(input[i], 2.0)) /
                            sqrt(105 * 9 * 11 * sqrt(M_PI / 2));
                break;
            case 7:
                output[i] = -16 * (-8 * pow(input[i], 7.0) + 84 * pow(input[i], 5.0) -
                                   210 * pow(input[i], 3.0) + 105 * (input[i])) *
                            exp(-pow(input[i], 2.0)) /
                            sqrt(105 * 9 * 11 * 13 * sqrt(M_PI / 2));
                break;
            case 8:
                output[i] = 16 * (16 * pow(input[i], 8.0) - 224 * pow(input[i], 6.0) +
                                  840 * pow(input[i], 4.0) - 840 * pow(input[i], 2.0) + 105) *
                            exp(-pow(input[i], 2.0)) /
                            sqrt(105 * 9 * 11 * 13 * 15 * sqrt(M_PI / 2));
                break;
        }
    }
}


template<typename T>
void mexh(const T *const input, T *const output, const size_t N) {
    size_t i = 0;
    for (i = 0; i < N; i++) {
        output[i] = (1 - pow(input[i], 2.0)) * exp(-pow(input[i], 2.0) / 2) * 2 /
                    (sqrt(3) * sqrt(sqrt(M_PI)));
    }
}

template<typename T>
void morl(const T *const input, T *const output, const size_t N) {
    size_t i = 0;
    for (i = 0; i < N; i++) {
        output[i] = cos(5 * input[i]) * exp(-pow(input[i], 2.0) / 2);
    }
}


template<typename T>
void cgau(const T *const input,
          T *const output_r, T *const output_i, const size_t N,
          const size_t number) {
    size_t i = 0;
    for (i = 0; i < N; i++) {
        switch (number) {
            case 1:
                output_r[i] = (-2 * input[i] * cos(input[i]) - sin(input[i])) *
                              exp(-pow(input[i], 2.0)) /
                              sqrt(2 * sqrt(M_PI / 2));
                output_i[i] = (2 * input[i] * sin(input[i]) - cos(input[i])) *
                              exp(-pow(input[i], 2.0)) /
                              sqrt(2 * sqrt(M_PI / 2));
                break;
            case 2:
                output_r[i] = (4 * pow(input[i], 2.0) * cos(input[i]) +
                               4 * input[i] * sin(input[i]) - 3 * cos(input[i])) *
                              exp(-pow(input[i], 2.0)) /
                              sqrt(10 * sqrt(M_PI / 2));
                output_i[i] = (-4 * pow(input[i], 2.0) * sin(input[i]) +
                               4 * input[i] * cos(input[i]) + 3 * sin(input[i])) *
                              exp(-pow(input[i], 2.0)) /
                              sqrt(10 * sqrt(M_PI / 2));
                break;
            case 3:
                output_r[i] = (-8 * pow(input[i], 3.0) * cos(input[i]) -
                               12 * pow(input[i], 2.0) * sin(input[i]) +
                               18 * input[i] * cos(input[i]) + 7 * sin(input[i])) *
                              exp(-pow(input[i], 2.0)) /
                              sqrt(76 * sqrt(M_PI / 2));
                output_i[i] = (8 * pow(input[i], 3.0) * sin(input[i]) -
                               12 * pow(input[i], 2.0) * cos(input[i]) -
                               18 * input[i] * sin(input[i]) + 7 * cos(input[i])) *
                              exp(-pow(input[i], 2.0)) /
                              sqrt(76 * sqrt(M_PI / 2));

                break;
            case 4:
                output_r[i] = (16 * pow(input[i], 4.0) * cos(input[i]) +
                               32 * pow(input[i], 3.0) * sin(input[i]) -
                               72 * pow(input[i], 2.0) * cos(input[i]) -
                               56 * input[i] * sin(input[i]) + 25 * cos(input[i])) *
                              exp(-pow(input[i], 2.0)) /
                              sqrt(764 * sqrt(M_PI / 2));;
                output_i[i] = (-16 * pow(input[i], 4.0) * sin(input[i]) +
                               32 * pow(input[i], 3.0) * cos(input[i]) +
                               72 * pow(input[i], 2.0) * sin(input[i]) -
                               56 * input[i] * cos(input[i]) - 25 * sin(input[i])) *
                              exp(-pow(input[i], 2.0)) /
                              sqrt(764 * sqrt(M_PI / 2));

                break;
            case 5:
                output_r[i] = (-32 * pow(input[i], 5.0) * cos(input[i]) -
                               80 * pow(input[i], 4.0) * sin(input[i]) +
                               240 * pow(input[i], 3.0) * cos(input[i]) +
                               280 * pow(input[i], 2.0) * sin(input[i]) -
                               250 * input[i] * cos(input[i]) - 81 * sin(input[i])) *
                              exp(-pow(input[i], 2.0)) /
                              sqrt(9496 * sqrt(M_PI / 2));
                output_i[i] = (32 * pow(input[i], 5.0) * sin(input[i]) -
                               80 * pow(input[i], 4.0) * cos(input[i]) -
                               240 * pow(input[i], 3.0) * sin(input[i]) +
                               280 * pow(input[i], 2.0) * cos(input[i]) +
                               250 * input[i] * sin(input[i]) - 81 * cos(input[i])) *
                              exp(-pow(input[i], 2.0)) /
                              sqrt(9496 * sqrt(M_PI / 2));

                break;
            case 6:
                output_r[i] = (64 * pow(input[i], 6.0) * cos(input[i]) +
                               192 * pow(input[i], 5.0) * sin(input[i]) -
                               720 * pow(input[i], 4.0) * cos(input[i]) -
                               1120 * pow(input[i], 3.0) * sin(input[i]) +
                               1500 * pow(input[i], 2.0) * cos(input[i]) +
                               972 * input[i] * sin(input[i]) - 331 * cos(input[i])) *
                              exp(-pow(input[i], 2.0)) /
                              sqrt(140152 * sqrt(M_PI / 2));
                output_i[i] = (-64 * pow(input[i], 6.0) * sin(input[i]) +
                               192 * pow(input[i], 5.0) * cos(input[i]) +
                               720 * pow(input[i], 4.0) * sin(input[i]) -
                               1120 * pow(input[i], 3.0) * cos(input[i]) -
                               1500 * pow(input[i], 2.0) * sin(input[i]) +
                               972 * input[i] * cos(input[i]) + 331 * sin(input[i])) *
                              exp(-pow(input[i], 2.0)) /
                              sqrt(140152 * sqrt(M_PI / 2));

                break;
            case 7:
                output_r[i] = (-128 * pow(input[i], 7.0) * cos(input[i]) -
                               448 * pow(input[i], 6.0) * sin(input[i]) +
                               2016 * pow(input[i], 5.0) * cos(input[i]) +
                               3920 * pow(input[i], 4.0) * sin(input[i]) -
                               7000 * pow(input[i], 3.0) * cos(input[i]) -
                               6804 * pow(input[i], 2.0) * sin(input[i]) +
                               4634 * input[i] * cos(input[i]) + 1303 * sin(input[i])) *
                              exp(-pow(input[i], 2.0)) /
                              sqrt(2390480 * sqrt(M_PI / 2));
                output_i[i] = (128 * pow(input[i], 7.0) * sin(input[i]) -
                               448 * pow(input[i], 6.0) * cos(input[i]) -
                               2016 * pow(input[i], 5.0) * sin(input[i]) +
                               3920 * pow(input[i], 4.0) * cos(input[i]) +
                               7000 * pow(input[i], 3.0) * sin(input[i]) -
                               6804 * pow(input[i], 2.0) * cos(input[i]) -
                               4634 * input[i] * sin(input[i]) + 1303 * cos(input[i])) *
                              exp(-pow(input[i], 2.0)) /
                              sqrt(2390480 * sqrt(M_PI / 2));

                break;
            case 8:
                output_r[i] = (256 * pow(input[i], 8.0) * cos(input[i]) +
                               1024 * pow(input[i], 7.0) * sin(input[i]) -
                               5376 * pow(input[i], 6.0) * cos(input[i]) -
                               12544 * pow(input[i], 5.0) * sin(input[i]) +
                               28000 * pow(input[i], 4.0) * cos(input[i]) +
                               36288 * pow(input[i], 3.0) * sin(input[i]) -
                               37072 * pow(input[i], 2.0) * cos(input[i]) -
                               20848 * input[i] * sin(input[i]) + 5937 * cos(input[i])) *
                              exp(-pow(input[i], 2.0)) /
                              sqrt(46206736 * sqrt(M_PI / 2));
                output_i[i] = (-256 * pow(input[i], 8.0) * sin(input[i]) +
                               1024 * pow(input[i], 7.0) * cos(input[i]) +
                               5376 * pow(input[i], 6.0) * sin(input[i]) -
                               12544 * pow(input[i], 5.0) * cos(input[i]) -
                               28000 * pow(input[i], 4.0) * sin(input[i]) +
                               36288 * pow(input[i], 3.0) * cos(input[i]) +
                               37072 * pow(input[i], 2.0) * sin(input[i]) -
                               20848 * input[i] * cos(input[i]) - 5937 * sin(input[i])) *
                              exp(-pow(input[i], 2.0)) /
                              sqrt(46206736 * sqrt(M_PI / 2));

                break;
        }
    }
}


template<typename T>
void shan(const T *const input, T *const output_r, T *const output_i, const size_t N,
          const T FB, const T FC) {
    size_t i = 0;
    for (i = 0; i < N; i++) {
        output_r[i] = cos(2 * M_PI * FC * input[i]) * sqrt(FB);
        output_i[i] = sin(2 * M_PI * FC * input[i]) * sqrt(FB);
        if (input[i] != 0) {
            output_r[i] *= sin(input[i] * FB * M_PI) / (input[i] * FB * M_PI);
            output_i[i] *= sin(input[i] * FB * M_PI) / (input[i] * FB * M_PI);
        }
    }
}

template<typename T>
void fbsp(const T *const input, T *const output_r, T *const output_i, const size_t N,
          const unsigned int M, const T FB, const T FC) {
    size_t i = 0;
    for (i = 0; i < N; i++) {
        if (input[i] != 0) {
            output_r[i] = cos(2 * M_PI * FC * input[i]) * sqrt(FB) * pow(
                    sin(M_PI * input[i] * FB / (T) M) / (M_PI * input[i] * FB / (T) M),
                    (T) M);
            output_i[i] = sin(2 * M_PI * FC * input[i]) * sqrt(FB) * pow(
                    sin(M_PI * input[i] * FB / (T) M) / (M_PI * input[i] * FB / (T) M),
                    (T) M);
        } else {
            output_r[i] = cos(2 * M_PI * FC * input[i]) * sqrt(FB);
            output_i[i] = sin(2 * M_PI * FC * input[i]) * sqrt(FB);
        }
    }
}


template<typename T>
void cmor(const T *const input, T *const output_r, T *const output_i, const size_t N,
          const T FB, const T FC) {
    size_t i = 0;
    for (i = 0; i < N; i++) {
        output_r[i] =
                cos(2 * M_PI * FC * input[i]) * exp(-pow(input[i], 2.0) / FB) /
                sqrt(M_PI * FB);
        output_i[i] =
                sin(2 * M_PI * FC * input[i]) * exp(-pow(input[i], 2.0) / FB) /
                sqrt(M_PI * FB);

    }
}
