#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "w97.h"

//define the CDF9/7 filter bank
#if USE_SQRT2_NORMALIZATION
static const FilterSetBior44 FilterCDF97 = {
  { -0.064538, -0.040688,  0.418091,  0.788485,  0.418091, -0.040688, -0.064538 }, //syn_low
  { -0.064538,  0.040688,  0.418091, -0.788485,  0.418091,  0.040688, -0.064538 }, //ana_high
  {  0.037827, -0.023849, -0.110624,  0.377403,  0.852699,  0.377403, -0.110624, -0.023849,  0.037827 }, //ana_low
  { -0.037827, -0.023849,  0.110624,  0.377403, -0.852699,  0.377403,  0.110624, -0.023849, -0.037827 } //syn_high
};
#else
static const FilterSetBior44 FilterCDF97 = {
  { -0.091270, -0.057542,  0.591270,  1.115086,  0.591270, -0.057542, -0.091270 },
  { -0.045635,  0.028771,  0.295635, -0.557543,  0.295635,  0.028771, -0.045635 },
  {  0.026748, -0.016864, -0.078223,  0.266864,  0.602949,  0.266864, -0.078223, -0.016864,  0.026748 },
  { -0.053496, -0.033728,  0.156446,  0.533728, -1.205898,  0.533728,  0.156446, -0.033728, -0.053496 }
};
#endif

void convWT(const double* ScalingFilter, int hLength, const double* WaveletFilter, int gLength, double* signal, int n, int stride)
{
  double* tempbank = (double*)malloc(n * sizeof(double));
  for (int i = 0; i < n; i++)
  {
    tempbank[i] = 0.0;
  }
  int filtershift_lo = (hLength-1) / 2;
  int filtershift_hi = (gLength-1) / 2 - 1;
  for (int j = 0; j < n/2; j++)
  {
    for (int k = 0; k < hLength; k++)
    {
      //symmetric padding before convolution
      int idx = 2 * j + k - filtershift_lo;
      if (idx < 0) idx = abs(idx) - 1;
      if (idx >= n) idx = 2 * n - 1 - idx;
      tempbank[j] += ScalingFilter[k] * signal[idx * stride];
    }
    for (int k = 0; k < gLength; k++)
    {
      int idx = 2 * j + k - filtershift_hi;
      if (idx < 0) idx = abs(idx) - 1;
      if (idx >= n) idx = 2 * n - 1 - idx;
      tempbank[j + n/2] += WaveletFilter[k] * signal[idx * stride];
    }
  }
  for (int i = 0; i < n; i++)
  {
    signal[i * stride] = tempbank[i];
  }
  free(tempbank);
}

void invconvWT(const double* ScalingFilter, int hLength, const double* WaveletFilter, int gLength, double* trafo, int n, int stride)
{
  double* tempbank = (double*)malloc(n * sizeof(double));
  for (int i = 0; i < n; i++)
  {
    tempbank[i] = 0.0;
  }
  int filtershift_lo = (hLength-1) / 2;
  int filtershift_hi = (gLength-1) / 2;
  for (int j = 0; j < n; j++)
  {
    for (int k = 0; k < hLength; k++)
    {
      int idx = j + k - filtershift_lo;
      if (idx % 2 == 0)
      {
        if (idx < 0) idx = abs(idx) - 1;
        if (idx >= n) idx = 2 * n - 1 - idx;
        tempbank[j] += ScalingFilter[k] * trafo[(idx/2) * stride];
      }
    }
    for (int k = 0; k < gLength; k++)
    {
      int idx = j + k - filtershift_hi;
      if (idx % 2 == 1)
      {
        if (idx < 0) idx = abs(idx) - 1;
        if (idx >= n) idx = 2 * n - 1 - idx;
        tempbank[j] += WaveletFilter[k] * trafo[(idx/2 + n/2) * stride];
      }
    }
  }
  for (int i = 0; i < n; i++)
  {
    trafo[i * stride] = tempbank[i];
  }
  free(tempbank);
}

void convWT_2d(double* x, int width, int height, int stride)
{
  const double* ScalingFilter = FilterCDF97.h_ana;
  int hLength = 9;
  const double* WaveletFilter = FilterCDF97.g_ana;
  int gLength = 7;
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    convWT(ScalingFilter, hLength, WaveletFilter, gLength, x + rowidx * stride, width, 1);
  }
  for (int colidx = 0; colidx < width; colidx++)
  {
    convWT(ScalingFilter, hLength, WaveletFilter, gLength, x + colidx, height, stride);
  }
}

void invconvWT_2d(double* x, int width, int height, int stride)
{
  const double* ScalingFilter = FilterCDF97.h_syn;
  int hLength = 7;
  const double* WaveletFilter = FilterCDF97.g_syn;
  int gLength = 9;
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    invconvWT(ScalingFilter, hLength, WaveletFilter, gLength, x + rowidx * stride, width, 1);
  }
  for (int colidx = 0; colidx < width; colidx++)
  {
    invconvWT(ScalingFilter, hLength, WaveletFilter, gLength, x + colidx, height, stride);
  }
}

void lwt97(double* x, int n, int stride)
{
  //predict 1
  double a = -1.586134342;
  for (int i = stride; i < (n - 2) * stride; i += 2*stride)
  {
    x[i] += a * (x[i - stride] + x[i + stride]);
  }
  x[(n - 1) * stride] += 2 * a * x[(n - 2) * stride];
  
  //update 1
  a = -0.05298011854;
  for (int i = 2 * stride; i < (n - 1) * stride; i += 2*stride)
  {
    x[i] += a * (x[i - stride] + x[i + stride]);
  }
  x[0] += 2 * a * x[stride];
  
  //predict 2
  a = 0.8829110762;
  for (int i = stride; i < (n - 2) * stride; i += 2*stride)
  {
    x[i] += a * (x[i - stride] + x[i + stride]);
  }
  x[(n - 1) * stride] += 2 * a * x[(n - 2) * stride];
    
  //update 2
  a = 0.4435068522;
  for (int i = 2 * stride; i < (n - 1) * stride; i += 2*stride)
  {
    x[i] += a * (x[i - stride] + x[i + stride]);
  }
  x[0] += 2 * a * x[stride];
  
  //scale
  a = 1/1.149604398;
  for (int i = 0; i < n ; i++)
  {
    if (i % 2 == 0) x[i * stride] /= a;
    if (i % 2 == 1) x[i * stride] *= -a;
#if !USE_SQRT2_NORMALIZATION
    x[i * stride] /= sqrt(2);
#endif
  }
  
  //pack
  double* tempbank = (double*)malloc(n * sizeof(double));
  for (int i = 0; i < n; i++)
  {
    if (i % 2 == 0) tempbank[i / 2] = x[i * stride];
    if (i % 2 == 1) tempbank[n/2 + i/2] = x[i * stride];
  }
  for (int i = 0; i < n; i++) x[i * stride] = tempbank[i];
  free(tempbank);
}

void ilwt97(double* x, int n, int stride)
{
  //unpack
  double* tempbank = (double*)malloc(n * sizeof(double));
  for (int i = 0; i < n/2; i++)
  {
    tempbank[i * 2] = x[i * stride];
    tempbank[i * 2 + 1] = x[(i + n/2) * stride];
  }
  for (int i = 0; i < n; i++) x[i * stride] = tempbank[i];
  free(tempbank);
  
  //undo scale
  double a = 1.149604398;
  for (int i = 0; i < n; i++)
  {
    if (i % 2 == 0) x[i * stride] /= a;
    if (i % 2 == 1) x[i * stride] *= -a;
#if !USE_SQRT2_NORMALIZATION
    x[i * stride] *= sqrt(2);
#endif
  }
  
  //undo update 2
  a = -0.4435068522;
  for (int i = 2 * stride; i < (n - 1) * stride; i += 2*stride)
  {
    x[i] += a * (x[i - stride] + x[i + stride]);
  }
  x[0] += 2 * a * x[stride];
  
  //undo predict 2
  a = -0.8829110762;
  for (int i = stride; i < (n - 2) * stride; i += 2*stride)
  {
    x[i] += a * (x[i - stride] + x[i + stride]);
  }
  x[(n - 1) * stride] += 2 * a * x[(n - 2) * stride];
  
  //undo update 1
  a = 0.05298011854;
  for (int i = 2 * stride; i < (n - 1) * stride; i += 2*stride)
  {
    x[i] += a * (x[i - stride] + x[i + stride]);
  }
  x[0] += 2 * a * x[stride];
  
  //undo predict 1
  a = 1.586134342;
  for (int i = stride; i < (n - 2) * stride; i += 2*stride)
  {
    x[i] += a * (x[i - stride] + x[i + stride]);
  } 
  x[(n - 1) * stride] += 2 * a * x[(n - 2) * stride];
}

void lwt97_2d(double* x, int width, int height, int stride)
{
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    lwt97(x + rowidx * stride, width, 1);
  }
  for (int colidx = 0; colidx < width; colidx++)
  {
    lwt97(x + colidx, height, stride);
  }
}

void ilwt97_2d(double* x, int width, int height, int stride)
{
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    ilwt97(x + rowidx * stride, width, 1);
  }
  for (int colidx = 0; colidx < width; colidx++)
  {
    ilwt97(x + colidx, height, stride);
  }
}