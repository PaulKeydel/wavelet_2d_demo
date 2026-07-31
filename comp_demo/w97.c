#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "w97.h"


void convWT(double* ScalingFilter, int hLength, double* WaveletFilter, int gLength, double* signal, int n, int stride)
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

void invconvWT(double* ScalingFilter, int hLength, double* WaveletFilter, int gLength, double* trafo, int n, int stride)
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

void convWT_2d(FilterSet* filters, double* x, int width, int height, int stride)
{
  double* ScalingFilter = filters->h_ana;
  int hLength = 9;
  double* WaveletFilter = filters->g_ana;
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

void invconvWT_2d(FilterSet* filters, double* x, int width, int height, int stride)
{
  double* ScalingFilter = filters->h_syn;
  int hLength = 7;
  double* WaveletFilter = filters->g_syn;
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
  double a;
  int i;
    
  //predict 1
  a = -1.586134342;
  for (i = 1*stride; i < (n-2)*stride; i += 2*stride)
  {
    x[i] += a*(x[i-stride]+x[i+stride]);
  }
  x[(n-1)*stride] += 2*a*x[(n-2)*stride];
  
  //update 1
  a = -0.05298011854;
  for (i = 2*stride; i < (n-1)*stride; i += 2*stride)
  {
    x[i] += a*(x[i-stride]+x[i+stride]);
  }
  x[0] += 2*a*x[stride];
  
  //predict 2
  a = 0.8829110762;
  for (i = 1*stride; i < (n-2)*stride; i += 2*stride)
  {
    x[i] += a*(x[i-stride]+x[i+stride]);
  }
  x[(n-1)*stride] += 2*a*x[(n-2)*stride];
    
  //update 2
  a = 0.4435068522;
  for (i = 2*stride; i < (n-1)*stride; i += 2*stride)
  {
    x[i] += a*(x[i-stride]+x[i+stride]);
  }
  x[0] += 2*a*x[stride];
  
  //scale
  //a=1.230174;
  a = 1/1.149604398; //sqrt2 normalization
  for (i = 0; i < n ; i++)
  {
      if (i%2) x[i*stride] *= -a; //divide by 2 for Taubmann
      else x[i*stride] /= a;
  }
  
  //pack
  double* tempbank = (double *)malloc(n*sizeof(double));
  for (i = 0; i < n; i++)
  {
    if (i%2 == 0) tempbank[i/2] = x[i*stride];
    else tempbank[n/2+i/2] = x[i*stride];
  }
  for (i = 0; i < n; i++) x[i*stride] = tempbank[i];
  free(tempbank);
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

void ilwt97(double* x, int n, int stride)
{
  double a;
  int i;
  
  //unpack
  double* tempbank = (double *)malloc(n*sizeof(double));
  for (i = 0; i < n/2; i++)
  {
    tempbank[i*2] = x[i*stride];
    tempbank[i*2 + 1] = x[(i+n/2)*stride];
  }
  for (i = 0; i < n; i++) x[i*stride] = tempbank[i];
  free(tempbank);
  
  //undo scale
  //a=1/1.230174;
  a = 1.149604398; //sqrt2 normalization
  for (i = 0; i < n; i++)
  {
    if (i%2) x[i*stride] *= -a; //multiply by 2a for Taubmann
    else x[i*stride] /= a;
  }
  
  //undo update 2
  a = -0.4435068522;
  for (i = 2*stride; i < (n-1)*stride; i += 2*stride)
  {
    x[i] += a*(x[i-stride]+x[i+stride]);
  }
  x[0] += 2*a*x[stride];
  
  //undo predict 2
  a = -0.8829110762;
  for (i = 1*stride; i < (n-2)*stride; i += 2*stride)
  {
    x[i] += a*(x[i-stride]+x[i+stride]);
  }
  x[(n-1)*stride] += 2*a*x[(n-2)*stride];
  
  //undo update 1
  a = 0.05298011854;
  for (i = 2*stride; i < (n-1)*stride; i += 2*stride)
  {
    x[i] += a*(x[i-stride]+x[i+stride]);
  }
  x[0] += 2*a*x[stride];
  
  //undo predict 1
  a = 1.586134342;
  for (i = 1*stride; i < (n-2)*stride; i += 2*stride)
  {
    x[i] += a*(x[i-stride]+x[i+stride]);
  } 
  x[(n-1)*stride] += 2*a*x[(n-2)*stride];
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