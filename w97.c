#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double* tempbank = 0;

void convWT(double* ScalingFilter, int hLength, double* WaveletFilter, int gLength, double* signal, int n, int stride)
{
  tempbank = (double*)malloc(n * sizeof(double));
  int i,idx,k;
  for (i = 0; i < n; i++)
  {
    tempbank[i] = 0;
  }
  for (idx = 0; idx < n/2; idx++)
  {
    for (i = 0; i <= hLength - 1; i++)
    {
      k = (abs(2*idx+i-(hLength-1)/2) < n) ? abs(2*idx+i-(hLength-1)/2) : 2*(n-1)-abs(2*idx+i-(hLength-1)/2);
      tempbank[idx] += ScalingFilter[i] * signal[k * stride];
    }
    for (i=0; i<=gLength-1; i++)
    {
      k = (abs(2*idx+i-(gLength-3)/2) < n) ? abs(2*idx+i-(gLength-3)/2) : 2*(n-1)-abs(2*idx+i-(gLength-3)/2);
      tempbank[idx+n/2] += WaveletFilter[i] * signal[k * stride];
    }
  }
  for (i = 0; i < n; i++)
  {
    signal[i * stride] = tempbank[i];
  }
  free(tempbank);
}

void convWT_2d(double* ScalingFilter, int hLength, double* WaveletFilter, int gLength, double* x, int width, int height, int stride)
{
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    convWT(ScalingFilter, hLength, WaveletFilter, gLength, x + rowidx * stride, width, 1);
  }
  for (int colidx = 0; colidx < width; colidx++)
  {
    convWT(ScalingFilter, hLength, WaveletFilter, gLength, x + colidx, height, stride);
  }
}

void invconvWT(double* ScalingFilter, int hLength, double* WaveletFilter, int gLength, double* trafo, int n, int stride)
{
  tempbank = (double*)malloc(n * sizeof(double));
  int i,idx,k;
  for (i = 0; i < n; i++)
  {
    tempbank[i] = 0;
  }
  for (idx = 0; idx < n; idx++)
  {
    for (i=idx/2-(hLength-1)/4; 2*i<=(idx+(hLength-1)/2); i++)
    {
      k = (abs(i) < n/2) ? abs(i) : n-1-abs(i);
      tempbank[idx]+=ScalingFilter[idx-2*i+(hLength-1)/2]*trafo[k*stride];
    }
    for (i=idx/2-(gLength+1)/4; 2*i<=(idx+(gLength-3)/2); i++) //oder <= ???
    {
      //k = (abs(i) < length/2) ? abs(i)-1 : length-1-abs(i);
      k=i;
      if (i<0) k=abs(i)-1;
      if (i>=n/2) k=n-i-2;
      tempbank[idx]+=WaveletFilter[idx-2*i+(gLength-3)/2]*trafo[(k+n/2)*stride];
    }
  }
  for (i = 0; i < n; i++)
  {
    trafo[i * stride] = tempbank[i];
  }
  free(tempbank);
}

void invconvWT_2d(double* ScalingFilter, int hLength, double* WaveletFilter, int gLength, double* x, int width, int height, int stride)
{
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
  tempbank = (double *)malloc(n*sizeof(double));
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
  tempbank = (double *)malloc(n*sizeof(double));
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