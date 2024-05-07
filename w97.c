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

int estQuantStepSize(double* x, int width, int height, int stride)
{
  double min = x[0];
  double max = x[0];
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      if (x[rowidx * stride + colidx] > max) max = x[rowidx * stride + colidx];
      if (x[rowidx * stride + colidx] < min) min = x[rowidx * stride + colidx];
    }
  }
  int bitdepth = (int)ceil(log2(max - min));
  return (int)pow(2, bitdepth / 2 + 1);
}

void quant(double* x, int width, int height, int stride, int stepsize)
{
  for (int rowidx = 0; rowidx < height; rowidx++)
  {
    for (int colidx = 0; colidx < width; colidx++)
    {
      int sgn = x[rowidx * stride + colidx] < 0 ? -1 : 1;
      x[rowidx * stride + colidx] = sgn * floor(fabs(x[rowidx * stride + colidx]) / stepsize) * stepsize;
    }
  }
}

void dbls_to_file(const char* fname, double* data, int len)
{
  int* int_res = (int*)malloc(len * sizeof(int));
  for (int i = 0; i < len; i++) 
  {
    int_res[i] = (int)round(data[i]);
  }
  FILE *fp = fopen(fname, "wb");
  if(fp == NULL)
  {
    printf("error creating file");
  }
  fwrite((const void*)int_res, sizeof(int), len, fp);
  fclose(fp);
  free(int_res);
}

void dbls_from_file(const char* fname, double* data, int len)
{
  int* int_res = (int*)malloc(len * sizeof(int));
  FILE *fp = fopen(fname, "rb");
  if(fp == NULL)
  {
    printf("error opening file");
  }
  fread((void*)int_res, sizeof(int), len, fp);
  fclose(fp);
  for (int i = 0; i < len; i++) 
  {
    data[i] = (double)int_res[i];
  }
  free(int_res);
}

int main(int argc, char **argv)
{
  if (argc != 4)
  {
    printf("Not enough parameter!\n");
    printf("Usage: w97 <image_file> <width> <height>\n");
    return -1;
  }

  //declarations
  int width = atoi(argv[2]);
  int height = atoi(argv[3]);
  double* x = (double*)malloc(width * height * sizeof(double));
  double* y = (double*)malloc(width * height * sizeof(double));
  int i,j,s;

  //load data and copy to orig.bin
  dbls_from_file(argv[1], x, width * height);
  dbls_from_file(argv[1], y, width * height);
  dbls_to_file("orig.bin", x, width * height);
    
  //filter set sqrt(2):sqrt(2) normalization
  double h_syn[7] = { -0.064538, -0.040688, 0.418091,  0.788485, 0.418091, -0.040688, -0.064538 };
  double g_ana[7] = { -0.064538,  0.040688, 0.418091, -0.788485, 0.418091,  0.040688, -0.064538 };
  double h_ana[9] = {  0.037827, -0.023849, -0.110624, 0.377403,  0.852699, 0.377403, -0.110624, -0.023849,  0.037827};
  double g_syn[9] = { -0.037827, -0.023849,  0.110624, 0.377403, -0.852699, 0.377403,  0.110624, -0.023849, -0.037827};
    
  //Taubmann:
  //h_syn: -0.091270, -0.057542,  0.591270, 1.115086,  0.591270, -0.057542, -0.091270
  //g_ana:  0.045635, -0.028771, -0.295635, 0.557543, -0.295635, -0.028771,  0.045635
  //h_ana: 0.026748, -0.016864, -0.078223,  0.266864, 0.602949,  0.266864, -0.078223, -0.016864, 0.026748
  //g_syn: 0.053496,  0.033728, -0.156446, -0.533728, 1.205898, -0.533728, -0.156446,  0.033728, 0.053496

  //read trafo levels
  printf("Trafo levels: ");
  scanf("%d", &s);
  printf("\n");

  //estimate step size for quantization depending on input bitdepth
  int stepsize0 = estQuantStepSize(x, width, height, width);
  printf("Quantization step size: %d\n", stepsize0);
  printf("\n");

  //do the forward 9/7 transform
  for(int i = 1; i <= s; i++)
  {
    lwt97_2d(x, width/pow(2,i-1), height/pow(2,i-1), width);
  
    convWT_2d(h_ana, 9, g_ana, 7, y, width/pow(2,i-1), height/pow(2,i-1), width);
  }
  
  //print the wavelet coefficients
  printf("Wavelets coefficients lifting scheme:\n");
  for (i = 0; i < 5; i++)
  {
    printf("%f ", x[i]);
  }
  printf(" ... ");
  for (i = width*height - 5; i < width*height; i++)
  {
    printf("%f ", x[i]);
  }
  printf("\n");
  printf("Wavelets coefficients convolution method:\n");
  for (i = 0; i < 5; i++)
  {
    printf("%f ", y[i]);
  }
  printf(" ... ");
  for (i = width*height - 5; i < width*height; i++)
  {
    printf("%f ", y[i]);
  }
  printf("\n");

  //quantization and safe coeffs to files
  quant(x, width, height, width, stepsize0);
  dbls_to_file("coeffs_quant.bin", x, width * height);
  dbls_to_file("coeffs.bin", y, width * height);

  //do the inverse 9/7 transform
  for(int i = s; i >= 1; i--)
  {
    ilwt97_2d(x, width/pow(2,i-1), height/pow(2,i-1), width);
      
    invconvWT_2d(h_syn, 7, g_syn, 9, y, width/pow(2,i-1), height/pow(2,i-1), width);
  }
  dbls_to_file("reco.bin", x, width * height);
  
  free(x);
  free(y);
  return 0;
}