//set scaling for coefficients
//if the flag is set the usual sqrt(2) normalization will be used in order to preserve the energy (in terms of l2 norm)
//if the flag is not set the analysis filters will be divided by sqrt(2) and the synthesis filters are multiplied by sqrt(2)
#define USE_SQRT2_NORMALIZATION 1

//define macros for transform scaling
#if USE_SQRT2_NORMALIZATION
#define TRAFO_SCALE    1.41421356 //sqrt(2)
#define TRAFO_SCALE_2D 2
#else
#define TRAFO_SCALE    1.0
#define TRAFO_SCALE_2D 1
#endif

//datatype for biorthogonal 4.4 wavelet filter banks
typedef struct
{
  double h_syn[7];
  double g_ana[7];
  double h_ana[9];
  double g_syn[9];
} FilterSetBior44;

//forward declaration of 2D CDF9/7 wavelet transformation methods
//via lifting scheme
void ilwt97_2d(double* x, int width, int height, int stride);
void lwt97_2d(double* x, int width, int height, int stride);

//via convolution
void invconvWT_2d(double* x, int width, int height, int stride);
void convWT_2d(double* x, int width, int height, int stride);