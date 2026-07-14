//Different sets of filters
typedef struct
{
  double h_syn[7];
  double g_ana[7];
  double h_ana[9];
  double g_syn[9];
} FilterSet;

static FilterSet FilterCDF97 = {
  { -0.064538, -0.040688, 0.418091,  0.788485, 0.418091, -0.040688, -0.064538 },
  { -0.064538,  0.040688, 0.418091, -0.788485, 0.418091,  0.040688, -0.064538 },
  {  0.037827, -0.023849, -0.110624, 0.377403,  0.852699, 0.377403, -0.110624, -0.023849,  0.037827 },
  { -0.037827, -0.023849,  0.110624, 0.377403, -0.852699, 0.377403,  0.110624, -0.023849, -0.037827 }
};

static FilterSet FilterTaubmann = {
  { -0.091270, -0.057542,  0.591270, 1.115086,  0.591270, -0.057542, -0.091270 },
  { 0.045635, -0.028771, -0.295635, 0.557543, -0.295635, -0.028771,  0.045635 },
  { 0.026748, -0.016864, -0.078223,  0.266864, 0.602949,  0.266864, -0.078223, -0.016864, 0.026748 },
  { 0.053496,  0.033728, -0.156446, -0.533728, 1.205898, -0.533728, -0.156446,  0.033728, 0.053496 }
};

//forward declaration of 2D CDF9/7 wavelet transformation methods

//via lifting scheme
void ilwt97_2d(double* x, int width, int height, int stride);
void lwt97_2d(double* x, int width, int height, int stride);

//via convolution
void invconvWT_2d(FilterSet* filters, double* x, int width, int height, int stride);
void convWT_2d(FilterSet* filters, double* x, int width, int height, int stride);