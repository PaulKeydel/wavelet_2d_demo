//forward declaration of 2D CDF9/7 wavelet transformation methods

//via lifting scheme
void ilwt97_2d(double* x, int width, int height, int stride);
void lwt97_2d(double* x, int width, int height, int stride);

//via convolution
void invconvWT_2d(double* ScalingFilter, int hLength, double* WaveletFilter, int gLength, double* x, int width, int height, int stride);
void convWT_2d(double* ScalingFilter, int hLength, double* WaveletFilter, int gLength, double* x, int width, int height, int stride);