
#include "mex.h"
#include "math.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "library.h"
#include "frot.h"
#include "fproj.h"


#define ABS(x)    (((x) > 0) ? (x) : (-(x)))


/* InitSigma gives the amount of smoothing applied to the image at the
first level of each octave.  In effect, this determines the sampling
needed in the image domain relative to amount of smoothing.  Good
values determined experimentally are in the range 1.2 to 1.8.
*/
/* float InitSigma_aa = 1.0;*/
static float InitSigma_aa = 1.6;

#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

/* Gaussian convolution kernels are truncated at this many sigmas from
the center.  While it is more efficient to keep this value small,
experiments show that for consistent scale-space analysis it needs
a value of about 3.0, at which point the Gaussian has fallen to
only 1% of its central value.  A value of 2.0 greatly reduces
keypoint consistency, and a value of 4.0 is better than 3.0.
*/
const float GaussTruncate1 = 4.0;


/* --------------------------- Blur image --------------------------- */


/* Same as ConvBuffer, but implemented with loop unrolling for increased
speed.  This is the most time intensive routine in keypoint detection,
so deserves careful attention to efficiency.  Loop unrolling simply
sums 5 multiplications at a time to allow the compiler to schedule
operations better and avoid loop overhead.  This almost triples
speed of previous version on a Pentium with gcc.
*/


void ConvBufferFast(float *buffer, float *kernel, int rsize, int ksize)
{
  int i;
  float *bp, *kp, *endkp;
  float sum;

  for (i = 0; i < rsize; i++) {
    sum = 0.0;
    bp = &buffer[i];
    kp = &kernel[0];
    endkp = &kernel[ksize];

    /* Loop unrolling: do 5 multiplications at a time. */
    //      while (kp + 4 < endkp) {
    //      sum += (double) bp[0] * (double) kp[0] + (double)  bp[1] * (double) kp[1] + (double) bp[2] * (double) kp[2] +
    //             (double) bp[3] * (double) kp[3] + (double)  bp[4] * (double) kp[4];
    //      bp += 5;
    //      kp += 5;
    //     }
    //      /* Do 2 multiplications at a time on remaining items. */
    //     while (kp + 1 < endkp) {
    //       sum += (double) bp[0] * (double) kp[0] + (double)  bp[1] * (double) kp[1];
    //	   bp += 2;
    //	   kp += 2;
    //	 }
    //      /* Finish last one if needed. */
    //		if (kp < endkp) {
    //		sum += (double) *bp * (double) *kp;
    //		}

    while (kp < endkp) {
      sum += *bp++ * *kp++;
    }

    buffer[i] = sum;
  }
}

/* Convolve image with the 1-D kernel vector along image rows.  This
is designed to be as efficient as possible.  Pixels outside the
image are set to the value of the closest image pixel.
*/
void ConvHorizontal(vector<float>& image, int width, int height, float *kernel, int ksize)
{
  int rows, cols, r, c, i, halfsize;
  float buffer[4000];
  vector<float> pixels(width*height);


  rows = height;
  cols = width;

  halfsize = ksize / 2;
  pixels = image;
  assert(cols + ksize < 4000);

  for (r = 0; r < rows; r++) {
    /* Copy the row into buffer with pixels at ends replicated for
    half the mask size.  This avoids need to check for ends
    within inner loop. */
    for (i = 0; i < halfsize; i++)
      buffer[i] = pixels[r*cols];
    for (i = 0; i < cols; i++)
      buffer[halfsize + i] = pixels[r*cols+i];
    for (i = 0; i < halfsize; i++)
      buffer[halfsize + cols + i] = pixels[r*cols+cols-1];

    ConvBufferFast(buffer, kernel, cols, ksize);
    for (c = 0; c < cols; c++)
      pixels[r*cols+c] = buffer[c];
  }
  image = pixels;
}


/* Same as ConvHorizontal, but apply to vertical columns of image.
*/
void ConvVertical(vector<float>& image, int width, int height, float *kernel, int ksize)
{
  int rows, cols, r, c, i, halfsize;
  float buffer[4000];
  vector<float> pixels(width*height);

  rows = height;
  cols = width;

  halfsize = ksize / 2;
  pixels = image;
  assert(rows + ksize < 4000);

  for (c = 0; c < cols; c++) {
    for (i = 0; i < halfsize; i++)
      buffer[i] = pixels[c];
    for (i = 0; i < rows; i++)
      buffer[halfsize + i] = pixels[i*cols+c];
    for (i = 0; i < halfsize; i++)
      buffer[halfsize + rows + i] = pixels[(rows - 1)*cols+c];

    ConvBufferFast(buffer, kernel, rows, ksize);
    for (r = 0; r < rows; r++)
      pixels[r*cols+c] = buffer[r];
  }

  image = pixels;
}



/* 1D Convolve image with a Gaussian of width sigma and store result back
in image.   This routine creates the Gaussian kernel, and then applies
it in horizontal (flag_dir=0) OR vertical directions (flag_dir!=0).
*/
void GaussianBlur1D(vector<float>& image, int width, int height, float sigma, int flag_dir)
{
  float x, kernel[100], sum = 0.0;
  int ksize, i;

  /* The Gaussian kernel is truncated at GaussTruncate sigmas from
  center.  The kernel size should be odd.
  */
  ksize = (int)(2.0 * GaussTruncate1 * sigma + 1.0);
  ksize = MAX(3, ksize);    /* Kernel must be at least 3. */
  if (ksize % 2 == 0)       /* Make kernel size odd. */
    ksize++;
  assert(ksize < 100);

  /* Fill in kernel values. */
  for (i = 0; i <= ksize; i++) {
    x = i - ksize / 2;
    kernel[i] = exp(- x * x / (2.0 * sigma * sigma));
    sum += kernel[i];
  }
  /* Normalize kernel values to sum to 1.0. */
  for (i = 0; i < ksize; i++)
    kernel[i] /= sum;

  if (flag_dir == 0)
  {
    ConvHorizontal(image, width, height, kernel, ksize);
  }
  else
  {
    ConvVertical(image, width, height, kernel, ksize);
  }
}




void perform_tilt_and_rotation( vector<float>& image, int width, int height, vector<float>& image_to_return, int& width_t, int& height_t, float theta, float t)
{

    int flag_dir = 1;
    int fproj_o;
    float fproj_p, fproj_bg;
    char fproj_i;
    float *fproj_x4, *fproj_y4;
    //  float frot_b=0;
    float frot_b=128;
    char *frot_k;

    frot_k = 0; // order of interpolation


    fproj_o = 3;
    fproj_p = 0;
    fproj_i = 0;
    fproj_bg = 0;
    fproj_x4 = 0;
    fproj_y4 = 0;




    //theta = theta * 180 / PI;
    float t1 = 1;
    float t2 = 1/t;

    vector<float> image_t;
    int width_r, height_r;

    // simulate a rotation: rotate the image with an angle theta. (the outside of the rotated image are padded with the value frot_b)
    frot(image, image_t, width, height, &width_r, &height_r, &theta, &frot_b , frot_k);

    /* Tilt */
    width_t = (int) (width_r * t1);
    height_t = (int) (height_r * t2);

    int fproj_sx = width_t;
    int fproj_sy = height_t;

    float fproj_x1 = 0;
    float fproj_y1 = 0;
    float fproj_x2 = width_t;
    float fproj_y2 = 0;
    float fproj_x3 = 0;
    float fproj_y3 = height_t;

    /* Anti-aliasing filtering along vertical direction */
    /* sigma_aa = InitSigma_aa * log2(t);*/
    float sigma_aa = InitSigma_aa * t / 2;
    GaussianBlur1D(image_t,width_r,height_r,sigma_aa,flag_dir);


    // simulate a tilt: subsample the image along the vertical axis by a factor of t.
    vector<float> image_tmp(width_t*height_t);
    fproj (image_t, image_tmp, width_r, height_r, &fproj_sx, &fproj_sy, &fproj_bg, &fproj_o, &fproj_p, &fproj_i , fproj_x1 , fproj_y1 , fproj_x2 , fproj_y2 , fproj_x3 , fproj_y3, fproj_x4, fproj_y4);

    image_to_return = image_tmp;

}










// Fonction principale (gère la liaison avec Matlab)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
const mxArray *prhs[])
{
  /* Vérification du nombre d'arguments */
    if ( nrhs != 4 ) {
        //mexPrintf("Wrong number of input arguments!");
    } else if ((nlhs != 2)) {
        //error_msg("Wrong number of output arguments!");
    }

  /* Vérification du type des arguments */
    if (!mxIsDouble(prhs[0])||mxIsComplex(prhs[0])||!mxIsDouble(prhs[1])||mxIsComplex(prhs[1])||!mxIsDouble(prhs[2])||mxIsComplex(prhs[2])||mxGetM(prhs[2])!=1||mxGetN(prhs[2])!=1)
    {
        //error_msg("Wrong type of input arguments!");
    }


  /* matrices d'entrée*/
    int w1,h1,w2,h2;  // W = number of rows, H = number of columns
    w1 = mxGetM(prhs[0]);
    h1 = mxGetN(prhs[0]);

    vector<float> ipixels1(mxGetPr(prhs[0]), mxGetPr(prhs[0]) + w1*h1);
    double* t = mxGetPr(prhs[1]);
    double* theta = mxGetPr(prhs[2]);

    vector<float> ipixels2; //output image

    perform_tilt_and_rotation(ipixels1, w1, h1, ipixels2, w2, h2, *theta, *t);

    //output image
    plhs[0] = mxCreateDoubleMatrix(w2, h2, mxREAL);
    double *data = mxGetPr(plhs[0]);
    for(int j = 0; j < (int) h2; j++)
            for(int i = 0; i < (int) w2; i++)  data[j*w2+i] = ipixels2[j*w2+i];
}

