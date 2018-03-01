
#include "mex.h"
#include "math.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>


#include "libNumerics/numerics.h"

#define ABS(x)    (((x) > 0) ? (x) : (-(x)))

/**
 * @brief Computes tilted coordinates \f$ (x^\prime,y^\prime) \in [1,m_x]\times[1,m_y]\f$ from image coordinates \f$(x,y) \in [1,n_x]\times[1,n_y]\f$ where
 * \f$ T_tR_\theta (x,y) = (x^\prime,y^\prime) \f$. This takes into account the fact that a digital rotation \f$R_\theta\f$ applies a translation to reframe the image.
 * @param (x,y) Enters image coordiantes and returns tilted coordinates
 * @param w Width from the original image (\f$ w = n_x \f$)
 * @param h Width from the original image (\f$ h = n_y \f$)
 * @param t Tilt applied in the y direction
 * @param theta Direction \f$\theta\f$ of the tilt
 * @author Mariano Rodríguez
 */
void imagecoor2tiltedcoor(float& x, float& y, const int& w, const int& h, const float& t, const float& theta)
{
    // Apply R_\theta

    // Prepare rotation matrix
    libNumerics::matrix<float> Rot(2,2), vec(2,1), corners(2,4);
    // Rot = [cos(Rtheta) -sin(Rtheta);sin(Rtheta) cos(Rtheta)];
    Rot(0,0) = cos(theta); Rot(0,1) = sin(theta);
    Rot(1,0) = -sin(theta); Rot(1,1) = cos(theta);


    vec(0,0) = x-1;
    vec(1,0) = y-1;

    // rotation -> [x1;y1] = Rot*[x1;y1]
    vec = (Rot*vec);
    x = vec(0,0);
    y = vec(1,0);


    // Translate so that image borders are in the first quadrant?
    float  x_ori, y_ori;
    // A = Rot*[ [0;h1] [w1;0] [w1;h1] [0;0] ];
    corners(0,0) = 0;   corners(0,1) = w-1;   corners(0,2) = w-1;   corners(0,3) = 0;
    corners(1,0) = 0;  corners(1,1) = 0;    corners(1,2) = h-1;   corners(1,3) = h-1;

    corners = Rot*corners;

    //x_ori = min(corners(1,:));
    //y_ori = min(corners(2,:));
    x_ori = corners(0,0); y_ori = corners(1,0);
    for(int i=1; i<4; i++)
    {
        if (x_ori>corners(0,i))
            x_ori = corners(0,i);

        if (y_ori>corners(1,i))
            y_ori = corners(1,i);
    }

    // translation and Tilt T_t
    x = x - x_ori + 1;
    y = (y - y_ori)/t + 1;
}


/**
 * @brief Computes image coordinates \f$(x,y) \in [1,n_x]\times[1,n_y]\f$ from tilted coordinates \f$T_tR_\theta(x,y) \in [1,m_x]\times[1,m_y]\f$
 * @param (x,y) Enters image coordinates and returns tilted coordinates
 * @param w Width from the original image (\f$ w = n_x \f$)
 * @param h Width from the original image (\f$ h = n_y \f$)
 * @param t Tilt that was applied in the y direction
 * @param theta Direction \f$\theta\f$ in which the tilt that was applied
 * @return True if the computed image-coordinates fall inside the true image boundaries.
 * @author Mariano Rodríguez
 */
bool tiltedcoor2imagecoor(float& x, float& y, const int& w, const int& h, const float& threshold , const float& t, const float& theta)
{

    // Get initial translation \tau_{x_ori,y_ori}

    // Prepare rotation matrix
    libNumerics::matrix<float> Rot(2,2), vec(2,1), corners(2,4);
    // Rot = [cos(Rtheta) -sin(Rtheta);sin(Rtheta) cos(Rtheta)];
    Rot(0,0) = cos(theta); Rot(0,1) = sin(theta);
    Rot(1,0) = -sin(theta); Rot(1,1) = cos(theta);

    // Translate so that image borders are in the first quadrant?
    float  x_ori, y_ori;
    // A = Rot*[ [0;h-1] [w-1;0] [w-1;h-1] [0;0] ];
    corners(0,0) = 0;   corners(0,1) = w-1;   corners(0,2) = w-1;   corners(0,3) = 0;
    corners(1,0) = 0;  corners(1,1) = 0;    corners(1,2) = h-1;   corners(1,3) = h-1;
    corners = Rot*corners;

    //x_ori = min(corners(1,:));
    //y_ori = min(corners(2,:));
    x_ori = corners(0,0); y_ori = corners(1,0);
    for(int i=1; i<4; i++)
    {
        if (x_ori>corners(0,i))
            x_ori = corners(0,i);

        if (y_ori>corners(1,i))
            y_ori = corners(1,i);
    }

    // Distance from the point to the true image borders
    // if less than threshold stop computations and return false
    float xvec[5], yvec[5], d;
    for(int i=0; i<4; i++)
    {// translation and Tilt T_t
        xvec[i] = corners(0,i) - x_ori + 1;
        yvec[i] = (corners(1,i) - y_ori)/t + 1;
        //cout<<"xvec="<<xvec[i]<< " yvec="<<yvec[i]<<endl;
    }
    xvec[4] = xvec[0];
    yvec[4] = yvec[0];
    //cout<<"xvec="<<xvec[4]<< " yvec="<<yvec[4]<<endl;

    for(int i=0; i<4; i++)
    {
        d = ABS( (yvec[i+1]-yvec[i])*x - (xvec[i+1]-xvec[i])*y + xvec[i+1]*yvec[i] - yvec[i+1]*xvec[i] ) / sqrt( pow(xvec[i+1]-xvec[i],2) + pow(yvec[i+1]-yvec[i],2) );
        //cout<<"d="<<d<<endl;
        if (d <= threshold)
            return(false);
    }


    // Inverse of Tilt T_t and translation
    x = (x-1)   + x_ori;
    y = (y-1)*t + y_ori;

    vec(0,0) = x;
    vec(1,0) = y;

    // rotation -> [x1;y1] = Rot^{-1}*[x1;y1]
    Rot(1,0) = -Rot(1,0); Rot(0,1) = -Rot(0,1);
    vec = (Rot*vec);
    x = vec(0,0) + 1;
    y = vec(1,0) + 1;

    if ( (x>w)||(y>h)||(x<1)||(y<1) )
        return(false);
    else
        return(true);
}


// Fonction principale (gère la liaison avec Matlab)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
const mxArray *prhs[])
{

    double* temp = mxGetPr(prhs[0]); float x = (float) *temp;
    temp = mxGetPr(prhs[1]); float y = (float) *temp;
    temp = mxGetPr(prhs[2]); int w = (int) *temp;
    temp = mxGetPr(prhs[3]); int h = (int) *temp;
    temp = mxGetPr(prhs[4]); float t = (float) *temp;
    temp = mxGetPr(prhs[5]); float theta = (float) *temp;

    temp =  mxGetPr(prhs[6]); int direction = (int) *temp;
    temp = mxGetPr(prhs[7]); float  threshold = (float) *temp;
    
    bool inside = true;


    switch (direction) {
      case 1:
      {
        imagecoor2tiltedcoor(x, y, w, h, t, theta);
        break;
      }
      case -1:
      {
        inside = tiltedcoor2imagecoor(x, y, w, h, threshold , t, theta);
        break;
      }
    }

    //output data
    plhs[0] = mxCreateDoubleMatrix(1, 3, mxREAL);
    double *data = mxGetPr(plhs[0]);
    data[0] = x;
    data[1] = y;
    if (inside)
      data[2] = 1.0;
    else
      data[2] = -1.0;


}
