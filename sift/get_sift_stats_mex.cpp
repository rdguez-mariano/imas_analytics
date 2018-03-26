#include "mex.h"
#include <string>
#include <sstream>
#include <iostream>
#include <omp.h>


#include <opencv2/core.hpp>
#include "opencv2/opencv.hpp"
#include "opencv2/opencv_modules.hpp"
#include "opencv2/core/core.hpp"
#include "opencv2/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include <opencv2/xfeatures2d.hpp>
#include "opencv2/flann/miniflann.hpp"

void floatarray2opencvimage(double *input_image,cv::Mat& output_image,int width, int height)
{
  output_image.create(height, width, CV_8UC1);//cv::imread(argv[2]);

  for (int i = 0;  i < output_image.rows; i++)
  {
    for (int j = 0; j < output_image.cols; j++)
    {
      output_image.data[((output_image.cols)*i)+j] = (uchar) floor(input_image[((output_image.cols)*i)+j]);
    }
  }
}


static inline void unpackOctave(const cv::KeyPoint& kpt, int& octave, int& layer, float& scale)
{
  octave = kpt.octave & 255;
  layer = (kpt.octave >> 8) & 255;
  octave = octave < 128 ? octave : (-128 | octave);
  scale = octave >= 0 ? 1.f/(1 << octave) : (float)(1 << -octave);
}


// Fonction principale (gère la liaison avec Matlab)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
  const mxArray *prhs[])
  {
    /* matrices d'entrée*/
    long w1,h1;  // W = number of rows, H = number of columns
    w1 = mxGetM(prhs[0]);
    h1 = mxGetN(prhs[0]);
    double* ipixels1 = mxGetPr(prhs[0]);
    const float groupratio = 0.8;

    cv::Mat queryImg;
    floatarray2opencvimage(ipixels1,queryImg,w1,h1);

    cv::Ptr<cv::FeatureDetector> detector = cv::xfeatures2d::SIFT::create();
    cv::Ptr<cv::DescriptorExtractor> extractor = cv::xfeatures2d::SIFT::create();

    cv::Mat dlist;
    std::vector<cv::KeyPoint> klist;
    detector->detect(queryImg, klist);
    extractor->compute(queryImg, klist, dlist);

    //RootSIFT
    if (false)
    {
      int rows;
      rows = dlist.rows;

      for (int ii=0;ii<rows;ii++)
      {
        cv::Mat temp,temp2;
        cv::normalize( dlist.row(ii), temp, 1, cv::NORM_L1);
        cv::sqrt(temp, temp2);
        temp2.row(0).copyTo(dlist.row(ii));
      }

    }


    float	dsq, distsq1, distsq2;
    float* mapratio = new float[klist.size()];
    int* groups = new int[klist.size()];
    for (int i=0; i< (int) klist.size(); i++)
      groups[i] = 0;
    int groupcont = 1;

    for (int i=0; i< (int) klist.size(); i++)
    {
      distsq1 = distsq2 = 99999999999;
      for (int j=0; j< (int) klist.size(); j++)
      {
        if(i!=j)
        {
          dsq = (float)cv::norm(dlist.row(i),dlist.row(j),cv::NORM_L2);

          if (dsq < distsq1) {
            distsq2 = distsq1;
            distsq1 = dsq;
          } else if (dsq < distsq2) {
            distsq2 = dsq;
          }
        }
      }
      mapratio[i] = distsq1/distsq2;

      if (groups[i]!=0)
        continue;

      std::vector<int> cgroup; //current group
      cgroup.clear();
      for (int j=0; j< (int) klist.size(); j++)
      {
        if((i!=j)&&(groups[j]==0))
        {
          dsq = (float)cv::norm(dlist.row(i),dlist.row(j),cv::NORM_L2);
          if (distsq1/dsq > groupratio)
          cgroup.push_back(j);
        }
      }

      if (cgroup.size()>=1)
      {
        groups[i] = groupcont;
        for (int j=0;j<(int)cgroup.size();j++)
          groups[cgroup[j]] = groupcont;
        groupcont++;
      }

    }

    // Matchings Output Matrix
    int wo = 9;
    plhs[0] = mxCreateDoubleMatrix(wo, klist.size(), mxREAL);
    double *data = mxGetPr(plhs[0]);

    for ( int i = 0; i < (int) klist.size(); i++)
    {
      cv::KeyPoint* kp = & klist[i];
      int octave, layer;
      float scale;
      unpackOctave(*kp, octave, layer, scale);

      data[i*wo+0]=(double)kp->pt.x;
      data[i*wo+1]=(double)kp->pt.y;
      data[i*wo+2]=(double)kp->angle;
      data[i*wo+3]= (double) octave;
      data[i*wo+4]= (double) layer;
      data[i*wo+5]= (double) scale;
      data[i*wo+6]= (double) kp->size;
      data[i*wo+7]= (double) mapratio[i];
      data[i*wo+8]= (double) groups[i];
    }
    delete[] mapratio;
  }
