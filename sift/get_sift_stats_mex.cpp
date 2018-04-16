#ifndef _NOMEX
#include "mex.h"
#endif

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

template<class T>
struct ClusterNode;

template<class T>
struct KpNode
{
  KpNode<T>* next;
  //KpNode* next;
  ClusterNode<T>* membership;
  T data; //this was my kp_id
  KpNode(const T &x);
  ~KpNode();
  static void UpdateMembership(KpNode<T>* k,ClusterNode<T>* c);
private:
  KpNode();
};

template<class T>
struct ClusterNode
{
  int number_of_KpNodes;
  KpNode<T>* first_Kp;
  KpNode<T>* last_Kp;
  ClusterNode<T>* next;
  ClusterNode<T>* prev;
  ClusterNode(KpNode<T>* Kp);
  ~ClusterNode();

  ClusterNode(); // don't use this constructor -> only for first cluster node
};

template<class T>
ClusterNode<T>::~ClusterNode()
{
  if (first_Kp!=0)
    delete first_Kp;
  if (next!=0)
    delete next;
}

template<class T>
ClusterNode<T>::ClusterNode()
{
  number_of_KpNodes = 0;
  next = 0;
  prev = 0;
  first_Kp = 0;
  last_Kp = 0;
}

template<class T>
ClusterNode<T>::ClusterNode(KpNode<T>* Kp)
:number_of_KpNodes(1)
{
  next = 0;
  prev = 0;
  first_Kp = Kp;
  last_Kp = Kp;
  Kp->membership = this;
}

template<class T>
KpNode<T>::KpNode(const T &x)
:data(x)
{
  next = 0;
}

template<class T>
KpNode<T>::~KpNode()
{
  if (next!=0)
  delete next;
}


template<class T>
void KpNode<T>::UpdateMembership(KpNode<T>* k, ClusterNode<T>* c)
{
  if(k->membership!=0)
  {
    k->membership->first_Kp = 0;
    k->membership->last_Kp = 0;
  }
  KpNode<T> * current = k;
  do
  {
    current->membership = c;
    current = current->next;
  }
  while(current!=0);
}


class irr_classes
{
public:
  int number_of_clusters;
  ClusterNode<int>* first_cluster;
  ClusterNode<int>* last_cluster;
  bool search_id(int kp_id, KpNode<int>** n);
  void add_singleton_cluster(int kp_id);
  void Merge(ClusterNode<int>* c1, ClusterNode<int>* c2);
  void PrintCluster();
  irr_classes();

  // linked by the same id
  std::vector<int> id_in;
  std::vector<KpNode<int>*> nodes;
};

void irr_classes::PrintCluster()
{
  int ClusterCount = 0;
  ClusterNode<int>* current_cluster = first_cluster;
  do {
    std::cout<<" Cluster id = "<<ClusterCount<<", N. Nodes = "<<current_cluster->number_of_KpNodes<<std::endl<<"{";
    KpNode<int>* current_node = current_cluster->first_Kp;

    while(current_node!=0)
    {
      std::cout<<"  "<<current_node->data;
      current_node = current_node->next;
    }
    std::cout<<"}"<<std::endl;
    current_cluster = current_cluster->next;
    ClusterCount++;
  }
  while(current_cluster!=0);
}

void irr_classes::Merge(ClusterNode<int>* c1, ClusterNode<int>* c2)
{
  if (c1!=c2)
  {
    // c1 is supposed to have more elements that c2
    if (c1->number_of_KpNodes<c2->number_of_KpNodes)
    {
      ClusterNode<int>* temp = c2;
      c2 = c1;
      c1 = temp;
    }

    if (c1->last_Kp==0)
    {
      std::cout<<"not supposed to happed this / number of inner nodes = "<<c1->number_of_KpNodes<<std::endl;
    }
    c1->last_Kp->next = c2->first_Kp;
    c1->last_Kp = c2->last_Kp;
    c1->number_of_KpNodes += c2->number_of_KpNodes;
    //std::cout<<"entering UpdateMembership "<<std::endl;
    if ((c2==last_cluster)&&(c2->prev!=0))
      last_cluster=c2->prev;
    KpNode<int>::UpdateMembership(c2->first_Kp,c1);
    if(c2->prev!=0)
    c2->prev->next = c2->next;
    if(c2->next!=0)
    c2->next->prev = c2->prev;
    c2->next = 0;
    c2->prev = 0;
    delete c2;
  }
}


irr_classes::irr_classes()
{
  number_of_clusters = 0;
  id_in.clear();
  nodes.clear();
  first_cluster = new ClusterNode<int>();
  last_cluster = first_cluster;
}

void irr_classes::add_singleton_cluster(int kp_id)
{
  KpNode<int>* newKp = new KpNode<int>(kp_id);
  //std::cout << kp_id<<std::endl;
  ClusterNode<int>* newCluster = new ClusterNode<int>(newKp);
  if (last_cluster != 0)
  {
    last_cluster->next = newCluster;
    newCluster->prev = last_cluster;
    last_cluster = newCluster;
  }
  else
  {
    std::cout<<"You're not supposed to use it without a the first cluster node"<<std::endl;
  }
  id_in.push_back(kp_id);
  nodes.push_back(newKp);
}

bool irr_classes::search_id(int kp_id, KpNode<int>** n)
{
  ClusterNode<int>* current_cluster = first_cluster;
  do {
    KpNode<int>* current_node = current_cluster->first_Kp;
    while(current_node!=0)
    {
      if (current_node->data ==kp_id)
      {
        *n = current_node;
        return(true);
      }
      current_node = current_node->next;
    }
    current_cluster = current_cluster->next;
  }
  while(current_cluster!=0);

  //if here, it means kp_id wasn't found
  return(false);
}





struct keypoint_simple {
  cv::Point2f pt;
  double	scale, size,
  angle, octave, layer;
  int id;
};
typedef std::pair<keypoint_simple,keypoint_simple> matching;
typedef std::vector<matching> matchingslist;

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

#ifdef _NOMEX
#include <string>
int main(int argc, char **argv)
  {
    /* matrices d'entrée*/
    std::string s1(argv[1]);
    std::string s2(argv[2]);
    cv::Mat queryImg = cv::imread(s1, CV_LOAD_IMAGE_GRAYSCALE );
    cv::Mat acImg = cv::imread(s2, CV_LOAD_IMAGE_GRAYSCALE );
#else
// Fonction principale (gère la liaison avec Matlab)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
  const mxArray *prhs[])
  {
    /* matrices d'entrée*/
    long w1,h1,wac,hac;  // W = number of rows, H = number of columns

    w1 = mxGetM(prhs[0]);
    h1 = mxGetN(prhs[0]);
    double* ipixels1 = mxGetPr(prhs[0]);
    wac = mxGetM(prhs[1]);
    hac = mxGetN(prhs[1]);
    double* ipixelsac = mxGetPr(prhs[1]);
    cv::Mat queryImg,acImg;
    floatarray2opencvimage(ipixels1,queryImg,w1,h1);
    floatarray2opencvimage(ipixelsac,acImg,wac,hac);
#endif


    int nthreads, maxthreads;
    /* Display info on OpenMP*/
    #pragma omp parallel
    {
      #pragma omp master
      {
        nthreads = omp_get_num_threads();
        maxthreads = omp_get_max_threads();
      }
    }
    std::cout<<"--> Using "<<nthreads<<" threads out of "<<maxthreads<<" <--"<<std::endl;

#ifndef _NOMEX
    double matchratio = (double) mxGetScalar(prhs[2]);
    int thr_same_kp = (int) mxGetScalar(prhs[3]);
#else
    double matchratio = atof(argv[3]);
    int thr_same_kp = atoi(argv[4]);
#endif
    cv::Ptr<cv::FeatureDetector> detector = cv::xfeatures2d::SIFT::create();
    cv::Ptr<cv::DescriptorExtractor> extractor = cv::xfeatures2d::SIFT::create();

    cv::Mat dlist,dlistac;
    std::vector<cv::KeyPoint> klist,klistac;

    #pragma omp parallel
    #pragma omp master
    {
      #pragma omp task shared(queryImg, klist, dlist)
      {
        detector->detect(queryImg, klist);
        extractor->compute(queryImg, klist, dlist);
      }
      #pragma omp task shared(acImg, klistac, dlistac)
      {
        detector->detect(acImg, klistac);
        extractor->compute(acImg, klistac, dlistac);
      }
    }

    std::cout<<klist.size()<<" query keypoints  /  "<<klistac.size()<<" ac keypoints"<<std::endl;

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

    //Group analysis
    double	dsq, distsq1;
    matchingslist matchings;
    std::vector<double> ac_distance;
    irr_classes * iclasses = new irr_classes();
    #pragma omp parallel for shared(matchings,ac_distance,iclasses,klist,dlist,klistac,dlistac) private(dsq,distsq1)  schedule(static)
    for (int i=0; i< (int) klist.size(); i++)
    {
      // minimal distance to the a-contrario descriptors
      distsq1 = (double)cv::norm(dlist.row(i),dlistac.row(0),cv::NORM_L2);
      for (int j=1; j< (int) klistac.size(); j++)
      {
        dsq = (double)cv::norm(dlist.row(i),dlistac.row(j),cv::NORM_L2);
        if (dsq < distsq1)
        distsq1 = dsq;
      }

      //see all others
      for (int j=0; j< (int) klist.size(); j++)
      {
        if((i!=j)&&(cv::norm(klist[i].pt-klist[j].pt)>thr_same_kp))
        {
          dsq = (double)cv::norm(dlist.row(i),dlist.row(j),cv::NORM_L2);
          if (dsq < matchratio*distsq1)
          {
            #pragma omp critical
            {
              bool already_in = false;
              for (int m =0;m<matchings.size();m++)
              {
                //if( ((matchings[m].first.id==i)&&((matchings[m].second.id==j))) || ((matchings[m].second.id==i)&&((matchings[m].first.id==j))) )
                if( ((cv::norm(klist[i].pt-matchings[m].first.pt)<thr_same_kp)&&(cv::norm(klist[j].pt-matchings[m].second.pt)<thr_same_kp)) ||
                ((cv::norm(klist[i].pt-matchings[m].second.pt)<thr_same_kp)&&(cv::norm(klist[j].pt-matchings[m].first.pt)<thr_same_kp)) )
                {
                  if(ac_distance[m]>dsq/distsq1)
                  ac_distance[m] = dsq/distsq1;
                  already_in = true;
                  break;// on est censé à trouver que un seul car la liste à un seul representant
                }
              }

              if (!already_in)
              {
                //add Match
                cv::KeyPoint* kp = & klist[i];
                int octave, layer;
                float scale;
                unpackOctave(*kp, octave, layer, scale);

                keypoint_simple k1, k2;

                k1.pt = kp->pt;
                k1.size = (double) kp->size;
                k1.angle = (double)kp->angle;
                k1.scale = (double) scale;
                k1.layer = (double) layer;
                k1.octave = (double) octave;


                kp = & klist[j];
                unpackOctave(*kp, octave, layer, scale);

                k2.pt = kp->pt;
                k2.size = (double) kp->size;
                k2.angle = (double)kp->angle;
                k2.scale = (double) scale;
                k2.layer = (double) layer;
                k2.octave = (double) octave;

                //add new id to id_in to build the similarity matrix

                int i_surr=i, j_surr=j;
                int ind_i = -1, ind_j = -1;
                bool i_id=false,j_id=false;
                for (int m=0;m<iclasses->id_in.size();m++)
                {
                  if(cv::norm(klist[iclasses->id_in[m]].pt-klist[i].pt)<thr_same_kp)
                  {
                    i_id = true;
                    i_surr = iclasses->id_in[m];
                    ind_i = m;
                  }
                  if(cv::norm(klist[iclasses->id_in[m]].pt-klist[j].pt)<thr_same_kp)
                  {
                    j_id = true;
                    j_surr = iclasses->id_in[m];
                    ind_j = m;
                  }
                }
                if (!i_id)
                {
                  iclasses->add_singleton_cluster(i_surr);
                  ind_i = iclasses->id_in.size()-1;
                }
                if (!j_id)
                {
                  iclasses->add_singleton_cluster(j_surr);
                  ind_j = iclasses->id_in.size()-1;
                }

                k1.id = i_surr;
                k2.id = j_surr;

                iclasses->Merge(iclasses->nodes[ind_i]->membership,iclasses->nodes[ind_j]->membership);

                matchings.push_back( matching(k1,k2) );
                ac_distance.push_back(dsq/distsq1);
              }
            }
          }
        }
      }
    }

#ifndef _NOMEX
    // dense keypoints Matrix
    int wo = 7;
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
    }

    // matches
    wo = 17;
    plhs[1] = mxCreateDoubleMatrix(wo, matchings.size(), mxREAL);
    data = mxGetPr(plhs[1]);
    for ( int i = 0; i < (int) matchings.size(); i++)
    {
      data[i*wo+0]=matchings[i].first.pt.x;
      data[i*wo+1]=matchings[i].first.pt.y;
      data[i*wo+2]=matchings[i].first.angle;
      data[i*wo+3]= matchings[i].first.octave;
      data[i*wo+4]= matchings[i].first.layer;
      data[i*wo+5]= matchings[i].first.scale;
      data[i*wo+6]= matchings[i].first.size;
      data[i*wo+7]=(double) matchings[i].first.id;

      data[i*wo+8]=matchings[i].second.pt.x;
      data[i*wo+9]=matchings[i].second.pt.y;
      data[i*wo+10]=matchings[i].second.angle;
      data[i*wo+11]= matchings[i].second.octave;
      data[i*wo+12]= matchings[i].second.layer;
      data[i*wo+13]= matchings[i].second.scale;
      data[i*wo+14]= matchings[i].second.size;
      data[i*wo+15]= (double)matchings[i].second.id;

      data[i*wo+16]= ac_distance[i];
    }

    // Similarity matrix and vector of representants
    wo = iclasses->id_in.size();
    int wo_kp = 3;
    plhs[2] = mxCreateDoubleMatrix(wo, wo, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(wo_kp, wo, mxREAL);
    data = mxGetPr(plhs[2]);
    double *data_kp = mxGetPr(plhs[3]);
    for ( int i = 0; i < (int) wo*wo; i++)
    data[0] = 0;
    for ( int i = 0; i < (int) wo*wo_kp; i++)
    data_kp[0] = 0;

    for ( int i = 0; i < (int) matchings.size(); i++)
    {
      int ind_1, ind_2;
      bool found1=false,found2=false;
      for(ind_1=0;ind_1<iclasses->id_in.size();ind_1++)
      {
        if(iclasses->id_in[ind_1]==matchings[i].first.id)
        {
          found1 = true;
          if(data_kp[ind_1*wo_kp+0]==0)
          {
            data_kp[ind_1*wo_kp+0]=matchings[i].first.pt.x;
            data_kp[ind_1*wo_kp+1]=matchings[i].first.pt.y;
            data_kp[ind_1*wo_kp+2]= (double)matchings[i].first.id;
          }
          break;
        }
      }
      for(ind_2=0;ind_2<iclasses->id_in.size();ind_2++)
      {
        if(iclasses->id_in[ind_2]==matchings[i].second.id)
        {
          found2 = true;
          if(data_kp[ind_2*wo_kp+0]==0)
          {
            data_kp[ind_2*wo_kp+0]=matchings[i].second.pt.x;
            data_kp[ind_2*wo_kp+1]=matchings[i].second.pt.y;
            data_kp[ind_2*wo_kp+2]=(double)matchings[i].second.id;
          }
          break;
        }
      }
      if(found1&&found2)
      {
        data[ind_1*wo+ind_2] = 1 - ac_distance[i];
        data[ind_2*wo+ind_1] = 1 - ac_distance[i];
        //std::cout<<"ok"<<std::endl;
      }
      else
      {
        std::cout<<"error: found1="<<found1<<", found2="<<found2<<std::endl;
        std::cout<<matchings[i].first.id<<" and "<<matchings[i].second.id<<" are not among the following"<<std::endl;
        for(ind_1=0;ind_1<iclasses->id_in.size();ind_1++)
        std::cout<<iclasses->id_in[ind_1]<<std::endl;

        return;
      }
    }
#else
      iclasses->PrintCluster();
#endif

  }
