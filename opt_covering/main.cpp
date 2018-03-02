#ifdef _OPENMP
#include <omp.h>
#endif
#include <sstream>
#include <string>
#include <iostream>
#include <stdio.h> /* printf */
#include <math.h> /* round, floor, ceil, trunc */
#include <vector>
#include <algorithm>

#ifdef _PNG
#include "io_png/io_png.h"
#endif
using namespace std;



typedef size_t idxtype;

float epsilon, logepsilon, beta_epsilon, r,logr, beta_r, region, logregion, beta_region;
int N;

#define ABS(x)    (((x) > 0) ? (x) : (-(x)))
const float point_radius = 2.0f; // for drawing points in the output image




void drawpoint(float * data, int width, int i, int j, float radius, float value)
{
    for (int x=-radius;x<=radius;x++)
        for (int y=-radius;y<=radius;y++)
        {
            if ((y+j>0)&&(x+i>0)&&(y+j<width)&&(x+i<width)&&(sqrt(x*x+y*y)<=radius))
            {
                data[(y+j)*width + x+i] = value;
            }
        }
}


void getdisks(float t,float psi, float r, float phistep, std::vector<float>& phivec, std::vector<float>& tupvec,std::vector<float>& tlowvec)
{
    float beta = ( pow(r,2)+1 )/(2*r);
    float phi = 0;
    while(phi <2*M_PI)
    {

        float Gphi = pow(cos(psi-phi),2);
        float dis = pow(beta,2) - ( Gphi/t + t*(1-Gphi) )*( (1-Gphi)/t + t*Gphi );
        float tupper, tlower;
        if( dis>0 )
        {
            tupper = ( beta + sqrt(dis) )/( Gphi/t + (1-Gphi)*t );
            tlower = ( beta - sqrt(dis) )/( Gphi/t + (1-Gphi)*t );
            if (tupper>1)
            {
                phivec.push_back(phi);
                tupvec.push_back(tupper);
                if (tlower<1)
                    tlowvec.push_back( 1.0f);
                else
                    tlowvec.push_back(tlower);
            }
        }
        else
            if (dis==0)
            {
                tupper = beta/( Gphi/t + (1-Gphi)*t );
                if (tupper>1)
                {
                    phivec.push_back(phi);
                    tupvec.push_back(tupper);
                    tlowvec.push_back(tupper);
                }
            }
        phi += phistep;
    }
}

void drawborder(float *data, std::vector<float> phivec, std::vector<float>& tvec,float seenregion, int logseenregionpixels)
{
    int i,j,wo = logseenregionpixels*2 +1; float x,y;
    for(int n=0;n<phivec.size();n++)
    {
        x = log(tvec[n])*logseenregionpixels/log(seenregion)*cos(phivec[n]);
        y = log(tvec[n])*logseenregionpixels/log(seenregion)*sin(phivec[n]);
        i = (int)round(x+logseenregionpixels);
        j = (int)round(y+logseenregionpixels);

        //drawpoint(data, wo, i, j, point_radius, -1);
        data[j*wo+i] = -2;
    }

}

std::vector<float> scalarmult(std::vector<float> vec, float scalar)
{
    for(int n=0;n<vec.size();n++)
        vec[n] = vec[n] + scalar;
    return(vec);
}

void drawdisks(float * data, const vector<float>& tilt_element, float r,float seenregion, int logseenregionpixels)
{
    const float phidiskstep = M_PI/100000;
    std::vector<float> phivec, tupvec, tlowvec;
    phivec.clear(); tupvec.clear(); tlowvec.clear();

    int i,j,wo = logseenregionpixels*2 +1; float x,y;
    i = 0 + logseenregionpixels;
    j = 0 + logseenregionpixels;
    drawpoint(data, wo, i, j, point_radius, -1);
    getdisks(1.0f,0, r, phidiskstep, phivec, tupvec, tlowvec);
    drawborder(data, phivec, tupvec,seenregion,logseenregionpixels);

    for(int n=1;n<=N;n++)
    {
        float phistep1 = tilt_element[2*n-1];
        float t1 = tilt_element[2*n-2];
        float phi1 = 0;
        phivec.clear(); tupvec.clear(); tlowvec.clear();
        getdisks(t1,phi1, r, phidiskstep, phivec, tupvec, tlowvec);
        while(phi1<M_PI)
        {
            x = log(t1)*logseenregionpixels/log(seenregion)*cos(phi1);
            y = log(t1)*logseenregionpixels/log(seenregion)*sin(phi1);
            i = (int)round(x+logseenregionpixels);
            j = (int)round(y+logseenregionpixels);

            drawpoint(data, wo, i, j, point_radius, -1);
            drawborder(data, scalarmult(phivec,phi1), tupvec,seenregion,logseenregionpixels);
            drawborder(data, scalarmult(phivec,phi1), tlowvec,seenregion,logseenregionpixels);

            phi1 = phi1+phistep1;
        }
    }

    // plot(log(t).*cos(psi),log(t).*sin(psi),DRAW_CENTER_DISKS_OPTS,'MarkerSize',12);
    // tempball.curvsup = plot(log(pointsup).*cos(phivect),log(pointsup).*sin(phivect),'-','Color',colorvec);
}


float transition_tilt(float t, float psi1, float s, float psi2);

void write_image_covering(const vector<float>& tilt_element, float r,float region,float seenregion, int logseenregionpixels)
{

    int wo =  2*logseenregionpixels + 1;
    int ho = wo;

    float * data = new float[wo*ho];
    for (int i=0;i<wo*ho;i++)
        data[i] = 0;


    // status for each pixel
    for (int x=-logseenregionpixels;x<=logseenregionpixels;x++)
        for (int y=-logseenregionpixels;y<=logseenregionpixels;y++)
        {

            float rho = sqrt(x*x +y*y);
            float phi = atan2(y,x);
            float t = exp( (log(seenregion)*rho/(float)logseenregionpixels) );

            int i = x+logseenregionpixels;
            int j = y+logseenregionpixels;

            if (transition_tilt(t,phi,1,0)<r)
                data[j*wo+i]++;

            for(int n=1;n<=N;n++)
            {
                float phistep1 = tilt_element[2*n-1];
                if (phistep1==0)
                    cout<<phistep1<<endl;
                float t1 = tilt_element[2*n-2];
                float phi1 = 0;
                while(phi1<M_PI)
                {
                    if (transition_tilt(t,phi,t1,phi1)<r)
                        data[j*wo+i]++;
                    phi1 = phi1+phistep1;
                }
            }
        }
    drawdisks(data, tilt_element, r,seenregion, logseenregionpixels);

    //image with orientations as in MATLAB
    reverse(data,data+wo*ho);
    for (int i=0;i<wo;i++)
        reverse(data+i*wo,data+(i+1)*wo-1);


    float* red = new float[3], *green = new float[3], *gray = new float[3], *white = new float[3],*blue = new float[3],*black = new float[3];

    red[0] = 204.0f;red[1] = 1.0f; red[2] = 51.0f;
    green[0] = 1.0f;green[1] = 150.0f; green[2] = 51.0f;
    black[0] = 0.0f; black[1] = 0.0f; black[2] = 0.0f;
    blue[0] = 255*0.871f; blue[1] = 255*0.922f; blue[2] = 255*0.98f;
    gray[0] = 150.0f;gray[1] = 150.0f; gray[2] = 150.0f;
    white[0] = 250.0f;white[1] = 250.0f; white[2] = 250.0f;

    float * rgb = new float[wo*ho*3];
    for(int c=0;c<3;c++)
        for(int j = 0; j < (int) ho; j++)
            for(int i = 0; i < (int) wo; i++)
            {
                if (data[j*wo+i]>1)
                    rgb[j*wo+i+c*(wo*ho)] = blue[c];
                if (data[j*wo+i]==1)
                    rgb[j*wo+i+c*(wo*ho)] = white[c];
                if (data[j*wo+i]==0)
                    rgb[j*wo+i+c*(wo*ho)] = gray[c];

                if (data[j*wo+i]==-1)
                    rgb[j*wo+i+c*(wo*ho)] = green[c];

                if (data[j*wo+i]==-2)
                    rgb[j*wo+i+c*(wo*ho)] = red[c];
            }

    //dashed line GAMMA in black
    float phi = 0;
    while(phi<2*M_PI)
    {
        int x = (int) round( log(region)*cos(phi)*logseenregionpixels/log(seenregion) );
        int y = (int) round( log(region)*sin(phi)*logseenregionpixels/log(seenregion) );
        int i = x+logseenregionpixels;
        int j = y+logseenregionpixels;

        if((i>0)&&(j>0)&&(j<ho)&&(i<ho))
            for(int c=0;c<3;c++)
            {
                rgb[j*wo+i+c*(wo*ho)] = black[c];
            }

        phi += 2*M_PI/(logseenregionpixels);
    }
#ifdef _PNG
    write_png_f32("covering.png", rgb, wo, ho, 3);
#endif
    delete[] rgb;
}



float phi0(float& t,float& beta)
{
    return (2*acos(sqrt( ( beta - (1/pow(t,2) + pow(t,2))/2 )/( 1 - (1/pow(t,2) + pow(t,2))/2 ) )));
}


float Fcout(const vector<float>& tilt_element)
{
    float ret = 1.0f;
    for(int n=1;n<=tilt_element.size()/2;n++)
        ret += (trunc(M_PI/tilt_element[2*n-1])+1)/tilt_element[2*n-2];
    return (ret);
}


float transition_tilt(float t, float psi1, float s, float psi2)
{
    float cos_2 = pow(cos(psi1-psi2),2);
    float g = ( pow(t/s,2) + 1 )*cos_2 + ( 1/pow(s,2) + pow(t,2) )*(1-cos_2);
    float G = (s/t)*g/2;
    return( G + sqrt(pow(G,2) - 1) );
}


void intersection_points(float t, float phi_step, float beta, vector<float>& tvec, float& phi)
{
    phi = phi_step/2;

    float Gphi = pow(cos(phi),2);
    float dis = pow(beta,2) - ( Gphi/t + t*(1-Gphi) )*( (1-Gphi)/t + t*Gphi );
    tvec.clear();

    if (dis==0)
    {
        tvec.push_back( beta/( Gphi/t + (1-Gphi)*t ) );
    }
    if (dis>0)
    {
        tvec.push_back(  ( beta - sqrt(dis) )/( Gphi/t + (1-Gphi)*t )  );
        tvec.push_back(  ( beta + sqrt(dis) )/( Gphi/t + (1-Gphi)*t )  );
    }
}



bool verify_cover(const vector<float>& tilt_element, float epsilon_local, float beta_epsilon_local)
{
    float covered_portion = r;

    vector<float> tvec_inter;
    float phi_inter;
    intersection_points(tilt_element[0], tilt_element[1], beta_r, tvec_inter, phi_inter);

    if(tvec_inter.size()==0)
        return(false);

    for(int n=1;n<tilt_element.size()/2;n++)
    {
        if (*std::min_element(tvec_inter.begin(),tvec_inter.end())>covered_portion)
            return(false);

        covered_portion = *std::max_element(tvec_inter.begin(),tvec_inter.end());

        if (covered_portion>region)
            return(true);

        // look for new intersection points
        intersection_points(tilt_element[2*(n+1)-2], tilt_element[2*(n+1)-1], beta_r, tvec_inter, phi_inter);
        if(tvec_inter.size()==0)
            return(false);


        float tstart = covered_portion, tstop = *std::min_element(tvec_inter.begin(),tvec_inter.end());

        float phase1 = tilt_element[2*n-1], phase2 = tilt_element[2*(n+1)-1];
        float t1     = tilt_element[2*n-2], t2     = tilt_element[2*(n+1)-2];

        float ti = tstart*epsilon_local;
        while (ti<tstop)
        {
            float phi_epsilon=phi0(ti,beta_epsilon_local);
            float phii=0;
            while (phii<=M_PI)
            {
                /*Checking 4 neigboor disks*/

                bool covered = false;
                if (transition_tilt(ti,phii,t1,floor(phii/phase1)*phase1)<r)
                    covered = true;
                if (transition_tilt(ti,phii,t1,ceil(phii/phase1)*phase1)<r)
                    covered = true;
                if (transition_tilt(ti,phii,t2,floor(phii/phase2)*phase2)<r)
                    covered = true;
                if (transition_tilt(ti,phii,t2,ceil(phii/phase2)*phase2)<r)
                    covered = true;
                if (covered==false)
                    return(false);
                phii=phii+phi_epsilon;
            }
            ti = ti*epsilon_local;
        }
    }

    if ((*std::max_element(tvec_inter.begin(),tvec_inter.end()))>region)
        return(true);
    else
        return(false);
}


template <class T, class idxer>
class supervector
{
public:
    class supertree
    {
    public:
        supertree(const vector<supertree*>& vec)
        {
            right = vec;
        }
        supertree(T val)
        {
            right.clear();
            value = val;
        }

        vector<supertree*> right; // way to the end... if equal to 0 this is the end
        T value;
    };

    //  supervector(const supervector<T,idxer>& vec)
    //  {
    //      ordered = vec.ordered;
    //      veclength = vec.veclength;
    //      total = vec.total;
    //      sizes.resize(vec.sizes.size());
    //      copy(vec.sizes.begin(), vec.sizes.end(), sizes.begin());
    //      Cinc.resize(vec.Cinc.size());
    //      copy(vec.Cinc.begin(), vec.Cinc.end(), Cinc.begin());
    //      if (vec.single_vector.empty())
    //      {
    //          set_of_vectors.resize(vec.set_of_vectors.size());
    //          for(int n=0;n<vec.set_of_vectors.size();n++)
    //          {
    //              set_of_vectors[n] = *(new supervector<T,idxer>(vec.set_of_vectors[n]));
    //          }
    //      }

    //  }

    supervector(const vector<T>& vec)
    {
        single_vector = vec;
        total = single_vector.size();
        set_of_vectors.clear();
        sizes.resize(1);
        Cinc.resize(1);
        Cinc[0] = 1;
        sizes[0] = single_vector.size();
        veclength = 1;
    }
    supervector(const vector< supervector<T,idxer>* >& supervec)
    {
        single_vector.clear();
        sizes.resize(supervec.size());
        Cinc.resize(supervec.size());
        total = 1;
        veclength = 0;
        for(idxer i=0;i<(idxer)supervec.size();i++)
        {
            total *= supervec[i]->total;
            sizes[i] = supervec[i]->total;
            set_of_vectors.push_back(*supervec[i]);
            veclength += supervec[i]->veclength;
        }

        Cinc[0] = 1;
        for(idxer j=1;j<(idxer)sizes.size();j++)
            Cinc[j] = Cinc[j-1]*sizes[j-1];
    }


    /**
     * @brief Translates a 1-dimensional index into an N-dimensional vector of indices
     * @param ind 1-dimensional index
     * @return Vector of indices \f$(i_1,...,i_N) \in \left\lbrace 1,C_1\rigth\rbrace \times \cdots \times \left\lbrace 1,C_N\rigth\rbrace\f$
     */
    vector<idxer> ind2prod(idxer ind)
    {
#ifdef _saveflag
        if (!((ind<total)&(ind>=0)))
            cerr<<"Error in this call -> ind2prod("<<ind<<")"<<" // not on allowed interval"<<endl;
#endif
        idxer N = Cinc.size();
        vector<idxer> prod(N);
        for(idxer j=1;j<=N;j++) //idxer might be only positive numbers
        {
            ldiv_t dv;// it might need dv{}
            dv = std::ldiv(ind,Cinc[N-j]);
            prod[N-j] = dv.quot;
            ind = dv.rem;
        }
        return(prod);
    }

    /**
     * @brief Translates an N-dimensional vector of indices into a 1-dimensional index
     * @param prod Vector of indices \f$(i_1,...,i_N) \in \left\lbrace 1,C_1\rigth\rbrace \times \cdots \times \left\lbrace 1,C_N\rigth\rbrace\f$
     * @return Index corresponding to the vector of indices
     */
    idxer prod2ind(const vector<idxer>& prod)
    {
#ifdef _saveflag
        if (prod.size()!=veclength)
            cerr<<"Error in this call -> prod2ind(prodvec) // prodvec.size()!=veclength="<<veclength<<endl;
        for(idxer i=0;i<veclength;i++)
            if (!((prod[i]<sizes[i])&(prod[i]>=0)))
                cerr<<"Error in this call -> prod2ind(prodvec) // prodvec["<<i<<"] not on allowed interval"<<endl;
#endif
        idxer N = (idxer) Cinc.size(), ind =0;
        for(idxer j=0;j<N;j++)
            ind += Cinc[j]*prod[j];

#ifdef _saveflag
        if (!((ind<total)&(ind>=0)))
            cerr<<"Error at the end of this call -> prod2ind(vector) // ind="<<ind<<" not on allowed interval"<<endl;
#endif
        return(ind);
    }

    supertree* get_from_ind(idxer i)
    {
        vector<idxer> prod = ind2prod(i);

        if (single_vector.empty())
        {
            vector<supertree*> trees(prod.size());
            for(idxer j = 0; j<prod.size(); j++)
                trees[j] = set_of_vectors[j].get_from_ind(prod[j]);
            supertree* linkingnode = new supertree(trees);
            return(linkingnode);
        }
        else
        {
            supertree* leaf = new supertree(single_vector[prod[0]]);
            return(leaf);
        }

    }

    T get_from_prod(const vector<idxer>& prod)
    {
        vector<T> values(veclength);
        if (single_vector.empty())
        {
            idxer valcount = 1;
            for (idxer j =0; j<prod.size();j++)
            {
                vector<T> subvalues = set_of_vectors[j].get_from_ind(prod[j]);
                for (idxer k =0; k<subvalues.size();k++)
                    values[valcount++] = subvalues[k];
            }
        }
        else
        {
            for (idxer j =0; j<prod.size();j++)
                values[j] = single_vector[prod[j]];
        }
        return(values);
    }

    idxer get_total()
    {
        return(total);
    }

    void vector_from_supertree(supertree* tree, vector<T>& values)
    {
        values.resize(veclength);
        idxer pos = 0;
        read_tree(tree,values,pos);
        delete tree;
    }

    void print_all_prod_values()
    {
        for(idxer j=0;j<total;j++)
        {
            vector<T> ret;
            get_vector_from_ind(j,ret);
            //vector_from_supertree(get_from_ind(j),ret);
            cout<<"("<<ret[0];
            for (idxer i=1;i<ret.size();i++)
                cout<<", "<<ret[i];
            cout<<")"<<endl;
        }
    }

    void print_prod_values(idxer j)
    {
        vector<T> ret;
        get_vector_from_ind(j,ret);
        //vector_from_supertree(get_from_ind(j),ret);
        cout<<"("<<ret[0];
        for (idxer i=1;i<ret.size();i++)
            cout<<", "<<ret[i];
        cout<<")"<<endl;
    }

    void get_vector_from_ind(idxer j,vector<T>& ret)
    {
        vector_from_supertree(get_from_ind(j),ret);
    }

    idxer get_veclength()
    {
        return(veclength);
    }

    T get_value_from_subvector(idxer i, idxer vec_ind)
    {
        vector<T> values(set_of_vectors[vec_ind].veclength);
        idxer pos = 0;
        read_tree(set_of_vectors[vec_ind].get_from_ind(i), values, pos);
        return(values[0]);
    }

    idxer get_size_from_subvector(idxer vec_ind)
    {
        return(sizes[vec_ind]);
    }

    idxer find_nearest_to_subvector(const vector<T>& values, idxer vec_ind)
    {
        return(set_of_vectors[vec_ind].find_nearest_to_vector(values));
    }

    idxer find_nearest_to_vector(const vector<T>& values)
    {
        if (single_vector.empty())
        {
            vector<idxer> prod(veclength);
            for (idxer j=0; j<veclength;j++)
            {
                vector<T> subvalue(1);
                subvalue[0] = values[j];
                prod[j] = set_of_vectors[j].find_nearest_to_vector(subvalue);
            }
            return(prod2ind(prod));
        }
        else
        {
            if (ordered)
                return( find_fast(values[0],0,single_vector.size()-1) );
            else
                return(find(values[0]));

        }
    }

private:
    idxer total;
    bool ordered = true; // if single_vector this property can be exploited for fast computation of an element

    // one of these two is always empty
    vector<T> single_vector;
    vector<supervector> set_of_vectors;

    //
    idxer veclength; // length for return values

    vector<idxer> sizes;

    /**
     * @brief Cinc A precomputed vector that is equal to \f$ Cinc = \left\lbrace 1, C_1, C_1*C_2, C_1*C_2*\cdots * C_{N-1} \rigth\rbrace\f$
     */
    vector<idxer> Cinc;

    idxer find_fast(T value, idxer i_min, idxer i_max)
    {

        if ((i_max-i_min)<=3)
        {
            idxer argmin = i_min;
            T min_dist = ABS(single_vector[i_min]-value), current_dist;
            for (idxer j=i_min;j<=i_max;j++)
            {
                current_dist = ABS(single_vector[j]-value);
                if (min_dist>current_dist)
                {
                    min_dist = current_dist;
                    argmin = j;
                }
            }
            return(argmin);
        }

        idxer middle = floor( (i_max+i_min)/2);

        if (single_vector[middle]>value)
            return(find_fast(value,i_min,middle+1));
        else
            return(find_fast(value,middle-1,i_max));
    }

    idxer find(T value)
    {
        idxer argmin = 0;
        T min_dist = ABS(single_vector[0]-value), current_dist;
        for (idxer j=1;j<single_vector.size();j++)
        {
            current_dist = ABS(single_vector[j]-value);
            if (min_dist>current_dist)
            {
                min_dist = current_dist;
                argmin = j;
            }
        }
        return(argmin);
    }


    //NOTE: it deletes branches after reading
    void read_tree(supertree* tree, vector<T>& values, idxer& pos)
    {
        if (tree->right.empty())
        {
            values[pos++] = tree->value;
        }
        else
            for(idxer i=0;i<(idxer)tree->right.size();i++)
            {
                read_tree(tree->right[i],values,pos);
                delete tree->right[i];
            }
    }
};


//n=1,2,3,...,N
void gen_next_phases_n(vector<idxtype>& feasiblesets_index, supervector<float,idxtype>& tilts_set, vector<idxtype>& tilts_index, int n)
{
    float tn = tilts_set.get_value_from_subvector( tilts_index[2*n-2], 2*n-2 );
    float tmax = tilts_set.get_value_from_subvector( tilts_set.get_size_from_subvector(2*n-2)-1, 2*n-2 );

    //float phase_epsilon = phi0(tn,beta_epsilon);
    //float phase_stop =2*phi0(tn,beta_r);// 2*phi0(tn,beta_r); // this is an heuristic, in theory it should be \pi

    float phase_epsilon = phi0(region,beta_epsilon);
    float phase_stop =2*phi0(tn,beta_r);// 2*phi0(tn,beta_r); // this is an heuristic, in theory it should be \pi


    if (phase_stop>M_PI)
        phase_stop = M_PI;

    //float tmaxphase_epsilon = phi0(tmax,beta_epsilon);
    vector<float> singlelem(1);


    if (n==N)
    {
        for(float phase = phase_epsilon; phase< phase_stop; phase+=phase_epsilon )
        {
            //tilts_index[2*n-1] = (idxtype) round(phase/tmaxphase_epsilon); //current index
            singlelem[0]=phase;
            tilts_index[2*n-1] = tilts_set.find_nearest_to_subvector(singlelem, 2*n-1); //current index
            feasiblesets_index.push_back(tilts_set.prod2ind(tilts_index));
        }
    }
    else
    {
        for(float phase = phase_epsilon; phase< phase_stop; phase+=phase_epsilon )
        {
            //tilts_index[2*n-1] = (idxtype) round(phase/tmaxphase_epsilon); // current index
            singlelem[0]=phase;
            tilts_index[2*n-1] = tilts_set.find_nearest_to_subvector(singlelem, 2*n-1); //current index
            gen_next_phases_n(feasiblesets_index,tilts_set,tilts_index,n+1);
        }
    }

}


//n=1,2,3,...,N
void gen_next_tilts_n(vector<idxtype>& feasiblesets_index, supervector<float,idxtype>& tilts_set, vector<idxtype>& tilts_index, int n)
{

    float tp;
    idxtype t_min, t_max;
    vector<float> singlelem(1);
    if (n>1)
    {
        tp= tilts_set.get_value_from_subvector( tilts_index[2*n-4], 2*n-4 );
        //t_min = floor(log( tp/pow(r,n-1) )/logepsilon);
        //t_max = floor(log( tp/pow(r,n-2) )/logepsilon);
        singlelem[0] = tp*r;
        t_min =  tilts_set.find_nearest_to_subvector(singlelem, 2*n-2);
        singlelem[0] = tp*pow(r,2);
        t_max =  tilts_set.find_nearest_to_subvector(singlelem, 2*n-2);
    }
    else
    {
        tp =1;
        t_min = 0;
        //t_max = floor( logr/logepsilon );
        singlelem[0] = tp*pow(r,2);
        t_max =  tilts_set.find_nearest_to_subvector(singlelem, 2*n-2);

    }

    if (n==N)
    {
        for (idxtype i=t_min; i<=t_max;i++)
        {
            tilts_index[2*n-2] = i;
            gen_next_phases_n(feasiblesets_index,tilts_set,tilts_index,1);
        }
    }
    else
    {
        for (idxtype i=t_min; i<=t_max;i++)
        {
            tilts_index[2*n-2] = i;
            gen_next_tilts_n(feasiblesets_index,tilts_set,tilts_index,n+1);
        }
    }
}

#include <algorithm>
#include <ctime>;
#include <cstdlib>
int myrandom (int i)
{
    return std::rand()%i;
}



double search_for_minimum(vector<float>& argmin,vector<idxtype>& feasiblesets_index, supervector<float,idxtype>& tilts_set, float fixed_Fcout)
{
    double start_time=0, run_time=0;
#ifdef _OPENMP
    start_time = omp_get_wtime();
    int maxthreads=omp_get_max_threads();
    omp_set_num_threads(maxthreads);
#endif

    float min_Fcout = fixed_Fcout;
    float internal_min_Fcout = min_Fcout;
    vector<float> tilt_element;
    //supervector<float,idxtype> tilts_set = tilts_set2;
    std::srand ( unsigned ( std::time(0) ) );
    std::random_shuffle( feasiblesets_index.begin(), feasiblesets_index.end(), myrandom);
#pragma omp parallel for shared(feasiblesets_index,min_Fcout,argmin) firstprivate(internal_min_Fcout,epsilon, beta_epsilon, r, beta_r, region,tilt_element) schedule(static)
    for(idxtype i=0;i<feasiblesets_index.size();i++)
    {
        tilts_set.get_vector_from_ind(feasiblesets_index[i],tilt_element);

        float thisFcout = Fcout(tilt_element);
        //cout<<thisFcout<<endl;
        if (internal_min_Fcout>thisFcout)
        {
#pragma omp atomic read
            internal_min_Fcout = min_Fcout;
            // verify cover only if really needed
            bool iscover = (internal_min_Fcout>thisFcout)&&(verify_cover(tilt_element, epsilon, beta_epsilon));
            if (iscover)
            {
#pragma omp critical
                {
                    if (min_Fcout>thisFcout)
                    {

                        min_Fcout = thisFcout;
                        argmin = tilt_element;
                    }
                }
            }
        }
    }

#ifdef _OPENMP
    run_time = omp_get_wtime() - start_time;
    return(run_time);
#endif
}



int main (int argc, char **argv)
{

    //INPUTS

    // number of concentric groups of tilts
    int N_local=2;

    // in order to uniformly discretize the space of tilts, i.e. distance between nearest elements
    float epsilon_local = 1.00854f;
    epsilon_local = 1.04f;
    if (argc>3)
    {
        epsilon_local =  atof(argv[3]);
        N_local = atoi(argv[4]);
    }
    float r_local = atof(argv[1]);
    float region_local = atof(argv[2]);


    // Setting up search parameters // TO BE USED GLOBALLY
    N=N_local;
    epsilon = epsilon_local;
    r = r_local;
    region = region_local;
    logepsilon = log(epsilon);
    logr = log(r);
    beta_epsilon = ( pow(epsilon,2)+1 )/(2*epsilon);
    beta_r = ( pow(r,2)+1 )/(2*r);
    logregion = log(region);
    beta_region = ( pow(region,2)+1 )/(2*region);



    // Creating data structure for handling discrete sets like [r^n,r^2n]_\epsilon
    vector< vector<float> > veci(2*N);
    float temp = logr/logepsilon;
    vector< supervector<float,idxtype>* > prodsvec(2*N);
    for(int n=1;n<=N;n++)
    {
        idxtype i_max;
        i_max = floor(n*temp);
        for (idxtype i=0;i<=i_max;i++)
            veci[2*n-2].push_back(pow(r,n)*pow(pow(epsilon,pow(sqrt(2),N-1)*n),i)); // varing steps like pow(epsilon,n)... lower tilts have greater costs!
        //cout<<"last = "<<veci[n-1][veci[n-1].size()-1] <<",  r^2n = " <<pow(r,2*n)<<",  last*epsilon = "<<veci[n-1][veci[n-1].size()-1]*epsilon<<endl;
        //cout<<veci[n-1].size()<<endl;
        prodsvec[2*n-2] = new supervector<float,idxtype>(veci[2*n-2]);

        //float lasttilt = veci[2*n-2][veci[2*n-2].size()-1];
        //float lastphi_epsilon = phi0(lasttilt,beta_epsilon);
        float lastphi_epsilon = phi0(region,beta_epsilon);
        idxtype max_phases = floor(M_PI/lastphi_epsilon);
        //cout<<max_phases<<endl;
        for (idxtype i=0;i<=max_phases;i++)
            veci[2*n-1].push_back(lastphi_epsilon*i);

        prodsvec[2*n-1] = new supervector<float,idxtype>(veci[2*n-1]);
    }


    supervector<float,idxtype> tilts_set(prodsvec);
    //cout<<tilts_set.get_veclength()<<" -> "<<tilts_set.get_total()<<endl;
    //tilts_set.print_prod_values(100);

    vector<idxtype> feasiblesets_index;
    vector<idxtype> tilts_index(2*N);

    gen_next_tilts_n(feasiblesets_index,tilts_set,tilts_index,1);


    cout<<"\nParameters:"<<endl;
    cout<<" r = "<<r<<", region = "<<region<<", epsilon = "<<epsilon<<", N = "<<N<<endl;
    cout<<" "<<feasiblesets_index.size()<<" feasible sets to test"<<endl<<endl;

    // SEARCH FOR MINIMUM
    cout<<" --> Starting a global search <--"<<endl;
    vector<float> argmin;
    double run_time = search_for_minimum( argmin, feasiblesets_index ,tilts_set,10000.0f );

    // Output
    cout<<" Argmin -> ";
    tilts_set.print_prod_values(tilts_set.find_nearest_to_vector(argmin));
    cout <<" Area Ratio = "<<Fcout(argmin)<<endl;
    printf("\n Work done in %lf seconds\n ",run_time);


    //    cout<<" Argmin (in MATLAB format) = (";
    //    for (int n=1;n<=argmin.size()/2;n++)
    //        cout<<" t"<<n<<"="<<argmin[2*n-2]<<"; phi"<<n<<"="<<argmin[2*n-1]<<";";
    //    cout<<" region="<<region<<"; radius="<<r<<"; )"<<endl;


    //  REFINED SEARCH
    cout<<"\n\n --> Refining argmin <--"<<endl;
    int numtilts = 10;
    run_time =0;
    cout <<"0. Area Ratio = "<<Fcout(argmin)<<" / ";
    tilts_set.print_prod_values(tilts_set.find_nearest_to_vector(argmin));
    idxtype rayon = floor((numtilts/(2*(N-1))));
    for(int iter=1;iter<=5;iter++)
    {
        vector< vector<float> > veci_refined(2*N);
        vector< supervector<float,idxtype>* > prodsvec_refined(2*N);
        //tilts_set.print_prod_values(tilts_set.find_nearest_to_vector(argmin));
        //vector<idxtype> argmin_prod_idx = tilts_set.ind2prod(tilts_set.find_nearest_to_vector(argmin));


        for(int n=1;n<=N;n++)
        {
            vector<float> value(1);
            value[0] = argmin[2*n-2];
            idxtype n_t = tilts_set.find_nearest_to_subvector(value,2*n-2);
            value[0] = argmin[2*n-1];
            idxtype phase_t = tilts_set.find_nearest_to_subvector(value,2*n-1);

            float tstart, tstop, phase_start, phase_stop;

            if (n_t>rayon)
                tstart = tilts_set.get_value_from_subvector(n_t-rayon,2*n-2);
            else
                tstart = tilts_set.get_value_from_subvector(0,2*n-2);

            if (n_t+rayon>=tilts_set.get_size_from_subvector(2*n-2))
                tstop = tilts_set.get_value_from_subvector(tilts_set.get_size_from_subvector(2*n-2)-1,2*n-2);
            else
                tstop = tilts_set.get_value_from_subvector(n_t+rayon,2*n-2);

            float epsilon_t = exp(log(tstop/tstart)/numtilts);//1.01f
            veci_refined[2*n-2].push_back(argmin[2*n-2]);
            for (float ti=argmin[2*n-2]/epsilon_t; ti>=tstart; ti/=epsilon_t)
                veci_refined[2*n-2].push_back(ti);
            for (float ti=argmin[2*n-2]*epsilon_t; ti<=tstop; ti*=epsilon_t)
                veci_refined[2*n-2].push_back(ti);

            std::sort(veci_refined[2*n-2].begin(), veci_refined[2*n-2].end());

            //cout<<"     "<<tstart<<"<"<<argmin[2*n-2]<<"<"<<tstop<<endl;

            prodsvec_refined[2*n-2] = new supervector<float,idxtype>(veci_refined[2*n-2]);


            if (phase_t>rayon)
                phase_start = tilts_set.get_value_from_subvector(phase_t-rayon,2*n-1);
            else
                phase_start = tilts_set.get_value_from_subvector(0,2*n-1);

            if (phase_t+rayon>=tilts_set.get_size_from_subvector(2*n-1))
                phase_stop = tilts_set.get_value_from_subvector(tilts_set.get_size_from_subvector(2*n-1)-1,2*n-1);
            else
                phase_stop = tilts_set.get_value_from_subvector(phase_t+rayon,2*n-1);

            float phase_epsilon = (phase_stop-phase_start)/numtilts;
            veci_refined[2*n-1].push_back(argmin[2*n-1]);
            for (float phasei=argmin[2*n-1]-phase_epsilon; phasei>=phase_start; phasei-=phase_epsilon)
                veci_refined[2*n-1].push_back(phasei);
            for (float phasei=argmin[2*n-1]+phase_epsilon; phasei<=phase_stop; phasei+=phase_epsilon)
                veci_refined[2*n-1].push_back(phasei);

            std::sort(veci_refined[2*n-1].begin(), veci_refined[2*n-1].end());
            //cout<<"     "<<phase_start<<"<"<<argmin[2*n-1]<<"<"<<phase_stop<<endl;
            prodsvec_refined[2*n-1] = new supervector<float,idxtype>(veci_refined[2*n-1]);
        }


        supervector<float,idxtype> tilts_set_refined(prodsvec_refined);
        vector<idxtype> feasiblesets_index_refined(tilts_set_refined.get_total());
        for (int n=0;n<tilts_set_refined.get_total();n++)
            feasiblesets_index_refined[n]=n;

        //tilts_set_refined.print_prod_values(tilts_set_refined.find_nearest_to_vector(argmin));
        //cout<<"searching in "<<tilts_set_refined.get_total()<<endl;
        vector<float> argmin_refined;
        if (feasiblesets_index_refined.size()!=0)
            run_time += search_for_minimum( argmin_refined, feasiblesets_index_refined ,tilts_set_refined, Fcout(argmin) );
        else
            break;

        //cout<<verify_cover(argmin,epsilon,beta_epsilon)<<" - "<<verify_cover(argmin_refined,epsilon,beta_epsilon)<<endl;

        if (argmin_refined.empty())
        {
            argmin_refined = argmin;
            tilts_set_refined = tilts_set;
        }
        else
            argmin = argmin_refined;

        cout <<iter<<". Area Ratio = "<<Fcout(argmin)<<" / ";
        tilts_set_refined.print_prod_values(tilts_set_refined.find_nearest_to_vector(argmin));
        tilts_set = tilts_set_refined;
        rayon = ceil( (numtilts+iter/(2*(N-1)))/iter);
        //cout<<rayon<<endl;
    }


    printf("\n Work done in %lf seconds\n ",run_time);
    cout<<"\n\n Argmin (in MATLAB format) = (";
    for (int n=1;n<=argmin.size()/2;n++)
        cout<<" t"<<n<<"="<<argmin[2*n-2]<<"; phi"<<n<<"="<<argmin[2*n-1]<<";";
    cout<<" region="<<region<<"; radius="<<r<<"; )"<<endl;


    cout<<"\n\n Argmin (in C++ format) = (";
    for (int n=1;n<=argmin.size()/2;n++)
        cout<<" tilts["<<n-1<<"]="<<argmin[2*n-2]<<"; phase["<<n-1<<"]="<<argmin[2*n-1]<<";";
    cout<<" region="<<region<<"; radius="<<r<<"; )"<<endl;


    cout<<"\n\n article table info \n";
    printf("$%.0f^\\circ$ & $%.0f^\\circ$ & $%.0f^\\circ$ & $%.3f$",180*acos(1/r)/M_PI, 180*acos(1/region)/M_PI, 180*acos(r/pow(region,2))/M_PI,Fcout(argmin));
    for (int n=1;n<=argmin.size()/2;n++)
        cout<<" & $"<<argmin[2*n-2]<<"$ & $"<<argmin[2*n-1]<<"$";
    cout<<"\\tabularnewline"<<endl;
    cout<<"\\hline"<<endl;

#ifdef _PNG
    write_image_covering(argmin, r, region, argmin[2*N-2]*r, 400);
    cout<<"\n Image covering.png created...\n";
#endif

    cout<<"Done..."<<endl;
}
