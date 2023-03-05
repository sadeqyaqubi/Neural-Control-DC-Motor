#include "mex.h"
#include <stdio.h>
#include <cmath>
#include <vector>
#include <process.h>
#include <cstdlib>
#define pi 3.1415927
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs [])

{
    double k1 = mxGetScalar(prhs [0]);
    double k2 = mxGetScalar(prhs [1]);
    double j1 = mxGetScalar(prhs [2]);
    double j2 = mxGetScalar(prhs [3]);
    double l1 = mxGetScalar(prhs [4]);
    double l2 = mxGetScalar(prhs [5]);
    double r1 = mxGetScalar(prhs [6]);
    double r2 = mxGetScalar(prhs [7]);
    double T = mxGetScalar(prhs [8]);
    double tf = mxGetScalar(prhs [9]);
    
    int stp = ceil(tf/T);;
    
    plhs[0] = mxCreateDoubleMatrix(stp, 13, mxREAL);
    double *out = mxGetPr( plhs [0]);
    
    double *t = out;
    double *x1 = t+stp;
    double *x2 = x1+stp;
    double *x3 = x2+stp;
    double *x4 = x3+stp;
    double *u1 = x4+stp;
    double *u2 = u1+stp;
    double *x1f = u2+stp;
    double *x2f = x1f+stp;
    double *x3f = x2f+stp;
    double *x4f = x3f+stp;
    double *debpar1 = x4f+stp;
    double *debpar2 = debpar1+stp;
    
//
    
    double x1_k;
    double x2_k;
    double x3_k;
    double x4_k;
    double u1_k;
    double u2_k;
    double x1_k1;
    double x2_k1;
    double x3_k1;
    double x4_k1;
    double u1_k1;
    double u2_k1;
    
    double a1_k;
    double a2_k;
    double alpha_k;
    double omegak1;
    double theta_k2;
    
    double b11;
    double b12;
    double b21;
    double b22;
    
    double phih1;
    double phih2;
    double phil1;
    double phil2;
    double phim1;
    double phim2;
    
    double detb;
    
    double lambda = 0.25;
    double alpha = 0.9;
    
    double x1hat;
    double x3hat;
    double x1hat1;
    double x3hat1;
    
    double xd3_0;
    double xd3_1;
    double xd2_1;
    double xd2_2;
    
    double ax3 = 2;
    double wx3 = 1;
    double ax2 = 2;
    double wx2 = 1;
    
    double stilde1;
    double stilde2;
    
    double a11;
    double a21;
    
    double randmin = 1.0;
    double randmax = 100.0;
    
    double m_mat11_1 = ((rand()/double (RAND_MAX))*randmax)+randmin ;
    double m_mat21_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat31_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat41_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat51_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat61_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat12_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat22_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat32_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat42_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat52_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat62_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    
    double sigma_mat11_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat21_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat31_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat41_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat51_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat61_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat12_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat22_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat32_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat42_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat52_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat62_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    
    double a_mat11_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat21_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat31_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat41_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat51_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat61_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat12_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat22_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat32_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat42_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat52_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat62_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    
    double b_mat11_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat21_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat31_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat41_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat51_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat61_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat12_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat22_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat32_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat42_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat52_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat62_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    
    double weight_mat1_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double weight_mat2_1 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    
    //
    double m_mat11_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat21_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat31_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat41_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat51_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat61_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat12_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat22_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat32_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat42_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat52_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat62_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    
    double sigma_mat11_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat21_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat31_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat41_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat51_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat61_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat12_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat22_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat32_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat42_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat52_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat62_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    
    double a_mat11_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat21_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat31_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat41_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat51_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat61_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat12_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat22_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat32_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat42_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat52_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat62_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    
    double b_mat11_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat21_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat31_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat41_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat51_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat61_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat12_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat22_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat32_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat42_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat52_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat62_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    
    double weight_mat1_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double weight_mat2_2 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
///
    double m_mat11_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat21_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat31_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat41_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat51_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat61_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat12_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat22_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat32_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat42_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat52_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat62_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    
    double sigma_mat11_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat21_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat31_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat41_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat51_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat61_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat12_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat22_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat32_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat42_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat52_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat62_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    
    double a_mat11_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat21_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat31_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat41_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat51_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat61_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat12_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat22_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat32_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat42_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat52_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat62_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    
    double b_mat11_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat21_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat31_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat41_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat51_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat61_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat12_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat22_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat32_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat42_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat52_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat62_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    
    double weight_mat1_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double weight_mat2_3 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
////
    
    double m_mat11_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat21_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat31_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat41_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat51_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat61_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat12_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat22_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat32_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat42_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat52_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double m_mat62_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    
    double sigma_mat11_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat21_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat31_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat41_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat51_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat61_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat12_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat22_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat32_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat42_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat52_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double sigma_mat62_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    
    double a_mat11_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat21_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat31_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat41_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat51_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat61_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat12_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat22_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat32_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat42_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat52_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double a_mat62_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    
    double b_mat11_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat21_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat31_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat41_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat51_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat61_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat12_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat22_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat32_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat42_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat52_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double b_mat62_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    
    double weight_mat1_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    double weight_mat2_4 = ((rand()/double (RAND_MAX))*randmax)+randmin  ;
    
    double o2mat11_1 ;
    double o2mat21_1 ;
    double o2mat31_1 ;
    double o2mat41_1 ;
    double o2mat51_1 ;
    double o2mat61_1;
    
    double o2mat12_1 ;
    double o2mat22_1;
    double o2mat32_1 ;
    double o2mat42_1;
    double o2mat52_1 ;
    double o2mat62_1 ;
    
    double o2mat11_2;
    double o2mat21_2 ;
    double o2mat31_2;
    double o2mat41_2;
    double o2mat51_2;
    double o2mat61_2;
    
    double o2mat12_2;
    double o2mat22_2;
    double o2mat32_2;
    double o2mat42_2;
    double o2mat52_2;
    double o2mat62_2;
    
    double o2mat11_3;
    double o2mat21_3;
    double o2mat31_3;
    double o2mat41_3;
    double o2mat51_3;
    double o2mat61_3;
    
    double o2mat12_3;
    double o2mat22_3;
    double o2mat32_3;
    double o2mat42_3;
    double o2mat52_3;
    double o2mat62_3;
    
    double o2mat11_4;
    double o2mat21_4;
    double o2mat31_4;
    double o2mat41_4;
    double o2mat51_4;
    double o2mat61_4;
    
    double o2mat12_4;
    double o2mat22_4;
    double o2mat32_4;
    double o2mat42_4;
    double o2mat52_4;
    double o2mat62_4;
    
    double o3mat1_1;
    double o3mat2_1;
    
    double o3mat1_2;
    double o3mat2_2;
    
    double o3mat1_3;
    double o3mat2_3;
    
    double o3mat1_4;
    double o3mat2_4;
    
    double o3sum_1 ;
    double o3sum_2;
    double o3sum_3;
    double o3sum_4;
    
    double o4mat1_1;
    double o4mat2_1;
    
    double o4mat1_2;
    double o4mat2_2;
    
    double o4mat1_3;
    double o4mat2_3;
    
    double o4mat1_4;
    double o4mat2_4;
    
    // WAVELET
    
    double nu_mat11_1;
    double nu_mat21_1;
    double nu_mat31_1;
    double nu_mat41_1;
    double nu_mat51_1;
    double nu_mat61_1;
    
    double nu_mat12_1;
    double nu_mat22_1;
    double nu_mat32_1;
    double nu_mat42_1;
    double nu_mat52_1;
    double nu_mat62_1;
    
    //
    double nu_mat11_2;
    double nu_mat21_2;
    double nu_mat31_2;
    double nu_mat41_2;
    double nu_mat51_2;
    double nu_mat61_2;
    
    double nu_mat12_2;
    double nu_mat22_2;
    double nu_mat32_2;
    double nu_mat42_2;
    double nu_mat52_2;
    double nu_mat62_2;
    ///
    double nu_mat11_3;
    double nu_mat21_3;
    double nu_mat31_3;
    double nu_mat41_3;
    double nu_mat51_3;
    double nu_mat61_3;
    
    double nu_mat12_3;
    double nu_mat22_3;
    double nu_mat32_3;
    double nu_mat42_3;
    double nu_mat52_3;
    double nu_mat62_3;
    ////
    double nu_mat11_4;
    double nu_mat21_4;
    double nu_mat31_4;
    double nu_mat41_4;
    double nu_mat51_4;
    double nu_mat61_4;
    
    double nu_mat12_4;
    double nu_mat22_4;
    double nu_mat32_4;
    double nu_mat42_4;
    double nu_mat52_4;
    double nu_mat62_4;
    
    
    double phi_mat11_1;
    double phi_mat21_1;
    double phi_mat31_1;
    double phi_mat41_1;
    double phi_mat51_1;
    double phi_mat61_1;
    
    double phi_mat12_1;
    double phi_mat22_1;
    double phi_mat32_1;
    double phi_mat42_1;
    double phi_mat52_1;
    double phi_mat62_1;
    
    //
    double phi_mat11_2;
    double phi_mat21_2;
    double phi_mat31_2;
    double phi_mat41_2;
    double phi_mat51_2;
    double phi_mat61_2;
    
    double phi_mat12_2;
    double phi_mat22_2;
    double phi_mat32_2;
    double phi_mat42_2;
    double phi_mat52_2;
    double phi_mat62_2;
    ///
    double phi_mat11_3;
    double phi_mat21_3;
    double phi_mat31_3;
    double phi_mat41_3;
    double phi_mat51_3;
    double phi_mat61_3;
    
    double phi_mat12_3;
    double phi_mat22_3;
    double phi_mat32_3;
    double phi_mat42_3;
    double phi_mat52_3;
    double phi_mat62_3;
    ////
    double phi_mat11_4;
    double phi_mat21_4;
    double phi_mat31_4;
    double phi_mat41_4;
    double phi_mat51_4;
    double phi_mat61_4;
    
    double phi_mat12_4;
    double phi_mat22_4;
    double phi_mat32_4;
    double phi_mat42_4;
    double phi_mat52_4;
    double phi_mat62_4;
    
    double omega_vec1_1;
    double omega_vec2_1;
    
    double omega_vec1_2;
    double omega_vec2_2;
    
    double omega_vec1_3;
    double omega_vec2_3;
    
    double omega_vec1_4;
    double omega_vec2_4;
    
    double xf_1;
    double xf_2;
    double xf_3;
    double xf_4;
    
    double o5mat1_1;
    double o5mat2_1;
    
    double o5mat1_2;
    double o5mat2_2;
    
    double o5mat1_3;
    double o5mat2_3;
    
    double o5mat1_4;
    double o5mat2_4;
    
    // FWNN UPDATE
    
    double pym11_1;
    double pym21_1;
    double pym31_1;
    double pym41_1;
    double pym51_1;
    double pym61_1;
    
    double pym12_1;
    double pym22_1;
    double pym32_1;
    double pym42_1;
    double pym52_1;
    double pym62_1;
    
//
    
    double pym11_2;
    double pym21_2;
    double pym31_2;
    double pym41_2;
    double pym51_2;
    double pym61_2;
    
    double pym12_2;
    double pym22_2;
    double pym32_2;
    double pym42_2;
    double pym52_2;
    double pym62_2;
    
///
    
    double pym11_3;
    double pym21_3;
    double pym31_3;
    double pym41_3;
    double pym51_3 ;
    double pym61_3 ;
    
    double pym12_3;
    double pym22_3;
    double pym32_3 ;
    double pym42_3;
    double pym52_3 ;
    double pym62_3 ;
    
    double pym11_4 ;
    double pym21_4 ;
    double pym31_4 ;
    double pym41_4 ;
    double pym51_4;
    double pym61_4 ;
    
    double pym12_4 ;
    double pym22_4 ;
    double pym32_4 ;
    double pym42_4;
    double pym52_4 ;
    double pym62_4 ;
    
    
    ////////////////////
    
    double pys11_1;
    double pys21_1;
    double pys31_1;
    double pys41_1;
    double pys51_1;
    double pys61_1;
    
    double pys12_1;
    double pys22_1;
    double pys32_1;
    double pys42_1;
    double pys52_1;
    double pys62_1;
    
//
    
    double pys11_2;
    double pys21_2;
    double pys31_2;
    double pys41_2;
    double pys51_2;
    double pys61_2;
    
    double pys12_2;
    double pys22_2;
    double pys32_2;
    double pys42_2;
    double pys52_2;
    double pys62_2;
    
///
    
    double pys11_3;
    double pys21_3;
    double pys31_3;
    double pys41_3;
    double pys51_3 ;
    double pys61_3 ;
    
    double pys12_3;
    double pys22_3;
    double pys32_3 ;
    double pys42_3;
    double pys52_3 ;
    double pys62_3 ;
    
    double pys11_4 ;
    double pys21_4 ;
    double pys31_4 ;
    double pys41_4 ;
    double pys51_4;
    double pys61_4 ;
    
    double pys12_4 ;
    double pys22_4 ;
    double pys32_4 ;
    double pys42_4;
    double pys52_4 ;
    double pys62_4 ;
    
    
    
    //////////////////////
    
    double pya11_1;
    double pya21_1;
    double pya31_1;
    double pya41_1;
    double pya51_1;
    double pya61_1;
    
    double pya12_1;
    double pya22_1;
    double pya32_1;
    double pya42_1;
    double pya52_1;
    double pya62_1;
    
//
    
    double pya11_2;
    double pya21_2;
    double pya31_2;
    double pya41_2;
    double pya51_2;
    double pya61_2;
    
    double pya12_2;
    double pya22_2;
    double pya32_2;
    double pya42_2;
    double pya52_2;
    double pya62_2;
    
///
    
    double pya11_3;
    double pya21_3;
    double pya31_3;
    double pya41_3;
    double pya51_3 ;
    double pya61_3 ;
    
    double pya12_3;
    double pya22_3;
    double pya32_3 ;
    double pya42_3;
    double pya52_3 ;
    double pya62_3 ;
    
    double pya11_4 ;
    double pya21_4 ;
    double pya31_4 ;
    double pya41_4 ;
    double pya51_4;
    double pya61_4 ;
    
    double pya12_4 ;
    double pya22_4 ;
    double pya32_4 ;
    double pya42_4;
    double pya52_4 ;
    double pya62_4 ;
    
    //////////////////////
    
    double pyb11_1;
    double pyb21_1;
    double pyb31_1;
    double pyb41_1;
    double pyb51_1;
    double pyb61_1;
    
    double pyb12_1;
    double pyb22_1;
    double pyb32_1;
    double pyb42_1;
    double pyb52_1;
    double pyb62_1;
    
//
    
    double pyb11_2;
    double pyb21_2;
    double pyb31_2;
    double pyb41_2;
    double pyb51_2;
    double pyb61_2;
    
    double pyb12_2;
    double pyb22_2;
    double pyb32_2;
    double pyb42_2;
    double pyb52_2;
    double pyb62_2;
    
///
    
    double pyb11_3;
    double pyb21_3;
    double pyb31_3;
    double pyb41_3;
    double pyb51_3 ;
    double pyb61_3 ;
    
    double pyb12_3;
    double pyb22_3;
    double pyb32_3 ;
    double pyb42_3;
    double pyb52_3 ;
    double pyb62_3 ;
    
    double pyb11_4 ;
    double pyb21_4 ;
    double pyb31_4 ;
    double pyb41_4 ;
    double pyb51_4;
    double pyb61_4 ;
    
    double pyb12_4 ;
    double pyb22_4 ;
    double pyb32_4 ;
    double pyb42_4;
    double pyb52_4 ;
    double pyb62_4 ;
    
    
    
    double pyw1_1 ;
    double pyw2_1 ;
    
    double pyw1_2;
    double pyw2_2 ;
    
    double pyw1_3;
    double pyw2_3;
    
    double pyw1_4 ;
    double pyw2_4 ;
    
    double eta1 = 1;
    double eta2 = 1;
    double eta3 = 1;
    double eta4 = 1;
    
    double alpha0 ;
    double alpha1 ;
    
    double rho;
    
    double xtilde_1 ;
    double xtilde_2 ;
    double xtilde_3 ;
    double xtilde_4;
    
    double thres = 20;
    
    
    
    t[0] = 0; x1[0] = 0; x2[0] = 0; x3[0] = 0; x4[0] = 0; u1[0] = 0; u2[0] = 0; x1f[0] = 0; x2f[0] = 0; x3f[0] = 0; x4f[0] = 0;
    t[1] = T; x1[1] = 0; x2[1] = 0; x3[1] = 0; x4[1] = 0; u1[1] = 0; u2[1] = 0; x1f[1] = 0; x2f[1] = 0; x3f[1] = 0; x4f[1] = 0;
    
    for (int k = 3 ; k<stp ; k++) {
        
        // SAMPLING
        
        t[k] = k*T;
        
        if (k==stp){
            t[k+1]= k*T;
        } else {
            t[k+1] = (k+1)*T;
        }
        
        xd3_0 = ax3 * cos(wx3 * t[k]);
        xd3_1 = ax3 * cos(wx3 * (t[k] + T));
        
        xd2_1 = ax2 * cos(wx2 * (t[k] + T));
        xd2_2 = ax2 * cos(wx2 * (t[k] + 2 * T));
        
        x1_k = x1[k-1]; x2_k = x2[k-1]; x3_k = x3[k-1]; x4_k = x4[k-1]; u1_k = u1[k-1]; u2_k = u2[k-1];
        x1_k1 = x1[k-2]; x2_k1 = x2[k-2]; x3_k1 = x3[k-2]; x4_k1 = x4[k-2]; u1_k1 = u1[k-2]; u2_k1 = u2[k-2];
        
        a1_k = x3_k + (3*T/2)*((-r1/l1)*(x3_k)+(-k1/l1)*x2_k) - (T/2)*((-r1/l1)*(x3_k1)+(-k1/l1)*x2_k1 + (1/l1)*u1_k1);
        a2_k = x4_k + (3*T/2)*((-r2/l2)*(x4_k)+(-k2/l2)*x2_k) - (T/2)*((-r2/l2)*(x4_k1)+(-k2/l2)*x2_k1 + (1/l2)*u2_k1);
        
        x1[k] = x1_k + (3*T/2)*x2_k - (T/2)*x2_k1;
        x2[k] = x2_k + (3*T/2)*((k1/j1)*x3_k + (k2/j2)*x4_k) - (T/2)*((k1/j1)*x3_k1 + (k2/j2)*x4_k1);
        x3[k] = a1_k + 1.5*(T/l1)*u1_k;
        x4[k] = a2_k + 1.5*(T/l2)*u2_k;
        
        
        if (t[k]>0.2) {
//  // FWNN ESTIMATION
            
            o2mat11_1 = exp(-(x1_k-m_mat11_1)*(x1_k-m_mat11_1))/(2*sigma_mat11_1*sigma_mat11_1);
            o2mat21_1 = exp(-(x2_k-m_mat21_1)*(x2_k-m_mat21_1))/(2*sigma_mat21_1*sigma_mat21_1);
            o2mat31_1 = exp(-(x3_k-m_mat31_1)*(x3_k-m_mat31_1))/(2*sigma_mat31_1*sigma_mat31_1);
            o2mat41_1 = exp(-(x4_k-m_mat41_1)*(x4_k-m_mat41_1))/(2*sigma_mat41_1*sigma_mat41_1);
            o2mat51_1 = exp(-(u1_k-m_mat51_1)*(u1_k-m_mat51_1))/(2*sigma_mat51_1*sigma_mat51_1);
            o2mat61_1 = exp(-(u2_k-m_mat61_1)*(u2_k-m_mat61_1))/(2*sigma_mat61_1*sigma_mat61_1);
            
            o2mat12_1 = exp(-(x1_k-m_mat12_1)*(x1_k-m_mat12_1))/(2*sigma_mat12_1*sigma_mat12_1);
            o2mat22_1 = exp(-(x2_k-m_mat22_1)*(x2_k-m_mat22_1))/(2*sigma_mat22_1*sigma_mat22_1);
            o2mat32_1 = exp(-(x3_k-m_mat32_1)*(x3_k-m_mat32_1))/(2*sigma_mat32_1*sigma_mat32_1);
            o2mat42_1 = exp(-(x4_k-m_mat42_1)*(x4_k-m_mat42_1))/(2*sigma_mat42_1*sigma_mat42_1);
            o2mat52_1 = exp(-(u1_k-m_mat52_1)*(u1_k-m_mat52_1))/(2*sigma_mat52_1*sigma_mat52_1);
            o2mat62_1 = exp(-(u2_k-m_mat62_1)*(u2_k-m_mat62_1))/(2*sigma_mat62_1*sigma_mat62_1);
            
            //
            o2mat11_2 = exp(-(x1_k-m_mat11_2)*(x1_k-m_mat11_2))/(2*sigma_mat11_2*sigma_mat11_3);
            o2mat21_2 = exp(-(x2_k-m_mat21_2)*(x2_k-m_mat21_2))/(2*sigma_mat21_2*sigma_mat21_3);
            o2mat31_2 = exp(-(x3_k-m_mat31_2)*(x3_k-m_mat31_2))/(2*sigma_mat31_2*sigma_mat31_3);
            o2mat41_2 = exp(-(x4_k-m_mat41_2)*(x4_k-m_mat41_2))/(2*sigma_mat41_2*sigma_mat41_3);
            o2mat51_2 = exp(-(u1_k-m_mat51_2)*(u1_k-m_mat51_2))/(2*sigma_mat51_2*sigma_mat51_3);
            o2mat61_2 = exp(-(u2_k-m_mat61_2)*(u2_k-m_mat61_2))/(2*sigma_mat61_2*sigma_mat61_3);
            
            o2mat12_2 = exp(-(x1_k-m_mat12_2)*(x1_k-m_mat12_2))/(2*sigma_mat12_2*sigma_mat12_2);
            o2mat22_2 = exp(-(x2_k-m_mat22_2)*(x2_k-m_mat22_2))/(2*sigma_mat22_2*sigma_mat22_2);
            o2mat32_2 = exp(-(x3_k-m_mat32_2)*(x3_k-m_mat32_2))/(2*sigma_mat32_2*sigma_mat32_2);
            o2mat42_2 = exp(-(x4_k-m_mat42_2)*(x4_k-m_mat42_2))/(2*sigma_mat42_2*sigma_mat42_2);
            o2mat52_2 = exp(-(u1_k-m_mat52_2)*(u1_k-m_mat52_2))/(2*sigma_mat52_2*sigma_mat52_2);
            o2mat62_2 = exp(-(u2_k-m_mat62_2)*(u2_k-m_mat62_2))/(2*sigma_mat62_2*sigma_mat62_2);
            ///
            o2mat11_3 = exp(-(x1_k-m_mat11_3)*(x1_k-m_mat11_3))/(2*sigma_mat11_3*sigma_mat11_3);
            o2mat21_3 = exp(-(x2_k-m_mat21_3)*(x2_k-m_mat21_3))/(2*sigma_mat21_3*sigma_mat21_3);
            o2mat31_3 = exp(-(x3_k-m_mat31_3)*(x3_k-m_mat31_3))/(2*sigma_mat31_3*sigma_mat31_3);
            o2mat41_3 = exp(-(x4_k-m_mat41_3)*(x4_k-m_mat41_3))/(2*sigma_mat41_3*sigma_mat41_3);
            o2mat51_3 = exp(-(u1_k-m_mat51_3)*(u1_k-m_mat51_3))/(2*sigma_mat51_3*sigma_mat51_3);
            o2mat61_3 = exp(-(u2_k-m_mat61_3)*(u2_k-m_mat61_3))/(2*sigma_mat61_3*sigma_mat61_3);
            
            o2mat12_3 = exp(-(x1_k-m_mat12_3)*(x1_k-m_mat12_3))/(2*sigma_mat12_3*sigma_mat12_3);
            o2mat22_3 = exp(-(x2_k-m_mat22_3)*(x2_k-m_mat22_3))/(2*sigma_mat22_3*sigma_mat22_3);
            o2mat32_3 = exp(-(x3_k-m_mat32_3)*(x3_k-m_mat32_3))/(2*sigma_mat32_3*sigma_mat32_3);
            o2mat42_3 = exp(-(x4_k-m_mat42_3)*(x4_k-m_mat42_3))/(2*sigma_mat42_3*sigma_mat42_3);
            o2mat52_3 = exp(-(u1_k-m_mat52_3)*(u1_k-m_mat52_3))/(2*sigma_mat52_3*sigma_mat52_3);
            o2mat62_3 = exp(-(u2_k-m_mat62_3)*(u2_k-m_mat62_3))/(2*sigma_mat62_3*sigma_mat62_3);
            
////
            o2mat11_4 = exp(-(x1_k-m_mat11_4)*(x1_k-m_mat11_4))/(2*sigma_mat11_4*sigma_mat11_4);
            o2mat21_4 = exp(-(x2_k-m_mat21_4)*(x2_k-m_mat21_4))/(2*sigma_mat21_4*sigma_mat21_4);
            o2mat31_4 = exp(-(x3_k-m_mat31_4)*(x3_k-m_mat31_4))/(2*sigma_mat31_4*sigma_mat31_4);
            o2mat41_4 = exp(-(x4_k-m_mat41_4)*(x4_k-m_mat41_4))/(2*sigma_mat41_4*sigma_mat41_4);
            o2mat51_4 = exp(-(u1_k-m_mat51_4)*(u1_k-m_mat51_4))/(2*sigma_mat51_4*sigma_mat51_4);
            o2mat61_4 = exp(-(u2_k-m_mat61_4)*(u2_k-m_mat61_4))/(2*sigma_mat61_4*sigma_mat61_4);
            
            o2mat12_4 = exp(-(x1_k-m_mat12_4)*(x1_k-m_mat12_4))/(2*sigma_mat12_4*sigma_mat12_4);
            o2mat22_4 = exp(-(x2_k-m_mat22_4)*(x2_k-m_mat22_4))/(2*sigma_mat22_4*sigma_mat22_4);
            o2mat32_4 = exp(-(x3_k-m_mat32_4)*(x3_k-m_mat32_4))/(2*sigma_mat32_4*sigma_mat32_4);
            o2mat42_4 = exp(-(x4_k-m_mat42_4)*(x4_k-m_mat42_4))/(2*sigma_mat42_4*sigma_mat42_4);
            o2mat52_4 = exp(-(u1_k-m_mat52_4)*(u1_k-m_mat52_4))/(2*sigma_mat52_4*sigma_mat52_4);
            o2mat62_4 = exp(-(u2_k-m_mat62_4)*(u2_k-m_mat62_4))/(2*sigma_mat62_4*sigma_mat62_4);
            
            
            o3mat1_1 = o2mat11_1*o2mat21_1*o2mat31_1*o2mat41_1*o2mat51_1*o2mat61_1;
            o3mat2_1 = o2mat12_1*o2mat22_1*o2mat32_1*o2mat42_1*o2mat52_1*o2mat62_1;
            
            o3mat1_2 = o2mat11_2*o2mat21_2*o2mat31_2*o2mat41_2*o2mat51_2*o2mat61_2;
            o3mat2_2 = o2mat12_2*o2mat22_2*o2mat32_2*o2mat42_2*o2mat52_2*o2mat62_2;
            
            o3mat1_3 = o2mat11_3*o2mat21_3*o2mat31_3*o2mat41_3*o2mat51_3*o2mat61_3;
            o3mat2_3 = o2mat12_3*o2mat22_3*o2mat32_3*o2mat42_3*o2mat52_3*o2mat62_3;
            
            o3mat1_4 = o2mat11_4*o2mat21_4*o2mat31_4*o2mat41_4*o2mat51_4*o2mat61_4;
            o3mat2_4 = o2mat12_4*o2mat22_4*o2mat32_4*o2mat42_4*o2mat52_4*o2mat62_4;
            
            o3sum_1 = o3mat1_1 + o3mat2_1;
            o3sum_2 = o3mat1_2 + o3mat2_2;
            o3sum_3 = o3mat1_3 + o3mat2_3;
            o3sum_4 = o3mat1_4 + o3mat2_4;
            
            o4mat1_1 = o3mat1_1/o3sum_1;
            o4mat2_1 = o3mat2_1/o3sum_1;
            
            o4mat1_2 = o3mat1_2/o3sum_2;
            o4mat2_2 = o3mat2_2/o3sum_2;
            
            o4mat1_3 = o3mat1_3/o3sum_3;
            o4mat2_3 = o3mat2_3/o3sum_3;
            
            o4mat1_4 = o3mat1_4/o3sum_4;
            o4mat2_4 = o3mat2_4/o3sum_4;
            
            if (o3sum_1==0){
                o4mat1_1 = 0.5;
                o4mat2_1 = 0.5;
            }
            
            if (o3sum_2==0){
                o4mat1_2 = 0.5;
                o4mat2_2 = 0.5;
            }
            
            if (o3sum_3==0){
                o4mat1_3 = 0.5;
                o4mat2_3 = 0.5;
            }
            
            if (o3sum_4==0){
                o4mat1_4 = 0.5;
                o4mat2_4 = 0.5;
            }
            
            
            // WAVELET
            
            nu_mat11_1 = (x1_k -b_mat11_1)/a_mat11_1;
            nu_mat21_1 = (x2_k -b_mat21_1)/a_mat21_1;
            nu_mat31_1 = (x3_k -b_mat31_1)/a_mat31_1;
            nu_mat41_1 = (x4_k -b_mat41_1)/a_mat41_1;
            nu_mat51_1 = (u1_k -b_mat51_1)/a_mat51_1;
            nu_mat61_1 = (u2_k -b_mat61_1)/a_mat61_1;
            
            nu_mat12_1 = (x1_k -b_mat12_1)/a_mat12_1;
            nu_mat22_1 = (x2_k -b_mat22_1)/a_mat22_1;
            nu_mat32_1 = (x3_k -b_mat32_1)/a_mat32_1;
            nu_mat42_1 = (x4_k -b_mat42_1)/a_mat42_1;
            nu_mat52_1 = (u1_k -b_mat52_1)/a_mat52_1;
            nu_mat62_1 = (u2_k -b_mat62_1)/a_mat62_1;
            
            //
            nu_mat11_2 = (x1_k -b_mat11_2)/a_mat11_2;
            nu_mat21_2 = (x2_k -b_mat21_2)/a_mat21_2;
            nu_mat31_2 = (x3_k -b_mat31_2)/a_mat31_2;
            nu_mat41_2 = (x4_k -b_mat41_2)/a_mat41_2;
            nu_mat51_2 = (u1_k -b_mat51_2)/a_mat51_2;
            nu_mat61_2 = (u2_k -b_mat61_2)/a_mat61_2;
            
            nu_mat12_2 = (x1_k -b_mat12_2)/a_mat12_2;
            nu_mat22_2 = (x2_k -b_mat22_2)/a_mat22_2;
            nu_mat32_2 = (x3_k -b_mat32_2)/a_mat32_2;
            nu_mat42_2 = (x4_k -b_mat42_2)/a_mat42_2;
            nu_mat52_2 = (u1_k -b_mat52_2)/a_mat52_2;
            nu_mat62_2 = (u2_k -b_mat62_2)/a_mat62_2;
            ///
            nu_mat11_3 = (x1_k -b_mat11_3)/a_mat11_3;
            nu_mat21_3 = (x2_k -b_mat21_3)/a_mat21_3;
            nu_mat31_3 = (x3_k -b_mat31_3)/a_mat31_3;
            nu_mat41_3 = (x4_k -b_mat41_3)/a_mat41_3;
            nu_mat51_3 = (u1_k -b_mat51_3)/a_mat51_3;
            nu_mat61_3 = (u2_k -b_mat61_3)/a_mat61_3;
            
            nu_mat12_3 = (x1_k -b_mat12_3)/a_mat12_3;
            nu_mat22_3 = (x2_k -b_mat22_3)/a_mat22_3;
            nu_mat32_3 = (x3_k -b_mat32_3)/a_mat32_3;
            nu_mat42_3 = (x4_k -b_mat42_3)/a_mat42_3;
            nu_mat52_3 = (u1_k -b_mat52_3)/a_mat52_3;
            nu_mat62_3 = (u2_k -b_mat62_3)/a_mat62_3;
            ////
            nu_mat11_4 = (x1_k -b_mat11_4)/a_mat11_4;
            nu_mat21_4 = (x2_k -b_mat21_4)/a_mat21_4;
            nu_mat31_4 = (x3_k -b_mat31_4)/a_mat31_4;
            nu_mat41_4 = (x4_k -b_mat41_4)/a_mat41_4;
            nu_mat51_4 = (u1_k -b_mat51_4)/a_mat51_4;
            nu_mat61_4 = (u2_k -b_mat61_4)/a_mat61_4;
            
            nu_mat12_4 = (x1_k -b_mat12_4)/a_mat12_4;
            nu_mat22_4 = (x2_k -b_mat22_4)/a_mat22_4;
            nu_mat32_4 = (x3_k -b_mat32_4)/a_mat32_4;
            nu_mat42_4 = (x4_k -b_mat42_4)/a_mat42_4;
            nu_mat52_4 = (u1_k -b_mat52_4)/a_mat52_4;
            nu_mat62_4 = (u2_k -b_mat62_4)/a_mat62_4;
            
            
            phi_mat11_1 = -nu_mat11_1*exp(-nu_mat11_1*nu_mat11_1/2.0);
            phi_mat21_1 = -nu_mat21_1*exp(-nu_mat21_1*nu_mat21_1/2.0);
            phi_mat31_1 = -nu_mat31_1*exp(-nu_mat31_1*nu_mat31_1/2.0);
            phi_mat41_1 = -nu_mat41_1*exp(-nu_mat41_1*nu_mat41_1/2.0);
            phi_mat51_1 = -nu_mat51_1*exp(-nu_mat51_1*nu_mat51_1/2.0);
            phi_mat61_1 = -nu_mat61_1*exp(-nu_mat61_1*nu_mat61_1/2.0);
            
            phi_mat12_1 = -nu_mat12_1*exp(-nu_mat12_1*nu_mat12_1/2.0);
            phi_mat22_1 = -nu_mat22_1*exp(-nu_mat22_1*nu_mat22_1/2.0);
            phi_mat32_1 = -nu_mat32_1*exp(-nu_mat32_1*nu_mat32_1/2.0);
            phi_mat42_1 = -nu_mat42_1*exp(-nu_mat42_1*nu_mat42_1/2.0);
            phi_mat52_1 = -nu_mat52_1*exp(-nu_mat52_1*nu_mat52_1/2.0);
            phi_mat62_1 = -nu_mat62_1*exp(-nu_mat62_1*nu_mat62_1/2.0);
            
            
            //
            phi_mat11_2 = -nu_mat11_2*exp(-nu_mat11_2*nu_mat11_2/2.0);
            phi_mat21_2 = -nu_mat21_2*exp(-nu_mat21_2*nu_mat21_2/2.0);
            phi_mat31_2 = -nu_mat31_2*exp(-nu_mat31_2*nu_mat31_2/2.0);
            phi_mat41_2 = -nu_mat41_2*exp(-nu_mat41_2*nu_mat41_2/2.0);
            phi_mat51_2 = -nu_mat51_2*exp(-nu_mat51_2*nu_mat51_2/2.0);
            phi_mat61_2 = -nu_mat61_2*exp(-nu_mat61_2*nu_mat61_2/2.0);
            
            phi_mat12_2 = -nu_mat12_2*exp(-nu_mat12_2*nu_mat12_2/2.0);
            phi_mat22_2 = -nu_mat22_2*exp(-nu_mat22_2*nu_mat22_2/2.0);
            phi_mat32_2 = -nu_mat32_2*exp(-nu_mat32_2*nu_mat32_2/2.0);
            phi_mat42_2 = -nu_mat42_2*exp(-nu_mat42_2*nu_mat42_2/2.0);
            phi_mat52_2 = -nu_mat52_2*exp(-nu_mat52_2*nu_mat52_2/2.0);
            phi_mat62_2 = -nu_mat62_2*exp(-nu_mat62_2*nu_mat62_2/2.0);
            
            ///
            phi_mat11_3 = -nu_mat11_3*exp(-nu_mat11_3*nu_mat11_3/2.0);
            phi_mat21_3 = -nu_mat21_3*exp(-nu_mat21_3*nu_mat21_3/2.0);
            phi_mat31_3 = -nu_mat31_3*exp(-nu_mat31_3*nu_mat31_3/2.0);
            phi_mat41_3 = -nu_mat41_3*exp(-nu_mat41_3*nu_mat41_3/2.0);
            phi_mat51_3 = -nu_mat51_3*exp(-nu_mat51_3*nu_mat51_3/2.0);
            phi_mat61_3 = -nu_mat61_3*exp(-nu_mat61_3*nu_mat61_3/2.0);
            
            phi_mat12_3 = -nu_mat12_3*exp(-nu_mat12_3*nu_mat12_3/2.0);
            phi_mat22_3 = -nu_mat22_3*exp(-nu_mat22_3*nu_mat22_3/2.0);
            phi_mat32_3 = -nu_mat32_3*exp(-nu_mat32_3*nu_mat32_3/2.0);
            phi_mat42_3 = -nu_mat42_3*exp(-nu_mat42_3*nu_mat42_3/2.0);
            phi_mat52_3 = -nu_mat52_3*exp(-nu_mat52_3*nu_mat52_3/2.0);
            phi_mat62_3 = -nu_mat62_3*exp(-nu_mat62_3*nu_mat62_3/2.0);
            
            ////
            phi_mat11_4 = -nu_mat11_4*exp(-nu_mat11_4*nu_mat11_4/2.0);
            phi_mat21_4 = -nu_mat21_4*exp(-nu_mat21_4*nu_mat21_4/2.0);
            phi_mat31_4 = -nu_mat31_4*exp(-nu_mat31_4*nu_mat31_4/2.0);
            phi_mat41_4 = -nu_mat41_4*exp(-nu_mat41_4*nu_mat41_4/2.0);
            phi_mat51_4 = -nu_mat51_4*exp(-nu_mat51_4*nu_mat51_4/2.0);
            phi_mat61_4 = -nu_mat61_4*exp(-nu_mat61_4*nu_mat61_4/2.0);
            
            phi_mat12_4 = -nu_mat12_4*exp(-nu_mat12_4*nu_mat12_4/2.0);
            phi_mat22_4 = -nu_mat22_4*exp(-nu_mat22_4*nu_mat22_4/2.0);
            phi_mat32_4 = -nu_mat32_4*exp(-nu_mat32_4*nu_mat32_4/2.0);
            phi_mat42_4 = -nu_mat42_4*exp(-nu_mat42_4*nu_mat42_4/2.0);
            phi_mat52_4 = -nu_mat52_4*exp(-nu_mat52_4*nu_mat52_4/2.0);
            phi_mat62_4 = -nu_mat62_4*exp(-nu_mat62_4*nu_mat62_4/2.0);
            
            
            
            
//
            omega_vec1_1 = (phi_mat11_1+phi_mat21_1+phi_mat31_1+phi_mat41_1+phi_mat51_1+phi_mat61_1)*weight_mat1_1;
            omega_vec2_1 = (phi_mat12_1+phi_mat22_1+phi_mat32_1+phi_mat42_1+phi_mat52_1+phi_mat62_1)*weight_mat2_1;
            
            omega_vec1_2 = (phi_mat11_2+phi_mat21_2+phi_mat31_2+phi_mat41_2+phi_mat51_2+phi_mat61_2)*weight_mat1_2;
            omega_vec2_2 = (phi_mat12_2+phi_mat22_2+phi_mat32_2+phi_mat42_2+phi_mat52_2+phi_mat62_2)*weight_mat2_2;
            
            omega_vec1_3 = (phi_mat11_3+phi_mat21_3+phi_mat31_3+phi_mat41_3+phi_mat51_3+phi_mat61_3)*weight_mat1_3;
            omega_vec2_3 = (phi_mat12_3+phi_mat22_3+phi_mat32_3+phi_mat42_3+phi_mat52_3+phi_mat62_3)*weight_mat2_3;
            
            omega_vec1_4 = (phi_mat11_4+phi_mat21_4+phi_mat31_4+phi_mat41_4+phi_mat51_4+phi_mat61_4)*weight_mat1_4;
            omega_vec2_4 = (phi_mat12_4+phi_mat22_4+phi_mat32_4+phi_mat42_4+phi_mat52_4+phi_mat62_4)*weight_mat2_4;
            
            xf_1 = omega_vec1_1*o4mat1_1 + omega_vec2_1*o4mat2_1;
            xf_2 = omega_vec1_2*o4mat1_2 + omega_vec2_2*o4mat2_2;
            xf_3 = omega_vec1_3*o4mat1_3 + omega_vec2_3*o4mat2_3;
            xf_4 = omega_vec1_4*o4mat1_4 + omega_vec2_4*o4mat2_4;
            
            x1f[k] = xf_1;
            x2f[k] = xf_2;
            x3f[k] = xf_3;
            x4f[k] = xf_4;
            
            
            
            o5mat1_1 = omega_vec1_1*o4mat1_1;
            o5mat2_1 = omega_vec2_1*o4mat2_1;
            
            o5mat1_2 = omega_vec1_2*o4mat1_2;
            o5mat2_2 = omega_vec2_2*o4mat2_2;
            
            o5mat1_3 = omega_vec1_3*o4mat1_3;
            o5mat2_3 = omega_vec2_3*o4mat2_3;
            
            o5mat1_4 = omega_vec1_4*o4mat1_4;
            o5mat2_4 = omega_vec2_4*o4mat2_4;
            
//  // FWNN UPDATE
            
            pym11_1 = (o5mat1_1 - xf_1)*o4mat1_1*(x1_k-m_mat11_1)/(sigma_mat11_1*sigma_mat11_1);
            pym21_1 = (o5mat1_1 - xf_1)*o4mat1_1*(x2_k-m_mat21_1)/(sigma_mat21_1*sigma_mat21_1);
            pym31_1 = (o5mat1_1 - xf_1)*o4mat1_1*(x3_k-m_mat31_1)/(sigma_mat31_1*sigma_mat31_1);
            pym41_1 = (o5mat1_1 - xf_1)*o4mat1_1*(x4_k-m_mat41_1)/(sigma_mat41_1*sigma_mat41_1);
            pym51_1 = (o5mat1_1 - xf_1)*o4mat1_1*(u1_k-m_mat51_1)/(sigma_mat51_1*sigma_mat51_1);
            pym61_1 = (o5mat1_1 - xf_1)*o4mat1_1*(u2_k-m_mat61_1)/(sigma_mat61_1*sigma_mat61_1);
            
            pym12_1 = (o5mat2_1 - xf_1)*o4mat2_1*(x1_k-m_mat12_1)/(sigma_mat12_1*sigma_mat12_1);
            pym22_1 = (o5mat2_1 - xf_1)*o4mat2_1*(x2_k-m_mat22_1)/(sigma_mat22_1*sigma_mat22_1);
            pym32_1 = (o5mat2_1 - xf_1)*o4mat2_1*(x3_k-m_mat32_1)/(sigma_mat32_1*sigma_mat32_1);
            pym42_1 = (o5mat2_1 - xf_1)*o4mat2_1*(x4_k-m_mat42_1)/(sigma_mat42_1*sigma_mat42_1);
            pym52_1 = (o5mat2_1 - xf_1)*o4mat2_1*(u1_k-m_mat52_1)/(sigma_mat52_1*sigma_mat52_1);
            pym62_1 = (o5mat2_1 - xf_1)*o4mat2_1*(u2_k-m_mat62_1)/(sigma_mat62_1*sigma_mat62_1);
            
            //
            pym11_2 = (o5mat1_2 - xf_2)*o4mat1_2*(x1_k-m_mat11_2)/(sigma_mat11_2*sigma_mat11_2);
            pym21_2 = (o5mat1_2 - xf_2)*o4mat1_2*(x2_k-m_mat21_2)/(sigma_mat21_2*sigma_mat21_2);
            pym31_2 = (o5mat1_2 - xf_2)*o4mat1_2*(x3_k-m_mat31_2)/(sigma_mat31_2*sigma_mat31_2);
            pym41_2 = (o5mat1_2 - xf_2)*o4mat1_2*(x4_k-m_mat41_2)/(sigma_mat41_2*sigma_mat41_2);
            pym51_2 = (o5mat1_2 - xf_2)*o4mat1_2*(u1_k-m_mat51_2)/(sigma_mat51_2*sigma_mat51_2);
            pym61_2 = (o5mat1_2 - xf_2)*o4mat1_2*(u2_k-m_mat61_2)/(sigma_mat61_2*sigma_mat61_2);
            
            pym12_2 = (o5mat2_2 - xf_2)*o4mat2_2*(x1_k-m_mat12_2)/(sigma_mat12_2*sigma_mat12_2);
            pym22_2 = (o5mat2_2 - xf_2)*o4mat2_2*(x2_k-m_mat22_2)/(sigma_mat22_2*sigma_mat22_2);
            pym32_2 = (o5mat2_2 - xf_2)*o4mat2_2*(x3_k-m_mat32_2)/(sigma_mat32_2*sigma_mat32_2);
            pym42_2 = (o5mat2_2 - xf_2)*o4mat2_2*(x4_k-m_mat42_2)/(sigma_mat42_2*sigma_mat42_2);
            pym52_2 = (o5mat2_2 - xf_2)*o4mat2_2*(u1_k-m_mat52_2)/(sigma_mat52_2*sigma_mat52_2);
            pym62_2 = (o5mat2_2 - xf_2)*o4mat2_2*(u2_k-m_mat62_2)/(sigma_mat62_2*sigma_mat62_2);
            ///
            pym11_3 = (o5mat1_3 - xf_3)*o4mat1_3*(x1_k-m_mat11_3)/(sigma_mat11_3*sigma_mat11_3);
            pym21_3 = (o5mat1_3 - xf_3)*o4mat1_3*(x2_k-m_mat21_3)/(sigma_mat21_3*sigma_mat21_3);
            pym31_3 = (o5mat1_3 - xf_3)*o4mat1_3*(x3_k-m_mat31_3)/(sigma_mat31_3*sigma_mat31_3);
            pym41_3 = (o5mat1_3 - xf_3)*o4mat1_3*(x4_k-m_mat41_3)/(sigma_mat41_3*sigma_mat41_3);
            pym51_3 = (o5mat1_3 - xf_3)*o4mat1_3*(u1_k-m_mat51_3)/(sigma_mat51_3*sigma_mat51_3);
            pym61_3 = (o5mat1_3 - xf_3)*o4mat1_3*(u2_k-m_mat61_3)/(sigma_mat61_3*sigma_mat61_3);
            
            pym12_3 = (o5mat2_3 - xf_3)*o4mat2_3*(x1_k-m_mat12_3)/(sigma_mat12_3*sigma_mat12_3);
            pym22_3 = (o5mat2_3 - xf_3)*o4mat2_3*(x2_k-m_mat22_3)/(sigma_mat22_3*sigma_mat22_3);
            pym32_3 = (o5mat2_3 - xf_3)*o4mat2_3*(x3_k-m_mat32_3)/(sigma_mat32_3*sigma_mat32_3);
            pym42_3 = (o5mat2_3 - xf_3)*o4mat2_3*(x4_k-m_mat42_3)/(sigma_mat42_3*sigma_mat42_3);
            pym52_3 = (o5mat2_3 - xf_3)*o4mat2_3*(u1_k-m_mat52_3)/(sigma_mat52_3*sigma_mat52_3);
            pym62_3 = (o5mat2_3 - xf_3)*o4mat2_3*(u2_k-m_mat62_3)/(sigma_mat62_3*sigma_mat62_3);
            ////
            pym11_4 = (o5mat1_4 - xf_4)*o4mat1_4*(x1_k-m_mat11_4)/(sigma_mat11_4*sigma_mat11_4);
            pym21_4 = (o5mat1_4 - xf_4)*o4mat1_4*(x2_k-m_mat21_4)/(sigma_mat21_4*sigma_mat21_4);
            pym31_4 = (o5mat1_4 - xf_4)*o4mat1_4*(x3_k-m_mat31_4)/(sigma_mat31_4*sigma_mat31_4);
            pym41_4 = (o5mat1_4 - xf_4)*o4mat1_4*(x4_k-m_mat41_4)/(sigma_mat41_4*sigma_mat41_4);
            pym51_4 = (o5mat1_4 - xf_4)*o4mat1_4*(u1_k-m_mat51_4)/(sigma_mat51_4*sigma_mat51_4);
            pym61_4 = (o5mat1_4 - xf_4)*o4mat1_4*(u2_k-m_mat61_4)/(sigma_mat61_4*sigma_mat61_4);
            
            pym12_4 = (o5mat2_4 - xf_4)*o4mat2_4*(x1_k-m_mat12_4)/(sigma_mat12_4*sigma_mat12_4);
            pym22_4 = (o5mat2_4 - xf_4)*o4mat2_4*(x2_k-m_mat22_4)/(sigma_mat22_4*sigma_mat22_4);
            pym32_4 = (o5mat2_4 - xf_4)*o4mat2_4*(x3_k-m_mat32_4)/(sigma_mat32_4*sigma_mat32_4);
            pym42_4 = (o5mat2_4 - xf_4)*o4mat2_4*(x4_k-m_mat42_4)/(sigma_mat42_4*sigma_mat42_4);
            pym52_4 = (o5mat2_4 - xf_4)*o4mat2_4*(u1_k-m_mat52_4)/(sigma_mat52_4*sigma_mat52_4);
            pym62_4 = (o5mat2_4 - xf_4)*o4mat2_4*(u2_k-m_mat62_4)/(sigma_mat62_4*sigma_mat62_4);
            
            
            
            
            ////////////////////
            pys11_1 = (o5mat1_1 - xf_1)*o4mat1_1*(x1_k-m_mat11_1)*(x1_k-m_mat11_1)/(sigma_mat11_1*sigma_mat11_1*sigma_mat11_1);
            pys21_1 = (o5mat1_1 - xf_1)*o4mat1_1*(x2_k-m_mat21_1)*(x2_k-m_mat21_1)/(sigma_mat21_1*sigma_mat21_1*sigma_mat21_1);
            pys31_1 = (o5mat1_1 - xf_1)*o4mat1_1*(x3_k-m_mat31_1)*(x3_k-m_mat31_1)/(sigma_mat31_1*sigma_mat31_1*sigma_mat31_1);
            pys41_1 = (o5mat1_1 - xf_1)*o4mat1_1*(x4_k-m_mat41_1)*(x4_k-m_mat41_1)/(sigma_mat41_1*sigma_mat41_1*sigma_mat41_1);
            pys51_1 = (o5mat1_1 - xf_1)*o4mat1_1*(u1_k-m_mat51_1)*(u1_k-m_mat51_1)/(sigma_mat51_1*sigma_mat51_1*sigma_mat51_1);
            pys61_1 = (o5mat1_1 - xf_1)*o4mat1_1*(u2_k-m_mat61_1)*(u2_k-m_mat61_1)/(sigma_mat61_1*sigma_mat61_1*sigma_mat61_1);
            
            pys12_1 = (o5mat2_1 - xf_1)*o4mat2_1*(x1_k-m_mat12_1)*(x1_k-m_mat12_1)/(sigma_mat12_1*sigma_mat12_1*sigma_mat12_1);
            pys22_1 = (o5mat2_1 - xf_1)*o4mat2_1*(x2_k-m_mat22_1)*(x2_k-m_mat22_1)/(sigma_mat22_1*sigma_mat22_1*sigma_mat22_1);
            pys32_1 = (o5mat2_1 - xf_1)*o4mat2_1*(x3_k-m_mat32_1)*(x3_k-m_mat32_1)/(sigma_mat32_1*sigma_mat32_1*sigma_mat32_1);
            pys42_1 = (o5mat2_1 - xf_1)*o4mat2_1*(x4_k-m_mat42_1)*(x4_k-m_mat42_1)/(sigma_mat42_1*sigma_mat42_1*sigma_mat42_1);
            pys52_1 = (o5mat2_1 - xf_1)*o4mat2_1*(u1_k-m_mat52_1)*(u1_k-m_mat52_1)/(sigma_mat52_1*sigma_mat52_1*sigma_mat52_1);
            pys62_1 = (o5mat2_1 - xf_1)*o4mat2_1*(u2_k-m_mat62_1)*(u2_k-m_mat62_1)/(sigma_mat62_1*sigma_mat62_1*sigma_mat62_1);
            
            //
            pys11_2 = (o5mat1_2 - xf_2)*o4mat1_2*(x1_k-m_mat11_2)*(x1_k-m_mat11_2)/(sigma_mat11_2*sigma_mat11_2*sigma_mat11_2);
            pys21_2 = (o5mat1_2 - xf_2)*o4mat1_2*(x2_k-m_mat21_2)*(x2_k-m_mat21_2)/(sigma_mat21_2*sigma_mat21_2*sigma_mat21_2);
            pys31_2 = (o5mat1_2 - xf_2)*o4mat1_2*(x3_k-m_mat31_2)*(x3_k-m_mat31_2)/(sigma_mat31_2*sigma_mat31_2*sigma_mat31_2);
            pys41_2 = (o5mat1_2 - xf_2)*o4mat1_2*(x4_k-m_mat41_2)*(x4_k-m_mat41_2)/(sigma_mat41_2*sigma_mat41_2*sigma_mat41_2);
            pys51_2 = (o5mat1_2 - xf_2)*o4mat1_2*(u1_k-m_mat51_2)*(u1_k-m_mat51_2)/(sigma_mat51_2*sigma_mat51_2*sigma_mat51_2);
            pys61_2 = (o5mat1_2 - xf_2)*o4mat1_2*(u2_k-m_mat61_2)*(u2_k-m_mat61_2)/(sigma_mat61_2*sigma_mat61_2*sigma_mat61_2);
            
            pys12_2 = (o5mat2_2 - xf_2)*o4mat2_2*(x1_k-m_mat12_2)*(x1_k-m_mat12_2)/(sigma_mat12_2*sigma_mat12_2*sigma_mat12_2);
            pys22_2 = (o5mat2_2 - xf_2)*o4mat2_2*(x2_k-m_mat22_2)*(x2_k-m_mat22_2)/(sigma_mat22_2*sigma_mat22_2*sigma_mat22_2);
            pys32_2 = (o5mat2_2 - xf_2)*o4mat2_2*(x3_k-m_mat32_2)*(x3_k-m_mat32_2)/(sigma_mat32_2*sigma_mat32_2*sigma_mat32_2);
            pys42_2 = (o5mat2_2 - xf_2)*o4mat2_2*(x4_k-m_mat42_2)*(x4_k-m_mat42_2)/(sigma_mat42_2*sigma_mat42_2*sigma_mat42_2);
            pys52_2 = (o5mat2_2 - xf_2)*o4mat2_2*(u1_k-m_mat52_2)*(u1_k-m_mat52_2)/(sigma_mat52_2*sigma_mat52_2*sigma_mat52_2);
            pys62_2 = (o5mat2_2 - xf_2)*o4mat2_2*(u2_k-m_mat62_2)*(u2_k-m_mat62_2)/(sigma_mat62_2*sigma_mat62_2*sigma_mat62_2);
            ///
            pys11_3 = (o5mat1_3 - xf_3)*o4mat1_3*(x1_k-m_mat11_3)*(x1_k-m_mat11_3)/(sigma_mat11_3*sigma_mat11_3*sigma_mat11_3);
            pys21_3 = (o5mat1_3 - xf_3)*o4mat1_3*(x2_k-m_mat21_3)*(x2_k-m_mat21_3)/(sigma_mat21_3*sigma_mat21_3*sigma_mat21_3);
            pys31_3 = (o5mat1_3 - xf_3)*o4mat1_3*(x3_k-m_mat31_3)*(x3_k-m_mat31_3)/(sigma_mat31_3*sigma_mat31_3*sigma_mat31_3);
            pys41_3 = (o5mat1_3 - xf_3)*o4mat1_3*(x4_k-m_mat41_3)*(x4_k-m_mat41_3)/(sigma_mat41_3*sigma_mat41_3*sigma_mat41_3);
            pys51_3 = (o5mat1_3 - xf_3)*o4mat1_3*(u1_k-m_mat51_3)*(u1_k-m_mat51_3)/(sigma_mat51_3*sigma_mat51_3*sigma_mat51_3);
            pys61_3 = (o5mat1_3 - xf_3)*o4mat1_3*(u2_k-m_mat61_3)*(u2_k-m_mat61_3)/(sigma_mat61_3*sigma_mat61_3*sigma_mat61_3);
            
            pys12_3 = (o5mat2_3 - xf_3)*o4mat2_3*(x1_k-m_mat12_3)*(x1_k-m_mat12_3)/(sigma_mat12_3*sigma_mat12_3*sigma_mat12_3);
            pys22_3 = (o5mat2_3 - xf_3)*o4mat2_3*(x2_k-m_mat22_3)*(x2_k-m_mat22_3)/(sigma_mat22_3*sigma_mat22_3*sigma_mat22_3);
            pys32_3 = (o5mat2_3 - xf_3)*o4mat2_3*(x3_k-m_mat32_3)*(x3_k-m_mat32_3)/(sigma_mat32_3*sigma_mat32_3*sigma_mat32_3);
            pys42_3 = (o5mat2_3 - xf_3)*o4mat2_3*(x4_k-m_mat42_3)*(x4_k-m_mat42_3)/(sigma_mat42_3*sigma_mat42_3*sigma_mat42_3);
            pys52_3 = (o5mat2_3 - xf_3)*o4mat2_3*(u1_k-m_mat52_3)*(u1_k-m_mat52_3)/(sigma_mat52_3*sigma_mat52_3*sigma_mat52_3);
            pys62_3 = (o5mat2_3 - xf_3)*o4mat2_3*(u2_k-m_mat62_3)*(u2_k-m_mat62_3)/(sigma_mat62_3*sigma_mat62_3*sigma_mat62_3);
            ////
            
            
            pys11_4 = (o5mat1_4 - xf_4)*o4mat1_4*(x1_k-m_mat11_4)*(x1_k-m_mat11_4)/(sigma_mat11_4*sigma_mat11_4*sigma_mat11_4);
            pys21_4 = (o5mat1_4 - xf_4)*o4mat1_4*(x2_k-m_mat21_4)*(x2_k-m_mat21_4)/(sigma_mat21_4*sigma_mat21_4*sigma_mat21_4);
            pys31_4 = (o5mat1_4 - xf_4)*o4mat1_4*(x3_k-m_mat31_4)*(x3_k-m_mat31_4)/(sigma_mat31_4*sigma_mat31_4*sigma_mat31_4);
            pys41_4 = (o5mat1_4 - xf_4)*o4mat1_4*(x4_k-m_mat41_4)*(x4_k-m_mat41_4)/(sigma_mat41_4*sigma_mat41_4*sigma_mat41_4);
            pys51_4 = (o5mat1_4 - xf_4)*o4mat1_4*(u1_k-m_mat51_4)*(u1_k-m_mat51_4)/(sigma_mat51_4*sigma_mat51_4*sigma_mat51_4);
            pys61_4 = (o5mat1_4 - xf_4)*o4mat1_4*(u2_k-m_mat61_4)*(u2_k-m_mat61_4)/(sigma_mat61_4*sigma_mat61_4*sigma_mat61_4);
            
            pys12_4 = (o5mat2_4 - xf_4)*o4mat2_4*(x1_k-m_mat12_4)*(x1_k-m_mat12_4)/(sigma_mat12_4*sigma_mat12_4*sigma_mat12_4);
            pys22_4 = (o5mat2_4 - xf_4)*o4mat2_4*(x2_k-m_mat22_4)*(x2_k-m_mat22_4)/(sigma_mat22_4*sigma_mat22_4*sigma_mat22_4);
            pys32_4 = (o5mat2_4 - xf_4)*o4mat2_4*(x3_k-m_mat32_4)*(x3_k-m_mat32_4)/(sigma_mat32_4*sigma_mat32_4*sigma_mat32_4);
            pys42_4 = (o5mat2_4 - xf_4)*o4mat2_4*(x4_k-m_mat42_4)*(x4_k-m_mat42_4)/(sigma_mat42_4*sigma_mat42_4*sigma_mat42_4);
            pys52_4 = (o5mat2_4 - xf_4)*o4mat2_4*(u1_k-m_mat52_4)*(u1_k-m_mat52_4)/(sigma_mat52_4*sigma_mat52_4*sigma_mat52_4);
            pys62_4 = (o5mat2_4 - xf_4)*o4mat2_4*(u2_k-m_mat62_4)*(u2_k-m_mat62_4)/(sigma_mat62_4*sigma_mat62_4*sigma_mat62_4);
            
            
            
            //////////////////////
            pya11_1 = o4mat1_1*weight_mat1_1*phi_mat11_1*(nu_mat11_1*nu_mat11_1 - 1)/(a_mat11_1);
            pya21_1 = o4mat1_1*weight_mat1_1*phi_mat21_1*(nu_mat21_1*nu_mat21_1 - 1)/(a_mat21_1);
            pya31_1 = o4mat1_1*weight_mat1_1*phi_mat31_1*(nu_mat31_1*nu_mat31_1 - 1)/(a_mat31_1);
            pya41_1 = o4mat1_1*weight_mat1_1*phi_mat41_1*(nu_mat41_1*nu_mat41_1 - 1)/(a_mat41_1);
            pya51_1 = o4mat1_1*weight_mat1_1*phi_mat51_1*(nu_mat51_1*nu_mat51_1 - 1)/(a_mat51_1);
            pya61_1 = o4mat1_1*weight_mat1_1*phi_mat61_1*(nu_mat61_1*nu_mat61_1 - 1)/(a_mat61_1);
            
            pya12_1 = o4mat2_1*weight_mat2_1*phi_mat12_1*(nu_mat12_1*nu_mat12_1 - 1)/(a_mat12_1);
            pya22_1 = o4mat2_1*weight_mat2_1*phi_mat22_1*(nu_mat22_1*nu_mat22_1 - 1)/(a_mat22_1);
            pya32_1 = o4mat2_1*weight_mat2_1*phi_mat32_1*(nu_mat32_1*nu_mat32_1 - 1)/(a_mat32_1);
            pya42_1 = o4mat2_1*weight_mat2_1*phi_mat42_1*(nu_mat42_1*nu_mat42_1 - 1)/(a_mat42_1);
            pya52_1 = o4mat2_1*weight_mat2_1*phi_mat52_1*(nu_mat52_1*nu_mat52_1 - 1)/(a_mat52_1);
            pya62_1 = o4mat2_1*weight_mat2_1*phi_mat62_1*(nu_mat62_1*nu_mat62_1 - 1)/(a_mat62_1);
            //
            pya11_2 = o4mat1_2*weight_mat1_2*phi_mat11_2*(nu_mat11_2*nu_mat11_2 - 1)/(a_mat11_2);
            pya21_2 = o4mat1_2*weight_mat1_2*phi_mat21_2*(nu_mat21_2*nu_mat21_2 - 1)/(a_mat21_2);
            pya31_2 = o4mat1_2*weight_mat1_2*phi_mat31_2*(nu_mat31_2*nu_mat31_2 - 1)/(a_mat31_2);
            pya41_2 = o4mat1_2*weight_mat1_2*phi_mat41_2*(nu_mat41_2*nu_mat41_2 - 1)/(a_mat41_2);
            pya51_2 = o4mat1_2*weight_mat1_2*phi_mat51_2*(nu_mat51_2*nu_mat51_2 - 1)/(a_mat51_2);
            pya61_2 = o4mat1_2*weight_mat1_2*phi_mat61_2*(nu_mat61_2*nu_mat61_2 - 1)/(a_mat61_2);
            
            pya12_2 = o4mat2_2*weight_mat2_2*phi_mat12_2*(nu_mat12_2*nu_mat12_2 - 1)/(a_mat12_2);
            pya22_2 = o4mat2_2*weight_mat2_2*phi_mat22_2*(nu_mat22_2*nu_mat22_2 - 1)/(a_mat22_2);
            pya32_2 = o4mat2_2*weight_mat2_2*phi_mat32_2*(nu_mat32_2*nu_mat32_2 - 1)/(a_mat32_2);
            pya42_2 = o4mat2_2*weight_mat2_2*phi_mat42_2*(nu_mat42_2*nu_mat42_2 - 1)/(a_mat42_2);
            pya52_2 = o4mat2_2*weight_mat2_2*phi_mat52_2*(nu_mat52_2*nu_mat52_2 - 1)/(a_mat52_2);
            pya62_2 = o4mat2_2*weight_mat2_2*phi_mat62_2*(nu_mat62_2*nu_mat62_2 - 1)/(a_mat62_2);
            ///
            pya11_3 = o4mat1_3*weight_mat1_3*phi_mat11_3*(nu_mat11_3*nu_mat11_3 - 1)/(a_mat11_3);
            pya21_3 = o4mat1_3*weight_mat1_3*phi_mat21_3*(nu_mat21_3*nu_mat21_3 - 1)/(a_mat21_3);
            pya31_3 = o4mat1_3*weight_mat1_3*phi_mat31_3*(nu_mat31_3*nu_mat31_3 - 1)/(a_mat31_3);
            pya41_3 = o4mat1_3*weight_mat1_3*phi_mat41_3*(nu_mat41_3*nu_mat41_3 - 1)/(a_mat41_3);
            pya51_3 = o4mat1_3*weight_mat1_3*phi_mat51_3*(nu_mat51_3*nu_mat51_3 - 1)/(a_mat51_3);
            pya61_3 = o4mat1_3*weight_mat1_3*phi_mat61_3*(nu_mat61_3*nu_mat61_3 - 1)/(a_mat61_3);
            
            pya12_3 = o4mat2_3*weight_mat2_3*phi_mat12_3*(nu_mat12_3*nu_mat12_3 - 1)/(a_mat12_3);
            pya22_3 = o4mat2_3*weight_mat2_3*phi_mat22_3*(nu_mat22_3*nu_mat22_3 - 1)/(a_mat22_3);
            pya32_3 = o4mat2_3*weight_mat2_3*phi_mat32_3*(nu_mat32_3*nu_mat32_3 - 1)/(a_mat32_3);
            pya42_3 = o4mat2_3*weight_mat2_3*phi_mat42_3*(nu_mat42_3*nu_mat42_3 - 1)/(a_mat42_3);
            pya52_3 = o4mat2_3*weight_mat2_3*phi_mat52_3*(nu_mat52_3*nu_mat52_3 - 1)/(a_mat52_3);
            pya62_3 = o4mat2_3*weight_mat2_3*phi_mat62_3*(nu_mat62_3*nu_mat62_3 - 1)/(a_mat62_3);
            ////
            pya11_4 = o4mat1_4*weight_mat1_4*phi_mat11_4*(nu_mat11_4*nu_mat11_4 - 1)/(a_mat11_4);
            pya21_4 = o4mat1_4*weight_mat1_4*phi_mat21_4*(nu_mat21_4*nu_mat21_4 - 1)/(a_mat21_4);
            pya31_4 = o4mat1_4*weight_mat1_4*phi_mat31_4*(nu_mat31_4*nu_mat31_4 - 1)/(a_mat31_4);
            pya41_4 = o4mat1_4*weight_mat1_4*phi_mat41_4*(nu_mat41_4*nu_mat41_4 - 1)/(a_mat41_4);
            pya51_4 = o4mat1_4*weight_mat1_4*phi_mat51_4*(nu_mat51_4*nu_mat51_4 - 1)/(a_mat51_4);
            pya61_4 = o4mat1_4*weight_mat1_4*phi_mat61_4*(nu_mat61_4*nu_mat61_4 - 1)/(a_mat61_4);
            
            pya12_4 = o4mat2_4*weight_mat2_4*phi_mat12_4*(nu_mat12_4*nu_mat12_4 - 1)/(a_mat12_4);
            pya22_4 = o4mat2_4*weight_mat2_4*phi_mat22_4*(nu_mat22_4*nu_mat22_4 - 1)/(a_mat22_4);
            pya32_4 = o4mat2_4*weight_mat2_4*phi_mat32_4*(nu_mat32_4*nu_mat32_4 - 1)/(a_mat32_4);
            pya42_4 = o4mat2_4*weight_mat2_4*phi_mat42_4*(nu_mat42_4*nu_mat42_4 - 1)/(a_mat42_4);
            pya52_4 = o4mat2_4*weight_mat2_4*phi_mat52_4*(nu_mat52_4*nu_mat52_4 - 1)/(a_mat52_4);
            pya62_4 = o4mat2_4*weight_mat2_4*phi_mat62_4*(nu_mat62_4*nu_mat62_4 - 1)/(a_mat62_4);
            
            //////////////////////
            
            pyb11_1 = o4mat1_1*weight_mat1_1*phi_mat11_1*(nu_mat11_1*nu_mat11_1 - 1)/(nu_mat11_1*a_mat11_1);
            pyb21_1 = o4mat1_1*weight_mat1_1*phi_mat21_1*(nu_mat21_1*nu_mat21_1 - 1)/(nu_mat21_1*a_mat21_1);
            pyb31_1 = o4mat1_1*weight_mat1_1*phi_mat31_1*(nu_mat31_1*nu_mat31_1 - 1)/(nu_mat31_1*a_mat31_1);
            pyb41_1 = o4mat1_1*weight_mat1_1*phi_mat41_1*(nu_mat41_1*nu_mat41_1 - 1)/(nu_mat41_1*a_mat41_1);
            pyb51_1 = o4mat1_1*weight_mat1_1*phi_mat51_1*(nu_mat51_1*nu_mat51_1 - 1)/(nu_mat51_1*a_mat51_1);
            pyb61_1 = o4mat1_1*weight_mat1_1*phi_mat61_1*(nu_mat61_1*nu_mat61_1 - 1)/(nu_mat61_1*a_mat61_1);
            
            pyb12_1 = o4mat2_1*weight_mat2_1*phi_mat12_1*(nu_mat12_1*nu_mat12_1 - 1)/(nu_mat12_1*a_mat12_1);
            pyb22_1 = o4mat2_1*weight_mat2_1*phi_mat22_1*(nu_mat22_1*nu_mat22_1 - 1)/(nu_mat22_1*a_mat22_1);
            pyb32_1 = o4mat2_1*weight_mat2_1*phi_mat32_1*(nu_mat32_1*nu_mat32_1 - 1)/(nu_mat32_1*a_mat32_1);
            pyb42_1 = o4mat2_1*weight_mat2_1*phi_mat42_1*(nu_mat42_1*nu_mat42_1 - 1)/(nu_mat42_1*a_mat42_1);
            pyb52_1 = o4mat2_1*weight_mat2_1*phi_mat52_1*(nu_mat52_1*nu_mat52_1 - 1)/(nu_mat52_1*a_mat52_1);
            pyb62_1 = o4mat2_1*weight_mat2_1*phi_mat62_1*(nu_mat62_1*nu_mat62_1 - 1)/(nu_mat62_1*a_mat62_1);
            
//
            
            pyb11_2 = o4mat1_2*weight_mat1_2*phi_mat11_2*(nu_mat11_2*nu_mat11_2 - 1)/(nu_mat11_2*a_mat11_2);
            pyb21_2 = o4mat1_2*weight_mat1_2*phi_mat21_2*(nu_mat21_2*nu_mat21_2 - 1)/(nu_mat21_2*a_mat21_2);
            pyb31_2 = o4mat1_2*weight_mat1_2*phi_mat31_2*(nu_mat31_2*nu_mat31_2 - 1)/(nu_mat31_2*a_mat31_2);
            pyb41_2 = o4mat1_2*weight_mat1_2*phi_mat41_2*(nu_mat41_2*nu_mat41_2 - 1)/(nu_mat41_2*a_mat41_2);
            pyb51_2 = o4mat1_2*weight_mat1_2*phi_mat51_2*(nu_mat51_2*nu_mat51_2 - 1)/(nu_mat51_2*a_mat51_2);
            pyb61_2 = o4mat1_2*weight_mat1_2*phi_mat61_2*(nu_mat61_2*nu_mat61_2 - 1)/(nu_mat61_2*a_mat61_2);
            
            pyb12_2 = o4mat2_2*weight_mat2_2*phi_mat12_2*(nu_mat12_2*nu_mat12_2 - 1)/(nu_mat12_2*a_mat12_2);
            pyb22_2 = o4mat2_2*weight_mat2_2*phi_mat22_2*(nu_mat22_2*nu_mat22_2 - 1)/(nu_mat22_2*a_mat22_2);
            pyb32_2 = o4mat2_2*weight_mat2_2*phi_mat32_2*(nu_mat32_2*nu_mat32_2 - 1)/(nu_mat32_2*a_mat32_2);
            pyb42_2 = o4mat2_2*weight_mat2_2*phi_mat42_2*(nu_mat42_2*nu_mat42_2 - 1)/(nu_mat42_2*a_mat42_2);
            pyb52_2 = o4mat2_2*weight_mat2_2*phi_mat52_2*(nu_mat52_2*nu_mat52_2 - 1)/(nu_mat52_2*a_mat52_2);
            pyb62_2 = o4mat2_2*weight_mat2_2*phi_mat62_2*(nu_mat62_2*nu_mat62_2 - 1)/(nu_mat62_2*a_mat62_2);
            
///
            
            pyb11_3 = o4mat1_3*weight_mat1_3*phi_mat11_3*(nu_mat11_3*nu_mat11_3 - 1)/(nu_mat11_3*a_mat11_3);
            pyb21_3 = o4mat1_3*weight_mat1_3*phi_mat21_3*(nu_mat21_3*nu_mat21_3 - 1)/(nu_mat21_3*a_mat21_3);
            pyb31_3 = o4mat1_3*weight_mat1_3*phi_mat31_3*(nu_mat31_3*nu_mat31_3 - 1)/(nu_mat31_3*a_mat31_3);
            pyb41_3 = o4mat1_3*weight_mat1_3*phi_mat41_3*(nu_mat41_3*nu_mat41_3 - 1)/(nu_mat41_3*a_mat41_3);
            pyb51_3 = o4mat1_3*weight_mat1_3*phi_mat51_3*(nu_mat51_3*nu_mat51_3 - 1)/(nu_mat51_3*a_mat51_3);
            pyb61_3 = o4mat1_3*weight_mat1_3*phi_mat61_3*(nu_mat61_3*nu_mat61_3 - 1)/(nu_mat61_3*a_mat61_3);
            
            pyb12_3 = o4mat2_3*weight_mat2_3*phi_mat12_3*(nu_mat12_3*nu_mat12_3 - 1)/(nu_mat12_3*a_mat12_3);
            pyb22_3 = o4mat2_3*weight_mat2_3*phi_mat22_3*(nu_mat22_3*nu_mat22_3 - 1)/(nu_mat22_3*a_mat22_3);
            pyb32_3 = o4mat2_3*weight_mat2_3*phi_mat32_3*(nu_mat32_3*nu_mat32_3 - 1)/(nu_mat32_3*a_mat32_3);
            pyb42_3 = o4mat2_3*weight_mat2_3*phi_mat42_3*(nu_mat42_3*nu_mat42_3 - 1)/(nu_mat42_3*a_mat42_3);
            pyb52_3 = o4mat2_3*weight_mat2_3*phi_mat52_3*(nu_mat52_3*nu_mat52_3 - 1)/(nu_mat52_3*a_mat52_3);
            pyb62_3 = o4mat2_3*weight_mat2_3*phi_mat62_3*(nu_mat62_3*nu_mat62_3 - 1)/(nu_mat62_3*a_mat62_3);
////
            pyb11_4 = o4mat1_4*weight_mat1_4*phi_mat11_4*(nu_mat11_4*nu_mat11_4 - 1)/(nu_mat11_4*a_mat11_4);
            pyb21_4 = o4mat1_4*weight_mat1_4*phi_mat21_4*(nu_mat21_4*nu_mat21_4 - 1)/(nu_mat21_4*a_mat21_4);
            pyb31_4 = o4mat1_4*weight_mat1_4*phi_mat31_4*(nu_mat31_4*nu_mat31_4 - 1)/(nu_mat31_4*a_mat31_4);
            pyb41_4 = o4mat1_4*weight_mat1_4*phi_mat41_4*(nu_mat41_4*nu_mat41_4 - 1)/(nu_mat41_4*a_mat41_4);
            pyb51_4 = o4mat1_4*weight_mat1_4*phi_mat51_4*(nu_mat51_4*nu_mat51_4 - 1)/(nu_mat51_4*a_mat51_4);
            pyb61_4 = o4mat1_4*weight_mat1_4*phi_mat61_4*(nu_mat61_4*nu_mat61_4 - 1)/(nu_mat61_4*a_mat61_4);
            
            pyb12_4 = o4mat2_4*weight_mat2_4*phi_mat12_4*(nu_mat12_4*nu_mat12_4 - 1)/(nu_mat12_4*a_mat12_4);
            pyb22_4 = o4mat2_4*weight_mat2_4*phi_mat22_4*(nu_mat22_4*nu_mat22_4 - 1)/(nu_mat22_4*a_mat22_4);
            pyb32_4 = o4mat2_4*weight_mat2_4*phi_mat32_4*(nu_mat32_4*nu_mat32_4 - 1)/(nu_mat32_4*a_mat32_4);
            pyb42_4 = o4mat2_4*weight_mat2_4*phi_mat42_4*(nu_mat42_4*nu_mat42_4 - 1)/(nu_mat42_4*a_mat42_4);
            pyb52_4 = o4mat2_4*weight_mat2_4*phi_mat52_4*(nu_mat52_4*nu_mat52_4 - 1)/(nu_mat52_4*a_mat52_4);
            pyb62_4 = o4mat2_4*weight_mat2_4*phi_mat62_4*(nu_mat62_4*nu_mat62_4 - 1)/(nu_mat62_4*a_mat62_4);
            
            
            
            pyw1_1 = o4mat1_1*(phi_mat11_1+phi_mat21_1+phi_mat31_1+phi_mat41_1+phi_mat51_1+phi_mat61_1);
            pyw2_1 = o4mat2_1*(phi_mat12_1+phi_mat22_1+phi_mat32_1+phi_mat42_1+phi_mat52_1+phi_mat62_1);
            
            pyw1_2 = o4mat1_2*(phi_mat11_2+phi_mat21_2+phi_mat31_2+phi_mat41_2+phi_mat51_2+phi_mat61_2);
            pyw2_2 = o4mat2_2*(phi_mat12_2+phi_mat22_2+phi_mat32_2+phi_mat42_2+phi_mat52_2+phi_mat62_2);
            
            pyw1_3 = o4mat1_3*(phi_mat11_3+phi_mat21_3+phi_mat31_3+phi_mat41_3+phi_mat51_3+phi_mat61_3);
            pyw2_3 = o4mat2_3*(phi_mat12_3+phi_mat22_3+phi_mat32_3+phi_mat42_3+phi_mat52_3+phi_mat62_3);
            
            pyw1_4 = o4mat1_4*(phi_mat11_4+phi_mat21_4+phi_mat31_4+phi_mat41_4+phi_mat51_4+phi_mat61_4);
            pyw2_4 = o4mat2_4*(phi_mat12_4+phi_mat22_4+phi_mat32_4+phi_mat42_4+phi_mat52_4+phi_mat62_4);
            
            eta1 = 1;
            eta2 = 1;
            eta3 = 1;
            eta4 = 1;
            
            alpha0 =
                    eta1*eta1*(pyb11_1*pyb11_1 + pyb21_1*pyb21_1 + pyb31_1*pyb31_1 + pyb41_1*pyb41_1 + pyb51_1*pyb51_1 + pyb61_1*pyb61_1 + pyb12_1*pyb12_1 + pyb22_1*pyb22_1 + pyb32_1*pyb32_1 + pyb42_1*pyb42_1 + pyb52_1*pyb52_1 + pyb62_1*pyb62_1)*(x1_k-xf_1)*(x1_k-xf_1) +
                    eta2*eta2*(pyb11_2*pyb11_2 + pyb21_2*pyb21_2 + pyb31_2*pyb31_2 + pyb41_2*pyb41_2 + pyb51_2*pyb51_2 + pyb61_2*pyb61_2 + pyb12_2*pyb12_2 + pyb22_2*pyb22_2 + pyb32_2*pyb32_2 + pyb42_2*pyb42_2 + pyb52_1*pyb52_2 + pyb62_1*pyb62_2)*(x2_k-xf_2)*(x2_k-xf_2) +
                    eta3*eta3*(pyb11_3*pyb11_3 + pyb21_3*pyb21_3 + pyb31_3*pyb31_3 + pyb41_3*pyb41_3 + pyb51_3*pyb51_3 + pyb61_3*pyb61_3 + pyb12_3*pyb12_3 + pyb22_3*pyb22_3 + pyb32_3*pyb32_3 + pyb42_3*pyb42_3 + pyb52_1*pyb52_3 + pyb62_1*pyb62_3)*(x3_k-xf_3)*(x3_k-xf_3) +
                    eta4*eta4*(pyb11_4*pyb11_4 + pyb21_4*pyb21_4 + pyb31_4*pyb31_4 + pyb41_4*pyb41_4 + pyb51_4*pyb51_4 + pyb61_4*pyb61_4 + pyb12_4*pyb12_4 + pyb22_4*pyb22_4 + pyb32_4*pyb32_4 + pyb42_4*pyb42_4 + pyb52_1*pyb52_4 + pyb62_1*pyb62_4)*(x4_k-xf_4)*(x4_k-xf_4)
                    +
                    eta1*eta1*(pya11_1*pya11_1 + pya21_1*pya21_1 + pya31_1*pya31_1 + pya41_1*pya41_1 + pya51_1*pya51_1 + pya61_1*pya61_1 + pya12_1*pya12_1 + pya22_1*pya22_1 + pya32_1*pya32_1 + pya42_1*pya42_1 + pya52_1*pya52_1 + pya62_1*pya62_1)*(x1_k-xf_1)*(x1_k-xf_1) +
                    eta2*eta2*(pya11_2*pya11_2 + pya21_2*pya21_2 + pya31_2*pya31_2 + pya41_2*pya41_2 + pya51_2*pya51_2 + pya61_2*pya61_2 + pya12_2*pya12_2 + pya22_2*pya22_2 + pya32_2*pya32_2 + pya42_2*pya42_2 + pya52_1*pya52_2 + pya62_1*pya62_2)*(x2_k-xf_2)*(x2_k-xf_2) +
                    eta3*eta3*(pya11_3*pya11_3 + pya21_3*pya21_3 + pya31_3*pya31_3 + pya41_3*pya41_3 + pya51_3*pya51_3 + pya61_3*pya61_3 + pya12_3*pya12_3 + pya22_3*pya22_3 + pya32_3*pya32_3 + pya42_3*pya42_3 + pya52_1*pya52_3 + pya62_1*pya62_3)*(x3_k-xf_3)*(x3_k-xf_3) +
                    eta4*eta4*(pya11_4*pya11_4 + pya21_4*pya21_4 + pya31_4*pya31_4 + pya41_4*pya41_4 + pya51_4*pya51_4 + pya61_4*pya61_4 + pya12_4*pya12_4 + pya22_4*pya22_4 + pya32_4*pya32_4 + pya42_4*pya42_4 + pya52_1*pya52_4 + pya62_1*pya62_4)*(x4_k-xf_4)*(x4_k-xf_4)
                    +
                    eta1*eta1*(pym11_1*pym11_1 + pym21_1*pym21_1 + pym31_1*pym31_1 + pym41_1*pym41_1 + pym51_1*pym51_1 + pym61_1*pym61_1 + pym12_1*pym12_1 + pym22_1*pym22_1 + pym32_1*pym32_1 + pym42_1*pym42_1 + pym52_1*pym52_1 + pym62_1*pym62_1)*(x1_k-xf_1)*(x1_k-xf_1) +
                    eta2*eta2*(pym11_2*pym11_2 + pym21_2*pym21_2 + pym31_2*pym31_2 + pym41_2*pym41_2 + pym51_2*pym51_2 + pym61_2*pym61_2 + pym12_2*pym12_2 + pym22_2*pym22_2 + pym32_2*pym32_2 + pym42_2*pym42_2 + pym52_1*pym52_2 + pym62_1*pym62_2)*(x2_k-xf_2)*(x2_k-xf_2) +
                    eta3*eta3*(pym11_3*pym11_3 + pym21_3*pym21_3 + pym31_3*pym31_3 + pym41_3*pym41_3 + pym51_3*pym51_3 + pym61_3*pym61_3 + pym12_3*pym12_3 + pym22_3*pym22_3 + pym32_3*pym32_3 + pym42_3*pym42_3 + pym52_1*pym52_3 + pym62_1*pym62_3)*(x3_k-xf_3)*(x3_k-xf_3) +
                    eta4*eta4*(pym11_4*pym11_4 + pym21_4*pym21_4 + pym31_4*pym31_4 + pym41_4*pym41_4 + pym51_4*pym51_4 + pym61_4*pym61_4 + pym12_4*pym12_4 + pym22_4*pym22_4 + pym32_4*pym32_4 + pym42_4*pym42_4 + pym52_1*pym52_4 + pym62_1*pym62_4)*(x4_k-xf_4)*(x4_k-xf_4)
                    +
                    eta1*eta1*(pys11_1*pys11_1 + pys21_1*pys21_1 + pys31_1*pys31_1 + pys41_1*pys41_1 + pys51_1*pys51_1 + pys61_1*pys61_1 + pys12_1*pys12_1 + pys22_1*pys22_1 + pys32_1*pys32_1 + pys42_1*pys42_1 + pys52_1*pys52_1 + pys62_1*pys62_1)*(x1_k-xf_1)*(x1_k-xf_1) +
                    eta2*eta2*(pys11_2*pys11_2 + pys21_2*pys21_2 + pys31_2*pys31_2 + pys41_2*pys41_2 + pys51_2*pys51_2 + pys61_2*pys61_2 + pys12_2*pys12_2 + pys22_2*pys22_2 + pys32_2*pys32_2 + pys42_2*pys42_2 + pys52_1*pys52_2 + pys62_1*pys62_2)*(x2_k-xf_2)*(x2_k-xf_2) +
                    eta3*eta3*(pys11_3*pys11_3 + pys21_3*pys21_3 + pys31_3*pys31_3 + pys41_3*pys41_3 + pys51_3*pys51_3 + pys61_3*pys61_3 + pys12_3*pys12_3 + pys22_3*pys22_3 + pys32_3*pys32_3 + pys42_3*pys42_3 + pys52_1*pys52_3 + pys62_1*pys62_3)*(x3_k-xf_3)*(x3_k-xf_3) +
                    eta4*eta4*(pys11_4*pys11_4 + pys21_4*pys21_4 + pys31_4*pys31_4 + pys41_4*pys41_4 + pys51_4*pys51_4 + pys61_4*pys61_4 + pys12_4*pys12_4 + pys22_4*pys22_4 + pys32_4*pys32_4 + pys42_4*pys42_4 + pys52_1*pys52_4 + pys62_1*pys62_4)*(x4_k-xf_4)*(x4_k-xf_4)
                    +
                    eta1*eta1*(pyw1_1*pyw1_1 + pyw2_1*pyw2_1)*(x1_k-xf_1)*(x1_k-xf_1) + eta2*eta2*(pyw1_2*pyw1_2 + pyw2_2*pyw2_2)*(x2_k-xf_2)*(x2_k-xf_2) + eta3*eta3*(pyw1_3*pyw1_3 + pyw2_3*pyw2_3)*(x3_k-xf_3)*(x3_k-xf_3) + eta4*eta4*(pyw1_4*pyw1_4 + pyw2_4*pyw2_4)*(x4_k-xf_4)*(x4_k-xf_4);
            
            alpha1 =
                    eta1*eta1*eta1*(pyb11_1*pyb11_1 + pyb21_1*pyb21_1 + pyb31_1*pyb31_1 + pyb41_1*pyb41_1 + pyb51_1*pyb51_1 + pyb61_1*pyb61_1 + pyb12_1*pyb12_1 + pyb22_1*pyb22_1 + pyb32_1*pyb32_1 + pyb42_1*pyb42_1 + pyb52_1*pyb52_1 + pyb62_1*pyb62_1)*(pyb11_1*pyb11_1 + pyb21_1*pyb21_1 + pyb31_1*pyb31_1 + pyb41_1*pyb41_1 + pyb51_1*pyb51_1 + pyb61_1*pyb61_1 + pyb12_1*pyb12_1 + pyb22_1*pyb22_1 + pyb32_1*pyb32_1 + pyb42_1*pyb42_1 + pyb52_1*pyb52_1 + pyb62_1*pyb62_1)*(x1_k-xf_1)*(x1_k-xf_1) +
                    eta2*eta2*eta2*(pyb11_2*pyb11_2 + pyb21_2*pyb21_2 + pyb31_2*pyb31_2 + pyb41_2*pyb41_2 + pyb51_2*pyb51_2 + pyb61_2*pyb61_2 + pyb12_2*pyb12_2 + pyb22_2*pyb22_2 + pyb32_2*pyb32_2 + pyb42_2*pyb42_2 + pyb52_1*pyb52_2 + pyb62_1*pyb62_2)*(pyb11_2*pyb11_2 + pyb21_2*pyb21_2 + pyb31_2*pyb31_2 + pyb41_2*pyb41_2 + pyb51_2*pyb51_2 + pyb61_2*pyb61_2 + pyb12_2*pyb12_2 + pyb22_2*pyb22_2 + pyb32_2*pyb32_2 + pyb42_2*pyb42_2 + pyb52_1*pyb52_2 + pyb62_1*pyb62_2)*(x2_k-xf_2)*(x2_k-xf_2) +
                    eta3*eta3*eta3*(pyb11_3*pyb11_3 + pyb21_3*pyb21_3 + pyb31_3*pyb31_3 + pyb41_3*pyb41_3 + pyb51_3*pyb51_3 + pyb61_3*pyb61_3 + pyb12_3*pyb12_3 + pyb22_3*pyb22_3 + pyb32_3*pyb32_3 + pyb42_3*pyb42_3 + pyb52_1*pyb52_3 + pyb62_1*pyb62_3)*(pyb11_3*pyb11_3 + pyb21_3*pyb21_3 + pyb31_3*pyb31_3 + pyb41_3*pyb41_3 + pyb51_3*pyb51_3 + pyb61_3*pyb61_3 + pyb12_3*pyb12_3 + pyb22_3*pyb22_3 + pyb32_3*pyb32_3 + pyb42_3*pyb42_3 + pyb52_1*pyb52_3 + pyb62_1*pyb62_3)*(x3_k-xf_3)*(x3_k-xf_3) +
                    eta4*eta4*eta4*(pyb11_4*pyb11_4 + pyb21_4*pyb21_4 + pyb31_4*pyb31_4 + pyb41_4*pyb41_4 + pyb51_4*pyb51_4 + pyb61_4*pyb61_4 + pyb12_4*pyb12_4 + pyb22_4*pyb22_4 + pyb32_4*pyb32_4 + pyb42_4*pyb42_4 + pyb52_1*pyb52_4 + pyb62_1*pyb62_4)*(pyb11_4*pyb11_4 + pyb21_4*pyb21_4 + pyb31_4*pyb31_4 + pyb41_4*pyb41_4 + pyb51_4*pyb51_4 + pyb61_4*pyb61_4 + pyb12_4*pyb12_4 + pyb22_4*pyb22_4 + pyb32_4*pyb32_4 + pyb42_4*pyb42_4 + pyb52_1*pyb52_4 + pyb62_1*pyb62_4)*(x4_k-xf_4)*(x4_k-xf_4)
                    +
                    eta1*eta1*eta1*(pya11_1*pya11_1 + pya21_1*pya21_1 + pya31_1*pya31_1 + pya41_1*pya41_1 + pya51_1*pya51_1 + pya61_1*pya61_1 + pya12_1*pya12_1 + pya22_1*pya22_1 + pya32_1*pya32_1 + pya42_1*pya42_1 + pya52_1*pya52_1 + pya62_1*pya62_1)*(pya11_1*pya11_1 + pya21_1*pya21_1 + pya31_1*pya31_1 + pya41_1*pya41_1 + pya51_1*pya51_1 + pya61_1*pya61_1 + pya12_1*pya12_1 + pya22_1*pya22_1 + pya32_1*pya32_1 + pya42_1*pya42_1 + pya52_1*pya52_1 + pya62_1*pya62_1)*(x1_k-xf_1)*(x1_k-xf_1) +
                    eta2*eta2*eta2*(pya11_2*pya11_2 + pya21_2*pya21_2 + pya31_2*pya31_2 + pya41_2*pya41_2 + pya51_2*pya51_2 + pya61_2*pya61_2 + pya12_2*pya12_2 + pya22_2*pya22_2 + pya32_2*pya32_2 + pya42_2*pya42_2 + pya52_1*pya52_2 + pya62_1*pya62_2)*(pya11_2*pya11_2 + pya21_2*pya21_2 + pya31_2*pya31_2 + pya41_2*pya41_2 + pya51_2*pya51_2 + pya61_2*pya61_2 + pya12_2*pya12_2 + pya22_2*pya22_2 + pya32_2*pya32_2 + pya42_2*pya42_2 + pya52_1*pya52_2 + pya62_1*pya62_2)*(x2_k-xf_2)*(x2_k-xf_2) +
                    eta3*eta3*eta3*(pya11_3*pya11_3 + pya21_3*pya21_3 + pya31_3*pya31_3 + pya41_3*pya41_3 + pya51_3*pya51_3 + pya61_3*pya61_3 + pya12_3*pya12_3 + pya22_3*pya22_3 + pya32_3*pya32_3 + pya42_3*pya42_3 + pya52_1*pya52_3 + pya62_1*pya62_3)*(pya11_3*pya11_3 + pya21_3*pya21_3 + pya31_3*pya31_3 + pya41_3*pya41_3 + pya51_3*pya51_3 + pya61_3*pya61_3 + pya12_3*pya12_3 + pya22_3*pya22_3 + pya32_3*pya32_3 + pya42_3*pya42_3 + pya52_1*pya52_3 + pya62_1*pya62_3)*(x3_k-xf_3)*(x3_k-xf_3) +
                    eta4*eta4*eta4*(pya11_4*pya11_4 + pya21_4*pya21_4 + pya31_4*pya31_4 + pya41_4*pya41_4 + pya51_4*pya51_4 + pya61_4*pya61_4 + pya12_4*pya12_4 + pya22_4*pya22_4 + pya32_4*pya32_4 + pya42_4*pya42_4 + pya52_1*pya52_4 + pya62_1*pya62_4)*(pya11_4*pya11_4 + pya21_4*pya21_4 + pya31_4*pya31_4 + pya41_4*pya41_4 + pya51_4*pya51_4 + pya61_4*pya61_4 + pya12_4*pya12_4 + pya22_4*pya22_4 + pya32_4*pya32_4 + pya42_4*pya42_4 + pya52_1*pya52_4 + pya62_1*pya62_4)*(x4_k-xf_4)*(x4_k-xf_4)
                    +
                    eta1*eta1*eta1*(pym11_1*pym11_1 + pym21_1*pym21_1 + pym31_1*pym31_1 + pym41_1*pym41_1 + pym51_1*pym51_1 + pym61_1*pym61_1 + pym12_1*pym12_1 + pym22_1*pym22_1 + pym32_1*pym32_1 + pym42_1*pym42_1 + pym52_1*pym52_1 + pym62_1*pym62_1)*(pym11_1*pym11_1 + pym21_1*pym21_1 + pym31_1*pym31_1 + pym41_1*pym41_1 + pym51_1*pym51_1 + pym61_1*pym61_1 + pym12_1*pym12_1 + pym22_1*pym22_1 + pym32_1*pym32_1 + pym42_1*pym42_1 + pym52_1*pym52_1 + pym62_1*pym62_1)*(x1_k-xf_1)*(x1_k-xf_1) +
                    eta2*eta2*eta2*(pym11_2*pym11_2 + pym21_2*pym21_2 + pym31_2*pym31_2 + pym41_2*pym41_2 + pym51_2*pym51_2 + pym61_2*pym61_2 + pym12_2*pym12_2 + pym22_2*pym22_2 + pym32_2*pym32_2 + pym42_2*pym42_2 + pym52_1*pym52_2 + pym62_1*pym62_2)*(pym11_2*pym11_2 + pym21_2*pym21_2 + pym31_2*pym31_2 + pym41_2*pym41_2 + pym51_2*pym51_2 + pym61_2*pym61_2 + pym12_2*pym12_2 + pym22_2*pym22_2 + pym32_2*pym32_2 + pym42_2*pym42_2 + pym52_1*pym52_2 + pym62_1*pym62_2)*(x2_k-xf_2)*(x2_k-xf_2) +
                    eta3*eta3*eta3*(pym11_3*pym11_3 + pym21_3*pym21_3 + pym31_3*pym31_3 + pym41_3*pym41_3 + pym51_3*pym51_3 + pym61_3*pym61_3 + pym12_3*pym12_3 + pym22_3*pym22_3 + pym32_3*pym32_3 + pym42_3*pym42_3 + pym52_1*pym52_3 + pym62_1*pym62_3)*(pym11_3*pym11_3 + pym21_3*pym21_3 + pym31_3*pym31_3 + pym41_3*pym41_3 + pym51_3*pym51_3 + pym61_3*pym61_3 + pym12_3*pym12_3 + pym22_3*pym22_3 + pym32_3*pym32_3 + pym42_3*pym42_3 + pym52_1*pym52_3 + pym62_1*pym62_3)*(x3_k-xf_3)*(x3_k-xf_3) +
                    eta4*eta4*eta4*(pym11_4*pym11_4 + pym21_4*pym21_4 + pym31_4*pym31_4 + pym41_4*pym41_4 + pym51_4*pym51_4 + pym61_4*pym61_4 + pym12_4*pym12_4 + pym22_4*pym22_4 + pym32_4*pym32_4 + pym42_4*pym42_4 + pym52_1*pym52_4 + pym62_1*pym62_4)*(pym11_4*pym11_4 + pym21_4*pym21_4 + pym31_4*pym31_4 + pym41_4*pym41_4 + pym51_4*pym51_4 + pym61_4*pym61_4 + pym12_4*pym12_4 + pym22_4*pym22_4 + pym32_4*pym32_4 + pym42_4*pym42_4 + pym52_1*pym52_4 + pym62_1*pym62_4)*(x4_k-xf_4)*(x4_k-xf_4)
                    +
                    eta1*eta1*eta1*(pys11_1*pys11_1 + pys21_1*pys21_1 + pys31_1*pys31_1 + pys41_1*pys41_1 + pys51_1*pys51_1 + pys61_1*pys61_1 + pys12_1*pys12_1 + pys22_1*pys22_1 + pys32_1*pys32_1 + pys42_1*pys42_1 + pys52_1*pys52_1 + pys62_1*pys62_1)*(pys11_1*pys11_1 + pys21_1*pys21_1 + pys31_1*pys31_1 + pys41_1*pys41_1 + pys51_1*pys51_1 + pys61_1*pys61_1 + pys12_1*pys12_1 + pys22_1*pys22_1 + pys32_1*pys32_1 + pys42_1*pys42_1 + pys52_1*pys52_1 + pys62_1*pys62_1)*(x1_k-xf_1)*(x1_k-xf_1)+
                    eta2*eta2*eta2*(pys11_2*pys11_2 + pys21_2*pys21_2 + pys31_2*pys31_2 + pys41_2*pys41_2 + pys51_2*pys51_2 + pys61_2*pys61_2 + pys12_2*pys12_2 + pys22_2*pys22_2 + pys32_2*pys32_2 + pys42_2*pys42_2 + pys52_1*pys52_2 + pys62_1*pys62_2)*(pys11_2*pys11_2 + pys21_2*pys21_2 + pys31_2*pys31_2 + pys41_2*pys41_2 + pys51_2*pys51_2 + pys61_2*pys61_2 + pys12_2*pys12_2 + pys22_2*pys22_2 + pys32_2*pys32_2 + pys42_2*pys42_2 + pys52_1*pys52_2 + pys62_1*pys62_2)*(x2_k-xf_2)*(x2_k-xf_2) +
                    eta3*eta3*eta3*(pys11_3*pys11_3 + pys21_3*pys21_3 + pys31_3*pys31_3 + pys41_3*pys41_3 + pys51_3*pys51_3 + pys61_3*pys61_3 + pys12_3*pys12_3 + pys22_3*pys22_3 + pys32_3*pys32_3 + pys42_3*pys42_3 + pys52_1*pys52_3 + pys62_1*pys62_3)*(pys11_3*pys11_3 + pys21_3*pys21_3 + pys31_3*pys31_3 + pys41_3*pys41_3 + pys51_3*pys51_3 + pys61_3*pys61_3 + pys12_3*pys12_3 + pys22_3*pys22_3 + pys32_3*pys32_3 + pys42_3*pys42_3 + pys52_1*pys52_3 + pys62_1*pys62_3)*(x3_k-xf_3)*(x3_k-xf_3) +
                    eta4*eta4*eta4*(pys11_4*pys11_4 + pys21_4*pys21_4 + pys31_4*pys31_4 + pys41_4*pys41_4 + pys51_4*pys51_4 + pys61_4*pys61_4 + pys12_4*pys12_4 + pys22_4*pys22_4 + pys32_4*pys32_4 + pys42_4*pys42_4 + pys52_1*pys52_4 + pys62_1*pys62_4)*(pys11_4*pys11_4 + pys21_4*pys21_4 + pys31_4*pys31_4 + pys41_4*pys41_4 + pys51_4*pys51_4 + pys61_4*pys61_4 + pys12_4*pys12_4 + pys22_4*pys22_4 + pys32_4*pys32_4 + pys42_4*pys42_4 + pys52_1*pys52_4 + pys62_1*pys62_4)*(x4_k-xf_4)*(x4_k-xf_4)
                    +
                    eta1*eta1*eta1*(pyw1_1*pyw1_1 + pyw2_1*pyw2_1)*(pyw1_1*pyw1_1 + pyw2_1*pyw2_1)*(x1_k-xf_1)*(x1_k-xf_1) + eta2*eta2*eta2*(pyw1_2*pyw1_2 + pyw2_2*pyw2_2)*(pyw1_2*pyw1_2 + pyw2_2*pyw2_2)*(x2_k-xf_2)*(x2_k-xf_2) + eta3*eta3*eta3*(pyw1_3*pyw1_3 + pyw2_3*pyw2_3)*(pyw1_3*pyw1_3 + pyw2_3*pyw2_3)*(x3_k-xf_3)*(x3_k-xf_3) + eta4*eta4*eta4*(pyw1_4*pyw1_4 + pyw2_4*pyw2_4)*(pyw1_4*pyw1_4 + pyw2_4*pyw2_4)*(x4_k-xf_4)*(x4_k-xf_4);
            
            debpar1[k] = alpha0;
            debpar2[k] =  alpha1;
            
            rho = alpha0/alpha1;
            if (abs(rho)<100){
                rho = rho;
            } else {
                rho = 0;
            }
            
            xtilde_1 = x1_k - xf_1;
            xtilde_2 = x2_k - xf_2;
            xtilde_3 = x3_k - xf_3;
            xtilde_4 = x4_k - xf_4;
            
            
            // UPDATING PARAMETERS
            
            m_mat11_1 += rho*xtilde_1*pym11_1;
            m_mat21_1 += rho*xtilde_1*pym21_1;
            m_mat31_1 += rho*xtilde_1*pym31_1;
            m_mat41_1 += rho*xtilde_1*pym41_1;
            m_mat51_1 += rho*xtilde_1*pym51_1;
            m_mat61_1 += rho*xtilde_1*pym61_1;
            
            m_mat12_1 += rho*xtilde_1*pym12_1;
            m_mat22_1 += rho*xtilde_1*pym22_1;
            m_mat32_1 += rho*xtilde_1*pym32_1;
            m_mat42_1 += rho*xtilde_1*pym42_1;
            m_mat52_1 += rho*xtilde_1*pym52_1;
            m_mat62_1 += rho*xtilde_1*pym62_1;
            //
            m_mat11_2 += rho*xtilde_2*pym11_2;
            m_mat21_2 += rho*xtilde_2*pym21_2;
            m_mat31_2 += rho*xtilde_2*pym31_2;
            m_mat41_2 += rho*xtilde_2*pym41_2;
            m_mat51_2 += rho*xtilde_2*pym51_2;
            m_mat61_2 += rho*xtilde_2*pym61_2;
            
            m_mat12_2 += rho*xtilde_2*pym12_2;
            m_mat22_2 += rho*xtilde_2*pym22_2;
            m_mat32_2 += rho*xtilde_2*pym32_2;
            m_mat42_2 += rho*xtilde_2*pym42_2;
            m_mat52_2 += rho*xtilde_2*pym52_2;
            m_mat62_2 += rho*xtilde_2*pym62_2;
            ///
            m_mat11_3 += rho*xtilde_3*pym11_3;
            m_mat21_3 += rho*xtilde_3*pym21_3;
            m_mat31_3 += rho*xtilde_3*pym31_3;
            m_mat41_3 += rho*xtilde_3*pym41_3;
            m_mat51_3 += rho*xtilde_3*pym51_3;
            m_mat61_3 += rho*xtilde_3*pym61_3;
            
            m_mat12_3 += rho*xtilde_3*pym12_3;
            m_mat22_3 += rho*xtilde_3*pym22_3;
            m_mat32_3 += rho*xtilde_3*pym32_3;
            m_mat42_3 += rho*xtilde_3*pym42_3;
            m_mat52_3 += rho*xtilde_3*pym52_3;
            m_mat62_3 += rho*xtilde_3*pym62_3;
            ////
            m_mat11_4 += rho*xtilde_4*pym11_4;
            m_mat21_4 += rho*xtilde_4*pym21_4;
            m_mat31_4 += rho*xtilde_4*pym31_4;
            m_mat41_4 += rho*xtilde_4*pym41_4;
            m_mat51_4 += rho*xtilde_4*pym51_4;
            m_mat61_4 += rho*xtilde_4*pym61_4;
            
            m_mat12_4 += rho*xtilde_4*pym12_4;
            m_mat22_4 += rho*xtilde_4*pym22_4;
            m_mat32_4 += rho*xtilde_4*pym32_4;
            m_mat42_4 += rho*xtilde_4*pym42_4;
            m_mat52_4 += rho*xtilde_4*pym52_4;
            m_mat62_4 += rho*xtilde_4*pym62_4;
            
            
            sigma_mat11_1 += rho*xtilde_1*pys11_1; // CORRECTED
            sigma_mat21_1 += rho*xtilde_1*pys21_1;
            sigma_mat31_1 += rho*xtilde_1*pys31_1;
            sigma_mat41_1 += rho*xtilde_1*pys41_1;
            sigma_mat51_1 += rho*xtilde_1*pys51_1;
            sigma_mat61_1 += rho*xtilde_1*pys61_1;
            
            sigma_mat12_1 += rho*xtilde_1*pys12_1;
            sigma_mat22_1 += rho*xtilde_1*pys22_1;
            sigma_mat32_1 += rho*xtilde_1*pys32_1;
            sigma_mat42_1 += rho*xtilde_1*pys42_1;
            sigma_mat52_1 += rho*xtilde_1*pys52_1;
            sigma_mat62_1 += rho*xtilde_1*pys62_1;
            //
            sigma_mat11_2 += rho*xtilde_2*pys11_2; // CORRECTED
            sigma_mat21_2 += rho*xtilde_2*pys21_2;
            sigma_mat31_2 += rho*xtilde_2*pys31_2;
            sigma_mat41_2 += rho*xtilde_2*pys41_2;
            sigma_mat51_2 += rho*xtilde_2*pys51_2;
            sigma_mat61_2 += rho*xtilde_2*pys61_2;
            
            sigma_mat12_2 += rho*xtilde_2*pys12_2;
            sigma_mat22_2 += rho*xtilde_2*pys22_2;
            sigma_mat32_2 += rho*xtilde_2*pys32_2;
            sigma_mat42_2 += rho*xtilde_2*pys42_2;
            sigma_mat52_2 += rho*xtilde_2*pys52_2;
            sigma_mat62_2 += rho*xtilde_2*pys62_2;
///
            sigma_mat11_3 += rho*xtilde_3*pys11_3;  // CORRECTED
            sigma_mat21_3 += rho*xtilde_3*pys21_3;
            sigma_mat31_3 += rho*xtilde_3*pys31_3;
            sigma_mat41_3 += rho*xtilde_3*pys41_3;
            sigma_mat51_3 += rho*xtilde_3*pys51_3;
            sigma_mat61_3 += rho*xtilde_3*pys61_3;
            
            sigma_mat12_3 += rho*xtilde_3*pys12_3;
            sigma_mat22_3 += rho*xtilde_3*pys22_3;
            sigma_mat32_3 += rho*xtilde_3*pys32_3;
            sigma_mat42_3 += rho*xtilde_3*pys42_3;
            sigma_mat52_3 += rho*xtilde_3*pys52_3;
            sigma_mat62_3 += rho*xtilde_3*pys62_3;
/////
            sigma_mat11_4 += rho*xtilde_4*pys11_4; // CORRECTED
            sigma_mat21_4 += rho*xtilde_4*pys21_4; // CORRECTED
            sigma_mat31_4 += rho*xtilde_4*pys31_4; // CORRECTED
            sigma_mat41_4 += rho*xtilde_4*pys41_4; // CORRECTED
            sigma_mat51_4 += rho*xtilde_4*pys51_4; // CORRECTED
            sigma_mat61_4 += rho*xtilde_4*pys61_4; // CORRECTED
            
            sigma_mat12_4 += rho*xtilde_4*pys12_4;
            sigma_mat22_4 += rho*xtilde_4*pys22_4;
            sigma_mat32_4 += rho*xtilde_4*pys32_4;
            sigma_mat42_4 += rho*xtilde_4*pys42_4;
            sigma_mat52_4 += rho*xtilde_4*pys52_4;
            sigma_mat62_4 += rho*xtilde_4*pys62_4;
            
            a_mat11_1 += rho*xtilde_1*pya11_1;// CORRECTED
            a_mat21_1 += rho*xtilde_1*pya21_1;
            a_mat31_1 += rho*xtilde_1*pya31_1;
            a_mat41_1 += rho*xtilde_1*pya41_1;
            a_mat51_1 += rho*xtilde_1*pya51_1;
            a_mat61_1 += rho*xtilde_1*pya61_1;
            
            a_mat12_1 += rho*xtilde_1*pya12_1;
            a_mat22_1 += rho*xtilde_1*pya22_1;
            a_mat32_1 += rho*xtilde_1*pya32_1;
            a_mat42_1 += rho*xtilde_1*pya42_1;
            a_mat52_1 += rho*xtilde_1*pya52_1;
            a_mat62_1 += rho*xtilde_1*pya62_1;
            //
            a_mat11_2 += rho*xtilde_2*pya11_2;// CORRECTED
            a_mat21_2 += rho*xtilde_2*pya21_2;
            a_mat31_2 += rho*xtilde_2*pya31_2;
            a_mat41_2 += rho*xtilde_2*pya41_2;
            a_mat51_2 += rho*xtilde_2*pya51_2;
            a_mat61_2 += rho*xtilde_2*pya61_2;
            
            a_mat12_2 += rho*xtilde_2*pya12_2;
            a_mat22_2 += rho*xtilde_2*pya22_2;
            a_mat32_2 += rho*xtilde_2*pya32_2;
            a_mat42_2 += rho*xtilde_2*pya42_2;
            a_mat52_2 += rho*xtilde_2*pya52_2;
            a_mat62_2 += rho*xtilde_2*pya62_2;
///
            a_mat11_3 += rho*xtilde_3*pya11_3;// CORRECTED
            a_mat21_3 += rho*xtilde_3*pya21_3;
            a_mat31_3 += rho*xtilde_3*pya31_3;
            a_mat41_3 += rho*xtilde_3*pya41_3;
            a_mat51_3 += rho*xtilde_3*pya51_3;
            a_mat61_3 += rho*xtilde_3*pya61_3;
            
            a_mat12_3 += rho*xtilde_3*pya12_3;
            a_mat22_3 += rho*xtilde_3*pya22_3;
            a_mat32_3 += rho*xtilde_3*pya32_3;
            a_mat42_3 += rho*xtilde_3*pya42_3;
            a_mat52_3 += rho*xtilde_3*pya52_3;
            a_mat62_3 += rho*xtilde_3*pya62_3;
/////
            a_mat11_4 += rho*xtilde_4*pya11_4;// CORRECTED
            a_mat21_4 += rho*xtilde_4*pya21_4;// CORRECTED
            a_mat31_4 += rho*xtilde_4*pya31_4;// CORRECTED
            a_mat41_4 += rho*xtilde_4*pya41_4;// CORRECTED
            a_mat51_4 += rho*xtilde_4*pya51_4;// CORRECTED
            a_mat61_4 += rho*xtilde_4*pya61_4;// CORRECTED
            
            a_mat12_4 += rho*xtilde_4*pya12_4;
            a_mat22_4 += rho*xtilde_4*pya22_4;
            a_mat32_4 += rho*xtilde_4*pya32_4;
            a_mat42_4 += rho*xtilde_4*pya42_4;
            a_mat52_4 += rho*xtilde_4*pya52_4;
            a_mat62_4 += rho*xtilde_4*pya62_4;
            
            
            b_mat11_1 += rho*xtilde_1*pyb11_1;// CORRECTED
            b_mat21_1 += rho*xtilde_1*pyb21_1;
            b_mat31_1 += rho*xtilde_1*pyb31_1;
            b_mat41_1 += rho*xtilde_1*pyb41_1;
            b_mat51_1 += rho*xtilde_1*pyb51_1;
            b_mat61_1 += rho*xtilde_1*pyb61_1;
            
            b_mat12_1 += rho*xtilde_1*pyb12_1;
            b_mat22_1 += rho*xtilde_1*pyb22_1;
            b_mat32_1 += rho*xtilde_1*pyb32_1;
            b_mat42_1 += rho*xtilde_1*pyb42_1;
            b_mat52_1 += rho*xtilde_1*pyb52_1;
            b_mat62_1 += rho*xtilde_1*pyb62_1;
            //
            b_mat11_2 += rho*xtilde_2*pyb11_2;// CORRECTED
            b_mat21_2 += rho*xtilde_2*pyb21_2;
            b_mat31_2 += rho*xtilde_2*pyb31_2;
            b_mat41_2 += rho*xtilde_2*pyb41_2;
            b_mat51_2 += rho*xtilde_2*pyb51_2;
            b_mat61_2 += rho*xtilde_2*pyb61_2;
            
            b_mat12_2 += rho*xtilde_2*pyb12_2;
            b_mat22_2 += rho*xtilde_2*pyb22_2;
            b_mat32_2 += rho*xtilde_2*pyb32_2;
            b_mat42_2 += rho*xtilde_2*pyb42_2;
            b_mat52_2 += rho*xtilde_2*pyb52_2;
            b_mat62_2 += rho*xtilde_2*pyb62_2;
///
            b_mat11_3 += rho*xtilde_3*pyb11_3;// CORRECTED
            b_mat21_3 += rho*xtilde_3*pyb21_3;
            b_mat31_3 += rho*xtilde_3*pyb31_3;
            b_mat41_3 += rho*xtilde_3*pyb41_3;
            b_mat51_3 += rho*xtilde_3*pyb51_3;
            b_mat61_3 += rho*xtilde_3*pyb61_3;
            
            b_mat12_3 += rho*xtilde_3*pyb12_3;
            b_mat22_3 += rho*xtilde_3*pyb22_3;
            b_mat32_3 += rho*xtilde_3*pyb32_3;
            b_mat42_3 += rho*xtilde_3*pyb42_3;
            b_mat52_3 += rho*xtilde_3*pyb52_3;
            b_mat62_3 += rho*xtilde_3*pyb62_3;
/////
            b_mat11_4 += rho*xtilde_4*pyb11_4;// CORRECTED
            b_mat21_4 += rho*xtilde_4*pyb21_4;// CORRECTED
            b_mat31_4 += rho*xtilde_4*pyb31_4;// CORRECTED
            b_mat41_4 += rho*xtilde_4*pyb41_4;// CORRECTED
            b_mat51_4 += rho*xtilde_4*pyb51_4;// CORRECTED
            b_mat61_4 += rho*xtilde_4*pyb61_4;// CORRECTED
            
            b_mat12_4 += rho*xtilde_4*pyb12_4;
            b_mat22_4 += rho*xtilde_4*pyb22_4;
            b_mat32_4 += rho*xtilde_4*pyb32_4;
            b_mat42_4 += rho*xtilde_4*pyb42_4;
            b_mat52_4 += rho*xtilde_4*pyb52_4;
            b_mat62_4 += rho*xtilde_4*pyb62_4;
            
            weight_mat1_1 += rho*xtilde_1*pyw1_1;
            weight_mat2_1 += rho*xtilde_1*pyw2_1;
            
            weight_mat1_2 += rho*xtilde_2*pyw1_2;
            weight_mat2_2 += rho*xtilde_2*pyw2_2;
            
            weight_mat1_3 += rho*xtilde_3*pyw1_3;
            weight_mat2_3 += rho*xtilde_3*pyw2_3;
            
            weight_mat1_4 += rho*xtilde_4*pyw1_4;
            weight_mat2_4 += rho*xtilde_4*pyw2_4;
        } else {
            x1f[k] = 0;
            x2f[k] = 0;
            x3f[k] = 0;
            x4f[k] = 0;
            debpar1[k] = 0;
            debpar2[k] =  0;
        }
        // CONTROL
        
        x1_k = x1[k];
        x2_k = x2[k];
        x3_k = x3[k];
        x4_k = x4[k];
        
        x1_k1 = x1[k-1];
        x2_k1 = x2[k-1];
        x3_k1 = x3[k-1];
        x4_k1 = x4[k-1];
        
        u1_k = u1[k-1];
        u2_k = u2[k-1];
        
        a1_k = x3_k + (3*T/2)*((-r1/l1)*(x3_k)+(-k1/l1)*x2_k) - (T/2)*((-r1/l1)*(x3_k1)+(-k1/l1)*x2_k1 + (1/l1)*u1_k);
        
        omegak1 = x2_k + (1.5*T)*((k1/j1)*x3_k + (k2/j2)*x4_k) - (0.5*T)*((k1/j1)*x3_k1 + (k2/j2)*x4_k1);
        theta_k2 = x1_k + (3*T/2)*omegak1 - (T/2)*x2_k1 + T*x2_k;
        
        alpha_k = x2_k + T*((k1/j1)*x3_k + (k2/j2)*x4_k) - 0.5*T*((k1/j1)*x3_k1 + (k2/j2)*x4_k1) +
                1.5*T*(k1/j1)*(x3_k + 1.5*T*((-r1/l1)*x3_k + (-k1/l1)*x2_k) - 0.5*T*((-r1/l1)*x3_k1 + (-k1/l1)*x2_k1 + (1/l1)*u1_k)) +
                1.5*T*(k2/j2)*(x4_k + 1.5*T*((-r2/l2)*x4_k + (-k2/l2)*x2_k) - 0.5*T*((-r2/l2)*x4_k1 + (-k2/l2)*x2_k1 + (1/l2)*u2_k));
        
        
        a11 = alpha_k - xd2_2 + lambda * (omegak1 - xd2_1);
        a21 = a1_k - xd3_1 + lambda * (x3_k - xd3_0);
        
        
        
        b11 = (9/4)*(T*T)*k1/(j1*l1);
        b12 = (9/4)*(T*T)*k2/(j2*l2);
        b21 = 1.5*(T/l1);
        b22 = 0;
        
        
        x3hat = a1_k + b21 * u1_k + b22 * u2_k;
        x3hat1 = x3_k;
        
        
        x1hat = alpha_k + b11 * u1_k + b12 * u2_k;
        x1hat1 = x2_k;
        
        
        stilde1 = x1hat - xd2_2 + lambda * (x1hat1 - xd2_1) ;
        stilde2 = x3hat - xd3_1 + lambda * (x3hat1 - xd3_0) ;
        
        phih1 = stilde1 - a11;
        phih2 = stilde2 - a21;
        phil1 = -stilde1 - a11;
        phil2 = -stilde2 - a21;
        
        phim1 = (1 - alpha) * phil1 + alpha * phih1;
        phim2 = (1 - alpha) * phil2 + alpha * phih2;
        
        detb = b11 * b22 - b12 * b21;
        
        u1[k] = (b22 * phim1 - b12 * phim2) / detb;
        u2[k] = (b11 * phim2 - b21 * phim1) / detb;
        
        
//  if (abs(u1[k])>thres){
//      u1[k] = thres*abs(u1[k])/u1[k];
//  }
//
//   if (abs(u2[k])>thres){
//      u2[k] = thres*abs(u2[k])/u2[k];
//   }
    }
}