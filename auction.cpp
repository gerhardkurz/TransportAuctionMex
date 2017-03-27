/*
 * MATLAB Wrapper for Transport Problem Auction Algorithm
 *
 * by Gerhard Kurz
 *
 * original auction code by Joseph D Walsh III <math@jdwalsh03.com>
 */

#include "mex.h"
#include <cmath>
#include <iostream>
#include <chrono>
#include <numeric>     //accumulate
#include "glob.hpp"    // mfloat, mfvec, mint, muint
#include "object.hpp"  // objlist, voblist
#include "gamap.hpp"   // GAmap

mfloat gEPS = std::sqrt(std::numeric_limits<mfloat>::epsilon());
mfloat gINF = std::numeric_limits<mfloat>::infinity();
muint  gVBS = 0;  // verbosity (0, 1, or 2)

/** --- GArun ------------------------------------------------------------------
 **/
void GArun (const mfvec& DWT, const mfvec& SWT, const voblist& A,
        const mfloat MX, const mfloat MN, const mfloat ST,
        objlist& T, mfvec& PR) {
    GAmap gslv (DWT, SWT, A, MX, MN, ST);
    gslv.Solve (T, PR);
    return;
}

/** --- Dual -------------------------------------------------------------------
 **/
mfloat Dual (const mfvec& DWT, const mfvec& SWT, const voblist& ARX,
        const mfvec& PR) {
    mfloat dcst = 0.0;
    mfvec  M;
    for (muint i = 0; i < DWT.size(); i++) {
        M.clear();
        for (muint j = 0; j < ARX[i].size(); j++) {
            M.push_back (ARX[i][j].c - PR[muint (ARX[i][j].j)]);
            if (gVBS > 1) {
                printf (" Exp %lu -> %lu : %f - %f = %f\n", i, j,
                        ARX[i][j].c, PR[muint (ARX[i][j].j)], M.back());
            }
        }
        dcst += DWT[i] * *std::max_element (M.begin(), M.end());
    }
    if (gVBS > 1) {
        printf("---------------\n");
    }
    for (muint it = 0; it < SWT.size(); it++) {
        dcst += SWT[it] * PR[it];
        if (gVBS > 1) {
            printf (" Price %lu : %f\n", it, PR[it]);
        }
    }
    if (gVBS > 1) {
        printf("---------------\n");
    }
    return dcst;
}

/** --- Primal -----------------------------------------------------------------
 **/
mfloat Primal (const char* aname, const voblist& ARX, const objlist& T) {
    mfloat pcst = 0.0;
    muint  tst;
    for (muint it = 0; it < T.size(); it++) {
        tst = 0;
        while (   (tst < ARX[muint (T[it].i)].size())
        && (ARX[muint (T[it].i)][tst].j != T[it].j)) {
            tst++;
        }
        if (tst < ARX[muint (T[it].i)].size()) {
            pcst += ARX[muint (T[it].i)][tst].c * T[it].c;
            if (gVBS > 1) {
                printf ("    %3ld -> %3ld : %f @ %f \n", T[it].j, T[it].i,
                        T[it].c, ARX[muint (T[it].i)][tst].c);
            }
        } else {
            printf("Overflow in %s cost calculation\n", aname);
        }
    }
    return pcst;
}

void mexFunction(int numOutputs, mxArray *outputs[],
        int numInputs, const mxArray *inputs[])
{
    try {
        /* Check for proper number of arguments */
        if (numInputs != 3) {
            mexErrMsgIdAndTxt("Auction:inputs",
                    "Three inputs are required.");
        }
        
        /*
         * if (numOutputs != 3) {
         * throw std::invalid_argument("Three outputs are required.");
         * }
         *
         * if (sampleWeights.cols() != numSamples) {
         * throw std::invalid_argument("Number of weights has to match the number of samples.");
         * }
         *
         * const unsigned int n = *mxGetPr(inputs[2]);
         *
         * if (numSamples < n) {
         * throw std::invalid_argument("Need more samples than Gaussian components.");
         * }        */
        
        if(mxGetN(inputs[0]) != 1) {
            mexErrMsgIdAndTxt("Auction:inputs",
                    "First Input must be a column  vector.");
        }
        
        if(mxGetN(inputs[1]) != 1) {
            mexErrMsgIdAndTxt("Auction:inputs",
                    "Second Input must be a column  vector.");
        }
        
        size_t nSources = mxGetM(inputs[0]);
        size_t nSinks = mxGetM(inputs[1]);
        
        if(mxGetM(inputs[2]) != nSources) {
            mexErrMsgIdAndTxt("Auction:inputs",
                    "Rows of cost matrix need to match number of sources.");
        }        
        
        if(mxGetN(inputs[2]) != nSinks) {
            mexErrMsgIdAndTxt("Auction:inputs",
                    "Columns of cost matrix need to match number of sinks.");
        }        
        
        //todo there is an assumption N>=M, otherwise soruces and sinks have to be swapped
                
        mfvec DWT; //demand
        DWT.assign(mxGetPr(inputs[0]), mxGetPr(inputs[0]) + nSources);
        mfvec SWT; //supply
        SWT.assign(mxGetPr(inputs[1]), mxGetPr(inputs[1]) + nSinks);
        
        if(!equal(std::accumulate(DWT.begin(), DWT.end(), 0), std::accumulate(SWT.begin(), SWT.end(), 0))){
            mexErrMsgIdAndTxt("Auction:inputs",
                    "Sum of demand and sum of supply nust match");
        }
        
        voblist ARX; //arcs in the graph
        ARX.resize(nSources);
        for(int i=0; i<nSources; i++){
            ARX[i].resize(nSinks);
            for(int j=0; j<nSinks; j++){
                Object o;
                o.c = *(mxGetPr(inputs[2]) + j*nSources + i);
                o.i=i;
                o.j=j;
                ARX[i][j] = o; //todo matrix may be transposed?
            }
        }
        
        mfloat eps = -1.0;
        mfloat stp = 0.25;
        mfloat min = -1.0;
        
        objlist T;     // transport plan
        mfvec   PR;    // price vector
        mfloat  pcst;  // primal cost
        mfloat  dcst;  // dual cost
        
        mfloat C = 0;  // maximum cost
        muint  ar = 0; // number of arcs in graph
        for (muint i = 0; i < ARX.size(); i++) {
            for (muint j = 0; j < ARX[i].size(); j++) {
                ar++;
                if (C < std::abs(ARX[i][j].c)) {
                    C = std::abs(ARX[i][j].c);
                }
            }
        }
        if (eps < gEPS) {
            eps = C / 5.0;
        }
        if (min < gEPS) {
            min = 1.0 / mfloat (SWT.size());
        }
        
        printf("  GRAPH: %lu sinks, %lu sources, %lu arcs\n", DWT.size(),
                SWT.size(), ar);
        printf("  EPS  : %f starting, %e minimum\n", eps, min);

        std::chrono::high_resolution_clock::time_point t1, t2;
        std::chrono::duration <mfloat> dur;        
        t1 = std::chrono::high_resolution_clock::now();
        
        GArun (DWT, SWT, ARX, eps, min, stp, T, PR);

        t2 = std::chrono::high_resolution_clock::now();
        dur = std::chrono::duration_cast <std::chrono::duration <mfloat> >
                (t2 - t1);
        if (T.empty()) {
            printf ("  General auction NOT SOLVED\n");
        } else {
            pcst = Primal ("GA", ARX, T);
            printf ("  General auction primal cost    : %25.15f\n", pcst);
            dcst = Dual (DWT, SWT, ARX, PR);
            printf ("  General auction dual cost      : %25.15f\n", dcst);
            printf ("  General auction diff           : %25.15f\n",
                    dcst-pcst);
        }
        printf ("  General auction time           : %13.3f sec\n",
                dur.count() );
        
        //print plan
        /*
        for(int i=0; i<T.size(); i++){
            printf("%f %i %i\n", T[i].c, T[i].i, T[i].j);
        }*/
        
        if(numOutputs >= 1){
            outputs[0] = mxCreateDoubleMatrix(nSources,nSinks,mxREAL);
            for(int i=0; i<T.size(); i++){
                *(mxGetPr(outputs[0]) + T[i].j*nSources +T[i].i) = T[i].c;
            }            
        }
        
    } catch (std::exception& ex) {
        //usage();
        mexErrMsgTxt(ex.what());
    }
}
