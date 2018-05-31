#include <stdio.h>
#include <math.h>
#include <iostream> 
using namespace std;
#define PI 3.1416
#define CHUNK 4096

inline double* han2()
{
    double *w = new double[CHUNK];
    for(int i =0;i<CHUNK;i++)
        w[i] = pow(sin(PI*(i/float(CHUNK))),2);
    return w;
}
//cache hann window
const double *han_buf = han2();


double frame1[CHUNK]={0};
double prev_x[CHUNK/2]={0};
double prev_han[CHUNK/2]={0};
// =======================
//  x-1    x
// .....|.....|
//   ^-----^
//    frame1
//   ^-----^
//     out
//      ^-----^
//    current frame
//   ^--^
//  prev_x
//   ^--^
//  prev_han = han(x)[CHUNK/2:end]
//   
//=======================
extern "C" void ola(double x[],double *out)
{
    //prepare two frames for processing
    for(int i =0;i<CHUNK/2;i++)
    {
        frame1[i] = prev_x[i];
        prev_x[i] = x[i+CHUNK/2];
    }
    for(int i =CHUNK/2;i<CHUNK;i++)
        frame1[i] = x[i-CHUNK/2];
    //Start frame1,x process
    
    //END frame1,x process


    for(int i =0;i<CHUNK;i++) 
        frame1[i] *= han_buf[i];
    //OLA
    for(int i =0;i<CHUNK/2;i++)
    {
        out[i] = frame1[i]+prev_han[i];

    }
    for(int i =0;i<CHUNK;i++) 
        x[i] *= han_buf[i];

    for(int i =CHUNK/2;i<CHUNK;i++)
        out[i] =frame1[i]+x[i-CHUNK/2];

    for(int i =0;i<CHUNK/2;i++) 
        prev_han[i] = x[i+CHUNK/2];

}

