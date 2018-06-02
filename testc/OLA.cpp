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
double prev_yhan[CHUNK/2]={0};
// =======================
//  x-1    x
// .....|.....|
//   ^-----^
//    frame1=cat(prev_x,x[:half])
//      ^-----^
//      frame2=x
//   ^--^
//  prev_x
//.................after processs............
//   ^-----^
//      y1
//      ^-----^
//         y2
//   ^-----^
//     out = cat(prev_y,y[:half])
//   ^--^
//  prev_yhan = han(x)[CHUNK/2:end]
//   
//=======================
extern "C" void ola(double x[],double *out)
{
    double y1[CHUNK],y2[CHUNK];
    //prepare two frames for processing
    for(int i =0;i<CHUNK/2;i++)
    {
        frame1[i] = prev_x[i];
        prev_x[i] = x[i+CHUNK/2];
    }
    for(int i =CHUNK/2;i<CHUNK;i++)
        frame1[i] = x[i-CHUNK/2];

    //process two frames: take frame1 ,x as input 
    for(int i =0;i<CHUNK;i++)
    {
        y1[i] = frame1[i]*2.0;
        y2[i] = x[i]*2.0; 
    }
    //END frame1,x process


    for(int i =0;i<CHUNK;i++) 
    {
        y1[i] *= han_buf[i];
        y2[i] *= han_buf[i];
    }
    //OLA
    for(int i =0;i<CHUNK/2;i++)
        out[i] = y1[i]+prev_yhan[i];
    for(int i =CHUNK/2;i<CHUNK;i++)
        out[i] = y1[i]+y2[i-CHUNK/2];

    for(int i =0;i<CHUNK/2;i++) 
        prev_yhan[i] = y2[i+CHUNK/2];

}

