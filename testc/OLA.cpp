#include <stdio.h>
#include <math.h>
#include <iostream> 
using namespace std;
#define PI 3.1416
#define CHUNK 4096
#define DELAY_BUF_SIZE int(CHUNK*2.5)
double delay_buf[DELAY_BUF_SIZE]={0};

inline double* han2()
{
    double *w = new double[CHUNK];
    for(int i =0;i<CHUNK;i++)
        w[i] = pow(sin(PI*(i/float(CHUNK))),2);
    return w;
}
//cache hann window
const double *han_buf = han2();



double prev_yhan[CHUNK/2]={0};

inline vector<double> get_frame_from_delay_buf(int from,int to)
{
    assert(to>from);
    to = max(to,0);
    from = min(from,DELAY_BUF_SIZE);
    vector<double> x( begin(delay_buf)+from, begin(delay_buf)+to);
    return x;
}
inline vector<double> my_proc(const int ana_end)
{
    vector<double> x= get_frame_from_delay_buf(ana_end-ceil(s*CHUNK),ana_end);
    vector<double> x_resamp= resamp(x,s,CHUNK);
    return x_resamp;

}
// =======================
//         x-1    x
//    ...|.....|.....|   delay buf
//      after processs
//    ...|.....|.....|    
//          ^-----^
//             y1
//             ^-----^
//                y2
//          ^-----^
//            out = cat(prev_y,y[:half])
//          ^--^
//         prev_yhan = han(x)[CHUNK/2:end]
//   
//=======================
extern "C" void ola(double x[],double *out)
{
    update_delay_buf(x);

    //prepare two frames for processing
    vector<double> y1 = my_proc(2*CHUNK);
    vector<double> y2 = my_proc(DELAY_BUF_SIZE);


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

