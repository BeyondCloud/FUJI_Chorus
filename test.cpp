#include <stdio.h>
#include <math.h>
//  g++ -Wall -O3 -std=c++11 test.cpp -shared -m64 -o chorus.dll
#define PI 3.1416
const double fs = 44100;
const double delay_sec = 0.013;
const double depth =0.003;
const double delay_samples = delay_sec * fs;
const double depth_samples  = depth * fs;
const double mod_rate = 3.2;
const double mod_arg = 2*3.1416*mod_rate/fs;
const double feedback = 0.3;
const int buf_len = delay_samples+depth_samples;
double prev_i = 0.0;
extern "C" void chorus(const double *x,double *out,double *delay_buf,int len)
{
    // printf("%d %lf",len,prev_i);

    for(int i=0;i<len;i++)
    {
        double modulated_sample = depth_samples*sin(mod_arg*i+prev_i);
        modulated_sample += (delay_samples-1);
        double interp1 =delay_buf[(int)floor(modulated_sample)];
        out[i] = interp1;
        double new_sample= x[i]+interp1*feedback;
        //roll 1 element

        for (int j = buf_len-1; j >= 1; j--)   
            delay_buf[j]=delay_buf[j-1];
        delay_buf[0] = new_sample;
    }
    prev_i = fmod((mod_arg*(len-1)+prev_i),(PI*2));
//     modulated_sample = depth_samples*math.sin(mod_arg*(i)+prev_i)
//     modulated_sample += (delay_samples-1)
// #       interpolate
//     interp1 =delay_buf[math.floor(modulated_sample)];
// #         interp2 =delay_buf[math.ceil(modulated_sample)];
// #         query_sample = math.modf(modulated_sample)[0] #get decimal part
//     delay_chnl[i] = interp1
//     new_sample= x[i]+delay_chnl[i]*feedback;
//     delay_buf[-1] = new_sample
//     delay_buf = np.roll(delay_buf,1)
}
extern "C" void test(double *x)
{
    x[0] = 1;
}