#include <iostream> 
#include <math.h> 
#include <complex>
#include <vector>

using namespace std; 
#define PI 3.1415926

inline vector<double> wsola(const vector<double> x,double s)
{
    const int ylen = ceil(s*x.size());
    vector<double> y(ylen, 0.0);
    for(int i = 0;i<ylen;i++)
        y[i] = x[floor(i/s)];

    return y;
}
inline vector<double> resamp(const vector<double> x,double s,int out_len)
{
    vector<double> x_resamp(out_len);
    for(double i =0;i<out_len;i++)
    {
        double interpx = i*s;
        double xl = floor(interpx);
        double dely = (xl<out_len-1)?(x[xl+1]-x[xl]):0;
        x_resamp[i] = x[xl]+(interpx-xl)*dely;
    }
    return x_resamp;

}
inline vector<double> pit_mod(const vector<double> x,double shift_note)
{
    const int xlen = x.size();
    const double s = pow(2.0,shift_note/12.0);
    vector<double> x_stretch= wsola(x,s);
    vector<double> x_resamp = resamp(x_stretch,s,xlen);
    return x_resamp;
    // linear interpolate

}
int main() { 

    vector<double> x{1,2,3,4,5,6,7,8};
    const double shift_note = 9.0;
    vector<double> y = pit_mod(x,shift_note);
    for(int i=0;i<int(y.size());i++)
        cout<<" "<<y[i];
    return 0; 
}
