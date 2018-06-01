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

//fix y move x
inline vector<double> crossCorr(vector<double> x,vector<double> y,int winLen)
{
    const int xlen = x.size();
    const int ylen = y.size();
    const int olen = xlen+ylen-winLen*2+1;
    vector<double> out(olen,0.0);

    for(int i=winLen-1;i<olen+winLen-1;i++)
    {
        for(int j=max(0,i-xlen+1);j<=min(i,ylen-1);j++)
        {
            out[i-winLen+1]+=(y[j]*x[xlen-i+j-1]);
            cout<<j<<" "<<xlen-i+j-1<<",";
        }
        cout<<endl;
    }
    return out;
}

int main() { 

    vector<double> x{1,2,3,4,5,6};
    vector<double> y{1,2,3,4};
    int winLen = 2;

    vector<double> ans = crossCorr(x,y,winLen); 

    for(int i=0;i<int(ans.size());i++)
        cout<<ans[i]<<" ";
    return 0; 
}
// inline vector<double> crossCorr(vector<double> x,vector<double> y,int winLen)
// {
//     const int xlen = x.size();
//     const int ylen = y.size();
//     const int olen = xlen+ylen-winLen*2+1;
//     vector<double> out(olen,0.0);

//     for(int i=0;i<olen;i++)
//     {
//         for(int j=max(0,i-xlen+1);j<=min(i,ylen-1);j++)
//         {
//             out[i]+=(y[j]*x[xlen-i+j-1]);
//             cout<<j<<" "<<xlen-i+j-1<<",";
//         }
//         cout<<endl;
//     }
//     return out;
// }

