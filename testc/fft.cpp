#include <iostream> 
#include <math.h> 
#include <complex.h>
#include <vector>

using namespace std; 
#define PI 3.1415926
#define cpx complex<double>
#define vcpx vector<cpx>
#define vvcpx vector<vcpx>

inline vcpx fft( vcpx x_cpx)
{

    const int len =  x_cpx.size();
    vcpx X(len);

    double x[len*2];
    for(int i=0;i<len;i++)
    {
        x[i*2] =real(x_cpx[i]);
        x[i*2+1] =imag(x_cpx[i]);
    }

    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;
    // reverse-binary reindexing
    n = len<<1;
    j=1;
    for (i=1; i<n; i+=2) {
        if (j>i) {
            swap(x[j-1], x[i-1]);
            swap(x[j], x[i]);
        }
        m = len;
        while (m>=2 && j>m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    };
        // here begins the Danielson-Lanczos section
    mmax=2;
    while (n>mmax) {
        istep = mmax<<1;
        theta = -(2*PI/mmax);
        wtemp = sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for (m=1; m < mmax; m += 2) {
            for (i=m; i <= n; i += istep) {
                j=i+mmax;
                tempr = wr*x[j-1] - wi*x[j];
                tempi = wr * x[j] + wi*x[j-1];
 
                x[j-1] = x[i-1] - tempr;
                x[j] = x[i] - tempi;
                x[i-1] += tempr;
                x[i] += tempi;
            }
            wtemp=wr;
            wr += wr*wpr - wi*wpi;
            wi += wi*wpr + wtemp*wpi;
        }
        mmax=istep;
    }
    for(int i=0;i<len;i++)
        X[i] = {x[i*2],x[i*2+1]};
        
    return X;
}
inline vector<double> ifft(vcpx x_cpx)
{
    const int xlen = x_cpx.size();
    //inner conjugate
    vector<double> x(xlen);
    for(int i=1;i<xlen;i++)
        x_cpx[i] =conj(x_cpx[i]);
    vcpx X = fft(x_cpx);

    for(int i=0;i<xlen;i++)
        x[i] = real(X[i])/xlen;
    return x;
}

inline vcpx d2cpx(vector<double> x)
{
    const int len  = x.size();
    vcpx X(len);
    for(int i=0;i<len;i++)
        X[i] ={x[i],0};
    return X;
}
// inline vvcpx stft(double x[])
// {   
//     const int anaHop = 64;
//     const int winLen = 1024;
//     const int winLenHalf = round(winLen/2);
//     // const int signalLength = x.size();
//     // const int numOfFrames = floor((signalLength - winLen)/anaHop + 1);
//     // double winPos[numOfFrames];
    
//     vvcpx A[2][2];
//     A[0][0] = vcpx(3.0,3.0);
//     return A;
//     // for(int i = 0;i<numOfFrames;i++)
//     //     winPos[i] = i*anaHop;

// }

inline vvcpx foo(double x[])
{
    vvcpx vec(2, vcpx(2));
    for(int i =0;i<2;i++)
    {
     for(int j =0;j<2;j++)
        {
            vec[i][j] = cpx(double(i), double(j));
        }       
    }
    return vec;
}
int main() { 

    double mat[4] = {1,2,3,4};

    vvcpx vec(2, vcpx(2));
    vec[0][0] =  cpx(1.2, 3.4);
    cout<<vec[0][0]<<endl;
    vvcpx v= foo(mat);
    for(int i =0;i<2;i++)
    {
     for(int j =0;j<2;j++)
        {
            cout<<" "<<v[i][j];
        }       
    }
    return 0; 
}
