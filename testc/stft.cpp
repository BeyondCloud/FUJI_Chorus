#include <iostream> 
#include <math.h> 
#include <complex.h>
#include <vector>

using namespace std; 
#define PI 3.1415926
#define cpx complex<double>
#define vcpx vector<cpx>
#define vvcpx vector<vcpx>
#define fs 22050
inline vector<double> han(int len,int order);
const vector<double> han1024 = han(1024,1);
inline vector<double> han(int len,int order)
{
    vector<double> w(len);
    for(int i =0;i<len;i++)
        w[i] = pow(sin(PI*(i/float(len))),order);
    return w;
}
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
inline vvcpx stft(vcpx x)
{   
    const int anaHop = 64;
    const int winLen = 1024;
    const int xlen = x.size();
    const int winLenHalf = round(winLen/2);
    const int fsAudio = fs;
    const int signalLength = x.size();
    const int numOfFrames = floor((anaHop+winLenHalf+xlen)/anaHop + 1);


    double xPadded[winLen+anaHop+winLenHalf+xlen];
    cout <<"pad:"<<winLen+anaHop+winLenHalf+xlen<<endl;
    cout <<"numOfFrames:"<<numOfFrames<<endl;
    int winPos[numOfFrames];
    for(int i =0;i<numOfFrames;i++)
        winPos[i] = i*anaHop;

    // cout <<"winPos:";
    // for(int i =0;i<numOfFrames;i++)
    //     cout<<winPos[i]<<" ";
    cout<<endl;

    for(int i =0;i<xlen;i++)
        xPadded[i+winLenHalf] = real(x[i]);
    cout<<real(x[0]);
    vvcpx spec(winLenHalf+1, vcpx(numOfFrames));
    vcpx xi(winLen);
    for(int i =0;i<numOfFrames;i++)
    {

        for(int j =0;j<winLen;j++)
            xi[j]= cpx(xPadded[winPos[i]+j]* han1024[j],0);
        vcpx Xi = fft(xi);
        // if(i == 30)
        // {
        //     cout<<xPadded[winPos[i]];
        //     for(int i =512;i<517;i++)
        //         cout<<Xi[i]<<" ";
        //     cout<<endl;
        // }
        for(int j =0;j<winLenHalf+1;j++)
            spec[j][i] = Xi[i];
    }
    return spec;
    // for(int i = 0;i<numOfFrames;i++)
    //     winPos[i] = i*anaHop;

}

// inline vvcpx foo(v,x)
// {
//     vvcpx vec(3, vcpx(2));
//     for(int i =0;i<2;i++)
//     {
//      for(int j =0;j<2;j++)
//         {
//             vec[i][j] = cpx(double(i), double(j));
//         }       
//     }
//     cout<<"s:"<<vec.size()<<" "<<vec[0].size();

//     return vec;
// }
int main() { 

    double mat[4] = {1,2,3,4};
    vcpx v(4096);
    for(double i=1;i<4097;i++)
        v[i-1] = cpx(sin(2.0*PI*i*40.0/4096.0),0);

    vvcpx spec= stft(v);

    // vvcpx vec(2, vcpx(2));
    // vec[0][0] =  cpx(1.2, 3.4);
    // cout<<vec[0][0]<<endl;
    
    // for(int i =0;i<2;i++)
    // {
    //  for(int j =0;j<2;j++)
    //     {
    //         cout<<" "<<v[i][j];
    //     }       
    // }
    return 0; 
}
