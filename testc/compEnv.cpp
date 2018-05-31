#include <iostream> 
#include <math.h>
#include <algorithm>  // std::min_element, std::max_element
#include <complex>
using namespace std; 
#define PI 3.1415926
inline double* han2(int len);
const double *han1024 = han2(1024);
const double *han8 = han2(8);

bool comp(int i, int j) { return i<j; }
inline complex<double>* d2cpx(double x[],int len)
{
    complex<double> *X = new complex<double>[len];
    for(int i=0;i<len;i++)
        X[i] ={x[i],0};
    return X;
}
inline complex<double>* fft(complex<double> x_cpx[],int len)
{


    complex<double> *X = new complex<double>[len];
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
inline double* ifft(complex<double> x_cpx[],int xlen)
{

    //inner conjugate
    double *x = new double[xlen];
    for(int i=1;i<xlen;i++)
        x_cpx[i] =conj(x_cpx[i]);
    complex<double> *X = fft(x_cpx,xlen);

    for(int i=0;i<xlen;i++)
        x[i] = real(X[i])/xlen;
    return x;
}

inline double* han2(int len)
{
    double *w = new double[len];
    for(int i =0;i<len;i++)
        w[i] = pow(sin(PI*(i/float(len))),2);
    return w;
}

inline double* conv(double *x,double *y,int xlen,int ylen)
{
    double *x_out = new double[xlen];
    const int y_ini = floor(ylen/2);
    for(int i =0;i<xlen;i++)
    {
        x_out[i] = 0;
        for(int j = min(y_ini+i,ylen-1);j>=max(0,i-xlen+1+y_ini);j--)
        {
            x_out[i] += x[i+y_ini-j]*y[j];
            // cout<<i+y_ini-j<<" "<<j<<",";
        }
        // cout<<endl;
    }
    return x_out;

}

inline double* compEnv( complex<double> X[],int Xlen,int filLen)
{
    double X_abs[Xlen];
    for(int i =0;i<Xlen;i++)
        X_abs[i] = abs(X[i]);
    double *win = han2(filLen);
    double *env = conv(X_abs,win,Xlen,filLen);
    const double *maxEnv = max_element(env,env+Xlen,comp);
    const double *maxX = max_element(X_abs,X_abs+Xlen,comp);
    
    for(int i =0;i<Xlen;i++)
    {
        env[i] /= (maxEnv[0]*maxX[0]);
        env[i] = max(env[i],0.02);
    }
    return env;
}

inline double* formantPres(double *ori,double *pit,int xlen)
{
    const int filLen = 24;
    for(int i =0;i<xlen;i++)
    {
        ori[i] *= han8[i];
        pit[i] *= han8[i];   
    }

    complex<double> *ORI = fft(d2cpx(ori,xlen),xlen);
    complex<double> *PIT = fft(d2cpx(pit,xlen),xlen);


    double *envO =  compEnv(ORI,xlen,filLen);
    double *envP =  compEnv(PIT,xlen,filLen);
    
    complex<double> fixed[xlen];
    double *out = new double[xlen];
    for(int i =0;i<xlen;i++)
        fixed[i] = (ORI[i]/envO[i])*envP[i];
    out = ifft(fixed,xlen);
    return out;
}

int main() { 
    double *win = han2(32);
    // for(int i =0;i<32;i++)
    //     cout<<win[i]<<" ";
    int x_len = 8;
    int y_len = 8;
    double x[x_len] = {1,2,3,4,5,6,7,8};
    double y[y_len] = {2,3,4,5,6,7,8,9};
    double *x_out = conv(x,y,x_len,y_len);
    double *ans = formantPres(x,y,x_len);

    for(int i =0;i<x_len;i++)
        cout<<ans[i]<<" ";

    // const double *maxX = max_element(x,x+x_len,comp);
    // cout<<maxX[0]<<" ";
    return 0; 
}

