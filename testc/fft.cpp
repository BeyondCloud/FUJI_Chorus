#include <iostream> 
#include <math.h> 
#include <complex>
using namespace std; 
#define PI 3.1415926
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
inline double* ifft(complex<double> x_cpx[],int len)
{

    //inner conjugate
    double *x = new double[len];
    for(int i=1;i<len;i++)
        x_cpx[i] =conj(x_cpx[i]);
    complex<double> *X = fft(x_cpx,len);

    for(int i=0;i<len;i++)
        x[i] = real(X[i])/len;
    return x;
}
inline complex<double>* d2cpx(double x[],int len)
{
    complex<double> *X = new complex<double>[len];
    for(int i=0;i<len;i++)
        X[i] ={x[i],0};
    return X;
}
int main() { 
    int len = 4;

    double x[len] = {1,2,3,4};
    // complex<double> x[len]={{1, 0},{2, 0},{3, 0},{4, 0}};
    // double x_cpx[len*2] = {0};
    // for(int i=0;i<len;i++)
    //     x_cpx[i*2] = x[i];
    complex<double> *x_cpx = d2cpx(x,len);


    complex<double> *X = fft(x_cpx,len);

    // for(int i=0;i<4;i++)
    //     cout<<real(X[i])<<endl;
    cout<<X[0]/3.0;
    double *out = ifft(X,len);
    for(int i=0;i<len;i++)
        cout<<out[i]<<endl;

    return 0; 
}
