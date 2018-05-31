#include <iostream> 
#include <math.h>
#include <algorithm>  // std::min_element, std::max_element
#include <complex>
#include <memory>
#include <vector>

using namespace std; 
#define PI 3.1415926
#define BUF_SIZE 8
inline vector<double> han2(int len);
const vector<double> buf_han = han2(BUF_SIZE);

inline vector<complex<double>> d2cpx(vector<double> x)
{
    const int len  = x.size();
    
    vector<complex<double>> X(len);
    for(int i=0;i<len;i++)
        X[i] ={x[i],0};
    return X;
}
inline vector<complex<double>> fft( vector<complex<double>> x_cpx)
{

    const int len =  x_cpx.size();
    vector<complex<double>> X(len);

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
inline vector<double> ifft(vector<complex<double>> x_cpx)
{
    const int xlen = x_cpx.size();
    //inner conjugate
    vector<double> x(xlen);
    for(int i=1;i<xlen;i++)
        x_cpx[i] =conj(x_cpx[i]);
    vector<complex<double>> X = fft(x_cpx);

    for(int i=0;i<xlen;i++)
        x[i] = real(X[i])/xlen;
    return x;
}

inline vector<double> han2(int len)
{
    vector<double> w(len);
    for(int i =0;i<len;i++)
        w[i] = pow(sin(PI*(i/float(len))),2);
    return w;
}

inline vector<double> conv(vector<double> x,vector<double> y)
{
    const int xlen  = x.size();
    const int ylen  = y.size();
    vector<double> x_out(xlen);
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

inline vector<double> compEnv( vector<complex<double>> X,int filLen)
{
    const int Xlen = X.size();
    vector<double> X_abs(Xlen);
    for(int i =0;i<Xlen;i++)
        X_abs[i] = abs(X[i]);
    vector<double> win = han2(filLen);
    vector<double> env = conv(X_abs,win);
    const double maxEnv = *max_element(begin(env), end(env));
    const double maxX = *max_element(begin(X_abs), end(X_abs));
    
    for(int i =0;i<Xlen;i++)
    {
        env[i] /= (maxEnv*maxX);
        env[i] = max(env[i],0.02);
    }
    return env;
}

inline vector<double> formantPres(vector<double> ori,vector<double> pit)
{
    const int filLen = 24;
    const int xlen = ori.size();
    vector<complex<double>> ORI = fft(d2cpx(ori));
    vector<complex<double>> PIT = fft(d2cpx(pit));


    vector<double> envO =  compEnv(ORI,filLen);
    vector<double> envP =  compEnv(PIT,filLen);
    
    vector<complex<double>> fixed(xlen);
    for(int i =0;i<xlen;i++)
        fixed[i] = (ORI[i]/envO[i])*envP[i];
    vector<double> out  = ifft(fixed);
    return out;
}

int main() { 

    cout<<"using buffer:"<<BUF_SIZE<<endl;
    vector<double> win = han2(32);
    // for(int i =0;i<32;i++)
    //     cout<<win[i]<<" ";
    int x_len = 8;
    int y_len = 8;
    vector<double> x = {1,2,3,4,5,6,7,8};
    vector<double> y = {2,3,4,5,6,7,8,9};
    vector<double> x_out = conv(x,y);
    //windowing
    for(int i =0;i<x_len;i++)
    {
        x[i] *= buf_han[i];
        y[i] *= buf_han[i];   
    }

    vector<double> ans = formantPres(x,y);

    for(int i =0;i<x_len;i++)
        cout<<ans[i]<<" ";

    // const double *maxX = max_element(x,x+x_len,comp);
    // cout<<maxX[0]<<" ";
    return 0; 
}

