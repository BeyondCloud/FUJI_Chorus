#include <iostream> 
#include <math.h> 
#include <complex.h>
#include <vector>
#include <algorithm>  // std::min_element, std::max_element
using namespace std; 
#define PI 3.1416
#define cpx complex<double>
#define vcpx vector<cpx>
#define vvcpx vector<vcpx>
#define fs 22050
#define CHUNK 4096
inline vector<double> han(int len,int order);
const vector<double> han1024 = han(1024,1);
const double ow[CHUNK]={0};
struct stftPar_t{
    int anaHop;
    int winLen;
    int xlen;
    int winLenHalf;
    int fsAudio;
    int signalLength;
    int numOfFrames;
    vector<int> winPos;
    int numOfIter;
    int origSigLen;
    int sigLen;
    vector<double> ow;
    stftPar_t()
    {
        anaHop = 64;
        winLen = 1024;
        xlen = CHUNK;
        winLenHalf = round(winLen/2);
        fsAudio = fs;
        signalLength = CHUNK;
        numOfFrames = floor((anaHop+winLenHalf+xlen)/anaHop + 1);
        winPos.resize(numOfFrames);
        for(int i =0;i<numOfFrames;i++)
            winPos[i] = i*anaHop;
        ow.resize(CHUNK);
        for(int i = 0 ;i<numOfFrames; i++)
        {
            for(int j = 0 ;j<winLen; j++)
            {
                int pos = winPos[i]+j-winLenHalf;
                if(pos>=0 && pos<CHUNK)
                  ow[pos] +=  pow(han1024[j],2);
            }  
        }

        numOfIter = 1;
        origSigLen = CHUNK;
        sigLen =  winPos[numOfFrames-1] +winLen -1;
    }
}stftPar; 

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
        // cout<<"r"<<x[i*2]<<" "<<x[i*2+1];
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
    {
        X[i] = cpx(x[i*2],x[i*2+1]);
    }
        
    return X;
}
inline vector<double> ifft(vcpx x_cpx)
{
    const int xlen = x_cpx.size();
    //inner conjugate
    vector<double> x(xlen);
    for(int i=1;i<xlen;i++)
    {
        x_cpx[i] =conj(x_cpx[i]);


    }
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


    vector<double> xPadded(stftPar.winLen+stftPar.anaHop+stftPar.winLenHalf+stftPar.xlen,0.0);


    for(int i =0;i<stftPar.xlen;i++)
        xPadded[i+stftPar.winLenHalf] = real(x[i]);


    vvcpx spec(stftPar.numOfFrames, vcpx(stftPar.winLenHalf+1));
    vcpx xi(stftPar.winLen);
    for(int i =0;i<stftPar.numOfFrames;i++)
    {

        for(int j =0;j<stftPar.winLen;j++)
            xi[j]= cpx(xPadded[stftPar.winPos[i]+j]* han1024[j],0);
        vcpx Xi = fft(xi);


        for(int j =0;j<stftPar.winLenHalf+1;j++)
        {

            spec[i][j] = Xi[j];
            // cout<<Xi;
        }

    }


    return spec;


}
inline vector<double> istft(vvcpx spec)
{   
    vector<double> out(stftPar.origSigLen,0.0);
   //restore other side of spec
    vvcpx X(stftPar.numOfFrames,vcpx(stftPar.winLen));

    for(int i=0;i<stftPar.numOfFrames;i++)
   {
        for(int j=0;j<stftPar.winLenHalf+1;j++)
            X[i][j] = spec[i][j];
        for(int j=stftPar.winLenHalf+1;j<stftPar.winLen;j++)
            X[i][j] = conj(spec[i][stftPar.winLen-j]);
   }
   for(int i=0;i<stftPar.numOfFrames;i++)
   {
        
        vector<double> xi = ifft(X[i]);
        for(int j=0;j<stftPar.winLen;j++)
        {
            int from = stftPar.winPos[i]+j;
            if(from >=0 && from <stftPar.origSigLen)
                out[from] += xi[j] * han1024[j];
        }
   }
   for(int i=0;i<stftPar.origSigLen;i++)
        out[i] /= stftPar.ow[i];
    return out;
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
            x_out[i] += x[i+y_ini-j]*y[j];
    }
    return x_out;

}
const vector<double> han24 = han(24,2);
inline vector<vector<double>> compEnv( vvcpx X)
{
    vector<vector<double>> env(stftPar.numOfFrames, vector<double>(stftPar.winLenHalf+1));
    vector<double> X_abs(stftPar.winLenHalf+1);
    for(int i =0;i<stftPar.numOfFrames;i++)
    {
        for(int j =0;j<stftPar.winLenHalf+1;j++)
            X_abs[j] = abs(X[i][j]);
       
        env[i] = conv(X_abs,han24);
        const double maxEnv = *max_element(begin(env[i]), end(env[i]));
        const double maxX = *max_element(begin(X_abs), end(X_abs));
        
        for(int j =0;j<stftPar.winLenHalf+1;j++)
        {
            env[i][j] /= (maxEnv*maxX);
            env[i][j] = max(env[i][j],0.02);
        } 
    }

    return env;
}
inline vector<double> formantPres(vector<double> pit,vector<double> ori)
{
    const int xlen = pit.size();
    vvcpx ORI = stft(d2cpx(ori));
    vvcpx PIT = stft(d2cpx(pit));


    vector<vector<double>>  envO =  compEnv(ORI);
    vector<vector<double>>  envP =  compEnv(PIT);
    
    for(int i =0;i<stftPar.numOfFrames;i++)
    {
        for(int j =0;j<stftPar.winLenHalf+1;j++)
            PIT[i][j] = (PIT[i][j]/envP[i][j])*envO[i][j];
    }

    vector<double> out  = istft(PIT);
    return out;
}
int main() { 

    vcpx v(4096);
    vector<double> vd(4096);
    
    for(double i=1.0;i<4097.0;i++)
    {
        v[i-1] = cpx(sin(2.0*PI*i*40.0/4096.0),0);
        vd[i-1] = sin(2.0*PI*i*40.0/4096.0);
    }
 
    vector<double> result  = formantPres(vd,vd);
    for(int i=0;i<CHUNK;i++)
        cout<<result[i]<<" ";
   
    return 0; 
}
