#include <iostream> 
#include <math.h> 
#include <complex.h>
#include <vector>

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
    cout <<"pad:"<<stftPar.winLen+stftPar.anaHop+stftPar.winLenHalf+stftPar.xlen<<endl;
    cout <<"numOfFrames:"<<stftPar.numOfFrames<<endl;
    // cout <<"winPos:";
    // for(int i =0;i<numOfFrames;i++)
    //     cout<<winPos[i]<<" ";
    cout<<endl;

    for(int i =0;i<stftPar.xlen;i++)
        xPadded[i+stftPar.winLenHalf] = real(x[i]);
    // cout<< xPadded[4604]<< " "<< xPadded[4605]<< " "<< xPadded[4606]<< " "<<endl;
   
    // for(int i =0;i<stftPar.xlen;i++)
    //  cout<< x[i]<< " ";
    

    vvcpx spec(stftPar.numOfFrames, vcpx(stftPar.winLen));
    vcpx xi(stftPar.winLen);
    for(int i =0;i<stftPar.numOfFrames;i++)
    {

        for(int j =0;j<stftPar.winLen;j++)
        {
            xi[j]= cpx(xPadded[stftPar.winPos[i]+j]* han1024[j],0);
            if(isnan(real(xi[j])))
            {
                cout<<stftPar.winPos[i]+j<<" ";
                cout<<xPadded[stftPar.winPos[i]+j]<<endl;
                // cout<<real(xi[j])<<" ||"<<imag(xi[j]);
            }
            // if (stftPar.winPos[i]+j >stftPar.winLen+stftPar.anaHop+stftPar.winLenHalf+stftPar.xlen)
            //     cout<<"errrr:";
        }
        vcpx Xi = fft(xi);
        // for(int j=0;j<10;j++)
        //     cout<<Xi[j]<<" ";

        // for(int j=0;j<10;j++)
        //     cout<<xPadded[stftPar.winPos[i]+j]<<" ";
        // cout<<endl;
        // if(i == 30)
        // {
        //     cout<<xPadded[winPos[i]];
        //     for(int i =512;i<517;i++)
        //         cout<<Xi[i]<<" ";
        //     cout<<endl;
        // }
        for(int j =0;j<stftPar.winLen;j++)
        {
            spec[i][j] = Xi[j];
            // cout<<Xi;
        }
    }

    // for(int i =0;i<stftPar.numOfFrames;i++)
    // {
    //     cout<<i<<":";
    //     for(int j =0;j<10;j++)
    //         cout<<real(spec[i][j])<<" ";
    //     cout<<endl;
    // }
    return spec;
    // for(int i = 0;i<numOfFrames;i++)
    //     winPos[i] = i*anaHop;

}
inline vector<double> istft(vvcpx spec)
{   
    vector<double> out(stftPar.origSigLen,0.0);

   for(int i=0;i<stftPar.numOfFrames;i++)
   {
        
        vector<double> xi = ifft(spec[i]);
       // for(int j=500;j<520;j++)
       //  {
       //      cout <<xi[j]<<" ";
       //  }
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
int main() { 

    vcpx v(4096);
    for(double i=1.0;i<4097.0;i++)
    {
        v[i-1] = cpx(sin(2.0*PI*i*40.0/4096.0),0);
        // cout<<v[i-1]<<" "<<i-1;
    }
    // for(int i=0;i<1024;i++)
    // {
    //     if(isnan(real(v[i])))
    //     {
    //         cout<<v[i];
    //         cout<<real(v[i])<<" ||"<<imag(v[i]);
    //     }

    //     // cout<<"r"<<x[i*2]<<" "<<x[i*2+1];
    // }
    // vcpx a= fft(v);
    
    // for(double i=38;i<45;i++)
    //     cout<<abs(a[i])<<" ";
    vvcpx spec= stft(v);
    for(double i=38;i<45;i++)
       cout<<abs(spec[0][i])<<" ";
    vector<double> out  =istft(spec);
    for(int i=0;i<4096;i++)
        cout<<out[i]<<" ";




    // cout<<"s size:";
    // cout<<spec.size()<<" "<<spec[0].size();

    // for(int i = CHUNK-1024;i<CHUNK; i++)
    // {
    //     cout<<stftPar.ow[i]<<" ";
    // }
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
