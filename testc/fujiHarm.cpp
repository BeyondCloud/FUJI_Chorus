//  g++ -Wall -std=c++11 fujiHarm.cpp -shared -m64 -o fujiHarm.dll

#include <iostream> 
#include <math.h>
#include <algorithm>  // std::min_element, std::max_element
#include <complex>
#include <memory>
#include <vector>
#include <assert.h>   
using namespace std; 
#define PI 3.1415926
#define CHUNK 4096
#define DELAY_BUF_SIZE int(CHUNK*2.5)

#define cpx complex<double>
#define vcpx vector<cpx>
#define vvcpx vector<vcpx>
#define fs 22050
const double ow[CHUNK]={0};


inline vvcpx stft(vcpx x);
inline vector<double> han(int len,int order);

const vector<double> han_buf = han(CHUNK,2);
inline vector<double> istft(vvcpx spec);
//global variable for OLA
const vector<double> han1024 = han(1024,1);
double prev_yhan[CHUNK/2]={0};
double delay_buf[DELAY_BUF_SIZE]={0};

//end OLA variable
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
        }
    }
    return out;
}

inline vector<double> wsola(const vector<double> x,double s)
{
    const int xlen = x.size();
    const double synHop = 512.0;
    const int winLen = 1024;
    const int winLenHalf = round(winLen/2);
    const int tol = 512;
    const vector<double> w = han(winLen,2);
    const double yLen = ceil(s*xlen);
    const int synPosLen =(yLen+winLenHalf)/synHop;
    vector<double> y(yLen);

    double synPos[synPosLen];
    vector<double> anaPos(synPosLen);
    vector<double> anaHop(synPosLen,0.0);
    double anaStep = (xlen)*(yLen/(yLen+winLenHalf))/(synPosLen-1.0);
    for(int i =0;i<synPosLen;i++)
    {
        synPos[i] = i*synHop;
        anaPos[i] = winLenHalf + round(anaStep*i);
    }

    for(int i =1;i<synPosLen;i++)
        anaHop[i] = anaPos[i]-anaPos[i-1];

    const double minFac = synHop/(*max_element(begin(anaHop), end(anaHop)));
    vector<double> x_pad(winLen+ceil(winLen*minFac)+tol*2+xlen,0.0);
    vector<double> yC(2*winLen+yLen,0.0);
    int del = 0;
    for(int i = winLenHalf+tol;i<winLenHalf+tol+xlen;i++)
        x_pad[i-1] = x[i-winLenHalf-tol];

    //start wsola
    for(int i= 0;i<synPosLen-1;i++)
    {
        int syn_from = synPos[i];
        int syn_to = synPos[i]+winLen;
        int ana_from = anaPos[i]+del;

        for(int j = syn_from;j<syn_to;j++)
        {
            yC[j+1] += x_pad[ana_from+j-syn_from]*w[j-syn_from];
        }
        vector<double> natProg(winLen);
        for(int j = 0;j<winLen;j++)
           natProg[j] = x_pad[ana_from+synHop+j];
        int next_ana_from = anaPos[i+1]-tol;
        int next_ana_to = anaPos[i+1]+winLen+tol;
        vector<double> xNextAnaWinRan(winLen+tol*2);
        for(int j = next_ana_from;j<next_ana_to;j++)
            xNextAnaWinRan[j-next_ana_from] = x_pad[j];

        vector<double> cc = crossCorr(xNextAnaWinRan,natProg,winLen);
        //argmax
        vector<double>::iterator maxit = max_element(cc.begin(), cc.end());
        int max_idx =distance(cc.begin(), maxit);

        del = tol - max_idx + 1;
    }

    for(int i = synPos[synPosLen-1];i<synPos[synPosLen-1]+winLen;i++)
        yC[i] += x_pad[anaPos[synPosLen-1]+del+i-synPos[synPosLen-1]]*w[i-synPos[synPosLen-1]];
    for(int i = 0;i<yLen;i++)
        y[i] = yC[i+winLenHalf+1];

    return y;
}

inline vector<complex<double>> d2cpx(vector<double> x)
{
    const int len  = x.size();
    
    vector<complex<double>> X(len);
    for(int i=0;i<len;i++)
        X[i] ={x[i],0};
    return X;
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
inline vector<double> resamp(const vector<double> x,double s,int out_len)
{
    vector<double> x_resamp(out_len);
    for(double i =0;i<out_len;i++)
    {
        double interpx = i*s;
        double xl = floor(interpx);
        double dely = (xl<(out_len-1))?(x[xl+1]-x[xl]):0;
        x_resamp[i] = x[xl]+(interpx-xl)*dely;
    }
    return x_resamp;

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
inline vector<double> get_frame_from_delay_buf(int from,int to)
{
    assert(to>from);
    to = max(to,0);
    from = min(from,DELAY_BUF_SIZE);
    vector<double> x( begin(delay_buf)+from, begin(delay_buf)+to);
    return x;
}
inline vector<double> fujiHarm(const int ana_end)
{
    const double shift_note = 4.0;
    const double s = pow(2.0,shift_note/12.0);
    
    vector<double> x= get_frame_from_delay_buf(ana_end-ceil(s*CHUNK),ana_end);
    vector<double> x_ori= get_frame_from_delay_buf(ana_end-CHUNK,ana_end);
    //add 0.1 for better result 
    // vector<double> x_stretch= wsola(x_vec,s+0.1);
    // vector<double> x_vec(x, x + xlen);
    vector<double> x_resamp= resamp(x,s,CHUNK);
    // return x_resamp;
    vector<double> x_pres = formantPres(x_resamp,x_ori);
    return x_pres;

    // for(int i = 0;i<int(x_stretch.size());i++)
    //     x[i] = x_stretch[i];

    // for(int i = 0;i<int(x_pres.size());i++)
    //     x[i] = x_pres[i];
}
inline void update_delay_buf(double *x)
{
    for(int i =0;i<DELAY_BUF_SIZE-CHUNK;i++)
        delay_buf[i] = delay_buf[i+CHUNK];
    
    for(int i =DELAY_BUF_SIZE-CHUNK;i<DELAY_BUF_SIZE;i++)
        delay_buf[i] = x[i-(DELAY_BUF_SIZE-CHUNK)];

}
extern "C" void ola(double *x,double *out)
{
    
    update_delay_buf(x);
    //prepare two frames for processing

    //process two frames: take frame1 ,x as input 
    vector<double> frame1 = get_frame_from_delay_buf(CHUNK,2*CHUNK);
    for(int i =0;i<CHUNK;i++)
    {
        //out[i] = frame1[i];
        out[i] = 0;
    }
    vector<double> y1 = fujiHarm(2*CHUNK);
    vector<double> y2 = fujiHarm(DELAY_BUF_SIZE);

    //END frame1,x process


    for(int i =0;i<CHUNK;i++) 
    {
        y1[i] *= han_buf[i];
        y2[i] *= han_buf[i];
    }
    //OLA
    for(int i =0;i<CHUNK/2;i++)
        out[i] += y1[i]+prev_yhan[i];
    for(int i =CHUNK/2;i<CHUNK;i++)
        out[i] += y1[i]+y2[i-CHUNK/2];

    for(int i =0;i<CHUNK/2;i++) 
        prev_yhan[i] = y2[i+CHUNK/2];

}
