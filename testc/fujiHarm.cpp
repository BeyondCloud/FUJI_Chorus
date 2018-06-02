//  g++ -Wall -std=c++11 fujiHarm.cpp -shared -m64 -o fujiHarm.dll

#include <iostream> 
#include <math.h>
#include <algorithm>  // std::min_element, std::max_element
#include <complex>
#include <memory>
#include <vector>

using namespace std; 
#define PI 3.1415926
#define CHUNK 4096
inline vector<double> han2(int len);
const vector<double> han_buf = han2(CHUNK);
//global variable for OLA
double frame1[CHUNK]={0};
double prev_x[CHUNK/2]={0};
double prev_yhan[CHUNK/2]={0};
//end OLA variable

inline vector<double> han2(int len)
{
    vector<double> w(len);
    for(int i =0;i<len;i++)
        w[i] = pow(sin(PI*(i/float(len))),2);
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
         //   cout<<j<<" "<<xlen-i+j-1<<",";
        }
       // cout<<endl;
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
    const vector<double> w = han2(winLen);
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
    for(int i =0;i<synPosLen;i++)
        cout<<" "<<synPos[i];
    cout<<endl;
    for(int i =0;i<synPosLen;i++)
        cout<<" "<<anaPos[i];
    cout<<endl;
    for(int i =1;i<synPosLen;i++)
        anaHop[i] = anaPos[i]-anaPos[i-1];

    const double minFac = synHop/(*max_element(begin(anaHop), end(anaHop)));
    vector<double> x_pad(winLen+ceil(winLen*minFac)+tol*2+xlen,0.0);
    vector<double> yC(2*winLen+yLen,0.0);
    int del = 0;
    for(int i = winLenHalf+tol;i<winLenHalf+tol+xlen;i++)
        x_pad[i-1] = x[i-winLenHalf-tol];

    //test x_pad
    // for(int i = winLenHalf+tol;i<winLenHalf+tol+10;i++)
    //     cout<<i<<" "<<x_pad[i]<<endl; 

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
        // for(int j = 0;j<300;j++)
        //     cout<<maxit[j]<<" ";

        // cout<<xNextAnaWinRan.size()<<""<<natProg.size()<<endl;
        cout<<max_idx<<endl;
        del = tol - max_idx + 1;
    }

    for(int i = synPos[synPosLen-1];i<synPos[synPosLen-1]+winLen;i++)
        yC[i] += x_pad[anaPos[synPosLen-1]+del+i-synPos[synPosLen-1]]*w[i-synPos[synPosLen-1]];
    cout<<synPos[synPosLen-1]<<" "<<anaPos[synPosLen-2];
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

inline vector<double> fujiHarm(double *x,int xlen)
{
    const double shift_note = 4.0;
    const double s = pow(2.0,shift_note/12.0);
    
    vector<double> x_vec(x, x + xlen);
    //add 0.1 for better result 
    vector<double> x_stretch= wsola(x_vec,s+0.1);

    vector<double> x_resamp= resamp(x_stretch,s,xlen);
    return x_resamp;

    //vector<double> x_pres = formantPres(x_vec,x_resamp);
    
    // for(int i = 0;i<int(x_stretch.size());i++)
    //     x[i] = x_stretch[i];

    // for(int i = 0;i<int(x_pres.size());i++)
    //     x[i] = x_pres[i];
}
extern "C" void ola(double *x,double *out)
{
    //prepare two frames for processing
    for(int i =0;i<CHUNK/2;i++)
    {
        frame1[i] = prev_x[i];
        prev_x[i] = x[i+CHUNK/2];
    }
    for(int i =CHUNK/2;i<CHUNK;i++)
        frame1[i] = x[i-CHUNK/2];

    //process two frames: take frame1 ,x as input 
    // for(int i =0;i<CHUNK;i++)
    // {
    //     y1[i] = frame1*2.0;
    //     y2[i] = x[i]*2.0; 
    // }
    vector<double> y1 = fujiHarm(frame1,CHUNK);
    vector<double> y2 = fujiHarm(x,CHUNK);
    
    //END frame1,x process


    for(int i =0;i<CHUNK;i++) 
    {
        y1[i] *= han_buf[i];
        y2[i] *= han_buf[i];
    }
    //OLA
    for(int i =0;i<CHUNK/2;i++)
        out[i] = y1[i]+prev_yhan[i];
    for(int i =CHUNK/2;i<CHUNK;i++)
        out[i] = y1[i]+y2[i-CHUNK/2];

    for(int i =0;i<CHUNK/2;i++) 
        prev_yhan[i] = y2[i+CHUNK/2];

}
