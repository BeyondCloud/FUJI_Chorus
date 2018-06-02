//  g++ -Wall -std=c++11 fujiHarm.cpp -shared -m64 -o fujiHarm.dll

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
        yC[i] += x_pad[anaPos[synPosLen-2]+del+i-synPos[synPosLen-1]]*w[i-synPos[synPosLen-1]];
    cout<<synPos[synPosLen-1]<<" "<<anaPos[synPosLen-2];
    for(int i = 0;i<yLen;i++)
        y[i] = yC[i+winLenHalf+1];

    return y;
}
extern "C" void fujiHarm(double *x,int xlen,double *out)
{
    const double shift_note = 9.0;
    const double s = pow(2.0,shift_note/12.0);
    vector<double> x_vec(x, x + xlen);

    vector<double> x_stretch= wsola(x_vec,s);
    for(int i = 0;i<int(x_stretch.size());i++)
        out[i] = x_stretch[i];
}
