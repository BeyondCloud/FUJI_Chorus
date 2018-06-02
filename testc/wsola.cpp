#include <iostream> 
#include <math.h> 
#include <complex>
#include <vector>
#include <algorithm>  // std::min_element, std::max_element
using namespace std; 
#define PI 3.1415926
inline vector<double> han2(int len)
{
    vector<double> w(len);
    for(int i =0;i<len;i++)
        w[i] = pow(sin(PI*(i/float(len))),2);
    return w;
}

inline vector<double> resamp(const vector<double> x,double s,int out_len)
{
    vector<double> x_resamp(out_len);
    for(double i =0;i<out_len;i++)
    {
        double interpx = i*s;
        double xl = floor(interpx);
        double dely = (xl<out_len-1)?(x[xl+1]-x[xl]):0;
        x_resamp[i] = x[xl]+(interpx-xl)*dely;
    }
    return x_resamp;

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
    // for(int i =0;i<synPosLen;i++)
    //     cout<<" "<<synPos[i];
    // cout<<endl;
    // for(int i =0;i<synPosLen;i++)
    //     cout<<" "<<anaPos[i];
    // cout<<endl;
    for(int i =1;i<synPosLen;i++)
        anaHop[i] = anaPos[i]-anaPos[i-1];

    const double minFac = synHop/(*max_element(begin(anaHop), end(anaHop)));
    vector<double> x_pad(winLen+tol*2+xlen,0.0);
    vector<double> yC(2*winLen+yLen,0.0);
    int del = 0;
    for(int i = winLenHalf+tol;i<winLenHalf+tol+xlen;i++)
        x_pad[i-1] = x[i-winLenHalf-tol];

    //test x_pad
    // for(int i = winLenHalf+tol;i<winLenHalf+tol+10;i++)
    //     cout<<i<<" "<<x_pad[i]<<endl; 

    //start wsola
    for(int i= 0;i<synPosLen;i++)
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
           natProg[j] = x_pad[ana_from]+synHop;
        int next_ana_from = anaPos[i+1]-tol;
        int next_ana_to = anaPos[i+1]+winLen-tol;
        vector<double> xNextAnaWinRan(winLen);
        for(int j = next_ana_from;j<next_ana_to;j++)
            xNextAnaWinRan[j-next_ana_from] = x[j];
        vector<double> cc = crossCorr(xNextAnaWinRan,natProg,winLen);
        //argmax
        vector<double>::iterator maxit = max_element(cc.begin(), cc.end());
        int max_idx =distance(cc.begin(), maxit);
        del = tol - max_idx + 1;
    }
    for(int i = synPos[synPosLen-1];i<synPos[synPosLen-1]+winLen;i++)
        yC[i] += x_pad[anaPos[synPosLen-2]+del+i-synPos[synPosLen-1]]*w[i-synPos[synPosLen-1]];
    cout<<synPos[synPosLen-1]<<" "<<anaPos[synPosLen-2];
    for(int i = 0;i<y.size();i++)
        y[i] = yC[i+winLenHalf+1];
    return y;
}
int main() { 

    vector<double> x{1,2,3,4,5,6};
    vector<double> y{1,2,3,4};

    int winLen = 2;
    //test window 
    // const vector<double> w = han2(1024);
    // for(int i=0;i<10;i++)
    //      cout<<w[i]<<" ";


    //test crossCorr
    // vector<double> ans = crossCorr(x,y,winLen); 
    // for(int i=0;i<int(ans.size());i++)
    //     cout<<ans[i]<<" ";


    vector<double> v(4096);
    const double s = 1.5;
    // t = [1/4096:1/4096:1];
    // x = sin(t*2*pi*40)';
    for(double i=1;i<4097;i++)
    {
        v[i] = sin(2.0*PI*i*40.0/4096.0);
    }
    vector<double> ans_w = wsola(v,s); 
    cout<<endl;
    for(int i=ans_w.size()-700;i<ans_w.size();i++)
        cout<<ans_w[i]<<" ";
    cout<<endl<<ans_w.size();
    return 0; 
}
