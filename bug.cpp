#include <stdio.h>
#include <math.h>
#include <iostream> 
#include <vector> 
using namespace std;
#define buf 3
#define PI 3.1415926
inline vector<double> han2(int len)
{
    vector<double> w(3);
    for(int i =0;i<len;i++)
        w[i] = pow(sin(PI*(i/float(len))),2);
    return w;
}
extern "C" void test()
{
    vector<double> han(3);
}
