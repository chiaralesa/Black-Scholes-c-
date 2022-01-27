// Black Scholes model.cpp 

#include <iostream>
#include <cmath>
using namespace std;


double normalCDF(double x)
{
    return 0.5 * (1.0 + erf(x / sqrt(2.)));
}

double Call(double S, double K, double T, double t, double sigma, double r)
{
    if (S < 0. )return 0.;
    if (sigma < 0. )
    {
        if (S < K * exp(-r * (T - t)))return 0.;
        else return S - K * exp(-r * (T - t));
    }
    if (fabs(T - t) < 0.) // maturity or not
    {
        if (S < K )return 0.;
        else return S - K;
    }
    // calculate option price
    double d1 = (log(S / K) + (r + sigma * sigma / 2) * (T - t)) / (sigma * sqrt(T - t));
    double d2 = (log(S / K) + (r - sigma * sigma / 2) * (T - t)) / (sigma * sqrt(T - t));
    return normalCDF(d1) * S - normalCDF(d2) * K * exp(-r * (T - t));
}

    
int main(){
    cout << "Call Option Price = " << Call(200, 150, 2, 1, 0.2, 0.1) << endl;
    return 0;
}
