
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>
#include <boost/math/constants/constants.hpp>
#include <boost/random.hpp>	
#include <iomanip>
using namespace std;

typedef vector<double> vec;

double normalCDF(double x)
{
    return 0.5 * (1.0 + erf(x / sqrt(2.)));
}

// Pricing European FX Option with Black Scholes
class FX_BS {
public:

    FX_BS(double S1
        , double rf
        , double rd
        , double sigma1
        , double TTM1
        , double t1
        , double K1) {

        S = S1;
        r_f = rf;
        r_d = rd;
        sigma = sigma1;
        TTM = TTM1;
        t = t1;
        K = K1;
    }

    void FX_BX_ds(double* d1, double* d2) {

        *d1 = (log(S / K) + (r_d - r_f + sigma * sigma / 2) * (TTM - t)) / (sigma * sqrt(TTM - t));
        *d2 = (log(S / K) + (r_d - r_f - sigma * sigma / 2) * (TTM - t)) / (sigma * sqrt(TTM - t));


    }

    double Call() {
        if (S < 0.) return 0.;
        if (sigma < 0.)
        {
            if (S < K * exp(-r_d * (TTM - t))) return 0.;
            else return S - K * exp(-r_d * (TTM - t));
        }
        if (fabs(TTM - t) < 0.) // maturity or not
        {
            if (S < K) return 0.;
            else return S - K;
        }

        double d1, d2;
        FX_BX_ds(&d1, &d2);

        return normalCDF(d1) * S * exp(-r_f * (TTM - t)) - normalCDF(d2) * K * exp(-r_d * (TTM - t));

    }



    double Put() {

        if (S < 0.) return 0.;
        if (sigma < 0.)
        {
            if (S > K * exp(-r_d * (TTM - t))) return 0.;
            else return K * exp(-r_d * (TTM - t)) - S;
        }
        if (fabs(TTM - t) < 0.) // maturity or not
        {
            if (S > K) return 0.;
            else return K - S;
        }

        double d1, d2;
        FX_BX_ds(&d1, &d2);

        return normalCDF(-d2) * K * exp(-r_d * (TTM - t)) - normalCDF(-d1) * S * exp(-r_f * (TTM - t));

    }

public:
    double S;
    double r_f;
    double r_d;
    double sigma;
    double TTM;
    double t;
    double K;



};

// pricing European FX Option with Montecarlo 
double MC_FX_BS(double S0
    , double r_f
    , double r_d
    , double sigma
    , double TTM
    , double K
    , int simulations
    , unsigned seed
    , string flag_Option_type) { //N is the number of steps frm the discretization

    double price = 0;
    double dt = TTM / simulations;

    default_random_engine generator(seed);

    // Main Simulation Loop

    for (int sim = 0; sim < simulations; sim++) {
        std::normal_distribution<double> normal(0, 1);
        double W = normal(generator);

        // Inner loop with discretization
        double St = S0 * exp((r_f - pow(sigma, 2) / 2) * TTM + sigma * sqrt(TTM) * W);

        if (flag_Option_type == "Call") {
            price += max((St * exp(-r_f * TTM)) - (K * exp(-r_d * TTM)), 0.0);
        }
        else if (flag_Option_type == "Put") {
            price += max((K * exp(-r_d * TTM)) - (St * exp(-r_f * TTM)), 0.0);


        }
    }

    double MC_price = (price / (double)simulations);
    return MC_price;


}

//Pricing FX European Option with Euler Discretization method
double MC_FX_BS_Euler(double S0
    , double r_f
    , double r_d
    , double sigma
    , double TTM
    , double K
    , int simulations
    , unsigned seed
    , int N
    , string flag_Option_type) { //N is the number of steps frm the discretization

    vec St(N + 1);
    St[0] = S0;
    double price = 0;
    double dt = TTM / N;

    default_random_engine generator(seed);

    // Main Simulation Loop

    for (int sim = 0; sim < simulations; sim++) {

        // Inner loop with discretization
        for (int i = 0; i < N; i++) {
            std::normal_distribution<double> normal(0, 1);
            double W = normal(generator);
            St[i + 1] = St[i] * (1 + r_f * dt + sigma * sqrt(dt) * W);
        }


        if (flag_Option_type == "Call") {
            price += max((St[N] * exp(-r_f * TTM)) - (K * exp(-r_d * TTM)), 0.0);
        }
        else if (flag_Option_type == "Put") {
            price += max((K * exp(-r_d * TTM)) - (St[N] * exp(-r_f * TTM)), 0.0);

        }
    }

    double MC_price = (price / (double)simulations);
    return MC_price;


}


// Pricing European Option under binomial pricing

class Binomial_FX_BS {




public:

    Binomial_FX_BS(double _S
        , double _r_f
        , double _r_d
        , double _sigma
        , double _U
        , double _TTM //TTM in fraction of year i.e. 6 months is 0.5
        , double _K
        , int _N
        , char _Call_Put) {

        S = _S;
        r_f = _r_f;
        r_d = _r_d;
        sigma = _sigma;
        U = _U;
        TTM = _TTM;
        K = _K;
        N = _N;
        Call_Put = _Call_Put;


    }


    vector<vector<double>> binomial_tree() {
        vector<vector<double>> binomial_matrix(N + 1, vector<double>(N + 1));
        double D = 1 / U;

        for (int i = 0; i < N + 1; i++) {
            for (int j = 0; j <= i; j++) {


                binomial_matrix[i][j] = K * pow(U, j) * pow(D, (i - j));

                //cout << binomial_matrix[i][j] << " ";

            }
            //cout << endl;
        }
        return binomial_matrix;
    }

    double Binomial_Pricer() {
        vector<vector<double>> S = binomial_tree();
        vector<vector<double>> Option_Price(N + 1, vector<double>(N + 1));
        vector<vector<double>> q(N + 1, vector<double>(N + 1)); // equivalent martingale measure

        for (int row = 0; row < N; row++) {
            for (int col = 0; col <= row; col++) {
                q[row][col] = ((S[row][col] * exp((r_d - r_f) * (TTM / N)) - S[row + 1][col])) / ((S[row + 1][col + 1] - S[row + 1][col]));
                // cout << q[row][col] << " ";
            }
            //  cout << endl;
        }
        double P = q[0][0]; // probability of going up

        for (int i = 0; i < N + 1; i++)
        {
            switch (Call_Put)
            {

            case 'C': Option_Price[N][i] = max((S[N][i] * exp(-r_f * TTM)) - (K * exp(-r_d * TTM)), 0.0); break;
            case 'P': Option_Price[N][i] = max((K * exp(-r_d * TTM)) - (S[N][i] * exp(-r_d * TTM)), 0.0); break;

            }
        }

        // go backward to discount the option price
        for (int row = N - 1; row >= 0; row--) {
            for (int col = 0; col <= row; col++) {
                Option_Price[row][col] = exp((-r_d) * (TTM / N)) * ((1 - P) * Option_Price[row + 1][col] + P * Option_Price[row + 1][col + 1]);
                //cout << Option_Price[row][col] << " ";
            }
            //cout << endl;
        }
        double binomial_result = Option_Price[0][0];
        // cout << binomial_result;
        return binomial_result;

    }

    ~Binomial_FX_BS() { };

private:
    double S;
    double r_f;
    double r_d;
    double sigma;
    double U;
    double TTM;
    double K;
    int N;
    char Call_Put;


};


int main() {
    double r_f = 0.0145; // foreign risk - free rate
    double r_d = 0.0228;// domestic risk - free rate
    double sigma = 0.077;// p.a.volatility
    double TTM = 0.5;//time to maturity
    double St = 0.6868; // 
    double K = 0.6868; //ATM K and St are equal
    int simulations = 10000;
    unsigned seed = 1234;
    int N = 20;
    double U = exp(sigma * sqrt(TTM / N));//upside stock price movement

    cout << U << endl;
    cout << 1 / U << endl;
    cout << "Pricing an EU Option Using Monte Carlo" << endl;
    cout << "Call Price is:" << MC_FX_BS(St, r_f, r_d, sigma, TTM, K, simulations, seed, "Call") << endl;
    cout << "Put Price is:" << MC_FX_BS(St, r_f, r_d, sigma, TTM, K, simulations, seed, "Put") << endl;

    cout << "Pricing an EU Option Using Euler Discretization" << endl;
    cout << "Call Price is:" << MC_FX_BS_Euler(St, r_f, r_d, sigma, TTM, K, simulations, seed, N, "Call") << endl;
    cout << "Put Price is:" << MC_FX_BS_Euler(St, r_f, r_d, sigma, TTM, K, simulations, seed, N, "Put") << endl;

    cout << "Pricing an EU Option Using Black Scholes" << endl;
    FX_BS pricer(St, r_f, r_d, sigma, TTM, 0, K);
    cout << "Call Price is:" << pricer.Call() << endl;
    cout << "Put Price is:" << pricer.Put() << endl;

    cout << "Pricing an EU Option Using Binomial Tree" << endl;
    Binomial_FX_BS binomial_call(St, r_f, r_d, sigma, U, TTM, K, N, 'C');
    Binomial_FX_BS binomial_put(St, r_f, r_d, sigma, U, TTM, K, N, 'P');
    cout << "Call Price is:" << binomial_call.Binomial_Pricer() << endl;
    cout << "Put Price is:" << binomial_put.Binomial_Pricer() << endl;

    return 0;


}



