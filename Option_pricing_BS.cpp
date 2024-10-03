#include <iostream>
#include <cmath>
using namespace std;

double EuroCall(double S, double K, double r, double q, double T, double sig);

double EuroPut(double S, double K, double r, double q, double T, double sig);

double NP(double x);
double N(double x);


int main()
{
 double S = 50;
 double K = 50;
 double r = 0.04;
 double q = 0.01;
 double T = 0.50;
 double sig = 0.40; 
 double callPrice, putPrice;

 callPrice = EuroCall(S,K,r,q,T,sig);
 putPrice =EuroPut(S,K,r,q,T,sig);

 cout << "Call price: " << callPrice << endl;
 cout << "Put price: " << putPrice << endl;

 cout << "Put price from Put-Call Parity: " << callPrice + K * exp(-r * T) - S*exp(-q*T) << endl;

cout << endl;

system("PAUSE");
return 0;
}



double EuroCall(
    double S,
    double K,
    double r,
    double q,
    double T,
    double sig
)

{
    double d1, d2;

    d1 = (log(S/K) +(r-q+(sig*sig)*0.5)*T) / (sig*sqrt(T));
    d2 = d1 - sig * sqrt(T);


return S * exp(-q*T) * N(d1) - K*exp(-r*T) * N(d2);
}

double EuroPut(
    double S,
    double K,
    double r,
    double q,
    double T,
    double sig
)

{
    double d1, d2;

    d1 = (log(S/K) +(r-q+(sig*sig)*0.5)*T) / (sig*sqrt(T));
    d2 = d1 - sig * sqrt(T);

return K*exp(-r*T) * N(d2) - S * exp(-q*T) * N(d1);
}


double NP(double x)
{
    return(1.0/sqrt(2.0 * 3.1415) * exp(-x * x * 0.5));
}

double N(double x)
{
    double a1 = 0.319;
    double a2 = -0.356;
    double a3 = 1.781;
    double a4 = - 1.8212;
    double a5 = 1.3302;
    double k;

    k = 1/(1+0.2316 * x);

    if ( x>= 0.0)
        {
        return (1-NP(x)*((a1*k) +(a2*k*k) +(a3*k*k*k) +(a4*k*k*k*k) + (a5*k*k*k*k*k)));
        }
    else
        {
            return (1-N(-x));
        }
}



