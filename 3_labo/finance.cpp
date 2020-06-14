#include "finance.h"

double Avanza(double S0, double mu, double sig, double dt, double Gauss) {
   return S0*exp((mu-pow(sig, 2.)/2.)*dt+sig*Gauss*sqrt(dt));
}

double Call(double r, double T, double S, double K) {
    return exp(-r*T)*max(0., S-K);
}
double Put(double r, double T, double S, double K) {
    return exp(-r*T)*max(0., K-S);
}