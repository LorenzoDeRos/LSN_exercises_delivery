#ifndef __finance_h__
#define __finance_h__

#include <cmath>
#include <cstdlib>
#include <algorithm>

using namespace std;

double Avanza(double S0, double mu, double sig, double dt, double Gauss);   //simula l'avanzamento del prezzo dell'asset
double Call(double r, double T, double S, double K);                        //stima il prezzo di un'opzione Call
double Put(double r, double T, double S, double K);                         //stima il prezzo di un'opzione Put

#endif //__finance_h__