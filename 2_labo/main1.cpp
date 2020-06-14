/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <numeric>
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   int M=10000;                     //totale di integrali calcolati
   int N=100;                       //numero di gruppi/blocchi
   int L=M/N;                       //numero di integrali calcolati in ogni blocco
   vector<double> Integral;         //Vector in cui sono state memorizzate le N stime dell'integrale (1 per ogni blocco)
   vector<double> Integral2;        //Vector in cui sono stati memorizati i quadrati di ogni valore di 'Integral'
   vector<double> Int_prog;         //Vector in cui sono state memorizzate le medie progressive dei valori di 'Integral'
   vector<double> Int2_prog;        //Vector in cui sono state memorizzate le medie progressive dei valori di 'Integral2'
   vector<double> err_prog;         //Vector in cui sono state memorizzate le deviazioni standard progressive dei valori di 'Int_prog'

   for (int i=0; i<N; i++) {        //Riempio l'i-esimo valore di 'Integral' con la media di 'L' punti estratti uniformemente tra 0 e 1
      Integral.push_back(0.);
      for (int j=0; j<L; j++) {
         Integral[i]+= M_PI/2.*cos(M_PI*rnd.Rannyu()/2.);
      }
      Integral[i]=Integral[i]/L;
      Integral2.push_back(pow(Integral[i], 2.));         //Riempio l'i-esimo valore di 'Integral2' con il quadrato di dell'i-esimo valore di 'Integral'
   }

   for (int i=0; i<N; i++){         
      Int_prog.push_back(accumulate(Integral.begin(), Integral.begin()+i+1, 0.)/(i+1));      //Riempio l'i-esimo valore di 'Int_prog' con la media dei valori di 'Integral' fino a i compreso
      Int2_prog.push_back(accumulate(Integral2.begin(), Integral2.begin()+i+1, 0.)/(i+1));   //Riempio l'i-esimo valore di 'Int2_prog' con la media dei valori di 'Integral2' fino a i compreso
    if (i==0)
        err_prog.push_back(0.);
    else
        err_prog.push_back(sqrt((Int2_prog[i] - Int_prog[i]*Int_prog[i])/i));                //Riempio l'i-esimo valore di 'err_prog' con la deviazione standard della media dei valori di 'Int_prog'
   }
   
   ofstream outIntegral;      //Stampo su file '# di estrazioni totali fatte per la stima' 'valore stimato' 'errore stimato'
   outIntegral.open("Integral.txt");
   for (int i=0; i<N; i++) {
      outIntegral << L*(i+1) << " " << Int_prog[i] << " " << err_prog[i] << endl;
   }

   outIntegral.close();

   //per fare l'integrale con importance sampling scelgo come distribuzione di probabilità
   //la retta che passa per (1, 0), la normalizzo e calcolo l'inversa della cumulativa
   //la funzione 'myProb' di random gen resistuisce la variabile richiesta

   for (int i=0; i<N; i++) {
      Integral[i]=0.;
      for (int j=0; j<L; j++) {
         double temp=rnd.myProb();                                //L'unica differenza rispetto al caso precedente è in questa riga ed è
         Integral[i]+= M_PI/2.*cos(M_PI*temp/2.)/(-2.*temp+2.);   //l'estrazione di punti tra 0 e 1 con diversa distribuzione di probabilità
      }                                                           //e la diversa funzione valutata nel punto estratto
      Integral[i]=Integral[i]/L;
      Integral2[i]=pow(Integral[i], 2.);
   }

   for (int i=0; i<N; i++){
      Int_prog[i]=accumulate(Integral.begin(), Integral.begin()+i+1, 0.)/(i+1);
      Int2_prog[i]=accumulate(Integral2.begin(), Integral2.begin()+i+1, 0.)/(i+1);
    if (i==0)
        err_prog[i]=0.;
    else
        err_prog[i]= sqrt((Int2_prog[i] - Int_prog[i]*Int_prog[i])/i);

   }
   
   ofstream outIntIS;
   outIntIS.open("IntIS.txt");
   for (int i=0; i<N; i++) {
      outIntIS << L*(i+1) << " " << Int_prog[i] << " " << err_prog[i] << endl;
   }

   outIntIS.close();

   rnd.SaveSeed();
   return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
