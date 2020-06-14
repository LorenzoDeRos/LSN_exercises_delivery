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

   int N=10000;

   ofstream outdie;
   outdie.open("die.txt");
   ofstream outexp;
   outexp.open("exp.txt");
   ofstream outlor;
   outlor.open("lor.txt");

   for (int i=0; i<N; i++) {           
      outdie << rnd.Die() << endl;                 //Stampa di una variabile estratta come un dado a 6 facce
      outexp << rnd.Exp(1.) << endl;               //Stampa di una variabile estratta come un'esponenziale
      outlor << rnd.Lorentz(0., 1.) << endl;       //Stampa di una variabile estratta come una lorentziana
   }

   outdie.close();
   outexp.close();
   outlor.close();



   ofstream outdie2;
   outdie2.open("die2.txt");
   ofstream outexp2;
   outexp2.open("exp2.txt");
   ofstream outlor2;
   outlor2.open("lor2.txt");

   for (int i=0; i<N; i++) {           
      outdie2 << (rnd.Die()+rnd.Die())/2. << endl;                         //Stampa della media di 2 variabili estratte come un dado a 6 facce
      outexp2 << (rnd.Exp(1.)+rnd.Exp(1.))/2. << endl;                     //Stampa della media di 2 variabili estratte come un'esponenziale
      outlor2 << (rnd.Lorentz(0., 1.)+rnd.Lorentz(0., 1.))/2. << endl;     //Stampa della media di 2 variabili estratte come una lorentziana
   }

   outdie2.close();
   outexp2.close();
   outlor2.close();



   ofstream outdie10;
   outdie10.open("die10.txt");
   ofstream outexp10;
   outexp10.open("exp10.txt");
   ofstream outlor10;
   outlor10.open("lor10.txt");

   int M=10;

   for (int i=0; i<N; i++) {           
      double temp1=0.;
      double temp2=0.;
      double temp3=0.;
      for (int j=0; j<M; j++) {
         temp1+=rnd.Die();
         temp2+=rnd.Exp(1.);
         temp3+=rnd.Lorentz(0., 1.);
      }
      outdie10 << temp1/M << endl;        //Stampa della media di 10 variabili estratte come un dado a 6 facce
      outexp10 << temp2/M << endl;        //Stampa della media di 10 variabili estratte come un'esponenziale
      outlor10 << temp3/M << endl;        //Stampa della media di 10 variabili estratte come una lorentziana
   }

   outdie10.close();
   outexp10.close();
   outlor10.close();



   ofstream outdie100;
   outdie100.open("die100.txt");
   ofstream outexp100;
   outexp100.open("exp100.txt");
   ofstream outlor100;
   outlor100.open("lor100.txt");

   M=100;

   for (int i=0; i<N; i++) {
      double temp1=0.;
      double temp2=0.;
      double temp3=0.;
      for (int j=0; j<M; j++) {
         temp1+=rnd.Die();
         temp2+=rnd.Exp(1.);
         temp3+=rnd.Lorentz(0., 1.);
      }
      outdie100 << temp1/M << endl;       //Stampa della media di 100 variabili estratte come un dado a 6 facce
      outexp100 << temp2/M << endl;       //Stampa della media di 100 variabili estratte come un'esponenziale
      outlor100 << temp3/M << endl;       //Stampa della media di 100 variabili estratte come una lorentziana
   }

   outdie100.close();
   outexp100.close();
   outlor100.close();



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
