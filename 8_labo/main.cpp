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
#include "posizione.h"

using namespace std;

double L=5.3;                  //lato del cubo in cui viene estratta uniformemente la posizione
double mu=0.79;
double sigma=0.62;
double norm=1/2/sqrt(M_PI)/sigma/(1+exp(-pow(mu/sigma, 2.)));

double wave_f(Posizione &P) {
   double x=P.getX();
   return exp(-pow( (x-mu)/sigma, 2.)/2.) + exp(-pow( (x+mu)/sigma, 2.)/2.);
}

double wave_f_second(Posizione &P) {
   double x=P.getX();
   return pow( (x-mu)/sigma/sigma, 2.)*exp(-pow( (x-mu)/sigma, 2.)/2.) + pow( (x+mu)/sigma/sigma, 2.)*exp(-pow( (x+mu)/sigma, 2.)/2.)-wave_f(P)/sigma/sigma;
}

double Potential(Posizione &P) {
   double x=P.getX();
   return pow(x, 4.)-5.*x*x/2.;
}

double Energy(Posizione &P) {
   return -wave_f_second(P)/2./wave_f(P) + Potential(P);
}

double Prob_density(Posizione &P) {
   return pow(wave_f(P), 2.);
}

double Accept(Posizione &P, Posizione &P_next) {
   return min(1., Prob_density(P_next)/Prob_density(P));
}
 
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

   
   double X_0=1.;
   Posizione P(X_0, 0., 0.);
   int M=pow(10, 6);             //total number of evaluations
   int Nblocks=100;              //number of groups
   int Nsteps=M/Nblocks;         //number of evaluations in every group (steps)
   
   double X;
   double E, E_prog, E2_prog, err_E_prog;
   double acc, acc_final;
   double alpha;
   vector<double> ave_E;         //si inseriscono le stime del valor medio della prima distribuzione di probabilità
                                 //una stima per ogni blocco
   vector<double> av2_E;
   vector<double> ave_acc;

  

   cout << "PRIMA FUNZIONE D'ONDA" << endl << endl;
   cout << "AVANZAMENTO UNIFORME IN CUBO" << endl;
   cout << "le prime 10 misure sono" << endl;

   ofstream measures_E;
   measures_E.open("prog_E.txt");
   ofstream function_plot;
   function_plot.open("probability.my");

   for(int iblock=0; iblock<Nblocks; iblock++) {
      E=0.;
      acc=0.;
      for(int j=0; j<Nsteps; j++) {
         X=rnd.Rannyu(P.getX()-L/2, P.getX()+L/2);
         Posizione P_next(X, 0., 0.);
         alpha=Accept(P, P_next);
         acc+=alpha;
         if(rnd.Rannyu()<alpha) {
            P=P_next;
         }
         function_plot << P.getX() << endl;
         E+=Energy(P);
      }
      E=E/Nsteps;
      acc=acc/Nsteps;
      ave_E.push_back(E);
      av2_E.push_back(E*E);
      ave_acc.push_back(acc);
   }

   acc_final=accumulate(ave_acc.begin(), ave_acc.end(), 0.)/Nblocks;

   for(int iblock=0; iblock<Nblocks; iblock++) {
      E_prog=accumulate(ave_E.begin(), ave_E.begin()+iblock+1, 0.)/(iblock+1);
      E2_prog=accumulate(av2_E.begin(), av2_E.begin()+iblock+1, 0.)/(iblock+1);
      if(iblock==0) {
         err_E_prog=0.;
      }
      else {
         err_E_prog=sqrt((E2_prog-pow(E_prog, 2.))/iblock);
      }
      measures_E << (iblock+1)*Nsteps << " " << E_prog << " " << err_E_prog << endl;
   }

   measures_E.close();

   cout << "Energia media finale" << endl;
   cout << E_prog << " pm " << err_E_prog << endl;
   cout << "valor medio finale della probabilità di accettazione" << endl;
   cout << acc_final << endl << endl;

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
