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

   int M=10000;                     //Numero totale di punti estratti nell'intervallo
   int N=100;                       //Numero di blocchi in cui questi punti sono stati raggruppati
   int L=M/N;                       //Numero di punti estratti per ogni blocco
   vector<double> ave;              //Vector in cui sono state memorizzate le N stime dell'integrale (1 per ogni blocco)
   vector<double> av2;              //Vector in cui sono stati memorizati i quadrati di ogni valore di 'ave'
   vector<double> ave_prog;         //Vector in cui sono state memorizzate le medie progressive dei valori di 'ave'
   vector<double> av2_prog;         //Vector in cui sono state memorizzate le medie progressive dei valori di 'av2'
   vector<double> err_prog;         //Vector in cui sono state memorizzate le deviazioni standard progressive dei valori di 'ave_prog'

   for (int i=0; i<N; i++) {  //Riempio l'i-esimo valore di 'ave' con la media di 'L' punti estratti uniformemente tra 0 e 1
      ave.push_back(0.);
      for (int j=0; j<L; j++) {
         ave[i]+=rnd.Rannyu();
      }
      ave[i]=ave[i]/L;
      av2.push_back(pow(ave[i], 2.));  //Riempio l'i-esimo valore di 'av2' con il quadrato di dell'i-esimo valore di 'ave'
   }

   for (int i=0; i<N; i++){
      ave_prog.push_back(accumulate(ave.begin(), ave.begin()+i+1, 0.)/(i+1)); //Riempio l'i-esimo valore di 'ave_prog' con la media dei valori di 'ave' fino a i compreso
      av2_prog.push_back(accumulate(av2.begin(), av2.begin()+i+1, 0.)/(i+1)); //Riempio l'i-esimo valore di 'av2_prog' con la media dei valori di 'av2' fino a i compreso
      if (i==0)
        err_prog.push_back(0.);
      else
        err_prog.push_back(sqrt((av2_prog[i]-pow(ave_prog[i], 2.))/i));       //Riempio l'i-esimo valore di 'err_prog' con la deviazione standard della media dei valori di 'ave_prog'
   }
   
   ofstream outave;  //Stampo su file '# di estrazioni totali fatte per la stima' 'valore stimato' 'errore stimato'
   outave.open("ave.txt");
   for (int i=0; i<N; i++) {
      outave << L*(i+1) << " " << ave_prog[i] << " " << err_prog[i] << endl;
   }

   outave.close();

   //Seguo lo stesso procedimento di prima, sono identici i vector e il loro contenuto

   for (int i=0; i<N; i++) {
      ave[i]=0.;
      for (int j=0; j<L; j++) {
         ave[i]+=pow(rnd.Rannyu()-0.5, 2);   //L'unica differenza, dovuta alla diversa funzione da integrare
      }
      ave[i]=ave[i]/L;
      av2[i]=pow(ave[i], 2.);
   }

   for (int i=0; i<N; i++){
      ave_prog[i]=accumulate(ave.begin(), ave.begin()+i+1, 0.)/(i+1);
      av2_prog[i]=accumulate(av2.begin(), av2.begin()+i+1, 0.)/(i+1);
      if (i==0)
        err_prog[i]=0.;
      else
        err_prog[i]= sqrt((av2_prog[i] - ave_prog[i]*ave_prog[i])/i);

   }

   ofstream outsigma;
   outsigma.open("sigma.txt");
   for (int i=0; i<N; i++) {
      outsigma << L*(i+1) << " " << ave_prog[i] << " " << err_prog[i] << endl;
   }

   outsigma.close();

   int Nesp = 100;                     //numero di esperimenti
   M = 100;                            //numero di sottointervalli in cui è stato diviso l'intervallo [0, 1]
   L = 10000;                          //numero di estrazioni per ogni sottointervallo considerato
   vector<double> chi2;                //vector contenente i valori di chi2 per ogni esperimento svolto

   for (int i=0; i<Nesp; i++) {
      vector<int> Nhit;                 //vector contenente il numero di volte in cui si estrae un valore all'interno del sottointervallo considerato ("colpito")
      for (int j=0; j<M; j++) {
         Nhit.push_back(0.);
         for (int k=0; k<L; k++) {
            double temp=rnd.Rannyu();
            if ((temp>double(j)/M)&&(temp<double(j+1)/M)) {
               Nhit[j]++;
            }
         }
      }
      chi2.push_back(0.);
      for (int j=0; j<M; j++) {        //Si calcola poi il 'chi2' di ogni esperimento confrontando valore atteso e osservato di 'colpi andati a segno'
         chi2[i]+=double(M)/L*pow(Nhit[j]-L/M, 2.);
      }
   }

   ofstream outchi2;
   outchi2.open("chi2.txt");
   for(int i=0; i<Nesp; i++) {   //Stampo su file '# dell'esperimento considerato' 'valore di chi2'
      outchi2 << i+1 << " " << chi2[i] << endl;
   }

   outchi2.close();

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
