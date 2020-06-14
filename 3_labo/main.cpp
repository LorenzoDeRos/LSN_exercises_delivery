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
#include "finance.h"
 
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

   double S0=100.;               // initial asset price
   double T=1.;                  // delivery time
   double K=100.;                // strike price
   double r=0.1;                 // risk-free interest rate
   double sig=0.25;              // volatility

   int M=10000;                  // numero di simulazioni
   int G=100;                    // numero di gruppi di simulazioni
   int N=M/G;                    // numero di simulazioni per ogni gruppo
   vector<double> C_ave;         // qui salvo le medie di Call e Put di ogni gruppo
   vector<double> C_av2;
   vector<double> P_ave;
   vector<double> P_av2;
   vector<double> C_prog;        // qui salvo le medie progressive delle medie di Call e Put
   vector<double> C2_prog;
   vector<double> P_prog;
   vector<double> P2_prog;
   vector<double> err_C;
   vector<double> err_P;
   double S;                     // variabili temporanee
   double C=0.;
   double P=0.;

   //PRIMA PARTE
   //L'evoluzione del prezzo dell'asset viene simulata in un unico avanzamento di un tempo T
   
   for (int i=0; i<G; i++) {                                //itero sul numero di gruppi di simulazioni
      for (int j=0; j<N; j++) {                             //itero sul numero di simulazioni in ogni gruppo
         S=Avanza(S0, r, sig, T, rnd.Gauss(0., 1.));        //simulo l'evoluzione del prezzo dell'asset
         C+=Call(r, T, S, K);                               //accumulo la stima per la j-esima simulazione del prezzo di un'opzione Call
         P+=Put(r, T, S, K);                                //accumulo la stima per la j-esima simulazione del prezzo di un'opzione Put
      }
      C=C/N;                                                //faccio la media delle stime per l'i-esimo gruppo del prezzo di un'opzione Call
      P=P/N;                                                //faccio la media delle stime per l'i-esimo gruppo del prezzo di un'opzione Put
      C_ave.push_back(C);                                   //memorizzo la media appena trovata per il prezzo dell'opzione Call
      C_av2.push_back(pow(C, 2.));                          //e memorizzo anche il suo quadrato
      P_ave.push_back(P);                                   //memorizzo la media appena trovata per il prezzo dell'opzione Put
      P_av2.push_back(pow(P, 2.));                          //e memorizzo anche il suo quadrato
   }
   
   for (int i=0; i<G; i++) {
         C_prog.push_back(accumulate(C_ave.begin(), C_ave.begin()+i+1, 0.)/(i+1));     //calcolo la media progressiva delle medie di ogni gruppi (Call)
         C2_prog.push_back(accumulate(C_av2.begin(), C_av2.begin()+i+1, 0.)/(i+1));    //calcolo la media progressiva delle medie di ogni gruppo (Put)
         P_prog.push_back(accumulate(P_ave.begin(), P_ave.begin()+i+1, 0.)/(i+1));     //calcolo la media progressiva dei quadrati delle medie di ogni gruppo (Call)
         P2_prog.push_back(accumulate(P_av2.begin(), P_av2.begin()+i+1, 0.)/(i+1));    //calcolo la media progressiva dei quadrati delle medie di ogni gruppo (Put)
         if (i==0) {
            err_C.push_back(0.);
            err_P.push_back(0.);
         }
         else {
            err_C.push_back(sqrt((C2_prog[i]-pow(C_prog[i], 2.))/i));                  //calcolo la deviazione standard progressiva delle medie del prezzo dell'opzione Call
            err_P.push_back(sqrt((P2_prog[i]-pow(P_prog[i], 2.))/i));                  //calcolo la deviazione standard progressiva delle medie del prezzo dell'opzione Put
         }
   }

   ofstream out1;
   out1.open("Call_1_passo.txt");
   ofstream out2;
   out2.open("Put_1_passo.txt");

   for(int i=0; i<G; i++) {      //stampo 'numero di simulazioni fatte' 'stima del prezzo dell'opzione' 'deviazione standard del prezzo dell'opzione'
      out1 << N*(i+1) << " " << C_prog[i] << " " << err_C[i] << endl;
      out2 << N*(i+1) << " " << P_prog[i] << " " << err_P[i] << endl;
   }

   out1.close();
   out2.close();


   //SECONDA PARTE
   //La struttura della simulazione è identica a quella della prima parte
   //L'unica differenza sta nel fatto che l'evoluzione del prezzo dell'asset viene simulata in
   //un numero passi=100 di avanzamenti, ciascuno di un tempo st=T/passi


   int passi=100;
   double dt=T/passi;
   C=0.;
   P=0.;

   for (int i=0; i<G; i++) {
      for (int j=0; j<N; j++) {
         S=S0;
         for (int k=0; k<passi; k++) {                      //l'unica differenza è in un questo cilo di avanzamenti successivi di un tempo dt ciascuno
            S=Avanza(S, r, sig, dt, rnd.Gauss(0., 1.));     //che prima invece erano un avanzamento unico di un tempo T
         }                                                  //il resto del codice è identico
         C+=Call(r, T, S, K);
         P+=Put(r, T, S, K);
      }
      C=C/N;
      P=P/N;
      C_ave[i]=C;
      C_av2[i]=pow(C, 2.);
      P_ave[i]=P;
      P_av2[i]=pow(P, 2.);
   }
   
   for (int i=0; i<G; i++) {
         C_prog[i]=accumulate(C_ave.begin(), C_ave.begin()+i+1, 0.)/(i+1);
         C2_prog[i]=accumulate(C_av2.begin(), C_av2.begin()+i+1, 0.)/(i+1);
         P_prog[i]=accumulate(P_ave.begin(), P_ave.begin()+i+1, 0.)/(i+1);
         P2_prog[i]=accumulate(P_av2.begin(), P_av2.begin()+i+1, 0.)/(i+1);
         if (i==0) {
            err_C[i]=0.;
            err_P[i]=0.;
         }
         else {
            err_C[i]=sqrt((C2_prog[i]-pow(C_prog[i], 2.))/i);
            err_P[i]=sqrt((P2_prog[i]-pow(P_prog[i], 2.))/i);
         }
   }

   out1.open("Call_100_passi.txt");
   out2.open("Put_100_passi.txt");

   for(int i=0; i<G; i++) {
      out1 << N*(i+1) << " " << C_prog[i] << " " << err_C[i] << endl;
      out2 << N*(i+1) << " " << P_prog[i] << " " << err_P[i] << endl;
   }

   out1.close();
   out2.close();



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
