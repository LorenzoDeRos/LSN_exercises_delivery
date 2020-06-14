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

   double d=3.;
   double L=1.;
   int M=100;                 //Numero di gruppi/blocchi di esperimenti
   int Nesp=100;              //Numero di esperimenti svolti per ogni blocco
   int N=10000;               //Numero di lanci per ogni esperimento
   double y1;                 //Coordinata y del primo estremo dello stuzzicadente
   double x2;                 //Coordinata x del secondo estremo dello stuzzicadente
   double y2;                 //Coordinata y del secondo estremo dello stuzzicadente
   double cos;
   vector<double> ave;        //Vector in cui sono state memorizzate le M stime di pi greco (1 per ogni esperimento)
   vector<double> av2;        //Vector in cui sono stati memorizati i quadrati di ogni valore di 'ave'
   vector<double> ave_prog;   //Vector in cui sono state memorizzate le medie progressive dei valori di 'ave'
   vector<double> av2_prog;   //Vector in cui sono state memorizzate le medie progressive dei valori di 'av2'
   vector<double> err_prog;   //Vector in cui sono state memorizzate le deviazioni standard progressive dei valori di 'ave_prog'
   
   for (int i=0; i<M; i++) {
      ave.push_back(0.);
      for (int j=0; j<Nesp; j++) {
         int temp=0;                                        //itero sul numero M di esperimenti
         for (int k=0; k<N; k++) {                          //itero sul numero N di lanci
            y1=rnd.Rannyu(0, d);                            //estrazione della y del primo estremo
            do {
            x2=rnd.Rannyu(-L, L);                           //estrazione delle coordinate di un punto nel quadrato
            y2=rnd.Rannyu(y1-L, y1+L);                      //di lato 2L centrato in y
            } while (sqrt(pow(x2, 2.)+pow(y2-y1, 2.))>L);   //continuo a estrarre punti nel quadrato finché non ne trovo uno dentro il cerchio di raggio L
            cos=x2/sqrt(pow(x2, 2.)+pow(y2-y1, 2.));        //una volta trovato calcolo il coseno dell'angolo di cui è inclinato lo stuzzicadente
            if (y2>y1) {                                    //calcolo quindi la coordinata y del secondo estremo
               y2=y1+L*abs(cos);
            }
            else {
               y2=y1-abs(cos);
            }
            if ((y2>d)||(y2<0.)) {           //se lo stuzzicadente è di traverso sulla griglia lo conteggio
               temp++;
            }
         }                                   //qui finisce un esperimento
         ave[i]+=2.*L*N/temp/d;              //calcolo la stima di pi greco per l'esperimento appena finito
      }
      ave[i]=ave[i]/Nesp;                    //calcolo la media delle stime di tutti gli esperimenti in un blocco
      av2.push_back(pow(ave[i], 2.));        //e il suo quadrato
   }

   for (int i=0; i<M; i++) {
      ave_prog.push_back(accumulate(ave.begin(), ave.begin()+i+1, 0.)/(i+1));    //Riempio l'i-esimo valore di 'ave_prog' con la media dei valori di 'ave' fino a i compreso
      av2_prog.push_back(accumulate(av2.begin(), av2.begin()+i+1, 0.)/(i+1));    //Riempio l'i-esimo valore di 'av2_prog' con la media dei valori di 'av2' fino a i compreso
      if (i==0)
        err_prog.push_back(0.);
      else
        err_prog.push_back(sqrt((av2_prog[i] - pow(ave_prog[i], 2.))/i));        //Riempio l'i-esimo valore di 'err_prog' con la deviazione standard della media dei valori di 'ave_prog'
   }

   ofstream out;
   out.open("pigreco.txt");

   for (int i=0; i<M; i++) {        //Stampo su file '# di lanci totali fatti per la stima' 'valore stimato' 'errore stimato'
      out << i*Nesp*N << " " << ave_prog[i] << " " << err_prog[i] << endl;
   }

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
