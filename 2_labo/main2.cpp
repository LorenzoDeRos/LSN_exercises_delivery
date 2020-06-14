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
#include "posizione.h"

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

   Posizione O(0., 0., 0.);
   double passo=1.;
   int M=10000;                        //numero totale di random walkers
   int G=100;                          //numero di gruppi/blocchi di random walkers
   int N=M/G;                          //numero di random walkers in ogni gruppo
   int P=100;                          //numero di passi di ogni random walker
   vector<double> dist;                //vector contenente il valor medio della distanza all'aumentare dei passi per ogni gruppo
   vector<double> dist_media;          //vector contenente il valor medio della distanza di tutti i gruppi
   vector<double> dist2_media;         //vector contenente il valor medio del quadrato della distanza di tutti i gruppi
   vector<double> err;                 //vector contenente la deviazione standard della media della distanza di tutti i gruppi

   for(int i=0; i<P; i++) {            //costruisco i vectors
      dist.push_back(0.);
      err.push_back(0.);
      dist_media.push_back(0.);
      dist2_media.push_back(0.);
   }

   //PRIMA PARTE
   //Random walkers su un reticolo discreto di passo=1

   for (int i=0; i<G; i++) {                                   //itero sul numero di gruppi di random walkers
      for(int j=0; j<N; j++) {                                 //itero sul numero di random walkers in ogni gruppo
         Posizione R(0., 0., 0.);                              //punto di partenza del random walker
         for(int k=0; k<P; k++) {                              //itero sul numero di passi compiuti da ogni random walker
            R.AvanzaReticolo(passo, rnd.Asse(), rnd.Verso());  //il random walker avanza di un passo lungo una direzione e verso qualunque lungo il reticolo
            dist[k]+=pow(R.Distanza(O), 2.);                   //accumulo nel k-esimo elemento di 'dist' la distanza dall'origine al quadrato alla quale ogni random walker di un gruppo giunge
         }
      }
      for(int l=0; l<P; l++) {                                 //itero sul numero di passi
         dist[l]=sqrt(dist[l]/N);                              //radice del valor medio della (distanza dall'origine al quadrato) per il gruppo appena considerato
         dist_media[l]+=dist[l];                               //accumulo il valor medio della distanza in funzione del numero di passi di ogni gruppo
         dist2_media[l]+=pow(dist[l], 2.);                     //accumulo il quadrato del valor medio della distanza in funzione del numero di passi di ogni gruppo
         dist[l]=0.;                                           //azzero il valor medio della distanza dall'origine del gruppo considerato
      }                                                        //per poi riempirlo con quello del gruppo successivo
   }

   for(int i=0; i<P; i++) {                                    //itero sul numero di passi
      dist_media[i]=dist_media[i]/G;                           //valor medio della media della distanza in funzione del numero di passi di ogni gruppo
      dist2_media[i]=dist2_media[i]/G;                         //valor medio della media al quadrato della distanza in funzione del numero di passi di ogni gruppo
      if (i==0) {
         err[i]=0.;
      }
      else {
         err[i]=sqrt((dist2_media[i]-pow(dist_media[i], 2.))/(G-1.));   //deviazione standard del valor medio della media della distanza in funzione del numero di passi di ogni gruppo
      }
   }

   ofstream outlattice;
   outlattice.open("lattice_walks.txt");
   for(int i=0; i<P; i++) {
         outlattice << i+1 << " " << dist_media[i] << " " << err[i] << endl;     //Stampo 'passi compiuti' 'stima finale della distanza dall'origine' 'deviazione standard sulla stima finale'
         dist_media[i]=0.;                //azzero per riutilizzare gli stessi vectors dopo
         dist2_media[i]=0.;
         err[i]=0.;
   }

   outlattice.close();


   //SECONDA PARTE
   //Random walkers in una direzione arbitraria in 3D che fa passi di lunghezza 1

   //L'implementazione è identica al caso precedente con l'unica differenza nella funzione avanzamento

   for (int i=0; i<G; i++) {
      for(int j=0; j<N; j++) {
         Posizione R(0., 0., 0.);
         for(int k=0; k<P; k++) {
            R.AvanzaDirezione(passo, rnd.Longitudine(), rnd.Latitudine());    //Questa è l'unica riga differente
            dist[k]+=pow(R.Distanza(O), 2.);
         }
      }
      for(int l=0; l<P; l++) {
         dist[l]=sqrt(dist[l]/N);
         dist_media[l]+=dist[l];
         dist2_media[l]+=pow(dist[l], 2.);
         dist[l]=0.;
      }      
   }

   for(int i=0; i<P; i++) {
      dist_media[i]=dist_media[i]/G;
      dist2_media[i]=dist2_media[i]/G;
      if (i==0) {
         err[i]=0.;
      }
      else {
         err[i]=sqrt((dist2_media[i]-pow(dist_media[i], 2.))/(G-1.));
      }
   }

   ofstream outdirection;
   outdirection.open("random_direction_walks.txt");
   for(int i=0; i<P; i++) {
         outdirection << i+1 << " " << dist_media[i] << " " << err[i] << endl;
   }

   outdirection.close();

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
