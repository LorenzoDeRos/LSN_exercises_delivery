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
#include "random.h"
#include "genetic.h"

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

   

   int ncity=32;
   int norg=100;
   int ngenerations=5000;
   double mutation_probability=0.10;
   double crossover_probability=0.8;

   int whichmodality;
   cout << "How do you want to build your starting organism?" << endl << endl;
   cout << "1) From a new completely random organism on a circumference" << endl;
   cout << "2) From a new completely random organism inside a square" << endl << endl;
   cout << "Please write the number of your answer" << endl;
   cin >> whichmodality;

   double x, y;
   int label;
   organism Adam;                   // organismo di partenza
   ofstream out_evolution;
   out_evolution.open("genetic_path_evolution.out");
   string file="genetic_best_path.out";

// Vengono estratte le coordinate delle 32 città su una circonferenza di raggio 1
// Senza perdità di generalità, la città di partenza e di arrivo ha le coordinate (1, 0) e label 1
// Ad ogni città viene assegnato un label in ordine di costruzione
   if(whichmodality==1) {
      x=1.;
      y=0.;
      label=1;

      Adam.newcity(x, y, label);   // prima città

      double radius=1.;
      double theta;

      for(label=2; label<=ncity; label++) {     // costruzione delle altre città
         theta=rnd.Rannyu(0., 2.*M_PI);
         x=radius*cos(theta);
         y=radius*sin(theta);
         Adam.newcity(x, y, label);
      }
      // starting_organism.print();
   }

// Vengono estratte le coordinate delle 32 città con distribuzione di probabilità uniforme
// in un quadrato di lato 1 e vertici (0, 0), (0, 1), (1, 1), (1, 0)
// senza perdità di generalità, la città di partenza e di arrivo non verrà modificata nell'algoritmo
// e ha label 1
   if(whichmodality==2) {
      x=rnd.Rannyu();
      y=rnd.Rannyu();
      label=1;

      Adam.newcity(x, y, label);      //città di partenza e di arrivo

      for(label=2; label<=ncity; label++) {        //costruzione altre città
         x=rnd.Rannyu();
         y=rnd.Rannyu();
         Adam.newcity(x, y, label);
      }
      // Adam.print();
   }

// Creo una popolazione iniziale di organismi tutti uguali ad Adam
   population old_pop;
   for(int iorg=0; iorg<norg; iorg++) {
      old_pop.neworganism(Adam);
   }

// Creo una biodiversità iniziale della popolazione
   old_pop.shuffle();

// Ordino gli organismi dal più adatto al meno adatto all'ambiente
   old_pop.popsort();

// Controllo di avere una popolazione secondo le ipotesi
   if(old_pop.check()) cout << "Initialized correctly" << endl;
   else cout << "BIG PROBLEM!" << endl;

// Questa sarà la popolazione in cui verrà salvata la nuova generazione
   population new_pop;

   for(int igeneration=0; igeneration<ngenerations; igeneration++) {
         if(igeneration%100==0) cout << igeneration << " generation" << endl;
         int index1, index2;                 // indice del primo e del secondo genitore all'interno della popolazione
         int cut_index;                      // indice del punto di taglio della sequenza genetica dei genitori
         for(int iorg=0; iorg<norg-1; iorg+=2) {         // vengono estratte in tutto 50 coppie di genitori
            index1=old_pop.selection(rnd.Rannyu());      // estrazione del primo genitore
//          cout << "first parent " << index1 << endl;
//          old_pop.print(index1);
            index2=old_pop.selection(rnd.Rannyu());      // estrazione del secondo genitore
//          cout << "second parent " << index2 << endl;
//          old_pop.print(index2);
            cut_index=int(rnd.Rannyu(2., ncity-1.));     // estrazione dell'indice di taglio della sequenza genetica
//          cout << "cut index " << cut_index << endl;
            if((index1!=index2)&&(rnd.Rannyu()<crossover_probability)) {
               old_pop.crossover(index1, index2, cut_index);   // crossover, i figli sono messi in degli organismi temporanei
//             cout << "crossover occured" << endl;
               new_pop.neworganism(old_pop.getson1());         // sono copiati i figli dagli organismi temporanei nella vecchia popolazione
               new_pop.neworganism(old_pop.getson2());         // all'interno della nuova popolazione
            }
            else {      // se il crossover non avviene, i figli sono i cloni dei genitori
               new_pop.neworganism(old_pop.getorganism(index1));
               new_pop.neworganism(old_pop.getorganism(index2));
            }
//          cout << "first son " << index1 << endl;
//          new_pop.print(iorg);
//          cout << "second son " << index2 << endl;
//          new_pop.print(iorg+1);
            if(rnd.Rannyu()<mutation_probability) {      // il primo figlio subisce la prima mutazione
//             cout << "first son mutation 1 occured" << endl;
//             new_pop.print(iorg);
               new_pop.mutation1(iorg, rnd.Rannyu(), rnd.Rannyu());
//             new_pop.print(iorg);
            }
            if(rnd.Rannyu()<mutation_probability) {      // il secondo figlio subisce la prima mutazione
//             cout << "second son mutation 1 occured" << endl;
//             new_pop.print(iorg+1);
               new_pop.mutation1(iorg+1, rnd.Rannyu(), rnd.Rannyu());
//             new_pop.print(iorg+1);
            }
            if(rnd.Rannyu()<mutation_probability) {      // il primo figlio subisce la seconda mutazione
//             cout << "first son mutation 2 occured" << endl;
//             new_pop.print(iorg);
               new_pop.mutation2(iorg, rnd.Rannyu(), rnd.Rannyu(), rnd.Rannyu());
//             new_pop.print(iorg);
            }
            if(rnd.Rannyu()<mutation_probability) {      // il secondo figlio subisce la seconda mutazione
//             cout << "second son mutation 2 occured" << endl;
//             new_pop.print(iorg+1);
               new_pop.mutation2(iorg+1, rnd.Rannyu(), rnd.Rannyu(), rnd.Rannyu());
//             new_pop.print(iorg+1);
            }
            if(rnd.Rannyu()<mutation_probability) {      // il primo figlio subisce la terza mutazione
//             cout << "first son mutation 3 occured" << endl;
//             new_pop.print(iorg);
               new_pop.mutation3(iorg, rnd.Rannyu(), rnd.Rannyu());
//             new_pop.print(iorg);
            }
            if(rnd.Rannyu()<mutation_probability) {      // il secondo figlio subisce la terza mutazione
//             cout << "second son mutation 3 occured" << endl;
//             new_pop.print(iorg+1);
               new_pop.mutation3(iorg+1, rnd.Rannyu(), rnd.Rannyu());
//             new_pop.print(iorg+1);
            }
            if(rnd.Rannyu()<mutation_probability) {      // il primo figlio subisce la quarta mutazione
//             cout << "first son mutation 4 occured" << endl;
//             new_pop.print(iorg);
               new_pop.mutation4(iorg, rnd.Rannyu(), rnd.Rannyu(), rnd.Rannyu());
//             new_pop.print(iorg);
            }
            if(rnd.Rannyu()<mutation_probability) {      // il secondo figlio subisce la quarta mutazione
//             cout << "second son mutation 4 occured" << endl;
//             new_pop.print(iorg+1);
               new_pop.mutation4(iorg+1, rnd.Rannyu(), rnd.Rannyu(), rnd.Rannyu());
//             new_pop.print(iorg+1);
            }
            if(rnd.Rannyu()<mutation_probability) {      // il primo figlio subisce la quinta mutazione
//             cout << "first son mutation 5 occured" << endl;
//             new_pop.print(iorg);
               new_pop.mutation5(iorg, rnd.Rannyu(), rnd.Rannyu(), rnd.Rannyu());
//             new_pop.print(iorg);
            }
            if(rnd.Rannyu()<mutation_probability) {      // il secondo figlio subisce la quinta mutazione
//             cout << "second son mutation 5 occured" << endl;
//             new_pop.print(iorg+1);
               new_pop.mutation5(iorg+1, rnd.Rannyu(), rnd.Rannyu(), rnd.Rannyu());
//             new_pop.print(iorg+1);
            }
            old_pop.clearsons();       // vengono liberati gli organismi temporanei della vecchia popolazione
         }
         if(!new_pop.check()) cout << "BIG PROBLEM! CROSSOVER" << endl;    // la nuova popolazione viene controllata
         old_pop.clear();              // la vecchia popolazione viene liberata
         new_pop.popsort();            // la nuova popolazione viene riordinata in base all'idoneità all'ambiente
         out_evolution << (new_pop.getorganism(0)).length() << endl;    // si stampa la lunghezza dell'organismo migliore
         for(int jorg=0; jorg<norg; jorg++) {   // la nuova popolazione diventa la vecchia popolazione
            old_pop.neworganism(new_pop.getorganism(jorg));
         }
         new_pop.clear();              // la nuova popolazione viene liberata
   }

   (old_pop.getorganism(0)).fileprint(file);       // viene stampato l'organismo migliore dell'ultima generazione

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
