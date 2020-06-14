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
   int nsteps=100000;
   int nmutations=5;
   int whichmutation;
   double starting_temperature=0.01;
   double final_temperature=0.001;
   double delta_temperature=0.001;
   double alpha;
   double accepted=0.;
   double attempted=0.;

   double x, y;
   int label;
   organism my_org;
   string file="metropolis_best_path.out";
   ofstream evolution;
   evolution.open("metropolis_path_evolution.out");

   int whichmodality;
   cout << "How do you want to build your starting organism?" << endl << endl;
   cout << "1) From the last organism I obtained" << endl;
   cout << "2) From a new completely random organism on a circumference" << endl;
   cout << "3) From a new completely random organism inside a square" << endl << endl;
   cout << "Please write the number of your answer" << endl;
   cin >> whichmodality;

// L'organismo viene inizializzato leggendo le città dal percorso finale dell'ultima simulazione
   if(whichmodality==1) {
      ifstream input;
      input.open(file);
      for(int icity=0; icity<ncity; icity++) {
         input >> label;
         input >> x;
         input >> y;
         my_org.newcity(x, y, label);
      }
      if(my_org.check()) cout << "Initialised correctly" << endl;
   }

// L'organismo viene inizializzato costruendo le città lungo una circonferenza di raggio 1
// La prima città (1, 0) è fissa
   if(whichmodality==2) {
      x=1.;
      y=0.;
      label=1;

      my_org.newcity(x, y, label);

      double radius=1.;
      double theta;

      for(label=2; label<=ncity; label++) {
         theta=rnd.Rannyu(0., 2.*M_PI);
         x=radius*cos(theta);
         y=radius*sin(theta);
         my_org.newcity(x, y, label);
      }

      if(my_org.check()) cout << "Initialised correctly" << endl;
   }

// L'organismo viene inizializzato costruendo le città all'interno di un quadrato di lato 1
// La prima città è fissa
   if(whichmodality==3) {
      x=rnd.Rannyu();
      y=rnd.Rannyu();
      label=1;

      my_org.newcity(x, y, label);

      for(label=2; label<=ncity; label++) {
         x=rnd.Rannyu();
         y=rnd.Rannyu();
         my_org.newcity(x, y, label);
      }

      if(my_org.check()) cout << "Initialised correctly" << endl;
   }

// La temperatura viene abbassata a mano a mano che il percorso si equilibra alla temperatura selezionata
for(double temperature=starting_temperature; temperature>=final_temperature; temperature-=delta_temperature) {
   cout << endl << "TEMPERATURE " << temperature << endl;
   for(int istep=0; istep<nsteps; istep++) {
      organism temp;
      for(int icity=0; icity<ncity; icity++) temp.newcity(my_org.getcity(icity));
      whichmutation=int(nmutations*rnd.Rannyu())+1;      // Estraggo la proposta di mutazione che l'organismo può subire
      if(whichmutation==1) temp.mutation1(rnd.Rannyu(), rnd.Rannyu());
      if(whichmutation==2) temp.mutation2(rnd.Rannyu(), rnd.Rannyu(), rnd.Rannyu());
      if(whichmutation==3) temp.mutation3(rnd.Rannyu(), rnd.Rannyu());
      if(whichmutation==4) temp.mutation4(rnd.Rannyu(), rnd.Rannyu(), rnd.Rannyu());
      if(whichmutation==5) temp.mutation5(rnd.Rannyu(), rnd.Rannyu(), rnd.Rannyu());
      alpha=min(1., exp(-(temp.length()-my_org.length())/temperature));    // Calcolo qual è la probabilità che questa mutazione venga accettata
      if(rnd.Rannyu()<alpha) {      // La mutazione viene accettata e l'organismo modificato
         my_org.clear();
         for(int icity=0; icity<ncity; icity++) my_org.newcity(temp.getcity(icity));
         accepted++;
      }
      if(!my_org.check()) cout << "BIG PROBLEM!" << endl;
      temp.clear();
      attempted++;
      if(istep%10==0) evolution << temperature << " " << my_org.length() << endl;
   }
   cout << "acceptance " << accepted/attempted << endl;
}
   my_org.fileprint(file);

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
