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
#include "mpi.h"

using namespace std;
 
int main (int argc, char *argv[]){

   int size, rank, itag=0;
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Status stat1, stat2, stat3;

   Random rnd;
   int seed[4];
   int p1, p2, temp1, temp2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      for(int irank=0; irank<size; irank++) {
         Primes >> temp1 >> temp2;
         if(rank==irank) {
            p1=temp1;
            p2=temp2;
         }
      }
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
   int ngenerations=1000;
   int nprocess=4;

   if(size!=nprocess) cout << "Wrong setting of the number of processes!" << endl;

   int imigration=50;
   double mutation_probability=0.10;
   double crossover_probability=0.8;
   int whichmodality=2;

   double x, y;
   int label;

   int labeltransf[ncity];
   int labelsend[ncity];
   int labelreceive[ncity];
   double xtransf[ncity];
   double xsend[ncity];
   double xreceive[ncity];
   double ytransf[ncity];
   double ysend[ncity];
   double yreceive[ncity];
   int order[nprocess+1];

   organism Adam;
   organism explorer;
   ofstream out_evolution;
   out_evolution.open("parallel_genetic_path_evolution" + to_string(rank) + ".out");
   string file="parallel_genetic_best_path" + to_string(rank) + ".out";

// L'organismo di partenza viene costruito nel processo con rank 0
if(rank==0) {
   
   if(whichmodality==1) {
      x=1.;
      y=0.;
      label=1;

      labeltransf[0]=label;
      xtransf[0]=x;
      ytransf[0]=y;

      double radius=1.;
      double theta;

      for(label=2; label<=ncity; label++) {
         theta=rnd.Rannyu(0., 2.*M_PI);
         x=radius*cos(theta);
         y=radius*sin(theta);
         labeltransf[label-1]=label;
         xtransf[label-1]=x;
         ytransf[label-1]=y;
      }

   }

   if(whichmodality==2) {
      x=rnd.Rannyu();
      y=rnd.Rannyu();
      label=1;
      labeltransf[0]=label;
      xtransf[0]=x;
      ytransf[0]=y;

      for(label=2; label<=ncity; label++) {
         x=rnd.Rannyu();
         y=rnd.Rannyu();
         labeltransf[label-1]=label;
         xtransf[label-1]=x;
         ytransf[label-1]=y;
      }
   }

}

// L'organismo di partenza viene poi condiviso agli altri processori,
// in modo che tutti i continenti abbiano organismi costituiti dalle stesse città
   MPI_Bcast(labeltransf,ncity,MPI_INT,0, MPI_COMM_WORLD);
   MPI_Bcast(xtransf,ncity,MPI_DOUBLE,0, MPI_COMM_WORLD);
   MPI_Bcast(ytransf,ncity,MPI_DOUBLE,0, MPI_COMM_WORLD);

   for(int icity=0; icity<ncity; icity++) {
      Adam.newcity(xtransf[icity], ytransf[icity], labeltransf[icity]);
   }

// Ogni continenti costruisce la propria popolazione iniziale
   population old_pop;
   for(int iorg=0; iorg<norg; iorg++) {
      old_pop.neworganism(Adam);
   }
// La popolazione iniziale viene diversificata
   old_pop.shuffle();
// e ordinata dall'organismo migliore al peggiore
   old_pop.popsort();
   if(old_pop.check()) cout << "initialized correctly" << endl;
   else cout << "BIG PROBLEM!" << endl;
// Questa popolazione conterrà la futura generazione
   population new_pop;



      for(int igeneration=0; igeneration<ngenerations; igeneration++) {

         int index1, index2, cut_index;
         if(igeneration%100==0) cout << igeneration << " generation" << endl;
         for(int iorg=0; iorg<norg-1; iorg+=2) {
            index1=old_pop.selection(rnd.Rannyu());
            index2=old_pop.selection(rnd.Rannyu());
            cut_index=int(rnd.Rannyu(2., ncity-1.));
            if((index1!=index2)&&(rnd.Rannyu()<crossover_probability)) {
               old_pop.crossover(index1, index2, cut_index);
               new_pop.neworganism(old_pop.getson1());
               new_pop.neworganism(old_pop.getson2());
            }
            else {
               new_pop.neworganism(old_pop.getorganism(index1));
               new_pop.neworganism(old_pop.getorganism(index2));
            }
            if(rnd.Rannyu()<mutation_probability) {
               new_pop.mutation1(iorg, rnd.Rannyu(), rnd.Rannyu());
            }
            if(rnd.Rannyu()<mutation_probability) {
               new_pop.mutation1(iorg+1, rnd.Rannyu(), rnd.Rannyu());
            }
            if(rnd.Rannyu()<mutation_probability) {
               new_pop.mutation2(iorg, rnd.Rannyu(), rnd.Rannyu(), rnd.Rannyu());
            }
            if(rnd.Rannyu()<mutation_probability) {
               new_pop.mutation2(iorg+1, rnd.Rannyu(), rnd.Rannyu(), rnd.Rannyu());
            }
            if(rnd.Rannyu()<mutation_probability) {
               new_pop.mutation3(iorg, rnd.Rannyu(), rnd.Rannyu());
            }
            if(rnd.Rannyu()<mutation_probability) {
               new_pop.mutation3(iorg+1, rnd.Rannyu(), rnd.Rannyu());
            }
            if(rnd.Rannyu()<mutation_probability) {
               new_pop.mutation4(iorg, rnd.Rannyu(), rnd.Rannyu(), rnd.Rannyu());
            }
            if(rnd.Rannyu()<mutation_probability) {
               new_pop.mutation4(iorg+1, rnd.Rannyu(), rnd.Rannyu(), rnd.Rannyu());
            }
            if(rnd.Rannyu()<mutation_probability) {
               new_pop.mutation5(iorg, rnd.Rannyu(), rnd.Rannyu(), rnd.Rannyu());
            }
            if(rnd.Rannyu()<mutation_probability) {
               new_pop.mutation5(iorg+1, rnd.Rannyu(), rnd.Rannyu(), rnd.Rannyu());
            }
            old_pop.clearsons();
         }
         if(!new_pop.check()) cout << "BIG PROBLEM! CROSSOVER" << endl;
         old_pop.clear();
         new_pop.popsort();
         out_evolution << (new_pop.getorganism(0)).length() << endl;
         for(int jorg=0; jorg<norg; jorg++) old_pop.neworganism(new_pop.getorganism(jorg));
         new_pop.clear();
//       MIGRAZIONE: avviene una volta ogni imigration
         if((igeneration%imigration==0)&&(igeneration!=0)) {
//          Il processore 0 costruisce in modo casuale l'ordine con cui verranno scambiati i migliori organismi
            if(rank==0) {
               vector<int> ordergen;
               for(int i=0; i<size; i++) ordergen.push_back(i);
               ordergen.push_back(0);
               random_shuffle(ordergen.begin()+1, ordergen.end()-1);
               for(int i=0; i<=size; i++) order[i]=ordergen[i];
            }
//          e l'ordine viene poi condiviso agli altri continenti
//          in modo che ogni continente sappia da chi lo riceverà, e a chi lo manderà
            MPI_Bcast(order,size+1,MPI_INT,0, MPI_COMM_WORLD);
//          Ogni continente prepara l'organismo da inviare
            for(int icity=0; icity<ncity; icity++) {
               labelsend[icity]=(old_pop.getorganism(0).getcity(icity)).getlabel();
               xsend[icity]=(old_pop.getorganism(0).getcity(icity)).getx();
               ysend[icity]=(old_pop.getorganism(0).getcity(icity)).gety();
            }
            for(int iorder=0; iorder<size; iorder++) {
//             Ogni continente effettua l'invio al continente successivo nell'ordine
               if(rank==order[iorder]) {
                  itag=1;
                  MPI_Send(labelsend, ncity, MPI_INT, order[iorder+1], itag, MPI_COMM_WORLD);
                  itag=2;
                  MPI_Send(xsend, ncity, MPI_DOUBLE, order[iorder+1], itag, MPI_COMM_WORLD);
                  itag=3;
                  MPI_Send(ysend, ncity, MPI_DOUBLE, order[iorder+1], itag, MPI_COMM_WORLD);
               }
//             Ogni continente successivo attende di ricevere dal precedente nell'ordine
               if(rank==order[iorder+1]) {
                  itag=1;
                  MPI_Recv(labelreceive, ncity, MPI_INT, order[iorder], itag, MPI_COMM_WORLD, &stat1);
                  itag=2;
                  MPI_Recv(xreceive, ncity, MPI_DOUBLE, order[iorder], itag, MPI_COMM_WORLD, &stat2);
                  itag=3;
                  MPI_Recv(yreceive, ncity, MPI_DOUBLE, order[iorder], itag, MPI_COMM_WORLD, &stat3);
               }
            }
//          Ogni continente salva quanto ricevuto
            for(int icity=0; icity<ncity; icity++) {
               explorer.newcity(xreceive[icity], yreceive[icity], labelreceive[icity]);
            }
//          e lo sostituisce al proprio miglior organismo
            old_pop.setorganism(0, explorer);
            explorer.clear();
         }

      }


   (old_pop.getorganism(0)).fileprint(file);





   MPI_Finalize();

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
