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

double Function(const Posizione &P, int whichfunction) {
   if(whichfunction==1) {
      return exp(-2*P.getR())/M_PI;
   }
   else {
      return pow(P.getZ(), 2.)*exp(-P.getR())/4/M_PI;
   }
}

double Accept(Posizione &P, Posizione &P_next, int whichfunction) {
   return min(1., Function(P_next, whichfunction)/Function(P, whichfunction));
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

   Posizione P1(0.75, 0., 0.);   //first function start
   //Posizione P1(75., 0., 0.);  //first function far away start
   Posizione P2(3.5, 0., 0.);    //second function start
   Posizione P;
   int M=pow(10, 6);             //total number of evaluations
   int Nblocks=100;              //number of groups
   int Nsteps=M/Nblocks;         //number of evaluations in every group (steps)
   double L, L1=2.5, L2=6.;      //lato del cubo in cui viene estratta uniformemente la posizione
   double sigma, sigma1=0.75, sigma2=1.85;
   double X, Y, Z;
   double r, r_prog, r2_prog, err_r_prog;
   double acc, acc_final;
   double alpha;
   vector<double> ave_r;         //si inseriscono le stime del valor medio della prima distribuzione di probabilità
                                 //una stima per ogni blocco
   vector<double> av2_r;
   vector<double> ave_acc;

   for(int i=0; i<Nblocks; i++) {
      ave_r.push_back(0.);
      av2_r.push_back(0.);
      ave_acc.push_back(0.);
   }

   

   cout << "PRIMA FUNZIONE D'ONDA" << endl << endl;
   cout << "AVANZAMENTO UNIFORME IN CUBO" << endl;
   cout << "le prime 10 misure sono" << endl;

   ofstream coordinates_f1_unif;
   coordinates_f1_unif.open("coordinates_f1_unif.txt");
   ofstream measures_f1_unif;
   measures_f1_unif.open("prog_r_f1_unif.txt");
   ofstream acceptance_f1_unif;
   acceptance_f1_unif.open("acc_f1_unif.txt");

   //ofstream out;
   //out.open("starting_far_away.txt");

   int whichfunction=1;
   P=P1;
   L=L1;

   for(int iblock=0; iblock<Nblocks; iblock++) {
      r=0.;
      acc=0.;
      for(int j=0; j<Nsteps; j++) {
         X=rnd.Rannyu(P.getX()-L/2, P.getX()+L/2);
         Y=rnd.Rannyu(P.getY()-L/2, P.getY()+L/2);
         Z=rnd.Rannyu(P.getZ()-L/2, P.getZ()+L/2);
         Posizione P_next(X, Y, Z);
         alpha=Accept(P, P_next, whichfunction);
         acc+=alpha;
         if(rnd.Rannyu()<alpha) {
            P=P_next;
         }
         if(iblock%10==0&&j%100==0) {
            coordinates_f1_unif << P.getZ() << " " << P.getY() << " " << P.getZ() << endl;
         }
         //if(iblock==0&&j<1000) {
         //   out << j+1 << " " << P.getR() << " " << alpha << endl;
         //}
         r+=P.getR();
      }
      r=r/Nsteps;
      acc=acc/Nsteps;
      ave_r[iblock]=r;
      av2_r[iblock]=r*r;
      ave_acc[iblock]=acc;
   }

   //out.close();

   acc_final=accumulate(ave_acc.begin(), ave_acc.end(), 0.)/Nblocks;

   for(int iblock=0; iblock<Nblocks; iblock++) {
      r_prog=accumulate(ave_r.begin(), ave_r.begin()+iblock+1, 0.)/(iblock+1);
      r2_prog=accumulate(av2_r.begin(), av2_r.begin()+iblock+1, 0.)/(iblock+1);
      if(iblock==0) {
         err_r_prog=0.;
      }
      else {
         err_r_prog=sqrt((r2_prog-pow(r_prog, 2.))/iblock);
      }
      measures_f1_unif << (iblock+1)*Nsteps << " " << r_prog << " " << err_r_prog << endl;
      acceptance_f1_unif << iblock+1 << " " << ave_acc[iblock] << endl;
   }

   coordinates_f1_unif.close();
   measures_f1_unif.close();
   acceptance_f1_unif.close();

   cout << "posizione media finale" << endl;
   cout << r_prog << " pm " << err_r_prog << endl;
   cout << "valor medio finale della probabilità di accettazione" << endl;
   cout << acc_final << endl << endl;

   

   cout << "AVANZAMENTO GAUSSIANO" << endl;
   cout << "le prime 10 misure sono" << endl;

   ofstream coordinates_f1_Gauss;
   coordinates_f1_Gauss.open("coordinates_f1_Gauss.txt");
   ofstream measures_f1_Gauss;
   measures_f1_Gauss.open("prog_r_f1_Gauss.txt");
   ofstream acceptance_f1_Gauss;
   acceptance_f1_Gauss.open("acc_f1_Gauss.txt");

   whichfunction=1;
   P=P1;
   sigma=sigma1;

   for(int iblock=0; iblock<Nblocks; iblock++) {
      r=0.;
      acc=0.;
      for(int j=0; j<Nsteps; j++) {
         X=rnd.Gauss(P.getX(), sigma);
         Y=rnd.Gauss(P.getY(), sigma);
         Z=rnd.Gauss(P.getZ(), sigma);
         Posizione P_next(X, Y, Z);
         alpha=Accept(P, P_next, whichfunction);
         acc+=alpha;
         if(rnd.Rannyu()<alpha) {
            P=P_next;
         }
         if(iblock%10==0&&j%100==0) {
            coordinates_f1_Gauss << P.getZ() << " " << P.getY() << " " << P.getZ() << endl;
         }
         if(iblock==0&&j<10) {
            cout << P.getR() << endl;
         }
         r+=P.getR();
      }
      r=r/Nsteps;
      acc=acc/Nsteps;
      ave_r[iblock]=r;
      av2_r[iblock]=r*r;
      ave_acc[iblock]=acc;
   }

   acc_final=accumulate(ave_acc.begin(), ave_acc.end(), 0.)/Nblocks;

   for(int iblock=0; iblock<Nblocks; iblock++) {
      r_prog=accumulate(ave_r.begin(), ave_r.begin()+iblock+1, 0.)/(iblock+1);
      r2_prog=accumulate(av2_r.begin(), av2_r.begin()+iblock+1, 0.)/(iblock+1);
      if(iblock==0) {
         err_r_prog=0.;
      }
      else {
         err_r_prog=sqrt((r2_prog-pow(r_prog, 2.))/iblock);
      }
      measures_f1_Gauss << (iblock+1)*Nsteps << " " << r_prog << " " << err_r_prog << endl;
      acceptance_f1_Gauss << iblock+1 << " " << ave_acc[iblock] << endl;
   }

   coordinates_f1_Gauss.close();
   measures_f1_Gauss.close();
   acceptance_f1_Gauss.close();

   cout << "posizione media finale" << endl;
   cout << r_prog << " pm " << err_r_prog << endl;
   cout << "valor medio finale della probabilità di accettazione" << endl;
   cout << acc_final << endl << endl;



   cout << endl << endl;
   cout << "SECONDA FUNZIONE D'ONDA" << endl << endl;
   cout << "AVANZAMENTO IN UN CUBO" << endl;
   cout << "le prime 10 misure sono" << endl;

   ofstream coordinates_f2_unif;
   coordinates_f2_unif.open("coordinates_f2_unif.txt");
   ofstream measures_f2_unif;
   measures_f2_unif.open("prog_r_f2_unif.txt");
   ofstream acceptance_f2_unif;
   acceptance_f2_unif.open("acc_f2_unif.txt");

   whichfunction=2;
   P=P2;
   L=L2;

   for(int iblock=0; iblock<Nblocks; iblock++) {
      r=0.;
      acc=0.;
      for(int j=0; j<Nsteps; j++) {
         X=rnd.Rannyu(P.getX()-L/2, P.getX()+L/2);
         Y=rnd.Rannyu(P.getY()-L/2, P.getY()+L/2);
         Z=rnd.Rannyu(P.getZ()-L/2, P.getZ()+L/2);
         Posizione P_next(X, Y, Z);
         alpha=Accept(P, P_next, whichfunction);
         acc+=alpha;
         if(rnd.Rannyu()<alpha) {
            P=P_next;
         }
         if(iblock%10==0&&j%100==0) {
            coordinates_f2_unif << P.getZ() << " " << P.getY() << " " << P.getZ() << endl;
         }
         if(iblock==0&&j<10) {
            cout << P.getR() << endl;
         }
         r+=P.getR();
      }
      r=r/Nsteps;
      acc=acc/Nsteps;
      ave_r[iblock]=r;
      av2_r[iblock]=r*r;
      ave_acc[iblock]=acc;
   }

   acc_final=accumulate(ave_acc.begin(), ave_acc.end(), 0.)/Nblocks;

   for(int iblock=0; iblock<Nblocks; iblock++) {
      r_prog=accumulate(ave_r.begin(), ave_r.begin()+iblock+1, 0.)/(iblock+1);
      r2_prog=accumulate(av2_r.begin(), av2_r.begin()+iblock+1, 0.)/(iblock+1);
      if(iblock==0) {
         err_r_prog=0.;
      }
      else {
         err_r_prog=sqrt((r2_prog-pow(r_prog, 2.))/iblock);
      }
      measures_f2_unif << (iblock+1)*Nsteps << " " << r_prog << " " << err_r_prog << endl;
      acceptance_f2_unif << iblock+1 << " " << ave_acc[iblock] << endl;
   }

   coordinates_f2_unif.close();
   measures_f2_unif.close();
   acceptance_f2_unif.close();

   cout << "posizione media finale" << endl;
   cout << r_prog << " pm " << err_r_prog << endl;
   cout << "valor medio finale della probabilità di accettazione" << endl;
   cout << acc_final << endl << endl;
   



   cout << "AVANZAMENTO GAUSSIANO" << endl;
   cout << "le prime 10 misure sono" << endl;

   ofstream coordinates_f2_Gauss;
   coordinates_f2_Gauss.open("coordinates_f2_Gauss.txt");
   ofstream measures_f2_Gauss;
   measures_f2_Gauss.open("prog_r_f2_Gauss.txt");
   ofstream acceptance_f2_Gauss;
   acceptance_f2_Gauss.open("acc_f2_Gauss.txt");

   whichfunction=2;
   P=P2;
   sigma=sigma2;

   for(int iblock=0; iblock<Nblocks; iblock++) {
      r=0.;
      acc=0.;
      for(int j=0; j<Nsteps; j++) {
         X=rnd.Gauss(P.getX(), sigma);
         Y=rnd.Gauss(P.getY(), sigma);
         Z=rnd.Gauss(P.getZ(), sigma);
         Posizione P_next(X, Y, Z);
         alpha=Accept(P, P_next, whichfunction);
         acc+=alpha;
         if(rnd.Rannyu()<alpha) {
            P=P_next;
         }
         if(iblock%10==0&&j%100==0) {
            coordinates_f2_Gauss << P.getZ() << " " << P.getY() << " " << P.getZ() << endl;
         }
         if(iblock==0&&j<10) {
            cout << P.getR() << endl;
         }
         r+=P.getR();
      }
      r=r/Nsteps;
      acc=acc/Nsteps;
      ave_r[iblock]=r;
      av2_r[iblock]=r*r;
      ave_acc[iblock]=acc;
   }

   acc_final=accumulate(ave_acc.begin(), ave_acc.end(), 0.);

   for(int iblock=0; iblock<Nblocks; iblock++) {
      r_prog=accumulate(ave_r.begin(), ave_r.begin()+iblock+1, 0.)/(iblock+1);
      r2_prog=accumulate(av2_r.begin(), av2_r.begin()+iblock+1, 0.)/(iblock+1);
      if(iblock==0) {
         err_r_prog=0.;
      }
      else {
         err_r_prog=sqrt((r2_prog-pow(r_prog, 2.))/iblock);
      }
      measures_f2_Gauss << (iblock+1)*Nsteps << " " << r_prog << " " << err_r_prog << endl;
      acceptance_f2_Gauss << iblock+1 << " " << ave_acc[iblock] << endl;
   }

   coordinates_f2_Gauss.close();
   measures_f2_Gauss.close();
   acceptance_f2_Gauss.close();

   cout << "posizione media finale" << endl;
   cout << r_prog << " pm " << err_r_prog << endl;
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
