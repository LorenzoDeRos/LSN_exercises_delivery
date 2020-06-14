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
#include <ostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 
  //ofstream outEvo;
  //outEvo.open("magnetization_evolution.dat");
  cout << "how do you want to use this program?" << endl << endl;
  cout << "please answer with the sentence number" << endl;
  cout << "1) starting from a random spin configuration" << endl;
  cout << "2) starting from an equilibrated configuration at the chosen temperature" << endl;
  cin >> answer;
  Input(); //Inizialization
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
      //if(iblk==1&&istep<500) {
      //  outEvo << istep << " " << walker[im]/(double)m_spin << endl;
      //}
    }
    Averages(iblk);   //Print results for current block
  }
  if(answer==1) ConfFinal(); //Write final configuration
  if(answer==2) DataFinal(); //Write the best estimates

  return 0;
}


void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
  if(answer==1) {
    for (int i=0; i<nspin; ++i)
    {
      if(rnd.Rannyu() >= 0.5) s[i] = 1;
      else s[i] = -1;
    }
  }
  else {
    ifstream ReadConfig;
    ReadConfig.open("config_" + to_string(temp) + ".dat");
    for(int i=0; i<nspin; ++i) {
      ReadConfig >> s[i];
    }
    
  }
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro)
{
  int o;
  double p, energy_old, energy_new, sm;
  double energy_up, energy_down;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);
    attempted++;

    if(metro==1) //Metropolis
    {
      energy_old=Boltzmann(s[o], o);
      sm=-s[o];
      energy_new=Boltzmann(sm, o);
      p=min(1., exp(-beta*(energy_new-energy_old)));
      if(rnd.Rannyu()<p) {
        s[o]=sm;
        accepted++;
      }
    }
    else //Gibbs sampling
    {
      sm=1;
      energy_up=Boltzmann(sm, o);
      energy_down=Boltzmann(-sm, o);
      p=1/(1+exp(beta*(energy_up - energy_down)));
      accepted++;
      if (rnd.Rannyu()<p) s[o]=sm;
      else s[o]=-sm;
    }
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     m+= s[i];
  }
  walker[iu] = u;
  walker[ic] = u*u;
  walker[im] = m;
  walker[ix] = m*m;
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Ene.open("output.ene.0",ios::app);
    Heat.open("output.heat.0",ios::app);
    Mag.open("output.mag.0",ios::app);
    Chi.open("output.chi.0",ios::app);
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    stima_c = beta*beta*(blk_av[ic]/blk_norm-pow(blk_av[iu]/blk_norm, 2.))/(double)nspin; //Heat capacity
    stima_m = blk_av[im]/blk_norm/(double)nspin; //Megnatization
    stima_x = beta*blk_av[ix]/blk_norm/(double)nspin; //Magnetic susceptibility
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    glob_av[ic] += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    glob_av[im] += stima_m;
    glob_av2[im] += stima_m*stima_m;
    glob_av[ix] += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    
    Ene << iblk <<  " " << stima_u << " " << glob_av[iu]/(double)iblk << " " << err_u << endl;
    Ene.close();
    Heat << iblk << " " << stima_c << " " << glob_av[ic]/(double)iblk << " " << err_c << endl;
    Heat.close();
    Mag << iblk << " " << stima_m << " " << glob_av[im]/(double)iblk << " " << err_m << endl;
    Mag.close();
    Chi << iblk << " " << stima_x << " " << glob_av[ix]/(double)iblk << " " << err_x << endl;
    Chi.close();

    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config_" + to_string(temp) + ".dat");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

void DataFinal(void) {
  ofstream WriteEne, WriteHeat, WriteMag, WriteChi;
  cout << "Print the best estimates for the quantity calculated" << endl << endl;
  if(metro==1) {
    if(h>0.01) {
      WriteMag.open("metro_mag.dat", ios::app);
      WriteMag << temp << " " << glob_av[im]/(double)nblk << " " << err_m << endl;
      WriteMag.close();
    }
    else {
      WriteEne.open("metro_ene.dat", ios::app);
      WriteEne << temp << " " << glob_av[iu]/(double)nblk << " " << err_u << endl;
      WriteHeat.open("metro_heat.dat", ios::app);
      WriteHeat << temp << " " << glob_av[ic]/(double)nblk << " " << err_c << endl;
      WriteChi.open("metro_chi.dat", ios::app);
      WriteChi << temp << " " << glob_av[ix]/(double)nblk << " " << err_x << endl;
      WriteEne.close();
      WriteHeat.close();
      WriteChi.close();
    }
  }
  else {
    if(h>0.01) {
      WriteMag.open("gibbs_mag.dat", ios::app);
      WriteMag << temp << " " << glob_av[im]/(double)nblk << " " << err_m << endl;
      WriteMag.close();
    }
    else {
      WriteEne.open("gibbs_ene.dat", ios::app);
      WriteEne << temp << " " << glob_av[iu]/(double)nblk << " " << err_u << endl;
      WriteHeat.open("gibbs_heat.dat", ios::app);
      WriteHeat << temp << " " << glob_av[ic]/(double)nblk << " " << err_c << endl;
      WriteChi.open("gibbs_chi.dat", ios::app);
      WriteChi << temp << " " << glob_av[ix]/(double)nblk << " " << err_x << endl;
      WriteEne.close();
      WriteHeat.close();
      WriteChi.close();
    }
  }
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
  if(iblk==1)
    return 0;
  else
    return sqrt(abs((sum2/(double)iblk - pow(sum/(double)iblk,2)))/(double)iblk);
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
