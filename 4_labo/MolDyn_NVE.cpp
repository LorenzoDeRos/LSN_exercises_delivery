/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <vector>
#include <numeric>
#include "MolDyn_NVE.h"

using namespace std;

vector<double> pot_block;
vector<double> pot2_block;
vector<double> kin_block;
vector<double> kin2_block;
vector<double> etot_block;
vector<double> etot2_block;
vector<double> temp_block;
vector<double> temp2_block;

int main(){ 
  cout << "How do you want to use this program?" << endl;
  cout << "You have 3 possibilities" << endl << endl;
  cout << "1) Start the simulation from a fcc configuration" << endl;
  cout << "2) Start the simulation where you left it, trying to correct the temperature" << endl;
  cout << "3) Start the simulation where you left it" << endl << endl;
  cout << "please, answer with the answer number" << endl;
  int answer;
  cin >> answer;
  InputParameters();             //Inizialization independent from the answer
  if(answer==1) InputStart();
  if (answer==2) InputCorrect();
  if(answer==3) InputFinal();

  double delta_vol;
  
  for (int i=0; i<nbins; i++) acc_gdir[i]=0.;
  double acc_pot=0.,acc_kin=0.,acc_etot=0.,acc_temp=0.;
 
  int nconf = 1;
  for(int iblock=1; iblock<=nblocks; ++iblock) {
    cout << "Block number: " << iblock << endl;
    for(int istep=1; istep <= nstep/nblocks; ++istep){
      Move();           //Move particles with Verlet algorithm
      if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
      if(istep%10 == 0){
        Measure();     //Properties measurement
        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
        acc_pot+=stima_pot;
        acc_kin+=stima_kin;
        acc_etot+=stima_etot;
        acc_temp+=stima_temp;
        for(int i=0; i<nbins; i++) {
          acc_gdir[i]+=gdir[i];
          gdir[i]=0.;
        }
      }
    }
    acc_pot=acc_pot/(nstep/nblocks)*10.;
    acc_kin=acc_kin/(nstep/nblocks)*10.;
    acc_etot=acc_etot/(nstep/nblocks)*10.;
    acc_temp=acc_temp/(nstep/nblocks)*10.;
    for(int i=0; i<nbins; i++) {
      delta_vol=4.*M_PI/3.*(pow((i+1.)*bin_size, 3.)-pow(i*bin_size, 3.));
      acc_gdir[i]=acc_gdir[i]/rho/m_part/delta_vol/(nstep/nblocks)*10.;
    }
    pot_block.push_back(acc_pot);
    pot2_block.push_back(acc_pot*acc_pot);
    kin_block.push_back(acc_kin);
    kin2_block.push_back(acc_kin*acc_kin);
    etot_block.push_back(acc_etot);
    etot2_block.push_back(acc_etot*acc_etot);
    temp_block.push_back(acc_temp);
    temp2_block.push_back(acc_temp*acc_temp);
    for(int i=0; i<nbins; i++) {
      gdir_prog[i]+=acc_gdir[i];
      gdir2_prog[i]+=acc_gdir[i]*acc_gdir[i];
      acc_gdir[i]=0.;
    }
  }

  Analyse();
  
  ofstream Gave;
  Gave.open("gdir_ave.out");

  for(int i=0; i<nbins; i++) {
    gdir_prog[i]=gdir_prog[i]/nblocks;
    gdir2_prog[i]=gdir2_prog[i]/nblocks;
    err_gdir[i]=sqrt(abs(gdir2_prog[i]-gdir_prog[i]*gdir_prog[i])/(double)(nblocks-1.));
    Gave << i*bin_size << " " << gdir_prog[i] << " " << err_gdir[i] << endl;
  }

  ConfFinal();         //Write final configuration to restart


  return 0;
}


void InputParameters(void){ //Prepare all stuff for the simulation
  ifstream ReadInput;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> nblocks;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

  nbins=100;
  bin_size=box/2./(double)nbins;

}

void InputStart(void) {
//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ifstream ReadConf;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//Prepare initial velocities
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i){
     vx[i] = rand()/double(RAND_MAX) - 0.5;
     vy[i] = rand()/double(RAND_MAX) - 0.5;
     vz[i] = rand()/double(RAND_MAX) - 0.5;

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }
   return;
}

void InputCorrect(void) {
//Read initial configuration
  cout << "Read initial configuration from file config.final " << endl << endl;
  cout << "Trying correcting the temperature" << endl << endl;
  double t=0.;
  double sf;
  ifstream ReadLastConf,ReadLastConfOld;
  ReadLastConf.open("config.final");
  ReadLastConfOld.open("config_old.final");
  for (int i=0; i<npart; ++i){
    ReadLastConf >> x[i] >> y[i] >> z[i];
    ReadLastConfOld >> xold[i] >> yold[i] >> zold[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
    xold[i] = xold[i] * box;
    yold[i] = yold[i] * box;
    zold[i] = zold[i] * box;
  }
  Move();
  for (int i=0; i<npart; ++i) {
    vx[i] = Pbc(x[i] - xold[i])/delta;
    vy[i] = Pbc(y[i] - yold[i])/delta;
    vz[i] = Pbc(z[i] - zold[i])/delta;
    t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
  }
  stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
  t=0.;
  cout << "La temperatura stimata in partenza è " << stima_temp << endl;
  cout << "La temperatura attesa è " << temp << endl;
  sf=sqrt(temp/stima_temp);
  cout << "quindi il fattore di riscalamento è  " << sf << endl;
  for (int i=0; i<npart; ++i) {
    vx[i]=vx[i]*sf;
    vy[i]=vy[i]*sf;
    vz[i]=vz[i]*sf;
    xold[i]=Pbc(x[i]-vx[i]*delta);
    yold[i]=Pbc(y[i]-vy[i]*delta);
    zold[i]=Pbc(z[i]-vz[i]*delta);
  }
  for (int i=0; i<npart; ++i) {
    t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
  }
  stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
  
  return;
}

void InputFinal(void) {
  //Read initial configuration
  cout << "Read initial configuration from file config.final " << endl << endl;
  ifstream ReadLastConf,ReadLastConfOld;
  ReadLastConf.open("config.final");
  ReadLastConfOld.open("config_old.final");
  for (int i=0; i<npart; ++i){
    ReadLastConf >> x[i] >> y[i] >> z[i];
    ReadLastConfOld >> xold[i] >> yold[i] >> zold[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
    xold[i] = xold[i] * box;
    yold[i] = yold[i] * box;
    zold[i] = zold[i] * box;
  }
  Move();
  for (int i=0; i<npart; ++i) {
    vx[i] = Pbc(x[i] - xold[i])/delta;
    vy[i] = Pbc(y[i] - yold[i])/delta;
    vz[i] = Pbc(z[i] - zold[i])/delta;
    xold[i]=Pbc(x[i]-vx[i]*delta);
    yold[i]=Pbc(y[i]-vy[i]*delta);
    zold[i]=Pbc(z[i]-vz[i]*delta);
  }
  return;
}

void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }

  return f;
}

void Measure(){ //Properties measurement
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;

  Epot.open("output_epot.out",ios::app);
  Ekin.open("output_ekin.out",ios::app);
  Temp.open("output_temp.out",ios::app);
  Etot.open("output_etot.out",ios::app);

  v = 0.0; //reset observables
  t = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] );
     dy = Pbc( yold[i] - yold[j] );
     dz = Pbc( zold[i] - zold[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     bin=int(dr/bin_size);
     gdir[bin]+=2;

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();

    return;
}

void Analyse(void) {
  ofstream ave_pot,ave_kin,ave_etot,ave_temp;
  ave_pot.open("ave_epot.out", ios::app);
  ave_kin.open("ave_ekin.out", ios::app);
  ave_etot.open("ave_etot.out", ios::app);
  ave_temp.open("ave_temp.out", ios::app);

  for(int iblock=0; iblock<nblocks; ++iblock) {
    pot_prog=accumulate(pot_block.begin(), pot_block.begin()+iblock+1, 0.)/(iblock+1);
    pot2_prog=accumulate(pot2_block.begin(), pot2_block.begin()+iblock+1, 0.)/(iblock+1);
    kin_prog=accumulate(kin_block.begin(), kin_block.begin()+iblock+1, 0.)/(iblock+1);
    kin2_prog=accumulate(kin2_block.begin(), kin2_block.begin()+iblock+1, 0.)/(iblock+1);
    etot_prog=accumulate(etot_block.begin(), etot_block.begin()+iblock+1, 0.)/(iblock+1);
    etot2_prog=accumulate(etot2_block.begin(), etot2_block.begin()+iblock+1, 0.)/(iblock+1);
    temp_prog=accumulate(temp_block.begin(), temp_block.begin()+iblock+1, 0.)/(iblock+1);
    temp2_prog=accumulate(temp2_block.begin(), temp2_block.begin()+iblock+1, 0.)/(iblock+1);

    if (iblock==0) {
      err_pot=0.;
      err_kin=0.;
      err_etot=0.;
      err_temp=0.;
    }
    else {
      err_pot=sqrt((pot2_prog-pow(pot_prog, 2.))/iblock);
      err_kin=sqrt((kin2_prog-pow(kin_prog, 2.))/iblock);
      err_etot=sqrt((etot2_prog-pow(etot_prog, 2.))/iblock);
      err_temp=sqrt((temp2_prog-pow(temp_prog, 2.))/iblock);
    }

    
    ave_pot << pot_prog << " " << err_pot << endl;
    ave_kin << kin_prog << " " << err_kin << endl;
    ave_etot << etot_prog << " " << err_etot << endl;
    ave_temp << temp_prog << " " << err_temp << endl;

  }

  ave_pot.close();
  ave_kin.close();
  ave_etot.close();
  ave_temp.close();
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;
  ofstream WriteConfOld;

  cout << "Print final configuration (time T) to file config.final " << endl;
  cout << "And print configuration (time T-dt) to file config.old " << endl << endl;
  WriteConf.open("config.final");
  WriteConfOld.open("config_old.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    WriteConfOld << xold[i]/box << "   " << yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf.close();
  WriteConfOld.close();
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
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
