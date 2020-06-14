/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//parameters, observables
const int m_props=4;
int n_props;
int iv,ik,it,ie;
const int gdir_size=1000;
double bin_size, nbins;
double gdir[gdir_size];
double stima_pot, stima_kin, stima_etot, stima_temp;
 
// averages
double acc,att;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

double acc_gdir[gdir_size];
double gdir_prog[gdir_size];
double gdir2_prog[gdir_size];
double pot_prog,pot2_prog,kin_prog,kin2_prog,etot_prog,etot2_prog,temp_prog,temp2_prog;
double err_gdir[gdir_size];
double err_pot,err_kin,err_etot,err_temp;

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, nblocks, iprint, seed;
double delta;

//functions
void InputParameters(void);
void InputStart(void);
void InputCorrect(void);
void InputFinal(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
void Analyse(void);
double Force(int, int);
double Pbc(double);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
