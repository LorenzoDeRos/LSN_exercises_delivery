#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

int main() {
    double temp;
    vector<double> epot;
    vector<double> pres;
    ifstream in_epot;
    in_epot.open("epot_allvalues.gas");
    ifstream in_pres;
    in_pres.open("pres_allvalues.gas");
    ofstream out_chi_epot;
    out_chi_epot.open("chi_epot.gas");
    ofstream out_chi_pres;
    out_chi_pres.open("chi_pres.gas");
    for(int i=0; i<500000; i++) {
        in_epot >> temp;
        epot.push_back(temp);
        in_pres >> temp;
        pres.push_back(temp);
    }
    int n=epot.size();
    double add1_epot;
    double add1_pres;
    double add2_epot;
    double add2_pres;
    double add3_epot;
    double add3_pres;
    double add4_epot=0.;
    double add4_pres=0.;
    double add5_epot=0.;
    double add5_pres=0.;
    double num_epot;
    double num_pres;
    double den_epot;
    double den_pres;
    for(int imeas=0; imeas<n; imeas++) {
        add4_epot+=epot[imeas]*epot[imeas];
        add4_pres+=pres[imeas]*pres[imeas];
        add5_epot+=epot[imeas];
        add5_pres+=pres[imeas];
    }
    add4_epot=add4_epot/(double)n;
    add4_pres=add4_pres/(double)n;
    add5_epot=add5_epot/(double)n;
    add5_pres=add5_pres/(double)n;
    for(int imeas=0; imeas<n; imeas++){
        if(imeas%1000==0) cout << imeas << " done" << endl;
        add1_epot=0.;
        add1_pres=0.;
        add2_epot=0.;
        add2_pres=0.;
        add3_epot=0.;
        add3_pres=0.;
        for(int jmeas=0; jmeas<n-imeas; jmeas++) {
            add1_epot+=epot[jmeas]*epot[imeas+jmeas];
            add1_pres+=pres[jmeas]*pres[imeas+jmeas];
            add2_epot+=epot[jmeas];
            add2_pres+=pres[jmeas];
            add3_epot+=epot[imeas+jmeas];
            add3_pres+=pres[imeas+jmeas];
        }
        num_epot=add1_epot/(double)(n-imeas)-add2_epot*add3_epot/(double)(n-imeas)/(double)(n-imeas);
        num_pres=add1_pres/(double)(n-imeas)-add2_pres*add3_pres/(double)(n-imeas)/(double)(n-imeas);
        den_epot=add4_epot-add5_epot*add5_epot;
        den_pres=add4_pres-add5_pres*add5_pres;
        out_chi_epot << num_epot/den_epot << endl;
        out_chi_pres << num_pres/den_pres << endl;
    }

}