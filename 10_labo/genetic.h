#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;

class city

{

    private:

        double m_x, m_y;                            // coordinate della città

        int m_label;                                // label identificativo della città

    public:

        city(double x, double y, int label) {       // creatore di una nuova città
            m_x=x;
            m_y=y;
            m_label=label;
        };

        double getx() {return m_x;};                // ritorna coordinata x

        double gety() {return m_y;};                // ritorna coordinata y

        int getlabel() {return m_label;};           // ritorna il label identificativo

        void setx(double x) {m_x=x;};               // modifica la coordinata x

        void sety(double y) {m_y=y;};               // modifica la coordinata y

        void setlabel(int label) {m_label=label;};  // modifica il label identificativo

//      misura la distanza tra questa città e un'altra città passata come parametro
        double dist(city c) {return sqrt(pow(m_x-c.getx(), 2.) + pow(m_y-c.gety(), 2.));};

};

class organism

{

    private:

        vector<city> m_org;                                             // un organismo è un vettore di città

    public:

        void newcity(double x, double y, int label) {                   // aggiunge una città all'organismo
            city temp(x, y, label);
            m_org.push_back(temp);
        };

        void newcity(city c) {                                          // aggiunge una città all'organismo
            m_org.push_back(c);
        }

        void setcity(int index, double x, double y, int label) {        // modifica una città all'organismo
            m_org[index].setx(x);
            m_org[index].sety(y);
            m_org[index].setlabel(label);
        };

        void setcity(int index, city c) {                               // modifica una città all'organismo
            m_org[index].setx(c.getx());
            m_org[index].sety(c.gety());
            m_org[index].setlabel(c.getlabel());
        }

        city getcity(int icity) {return m_org[icity];};                 // ritorna una città dell'organismo

        long unsigned int getncity() {return m_org.size();};            // ritorna il numero di città dell'organismo

//      mischia l'ordine delle città di un organismo, senza mai spostare la città di partenza
        void shuffle() {random_shuffle(m_org.begin()+1, m_org.end());};

//      ritorna la lunghezza dell'organismo ricordandosi di tornare alla partenza dopo l'ultima città
        double length() {
            double ris=0.;
            for(long unsigned int icity=0; icity<m_org.size()-1; icity++) {
                ris+=m_org[icity].dist(m_org[icity+1]);
            }
            ris+=m_org[0].dist(m_org[m_org.size()-1]);
            return ris;
        };

        void clear() {m_org.clear();};      // libera l'organismo da tutte le città

        void print() {                      // stampa label e coordinate di tutte le città dell'organismo
            for(long unsigned int icity=0; icity<m_org.size(); icity++) {
                cout << m_org[icity].getlabel() << " " << m_org[icity].getx() << " " << m_org[icity].gety() << endl;
            }
        };

        void fileprint(string s) {          // stampa su file label e coordinate di tutte le città dell'organismo
            ofstream out;
            out.open(s);
            for(long unsigned int icity=0; icity<m_org.size(); icity++) {
                out << m_org[icity].getlabel() << " " << m_org[icity].getx() << " " << m_org[icity].gety() << endl;
            }
        };

        void stringprint() {                // stampa i label di tutte le città dell'organismo (sequenza genica identificativa)
            cout << endl;
            for(long unsigned int icity=0; icity<m_org.size(); icity++) {
                cout << m_org[icity].getlabel() << " ";
            }
            cout << endl;
        };

//      controlla che l'organismo abbia la città 1 come partenza e che possieda tutte le città
        bool check() {
            bool ris=true;
            vector<int> labels;
            for(long unsigned int icity=0; icity<m_org.size(); icity++){
                labels.push_back(m_org[icity].getlabel());
            }
            ris=ris&&(m_org[0].getlabel()==1);
            for(long unsigned int icity=1; icity<=m_org.size(); icity++){
                ris=ris&&(find(labels.begin(), labels.end(), icity) != labels.end());
            }
            return ris;
        };

//      cerca una certa città all'interno della sequenza genica compresa tra due indici
        bool found(int index1, int index2, int label) {
            vector<int> labels;
            for(int icity=index1; icity<index2; icity++) {
                labels.push_back(m_org[icity].getlabel());
            }
            return (find(labels.begin(), labels.end(), label) != labels.end());
        };

//      scambia due città
        void mutation1(double rnd1, double rnd2) {
            int index1=int((m_org.size()-1.)*rnd1+1.);                      //indexes must be in [1, m_org.size()-1]
            int index2=int((m_org.size()-1.)*rnd2+1.);
            iter_swap(m_org.begin()+index1, m_org.begin()+index2);
        };

//      cicla una certa sequenza genica
        void mutation2(double rnd1, double rnd2, double rnd3) {
            int index=int((m_org.size()-3.)*rnd1+1.);                       //index must be in [1, m_org.size()-3]
            int dim=int((m_org.size()-index-2)*rnd2+3.);                    //dim must be in [3, m_org.size()-index]
            int nshift=int((dim-1.)*rnd3+1.);                               //nshift must be in [1, dim-1]
            for(int ishift=0; ishift<nshift; ishift++) {
                m_org.emplace(m_org.begin()+index, m_org[index+dim-1]);
                m_org.erase(m_org.begin()+index+dim);
            }
        };

//      specchia una sequenza di città
        void mutation3(double rnd1, double rnd2) {                          //index1 must be in [1, m_org.size()-3]
            int index1=int((m_org.size()-3.)*rnd1+1.);
            int index2=int((m_org.size()-index1-2.)*rnd2+index1+2.);
            reverse(m_org.begin()+index1, m_org.begin()+index2);            //index2 must be in [index1+2, m_org.size()-1]
        };

//      scambia sequenze di città
        void mutation4(double rnd1, double rnd2, double rnd3) {
            int dim=int((m_org.size()-3.)*rnd1)+2;                          //dim must be in [2, m_org.size()-2]
            int index1=int((m_org.size()-dim)*rnd2)+1;                      //indexes must be in [1, m_org.size()-dim]
            int index2=int((m_org.size()-dim)*rnd3)+1;
            for(int ichange=0; ichange<dim; ichange++) {
                iter_swap(m_org.begin()+index1+ichange, m_org.begin()+index2+ichange);
            }
        };

//      scambia sequenze di città in ordine inverso
        void mutation5(double rnd1, double rnd2, double rnd3) {
            int dim=int((m_org.size()-3.)*rnd1)+2;                          //dim must be in [2, m_org.size()-2]
            int index1=int((m_org.size()-dim)*rnd2)+1;                      //indexes must be in [1, m_org.size()-dim]
            int index2=int((m_org.size()-dim)*rnd3)+1;
            for(int ichange=0; ichange<dim; ichange++) {
                iter_swap(m_org.begin()+index1+ichange, m_org.begin()+index2+dim-1-ichange);
            }
        };

};

class population

{

    private:

        vector<organism> m_pop;     // una popolazione è un vettore di organismi
        
        organism m_son1, m_son2;    // organismi in cui vengono salvati temporaneamente i figli

    public:

        void neworganism(organism org) {m_pop.push_back(org);};     // aggiunge un organismo alla popolazione

        organism getorganism(int iorg) {return m_pop[iorg];};       // restituisce un organismo

        void setorganism(int iorg, organism org) {                  // sostituisce un organismo con quello passato come parametro
            m_pop[iorg].clear();
            for(long unsigned int icity=0; icity<org.getncity(); icity++) {
                m_pop[iorg].newcity(org.getcity(icity));
            }
        };

        organism getson1() {return m_son1;};                        // restituisce il primo figlio

        organism getson2() {return m_son2;};                        // restituisce il secondo figlio

        void clear() {m_pop.clear();};                              // libera tutti gli organismi della popolazione
        
        void clearsons() {                                          // libera gli organismi contenenti i figli
            m_son1.clear();
            m_son2.clear();
        };

        void shuffle() {                                            // mischia le città di ogni individuo
            for(long unsigned int iorg=0; iorg<m_pop.size(); iorg++) {
                m_pop[iorg].shuffle();
            }
        };

        void print() {                                              // stampa ogni organismo e la sua lunghezza
            for(long unsigned int iorg=0; iorg<m_pop.size(); iorg++) {
                cout << endl;
                cout << iorg+1 << " organism" << endl << endl;
                m_pop[iorg].print();
                cout << endl;
                cout << m_pop[iorg].length() << " length" << endl << endl;
            }
        };

        void print(int index) {                                     // stampa un organismo e la sua lunghezza
            m_pop[index].print();
            cout << m_pop[index].length() << endl;
        };

        void printlengths() {                                       // stampa la lunghezza di ogni organismo
            for(long unsigned int iorg=0; iorg<m_pop.size(); iorg++) {
                cout << endl;
                cout << iorg+1 << " organism " << m_pop[iorg].length() << endl;
            }
        };

        void stringprint(int iorg) {                                // stampa la sequenza genica di ogni organismo
            m_pop[iorg].stringprint();
        };

        bool check() {                                              // controlla la popolazione
            bool ris=true;
            for(long unsigned int iorg=0; iorg<m_pop.size(); iorg++) {
                ris=ris&&m_pop[iorg].check();
            }
            return ris;
        };

//      ordina la popolazione dall'organismo più corto a quello più lungo
        void popsort(){
            sort(m_pop.begin(),
            m_pop.end(),
            [](organism one, organism two){
                return one.length() < two.length();
            });
        };

        int selection(double rnd) {                                 // selezione un genitore
            return int(m_pop.size()*rnd*rnd);
        };

//      Fanno una mutazione ad un organismo
        void mutation1(int iorg, double rnd1, double rnd2) {
            m_pop[iorg].mutation1(rnd1, rnd2);
        };
        void mutation2(int iorg, double rnd1, double rnd2, double rnd3) {
            m_pop[iorg].mutation2(rnd1, rnd2, rnd3);
        };
        void mutation3(int iorg, double rnd1, double rnd2) {
            m_pop[iorg].mutation3(rnd1, rnd2);
        };
        void mutation4(int iorg, double rnd1, double rnd2, double rnd3) {
            m_pop[iorg].mutation4(rnd1, rnd2, rnd3);
        };
        void mutation5(int iorg, double rnd1, double rnd2, double rnd3) {
            m_pop[iorg].mutation5(rnd1, rnd2, rnd3);
        };

//      genera 2 figli a partire da due genitori
        void crossover(int iorg, int jorg, int index){                  //index must be in[2, ncity-2]
            for(int icity=0; icity<index; icity++) {                    // copia la prima parte della sequenza genica
                m_son1.newcity(m_pop[iorg].getcity(icity));
                m_son2.newcity(m_pop[jorg].getcity(icity));
            }
            for(long unsigned int icity=0; icity<m_pop[iorg].getncity(); icity++) {       // completa la seconda parte con le città mancanti
                if(!m_son2.found(0, index, (m_pop[iorg].getcity(icity)).getlabel())) {    // nell'ordine in cui compaiono nell'altro genitore
                    m_son2.newcity(m_pop[iorg].getcity(icity));
                }
                if(!m_son1.found(0, index, (m_pop[jorg].getcity(icity)).getlabel())) {
                    m_son1.newcity(m_pop[jorg].getcity(icity));
                }
            }
        };

};