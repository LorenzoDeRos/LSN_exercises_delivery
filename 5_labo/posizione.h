#ifndef __posizione_h__
#define __posizione_h__

class Posizione {
public:
  // costruttori
  Posizione();
  Posizione(double x, double y, double z);
  // distruttore
  ~Posizione();
  // metodi
  double getX() const; // coordinate cartesiane
  double getY() const;
  double getZ() const;
  double getR() const;
  double Distanza(const Posizione &) const;                                     // distanza da un altro punto
  void AvanzaReticolo(double passo, int direzione, int verso);                  //traslazione discreta lungo un reticolo cubico 3D da effettuare sulla posizione
  void AvanzaDirezione(double passo, double longitudine, double latitudine);    //traslazione discreta in una direzione qualunque da effettuare sulla posizione
  void AvanzaXYZ(double X, double Y, double Z);					                        //traslazione delle coordinate in ingresso
  Posizione& operator=(const Posizione&);
private:
  double m_x, m_y, m_z;
};

#endif // __posizione_h__
