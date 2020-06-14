#include "posizione.h"
#include <cmath>

// costruttore di default
Posizione::Posizione() {
  m_x = 0.;
  m_y = 0.;
  m_z = 0.;
}

// costruttore a partire da una terna cartesiana
Posizione::Posizione(double x, double y, double z) {
  m_x = x;
  m_y = y;
  m_z = z;
}

//overloading operatore di assegnazione
Posizione& Posizione::operator=(const Posizione & P) {
  m_x=P.getX();
  m_y=P.getY();
  m_z=P.getZ();
  return *this;
}

// distruttore (puo` essere vuoto)
Posizione::~Posizione() {}

// Coordinate cartesiane
double Posizione::getX() const { return m_x; }

double Posizione::getY() const { return m_y; }

double Posizione::getZ() const { return m_z; }

// Coordinate polari

double Posizione::getR() const { return sqrt(m_x*m_x+m_y*m_y+m_z*m_z);}

// distanza da un altro punto

double Posizione::Distanza(const Posizione &b) const {
  return sqrt(pow(getX() - b.getX(), 2) + pow(getY() - b.getY(), 2) +
              pow(getZ() - b.getZ(), 2));
}

void Posizione::AvanzaXYZ(double X, double Y, double Z) {
  m_x+=X;
  m_y+=Y;
  m_z+=Z;
}

void Posizione::AvanzaReticolo(double passo, int direzione, int verso) {
  if (direzione==1) {
    m_x+=passo*verso;
  }
  if (direzione==2) {
    m_y+=passo*verso;
  }
  if (direzione==3) {
    m_z+=passo*verso;
  }
}

void Posizione::AvanzaDirezione(double passo, double longitudine, double latitudine) {
  m_x+=passo*sin(latitudine)*cos(longitudine);
  m_y+=passo*sin(latitudine)*sin(longitudine);
  m_z+=passo*cos(latitudine);
}