//  Defintion of the CKeplerIntegrator class
//  Integrators for the Kepler Problem, Chapter 4 and Section 5.5

#pragma once

#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>

class CKeplerIntegrator
{
public:
   CKeplerIntegrator(void);
   CKeplerIntegrator(double,double,double);
   ~CKeplerIntegrator(void);
   void assign(double,double,double);
   double Hamiltonean(void);
   void ExplicitEuler(void);
   void SymplecticEuler1(void);
   void SymplecticEuler2(void);
   void ImplicitEuler(void);

private:
   void InitialValues(void);

private:
   double m_dP10,
      m_dP20,
      m_dQ10,
      m_dQ20,
      m_dE,
      m_dDeltaT,
      m_dMaxT;
   double m_dP1,  // generalized momentum, component 1
      m_dP2,      // generalized momentum, component 2
      m_dQ1,  // generalized space coordinate, component 1
      m_dQ2,  // generalized space coordinate, component 2
      m_dT,   // time instance t
      m_dH;   // energy H(t)

// access for results

public:
   int m_nMaxElem;  //  # of elements-1 in the vectors:
   std::vector<double> m_vdT,     // time instances
      m_vdP1,      // results for generalized momenta, component 1
      m_vdP2,      // results for generalized momenta, component 2
      m_vdQ1,  // results for generalized space coordinates, component 1
      m_vdQ2,  // results for generalized space coordinates, component 2
      m_vdH;   // energy H(t)
};

