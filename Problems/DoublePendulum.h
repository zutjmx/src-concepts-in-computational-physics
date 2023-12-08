//  Definition of the CDoublePendulum class
//  Chapter 6

#pragma once

#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>

class CDoublePendulum
{
public:
   CDoublePendulum(void);
   CDoublePendulum(int,double,double,double,double);
   ~CDoublePendulum(void);
   void assign(int,double,double,double,double);
   void SetInitialValues(double phi1, double phi2, double p1, double p2) {
//  input parameters:
//  phi1 ... initial angle \varphi_1 of mass 1
//  phi2 ... initial angle \varphi_2 of mass 2
//  p1 ..... initial momentum of mass 1
//  p2 ..... initial momentum of mass 2
      m_dPhi10 = phi1;
      m_dPhi20 = phi2;
      m_dP10 = p1;
      m_dP20 = p2;
   }
   void RungeKutta(void);
   void Poincare(void);

private:
   double DotPhi1(double,double,double,double);
   double DotPhi2(double,double,double,double);
   double DotP1(double,double,double,double);
   double DotP2(double,double,double,double);

//  Variables

private:
   double m_dMass,         //  point mass
      m_dG,                //  gravitational acceleration
      m_dL,                //  angular momentum
      m_dMl2,              //  1.0/m\ell^2
      m_dMgl,              //  m*g*\ell
      m_dT,                //  time
      m_dDeltaT,           //  time step
      m_dPhi10,            //  initial values
      m_dPhi20,
      m_dP10,
      m_dP20;
   int m_nMaxSteps;        //  Max. # of steps

public:                    //  give acces to the results
   std::vector<double> m_vdT,  //  stored time instances
      m_vdPhi1,                //  stored angles \varphi_1
      m_vdPhi2,                //  stored angles \varphi_2
      m_vdP1,                  //  stored momenta p_1
      m_vdP2;                  //  stored momenta p_2
};

