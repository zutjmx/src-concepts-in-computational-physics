//  Definition of the CTdSchrodEq class, time-dependent Schroedinger equation
//  Section 11.5

#pragma once

#include <complex>

class CTdSchrodEq
{
public:
   CTdSchrodEq(void);
   CTdSchrodEq(int,int,double,double,double,double,double,double);
   ~CTdSchrodEq(void);
   void assign(int,int,double,double,double,double,double,double);
   void SetPotential(int);
   void InitPsi(void);
   void SolveEquation(int);

private:
   int m_nN,             //  number of grid points
      m_nM;              //  max. number of time steps
   double m_dL,          //  over all length
      m_dDeltaX,         //  space discretization
      m_dDeltaX2,        //  space discretization squared
      m_dQ0,             //  momentum q of the wave packet, Eq. (11.77)
      m_dX0,             //  center position of the wave packet
      m_dMass,           //  particle's mass
      m_dHbar,           //  Planck's constant
      m_dHbar2,          //  Planck's constant squared
      m_dSigma2,         //  \sigma^2, width of the wave packet squared
      m_dDeltaT;         //  time step
   std::complex<double> m_cZero,
      m_cI;
   std::vector<double> m_vdV;   //  vector to save potential
   std::vector<std::complex<double> > m_vcA,
      m_vcB,
      m_vcOmega;

public:  //  allow access to the results
   std::vector<double> m_vdX;    //  grid-points
   matrix2d<std::complex<double> > m_mcPsi;  //  \Psi(x,t)
};

