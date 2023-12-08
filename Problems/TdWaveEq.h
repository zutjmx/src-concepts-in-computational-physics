//  Definition of the CTdWaveEq class, wave equation, Section 11.4

#pragma once

class CTdWaveEq
{
public:
   CTdWaveEq(void);
   CTdWaveEq(int,int,double,double,double,double);
   ~CTdWaveEq(void);
   void assign(int,int,double,double,double,double);
   void SolveEquation(int);
   void InitBoundCond(void);

private:
   int m_nN,         //  number of grid-points
      m_nK;          //  max. number of time steps
   double m_dC,      //  speed of wave propagation
      m_dL,          //  length of the string
      m_dLambda,     //  CFL stability criterion, Eq. (11.31)
      m_dA,          //  amplitude
      m_dH,          //  space discretization
      m_dDeltaT,     //  time step \Delta t
      m_dPi;         //  number Pi in machine accuracy

public:   //  allow access to results
   matrix2d<double> m_mU;     //  result u(x,t)
   std::vector<double> m_vX;  //  grid-points
};

