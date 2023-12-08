//  Defintion of the CTdHeatEq class, time dependent heat equation, Section 11.3

#pragma once

class CTdHeatEq
{
public:
   CTdHeatEq(void);
   CTdHeatEq(int,int,int,double,double);
   ~CTdHeatEq(void);
   void assign(int n, int tmax, int k, double l, double kappa);
   void ExplicitEuler(int k);
   int ImplicitEuler(int k);

private:
   int TriDiag();
   void InitBoundCond(void);

private:
   int m_nN,         //  space discretization
      m_nTmax,       //  max. # of time steps
      m_nK;          //  array size
   double m_dL,      //  length of the rod
      m_dKappa,      //  thermal diffusivity
      m_dDeltaT,     //  time step
      m_dH,          //  space step h
      m_dH2,         //  h^2
      m_dKapDt;      //  kappa*deltat
   std::vector<double>m_vDiag,  //  diagonal of the A-Matrix
      m_vSuperd,                //  super-diagonal of the A-Matrix
      m_vSubd,                  //  sub-diagonal of the A-Matrix
      m_vGamma,                 //  Heat Sink/Source
      m_vRhs,                   //  right hand side
      m_vTemp;                  //  temperatures

public:   //  enable access to the results
   matrix2d<double> m_mT;    //  Temperature distribution T(x,t)
   std::vector<double> m_vdX;   //  space grid
};
