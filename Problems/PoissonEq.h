//  Definition of the CPoissonEq class, Section 11.2

#pragma once
class CPoissonEq
{
public:
   CPoissonEq(void);
   CPoissonEq(int,int,double,double);
   ~CPoissonEq(void);
   void assign(int,int,double,double);
   void Rho(int);
   void SolveEquation(double eta);
   void SetDirichlet(double x0, double xL, double y0, double yL);

private:
   int m_nN,        //  # of grid-points x-direction
      m_nM,         //  # of grid-points y-direction
      m_nNh,        //  m_nN/2
      m_nMh,        //  m_nM/2
      m_nNstep,     //  step size for \rho(x,y), x-direction
      m_nMstep;     //  step size for \rho(x,y), y-direction
   double m_dLx,    //  length in x-direction
      m_dLy,        //  length in y-direction
      m_dH,         //  distance between grid-points, x-direction
      m_dK,         //  distance between grid-points, y-direction
      m_dH2,        //  m_dH*m_dH
      m_dK2,        //  m_dK*m_dK
      m_dFact;      //  0.5*(m_dH2+m_dK2)
      matrix2d<double> m_mPhit;  //  intermediate storage

public:  //  open access to results 
   matrix2d<double> m_mPhi;  // \phi(x,y)
   matrix2d<double> m_mRho;  // \rho(x,y)
};

