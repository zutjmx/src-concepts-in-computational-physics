//  Realization of the CTdWaveEq class, wave equation, Section 11.4

#include "stdafx.h"
#include <cmath>
#include <vector>
#include "matrix2d.h"
#include "TdWaveEq.h"

using namespace std;

/*****************************************************************************/
CTdWaveEq::CTdWaveEq(void)
/*****************************************************************************/
//  empty constructor, please use the property assign to set the
//  wave equation's parameters
{
}
/*****************************************************************************/
CTdWaveEq::CTdWaveEq(int n, int k, double c, double l, double lambda, double a)
/*****************************************************************************/
//  input parameters:
//  n .......... number of grid-points
//  k .......... max. number of time steps
//  c .......... speed of wave propagation
//  l .......... length of the string
//  lambda ..... CFL stability criterion, Eq. (11.31)
//  a .......... amplitude
{
   double dStep;

   m_nN = n;
   m_nK = k;
   m_dC = c;
   m_dL = l;
   m_dLambda = lambda;
   m_dLambda *= m_dLambda;
   m_dA = a;
   m_dH = m_dL/(m_nN-1);            //  space discretization
   m_dDeltaT = m_dLambda*m_dH/m_dC; //  calculate the time step, Eq. (11.31)
   m_dPi = 4.0*atan(1.0);           //  number Pi in machine accuracy
   m_mU.assign(m_nN,m_nK,0.0);      //  array for solution u(x,t)
   m_vX.assign(m_nN,0.0);           //  vector for grid-points
   dStep = m_dL/(m_nN-1);
   for (register int i=0; i<m_nN; ++i)
      m_vX[i] = i*dStep;            //  load gridpoints
}
/*****************************************************************************/
CTdWaveEq::~CTdWaveEq(void)
/*****************************************************************************/
{
}
/*****************************************************************************/
void CTdWaveEq::assign(int n, int k, double c, double l, double lambda, double a)
/*****************************************************************************/
//  use this property to initialze the heat equations's parameters if an
//  empty constructor has been used to instantiate the CTdHeatEq class

//  input parameters:
//  n .......... number of grid-points
//  k .......... max. number of time steps
//  c .......... speed of wave propagation
//  l .......... length of the string
//  lambda ..... CFL stability criterion, Eq. (11.31)
//  a .......... amplitude
{
   double dStep;

   m_nN = n;
   m_nK = k;
   m_dC = c;
   m_dL = l;
   m_dLambda = lambda;
   m_dLambda *= m_dLambda;
   m_dA = a;
   m_dH = m_dL/(m_nN-1);            //  space discretization
   m_dDeltaT = m_dLambda*m_dH/m_dC; //  calculate the time step, Eq. (11.31)
   m_dPi = 4.0*atan(1.0);           //  number Pi in machine accuracy
   m_mU.assign(m_nN,m_nK,0.0);      //  array for solution u(x,t)
   m_vX.assign(m_nN,0.0);           //  vector for grid-points
   dStep = m_dL/(m_nN-1);
   for (register int i=0; i<m_nN; ++i)
      m_vX[i] = i*dStep;            //  load gridpoints
}
/*****************************************************************************/
void CTdWaveEq::InitBoundCond(void)
/*****************************************************************************/
//  Load initial and boundary conditions.
{
   int nNo2=m_nN/2;

   //  Initial conditions, Eq. (11.39)

   for (register int i=nNo2; i<m_nN; ++i)
      m_mU[i][0] = m_dA*sin(2.0*m_dPi*m_vX[i]/m_dL);
   for (register int i=1; i<m_nN-1; ++i)  //  Eq. (11.37)
      m_mU[i][1] = (1.0-m_dLambda)*m_mU[i][0]+0.5*m_dLambda*(m_mU[i-1][0]+
         m_mU[i+1][0]);

   //  Boundary conditions, the string is fixed at both ends

   m_mU.SetRows(0,0.0);
   m_mU.SetRows(m_nN-1,0.0);
}
/*****************************************************************************/
void CTdWaveEq::SolveEquation(int k)
/*****************************************************************************/
//  solve Eq. (11.30) for k time steps
{
   for (register int i=1; i<k; ++i)   //  loop over all k time steps
      for (register int j=1; j<m_nN-1; ++j)  //  loop over the grid-points
         m_mU[j][i+1] = 2.0*(1.0-m_dLambda)*m_mU[j][i]-m_mU[j][i-1]+
            m_dLambda*(m_mU[j-1][i]+m_mU[j+1][i]);
}
