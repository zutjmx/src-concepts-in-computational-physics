//  Realization of the CPoissonEq class, Section 11.2

#include "stdafx.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include "matrix2d.h"
#include "PoissonEq.h"

using namespace std;

inline double max(const double a, const double b) {
   return a >= b ? a : b;
}
inline double min(const double a, const double b) {
   return a <= b ? a : b;
}

/*****************************************************************************/
CPoissonEq::CPoissonEq(void)
/*****************************************************************************/
//  empty constructor, please use the property assign to set the
//  Poisson equation's parameters
{
}
/*****************************************************************************/
CPoissonEq::CPoissonEq(int n, int m, double dLx, double dLy)
/*****************************************************************************/
//  input parameters:
//  n ..... # of grid-points in x-direction
//  m ..... # of grid-points in y-direction
//  dLx ... length in x-direction
//  dLy ... length in y-direction
{
   m_nN = n;
   m_nM = m;
   m_nNh = m_nN/2;
   m_nMh = m_nM/2;
   m_nNstep = m_nN/10;
   m_nMstep = m_nM/10;
   m_dLx = dLx;
   m_dLy = dLy;
   m_dH = dLx/(n-1);  //  distance between grid-points, x-direction
   m_dK = dLy/(m-1);  //  distance between grid-points, y-direction
   m_dH2 = m_dH*m_dH;
   m_dK2 = m_dK*m_dK;
   m_dFact = 0.5/(m_dH2+m_dK2);
   m_mRho.assign(n,m,0.0);    //  \rho(x,y)
   m_mPhi.assign(n,m,0.0);    //  \phi(x,y)
   m_mPhit.assign(n,m,0.0);
}
/*****************************************************************************/
CPoissonEq::~CPoissonEq(void)
/*****************************************************************************/
{
}
/*****************************************************************************/
void CPoissonEq::assign(int n, int m, double dLx, double dLy)
/*****************************************************************************/
//  use this property to initialze the Poisson equations's parameters if an
//  empty constructor has been used to instantiate the CPoissonEq class

//  input parameters:
//  n ..... # of grid-points in x-direction
//  m ..... # of grid-points in y-direction
//  dLx ... length in x-direction
//  dLy ... length in y-direction
{
   m_nN = n;
   m_nM = m;
   m_nNh = m_nN/2;
   m_nMh = m_nM/2;
   m_nNstep = m_nN/10;
   m_nMstep = m_nM/10;
   m_dLx = dLx;
   m_dLy = dLy;
   m_dH = dLx/n;  //  distance between grid-points, x-direction
   m_dK = dLy/m;  //  distance between grid-points, y-direction
   m_dH2 = m_dH*m_dH;
   m_dK2 = m_dK*m_dK;
   m_dFact = 1.0/(2.0*(m_dH2+m_dK2));
   m_mRho.assign(n,m,0.0);    //  \rho(x,y)
   m_mPhi.assign(n,m,0.0);    //  \phi(x,y)
   m_mPhit.assign(n,m,0.0);
}
/*****************************************************************************/
void CPoissonEq::Rho(int sw)
/*****************************************************************************/
//  generate the charge density distribution \rho(x,y) according to
//  Eqs. (11.14)
{
   switch (sw) {
   case 1:                 // \rho(x,y), Eq. (11.14a)
      for (register int i=m_nNh-m_nNstep-1; i<m_nNh+m_nNstep; ++i)
         for (register int j=m_nMh-m_nMstep-1; j<m_nMh+m_nMstep; ++j)
            m_mRho[i][j] = 50.0;
      break;
   case 2:                 // \rho(x,y), Eq. (11.14b)
      for (register int i=m_nNh-m_nNstep-1; i<m_nNh+m_nNstep; ++i)      //  Omega_1,3
         for (register int j=m_nMh-m_nMstep; j<m_nMh; ++j)
            m_mRho[i][j] = 50.0;
      for (register int i=m_nNh-m_nNstep-1; i<m_nNh+m_nNstep; ++i)      //  Omega_2,4
         for (register int j=m_nMh; j<m_nMh+m_nMstep; ++j)
            m_mRho[i][j] = -50.0;
      break;
   case 3:                 // \rho(x,y), Eq. (11.14c)
      for (register int i=m_nNh-m_nNstep-1; i<m_nNh; ++i)      //  Omega_1
         for (register int j=m_nMh-m_nMstep; j<m_nMh; ++j)
            m_mRho[i][j] = 50.0;
      for (register int i=m_nNh; i<m_nNh+m_nNstep; ++i)      //  Omega_2
         for (register int j=m_nMh; j<m_nMh+m_nMstep; ++j)
            m_mRho[i][j] = 50.0;
      for (register int i=m_nNh-m_nNstep-1; i<m_nNh; ++i)      //  Omega_3
         for (register int j=m_nMh; j<m_nMh+m_nMstep; ++j)
            m_mRho[i][j] = -50.0;
      for (register int i=m_nNh; i<m_nNh+m_nNstep; ++i)      //  Omega_4
         for (register int j=m_nMh-m_nMstep; j<m_nMh; ++j)
            m_mRho[i][j] = -50.0;
   default:
      break;
   }
}
/*****************************************************************************/
void CPoissonEq::SolveEquation(double eta)
/*****************************************************************************/
//  solve Poisson equation according to Eq. (11.12)
//  input parameter:
//  eta .... accuracy
{
   double dPhijmo,dPhijpo,xx,dMax;

   while (1) {  
      for (register int i=1; i<m_nN-1; ++i)  //  Eq. (11.12)
         for (register int j=1; j<m_nM-1; ++j) {
            dPhijmo = m_mPhi[i][j-1];
            dPhijpo = m_mPhi[i][j+1];
            m_mPhi[i][j] = m_dFact*m_dH2*m_dK2*m_mRho[i][j];
            xx = m_dK2*(m_mPhi[i-1][j]+m_mPhi[i+1][j]);
            xx += m_dH2*(dPhijpo+dPhijmo);
            m_mPhi[i][j] += m_dFact*xx;
         }
      dMax = 0.0;   //  find max. error!
      for (register int i=0; i<m_nN; ++i)
         for (register int j=0; j<m_nM; ++j)
            dMax = max(dMax,fabs(m_mPhi[i][j]-m_mPhit[i][j]));
      if (fabs(dMax) < eta)  // accuracy reached?
         break;
      m_mPhit = m_mPhi;  // save results for next step
   }
}
/*****************************************************************************/
void CPoissonEq::SetDirichlet(double x0, double xL, double y0, double yL)
/*****************************************************************************/
//  set Dirichlet boundary conditions, Eq. (11.9)
//  input parameters:
//  x0 ... left hand x-direction boundary condition
//  xL ... right hand x-direction boundary condition
//  y0 ... left hand y-direction boundary condition
//  yL ... right hand y-direction boundary condition
{
   m_mPhi.SetCols(0,x0);
   m_mPhi.SetCols(m_nN-1,xL);
   m_mPhi.SetRows(0,y0);
   m_mPhi.SetRows(m_nM-1,yL);
}
