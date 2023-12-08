//  Realization of the CTdHeatEq class, time-dependent heat equation, Section 11.3

#include "stdafx.h"
#include <vector>
#include "matrix2d.h"
#include "TdHeatEq.h"

using namespace std;

/*****************************************************************************/
CTdHeatEq::CTdHeatEq(void)
/*****************************************************************************/
//  empty constructor, please use the property assign to set the
//  heat equation's parameters
{
}
/*****************************************************************************/
CTdHeatEq::CTdHeatEq(int n, int tmax, int k, double l, double kappa)
/*****************************************************************************/
//  input parameters:
//  n ....... number of grid-points
//  tmax .... max. number of time steps
//  k ....... max. array size
//  l ....... length of the rod
//  kappa ... thermal diffusivity
{
   m_nN = n;
   m_nTmax = tmax;
   m_nK = k;
   m_dL = l;
   m_dKappa = kappa;
   m_dDeltaT = static_cast<double>(m_nTmax)/static_cast<double>(k-1);
       //  \Delta t, time step
   m_dKapDt = m_dKappa*m_dDeltaT;
   m_dH = m_dL/(n-1);    //  distance between two grid-points
   m_dH2 = 1.0/(m_dH*m_dH);
   m_mT.assign(n,k,0.0); //  array for temperature profile T(x,t)
   m_vdX.assign(n,0.0);  //  vector for grid-points
   for (register int i=0; i<n; ++i)
      m_vdX[i] = i*m_dH; //  load the grid-points
}
/*****************************************************************************/
CTdHeatEq::~CTdHeatEq(void)
/*****************************************************************************/
{
}
/*****************************************************************************/
void CTdHeatEq::assign(int n, int tmax, int k, double l, double kappa)
/*****************************************************************************/
//  use this property to initialze the heat equations's parameters if an
//  empty constructor has been used to instantiate the CTdHeatEq class

//  input parameters:
//  n ....... number of grid-points
//  tmax .... max. number of time steps
//  k ....... max. array size
//  l ....... length of the rod
//  kappa ... thermal diffusivity
{
   m_nN = n;
   m_nTmax = tmax;
   m_nK = k;
   m_dL = l;
   m_dKappa = kappa;
   m_dDeltaT = static_cast<double>(m_nTmax)/static_cast<double>(k-1);
       //  \Delta t, time step
   m_dKapDt = m_dKappa*m_dDeltaT;
   m_dH = m_dL/(n-1);    //  distance between two grid-points
   m_dH2 = 1.0/(m_dH*m_dH);
   m_mT.assign(n,k,0.0); //  array for temperature profile T(x,t)
   m_vdX.assign(n,0.0);  //  vector for grid-points
   for (register int i=0; i<n; ++i)
      m_vdX[i] = i*m_dH; //  load the grid-points
}
/*****************************************************************************/
void CTdHeatEq::InitBoundCond(void)
/*****************************************************************************/
//  set boundary and initial conditions, Eqs. (11.24) and (11.25)
{
   m_mT.SetRows(0,0.0);           //  boundary conditions (11.24)
   m_mT.SetRows(m_nN-1,2.0);
   for (register int i=1; i<m_nN-1; ++i)
      m_mT[1][i] = 0.0;           //  initial condition (11.25)
}
/*****************************************************************************/
void CTdHeatEq::ExplicitEuler(int k)
/*****************************************************************************/
{
   double x;

   InitBoundCond();
   for (register int i=1; i<=k; ++i)    //  solve Eq. (11.19)
      for (register int j=1; j<m_nN-1; ++j) {
         x = (m_mT[j-1][i-1]-2.0*m_mT[j][i-1]+m_mT[j+1][i-1])*m_dH2;
         m_mT[j][i] = m_mT[j][i-1]+m_dKapDt*x;
      }
}
/*****************************************************************************/
int CTdHeatEq::TriDiag()
/*****************************************************************************/
//  solve Eq. (11.22) with tri-diagonal matrix A
//  on return:
//  = 1 ... first element of the diagonal is zero
//  = 2 ... the tri-diagonal matrix is ill conditioned
{
   int j;
   double bet;

   int n=m_vSubd.size();
   vector<double> gam(n);
   if (m_vDiag[0] == 0.0) {
      return 1;          //   Error 1
   }
   m_vTemp[0] = m_vRhs[0]/(bet = m_vDiag[0]);
   for (j=1; j<n; ++j) {
      gam[j] = m_vSuperd[j-1]/bet;
      bet = m_vDiag[j]-m_vSubd[j]*gam[j];
      if (bet == 0.0) {
         return 2;    //  Error 2
      }
      m_vTemp[j] = (m_vRhs[j]-m_vSubd[j]*m_vTemp[j-1])/bet;
   }
   for (j=(n-2); j>=0; --j)
      m_vTemp[j] -= gam[j+1]*m_vTemp[j+1];
   return 0;
}
/*****************************************************************************/
int CTdHeatEq::ImplicitEuler(int k)
/*****************************************************************************/
//  Solve Eq. (11.22) for k time steps
{
   int ret=0,nN=m_nN-2;
   double dDiag=1.0+2.0*m_dKappa*m_dDeltaT*m_dH2,
      dOffDiag=-m_dKappa*m_dDeltaT*m_dH2;
   vector<double> vF(nN,0.0);

   InitBoundCond();
   m_vDiag.assign(nN,dDiag);  //  define the diagonal elements of matrix A
                              //  Eq. (11.23)
   m_vSuperd.assign(nN,dOffDiag);  //  define the super-diagonal elements of matrix A
   m_vSubd.assign(nN,dOffDiag);    //  define the sub-diagonal elements of matrix A
   m_vRhs.assign(nN,0.0);          //  right hand side of Eq. (11.22)
   m_vTemp.assign(nN,0.0);
   vF[nN-1] = m_dKappa*m_dDeltaT*m_dH2*m_mT[m_nN-1][0];
                   //  last element of vector F
   for (register int i=1; i<=k; ++i) {   //  run over all k time steps
      for (register int j=0; j<nN; ++j)  //  generate right hand side elements of
                                         //  Eq. (11.22)
         m_vRhs[j] = m_mT[j+1][i-1]+vF[j];
      ret = TriDiag();                   //  solve Eq. (11.22)
      for (register int j=1; j<m_nN-1; ++j)
         m_mT[j][i] = m_vTemp[j-1];      //  save temeprature profile
   }
   return ret;
}
