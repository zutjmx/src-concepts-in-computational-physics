//  Definition of the CFractLevy class, Section 17.4, fractal time Lévy flight

#include "stdafx.h"
#include <stdlib.h>
#include <cmath>
#include <vector>
#include "FractLevy.h"

using namespace std;

/*****************************************************************************/
CFractLevy::CFractLevy(void)
/*****************************************************************************/
//  empty constructor, please use the property assign to define the
//  system's parameters
{
}
/*****************************************************************************/
CFractLevy::CFractLevy(int n, double a, double l, double b, double t)
/*****************************************************************************/
//  input parameters:
//  n ... number of time steps
//  a ... Lévy index \alpha, Eq. (17.78)
//  l ... minimal flight length, Eq. (17.78)
//  b ... parameter \beta, Eq. (17.81)
//  t ... minimal waighting time, Eq. (17.81)
{
   m_nN = n;
   m_dAlpha = a;
   m_dInvAlpha = 1.0/m_dAlpha;
   m_dL = l;
   m_dBeta = b;
   m_dInvBeta = 1.0/m_dBeta;
   m_dTau = t;
   m_vdT.assign(m_nN,0.0);
   m_vdX.assign(m_nN,0.0);
}
/*****************************************************************************/
CFractLevy::~CFractLevy(void)
/*****************************************************************************/
{
}
/*****************************************************************************/
void CFractLevy::assign(int n, double a, double l, double b, double t)
/*****************************************************************************/
//  use this property to define the parameters of the Lévy flight if the empty
//  constructor has been used to instantiate the CFractLevy class.

//  input parameters:
//  n ... number of time steps
//  a ... Lévy index \alpha, Eq. (17.78)
//  l ... minimal flight length, Eq. (17.78)
//  b ... parameter \beta, Eq. (17.81)
//  t ... minimal waighting time, Eq. (17.81)
{
   m_nN = n;
   m_dAlpha = a;
   m_dInvAlpha = 1.0/m_dAlpha;
   m_dL = l;
   m_dBeta = b;
   m_dInvBeta = 1.0/m_dBeta;
   m_dTau = t;
   m_vdT.assign(m_nN,0.0);
   m_vdX.assign(m_nN,0.0);
}
/*****************************************************************************/
void CFractLevy::RunFractLevy(void)
/*****************************************************************************/
//  generate the curves of Fig. 17.10
{
   int nSign;
   double dDeltaX,dDeltaT,r;

   for (register int i=1; i<m_nN; ++i) {
      nSign = RandInt() < 0.5 ? -1 : 1;
      r = RandInt();
      if (r > 0.9999)
         r = RandInt();
      dDeltaX = nSign*m_dL/pow(1.0-r,m_dInvAlpha);  //  new x-step
      r = RandInt();
      if (r > 0.9999)
         r = RandInt();
      dDeltaT = m_dTau/pow(1.0-r,m_dInvBeta);       //  new time step
      m_vdX[i] = m_vdX[i-1]+dDeltaX;                //  x(t)
      m_vdT[i] = m_vdT[i-1]+dDeltaT;                //  t
   }
}
