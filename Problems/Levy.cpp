//  Realization of the CLevy class, Section 17.4, Lévy flight

#include "stdafx.h"
#include <stdlib.h>
#include <cmath>
#include <vector>
#include "Levy.h"


/*****************************************************************************/
CLevy::CLevy(void)
/*****************************************************************************/
//  empty constructor, please use the property assign to define the
//  system's parameters
{
}
/*****************************************************************************/
CLevy::CLevy(int n, double a, double l)
/*****************************************************************************/
//  input parameters:
//  n ... number of time steps
//  a ... Lévy index \alpha, Eq. (17.78)
//  l ... minimal flight length, Eq. (17.78)
{
   m_nN = n;
   m_dAlpha = a;
   m_dInvAlpha = 1.0/m_dAlpha;
   m_dL = l;
   m_dPi = 4.0*atan(1.0);  //  number Pi in machine accuracy
   m_vdT.assign(m_nN,0.0);
   m_vdX.assign(m_nN,0.0);
}
/*****************************************************************************/
CLevy::~CLevy(void)
/*****************************************************************************/
{
}
/*****************************************************************************/
void CLevy::assign(int n, double a, double l)
/*****************************************************************************/
//  use this property to define the parameters of the Lévy flight if the empty
//  constructor has been used to instantiate the CLevy class.

//  input parameters:
//  n ... number of time steps
//  a ... Lévy index \alpha, Eq. (17.78)
//  l ... minimal flight length, Eq. (17.78)
{
   m_nN = n;
   m_dAlpha = a;
   m_dInvAlpha = 1.0/m_dAlpha;
   m_dL = l;
   m_dPi = 4.0*atan(1.0);  //  number Pi in machine accuracy
   m_vdT.assign(m_nN,0.0);
   m_vdX.assign(m_nN,0.0);
}
/*****************************************************************************/
void CLevy::RunLevy(void)
/*****************************************************************************/
//  execute the Lévy flight to generate the curves of Fig. 17.7
{
   int nSign;
   double dDeltaX,dDeltaT,r;

   for (register int i=1; i<m_nN; ++i) {
      nSign = RandInt() < 0.5 ? -1 : 1;
      r = RandInt();
      if (r > 0.9999)
         r = RandInt();
      dDeltaX = nSign*m_dL/pow(1.0-r,m_dInvAlpha);  //  step in x-direction
      r = RandInt();
      dDeltaT = -log(1.0-r);    //  next time step from exponential distribution
      m_vdX[i] = m_vdX[i-1]+dDeltaX;
      m_vdT[i] = m_vdT[i-1]+dDeltaT;
   }
}
/*****************************************************************************/
void CLevy::RunLevy2d(void)
/*****************************************************************************/
//  execute the 2D Lévy flight to generate the curves of Fig. 17.8
{
   int nSign;
   double dDeltaL,dPhi,r;

   for (register int i=1; i<m_nN; ++i) {
      nSign = RandInt() < 0.5 ? -1 : 1;
      r = RandInt();
      if (r > 0.9999)
         r = RandInt();
      dDeltaL = nSign*m_dL/pow(1.0-RandInt(),m_dInvAlpha);
      dPhi = RandInt()*2.0*m_dPi;
      m_vdX[i] = m_vdX[i-1]+dDeltaL*cos(dPhi);   //  x-coordinate
      m_vdT[i] = m_vdT[i-1]+dDeltaL*sin(dPhi);   //  y-coordinate
   }
}
/*****************************************************************************/
void CLevy::RunWiener2d(void)
/*****************************************************************************/
//  execute the 2D random walk to generate the curves of Fig. 17.8
{
   double dDeltaL,dPhi;

   m_vdX[0] = m_vdT[0] = 0.0;
   for (register int i=1; i<m_nN; ++i) {
      dDeltaL = RandInt()*sqrt(m_dL);
      dPhi = RandInt()*2.0*m_dPi;
      m_vdX[i] = m_vdX[i-1]+dDeltaL*cos(dPhi);   //  x-coordinate
      m_vdT[i] = m_vdT[i-1]+dDeltaL*sin(dPhi);   //  y-coordinate
   }
}
