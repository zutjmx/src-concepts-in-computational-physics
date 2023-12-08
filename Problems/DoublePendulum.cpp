// Realization of the CDoublePendulum class
// Chapter 6, Double Pendulum

#include "stdafx.h"
#include <cmath>
#include <vector>
#include "DoublePendulum.h"

using namespace std;

inline int sign(double a) {
   return a >= 0.0 ? 1 : -1;
}

/*****************************************************************************/
CDoublePendulum::CDoublePendulum(void)
/*****************************************************************************/
//  empty constructor, please use the property assign to set the
//  double pendulum's parameters
{
}
/*****************************************************************************/
CDoublePendulum::CDoublePendulum(int n, double dt, double m, double l, double g)
/*****************************************************************************/
//  input parameters:
//  n .... max. number of time steps
//  dt ... time step
//  m .... mass
//  l .... length of pendulum
//  g .... acceleration due to gravity
{
   m_nMaxSteps = n;
   m_dDeltaT = dt;
   m_dMass = m;
   m_dL = l;
   m_dG = g;
   m_dMl2 = 1.0/(m_dMass*m_dL*m_dL);
   m_dMgl = m_dMass*m_dL*m_dG;
   m_dT = 0.0;   //  start at time instance zero
}
/*****************************************************************************/
CDoublePendulum::~CDoublePendulum(void)
/*****************************************************************************/
{
}
/*****************************************************************************/
void CDoublePendulum::assign(int n, double dt, double m, double l, double g)
/*****************************************************************************/
//  use this property to initialze the double pendulum's parameters if an
//  empty constructor has been used to instantiate the CDoublePendulum class

//  input parameters:
//  n .... max. number of time steps
//  dt ... time step
//  m .... mass
//  l .... length of pendulum
//  g .... acceleration due to gravity
{
   m_nMaxSteps = n;
   m_dDeltaT = dt;
   m_dMass = m;
   m_dL = l;
   m_dG = g;
   m_dMl2 = 1.0/(m_dMass*m_dL*m_dL);
   m_dMgl = m_dMass*m_dL*m_dG;
   m_dT = 0.0;  // start at time instance zero
}
/*****************************************************************************/
double CDoublePendulum::DotPhi1(double phi1, double phi2, double p1, double p2)
/*****************************************************************************/
//  calculate derivative \dot{\varphi}_1 according to Eq. (6.22a)
{
   double sinp = sin(phi1-phi2);

   sinp *= sinp;
   sinp = 1.0+sinp;
   return m_dMl2*(p1-p2*cos(phi1-phi2))/sinp;
}
/*****************************************************************************/
double CDoublePendulum::DotPhi2(double phi1, double phi2, double p1, double p2)
/*****************************************************************************/
//  calculate derivative \dot{\varphi}_1 according to Eq. (6.22b)
{
   double sinp = sin(phi1-phi2);

   sinp *= sinp;
   sinp = 1.0+sinp;
   return m_dMl2*(2.0*p2-p1*cos(phi1-phi2))/sinp;
}
/*****************************************************************************/
double CDoublePendulum::DotP1(double phi1, double phi2, double p1, double p2)
/*****************************************************************************/
//  calculate derivative \dot{p}_1 according to Eq. (6.22c)
{
   double sinp = sin(phi1-phi2),
      cosp = cos(phi1-phi2),
      zw,den;

   den = 1.0/(1.0+sinp*sinp);   //  1.0/(1.0+\sin^2(\phi_1-\phi_2)
   zw = -p1*p2*sinp;
   zw += (p1*p1+2.0*p2*p2-2.0*p1*p2*cosp)*sinp*cosp*den;
   zw *= m_dMl2*den;
   zw -= 2.0*m_dMgl*sin(phi1);
   return zw;
}
/*****************************************************************************/
double CDoublePendulum::DotP2(double phi1, double phi2, double p1, double p2)
/*****************************************************************************/
//  calculate derivative \dot{\varphi}_1 according to Eq. (6.22d)
{
   double sinp = sin(phi1-phi2),
      cosp = cos(phi1-phi2),
      zw,den;

   den = 1.0/(1.0+sinp*sinp);   //  1.0/(1.0+\sin^2(\phi_1-\phi_2)
   zw = p1*p2*sinp;
   zw -= (p1*p1+2.0*p2*p2-2.0*p1*p2*cosp)*sinp*cosp*den;
   zw *= m_dMl2*den;
   zw -= m_dMgl*sin(phi2);
   return zw;
}
/*****************************************************************************/
void CDoublePendulum::RungeKutta(void)
/*****************************************************************************/
{
   double dPhi1=m_dPhi10,dPhi2=m_dPhi20,dP1=m_dP10,dP2=m_dP20,
      y11,y12,y13,y14,y21,y22,y23,y24,y31,y32,y33,y34,y41,y42,y43,y44,
      dPhi1n,dPhi2n,dP1n,dP2n,fy1[4],fy2[4],fy3[4];

   m_vdT.assign(m_nMaxSteps+2,0.0);   //  initialize result vectors
   m_vdPhi1.assign(m_nMaxSteps+2,0.0);
   m_vdPhi2.assign(m_nMaxSteps+2,0.0);
   m_vdP1.assign(m_nMaxSteps+2,0.0);
   m_vdP2.assign(m_nMaxSteps+2,0.0);
   m_vdT[0] = m_dT;
   m_vdPhi1[0] = m_dPhi10;
   m_vdPhi2[0] = m_dPhi20;
   m_vdP1[0] = m_dP10;
   m_vdP2[0] = m_dP20;

   for (register int j=1; j<=m_nMaxSteps; ++j) {
            //  Runge-Kutta scheme according to Eq. (6.26) for the
            //  vectors y, Eq. (6.24) and F(y), Eq. (6.25)
      y11 = dPhi1;   //  first equation of (6.26), Y_1 = y_n
      y12 = dPhi2;
      y13 = dP1;
      y14 = dP2;
                     //  second equation of (6.26), Y_2 = y_n+...
      y21 = dPhi1+m_dDeltaT*0.5*(fy1[0] = DotPhi1(y11,y12,y13,y14));
      y22 = dPhi2+m_dDeltaT*0.5*(fy1[1] = DotPhi2(y11,y12,y13,y14));
      y23 = dP1+m_dDeltaT*0.5*(fy1[2] = DotP1(y11,y12,y13,y14));
      y24 = dP2+m_dDeltaT*0.5*(fy1[3] = DotP2(y11,y12,y13,y14));
                     //  third equation of (6.26), Y_3 = y_n+...
      y31 = dPhi1+m_dDeltaT*0.5*(fy2[0] = DotPhi1(y21,y22,y23,y24));
      y32 = dPhi2+m_dDeltaT*0.5*(fy2[1] = DotPhi2(y21,y22,y23,y24));
      y33 = dP1+m_dDeltaT*0.5*(fy2[2] = DotP1(y21,y22,y23,y24));
      y34 = dP2+m_dDeltaT*0.5*(fy2[3] = DotP2(y21,y22,y23,y24));
                     //  fourth equation of (6.26), Y_4 = y_n+...
      y41 = dPhi1+m_dDeltaT*(fy3[0] = DotPhi1(y31,y32,y33,y34));
      y42 = dPhi2+m_dDeltaT*(fy3[1] = DotPhi2(y31,y32,y33,y34));
      y43 = dP1+m_dDeltaT*(fy3[2] = DotP1(y31,y32,y33,y34));
      y44 = dP2+m_dDeltaT*(fy3[3] = DotP2(y31,y32,y33,y34));
                     //  last equation of (6.26), y_{n+1} = y_n+...
      m_vdPhi1[j] = dPhi1n = dPhi1+m_dDeltaT*(fy1[0]+2.0*fy2[0]+2.0*fy3[0]+
         DotPhi1(y41,y42,y43,y44))/6.0;
      m_vdPhi2[j] = dPhi2n = dPhi2+m_dDeltaT*(fy1[1]+2.0*fy2[1]+2.0*fy3[1]+
         DotPhi2(y41,y42,y43,y44))/6.0;
      m_vdP1[j] = dP1n = dP1+m_dDeltaT*(fy1[2]+2.0*fy2[2]+2.0*fy3[2]+
         DotP1(y41,y42,y43,y44))/6.0;
      m_vdP2[j] = dP2n = dP2+m_dDeltaT*(fy1[3]+2.0*fy2[3]+2.0*fy3[3]+
         DotP2(y41,y42,y43,y44))/6.0;
      m_vdT[j] = m_dT += m_dDeltaT;
      dPhi1 = dPhi1n;   //  save results for next time step
      dPhi2 = dPhi2n;
      dP1 = dP1n;
      dP2 = dP2n;
   }
}
/*****************************************************************************/
void CDoublePendulum::Poincare(void)
/*****************************************************************************/
//  generate Poincare plot
{
   double dPhi1=m_dPhi10,dPhi2=m_dPhi20,dP1=m_dP10,dP2=m_dP20,
      y11,y12,y13,y14,y21,y22,y23,y24,y31,y32,y33,y34,y41,y42,y43,y44,
      dPhi1n,dPhi2n,dP1n,dP2n,fy1[4],fy2[4],fy3[4];

   for (register int j=1; j<=m_nMaxSteps; ++j) {
            //  Runge-Kutta scheme according to Eq. (6.26) for the
            //  vectors y, Eq. (6.24) and F(y), Eq. (6.25)
      y11 = dPhi1;
      y12 = dPhi2;
      y13 = dP1;
      y14 = dP2;
      y21 = dPhi1+m_dDeltaT*0.5*(fy1[0] = DotPhi1(y11,y12,y13,y14));
      y22 = dPhi2+m_dDeltaT*0.5*(fy1[1] = DotPhi2(y11,y12,y13,y14));
      y23 = dP1+m_dDeltaT*0.5*(fy1[2] = DotP1(y11,y12,y13,y14));
      y24 = dP2+m_dDeltaT*0.5*(fy1[3] = DotP2(y11,y12,y13,y14));
      y31 = dPhi1+m_dDeltaT*0.5*(fy2[0] = DotPhi1(y21,y22,y23,y24));
      y32 = dPhi2+m_dDeltaT*0.5*(fy2[1] = DotPhi2(y21,y22,y23,y24));
      y33 = dP1+m_dDeltaT*0.5*(fy2[2] = DotP1(y21,y22,y23,y24));
      y34 = dP2+m_dDeltaT*0.5*(fy2[3] = DotP2(y21,y22,y23,y24));
      y41 = dPhi1+m_dDeltaT*(fy3[0] = DotPhi1(y31,y32,y33,y34));
      y42 = dPhi2+m_dDeltaT*(fy3[1] = DotPhi2(y31,y32,y33,y34));
      y43 = dP1+m_dDeltaT*(fy3[2] = DotP1(y31,y32,y33,y34));
      y44 = dP2+m_dDeltaT*(fy3[3] = DotP2(y31,y32,y33,y34));
      dPhi1n = dPhi1+m_dDeltaT*(fy1[0]+2.0*fy2[0]+2.0*fy3[0]+
         DotPhi1(y41,y42,y43,y44))/6.0;
      dPhi2n = dPhi2+m_dDeltaT*(fy1[1]+2.0*fy2[1]+2.0*fy3[1]+
         DotPhi2(y41,y42,y43,y44))/6.0;
      dP1n = dP1+m_dDeltaT*(fy1[2]+2.0*fy2[2]+2.0*fy3[2]+
         DotP1(y41,y42,y43,y44))/6.0;
      dP2n = dP2+m_dDeltaT*(fy1[3]+2.0*fy2[3]+2.0*fy3[3]+
         DotP2(y41,y42,y43,y44))/6.0;

//       check for a sign change in angle \varphi_2! We do this
//       instead of looking for \varphi_2(t) = 0 to see when the
//       second pendulum crosses the vertical plane from the left
//       hand side. See page 84

      if (sign(dPhi2)!=sign(dPhi2n) && dP2n>0.0) {
         m_vdPhi1.push_back(dPhi1n);  // save the coordinates of
         m_vdP1.push_back(dP1n);      // the crossing
      }
      dPhi1 = dPhi1n;
      dPhi2 = dPhi2n;
      dP1 = dP1n;
      dP2 = dP2n;
   }
   m_nMaxSteps = m_vdP1.size();  //  number of Poincare dots
}
