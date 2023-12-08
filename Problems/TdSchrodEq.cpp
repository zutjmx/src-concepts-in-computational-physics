//  Realization of the CTdSchrodEq class, time-dependent Schroedinger equation
//  Section 11.5

#include "stdafx.h"
#include <cmath>
#include <complex>
#include <vector>
#include "matrix2d.h"
#include "TdSchrodEq.h"

using namespace std;

/*****************************************************************************/
CTdSchrodEq::CTdSchrodEq(void)
/*****************************************************************************/
//  empty constructor, please use the property assign to set the
//  Schroedinger equation's parameters
{
   m_cZero = complex<double>(0.0,0.0);
   m_cI = complex<double>(0.0,1.0);
}
/*****************************************************************************/
CTdSchrodEq::CTdSchrodEq(int n, int m, double q, double x, double ma, double h,
   double s, double dt)
/*****************************************************************************/
//  input parameters:
//  n .... number of grid points
//  m .... max. number of time steps
//  q .... momentum q of the wave packet, Eq. (11.77)
//  x .... center position of the wave packet
//  ma ... particle's mass
//  h .... Planck's constant
//  s .... \sigma^2, width of the wave packet squared
//  dt ... time step
{
   double dH;

   m_nN = n;
   m_nM = m;
   m_dQ0 = q;
   m_dX0 = x;
   m_dMass = ma;
   m_dHbar = h;
   m_dHbar2 = h*h;
   m_dSigma2 = s;
   m_dDeltaT = dt;
   m_cZero = complex<double>(0.0,0.0); 
   m_cI = complex<double>(0.0,1.0);
   m_mcPsi.assign(m_nM,m_nN,m_cZero);  // array \Psi(x,t)
   m_dL = m_nN-1.0;                    // length
   m_dDeltaX = m_dL/(m_nN-1);          // space discretization
   m_dDeltaX2 = m_dDeltaX*m_dDeltaX;
   dH = m_dL/(m_nN-1);
   m_vdX.assign(m_nN,0.0);             // vector for grid points
   for (register int i=1; i<m_nN; ++i)
      m_vdX[i] = i*dH;                 // set grid points
}
/*****************************************************************************/
CTdSchrodEq::~CTdSchrodEq(void)
/*****************************************************************************/
{
}
/*****************************************************************************/
void CTdSchrodEq::assign(int n, int m, double q, double x, double ma, double h,
                         double s, double dt)
/*****************************************************************************/
//  use this property to initialze the heat equations's parameters if an
//  empty constructor has been used to instantiate the CTdSchrodEq class

//  input parameters:
//  n .... number of grid points
//  m .... max. number of time steps
//  q .... momentum q of the wave packet, Eq. (11.77)
//  x .... center position of the wave packet
//  ma ... particle's mass
//  h .... Planck's constant
//  s .... \sigma^2, width of the wave packet squared
//  dt ... time step
{
   double dH;

   m_nN = n;
   m_nM = m;
   m_dQ0 = q;
   m_dX0 = x;
   m_dMass = ma;
   m_dHbar = h;
   m_dHbar2 = h*h;
   m_dSigma2 = s;
   m_dDeltaT = dt;
   m_cZero = complex<double>(0.0,0.0);
   m_cI = complex<double>(0.0,1.0);
   m_mcPsi.assign(m_nM,m_nN,m_cZero);  // array \Psi(x,t)
   m_dL = m_nN-1.0;                    // length
   m_dDeltaX = m_dL/(m_nN-1);          // space discretization
   m_dDeltaX2 = m_dDeltaX*m_dDeltaX;
   dH = m_dL/(m_nN-1);
   m_vdX.assign(m_nN,0.0);             // vector for grid points
   for (register int i=1; i<m_nN; ++i)
      m_vdX[i] = i*dH;                 // set grid points
}
/*****************************************************************************/
void CTdSchrodEq::SetPotential(int k)
/*****************************************************************************/
{
   m_vdV.assign(m_nN,0.0);
   switch (k) {
   case 1:        //  potential V_1(x), Eq. (11.78)
      for (register int i=249; i<260; ++i)
         m_vdV[i] = 0.7;
      break;
   case 2:        //  potential V_2(x), Eq. (11.79)
      for (register int i=249; i<260; ++i)
         m_vdV[i] = 0.7;
      for (register int i=299; i<310; ++i)
         m_vdV[i] = 0.7;
      break;
   default:
      break;
   }
}
/*****************************************************************************/
void CTdSchrodEq::InitPsi(void)
/*****************************************************************************/
//  Initialize the wave function
{
   double x,y,re,im;
   complex<double> cX;

   //  initial condition Eq. (11.77), wave packet, step 1 on page 153
   for (register int i=0; i<m_nN; ++i) {
      x = m_vdX[i]-m_dX0;
      x *= x;
      m_mcPsi[0][i] = exp(m_cI*m_dQ0*m_vdX[i])*exp(-x/m_dSigma2);
   }
   //  boundary condition of step 1 in page 153
   m_mcPsi[0][m_nN-1] = m_mcPsi[0][0] = m_cZero;
   x = 0;
   for (register int i=0; i<m_nN; ++i) {
      y = abs(m_mcPsi[0][i]);
      x += y*y;
   }
   x *= m_dDeltaX;
   x = 1.0/sqrt(x);
   for (register int i=0; i<m_nN; ++i)
      m_mcPsi[0][i] *= x;
   m_vcA.assign(m_nN,m_cZero);  //  vector a_k
   m_vcB.assign(m_nN,m_cZero);  //  vector b_k
   m_vcOmega.assign(m_nN,m_cZero);  //  vector \Omega_k
   cX = m_vcA[0];
   re = 2.0*(1.0+m_dMass*m_dDeltaX2*m_vdV[1]/m_dHbar2);
   im = 4.0*m_dMass*m_dDeltaX2/(m_dDeltaT*m_dHbar);
   m_vcA[1] = complex<double>(re,-im);  // Eq. (11.71), step 2
   cX = m_vcA[1];
   for (register int i=2; i<m_nN-1; ++i) {  //  Eq. (11.72), step 2
      re = 2.0*(1.0+m_dMass*m_dDeltaX2*m_vdV[i]/m_dHbar2);
      m_vcA[i] = complex<double>(re,-im)-1.0/m_vcA[i-1];
   }
}
/*****************************************************************************/
void CTdSchrodEq::SolveEquation(int k)
/*****************************************************************************/
//  solve time-dependent Schroedinger equation following the step outlined
//  on page 153. k time steps are executed
{
   double x,y,re,im;

   im = 4.0*m_dMass*m_dDeltaX2/(m_dDeltaT*m_dHbar);
   for (register int n=1; n<=k; ++n) {  //  loop over all time steps, step 3
      for (register int i=1; i<m_nN-1; ++i) {  //  step 4
         re = 2.0*(1.0+m_dMass*m_dDeltaX2*m_vdV[i]/m_dHbar2);
         m_vcOmega[i] = -m_mcPsi[n-1][i-1]+complex<double>(re,im)*m_mcPsi[n-1][i]-
                        m_mcPsi[n-1][i+1]; //  Eq. (11.73)
      }
      m_vcB[1] = m_vcOmega[1];    // step 5, Eq. (11.74)
      for (register int i=2; i<m_nN-1; ++i)  // step 5, Eq. (11.75)
         m_vcB[i] = m_vcB[i-1]/m_vcA[i-1]+m_vcOmega[i];
      m_mcPsi[n][m_nN-1] = m_cZero;
      for (register int i=m_nN-2; i>=1; --i)  // step 6, Eq. (11.76)
         m_mcPsi[n][i] = (m_mcPsi[n][i+1]-m_vcB[i])/m_vcA[i];
      m_mcPsi[n][0] = m_cZero;
      x = 0.0;
      for (register int i=0; i<m_nN; ++i) {  //  normalize wave function
         y = abs(m_mcPsi[n][i]);
         x += y*y;
      }
      x *= m_dDeltaX;
      x = 1.0/sqrt(x);
      for (register int i=0; i<m_nN; ++i)
         m_mcPsi[n][i] *= x;
   }
}
