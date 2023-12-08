//  Realization of the CKeplerIntegrator class
//  Integrators for the Kepler Problem, Chapter 4 and Section 5.5

#include "stdafx.h"
#include <cmath>
#include <vector>
#include "Integrator.h"

using namespace std;

/*****************************************************************************/
CKeplerIntegrator::CKeplerIntegrator(void) //  empty constructor
/*****************************************************************************/
//  please use the assign property to assign the relevant parameters
{
}
/*****************************************************************************/
CKeplerIntegrator::CKeplerIntegrator(double e, double dt, double mt)
/*****************************************************************************/
//  standard constructor
//  e .... excentricity
//  dt ... time step
//  mt ... maximum time
{
   m_dE = e;
   m_dDeltaT = dt;
   m_dMaxT = mt;
   InitialValues();
}
/*****************************************************************************/
CKeplerIntegrator::~CKeplerIntegrator(void)
/*****************************************************************************/
{
}
/*****************************************************************************/
void CKeplerIntegrator::assign(double e, double dt, double mt)
/*****************************************************************************/
//  required if you initialize the CKeplerIntegrator class with an empty constructor
//  e .... excentricity
//  dt ... time step
//  mt ... maximum time
{
   m_dE = e;
   m_dDeltaT = dt;
   m_dMaxT = mt;
   InitialValues();
}
/*****************************************************************************/
void CKeplerIntegrator::InitialValues()
/*****************************************************************************/
{
   m_nMaxElem = static_cast<int>(m_dMaxT/m_dDeltaT)+1; // maximum number of
                                                       // elements
   m_vdT.assign(m_nMaxElem,0.0);     //  store time steps
   m_vdP1.assign(m_nMaxElem,0.0);    //  store generalized momentum, component 1
   m_vdP2.assign(m_nMaxElem,0.0);    //  store generalized momentum, component 2
   m_vdQ1.assign(m_nMaxElem,0.0);    //  store generalized space coordinate,
                                     //  component 1
   m_vdQ2.assign(m_nMaxElem,0.0);    //  store generalized space coordinate,
                                     //  component 2
   m_vdH.assign(m_nMaxElem,0.0);     //  store total energy
   m_dP1 = m_dP10 = m_dQ2 = m_dQ20 = 0.0;  //  initial values Eq. (5.70) and
   m_dP2 = m_dP20 = sqrt((1.0+m_dE)/(1.0-m_dE));  //  Eq. (5.71)
   m_dQ1 = m_dQ10 = 1.0-m_dE;
   m_vdH[0] = m_dH = Hamiltonean();
   m_vdT[0] = m_dT = 0.0;            //  time initialization
   m_vdP1[0] = m_dP1;                //  store starting point
   m_vdP2[0] = m_dP2;
   m_vdQ1[0] = m_dQ1;
   m_vdQ2[0] = m_dQ2;
}
/*****************************************************************************/
double CKeplerIntegrator::Hamiltonean()
/*****************************************************************************/
//  calculate the Hamilton function, Eq. (5.63)
{
   return 0.5*(m_dP1*m_dP1+m_dP2*m_dP2)-1.0/sqrt(m_dQ1*m_dQ1+m_dQ2*m_dQ2);
}
/*****************************************************************************/
void CKeplerIntegrator::ExplicitEuler(void)
/*****************************************************************************/
//  Explicit Euler integrator according to Eqs. (5.65)
{
   double arg,sqarg,P11,P21,Q11,Q21;
   register int i;

   m_dT += m_dDeltaT;      //  define first time instance after initialization
   for (i=1; m_dT <= m_dMaxT; m_dT+=m_dDeltaT, ++i) {
                           //  loop over all time instances
      arg = m_dQ1*m_dQ1+m_dQ2*m_dQ2;
      sqarg = arg*sqrt(arg);       // denominator of Eqs. (5.65a) and (5.65b)
      P11 = m_dP1-m_dQ1*m_dDeltaT/sqarg; //  Eq. (5.65a)
      P21 = m_dP2-m_dQ2*m_dDeltaT/sqarg; //  Eq. (5.65b)
      Q11 = m_dQ1+m_dP1*m_dDeltaT;       //  Eq. (5.65c)
      Q21 = m_dQ2+m_dP2*m_dDeltaT;       //  Eq. (6.65d)
      m_vdP1[i] = m_dP1 = P11;           //  store results
      m_vdP2[i] = m_dP2 = P21;
      m_vdQ1[i] = m_dQ1 = Q11;
      m_vdQ2[i] = m_dQ2 = Q21;
      m_vdT[i] = m_dT;
      m_vdH[i] = m_dH = Hamiltonean();   //  calculate energy
   }
   m_nMaxElem = i-1;
}
/*****************************************************************************/
void CKeplerIntegrator::SymplecticEuler1(void)
/*****************************************************************************/
//  Symplectic Euler integrator according to Eqs. (5.68)
{
   double arg,sqarg,P11,P21,Q11,Q21;
   register int i;

   m_dT += m_dDeltaT;           //  define first time instance after initialization
   for (i=1; m_dT <= m_dMaxT; m_dT+=m_dDeltaT, ++i) {
                                //  loop over all time instances
      arg = m_dQ1*m_dQ1+m_dQ2*m_dQ2;
      sqarg = arg*sqrt(arg);    //  denominator of Eqs. (5.68a) and (5.68b)
      sqarg = m_dDeltaT/sqarg;
      P11 = m_dP1-m_dQ1*sqarg;  //  Eq. (5.68a)
      P21 = m_dP2-m_dQ2*sqarg;  //  Eq. (5.68b)
      Q11 = m_dQ1+m_dP1*m_dDeltaT-m_dQ1*m_dDeltaT*sqarg;
                                //  Eq. (5.68c)
      Q21 = m_dQ2+m_dP2*m_dDeltaT-m_dQ2*m_dDeltaT*sqarg;
                                //  Eq. (5.68d)
      m_vdP1[i] = m_dP1 = P11;  //  store results
      m_vdP2[i] = m_dP2 = P21;
      m_vdQ1[i] = m_dQ1 = Q11;
      m_vdQ2[i] = m_dQ2 = Q21;
      m_vdH[i] = m_dH = Hamiltonean(); //  calculate energy
      m_vdT[i] = m_dT;          //  store time instance
   }
   m_nMaxElem = i-1;
}
/*****************************************************************************/
void CKeplerIntegrator::SymplecticEuler2(void)
/*****************************************************************************/
//  Symplectic Euler integrator according to Eqs. (5.69)
{
   double arg,sqarg,P11,P21,Q11,Q21;
   register int i;

   m_dT += m_dDeltaT;           //  define first time instance after initialization
   for (i=1; m_dT <= m_dMaxT; m_dT+=m_dDeltaT, ++i) {
                                //  loop over all time instances
      Q11 = m_dQ1+m_dP1*m_dDeltaT;  // Eq. (5.69c)
      Q21 = m_dQ2+m_dP2*m_dDeltaT;  // Eq. (5.69d)
      arg = Q11*Q11+Q21*Q21;
      sqarg = arg*sqrt(arg);
      sqarg = m_dDeltaT/sqarg;  // denominator of Eqs. (5.69a) and (5.69b)
      P11 = m_dP1-Q11*sqarg;    // Eq. (5.69a)
      P21 = m_dP2-Q21*sqarg;    // Eq. (5.69b)
      m_vdP1[i] = m_dP1 = P11;  // store results
      m_vdP2[i] = m_dP2 = P21;
      m_vdQ1[i] = m_dQ1 = Q11;
      m_vdQ2[i] = m_dQ2 = Q21;
      m_vdH[i] = m_dH = Hamiltonean();  //  calculate energies
      m_vdT[i] = m_dT;          //  store time instance
   }
   m_nMaxElem = i-1;
}
/*****************************************************************************/
void CKeplerIntegrator::ImplicitEuler(void)
/*****************************************************************************/
//  Implicit Euler integrator according to Eqs. (5.66)
{
   double arg,sqarg,P11,P12,P21,P22,Q11,Q12,Q21,Q22;
   register int i,j;
   bool bFlag;

   m_dT += m_dDeltaT;      //  define first time instance after initialization
   P12 = P11 = m_dP1;
   P22 = P21 = m_dP2;
   Q12 = m_dQ1;
   Q22 = m_dQ2;
   for (i=1; m_dT <= m_dMaxT; m_dT+=m_dDeltaT, ++i) {
                           //  loop over all time instances
      for (j=1; j<=100; ++j) {
                           //  predictor - corrector loop
         Q11 = m_dQ1+P11*m_dDeltaT;    //  predictor, Eq. (5.66c)
         Q21 = m_dQ2+P21*m_dDeltaT;    //  Eq. (5.66d)
         arg = Q11*Q11+Q21*Q21;
         sqarg = arg*sqrt(arg);   //  denominator of Eqs. (5.66a) and (5.66b)
         sqarg = m_dDeltaT/sqarg;
         P11 = m_dP1-Q11*sqarg;        //  corrector, Eq. (5.66a)
         P21 = m_dP2-Q21*sqarg;        //  Eq. (5.66b)
         bFlag = fabs(P11-P12)>1.0e-9 || fabs(P21-P22)>1.0e-9 ||
                 fabs(Q11-Q12)>1.0e-9 || fabs(Q21-Q22)>1.0e-9;
                                       //  check accuracy
         if (bFlag) {                  //  accuracy not reached
            P12 = P11;
            P22 = P21;
            Q12 = Q11;
            Q22 = Q21;
         } else                        //  everyting is fine now
            break;                     //  leave the loop
      }
      if (!bFlag) {
         m_vdP1[i] = m_dP1 = P11;      //  store the results
         m_vdP2[i] = m_dP2 = P21;
         m_vdQ1[i] = m_dQ1 = Q11;
         m_vdQ2[i] = m_dQ2 = Q21;
         m_vdH[i] = m_dH = Hamiltonean();  //  calculate energy
         m_vdT[i] = m_dT;               //  store time instance
      } else {
         cout << "No convergence" << endl; // error exit
         break;
      }
   }
   m_nMaxElem = i-1;
}
