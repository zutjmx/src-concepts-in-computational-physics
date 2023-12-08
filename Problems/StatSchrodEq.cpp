//  Realization of the CStatSchrodEq class, Chapter 10

#include "stdafx.h"
#include <vector>
#include "StatSchrodEq.h"

using namespace std;

/*****************************************************************************/
CStatSchrodEq::CStatSchrodEq(void)
/*****************************************************************************/
//  empty constructor, please use the property assign to set the
//  Schroedinger equation's parameters
{
   m_nMax = m_nMaxIter = 0;
   m_dPi = 4.0*atan(1.0);   // number Pi in machine accuracy
}
/*****************************************************************************/
CStatSchrodEq::CStatSchrodEq(int nMax, int nMaxIter)
/*****************************************************************************/
//  input parameters:
//  nMax ............ number of grid-points
//  nMaxIter ........ mx. # of iterations
{
   m_nMax = nMax;
   m_nMaxIter = nMaxIter;
   m_dPi = 4.0*atan(1.0);   // number Pi in machine accuracy
   m_dH = 1.0/nMax;         // step size
   m_dH2 = m_dH*m_dH;       // step size squared
   m_d5ov6 = 5.0/6.0;
   m_d1ov6 = 1.0/6.0;
   m_vVl.assign(m_nMax+1,0.0);    // store potential 
   m_vPhil.assign(m_nMax+1,0.0);  //  \phi_n(s_l)
   m_vPhil[1] = 2.0*m_dH;         //  Eq. (10.56)
   m_vSl.assign(m_nMax+1,0.0);
   m_vSl2.assign(m_nMax+1,0.0);
   for (register int i=1; i<=m_nMax; ++i) {
      m_vSl[i] = i*m_dH;        //  grid-points  s_l
      m_vSl2[i] = m_vSl[i]*m_vSl[i];
   }
}
/*****************************************************************************/
CStatSchrodEq::~CStatSchrodEq(void)
/*****************************************************************************/
{
}
/*****************************************************************************/
void CStatSchrodEq::assign(int nMax, int nMaxIter)
/*****************************************************************************/
//  use this property to initialze the Schroedinger equations's parameters if an
//  empty constructor has been used to instantiate the CStatSchrodEq class

//  input parameters:
//  nMax ............ number of grid-points
//  nMaxIter ........ mx. # of iterations
{
   m_nMax = nMax;
   m_nMaxIter = nMaxIter;
   m_dH = 1.0/nMax;         // step size
   m_dH2 = m_dH*m_dH;       // step size squared
   m_d5ov6 = 5.0/6.0;
   m_d1ov6 = 1.0/6.0;
   m_vVl.assign(m_nMax+1,0.0);     //  store potential
   m_vPhil.assign(m_nMax+1,0.0);   //  \phi_n(s_l)
   m_vPhil[1] = 2.0*m_dH;          //  Eq. (10.56)
   m_vSl.assign(m_nMax+1,0.0);
   m_vSl2.assign(m_nMax+1,0.0);
   for (register int i=1; i<=m_nMax; ++i) {
      m_vSl[i] = i*m_dH;    //  grid-points s_l
      m_vSl2[i] = m_vSl[i]*m_vSl[i];
   }
}
/*****************************************************************************/
double CStatSchrodEq::SolveEquation(double eps)
/*****************************************************************************/
//  calulate the \phi_n(s_{l+1}) following Eq. (10.55)
{
   for (register int l=2; l<=m_nMax; ++l) {
      m_vPhil[l] = (2.0*(1.0-m_d5ov6*m_dH2*(eps-m_vVl[l-1]))*m_vPhil[l-1]-
                    (1.0+m_d1ov6*m_dH2*(eps-m_vVl[l-2]))*m_vPhil[l-2])/
                   (1.0+m_d1ov6*m_dH2*(eps-m_vVl[l]));
   }
   return m_vPhil[m_nMax];
}
/*****************************************************************************/
bool CStatSchrodEq::Numerov(double& eps, double step, double eta)
/*****************************************************************************/
//  Execute the Numerov algorithm of page 129
//  input parameter:
//  eps .... energy \epsilon_a of step 1 on page 129
//  step ... energy \epsilon_b = \epsilon_a+step
//  accuracy
//
//  on return = true: successful
//  eps .... eigenvalue
//  on return = false: iteration failed.
{
   bool bFlag = false;
   register int k = 1;
   double dEpsa = eps,       //  step 1
      dEpsb = dEpsa+step,    //  step 1
      dEpsc,
      dPhiNa,
      dPhiNb,
      dPhiNc,
      dEta = eta;

   while (k <= m_nMaxIter) {
      dPhiNa = SolveEquation(dEpsa);  // step 2
      dPhiNb = SolveEquation(dEpsb);  // step 2
      if (dPhiNa*dPhiNb > 0.0) {      // step 3
         dEpsa = dEpsb;
         dEpsb = dEpsa+1.0;
      } else {                        // step 4
         dEpsc = 0.5*(dEpsa+dEpsb);
         dPhiNc = SolveEquation(dEpsc);
         if (dPhiNa*dPhiNc < 0.0)     // step 5
            dEpsb = dEpsc;
         else if (dPhiNb*dPhiNc < 0.0) { // step 6
            bFlag = true;
            dEpsa = dEpsc;
         }
      }
      ++k;
      if (fabs(dEpsa-dEpsb) < dEta) {  //  step 7
         Normalize();        //  normalize wave function
         eps = dEpsa;        //  return eigenvalue
         break;
      }
   }
   return bFlag;             //  return status
}
/*****************************************************************************/
void CStatSchrodEq::Normalize(void)
/*****************************************************************************/
//  normalize wave function, Eq. (10.58)
{
   double sum = 0.0;

   for (register int i=0; i<=m_nMax; ++i)
      sum += m_vPhil[i]*m_vPhil[i];
   sum *= m_dH;
   sum = 1.0/sqrt(sum);
   for (register int i=0; i<=m_nMax; ++i)
      m_vPhil[i] *= sum;
}
/*****************************************************************************/
double CStatSchrodEq::ExpValue(vector<double>& v)
/*****************************************************************************/
//  calculate expectation value <v> from grid values in vector v
{
   double sum = 0.0;
   for (register int i=0; i<(int)v.size(); ++i)
      sum += m_vPhil[i]*m_vPhil[i]*v[i];
   return sum*m_dH;
}
/*****************************************************************************/
void CStatSchrodEq::Potential1(void)
/*****************************************************************************/
//  define potential v_1 of Eq. (10.59)
{
   for (register int i=0; i<=m_nMax; ++i)
      m_vVl[i] = 50.0*cos(m_dPi*m_vSl[i]);
}
/*****************************************************************************/
void CStatSchrodEq::Potential2(void)
/*****************************************************************************/
//  define potential v_2 of Eq. (10.59)
{
   for (register int i=0; i<=m_nMax; ++i) {
      double x = m_vSl[i]-0.5;
      m_vVl[i] = 50.0*exp(-x*x/0.08);
   }
}
/*****************************************************************************/
void CStatSchrodEq::Potential3(void)
/*****************************************************************************/
//  define potential v_3 of Eq. (10.59)
{
   for (register int i=0; i<=m_nMax; ++i)
      m_vVl[i] = 50.0*m_vSl[i];
}
