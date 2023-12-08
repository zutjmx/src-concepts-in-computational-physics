//  ´Realization of the CRejectSamp class, Section 13.3 Rejection metod

#include "stdafx.h"
#include <stdlib.h>
#include <cmath>
#include <vector>
#include "RejectSamp.h"

using namespace std;

/*****************************************************************************/
CRejectSamp::CRejectSamp(void)
/*****************************************************************************/
//  empty constructor, please use the property assign to define the
//  function q(x), h(x), and ih(x) used to generate the random numbers, Eqs. (13.33)
//  (13.24), and (13.26) or (13.27)
{
   m_Qx = 0;
   m_Hx = 0;
}
/*****************************************************************************/
CRejectSamp::CRejectSamp(double l, double s, double c, double(*q)(double),
   double(*h)(double),double(*ih)(double))
/*****************************************************************************/
//  input parameters:
//  l .... parameter \lambda, Eq. (13.24)
//  s .... parameter \sigma, Eq. (13.33)
//  c .... c_min, Eq. (13.41)
//  q ..... pointer to the pdf q(x), Eq. (13.33)
//  h ..... pointer to the envelop function f(x), Eq. (13.24)
//  ih .... pointer to the inverse envelop function, Eq. (13.26) or (13.27)
{
   m_dLambda = l;
   m_dSigma = s;
   m_dCmin = c;
   m_Qx = q;
   m_Hx = h;
   m_InvHx = ih;
}
/*****************************************************************************/
CRejectSamp::~CRejectSamp(void)
/*****************************************************************************/
{
}
/*****************************************************************************/
void CRejectSamp::assign(double l, double s, double c, double(*q)(double),
   double(*h)(double),double(*ih)(double))
/*****************************************************************************/
//  use this property to initialze the to define the functions q(x), h(x), and
//  ih(x) if an empty constructor has been used to instantiate the CRejectSamp
//  class.

//  input parameters:
//  l .... parameter \lambda, Eq. (13.24)
//  s .... parameter \sigma, Eq. (13.33)
//  c .... c_min, Eq. (13.41)
//  q ..... pointer to the pdf q(x), Eq. (13.33)
//  h ..... pointer to the envelop function f(x), Eq. (13.24)
//  ih .... pointer to the inverse envelop function, Eq. (13.26) or (13.27)
{
   m_dLambda = l;
   m_dSigma = s;
   m_dCmin = c;
   m_Qx = q;
   m_Hx = h;
   m_InvHx = ih;
}
/*****************************************************************************/
void CRejectSamp::Histogram(int m, int n)
/*****************************************************************************/
//  generate histogram, Fig. 13.3
//  the random numbers are generated following the scheme on page 182
{
   double r,xt,dStep;
   register int i=0;

   m_vdX.assign(n,0.0);
   while (i < n) {
      r = RandInt();       //  step 1
      xt = (*m_InvHx)(r);  //  step 2
      r = RandInt();
      if (r < m_Qx(xt)/(m_dCmin*m_Hx(xt))) {  //  accept ? step 3
         r = RandInt();    //  step 4, get sign
         if (r < 0.5)
            m_vdX[i++] = -xt;  //  save random numbers
         else
            m_vdX[i++] = xt;
      }
   }
   xt = m_vdX[0];
   for (register int i=1; i<n; ++i)
      xt = xt >= m_vdX[i] ? xt : m_vdX[i];
   xt += 1.0e-7;
   dStep = 2.0*xt/m;
   m_vdRange.assign(m+1,0.0);
   m_vnCount.assign(m,0);

   //   generate histogram

   m_dHistScale = static_cast<double>(m)/static_cast<double>(n)/(2.0*xt);
   m_vdRange[0] = -dStep*0.5*m;
   for (register int i=1; i<=m; ++i)
      m_vdRange[i] = m_vdRange[i-1]+dStep;
   for (register int i=0; i<n; ++i) {
      for (register int j=0; j<m; ++j) {
         xt = m_vdX[i];
         if (xt >= m_vdRange[j] && xt < m_vdRange[j+1]) {
            ++m_vnCount[j];
            break;
         }
      }
   }
}
