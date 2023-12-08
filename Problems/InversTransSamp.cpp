//  Realization of the CInversTransSamp class, Section 13.2

#include "stdafx.h"
#include <stdlib.h>
#include <cmath>
#include <vector>
#include "InversTransSamp.h"

using namespace std;

/*****************************************************************************/
CInversTransSamp::CInversTransSamp(void)
/*****************************************************************************/
//  empty constructor, please use the property assign to define the
//  function f(x) used to generate the random numbers, Eq. (13.26) or (13.27)
{
   m_Fx = 0;
}
/*****************************************************************************/
CInversTransSamp::CInversTransSamp(double l, double(*f)(double))
/*****************************************************************************/
//  input parameters:
//  l ..... parameter \lambda, Eq. (13.24)
//  f ..... pointer to the function f(x), Eq. (13.26) or (13.27)
{
   m_Fx = f;
   m_dLambda = l;
}
/*****************************************************************************/
CInversTransSamp::~CInversTransSamp(void)
/*****************************************************************************/
{
}
/*****************************************************************************/
void CInversTransSamp::assign(double l, double(*f)(double))
/*****************************************************************************/
//  use this property to define the function f(x) if an empty
//  constructor has been used to instantiate the CInversTransSamp class.

//  input parameters:
//  l ..... parameter \lambda, Eq. (13.24)
//  f ..... pointer to the function f(x), Eq. (13.26) or (13.27)
{
   m_Fx = f;
   m_dLambda = l;
}
/*****************************************************************************/
void CInversTransSamp::Histogram(int m, int n)
/*****************************************************************************/
//  generate the histograms, Fig. 13.1
{
   double xx,dStep;

   m_vdX.assign(n,0.0);
   for (register int i=0; i<n; ++i) {
      xx = RandInt();  //  get random number x \in [0,1)
      m_vdX[i] = (*m_Fx)(xx);  //  Eq. (13.26) or (13.27)
   }
   xx = m_vdX[0];
   for (register int i=1; i<n; ++i)
      xx = xx >= m_vdX[i] ? xx : m_vdX[i];
   xx += 1.0e-7;           //  determine range ...
   dStep = xx/m;
   m_vdRange.assign(m+1,0.0);
   m_vnCount.assign(m,0);
   m_dHistScale = static_cast<double>(m)/static_cast<double>(n)/xx;
   for (register int i=1; i<=m; ++i)
      m_vdRange[i] = i*dStep;
   for (register int i=0; i<n; ++i) { //  generate histogram
      for (register int j=0; j<m; ++j) {
         xx = m_vdX[i];
         if (xx >= m_vdRange[j] && xx < m_vdRange[j+1]) {
            ++m_vnCount[j];
            break;
         }
      }
   }
}
