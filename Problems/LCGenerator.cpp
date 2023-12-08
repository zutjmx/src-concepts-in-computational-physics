// Realization of the CLCGenerator class, linear congruential generator
// Section 12.2

#include "stdafx.h"
#include <cmath>
#include <vector>
#include "LCGenerator.h"

using namespace std;

/*****************************************************************************/
CLCGenerator::CLCGenerator(void)
/*****************************************************************************/
//  empty constructor. set the Park-Miller parameters Eq. (12.10)
{
   m_nA = static_cast<unsigned long long>(pow(7.0,5.0));
   m_nC = 0;
   m_nM = (static_cast<unsigned long long>(2<<30))-1;
   m_dM = static_cast<double>(m_nM);
   m_nX0 = 281;  //  seed
}
/*****************************************************************************/
CLCGenerator::CLCGenerator(unsigned long long a, unsigned long long c,
   unsigned long long m, unsigned long long x0)
/*****************************************************************************/
//  input parameters:
//  a .... Park-Miller parameter a
//  c .... Park-Miller parameter c
//  m .... Park-Miller parameter m
//  x0 ... seed
{
   m_nA = a;
   m_nC = c;
   m_dM = static_cast<double>(m_nM = m);
   m_nX0 = x0;
}
/*****************************************************************************/
CLCGenerator::~CLCGenerator(void)
/*****************************************************************************/
{
}
/*****************************************************************************/
void CLCGenerator::assign(unsigned long long a, unsigned long long c,
   unsigned long long m, unsigned long long x0)
/*****************************************************************************/
//  use this property to change the Park-Miller parameters

//  input parameters:
//  a .... Park-Miller parameter a
//  c .... Park-Miller parameter c
//  m .... Park-Miller parameter m
//  x0 ... seed
{
   m_nA = a;
   m_nC = c;
   m_dM = static_cast<double>(m_nM = m);
   m_nX0 = x0;
}
/*****************************************************************************/
double CLCGenerator::rand()
/*****************************************************************************/
{
   m_nX0 = (m_nA*m_nX0+m_nC)%m_nM;   //  Eq. (12.9)
   return static_cast<double>(m_nX0)/(m_dM+1.0);  //  return randum number
                                                  //  r \in [0,1)
}
/*****************************************************************************/
void CLCGenerator::SpectralTest(int n)
/*****************************************************************************/
//  Perform spectral test for 2n random numbers, Fig. 12.1
{
   m_vdX1.assign(n,0.0);
   m_vdX2.assign(n,0.0);
   for (register int i=0; i<n; ++i) {
      m_vdX1[i] = rand();
      m_vdX2[i] = rand();
   }
}
/*****************************************************************************/
void CLCGenerator::Histogram(int m, int n)
/*****************************************************************************/
//  Perform histogram test for m bins and n random numbers, Fig. 12.2 and 12.3
{
   double dStep=1.0/m,
      dXn;

   m_dHistScale = static_cast<double>(m)/static_cast<double>(n);
   m_vdRange.assign(m+1,0.0);
   m_vnCount.assign(m,0);
   for (register int i=1; i<=m; ++i)
      m_vdRange[i] = i*dStep;   //  define histogram range
   for (register int i=0; i<n; ++i) {
      dXn = rand();
      for (register int j=0; j<m; ++j) {
         if (dXn >= m_vdRange[j] && dXn < m_vdRange[j+1]) {
            ++m_vnCount[j];    //  bin count
            break;
         }
      }
   }
}
