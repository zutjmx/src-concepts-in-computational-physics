//  Realization of the CNormDist class
//  generate normaly distributed random numbers by direct sampling

#include "stdafx.h"
#include <stdlib.h>
#include <cmath>
#include "NormDist.h"

/*****************************************************************************/
CNormDist::CNormDist(void)
/*****************************************************************************/
{
   m_dPi = 4.0*atan(1.0);   //  number Pi in system accuracy
}
/*****************************************************************************/
CNormDist::~CNormDist(void)
/*****************************************************************************/
{
}
/*****************************************************************************/
void CNormDist::sample(double& z1, double& z2)
/*****************************************************************************/
//  sample the normal distribution with mean zero and variance 1 using
//  direct sampling as described in Sec. 13.1
{
   double r1,r2,x,y;

   r1 = RandInt();
   if (r1 == 0.0)       //  don't allow zero!
      r1 = RandInt();
   r2 = RandInt();
   x = -1.0*log(r1);
   y = 2.0*m_dPi*r2;
   z1 = sqrt(x)*cos(y);  // apply Eq. (13.9)
   z2 = sqrt(x)*sin(y);
}
