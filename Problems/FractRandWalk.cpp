//  Realization of the CFractRandWalk class, Section 17.4, fractal random walk

#include "stdafx.h"
#include <stdlib.h>
#include <cmath>
#include <random>
#include <vector>
#include "FractRandWalk.h"

using namespace std;

/*****************************************************************************/
CFractRandWalk::CFractRandWalk(void)
/*****************************************************************************/
//  empty constructor, please use the property assign to define the
//  system's parameters
{
}
/*****************************************************************************/
CFractRandWalk::CFractRandWalk(int n, double b, double t)
/*****************************************************************************/
//  input parameters:
//  n ... number of time steps
//  b ... parameter \beta, Eq. (17.81)
//  t ... minimal waighting time, Eq. (17.81)
{
   m_nN = n;
   m_dBeta = b;
   m_dInvBeta = 1.0/m_dBeta;
   m_dTau = t;
   m_vdX.assign(m_nN,0.0);
   m_vdT.assign(m_nN,0.0);
}
/*****************************************************************************/
CFractRandWalk::~CFractRandWalk(void)
/*****************************************************************************/
{
}
/*****************************************************************************/
void CFractRandWalk::assign(int n, double b, double t)
/*****************************************************************************/
//  use this property to define the parameters of the Lévy flight if the empty
//  constructor has been used to instantiate the CFractRandWalk class.

//  input parameters:
//  n ... number of time steps
//  b ... parameter \beta, Eq. (17.81)
//  t ... minimal waighting time, Eq. (17.81)
{
   m_nN = n;
   m_dBeta = b;
   m_dInvBeta = 1.0/m_dBeta;
   m_dTau = t;
   m_vdX.assign(m_nN,0.0);
   m_vdT.assign(m_nN,0.0);
}
/*****************************************************************************/
void CFractRandWalk::RunFractRandWalk(void)
/*****************************************************************************/
//  generate the curves of Fig. 17.9
{
   double dDeltaT,r;

//  See for explanation:
//  Nicolai M. Josuttis: The C++ Standard Library, Second Edition
//  Addison Wesley (2012)
//  ISBN: 978-0-321-62321-8
//  Chapter 17
//
//  If this doesn't work on UNIX (LINUX) please use the CNormDist class to
//  generate random numbers which follow a normal distribution with mean zero
//  and variance one

#ifdef _MSC_VER
   typedef std::ranlux64_base_01 Myeng; 
#else
   typedef std::default_random_engine Myeng;
#endif
   typedef std::normal_distribution<double> Mydist; 
   static Myeng eng; 
   static Mydist dist(0.0, 1.0);      //  normal distribution, mean zero and variance one 
   static Mydist::input_type engval = eng(); 
   static Mydist::result_type distval = dist(eng);

   distval = distval;  // to quiet "unused" warnings 
   engval = engval; 
   
   dist.reset(); // discard any cached values 

   for (register int i=1; i<m_nN; ++i) {
      m_vdX[i] = m_vdX[i-1]+dist(eng);  //  new x-coordinate
      r = RandInt();
      if (r > 0.9999)
         r = RandInt();
      dDeltaT = m_dTau/pow(1.0-r,m_dInvBeta);
      m_vdT[i] = m_vdT[i-1]+dDeltaT;    //  new time instance
   }
}
