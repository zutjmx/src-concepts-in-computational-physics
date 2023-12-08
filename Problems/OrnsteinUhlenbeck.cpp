//  Realization of the COrnsteinUhlenbeck class, Section 17.3

#include "stdafx.h"
#include <cmath>
#include <vector>
#include <random>
#include "OrnsteinUhlenbeck.h"

using namespace std;

/*****************************************************************************/
COrnsteinUhlenbeck::COrnsteinUhlenbeck(void)
/*****************************************************************************/
//  empty constructor, please use the property assign to define the
//  system's parameters
{
}
/*****************************************************************************/
COrnsteinUhlenbeck::COrnsteinUhlenbeck(int n, double beta, double a, double step)
/*****************************************************************************/
//  input parameters:
//  n ...... number of time steps
//  beta ... parameter \beta of Eq. (17.57)
//  a ...... parameter A of Eq. (17.57)
//  step ... time step
{
   m_nSteps = n;
   m_dBeta = beta;
   m_dA = a;
   m_dTstep = step;
   m_dVar = m_dA/(2.0*m_dBeta)*(1.0 - exp(-2.0*m_dBeta*m_dTstep));
      //  variance of the normal distribution, Eq. (17.63)
   m_vdV.assign(m_nSteps,1.0);  //  assign vectors
   m_vdX.assign(m_nSteps,1.0);
   m_vdTime.assign(m_nSteps,0.0);
}
/*****************************************************************************/
COrnsteinUhlenbeck::~COrnsteinUhlenbeck(void)
/*****************************************************************************/
{
}
/*****************************************************************************/
void COrnsteinUhlenbeck::assign(int n, double beta, double a, double step)
/*****************************************************************************/
//  use this property to define the parameters of the Ornstein-Uhlenbeck process
//  if the empty constructor has been used to instantiate the COrnsteinUhlenbeck
//  class.

//  input parameters:
//  n ...... number of time steps
//  beta ... parameter \beta of Eq. (17.57)
//  a ...... parameter A of Eq. (17.57)
//  step ... time step
{
   m_nSteps = n;
   m_dBeta = beta;
   m_dA = a;
   m_dTstep = step;
   m_dVar = m_dA/(2.0*m_dBeta)*(1.0 - exp(-2.0*m_dBeta*m_dTstep));
      //  variance of the normal distribution, Eq. (17.63)
   m_vdV.assign(m_nSteps,1.0);
   m_vdX.assign(m_nSteps,1.0);
   m_vdTime.assign(m_nSteps,0.0);
}
/*****************************************************************************/
void COrnsteinUhlenbeck::RunOrnUhl(double v, double x)
/*****************************************************************************/
//  input parameters:
//  v ... start velocity v_0
//  x ... start coordinate x_0
{
   double dZ;

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
   m_vdV[0] = v; // save intial values
   m_vdX[0] = x;
   for (register int n=1; n<m_nSteps; ++n) {
      dZ = sqrt(m_dVar)*dist(eng);  //  get the correct variance
      m_vdV[n] = m_vdV[n-1]*exp(-m_dBeta*m_dTstep)+dZ;
                    //  velocities v(t), Eq. (17.61)
      m_vdX[n] = m_vdX[n-1]+m_dTstep*m_vdV[n];
                    //  trajectory x(t)
   }
}
