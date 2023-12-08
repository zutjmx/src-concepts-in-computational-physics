//  Realization of the CSteepDesc class, Appendix H: Deterministic Optimization

#include "stdafx.h"
#include <cmath>
#include <vector>
#include "SteepDesc.h"

using namespace std;

/*****************************************************************************/
CSteepDesc::CSteepDesc(void)
/*****************************************************************************/
//  empty constructor, please use the property assign to define the
//  function f(x,y) and its first derivatives for search of the global minimum
{
   m_Fxy = m_GradFx = m_GradFy = 0;
   m_vdXiter.assign(100,1.0);
   m_vdYiter.assign(100,1.0);
}
/*****************************************************************************/
CSteepDesc::CSteepDesc(double (*f)(double,double), double (*gfx)(double,double),
   double (*gfy)(double,double))
/*****************************************************************************/
//  input parameters:
//  f ..... pointer to the function f(x,y)
//  gfx ... pointer to the function which returns the first partial derivative of
//          f(x,y) with respect to x
//  gfy ... pointer to the function which returns the first partial derivative of
//          f(x,y) with respect to y
{
   m_Fxy = f;
   m_GradFx = gfx;
   m_GradFy = gfy;
   m_vdXiter.assign(100,1.0);
   m_vdYiter.assign(100,1.0);
}
/*****************************************************************************/
CSteepDesc::~CSteepDesc(void)
/*****************************************************************************/
{
}
/*****************************************************************************/
void CSteepDesc::assign(double (*f)(double,double), double (*gfx)(double,double),
   double (*gfy)(double,double))
/*****************************************************************************/
//  use this property to initialze the to define the function f(x,y) and its
//  first partial derivatives for the search of the global minimum if an empty
//  constructor has been used to instantiate the CSteepDesc class.

//  input parameters:
//  f ..... pointer to the function f(x,y)
//  gfx ... pointer to the function which returns the first partial derivative of
//          f(x,y) with respect to x
//  gfy ... pointer to the function which returns the first partial derivative of
//          f(x,y) with respect to y
{
   m_Fxy = f;
   m_GradFx = gfx;
   m_GradFy = gfy;
   m_vdXiter.assign(100,1.0);
   m_vdYiter.assign(100,1.0);
}
/*****************************************************************************/
int CSteepDesc::Iterate(double x, double y)
/*****************************************************************************/
//  Execute iteration scheme, page 343
//  input parameters:
//  x ... x-coordinate of the starting point
//  y ... y-coordinate of the starting point
{
   double dXalt,dYalt,dGradFxAlt,dGradFyAlt,dEps,dAlphaa,dAlphab,
      dAlphac,dXneua,dYneua,dXneub,dYneub,dXneuc,dYneuc,dGradFxNeua,
      dGradFyNeua,dGradFxNeub,dGradFyNeub,dSkalarpa,dSkalarpb,dSkalarpc,
      dGradFxNeuc,dGradFyNeuc,dAlpha,dEta;
   int i=1,nIter;

   m_vdXiter[0] = x;  //  save starting point
   m_vdYiter[0] = y;
   dEta = 1.0;        //  make sure iteration starts
   while (dEta > 1.0e-6) {  //  start iteration, Page 343
      dXalt = m_vdXiter[i-1];  //  step 1
      dYalt = m_vdYiter[i-1];
      dGradFxAlt = (*m_GradFx)(dXalt,dYalt);  //  step 2
      dGradFyAlt = (*m_GradFy)(dXalt,dYalt);
      dEps = 1.0;   //  make sure bisection starts
      dAlphaa = 0.0;  //  step a
      dAlphab = 0.04;
      nIter = 1;
      while (dEps > 1.0e-5) {  //  bisection
         dXneua = dXalt-dAlphaa*dGradFxAlt;
         dYneua = dYalt-dAlphaa*dGradFyAlt;
         dGradFxNeua = (*m_GradFx)(dXneua,dYneua); //  Eq. H.10
         dGradFyNeua = (*m_GradFy)(dXneua,dYneua);
         dXneub = dXalt-dAlphab*dGradFxAlt;
         dYneub = dYalt-dAlphab*dGradFyAlt;
         dGradFxNeub = (*m_GradFx)(dXneub,dYneub);
         dGradFyNeub = (*m_GradFy)(dXneub,dYneub);
         dSkalarpa = dGradFxNeua*dGradFxAlt+dGradFyNeua*dGradFyAlt;
         dSkalarpb = dGradFxNeub*dGradFxAlt+dGradFyNeub*dGradFyAlt;
         dEps = fabs(dSkalarpa);       //  new accuracy
         if (dSkalarpa*dSkalarpb > 0)  //  bad choice, increase \alpha^b_n
            dAlphab = dAlphab*2.0;     //  step b
         else if (dSkalarpa*dSkalarpb < 0) {  //  step d
            dAlphac = (dAlphaa+dAlphab)*0.5;  //  Eq. H.11
            dXneuc = dXalt-dAlphac*dGradFxAlt;
            dYneuc = dYalt-dAlphac*dGradFyAlt;
            dGradFxNeuc = (*m_GradFx)(dXneuc,dYneuc);
            dGradFyNeuc = (*m_GradFy)(dXneuc,dYneuc);
            dSkalarpc = dGradFxNeuc*dGradFxAlt+dGradFyNeuc*dGradFyAlt;
            if (dSkalarpc*dSkalarpa < 0)
                dAlphab = dAlphac;
            else if (dSkalarpc*dSkalarpb < 0)
                dAlphaa = dAlphac;
            else
                return 1;  //  there is a problem!
         }
         if (++nIter > 10000)
            return 2;      //  bisection doesn't converge
      }
      dAlpha = dAlphaa;    //  step 4
      m_vdXiter[i] = m_vdXiter[i-1]-dAlpha*(*m_GradFx)(m_vdXiter[i-1],m_vdYiter[i-1]);
      m_vdYiter[i] = m_vdYiter[i-1]-dAlpha*(*m_GradFy)(m_vdXiter[i-1],m_vdYiter[i-1]);
           //  new iteration points
      double x=fabs(m_vdXiter[i]-dXalt),
         y=fabs(m_vdYiter[i]-dYalt);
      dEta = x*x+y*y;  //  calculate new iteration accuracy
      ++i;
   }
   m_nImax = i-1;
   return 0;
}
