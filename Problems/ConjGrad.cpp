//  Realization of the CConjGrad class, Appendix H: Deterministic Optimization

#include "stdafx.h"
#include <cmath>
#include <vector>
#include "ConjGrad.h"

using namespace std;

/*****************************************************************************/
CConjGrad::CConjGrad(void)
/*****************************************************************************/
//  empty constructor, please use the property assign to define the
//  function f(x,y) and its first derivatives for search of the global minimum
{
   m_Fxy = m_GradFx = m_GradFy = 0;
   m_vdXiter.assign(100,1.0);
   m_vdYiter.assign(100,1.0);
}
/*****************************************************************************/
CConjGrad::CConjGrad(double (*f)(double,double), double (*gfx)(double,double),
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
CConjGrad::~CConjGrad(void)
/*****************************************************************************/
{
}
/*****************************************************************************/
void CConjGrad::assign(double (*f)(double,double), double (*gfx)(double,double),
   double (*gfy)(double,double))
/*****************************************************************************/
//  use this property to initialze the to define the function f(x,y) and its
//  first partial derivatives for the search of the global minimum if an empty
//  constructor has been used to instantiate the CConjGrad class.

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
int CConjGrad::Iterate(double x, double y)
/*****************************************************************************/
//  Execute iteration scheme, page 349. Works for convex functions f(x,y)
//  as well as function with more than one minimum
//  input parameters:
//  x ... x-coordinate of the starting point
//  y ... y-coordinate of the starting point
{
   int i=1,nIter;
   double dEta,dXalt,dYalt,dXaltAlt,dYaltAlt,dpPsik[2],dpGrad[2],dLambdaa,
      dLambdab,dLambdac,dEps,dGradFxNeua,dGradFxNeub,dGradFxNeuc,dLambdak,
      dpGradAlt[2],dBeta,z,zz,dXneu,dYneu;

   m_vdXiter[0] = x;   //  save starting point
   m_vdYiter[0] = y;
   dEta = 1.0;         //  make sure the iteration starts
   while (dEta > 1.0e-6) {
      dXaltAlt = dXalt = m_vdXiter[i-1];
      dYaltAlt = dYalt = m_vdYiter[i-1];
      dpGrad[0] = dpPsik[0] = -(*m_GradFx)(dXalt,dYalt);  //  Eq. H.49
      dpGrad[1] = dpPsik[1] = -(*m_GradFy)(dXalt,dYalt);
      dLambdaa = 0.0;
      dLambdab = 1.0e-4;
      dEps = 1.0;
      nIter = 1;
      while (dEps > 1.0e-6) {
         dGradFxNeua = (*m_GradFx)(dXalt+dLambdaa*dpPsik[0],dYalt+dLambdaa*dpPsik[1])*dpPsik[0]+
                       (*m_GradFy)(dXalt+dLambdaa*dpPsik[0],dYalt+dLambdaa*dpPsik[1])*dpPsik[1];
         dGradFxNeub = (*m_GradFx)(dXalt+dLambdab*dpPsik[0],dYalt+dLambdab*dpPsik[1])*dpPsik[0]+
                       (*m_GradFy)(dXalt+dLambdab*dpPsik[0],dYalt+dLambdab*dpPsik[1])*dpPsik[1];
         if (dGradFxNeua*dGradFxNeub > 0)
            dLambdab *= 2.0;   //  increase \lambda_b
         else if (dGradFxNeua*dGradFxNeub < 0) {  //  bisection routine
            dLambdac = 0.5*(dLambdaa+dLambdab);
            dGradFxNeuc = (*m_GradFx)(dXalt+dLambdac*dpPsik[0],dYalt+dLambdac*dpPsik[1])*dpPsik[0]+
                       (*m_GradFy)(dXalt+dLambdac*dpPsik[0],dYalt+dLambdac*dpPsik[1])*dpPsik[1];
            if (dGradFxNeua*dGradFxNeuc < 0)
               dLambdab = dLambdac;
            else
               dLambdaa = dLambdac;
         }
         dGradFxNeua = (*m_GradFx)(dXalt+dLambdaa*dpPsik[0],dYalt+dLambdaa*dpPsik[1])*dpPsik[0]+
                       (*m_GradFy)(dXalt+dLambdaa*dpPsik[0],dYalt+dLambdaa*dpPsik[1])*dpPsik[1];
         dEps = fabs(dGradFxNeua);
         if (++nIter > 50) {
            dEps = 0.0;
            dLambdaa = 0.0;
            return 2;    //  Error 2
         }
      }
      dLambdak = dLambdaa;
      dXalt += dLambdak*dpPsik[0];   //  Eq. H.47
      dYalt += dLambdak*dpPsik[1];
      if ((*m_Fxy)(dXaltAlt,dYaltAlt) < (*m_Fxy)(dXalt,dYalt))
         return 1;    //  Error 1
      dpGradAlt[0] = dpGrad[0];
      dpGradAlt[1] = dpGrad[1];
      dpGrad[0] = -(*m_GradFx)(dXalt,dYalt);
      dpGrad[1] = -(*m_GradFy)(dXalt,dYalt);
      z = dpGrad[0]*dpGrad[0]+dpGrad[1]*dpGrad[1];
      zz =dpGradAlt[0]*dpGradAlt[0]+dpGradAlt[1]*dpGradAlt[1];
      dBeta = z/zz;
      dpPsik[0] = dpGrad[0]+dBeta*dpPsik[0];
      dpPsik[1] = dpGrad[1]+dBeta*dpPsik[1];
      dLambdaa = 0.0;
      dLambdab = 1.0e-5;
      dEps = 1.0;
      nIter = 1;
      while (dEps > 1.0e-6) {
         dGradFxNeua = (*m_GradFx)(dXalt+dLambdaa*dpPsik[0],dYalt+dLambdaa*dpPsik[1])*dpPsik[0]+
                       (*m_GradFy)(dXalt+dLambdaa*dpPsik[0],dYalt+dLambdaa*dpPsik[1])*dpPsik[1];
         dGradFxNeub = (*m_GradFx)(dXalt+dLambdab*dpPsik[0],dYalt+dLambdab*dpPsik[1])*dpPsik[0]+
                       (*m_GradFy)(dXalt+dLambdab*dpPsik[0],dYalt+dLambdab*dpPsik[1])*dpPsik[1];
         if (dGradFxNeua*dGradFxNeub > 0)
            dLambdab *= 2.0;   //  increase \lambda_b
         else if (dGradFxNeua*dGradFxNeub < 0) {  //  bisection routine
            dLambdac = 0.5*(dLambdaa+dLambdab);
            dGradFxNeuc = (*m_GradFx)(dXalt+dLambdac*dpPsik[0],dYalt+dLambdac*dpPsik[1])*dpPsik[0]+
                       (*m_GradFy)(dXalt+dLambdac*dpPsik[0],dYalt+dLambdac*dpPsik[1])*dpPsik[1];
            if (dGradFxNeua*dGradFxNeuc < 0)
               dLambdab = dLambdac;
            else
               dLambdaa = dLambdac;
         }
         dGradFxNeua = (*m_GradFx)(dXalt+dLambdaa*dpPsik[0],dYalt+dLambdaa*dpPsik[1])*dpPsik[0]+
                       (*m_GradFy)(dXalt+dLambdaa*dpPsik[0],dYalt+dLambdaa*dpPsik[1])*dpPsik[1];
         dEps = fabs(dGradFxNeua);
         if (++nIter > 50) {
            dEps = 0.0;
            dLambdaa = 0.0;
         }
      }
      dLambdak = dLambdaa;
      dXneu = dXalt+dLambdak*dpPsik[0];  //  Eq. H.47
      dYneu = dYalt+dLambdak*dpPsik[1];
      if ((*m_Fxy)(dXneu,dYneu) > (*m_Fxy)(dXalt,dYalt)) {
         dXneu = dXalt;
         dYneu = dYalt;
      }
      m_vdXiter[i] = dXneu;
      m_vdYiter[i++] = dYneu;
      dEta = fabs((*m_Fxy)(dXneu,dYneu)-(*m_Fxy)(dXaltAlt,dYaltAlt));
   }
   m_nImax = i-1;
   return 0;
}
