//  Realization of the CHeatEquation class, Chapter 9

#include "stdafx.h"
#include <cmath>
#include "HeatEquation.h"

using namespace std;

/*****************************************************************************/
CHeatEquation::CHeatEquation(void)
/*****************************************************************************/
//  empty constructor, please use the property assign to set the
//  heat equation's parameters
{
}
/*****************************************************************************/
CHeatEquation::~CHeatEquation(void)
/*****************************************************************************/
{
}
/*****************************************************************************/
CHeatEquation::CHeatEquation(int n, double t0, double tn, double l,
   double kap)
/*****************************************************************************/
//  input parameters:
//  n .... # of steps
//  t0 ... temperature T_0, see Fig. 9.1
//  tn ... temperature T_N, see Fig. 9.1
//  l .... length of the rod
//  kap .. thermal diffusivity
{
   m_nN = n;
   m_dT0 = t0;
   m_dTN = tn;
   m_dKappa = kap;
   m_dL = l;
   m_dH = m_dL/m_nN;
   m_dH2 = m_dH*m_dH;
   m_vGridPoint.assign(n+1,0.0);  //  define the grid points
   for (register int i=1; i<=n; ++i)
      m_vGridPoint[i] = m_vGridPoint[i-1]+m_dH;
   m_vDiag.assign(n,-2.0);        //  define the diagonal of matrix A,
                                  //  Eq. (9.8)
   m_vSuperd.assign(n,1.0);       //  define the super-diagonal and the
   m_vSubd.assign(n,1.0);         //  sub-diagonal of matrix , Eq. (9.8)
   m_vTemp.assign(n,0.0);
}
/*****************************************************************************/
void CHeatEquation::assign(int n, double t0, double tn, double l,
   double kap)
/*****************************************************************************/
//  use this property to initialze the heat equation's parameters if an
//  empty constructor has been used to instantiate the CHeatEquation class

//  input parameters:
//  n .... # of steps
//  t0 ... temperature T_0, see Fig. 9.1
//  tn ... temperature T_N, see Fig. 9.1
//  l .... length of the rod
//  kap .. thermal diffusivity
{
   m_nN = n;
   m_dT0 = t0;
   m_dTN = tn;
   m_dKappa = kap;
   m_dL = l;
   m_dH = m_dL/m_nN;              //  step size
   m_dH2 = m_dH*m_dH;
   m_vGridPoint.assign(n+1,0.0);  //  define the grid points
   for (register int i=1; i<=n; ++i)
      m_vGridPoint[i] = m_vGridPoint[i-1]+m_dH;
   m_vDiag.assign(n,-2.0);        //  define the diagonal of matrix A,
                                  //  Eq. (9.8)
   m_vSuperd.assign(n,1.0);       //  define the super-diagonal and the
   m_vSubd.assign(n,1.0);         //  sub-diagonal of matrix , Eq. (9.8)
   m_vTemp.assign(n,0.0);
}
/*****************************************************************************/
void CHeatEquation::SetHeatSource(double theta, double ell, double center)
/*****************************************************************************/
//  define the heat source/drain according to Eq. (9.21)
//  input parameters:
//  theta ... maximum height
//  ell ..... width of the heat source/drain
//  center .. position of the heat source/drain
{
   double fact,ell2,fact1;

   m_dTheta = theta;
   m_dEll = ell;
   ell2 = ell*ell;
   fact = theta/ell;
   fact1 = m_dH2/m_dKappa;
   m_vGamma.assign(m_nN,0.0);
   m_vRhs.assign(m_nN,0.0);
   for (register int i=0; i<m_nN; ++i) {
      double x = m_vGridPoint[i+1]-center;
      x = x*x;
      m_vGamma[i] = fact*exp(-x/ell);    //  Eq. (9.21)
   }
   m_vRhs[0] = fact1*m_vGamma[0]-m_dT0;  //  calculate the vector F, Eq. (9.24)
   for (register int i=1; i<m_nN-1; ++i)
      m_vRhs[i] = fact1*m_vGamma[i];
   m_vRhs[m_nN-1] = fact1*m_vGamma[m_nN-1]-m_dTN;
}
/*****************************************************************************/
int CHeatEquation::TriDiag()
/*****************************************************************************/
//  solve Eq. (9.7) with tri-diagonal matrix A
//  on return:
//  = 1 ... first element of the diagonal is zero
//  = 2 ... the tri-diagonal matrix is ill conditioned
{
   int j;
   double bet;

   int n=m_vSubd.size();
   vector<double> gam(n);
   if (m_vDiag[0] == 0.0) {
      return 1;          //   Error 1
   }
   m_vTemp[0] = m_vRhs[0]/(bet = m_vDiag[0]);
   for (j=1; j<n; ++j) {
      gam[j] = m_vSuperd[j-1]/bet;
      bet = m_vDiag[j]-m_vSubd[j]*gam[j];
      if (bet == 0.0) {
         return 2;    //  Error 2
      }
      m_vTemp[j] = (m_vRhs[j]-m_vSubd[j]*m_vTemp[j-1])/bet;
   }
   for (j=(n-2); j>=0; --j)
      m_vTemp[j] -= gam[j+1]*m_vTemp[j+1];
   return 0;
}
/*****************************************************************************/
int CHeatEquation::SolveEquation(void)
/*****************************************************************************/
{
   int ret = 0;

   ret = TriDiag();
   return ret;
}
