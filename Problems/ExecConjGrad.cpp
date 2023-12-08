//  Execute the CConjGrad class, Appendix H: Deterministic Optimization

#include "stdafx.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <cmath>
#include "matrix2d.h"
#include "ConjGrad.h"

using namespace std;

/*****************************************************************************/
static double f(double x, double y)
/*****************************************************************************/
//  Definition of the function f(x,y) according to Eq. H.6
//  input parameters:
//  x ... x-coordinate
//  y ... y-coordinate
//  on return:
//  the function value f(x,y)
{
   return cos(2.0*x)+sin(4.0*y)+exp(1.5*x*x+0.7*y*y)+2.0*x;
}
/*****************************************************************************/
static double GradFx(double x, double y)
/*****************************************************************************/
//  First partial derivative of the function f(x,y) (Eq. H.6) with respect to x,
//  Eq. H.7
//  input parameters:
//  x ... x-coordinate
//  y ... y-coordinate
//  on return:
//  the function value \frac{\partial f(x,y)}{\partial x}
{
   return -2.0*sin(2.0*x)+3.0*x*exp(1.5*x*x+0.7*y*y)+2.0;
}
/*****************************************************************************/
static double GradFy(double x, double y)
/*****************************************************************************/
//  First partial derivative of the function f(x,y) (Eq. H.6) with respect to y,
//  Eq. H.8
//  input parameters:
//  x ... x-coordinate
//  y ... y-coordinate
//  on return:
//  the function value \frac{\partial f(x,y)}{\partial y}
{
   return 4.0*cos(4.0*y)+1.4*y*exp(1.5*x*x+0.7*y*y);
}
/*****************************************************************************/
static double f1(double x, double y)
/*****************************************************************************/
//  Definitionh of the function f(x,y) according to Eq. H.51
//  input parameters:
//  x ... x-coordinate
//  y ... y-coordinate
//  on return:
//  the function value f(x,y)
{
   return x*x + 10.0*y*y;
}
/*****************************************************************************/
static double GradFx1(double x, double y)
/*****************************************************************************/
//  First partial derivative of the function f(x,y) (Eq. H.51) with respect to x
//  input parameters:
//  x ... x-coordinate
//  y ... y-coordinate
//  on return:
//  the function value \frac{\partial f(x,y)}{\partial x}
{
   return 2.0*x;
}
/*****************************************************************************/
static double GradFy1(double x, double y)
/*****************************************************************************/
//  First partial derivative of the function f(x,y) (Eq. H.51) with respect to y
//  input parameters:
//  x ... x-coordinate
//  y ... y-coordinate
//  on return:
//  the function value \frac{\partial f(x,y)}{\partial y}
{
   return 20.0*y;
}
/*****************************************************************************/
int ExecuteConjGrad()
/*****************************************************************************/
{
   ofstream ofLst;
   int ret=0;
   double dXstart=0.8,  //  define x-ccordinate of the starting point
      dYstart=1.05;     //  define y-ccordinate of the starting point
//   double dXstart=0.8,
//      dYstart=-0.75;
//   double dXstart=-1.05,
//      dYstart=1.05;
//   double dXstart=1.9,
//      dYstart=0.4;
   CConjGrad C(f,GradFx,GradFy);  //  work with Eq. H.6
//   ConjGrad C(f1,GradFx1,GradFy1);  //  work with Eq. H.51
   string sListFile("D:\\CompPhys\\cg1.dat");
        //  data file for Fig. H.2, first path
//   string sListFile("D:\\CompPhys\\cg2.dat");
        //  data file for Fig. H.2, second path
//   string sListFile("D:\\CompPhys\\cg3.dat");
        //  data file for Fig. H.2, third path
//   string sListFile("D:\\CompPhys\\cg4.dat");
        //  data file for Fig. H.3

//   Appendix H, conjugate gradients

   if (ret = C.Iterate(dXstart,dYstart)) {
      if (ret == 1)
         cout << "Error1 in bisection routine" << endl;
      else if (ret == 2)
         cout << "Error2" << endl;
   }
   ofLst.open(sListFile,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile << endl;
      return 1;
   }
   ofLst << fixed << setprecision(5);
   for (register int i=0; i<C.m_nImax; ++i)
      ofLst << setw(12) << C.m_vdXiter[i] << setw(12) << C.m_vdYiter[i] << endl;
   ofLst.close();

   return ret;
}
