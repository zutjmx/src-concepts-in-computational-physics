//  Execution of the CRejectSamp class, Section 13.3 Rejection metod

#include "stdafx.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <cmath>
#include "matrix2d.h"
#include "RejectSamp.h"

using namespace std;

//  use these parameters to define constants used in the functions
//  qx(x), hx(x), ihx(x)
static double x_dLambda,  //  parameter \lambda, Eq. (13.24)
   x_dSigma,              //  parameter \sigma, Eq. (13.33)
   x_dPi;                 //  number Pi in machine accuracy

/*****************************************************************************/
static double qx(double x)
/*****************************************************************************/
//  function q(x), Eq. (13.33)
{
   extern double x_dSigma,x_dPi;
   double sig2=x_dSigma*x_dSigma;
   return sqrt(2.0/(x_dPi*sig2))*exp(-x*x/(2.0*sig2));
}
/*****************************************************************************/
static double hx(double x)
/*****************************************************************************/
//  envelop function, exponential distribution Eq. (13.24)
{
   extern double x_dLambda;
   return 1.0/x_dLambda*exp(-x/x_dLambda);
}
/*****************************************************************************/
static double ihx(double x)
/*****************************************************************************/
//  inverse transformation of exp. function, Eq. (13.26)
{
   extern double x_dLambda;
   return -x_dLambda*log(1.0-x); //  Eq. (13.27)
}
/*****************************************************************************/
int ExecuteRejectSamp()
/*****************************************************************************/
{
   ofstream ofLst;
   double dLambda=1.0,  //  optimum parameter \lambda, Eq. (13.40)
      dSigma=1.0,       //  is equal to the parameter \sigma, Eq. (13.33)
      dCmin;            //  c_min, Eq. (13.41)
//   int nNum=1000,     //  number of random numbers generated
//   int nNum=10000,
   int nNum=100000,
      nMmax=50;
//   string sListFile("D:\\CompPhys\\rj1.dat");
        //  data file for Fig. 13.3a, nNum = 1000
//   string sListFile("D:\\CompPhys\\rj2.dat");
        //  data file for Fig. 13.3b, nNum = 10000
   string sListFile("D:\\CompPhys\\rj3.dat");
        //  data file for Fig. 13.3c, nNum = 100000

//   Chapter 13, Rejection Method

   x_dLambda = dLambda;        //  set external variables
   x_dSigma = dSigma;
   x_dPi = 4.0*atan(1.0);       //  number Pi
   dCmin = sqrt(2.0*exp(1.0)/x_dPi);  // c_min, Eq. (13.41)
   CRejectSamp R(dLambda,dSigma,dCmin,qx,hx,ihx);
   R.seed(131);                 //  set seed
   R.Histogram(nMmax,nNum);     //  generate the histogram Fig. 13.3
   ofLst.open(sListFile,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile << endl;
      return 1;
   }
   ofLst << fixed << setprecision(6);
   for (register int i=0; i<nMmax; ++i) {
      ofLst << setw(12) << R.m_vdRange[i] << setw(12)
            << R.m_vnCount[i]*R.m_dHistScale << endl;
      ofLst << setw(12) << R.m_vdRange[i+1] << setw(12)
            << R.m_vnCount[i]*R.m_dHistScale << endl;
      ofLst << setw(12) << R.m_vdRange[i+1] << setw(12) << 0.0 << endl;
   }
   ofLst.close();

   return 0;
}
