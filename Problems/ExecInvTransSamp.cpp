//  Execution of the CInversTransSamp class, Section 13.2

#include "stdafx.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <cmath>
#include "matrix2d.h"
#include "InversTransSamp.h"

using namespace std;

static double x_dLambda;  //  \lambda of Eq. (13.24)

/*****************************************************************************/
static double InvPx(double x)
/*****************************************************************************/
//  Inverse function to sample, (Eq. 13.26)
//  input parameter:
//  x ... random number x \in [0,1)
{
   return -x_dLambda*log(1.0-x);
//  if x \in (0,1] please use
//   return -x_dLambda*log(x)
}
/*****************************************************************************/
int ExecuteInversTransSamp()
/*****************************************************************************/
{
   ofstream ofLst;
//   double dLambda=1.0;  //  parameter \lambda of Eq. 13.24
   double dLambda=5.0;
   int nNum=100000,       //  number of random numbers to generate
      nMmax=50;           //  number of histogram bins
//   string sListFile("D:\\CompPhys\\rs1.dat");
        //  data file for Fig. 13.1a, \lambda = 1.0
   string sListFile("D:\\CompPhys\\rs2.dat");
        //  data file for Fig. 13.1b, \lambda = 5.0
   CInversTransSamp I(dLambda,InvPx);

   x_dLambda = dLambda;   //  set static var. for \lambda
   I.seed(131);           //  set seed of system random number generator
   I.Histogram(nMmax,nNum);  //  create histogram, Fig. 13.1
   ofLst.open(sListFile,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile;
      return 1;
   }
   ofLst << fixed << setprecision(6);
   for (register int i=0; i<nMmax; ++i) {
      ofLst << setw(12) << I.m_vdRange[i] << setw(12)
            << I.m_vnCount[i]*I.m_dHistScale << endl;
      ofLst << setw(12) << I.m_vdRange[i+1] << setw(12)
            << I.m_vnCount[i]*I.m_dHistScale << endl;
      ofLst << setw(12) << I.m_vdRange[i+1] << setw(12) << 0.0 << endl;
   }
   ofLst.close();

   return 0;
}
