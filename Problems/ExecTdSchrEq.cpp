//  Solve the time-dependent Schroedinger equation, Section 11.5

#include "stdafx.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <cmath>
#include "matrix2d.h"
#include "TdSchrodEq.h"

using namespace std;

/*****************************************************************************/
int ExecuteTdSchrodEq()
/*****************************************************************************/
{
   ofstream ofLst;
   int nN= 500,           //  number of grid points
      nM= 2701,           //  max. number of time steps
//      nStep = 500;      //  number of time steps  
//      nStep = 1000;
//      nStep=1500;
//      nStep=2000;
      nStep=2500;
   double dQ0=2.0,       //  momentum q of the wave packet, Eq. (11.77)
      dX0=200.0,         //  center position of the wave packet
      dMass=1.0,         //  particle's mass
      dHbar=1.0,         //  Planck's constant
      dSigma2=800,       //  \sigma^2, width of the wave packet squared
      dDeltaT=0.1,       //  time step
      x;
   vector<double> vdAbsPsi(nN,0.0);
   CTdSchrodEq T(nN,nM,dQ0,dX0,dMass,dHbar,dSigma2,dDeltaT);
//   string sListFile("D:\\CompPhys\\TdSchrod0.dat");
         //  save results for nStep = 0, Fig. 11.8
//   string sListFile("D:\\CompPhys\\TdSchrod11.dat");
         //  save results for nStep = 500, Fig. 11.8
//   string sListFile("D:\\CompPhys\\TdSchrod12.dat");
         //  save results for nStep = 1000, Fig. 11.8
//   string sListFile("D:\\CompPhys\\TdSchrod13.dat");
         //  save results for nStep = 1500, Fig. 11.8
//   string sListFile("D:\\CompPhys\\TdSchrod21.dat");
         //  save results for nStep = 500, Fig. 11.9
//   string sListFile("D:\\CompPhys\\TdSchrod22.dat");
         //  save results for nStep = 1000, Fig. 11.9
//   string sListFile("D:\\CompPhys\\TdSchrod23.dat");
         //  save results for nStep = 1500, Fig. 11.9
//   string sListFile("D:\\CompPhys\\TdSchrod24.dat");
         //  save results for nStep = 2000, Fig. 11.9
   string sListFile("D:\\CompPhys\\TdSchrod25.dat");
         //  save results for nStep = 2500, Fig. 11.9

//   T.SetPotential(1);   //  define potential
   T.SetPotential(2);
   T.InitPsi();           //  initial conditions
   T.SolveEquation(nStep);
   for (register int i=0; i<nN; ++i) {
      x = abs(T.m_mcPsi[nStep][i]);
      vdAbsPsi[i] = x*x;
   }
   ofLst.open(sListFile,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile << endl;
      return 1;
   }
   ofLst << scientific << setprecision(5);
   for (register int i=0; i<nN; ++i)
      ofLst << setw(15) << T.m_vdX[i] << setw(15) << vdAbsPsi[i] << endl;
   ofLst.close();
   return 0;
}
