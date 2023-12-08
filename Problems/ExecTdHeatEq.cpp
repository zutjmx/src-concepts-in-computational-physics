//  Solve the time-dependent heat equation, Section 11.3

#include "stdafx.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <cmath>
#include "matrix2d.h"
#include "TdHeatEq.h"

using namespace std;

/*****************************************************************************/
int ExecuteTdHeatEq()
/*****************************************************************************/
{
   ofstream ofLst;
   int nN=20,        //  number of grid-points
      nK=2000,       //  max. array size
//      nTmax=1399,     //  max. number of time steps
      nTmax=950,
//      nStep=25;    //  number of time steps, Figs. 11.3, 11.4, 11.5
//      nStep=50;
//      nStep=100;
//      nStep=150;
//      nStep=200;
      nStep=300;
   double dL=20.0,   //  length of the rod
      dKappa=1.0,    //  thermal diffusivity
      dDeltaT,       //  \Delta t
      dClf,          //  stability criterion, Eq. (11.20)
      dH;            //  distance between grid-points squared
//   string sListFile("D:\\CompPhys\\TdHeatEqEx12.dat");
         //  Explicit Euler: save results for nStep = 25, Fig. 11.3
//   string sListFile("D:\\CompPhys\\TdHeatEqEx13.dat");
         //  Explicit Euler: save results for nStep = 50, Fig. 11.3
//   string sListFile("D:\\CompPhys\\TdHeatEqEx14.dat");
         //  Explicit Euler: save results for nStep = 100, Fig. 11.3
//   string sListFile("D:\\CompPhys\\TdHeatEqEx15.dat");
         //  Explicit Euler: save results for nStep = 150, Fig. 11.3
   string sListFile("D:\\CompPhys\\TdHeatEqEx16.dat");
         //  Explicit Euler: save results for nStep = 300, Fig. 11.3
//   string sListFile("D:\\CompPhys\\TdHeatEqEx22.dat");
         //  Explicit Euler: save results for nStep = 25, Fig. 11.4
//   string sListFile("D:\\CompPhys\\TdHeatEqEx23.dat");
         //  Explicit Euler: save results for nStep = 50, Fig. 11.4
//   string sListFile("D:\\CompPhys\\TdHeatEqEx24.dat");
         //  Explicit Euler: save results for nStep = 100, Fig. 11.4
//   string sListFile("D:\\CompPhys\\TdHeatEqEx25.dat");
         //  Explicit Euler: save results for nStep = 150, Fig. 11.4
//   string sListFile("D:\\CompPhys\\TdHeatEqEx26.dat");
         //  Explicit Euler: save results for nStep = 200, Fig. 11.4
//   string sListFile("D:\\CompPhys\\TdHeatEqIm12.dat");
         //  Implicit Euler: save results for nStep = 25, Fig. 11.5
//   string sListFile("D:\\CompPhys\\TdHeatEqIm13.dat");
         //  Implicit Euler: save results for nStep = 50, Fig. 11.5
//   string sListFile("D:\\CompPhys\\TdHeatEqIm14.dat");
         //  Implicit Euler: save results for nStep = 100, Fig. 11.5
//   string sListFile("D:\\CompPhys\\TdHeatEqIm15.dat");
         //  Implicit Euler: save results for nStep = 150, Fig. 11.5
//   string sListFile("D:\\CompPhys\\TdHeatEqIm16.dat");
         //  Implicit Euler: save results for nStep = 300, Fig. 11.5

//  Chapter 11,  Time Dependent Heat Equation

   dDeltaT = static_cast<double>(nTmax)/static_cast<double>(nK-1);
             //  calculate time step \Delta t
//   cout << fixed << setprecision(5) << setw(12) << dDeltaT << endl;
   dH = dL/(nN-1);            //  distance between grid-points
   dH *= dH;
   dClf = dKappa*dDeltaT/dH;  //  stability criterion Eq. (11.20)
//   cout << fixed << setprecision(5) << setw(12) << dClf << endl;

   CTdHeatEq T(nN,nTmax,nK,dL,dKappa);
   T.ExplicitEuler(nStep);    //  execute explicit Euler integrator
//   T.ImplicitEuler(nStep);    //  execute implicit Euler integrator

//      save results on file

   ofLst.open(sListFile,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile << endl;
      return 1;
   }
   ofLst << scientific << setprecision(5);
   for(register int i=0; i<nN; ++i)
      ofLst << setw(12) << T.m_vdX[i] << setw(15) << T.m_mT[i][nStep] << endl;
   ofLst.close();
   return 0;
}
