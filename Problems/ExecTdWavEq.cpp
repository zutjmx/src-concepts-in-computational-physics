//  Solve the wave equation, Section 11.4

#include "stdafx.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <cmath>
#include "matrix2d.h"
#include "TdWaveEq.h"

using namespace std;

/*****************************************************************************/
int ExecuteTdWaveEq()
/*****************************************************************************/
{
   ofstream ofLst;
   int nN=100,       //  number of grid-points
      nK=301,        //  max. number of time steps
//      nStep = 0;   //  number of time steps
//      nStep=25;
//      nStep=50;
//      nStep=100;
//      nStep=150;
      nStep=200;
   double dC=2.0,    //  speed of wave propagation
      dL=1.0,        //  length of the string
//      dLambda=0.5, //  CFL stability criterion, Eq. (11.31)
      dLambda=1.01,
      dA=1.0;        //  amplitude
   CTdWaveEq T(nN,nK,dC,dL,dLambda,dA);
//   string sListFile("D:\\CompPhys\\TdWave11.dat");
         //  save results for nStep = 0, dLambda=0.5, Fig. 11.6
//   string sListFile("D:\\CompPhys\\TdWave12.dat");
         //  save results for nStep = 25, dLambda=0.5, Fig. 11.6
//   string sListFile("D:\\CompPhys\\TdWave13.dat");
         //  save results for nStep = 50, dLambda=0.5, Fig. 11.6
//   string sListFile("D:\\CompPhys\\TdWave14.dat");
         //  save results for nStep = 100, dLambda=0.5, Fig. 11.6
//   string sListFile("D:\\CompPhys\\TdWave15.dat");
         //  save results for nStep = 150, dLambda=0.5, Fig. 11.6
//   string sListFile("D:\\CompPhys\\TdWave16.dat");
         //  save results for nStep = 200, dLambda=0.5, Fig. 11.6
//   string sListFile("D:\\CompPhys\\TdWave22.dat");
         //  save results for nStep = 25, dLambda=1.01, Fig. 11.7
//   string sListFile("D:\\CompPhys\\TdWave23.dat");
         //  save results for nStep = 50, dLambda=1.01, Fig. 11.7
//   string sListFile("D:\\CompPhys\\TdWave24.dat");
         //  save results for nStep = 100, dLambda=1.01, Fig. 11.7
//   string sListFile("D:\\CompPhys\\TdWave25.dat");
         //  save results for nStep = 150, dLambda=1.01, Fig. 11.7
   string sListFile("D:\\CompPhys\\TdWave26.dat");
         //  save results for nStep = 200, dLambda=1.01, Fig. 11.7

//  Chapter 11,  Time Dependent Wave Equation

   T.InitBoundCond();       //  set initial and boundary conditions
   T.SolveEquation(nStep);  //  solve for nStep time steps

//     save results on file

   ofLst.open(sListFile,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile << endl;
      return 1;
   }
   ofLst << scientific << setprecision(5);
   for(register int i=0; i<nN; ++i) {
      ofLst << setw(15) << T.m_vX[i] << setw(15) << T.m_mU[i][nStep] << endl;
   }
   ofLst.close();
   return 0;
}
