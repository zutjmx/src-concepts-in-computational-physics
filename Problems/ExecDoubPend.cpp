//  Execute the double pendulum problem, Chapter 6

#include "stdafx.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <cmath>
#include "DoublePendulum.h"

using namespace std;

/*****************************************************************************/
int ExecuteDoublePend()
/*****************************************************************************/
{
   ofstream ofLst;
   int nMaxSteps=60000;  //  maximum # of steps
//   int nMaxSteps=24000;
//   int nMaxSteps=36000;
//   int nMaxSteps=480000;  //  for Poincare plots
   double dDeltaT=0.001, //  time step
      dG=9.8067,         //  acceleration due to gravity
      dMass=1.0,         //  pendulum's mass
      dL=1.0;            //  length

//  Chapter 6, Double Pendulum

   CDoublePendulum P(nMaxSteps,dDeltaT,dMass,dL,dG);
//   string sListFile("D:\\CompPhys\\dp1.dat");
         //  file for results Fig. 6.2
//   string sListFile("D:\\CompPhys\\dp2.dat");
         //  file for results Fig. 6.3
//   string sListFile("D:\\CompPhys\\dp3.dat");
         //  file for results Fig. 6.4
//   string sListFile("D:\\CompPhys\\dp4.dat");
   string sListFile("D:\\CompPhys\\dp5.dat");
         //  file for results Fig. 6.5
//   string sListFile("D:\\CompPhys\\dp6.dat");
         //  file for results Fig. 6.6
//   string sListFile("D:\\CompPhys\\pdp1.dat");
         //  file for results Fig. 6.9
//   string sListFile("D:\\CompPhys\\pdp2.dat");
         //  file for results Fig. 6.10
//   string sListFile("D:\\CompPhys\\pdp5.dat");
         //  file for results Fig. 6.11

//  set initial values

//   P.SetInitialValues(0.0,0.0,4.0,2.0);  // initial values Fig. 6.2, Fig. 6.9
//   P.SetInitialValues(1.0,0.0,0.0,3.0);  // initial values Fig. 6.3
//   P.SetInitialValues(0.0,0.0,0.0,4.0);  // initial values Fig. 6.4, Fig. 6.11
//   P.SetInitialValues(0.0,0.0,0.0,5.0);  // initial values Fig. 6.5
//   P.SetInitialValues(0.0,0.0,0.0,6.5);  // initial values Fig. 6.6
//   P.SetInitialValues(1.0,0.0,0.0,3.0);  // initial values Fig. 6.10

//   Execute Runge Kutta

   P.RungeKutta();

//   save results on file

   ofLst.open(sListFile,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile << endl;
      return 1;
   }
   for (register int i=0; i<=nMaxSteps; ++i)
      ofLst << fixed << setw(12) << setprecision(5) << P.m_vdT[i]
            << setw(12) << setprecision(5) << P.m_vdPhi1[i]
            << setw(12) << setprecision(5) << P.m_vdPhi2[i]
            << setw(12) << setprecision(5) << P.m_vdP1[i]
            << setw(12) << setprecision(5) << P.m_vdP2[i] << endl;
   ofLst.close();

//  Execute Poincare: create a Poincare plot

   //P.Poincare();

//   save results on file

   //ofLst.open(sListFile,ios::out);
   //if (!ofLst) {
   //   cout << "Can't open file: " << sListFile << endl;
   //   return 1;
   //}
   //for (register int i=0; i<static_cast<int>(P.m_vdPhi1.size()); ++i)
   //   ofLst << fixed << setw(12) << setprecision(5) << P.m_vdPhi1[i]
   //         << setw(12) << setprecision(5) << P.m_vdP1[i] << endl;
   //ofLst.close();

   return 0;
}
