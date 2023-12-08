
#include "stdafx.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include "Integrator.h"

using namespace std;

/*****************************************************************************/
int ExecuteKepler()
/*****************************************************************************/
{
   ofstream ofLst;
   double dE=0.6,     //  excentrizity
//      dDeltaT=0.1,
      dDeltaT=0.005,  //  time step
      dMaxT=10.0;     //  maximum time

//  Chapter 5, Kepler Problem

//  File names to store the results:

//   string sListFile("D:\\CompPhys\\ExplicitEuler.dat");
               //  save data for Fig. 5.1a and Fig. 5.3
   string sListFile("D:\\CompPhys\\ImplicitEuler.dat");
               //  save data for Fig. 5.1b and Fig. 5.3
//   string sListFile("D:\\CompPhys\\SymplecticEuler1.dat");
               //  save data for Fig. 5.1c and Fig. 5.3
//   string sListFile("D:\\CompPhys\\SymplecticEuler2.dat");
               //  save data for Fig. 5.1d and Fig. 5.3

   CKeplerIntegrator I(dE,dDeltaT,dMaxT); //  Construct the CKeplerIntegrator class

//   I.ExplicitEuler();       //  execute the Explicit Euler Integrator
   I.ImplicitEuler();       //  execute the Implicit Euler Integrator
//   I.SymplecticEuler1();       //  execute the Symplectic Euler Integrator 1
//   I.SymplecticEuler2();       //  execute the Symplectic Euler Integrator 2
   ofLst.open(sListFile,ios::out);  // create output file
   if (!ofLst) {
      cout << "Can't open file: " << sListFile << endl;
      return 1;
   }
   for (register int i=0; i<=I.m_nMaxElem; ++i) //  write results on output file
      ofLst << fixed << setw(12) << setprecision(5) << I.m_vdT[i]
            << setw(12) << setprecision(5) << I.m_vdP1[i]
            << setw(12) << setprecision(5) << I.m_vdP2[i]
            << setw(12) << setprecision(5) << I.m_vdQ1[i]
            << setw(12) << setprecision(5) << I.m_vdQ2[i]
            << setw(12) << setprecision(5) << I.m_vdH[i] << endl;
   ofLst.close();

   return 0;
}
