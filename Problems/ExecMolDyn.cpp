//  Execute the Molecular Dynamics Problem of Chapter 7

#include "stdafx.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <cmath>
#include "MolDyn.h"

using namespace std;

/*****************************************************************************/
int ExecuteMolDyn()
/*****************************************************************************/
{
   ofstream ofLst;
   int nMax=100,                //  # of particles
      nNum=10,                  //  # of particles per row in a
                                //  square lattice
      nSteps=3000;              //  # of time steps
   double dMass=1.0,            //  particle mass
      dG=9.81,                  //  acceleration due to gravity
      dSigma = 1.0,             //  Lennard-Jones potential \sigma, Eq. (7.4)
      dEps = 1.0,               //  Lennard-Jones potential \epsilon, Eq. (7.4)
      dLength = 30.0,           //  box-length
      dH0 = 10.0,               //  initial height
      dT0 = 0.0,                //  starting time
      dT = 0.0,                 //  time
      dTau = 1.0e-3,            //  time step
      dXstart = 10.5,           //  x-position of the first particle to the left
      dDeltaX = 1.0,            //  x-distance between particels, equal dEps
      dDeltaY = 1.0;            //  y-distance between particels, equal dEps
//   string sListFile("D:\\CompPhys\\ExMolDyn0.dat");
         //  data file: initial configuration, nSteps=0, Fig. 7.3a
//   string sListFile("D:\\CompPhys\\ExMolDyn1.dat");
         //  data file: configuration after nSteps=1200 time steps, Fig. 7.3b
//   string sListFile("D:\\CompPhys\\ExMolDyn2.dat");
         //  data file: configuration after nSteps=1800 time steps, Fig. 7.3c
   string sListFile("D:\\CompPhys\\ExMolDyn3.dat");
         //  data file: configuration after nSteps=3000 time steps, Fig. 7.3d

//  Chapter 7, Molecular Dynamics;

   CMolDyn M(nMax,nNum);
   M.SetParameter(dMass,dG,dSigma,dEps,dLength);
   M.SetTimeParam(dT0,dTau);
   M.InitSqareLattice(dXstart,dH0,dDeltaX,dDeltaY);
   dT = dT0;    //  set start time
   for (register int k=1; k<=nSteps; ++k) {  //  execute the leap-frog algorithm
      dT = M.LeapFrog(dT);
      if (!(k%100))      //  just a control so you know the program is working
         cout << setw(5) << k << endl;
   }

//  save results on file

   ofLst.open(sListFile,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile << endl;
      return 1;
   }
   ofLst << fixed << setprecision(6);
   for (register int i=0; i<nMax; ++i)
      ofLst << setw(15) << M.m_vXi[i]
            << setw(15) << M.m_vYi[i]
            << setw(15) << dT    //  just for control; time elapsed
            << setw(15) << M.m_vV[i] << endl;
   ofLst.close();

   return 0;
}
