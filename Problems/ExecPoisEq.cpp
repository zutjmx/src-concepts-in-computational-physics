//  Execution of the CPoissonEq class, Section 11.2

#include "stdafx.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <cmath>
#include "matrix2d.h"
#include "PoissonEq.h"

using namespace std;

/*****************************************************************************/
int ExecutePoissonEq()
/*****************************************************************************/
{
   ofstream ofLst;
   int nMaxX = 100,    //  number of grid points in x-direction
      nMaxY=100;       //  number of grid points in y-direction
   double dLx=10.0,    //  length in x-direction
      dLy=10.0,        //  length in y-direction
      dEta=1.0e-6,     //  accuracy
      dVx=0.0,         //  boundary condition x-direction
      dVy=0.0;         //  boundary condition y-direction
//   string sListFile("D:\\CompPhys\\rho1.dat");
        //  data file for potential, Fig. 11.1a
//   string sListFile("D:\\CompPhys\\rho2.dat");
        //  data file for potential, Fig. 11.1b
//   string sListFile("D:\\CompPhys\\rho3.dat");
        //  data file for potential, Fig. 11.1c
//   string sListFile("D:\\CompPhys\\phi1.dat");
        //  data file for solution, Fig. 11.2a
//   string sListFile("D:\\CompPhys\\phi2.dat");
        //  data file for solution, Fig. 11.2b
   string sListFile("D:\\CompPhys\\phi3.dat");
        //  data file for solution, Fig. 11.2c

//  Chapter 11,  Poisson Equation

   CPoissonEq P(nMaxX,nMaxY,dLx,dLy);
//   P.Rho(1);   //  charge density \rho_1, Eq. (11.14a)
//   P.Rho(2);   //  charge density \rho_1, Eq. (11.14b)
   P.Rho(3);   //  charge density \rho_1, Eq. (11.14c)
   P.SetDirichlet(0.0,dVx,0.0,dVy);  //  set Dirichlet boundary conditions, Eqs. (11.9)
   P.SolveEquation(dEta);            //  solve Poisson equation

//  dave results on file

   ofLst.open(sListFile,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile << endl;
      return 1;
   }
   ofLst << fixed << setprecision(4);
   for (register int i=0; i<nMaxX; ++i) {
      for (register int j=0; j<nMaxY; ++j)
         ofLst << setw(12) << P.m_mPhi[i][j];  //  store solutions
//         ofLst << setw(12) << P.m_mRho[i][j];  //  store charge densities
      ofLst << endl;
   }
	return 0;
}
