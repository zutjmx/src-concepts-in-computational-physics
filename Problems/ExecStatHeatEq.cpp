#include "stdafx.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <cmath>
#include "HeatEquation.h"

using namespace std;

/*****************************************************************************/
int ExecuteStatHeatEq()
/*****************************************************************************/
{
   ofstream ofLst;
   int nN=100,        //  # of steps
      ret=0;          //  error return
   double dT0 = 0.0,  //  temperature T_0, see Fig. 9.1
      dTN = 2.0,      //  temperature T_N, see Fig. 9.1
      dLength=10.0,   //  length of the rod
      dKappa=1.0,     //  thermal diffusivity
      dTheta=-0.4,    //  max. height of the heat source/drain
      dL=1.0,         //  width of the heat source/drain
      dX=5.0;         //  position of the heat source/drain

//  Chapter 9, stationary heat equation

   CHeatEquation H(nN,dT0,dTN,dLength,dKappa);
//   string sListFile("D:\\CompPhys\\HeatN5.dat");
         //  save results for nN = 5, Fig. 9.2
//   string sListFile("D:\\CompPhys\\HeatN10.dat");
         //  save results for nN = 10, Fig. 9.3
   string sListFile("D:\\CompPhys\\HeatN100T.dat");
         //  save results for nN = 100, Fig. 9.4

   H.SetHeatSource(dTheta,dL,dX);  //  define heat source/drain, Eq. (9.21)
   if (!(ret=H.SolveEquation())) { //  solve stat. heat equation
      ofLst.open(sListFile,ios::out);

//      save data on file

      if (!ofLst) {
         cout << "Can't open file: " << sListFile << endl;
         return 1;
      }
      ofLst << fixed << setprecision(6);
      ofLst << setw(12) << H.m_vGridPoint[0] << setw(12) << dT0 << endl;
            //  save left hand boundary condition
      for (register int k=0; k<nN-1; k++)
         ofLst << setw(12) << H.m_vGridPoint[k+1] << setw(12) << H.m_vTemp[k] << endl;
      ofLst << setw(12) << H.m_vGridPoint[nN] << setw(12) << dTN << endl;
            //  save right hand boundary condition
      ofLst.close();
   } else
      cout << "Error: " << ret << endl;
   return ret;
}
