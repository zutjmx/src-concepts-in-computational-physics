//  Execution of the COrnsteinUhlenbeck class, Section 17.3

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <vector>
#include <string>
#include <random>
#include "OrnsteinUhlenbeck.h"

using namespace std;

/*****************************************************************************/
int ExecuteOrnUhl(void)
/*****************************************************************************/
{
   ofstream ofLst;
//   int nNstep=1000; //  number of time steps for Fig. 17.4
   int nNstep=100000; //  number of time steps for Fig. 17.5
   double dBeta=1.0,  //  parameter \beta of Eq. (17.57)
      dA=5.0,         //  parameter A of Eq. (17.57)
      dTstep=0.01,    //  time step
      dT=0.0;         //  starting time
   //string sListFile("D:\\CompPhys\\OrnUhl1.dat");
        //  data file for Fig. 17.4, first curve
   //string sListFile1("D:\\CompPhys\\OrnUhl2.dat");
        //  data file for Fig. 17.4, second curve
   //string sListFile2("D:\\CompPhys\\OrnUhl3.dat");
        //  data file for Fig. 17.4, third curve
   string sListFile("D:\\CompPhys\\OrnUhl4.dat");
        //  data file for Fig. 17.5, first curve
   string sListFile1("D:\\CompPhys\\OrnUhl5.dat");
        //  data file for Fig. 17.5, second curve
   string sListFile2("D:\\CompPhys\\OrnUhl6.dat");
        //  data file for Fig. 17.5, third curve
   COrnsteinUhlenbeck O(nNstep,dBeta,dA,dTstep);

   O.RunOrnUhl(5.0,0.0);    //  generate first curve, v_0 = 5
   ofLst.open(sListFile,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile << endl;
      return 1;
   }
   ofLst << fixed << setprecision(6);
   dT = 0.0;
   for (register int i=0; i<nNstep; ++i, dT+=dTstep)
      ofLst << setw(10) << dT << setw(15) << O.m_vdV[i]
            << setw(15) << O.m_vdX[i] << endl;
   ofLst.close();

   O.RunOrnUhl(0.0,0.0);    //  generate second curve, v_0 = 0
   ofLst.open(sListFile1,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile1 << endl;
      return 1;
   }
   ofLst << fixed << setprecision(6);
   dT = 0.0;
   for (register int i=0; i<nNstep; ++i, dT+=dTstep)
      ofLst << setw(10) << dT << setw(15) << O.m_vdV[i]
            << setw(15) << O.m_vdX[i] << endl;
   ofLst.close();

   O.RunOrnUhl(10.0,0.0);    //  generate second curve, v_0 = 10
   ofLst.open(sListFile2,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile2 << endl;
      return 1;
   }
   ofLst << fixed << setprecision(6);
   dT = 0.0;
   for (register int i=0; i<nNstep; ++i, dT+=dTstep)
      ofLst << setw(10) << dT << setw(15) << O.m_vdV[i]
            << setw(15) << O.m_vdX[i] << endl;
   ofLst.close();

   return 0;
}
