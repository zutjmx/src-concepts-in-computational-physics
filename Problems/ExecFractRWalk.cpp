//  Execution of the CFractRandWalk class, Section 17.4, fractal random walk

#include "stdafx.h"
#include "stdlib.h"
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <vector>
#include <string>
#include "FractRandWalk.h"

using namespace std;

/*****************************************************************************/
int ExecuteFractRandWalk(void)
/*****************************************************************************/
{
   ofstream ofLst;
   int nN=1000;        //  number of time steps
   double dBeta=0.8,   //  parameter \beta, Eq. (17.81)
      dTau=0.1;        //  minimal waiting time, Eq. (17.81)
   string sListFile("D:\\CompPhys\\ftrw1.dat");
        //  data file for Fig. 17.9, first curve
   string sListFile1("D:\\CompPhys\\ftrw2.dat");
        //  data file for Fig. 17.9, second curve
   string sListFile2("D:\\CompPhys\\ftrw3.dat");
        //  data file for Fig. 17.9, third curve
   CFractRandWalk F(nN,dBeta,dTau);

   F.seed(1163);
   F.RunFractRandWalk();  //  run fractal random walk for Fig. 17.9, first curve
   ofLst.open(sListFile,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile << endl;
      return 1;
   }
   ofLst << scientific << setprecision(6);
   for (register int i=0; i<nN-1; ++i) {  //  generate steps
      ofLst << setw(15) << F.m_vdT[i] << setw(15) << F.m_vdX[i] << endl;
      ofLst << setw(15) << F.m_vdT[i+1] << setw(15) << F.m_vdX[i] << endl;
      ofLst << setw(15) << F.m_vdT[i+1] << setw(15) << F.m_vdX[i+1] << endl;
   }
   ofLst.close();

   F.seed(14151);
   F.RunFractRandWalk();  //  run fractal random walk for Fig. 17.9, second curve
   ofLst.open(sListFile1,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile1 << endl;
      return 1;
   }
   ofLst << scientific << setprecision(6);
   for (register int i=0; i<nN-1; ++i) {  //  generate steps
      ofLst << setw(15) << F.m_vdT[i] << setw(15) << F.m_vdX[i] << endl;
      ofLst << setw(10) << F.m_vdT[i+1] << setw(15) << F.m_vdX[i] << endl;
      ofLst << setw(10) << F.m_vdT[i+1] << setw(15) << F.m_vdX[i+1] << endl;
   }
   ofLst.close();

   F.seed(11679);
   F.RunFractRandWalk();  //  run fractal random walk for Fig. 17.9, third curve
   ofLst.open(sListFile2,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile2 << endl;
      return 1;
   }
   ofLst << scientific << setprecision(6);
   for (register int i=0; i<nN-1; ++i) {  //  generate steps
      ofLst << setw(15) << F.m_vdT[i] << setw(15) << F.m_vdX[i] << endl;
      ofLst << setw(10) << F.m_vdT[i+1] << setw(15) << F.m_vdX[i] << endl;
      ofLst << setw(10) << F.m_vdT[i+1] << setw(15) << F.m_vdX[i+1] << endl;
   }
   ofLst.close();

   return 0;
}
