//  Execution of the CRandWalk class, Chapter 17, Random Walk

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <vector>
#include <string>
#include "RandWalk.h"

using namespace std;

/*****************************************************************************/
int ExecuteRandWalk(void)
/*****************************************************************************/
{
   ofstream ofLst;
//   int nSteps=50;   //  number of time steps, Fig. 17.1a
//   int nSteps=100;   //  number of time steps, Fig. 17.1b
   int nSteps=1000;   //  number of time steps, Fig. 17.1c
   double dP=0.5;     //  probability p, step 1, page 245
//   string sListFile("D:\\CompPhys\\randw.dat");
        //  data file for Fig. 17.1a, first curve
//   string sListFile1("D:\\CompPhys\\randw1.dat");
        //  data file for Fig. 17.1a, second curve
//   string sListFile2("D:\\CompPhys\\randw2.dat");
        //  data file for Fig. 17.1a, third curve
//   string sListFile("D:\\CompPhys\\randw3.dat");
        //  data file for Fig. 17.1b, first curve
//   string sListFile1("D:\\CompPhys\\randw4.dat");
        //  data file for Fig. 17.1b, second curve
//   string sListFile2("D:\\CompPhys\\randw5.dat");
        //  data file for Fig. 17.1b, third curve
   string sListFile("D:\\CompPhys\\randw6.dat");
        //  data file for Fig. 17.1c, first curve
   string sListFile1("D:\\CompPhys\\randw7.dat");
        //  data file for Fig. 17.1c, second curve
   string sListFile2("D:\\CompPhys\\randw8.dat");
        //  data file for Fig. 17.1c, third curve
   CRandWalk R(nSteps,dP);

   R.seed(1731);
   R.RunRandWalk();       //  generate first curve
   ofLst.open(sListFile,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile << endl;
      return 1;
   }
   ofLst << fixed << setprecision(5);
   for (register int i=0; i<nSteps-1; ++i) {  //  generate steps
      ofLst << setw(10) << i << setw(15) << R.m_vnX[i] << endl;
      ofLst << setw(10) << (i+1) << setw(15) << R.m_vnX[i] << endl;
      ofLst << setw(10) << (i+1) << setw(15) << R.m_vnX[i+1] << endl;
   }
   ofLst.close();

   R.RunRandWalk();       //  generate second curve
   ofLst.open(sListFile1,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile1 << endl;
      return 1;
   }
   ofLst << fixed << setprecision(5);
   for (register int i=0; i<nSteps-1; ++i) {  //  generate steps
      ofLst << setw(10) << i << setw(15) << R.m_vnX[i] << endl;
      ofLst << setw(10) << (i+1) << setw(15) << R.m_vnX[i] << endl;
      ofLst << setw(10) << (i+1) << setw(15) << R.m_vnX[i+1] << endl;
   }
   ofLst.close();

   R.RunRandWalk();       //  generate third curve
   ofLst.open(sListFile2,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile2 << endl;
      return 1;
   }
   ofLst << fixed << setprecision(5);
   for (register int i=0; i<nSteps-1; ++i) {  //  generate steps
      ofLst << setw(10) << i << setw(15) << R.m_vnX[i] << endl;
      ofLst << setw(10) << (i+1) << setw(15) << R.m_vnX[i] << endl;
      ofLst << setw(10) << (i+1) << setw(15) << R.m_vnX[i+1] << endl;
   }
   ofLst.close();

   return 0;
}