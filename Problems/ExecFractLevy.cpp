//  Execution of the CFractLevy class, Section 17.4, fractal time Lévy flight

#include "stdafx.h"
#include "stdlib.h"
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <vector>
#include <string>
#include "FractLevy.h"

using namespace std;

/*****************************************************************************/
int ExecuteFractLevy(void)
/*****************************************************************************/
{
   ofstream ofLst;
   int nN=1000;
   double dAlpha=1.3,
      dL=0.01,
      dBeta=0.8,
      dTau=0.1;
   string sListFile("D:\\CompPhys\\ftlf1.dat");
        //  data file for Fig. 17.10, first curve
   string sListFile1("D:\\CompPhys\\ftlf2.dat");
        //  data file for Fig. 17.10, second curve
   string sListFile2("D:\\CompPhys\\ftlf3.dat");
        //  data file for Fig. 17.10, third curve
   CFractLevy F(nN,dAlpha,dL,dBeta,dTau);

   F.seed(17563);
   F.RunFractLevy();  //  run fractal time Lévy flight for Fig. 17.10, first curve
   ofLst.open(sListFile,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile << endl;
      return 1;
   }
   ofLst << scientific << setprecision(6);
   for (register int i=0; i<nN-1; ++i) {  //  generate steps
      ofLst << setw(16) << F.m_vdT[i] << setw(16) << F.m_vdX[i] << endl;
      ofLst << setw(16) << F.m_vdT[i+1] << setw(16) << F.m_vdX[i] << endl;
      ofLst << setw(16) << F.m_vdT[i+1] << setw(16) << F.m_vdX[i+1] << endl;
   }
   ofLst.close();

   F.seed(4151);
   F.RunFractLevy();  //  run fractal time Lévy flight for Fig. 17.10, second curve
   ofLst.open(sListFile1,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile1 << endl;
      return 1;
   }
   ofLst << scientific << setprecision(6);
   for (register int i=0; i<nN-1; ++i) {  //  generate steps
      ofLst << setw(16) << F.m_vdT[i] << setw(16) << F.m_vdX[i] << endl;
      ofLst << setw(16) << F.m_vdT[i+1] << setw(16) << F.m_vdX[i] << endl;
      ofLst << setw(16) << F.m_vdT[i+1] << setw(16) << F.m_vdX[i+1] << endl;
   }
   ofLst.close();

   F.seed(11679);
   F.RunFractLevy();  //  run fractal time Lévy flight for Fig. 17.10, third curve
   ofLst.open(sListFile2,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile2 << endl;
      return 1;
   }
   ofLst << scientific << setprecision(6);
   for (register int i=0; i<nN-1; ++i) {  //  generate steps
      ofLst << setw(16) << F.m_vdT[i] << setw(16) << F.m_vdX[i] << endl;
      ofLst << setw(16) << F.m_vdT[i+1] << setw(16) << F.m_vdX[i] << endl;
      ofLst << setw(16) << F.m_vdT[i+1] << setw(16) << F.m_vdX[i+1] << endl;
   }
   ofLst.close();

   return 0;
}
