//  Execution of the CLevy class, Section 17.4, Lévy flight

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <vector>
#include <string>
#include "Levy.h"

using namespace std;

/*****************************************************************************/
int ExecuteLevy(void)
/*****************************************************************************/
{
   ofstream ofLst;
//   int nNsteps=1000;  //  number of time steps for Fig. 17.7
   int nNsteps=100;  //  number of time steps for Fig. 17.8
   double dAlpha=1.3,  //  Lévy index \alpha, Eq. (17.78)
//      dL=0.001;        //  minimum flight length, Fig. 17.7
      dL=0.1;        //  minimum flight length, Fig. 17.8
   //string sListFile("D:\\CompPhys\\levya1.dat");
        //  data file for Fig. 17.7, first curve
   //string sListFile1("D:\\CompPhys\\levya2.dat");
        //  data file for Fig. 17.7, second curve
   //string sListFile2("D:\\CompPhys\\levya3.dat");
        //  data file for Fig. 17.7, third curve
   string sListFile("D:\\CompPhys\\levyb1.dat");
        //  data file for Fig. 17.8, Lévy flight
   string sListFile1("D:\\CompPhys\\levyb2.dat");
        //  data file for Fig. 17.8, random walk
   CLevy L(nNsteps,dAlpha,dL);

   //L.RunLevy();   //  run Lévy flight for Fig. 17.7, first curve
   //ofLst.open(sListFile,ios::out);
   //if (!ofLst) {
   //   cout << "Can't open file: " << sListFile << endl;
   //   return 1;
   //}
   //ofLst << fixed << setprecision(6);
   //for (register int i=0; i<nNsteps-1; ++i) {
   //   ofLst << setw(15) << L.m_vdT[i] << setw(15) << L.m_vdX[i] << endl;
   //   ofLst << setw(10) << L.m_vdT[i+1] << setw(15) << L.m_vdX[i] << endl;
   //   ofLst << setw(10) << L.m_vdT[i+1] << setw(15) << L.m_vdX[i+1] << endl;
   //}
   //ofLst.close();

//   L.RunLevy();   //  run Lévy flight for Fig. 17.7, second curve
//   ofLst.open(sListFile1,ios::out);
//   if (!ofLst) {
//      cout << "Can't open file: " << sListFile1 << endl;
//      return 1;
//   }
//   ofLst << fixed << setprecision(6);
//   for (register int i=0; i<nNsteps-1; ++i) {
//      ofLst << setw(15) << L.m_vdT[i] << setw(15) << L.m_vdX[i] << endl;
//      ofLst << setw(10) << L.m_vdT[i+1] << setw(15) << L.m_vdX[i] << endl;
//      ofLst << setw(10) << L.m_vdT[i+1] << setw(15) << L.m_vdX[i+1] << endl;
//   }
//   ofLst.close();
//
//   L.RunLevy();   //  run Lévy flight for Fig. 17.7, third curve
//   ofLst.open(sListFile2,ios::out);
//   if (!ofLst) {
//      cout << "Can't open file: " << sListFile2 << endl;
//      return 1;
//   }
//   ofLst << fixed << setprecision(6);
//   for (register int i=0; i<nNsteps-1; ++i) {
//      ofLst << setw(15) << L.m_vdT[i] << setw(15) << L.m_vdX[i] << endl;
//      ofLst << setw(10) << L.m_vdT[i+1] << setw(15) << L.m_vdX[i] << endl;
//      ofLst << setw(10) << L.m_vdT[i+1] << setw(15) << L.m_vdX[i+1] << endl;
//   }
//   ofLst.close();

//   L.seed(171);
   L.RunLevy2d();   //  run Lévy flight for Fig. 17.8
   ofLst.open(sListFile,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile << endl;
      return 1;
   }
   ofLst << fixed << setprecision(6);
   for (register int i=0; i<nNsteps-1; ++i)
      ofLst << setw(15) << L.m_vdX[i] << setw(15) << L.m_vdT[i] << endl;
   ofLst.close();

   L.RunWiener2d();   //  run random walk for Fig. 17.8
   ofLst.open(sListFile1,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile1 << endl;
      return 1;
   }
   ofLst << fixed << setprecision(6);
   for (register int i=0; i<nNsteps-1; ++i)
      ofLst << setw(15) << L.m_vdX[i] << setw(15) << L.m_vdT[i] << endl;
   ofLst.close();

   return 0;
}
