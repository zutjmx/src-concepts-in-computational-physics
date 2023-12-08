
#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <vector>
#include <string>

using namespace std;

/*****************************************************************************/
int ExecutePoissonProc(void)
/*****************************************************************************/
{
   ofstream ofLst;
   int nN=26;
   double r,     //  uniformly distributed random number in [0,1)
      dTau;      //  waiting time
   vector <int> vnN(nN,0);
   vector <double> vdT(nN,0.0);
   string sListFile("D:\\CompPhys\\poisson.dat");
        //  data file for Fig. 16.3, first curve
   string sListFile1("D:\\CompPhys\\poisson1.dat");
        //  data file for Fig. 16.3, second curve
   string sListFile2("D:\\CompPhys\\poisson2.dat");
        //  data file for Fig. 16.3, third curve

   for (register int i=0; i<nN; ++i)
      vnN[i] = i;
   for (register int i=1; i<nN; ++i) {
      r = static_cast<double>(rand())/(static_cast<double>(RAND_MAX)+1.0);
      dTau = -log(1.0-r);      //  Eq. (13.26), \lambda = 1.0
      vdT[i] = vdT[i-1]+dTau;  //  next time step
   }
   ofLst.open(sListFile,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile << endl;
      return 1;
   }
   ofLst << fixed << setprecision(5);
   for (register int i=0; i<nN-1; ++i) {  //  generate steps
      ofLst << setw(15) << vdT[i] << setw(10) << vnN[i] << endl;
      ofLst << setw(15) << vdT[i+1] << setw(10) << vnN[i] << endl;
      ofLst << setw(15) << vdT[i+1] << setw(10) << vnN[i+1] << endl;
   }
   ofLst.close();

   for (register int i=0; i<nN; ++i)
      vnN[i] = i;
   for (register int i=1; i<nN; ++i) {
      r = static_cast<double>(rand())/(static_cast<double>(RAND_MAX)+1.0);
      dTau = -log(1.0-r);      //  Eq. (13.26), \lambda = 1.0
      vdT[i] = vdT[i-1]+dTau;  // next time step
   }
   ofLst.open(sListFile1,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile;
      return 1;
   }
   ofLst << fixed << setprecision(5);
   for (register int i=0; i<nN-1; ++i) {  //  generate steps
      ofLst << setw(15) << vdT[i] << setw(10) << vnN[i] << endl;
      ofLst << setw(15) << vdT[i+1] << setw(10) << vnN[i] << endl;
      ofLst << setw(15) << vdT[i+1] << setw(10) << vnN[i+1] << endl;
   }
   ofLst.close();

   for (register int i=0; i<nN; ++i)
      vnN[i] = i;
   for (register int i=1; i<nN; ++i) {
      r = static_cast<double>(rand())/(static_cast<double>(RAND_MAX)+1.0);
      dTau = -log(1.0-r);      //  Eq. (13.26), \lambda = 1.0
      vdT[i] = vdT[i-1]+dTau;  //  next time step
   }
   ofLst.open(sListFile2,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile;
      return 1;
   }
   ofLst << fixed << setprecision(5);
   for (register int i=0; i<nN-1; ++i) {  //  generate steps
      ofLst << setw(15) << vdT[i] << setw(10) << vnN[i] << endl;
      ofLst << setw(15) << vdT[i+1] << setw(10) << vnN[i] << endl;
      ofLst << setw(15) << vdT[i+1] << setw(10) << vnN[i+1] << endl;
   }
   ofLst.close();

   return 0;
}
