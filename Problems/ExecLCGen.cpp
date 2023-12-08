//  Execute test for the linear congruential generator, Section 12.2

#include "stdafx.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <cmath>
#include "matrix2d.h"
#include "LCGenerator.h"

using namespace std;

/*****************************************************************************/
int ExecuteLCGenerator()
/*****************************************************************************/
{
   ofstream ofLst;

//   Park-Miller parameters of the 'good' generator

   //unsigned long long nA = static_cast<unsigned long long>(pow(7.0,5.0)),
   //   nC = 0,
   //   nM = (static_cast<unsigned long long>(2<<30))-1;
   //   nX0 = 281;

//   Park-Miller parameters of the 'bad' generator

   unsigned long long nA = 137,     //  Park-Miller parameter a
      nC = 0,                       //  Park-Miller parameter c
      nM = 2<<10,     //  Park-Miller parameter m for Fig. 14.1b
//      nM = 2<<9,     //  Park-Miller parameter m for Fig. 12.1b
      nX0 = 281;                    //  seed
   int nNmax = 500,  //  number of random numbers, for Fig. 14.1b
//   int nNmax = 1000,  //  number of random numbers for Figs. 12.1a, 12.1b, 14.1a
//   int nNmax = 10000,
//   int nNmax = 100000,  //  number of random numbers for Fig. 12.2a
//   int nNmax = 1000000,  //  number of random numbers for Fig. 12.2b
//   int nNmax = 10000000,  //  number of random numbers for Fig. 12.2c
      nMmax= 100;    //  number of bins for histogram test
//   string sListFile("D:\\CompPhys\\ran1.dat");
        //  data file for Fig. 12.1a
//   string sListFile("D:\\CompPhys\\ran2.dat");
        //  data file for Fig. 12.1b
//   string sListFile("D:\\CompPhys\\ran3.dat");
        //  data file for Fig. 14.1a
   string sListFile("D:\\CompPhys\\ran5.dat");
        //  data file for Fig. 14.1b
//   string sListFile("D:\\CompPhys\\ranh1.dat");
        //  data file for Fig. 12.2a
//   string sListFile("D:\\CompPhys\\ranh2.dat");
        //  data file for Fig. 12.2b
//   string sListFile("D:\\CompPhys\\ranh3.dat");
        //  data file for Fig. 12.2c
//   string sListFile("D:\\CompPhys\\ranh4.dat");
        //  data file for Fig. 12.3
   CLCGenerator L(nA,nC,nM,nX0);

//   Chapter 12, Pseudo Random Number

//   Spectral Test

   L.SpectralTest(nNmax);
   ofLst.open(sListFile,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile << endl;
      return 1;
   }
   ofLst << fixed << setprecision(6);
   for (register int i=0; i<nNmax; ++i)
      ofLst << setw(12) << L.m_vdX1[i] << setw(12) << L.m_vdX2[i] << endl;
   ofLst.close();

   //  Histogram

   //L.Histogram(nMmax,nNmax);
   //ofLst.open(sListFile,ios::out);
   //if (!ofLst) {
   //   cout << "Can't open file: " << sListFile << endl;
   //   return 1;
   //}
   //ofLst << fixed << setprecision(6);
   //for (register int i=0; i<nMmax; ++i) {
   //   ofLst << setw(12) << L.m_vdRange[i] << setw(12)
   //         << L.m_vnCount[i]*L.m_dHistScale << endl;
   //   ofLst << setw(12) << L.m_vdRange[i+1] << setw(12)
   //         << L.m_vnCount[i]*L.m_dHistScale << endl;
   //   ofLst << setw(12) << L.m_vdRange[i+1] << setw(12) << 0.0 << endl;
   //}
   //ofLst.close();

   return 0;
}
