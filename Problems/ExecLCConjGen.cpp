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
   //unsigned __int64 nA = static_cast<unsigned __int64>(pow(7.0,5.0)),
   //   nC = 0,
   //   nM = static_cast<unsigned __int64>(pow(2.0,31.0))-1,
   //   nX0 = 281;
   unsigned __int64 nA = 137,
      nC = 0,
//      nM = 2<<10,
      nM = 2<<9,
      nX0 = 281;
   int nNmax = 500,
//   int nNmax = static_cast<int>(1.0e+3),
//   int nNmax = static_cast<int>(1.0e+7),
      nMmax= 100;
//   string sListFile("D:\\CompPhys\\ran1.dat");
//   string sListFile("D:\\CompPhys\\ran2.dat");
//   string sListFile("D:\\CompPhys\\ran3.dat");
//   string sListFile("D:\\CompPhys\\ran4.dat");
   string sListFile("D:\\CompPhys\\ran5.dat");
//   string sListFile("D:\\CompPhys\\ranh1.dat");
//   string sListFile("D:\\CompPhys\\ranh2.dat");
//   string sListFile("D:\\CompPhys\\ranh3.dat");
//   string sListFile("D:\\CompPhys\\ranh4.dat");
//   string sListFile("D:\\CompPhys\\ranhtest.dat");
   LCGenerator L(nA,nC,nM,nX0);

//   Chapter 12, Pseudo Random Number

//   Spectral Test

   L.SpectralTest(nNmax);
   ofLst.open(sListFile,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile;
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
   //   cout << "Can't open file: " << sListFile;
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
