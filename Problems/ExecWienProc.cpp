//  Execute the Wiener Process, Section 16.3

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <vector>
#include <string>
#include <random>
#include "matrix2d.h"

using namespace std;

/*****************************************************************************/
int ExecuteWienerProc(void)
/*****************************************************************************/
{
   ofstream ofLst;
   int nNum=1000;       //  number of time steps
   double dTstep=0.01,  //  time step
//      dMu=0.0;        //  no drift,  Fig. 16.2
      dMu=1.0;          //  with drift, Fig. 17.3

//  See for explanation:
//  Nicolai M. Josuttis: The C++ Standard Library, Second Edition
//  Addison Wesley (2012)
//  ISBN: 978-0-321-62321-8
//  Chapter 17
//
//  If this doesn't work on UNIX (LINUX) please use the CNormDist class to
//  generate random numbers which follow a normal distribution with mean zero
//  and variance one

#ifdef _MSC_VER
   typedef std::ranlux64_base_01 Myeng;
#else
   typedef std::default_random_engine Myeng;
#endif
   typedef std::normal_distribution<double> Mydist; 
   Myeng eng; 
   Mydist dist(0.0, 1.0);      //  normal distribution, mean zero and variance one
   Mydist::input_type engval = eng(); 
   Mydist::result_type distval = dist(eng);
   vector <double> vdDistr(nNum,0.0),
      vdTime(nNum,0.0);
//   string sListFile("D:\\CompPhys\\wiener.dat");
        //  data file for Fig. 16.2, first curve
//   string sListFile1("D:\\CompPhys\\wiener1.dat");
        //  data file for Fig. 16.2, second curve
//   string sListFile2("D:\\CompPhys\\wiener2.dat");
        //  data file for Fig. 16.2, third curve
   string sListFile("D:\\CompPhys\\wiener3.dat");
        //  data file for Fig. 17.3, first curve
   string sListFile1("D:\\CompPhys\\wiener4.dat");
        //  data file for Fig. 17.3, second curve
   string sListFile2("D:\\CompPhys\\wiener5.dat");
        //  data file for Fig. 17.3, third curve
   
   distval = distval;  // to quiet "unused" warnings 
   engval = engval; 
   
   dist.reset(); // discard any cached values 
   vdDistr[0] = vdTime[0] = 0.0;    //  generate first curve
   for (register int i=1; i<nNum; ++i) {
      vdDistr[i] = vdDistr[i-1]+sqrt(dTstep)*dist(eng)+dMu*dTstep; 
          //  Eq. (16.31) plus drift term Eq. (17.42).
          //  See also footnote 2 on page 225
      vdTime[i] = i*dTstep;
   }
   ofLst.open(sListFile,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile << endl;
      return 1;
   }
   ofLst << fixed << setprecision(6);
   for (register int i=0; i<nNum; ++i)
      ofLst << setw(12) << vdTime[i] << setw(12)
            << vdDistr[i] << endl;
   ofLst.close();

   vdDistr[0] = vdTime[0] = 0.0;    //  generate second curve
   for (register int i=1; i<nNum; ++i) {
      vdDistr[i] = vdDistr[i-1]+sqrt(dTstep)*dist(eng)+dMu*dTstep;
      vdTime[i] = i*dTstep;
   }
   ofLst.open(sListFile1,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile1 << endl;
      return 1;
   }
   ofLst << fixed << setprecision(6);
   for (register int i=0; i<nNum; ++i)
      ofLst << setw(12) << vdTime[i] << setw(12)
            << vdDistr[i] << endl;
   ofLst.close();

   vdDistr[0] = vdTime[0] = 0.0;    //  generate third curve
   for (register int i=1; i<nNum; ++i) {
      vdDistr[i] = vdDistr[i-1]+sqrt(dTstep)*dist(eng)+dMu*dTstep;
      vdTime[i] = i*dTstep;
   }
   ofLst.open(sListFile2,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile2 << endl;
      return 1;
   }
   ofLst << fixed << setprecision(6);
   for (register int i=0; i<nNum; ++i)
      ofLst << setw(12) << vdTime[i] << setw(12)
            << vdDistr[i] << endl;
   ofLst.close();

   return 0;
}
