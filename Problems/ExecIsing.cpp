//  Execute the CIsing class, Chapter 15: Ising Model

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <vector>
#include <string>
#include "matrix2d.h"
#include "Ising.h"

using namespace std;

/*****************************************************************************/
int ExecuteIsing(void)
/*****************************************************************************/
{
   ofstream ofLst;
//   int nN=5,        //  system size
//   int nN=20,
   int nN=50,
//   int nN=100,
      nNsweep=30,     //  number of sweeps
//      nNsweep=3000,
      nTemp=50,       //  number of temperature steps
//      nTemp=1,
//      nMeasure=75,
//      nMeasure=100,
      nMeasure=150,   //  number of measurements
//      nMeasure=200,
//      nMeasure=500,
//      nMeasure=5000,
//      nMeasure=7500,
//      nMeasure=10000,
//      nMeasure=15000,
//      nTherm=500,
      nTherm=1000,    //  number of thermalization steps
      nStep=10;       //  number of sweeps between measurements
   double dKbTstart=3.0,  //  starting temperature
//   double dKbTstart=1.175,
      dKbTend=0.0,        //  final temperature
      dJ=0.5,             //  exchange parameter
      dH=0.0;             //  external field
   bool bColdWarm=true,     //  = true, warm start
//      bRandSampl=true;      //  = true, random sweep
      bRandSampl=false;      //  = false, sequential sweep
//   string sListFile("D:\\CompPhys\\isspin1.dat");
        //  data file for Fig. 15.4
//   string sListFile("D:\\CompPhys\\isspin2.dat");
        //  data file for Fig. 15.7
//   string sListFile("D:\\CompPhys\\iscpw1.dat");
        //  data file for Fig. 15.3
//   string sListFile("D:\\CompPhys\\iscpN5.dat");
        //  data file for Fig. 15.5, nN = 5
//   string sListFile("D:\\CompPhys\\iscpN20.dat");
        //  data file for Fig. 15.5, nN = 20
   string sListFile("D:\\CompPhys\\iscpN50.dat");
        //  data file for Fig. 15.5, nN = 50
//   string sListFile("D:\\CompPhys\\iscpN100.dat");
        //  data file for Fig. 15.5, nN = 100
   CIsing I(nN,dJ,dH,dKbTstart,dKbTend,nTemp);

//   I.seed(131);
//   I.RunIsingSweeps(bColdWarm,nNsweep);   //  use this command to generate the
                                            //  data file for Fig. 15.3
   I.RunIsingModel(bColdWarm,bRandSampl,nMeasure,nTherm,nStep);
   ofLst.open(sListFile,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile << endl;
      return 1;
   }

   //   generate matrix data file for Figs. 15.4 and 15.7

   //ofLst << fixed << setprecision(5);
   //for (register int i=0; i<static_cast<int>(I.m_vnSigma.size()); ++i) {
   //   ofLst << setw(4) << I.m_vnSigma[i];
   //   if (!((i+1)%nN))
   //      ofLst << endl;
   //}

   //   generate data file for Fig. 15.3

   //ofLst << fixed << setprecision(5);
   //for (register int i=0; i<I.m_nKindex; ++i) {
   //   if (!(i%100))
   //      ofLst << setw(8) << i << setw(12) << I.m_vdAvgM[i] << setw(12) << I.m_vdAvgE[i]
   //            <<  endl;
   //}

   //   generate standard data file with error information

   ofLst << scientific << setprecision(5);
   for (register int i=0; i<nTemp; ++i) {
      double dKbT = I.m_vdTemp[i],
         dErrE = sqrt(fabs(I.m_vdVarE[i])*nN),
         dErrM = sqrt(fabs(I.m_vdVarM[i])*nN),
         dCv = I.m_vdVarEn[i]*nN*nN/(dKbT*dKbT),  // Eq. (15.20)
         dErrCv = I.m_vdErrVarEn[i]*nN*nN,
         dChi = I.m_vdVarMn[i]*nN*nN/dKbT,        // Eq. (15.24)
         dErrChi = I.m_vdErrVarMn[i]*nN*nN;
      ofLst << setw(15) << dKbT                 //  temperature
            << setw(15) << fabs(I.m_vdM[i])     //  magnetization
            << setw(15) << dErrM                //  error on magnetization
            << setw(15) << I.m_vdE[i]           //  energy
            << setw(15) << dErrE                //  error on energy
            << setw(15) << dCv                  //  heat capacity
            << setw(15) << dErrCv               //  error on heat capacity
            << setw(15) << dChi                 //  susceptibility
            << setw(15) << dErrChi << endl;     //  error on susceptibility
   }
   ofLst.close();
   return 0;
}
