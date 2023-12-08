//  Execute the CPottsModel class, Chapter 15: Ising Model

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <vector>
#include <string>
#include "matrix2d.h"
#include "PottsModel.h"

using namespace std;

/*****************************************************************************/
int ExecutePottsModel(void)
/*****************************************************************************/
{
   ofstream ofLst;
   int nN=40,            //  system size
      nTemp=50,          //  number of temperatures
//      nMeasure=2000,     //  number of measurements
      nMeasure=10000,
      nTherm=1000,       
      nSteps=10,
//      nQ=2;            //  number of q-states
//      nQ=3;
//      nQ=4;
//      nQ=5;
//      nQ=6;
//      nQ=7;
      nQ=8;
   double dJ=0.5,        //  exchange parameter
      dH=0.0,            //  external field
//      dKbTstart=1.0,   //  starting temperature
//      dKbTend=0.0;     //  final temperature
//      dKbTstart=0.47,  //  starting temperature, Fig. 18.6
//      dKbTend=0.47;    //  final temperature, Fig. 18.6
//      dKbTstart = 0.56,  //  starting temperature, Fig. 18.7a
//      dKbTend = 0.56;    //  final temperature, Fig. 18.7a
      dKbTstart=0.372,  //  starting temperature, Fig. 18.7b
      dKbTend=0.372;    //  starting temperature, Fig. 18.7b
//   bool bRandSweep=false;  //  sequential sweep
   bool bRandSweep=true;     //  random sweep
//   string sListFile("D:\\CompPhys\\pmspin.dat");
        //  data file for Fig. 18.6
//   string sListFile("D:\\CompPhys\\PottsQ2.dat");
        //  data file for Figs. 18.2, 18.3, 18.4, and 18.5; q = 2
//   string sListFile("D:\\CompPhys\\PottsQ3.dat");
        //  data file for Figs. 18.2, 18.3, 18.4, and 18.5; q = 3
//   string sListFile("D:\\CompPhys\\PottsQ4.dat");
        //  data file for Figs. 18.2, 18.3, 18.4, and 18.5; q = 4
//   string sListFile("D:\\CompPhys\\PottsQ5.dat");
        //  data file for Figs. 18.2, 18.3, 18.4, and 18.5; q = 5
//   string sListFile("D:\\CompPhys\\PottsQ6.dat");
        //  data file for Figs. 18.2, 18.3, 18.4, and 18.5; q = 6
//   string sListFile("D:\\CompPhys\\PottsQ7.dat");
        //  data file for Figs. 18.2, 18.3, 18.4, and 18.5; q = 7
//   string sListFile("D:\\CompPhys\\PottsQ8.dat");
        //  data file for Figs. 18.2, 18.3, 18.4, and 18.5; q = 8
//   string sListFile("D:\\CompPhys\\Enpottsq2.dat");
        //  data file for Figs. 18.7a
   string sListFile("D:\\CompPhys\\Enpottsq8.dat");
        //  data file for Figs. 18.7b
   CPottsModel P(nN,dJ,dH,dKbTstart,dKbTend,nTemp);
   
   P.RunPottsModel(bRandSweep,nMeasure,nTherm,nSteps,nQ);
   ofLst.open(sListFile,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile << endl;
      return 1;
   }

   //   generate 2D matrix for Fig. 18.6

   //for (register int i=0; i<static_cast<int>(P.m_vnSigma.size()); ++i) {
   //   ofLst << setw(4) << P.m_vnSigma[i];
   //   if (!((i+1)%nN))
   //      ofLst << endl;
   //}
//   ofLst << scientific << setprecision(5);
   //for (register int i=0; i<nTemp; ++i) {
   //   double dKbT = P.m_vdTemp[i],
   //      dErrE = sqrt(fabs(I.m_vdVarE[i])*nN),
   //      dErrM = sqrt(fabs(I.m_vdVarM[i])*nN),
   //      dCv = I.m_vdVarEn[i]*nN*nN/(dKbT*dKbT),  // Eq. (15.20)
   //      dErrCv = I.m_vdErrVarEn[i]*nN*nN,
   //      dChi = I.m_vdVarMn[i]*nN*nN/dKbT,        // Eq. (15.24)
   //      dErrChi = I.m_vdErrVarMn[i]*nN*nN;
   //   ofLst << setw(15) << dKbT                 //  temperature
   //         << setw(15) << fabs(P.m_vdM[i])     //  magnetization
   //         << setw(15) << dErrM                //  error on magnetization
   //         << setw(15) << P.m_vdE[i]           //  energy
   //         << setw(15) << dErrE                //  error on energy
   //         << setw(15) << dCv                  //  heat capacity
   //         << setw(15) << dErrCv               //  error on heat capacity
   //         << setw(15) << dChi                 //  susceptibility
   //         << setw(15) << dErrChi << endl;     //  error on susceptibility
   //}

   //  generate data file for Figs. 18.7

   ofLst << scientific << setprecision(5);
   for (register int i=0; i<nMeasure; ++i)
      ofLst << setw(15) << P.m_vdEn[i] << endl;
   ofLst.close();

   return 0;
}
