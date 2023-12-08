//  Execution of the CStatSchrodEq class. Solve the particle in a box
//  problem. Chapter 10

#include "stdafx.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <cmath>
#include "StatSchrodEq.h"

using namespace std;

/*****************************************************************************/
int ExecuteStatSchrodEq()
/*****************************************************************************/
{
   ofstream ofLst;
   bool bFlag = false;
   int nMax = 100,          //  # of grid-points
      nMaxIter = 100,       //  max # of iterations
      k = 1;
               //  set energy \epsilon_a
//   double dEpsa = 3.0,       //  n = 1, eigenvalue 1, plain box 
//   double dEpsa = 15.0,      //  n = 2
//   double dEpsa = 40.0,      //  n = 3
//   double dEpsa = 72.0,      //  n = 4
//   double dEpsa = 120.0,      //  n = 5
//   double dEpsa = 15.0,      //  n = 1, eigenvalue 1, pox with potential v_1
                               //  according to Eq. (10.59)
//   double dEpsa = 40.0,      //  n = 2, v_1
//   double dEpsa = 80.0,      //  n = 3, v_1
//   double dEpsa = 120.0,      //  n = 4, v_1
//   double dEpsa = 170.0,      //  n = 5, v_1
//   double dEpsa = 30.0,      //  n = 1, eigenvalue 1, pox with potential v_2
                               //  according to Eq. (10.59)
//   double dEpsa = 40.0,      //  n = 2, v_2
//   double dEpsa = 65.0,      //  n = 3, v_2
//   double dEpsa = 100.0,      //  n = 4, v_2
//   double dEpsa = 145.0,      //  n = 5, v_2
//   double dEpsa = 20.0,      //  n = 1, eigenvalue 1, pox with potential v_3
                               //  according to Eq. (10.59)
//   double dEpsa = 40.0,      //  n = 2, v_3
//   double dEpsa = 70.0,      //  n = 3, v_3
//   double dEpsa = 100.0,      //  n = 4, v_3
   double dEpsa = 145.0,       //  n = 5, v_3
      dStep = 1.0,       //  Numerov: dEpsb = dEpsa+dStep,
      dSexp,
      dSexp2,
      dVarS,
      dEta = 1.0e-10;    //  accuracy
   CStatSchrodEq S(nMax,nMaxIter);
//   string sListFile("D:\\CompPhys\\StatSchr1.dat");
        //  data file for eigenvalue 1 no potential, Fig. 10.1
//   string sListFile("D:\\CompPhys\\StatSchr2.dat");
        //  data file for eigenvalue 2 no potential, Fig. 10.1
//   string sListFile("D:\\CompPhys\\StatSchr3.dat");
        //  data file for eigenvalue 3 no potential, Fig. 10.1
//   string sListFile("D:\\CompPhys\\StatSchr4.dat");
        //  data file for eigenvalue 4 no potential, Fig. 10.1
//   string sListFile("D:\\CompPhys\\StatSchr5.dat");
        //  data file for eigenvalue 5 no potential, Fig. 10.1
//   string sListFile("D:\\CompPhys\\StatSchrPot1.dat");
        //  data file, potential 1, Eq. (10.59), Fig. 10.2
//   string sListFile("D:\\CompPhys\\StatSchrPot2.dat");
        //  data file, potential 2, Eq. (10.59), Fig. 10.2
//   string sListFile("D:\\CompPhys\\StatSchrPot3.dat");
        //  data file, potential 3, Eq. (10.59), Fig. 10.2
//   string sListFile("D:\\CompPhys\\StatSchr11.dat");
        //  data file for eigenvalue 1, potential v_1, Fig. 10.3
//   string sListFile("D:\\CompPhys\\StatSchr12.dat");
        //  data file for eigenvalue 2, potential v_1, Fig. 10.3
//   string sListFile("D:\\CompPhys\\StatSchr13.dat");
        //  data file for eigenvalue 3, potential v_1, Fig. 10.3
//   string sListFile("D:\\CompPhys\\StatSchr14.dat");
        //  data file for eigenvalue 4, potential v_1, Fig. 10.3
//   string sListFile("D:\\CompPhys\\StatSchr15.dat");
        //  data file for eigenvalue 5, potential v_1, Fig. 10.3
//   string sListFile("D:\\CompPhys\\StatSchr21.dat");
        //  data file for eigenvalue 1, potential v_2, Fig. 10.4
//   string sListFile("D:\\CompPhys\\StatSchr22.dat");
        //  data file for eigenvalue 2, potential v_2, Fig. 10.4
//   string sListFile("D:\\CompPhys\\StatSchr23.dat");
        //  data file for eigenvalue 3, potential v_2, Fig. 10.4
//   string sListFile("D:\\CompPhys\\StatSchr24.dat");
        //  data file for eigenvalue 4, potential v_2, Fig. 10.4
//   string sListFile("D:\\CompPhys\\StatSchr25.dat");
        //  data file for eigenvalue 5, potential v_2, Fig. 10.4
//   string sListFile("D:\\CompPhys\\StatSchr31.dat");
        //  data file for eigenvalue 1, potential v_3, Fig. 10.5
//   string sListFile("D:\\CompPhys\\StatSchr32.dat");
        //  data file for eigenvalue 2, potential v_3, Fig. 10.5
//   string sListFile("D:\\CompPhys\\StatSchr33.dat");
        //  data file for eigenvalue 3, potential v_3, Fig. 10.5
//   string sListFile("D:\\CompPhys\\StatSchr34.dat");
        //  data file for eigenvalue 4, potential v_3, Fig. 10.5
   string sListFile("D:\\CompPhys\\StatSchr35.dat");
        //  data file for eigenvalue 1, potential v_3, Fig. 10.5

//  Chapter 10, Stationary Schroedinger Equation

//   S.Potential1();      //  define potential v_1
//   S.Potential2();      //  define potential v_2
   S.Potential3();      //  define potential v_3

//  save potential for plot, if required

      //ofLst.open(sListFile,ios::out);
      //if (!ofLst) {
      //   cout << "Can't open file: " << sListFile << endl;
      //   return 1;
      //}
      //ofLst << setprecision(4);
      //for (register int i=0; i<=nMax; ++i)
      //   ofLst << scientific << setw(15) << S.m_vSl[i]
      //         << setw(15) << S.m_vVl[i] << endl;
      //ofLst.close();
   if (S.Numerov(dEpsa,dStep,dEta)) {
      //  if required:
      dSexp = S.ExpValue(S.m_vSl);   // expectation value <x>, Eq. (10.49)
      dSexp2 = S.ExpValue(S.m_vSl2); // expectation value <x^2>, Eq. (10.50)
      dVarS = dSexp2-dSexp*dSexp;    // var(x), Eq. (10.51)

//   save results on file

      ofLst.open(sListFile,ios::out);
      if (!ofLst) {
         cout << "Can't open file: " << sListFile << endl;
         return 1;
      }
      ofLst << setprecision(4);
      for (register int i=0; i<=nMax; ++i)
         ofLst << scientific << setw(15) << S.m_vSl[i]
               << setw(15) << S.m_vPhil[i] << endl;
      ofLst.close();
   } else {
      cout << "Not converged" << endl;
      return 1;
   }

   return 0;
}
