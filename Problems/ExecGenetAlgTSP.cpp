//  Execution of the CGeneticAlgTSP class, Section 20.4: genetic algorithms
//  traveling salesperson problem

#include "stdafx.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <cmath>
#include "matrix2d.h"
#include "GeneticAlgTSP.h"

using namespace std;

/*****************************************************************************/
int ExecuteGenetAlgTSP()
/*****************************************************************************/
{
   ofstream ofLst;
   int nCitiesX=5,     //  number of lattice points in x-direction
      nCitiesY=6,      //  number of lattice points in y-direction
      nNumberCities=nCitiesX*nCitiesY,  //  total # of cities
      nPopulations=5000,   //  number of populations, must be even!
      nGenerations=5000;   //  number of generations
   double dMutp=1.0;       //  mutation probability
   string sListFile("D:\\CompPhys\\tspbas.dat");     //  data files for reg. grid
   string sListFile1="D:\\CompPhys\\tspgastart1.dat";
        //  data file for Fig. 20.5a
   string sListFile2="D:\\CompPhys\\tspgaend.dat";
        //  data file for Fig. 20.5a
   CGeneticAlgTSP G(nCitiesX,nCitiesY,nPopulations,nGenerations,dMutp);

   G.seed(13179);
   G.RegularGrid();    //  create regular lattice, Fig. 20.5

   G.Initialize();

/*   dLength = G.GetLength(G.m_mnIndividuals[1001]);

   ofLst.open(sListFile1,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile1 << endl;
      return 1;
   }
   ofLst << fixed << setprecision(6);
   int xx,yy;
   for (register int i=0; i<nNumberCities; ++i) {
      int j = G.m_vnRoute[i];
      if (!i) {
         xx = G.m_vnX[j]+1;
         yy = G.m_vnY[j]+1;
      }
      ofLst << setw(10) << G.m_vnX[j]+1 << setw(10) << G.m_vnY[j]+1 << endl;
   }
   ofLst << setw(10) << xx << setw(10) << yy << endl;
   ofLst.close();
*/

   G.RunGeneticAlgTSP();
   cout << "Final length = " << fixed << setw(15) << setprecision(5) << G.m_dLength
        << ", # of generations = " << setw(7) << G.m_nFinGeneration << endl;
   //  display the final result

   ofLst.open(sListFile2,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile2 << endl;
      return 1;
   }
   ofLst << fixed << setprecision(6);
   int xx,yy;
   for (register int i=0; i<nNumberCities; ++i) {
      int j = G.m_vnRoute[i];
      if (!i) {
         xx = G.m_vnX[j]+1;
         yy = G.m_vnY[j]+1;
      }
      ofLst << setw(10) << G.m_vnX[j]+1 << setw(10) << G.m_vnY[j]+1 << endl;
   }
   ofLst << setw(10) << xx << setw(10) << yy << endl;
   ofLst.close();

   return 0;
}
