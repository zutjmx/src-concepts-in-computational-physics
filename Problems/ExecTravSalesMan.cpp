#include "stdafx.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <cmath>
#include "TravSalesMan.h"

using namespace std;

/*****************************************************************************/
int ExecuteTravSalesMan()
/*****************************************************************************/
{
   ofstream ofLst;
   int nCitiesX=6,     //  number of lattice points in x-direction
      nCitiesY=6,      //  number of lattice points in y-direction
      nNumberCities=nCitiesX*nCitiesY,  //  total # of cities
      nNumberRoutes=5000;   //  max. number of routes
   double dQ=0.99,          //  geometric cooling schedule
          dEps=1.0e-6;
   string sListFile("D:\\CompPhys\\tspbas.dat");     //  data files for reg. grid
   string sListFile1("D:\\CompPhys\\tspstart.dat");
        //  data file for Fig. 20.1a
   string sListFile2("D:\\CompPhys\\tspend.dat");
        //  data file for Fig. 20.1b
   string sListFile3("D:\\CompPhys\\tsptemp.dat");
        //  data file for Fig. 20.3
//   string sListFile("D:\\CompPhys\\tspclbas.dat");     //  data files for reg. grid
//   string sListFile1("D:\\CompPhys\\tspclstart.dat");
        //  data file for Fig. 20.2a
//   string sListFile2("D:\\CompPhys\\tspclend.dat");
        //  data file for Fig. 20.2b
//   string sListFile3("D:\\CompPhys\\tspcltemp.dat");
        //  data file for Fig. 20.4
   CTravSalesMan T(nCitiesX,nCitiesY);

   T.seed(131);
   T.RegularLattice();    //  create regular lattice, Fig. 20.1
//   T.ClusterCities(5);    //  create city cluster, Fig. 20.2

   T.Initialize();
   cout << "T0 = " << setprecision(6) << setw(15) << T.m_dT0 << endl;
      //  just to inform you the program is working
   T.RunTravSalesMan(nNumberRoutes,dEps,dQ);

/*   ofLst.open(sListFile2,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile2 << endl;
      return 1;
   }
   ofLst << fixed << setprecision(6);
   int xx,yy;
   for (register int i=0; i<nNumberCities; ++i) {
      int j = T.m_vnRoute[i];
      if (i == 1) {
         xx = T.m_vnX[j]+1;
         yy = T.m_vnY[j]+1;
      }
      ofLst << setw(10) << T.m_vnX[j]+1 << setw(10) << T.m_vnY[j]+1 << endl;
   }
   ofLst << setw(10) << xx << setw(10) << yy << endl;
   ofLst.close();  */

/*   ofLst.open(sListFile3,ios::out);
   if (!ofLst) {
      cout << "Can't open file: " << sListFile3 << endl;
      return 1;
   }
   ofLst << fixed << setprecision(6);
   for (register int i=0; i<static_cast<int>(T.m_vdTemp.size()); ++i)
      ofLst << setw(15) << T.m_vdTemp[i] << setw(15) << T.m_vdMitl[i]
            << setw(15) << T.m_vdVar[i] << endl;
   ofLst.close();  */

   return 0;
}
