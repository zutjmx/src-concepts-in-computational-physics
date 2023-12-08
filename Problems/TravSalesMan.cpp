//  Realization the CTravSalesMan class, Section 20.3: Simulated Annealing

#include "stdafx.h"
#include "stdlib.h"
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include "TravSalesMan.h"

using namespace std;

/*****************************************************************************/
CTravSalesMan::CTravSalesMan(void)
/*****************************************************************************/
//  empty constructor, please use the property assign to define the
//  system's parameters
{
   m_nCitiesX = m_nCitiesY = m_nNumberCities = 0;
}
/*****************************************************************************/
CTravSalesMan::CTravSalesMan(int n, int m)
/*****************************************************************************/
//  input parameters:
//  n ... number of lattice points in x-directions
//  m ... number of lattice points in y-directions
{
   m_nCitiesX = n;
   m_nCitiesY = m;
   m_nNumberCities = n*m;
   m_vnX1.assign(n,0);
   m_vnY1.assign(m,0);
   m_vnX.assign(m_nNumberCities,0);
   m_vnY.assign(m_nNumberCities,0);
   m_vnInArr.assign(m_nNumberCities,0);
   m_vnRoute.assign(m_nNumberCities,0);
   for (register int i=0; i<m_nNumberCities; ++i)
       m_vnInArr[i] = i;
}
/*****************************************************************************/
CTravSalesMan::~CTravSalesMan(void)
/*****************************************************************************/
{
}
/*****************************************************************************/
void CTravSalesMan::assign(int n, int m)
/*****************************************************************************/
//  use this property to define the parameters of the traveling salesperson problem
//  if the empty constructor has been used to instantiate the CTravSalesMan class.

//  input parameters:
//  n ... number of lattice points in x-directions
//  m ... number of lattice points in y-directions
{
   m_nCitiesX = n;
   m_nCitiesY = m;
   m_nNumberCities = n*m;
   m_vnX1.assign(n,0);
   m_vnY1.assign(m,0);
   m_vnX.assign(m_nNumberCities,0);
   m_vnY.assign(m_nNumberCities,0);
   m_vnInArr.assign(m_nNumberCities,0);
   m_vnRoute.assign(m_nNumberCities,0);
   for (register int i=0; i<m_nNumberCities; ++i)
       m_vnInArr[i] = i;
}
/*****************************************************************************/
void CTravSalesMan::RegularLattice()
/*****************************************************************************/
//  arrange cities within a regular lattice as in Fig. 20.1
{
   for (register int i=0; i<m_nCitiesX; ++i)
      m_vnX1[i] = i;
   for (register int i=0; i<m_nCitiesY; ++i)
      m_vnY1[i] = i;
   for (register int i=0; i<m_nNumberCities; i+=m_nCitiesX) 
      for (register int j=0; j<m_nCitiesX; ++j)
         m_vnX[i+j] = m_vnX1[j];
   for (register int i=0, j=0; i<m_nNumberCities; i+=m_nCitiesX, ++j)
      for (register int k=0; k<m_nCitiesX; ++k)
         m_vnY[i+k] = j;
}
/*****************************************************************************/
void CTravSalesMan::ClusterCities(int offset)
/*****************************************************************************/
//  arrange cities within a cluster configuration as in Fig. 20.2
//  input parameter:
//  offset ... number of lattice constants to offset the clusters
{
   for (register int i=0; i<m_nCitiesX; ++i)
      m_vnX1[i] = i;
   for (register int i=0; i<m_nCitiesY; ++i)
      m_vnY1[i] = i;
   for (register int i=0; i<m_nNumberCities; i+=m_nCitiesX) 
      for (register int j=0; j<m_nCitiesX; ++j)
         m_vnX[i+j] = m_vnX1[j];
   for (register int i=0, j=0; i<m_nNumberCities; i+=m_nCitiesX, ++j)
      for (register int k=0; k<m_nCitiesX; ++k)
         m_vnY[i+k] = j;

   for (register int i=0; i<m_nNumberCities; ++i) {
      if (m_vnX[i] > 2)
         m_vnX[i] += offset;
      if (m_vnY[i] > 2)
         m_vnY[i] += offset;
   }
}
/*****************************************************************************/
double CTravSalesMan::CalcT0()
/*****************************************************************************/
//  calculate starting temperature according to Eq. (20.13)
{
   int nLength=m_vnX.size()-1,
      nNumInt=10000,
      r;
   double dT0=8.0,dL,dMitL=0.0,dMitL2=0.0,dVarL;
   vector<int> vnInArr;

   for (register int j=1; j<=nNumInt; ++j) {
      vnInArr = m_vnInArr;
      for (register int i=0; i<m_nNumberCities; ++i) {
         r = static_cast<int>(floor(RandInt()*vnInArr.size()));
         m_vnRoute[i] = vnInArr[r];        //  generate route
         vnInArr.erase(vnInArr.begin()+r); //  eliminate the city so that
                                           //  each city is visited only once
      }
      dL = LengthT(m_vnRoute);             //  calculate the length of the route
      dMitL += dL;                         //  accumulate the length
      dMitL2 += dL*dL;
   }
   dMitL /= nNumInt;
   dMitL2 /= nNumInt;
   dVarL = dMitL2-dMitL*dMitL;             //  variance
   return sqrt(dVarL);                     //  Eq. (20.13)
}
/*****************************************************************************/
double CTravSalesMan::LengthT(vector<int>& vnRoute)
/*****************************************************************************/
//  calculate the length of a particular route
//  input parameter;
//  vnRoute ... vector which defines the route through the cities
{
   int nLength=vnRoute.size(),
      r1,r2;
   double dZ=0.0,dDiffX,dDiffY;

   for (register int i=0; i<nLength-1; ++i) {
      r1 = vnRoute[i];
      r2 = vnRoute[i+1];
      dDiffX = m_vnX[r1]-m_vnX[r2];
      dDiffY = m_vnY[r1]-m_vnY[r2];
      dZ += sqrt(dDiffX*dDiffX+dDiffY*dDiffY);
   }
   r1 = vnRoute[0];
   r2 = vnRoute[nLength-1];
   dDiffX = m_vnX[r1]-m_vnX[r2];
   dDiffY = m_vnY[r1]-m_vnY[r2];
   dZ += sqrt(dDiffX*dDiffX+dDiffY*dDiffY);
   return dZ;
}
/*****************************************************************************/
void CTravSalesMan::Initialize()
/*****************************************************************************/
{
   int r;
   vector<int> vnInArr(m_vnInArr);  //  make a copy of the original vector
                                    //  of city numbers

   m_dT0 = CalcT0();      //  calculate the starting temperature
   for (register int i=0; i<m_nNumberCities; ++i) {  //  define a first route
      r = static_cast<int>(floor(RandInt()*vnInArr.size()));
      m_vnRoute[i] = vnInArr[r];        //  by randomly choosing cities
      vnInArr.erase(vnInArr.begin()+r); //  eliminate the city so that
                                        //  each city is visited only once
   }
}
/*****************************************************************************/
void CTravSalesMan::RunTravSalesMan(int nNum, double dEps, double dQ)
/*****************************************************************************/
//  run the simulated annealing process to solve the traveling salesperson
//  problem
//  input parameters:
//  nNum ... number of modifications
//  dEps ... accuracy
//  dQ ..... geometric cooling schedule, Eq. (20.14)
{
   int nCount=0,r,r1,
      nPower=1;  
   double dT,dMitll=0.0,dMitl2=0.0,dMitlSave,dCrit,dL,dE,dPr,dVar;
   vector<int> vnRouteTemp;

   while (nCount < 10) {          //  stop iteration after 10 consecutive
                                  //  steps with required accuracy
      dT = m_dT0*pow(dQ,nPower);  //  geometric cooling schedule, Eq. (20.14)
      dMitlSave = dMitll;
      dMitll = dMitl2 = 0.0;
      for (register int i=1; i<=nNum; ++i) {
         r = r1 = 1;
         while (r == r1) {  //  get two cities randomly
            r = static_cast<int>(floor(RandInt()*m_nNumberCities));
            r1 = static_cast<int>(floor(RandInt()*m_nNumberCities));
         }
         vnRouteTemp = m_vnRoute;
         vnRouteTemp[r] = m_vnRoute[r1]; //  interchange the two cities
         vnRouteTemp[r1] = m_vnRoute[r];
         dE = LengthT(vnRouteTemp)-LengthT(m_vnRoute); //  new length
         dE = exp(-dE/dT);
         dPr = min(dE,1.0);                 //  Eq. (20.8)
         if (dPr == 1.0 || RandInt() < dE)  //  Metropolis
            m_vnRoute = vnRouteTemp;        //  accept the new route
         dL = LengthT(m_vnRoute);           //  length of the new route
         dMitll += dL;                      //  accumulate the lengths
         dMitl2 += dL*dL;                   //  accumulate the squares of the lengths
      }
      dMitll /= nNum;                       //  avg length
      dMitl2 /= nNum;                       //  avg square of length
      dVar = dMitl2-dMitll*dMitll;          //  variance
      if (!((nPower++)%10)) {               //  increase nPower for geometric cooling
                                            //  schedule
         m_vdTemp.push_back(dT);            //  save results every tenth step
         m_vdMitl.push_back(dMitll);
         m_vdVar.push_back(dVar);
      }
      dCrit = fabs(dMitll-dMitlSave);       //  calculate error between steps
      cout << "T = " << fixed << setprecision(5) << setw(10) << dT << "    "
           << scientific << setw(15) << setprecision(5) << dCrit << endl;
      //  control output
      if (dCrit <= dEps)                    //  accuracy reached
         ++nCount;                          //  increase accuracy steps count
      else
         nCount = 0;
   }
   dL = LengthT(m_vnRoute);                 //  final length
   cout << "L = " << fixed << setprecision(5) << setw(15) << dL << endl;
   //  report final result
}
