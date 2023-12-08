//  Realization of the CGeneticAlgTSP class, Section 20.4: genetic algorithms
//  traveling salesperson problem

#include "stdafx.h"
#include "stdlib.h"
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <vector>
#include "matrix2d.h"
#include "GeneticAlgTSP.h"

using namespace std;

/*****************************************************************************/
CGeneticAlgTSP::CGeneticAlgTSP(void)
/*****************************************************************************/
//  empty constructor, please use the property assign to define the
//  system's parameters
{
   m_nCitiesX = m_nCitiesY = m_nNumberCities = 0;
}
/*****************************************************************************/
CGeneticAlgTSP::CGeneticAlgTSP(int n, int m, int p, int g, double mut)
/*****************************************************************************/
//  input parameters:
//  n ..... number of lattice points in x-directions
//  m ..... number of lattice points in y-directions
//  p ..... number of populations
//  g ..... number of generations
//  mut ... mutation probability
{
   m_nCitiesX = n;
   m_nCitiesY = m;
   m_nNumberCities = n*m;
   m_nPopulations = p;
   m_nGenerations = g;
   m_dMutp = mut;
   m_vnX1.assign(n,0);
   m_vnY1.assign(m,0);
   m_vnRoute.assign(m_nNumberCities,0);
   m_vnX.assign(m_nNumberCities,0);
   m_vnY.assign(m_nNumberCities,0);
   m_mnIndividuals.assign(m_nPopulations,m_nNumberCities,0);
}
/*****************************************************************************/
CGeneticAlgTSP::~CGeneticAlgTSP(void)
/*****************************************************************************/
{
}
/*****************************************************************************/
void CGeneticAlgTSP::assign(int n, int m, int p, int g, double mut)
/*****************************************************************************/
//  use this property to define the parameters of the traveling salesperson problem
//  if the empty constructor has been used to instantiate the CGeneticAlgTSP class.

//  input parameters:
//  n ..... number of lattice points in x-directions
//  m ..... number of lattice points in y-directions
//  p ..... number of populations
//  g ..... number of generations
//  mut ... mutation probability
{
   m_nCitiesX = n;
   m_nCitiesY = m;
   m_nNumberCities = n*m;
   m_nPopulations = p;
   m_nGenerations = g;
   m_dMutp = mut;
   m_vnX1.assign(n,0);
   m_vnY1.assign(m,0);
   m_vnRoute.assign(m_nNumberCities,0);
   m_vnX.assign(m_nNumberCities,0);
   m_vnY.assign(m_nNumberCities,0);
   m_mnIndividuals.assign(m_nPopulations,m_nNumberCities,0);
}
/*****************************************************************************/
void CGeneticAlgTSP::RegularGrid()
/*****************************************************************************/
//  arrange cities within a regular lattice as in Fig. 20.5
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
void CGeneticAlgTSP::Initialize()
/*****************************************************************************/
//  initialize the various individuals of all populations with a certain
//  (random) route.
{
   int r1;

   for (register int i=0; i<m_nPopulations; ++i)
      m_mnIndividuals[i][0] = 0;
   for (register int i=0; i<m_nPopulations; ++i) {
      for (register int j=1; j<m_nNumberCities; ++j) {
         r1 = static_cast<int>(floor(RandInt()*(j+1)));
         m_mnIndividuals[i][j] = r1;
      }
   }
}
/*****************************************************************************/
vector<int>& CGeneticAlgTSP::FlipLR(int* n, int length)
/*****************************************************************************/
//  flip a vector left to right
//  input parameters:
//  n ......... pointer to a one-dimensional array
//  length .... number of elements in the array
//  on return:
//  the vector containing the flipped elements is returned
{
   static vector<int> ind(length,0);

   for (register int i=length-1, k=0; i>=0; --i, ++k)
      ind[k] = n[i];
   return ind;
}
/*****************************************************************************/
void CGeneticAlgTSP::GetTour(vector<int>& ind)
/*****************************************************************************/
//  create a tour in vector m_vnRoute which touches every city only once
{
   int size=ind.size();
   vector<int>vHelp(size,0);

   for (register int i=0; i<size; ++i)
      vHelp[i] = i;
   for (register int i=0; i<size; ++i) {
      int r = ind[i];
      m_vnRoute[i] = vHelp[r];
      vHelp.erase(vHelp.begin()+r);
   }
}
/*****************************************************************************/
double CGeneticAlgTSP::LengthT(vector<int>& vnR)
/*****************************************************************************/
//  calculate the length of a particular route
//  input parameter;
//  vnR ... vector which defines the route through the cities
{
   int nLength=m_vnX.size(),r1,r2;
   double dZ=0.0,dDiffX,dDiffY;

   for (register int i=0; i<nLength-1; ++i) {
      r1 = vnR[i];
      r2 = vnR[i+1];
      dDiffX = m_vnX[r1]-m_vnX[r2];
      dDiffY = m_vnY[r1]-m_vnY[r2];
      dZ += sqrt(dDiffX*dDiffX+dDiffY*dDiffY);
   }
   r1 = vnR[0];
   r2 = vnR[nLength-1];
   dDiffX = m_vnX[r1]-m_vnX[r2];
   dDiffY = m_vnY[r1]-m_vnY[r2];
   dZ += sqrt(dDiffX*dDiffX+dDiffY*dDiffY);
   return dZ;
}
/*****************************************************************************/
double CGeneticAlgTSP::GetLength(int* n)
/*****************************************************************************/
//  calculate the length of a particular tour
//  input parameter:
//  n .... array with city numbers of the tour
//  on return the length of the tour is provided
{
   double dLength=0.0;
   int size = m_vnX.size();
   vector<int> ind(size,0);

   ind = FlipLR(n,size);
   GetTour(ind);
   dLength = LengthT(m_vnRoute);
   return dLength;
}
/*****************************************************************************/
void CGeneticAlgTSP::SortLength(vector<double>& vdLength, vector<int>& vnIndex)
/*****************************************************************************/
//  sort individuals according to their tour length
//  input parameters:
//  vdLength ... vector with tour lengths
//  vdIndex .... vector with corresponding city indices
//  on return the vectors vdLength and vnIndex contain the lengths
//  in ascending order with their corresponding city indices
{
   int nSize=vdLength.size(),nb;
   register int i,j;
   double a;

   for (j=1; j<nSize; ++j) {
      a = vdLength[j];
      nb = vnIndex[j];
      i = j;
      while (i > 0 && vdLength[i-1] > a) {
	 vdLength[i] = vdLength[i-1];
         vnIndex[i] = vnIndex[i-1];
         --i;
      }
      vdLength[i] = a;
      vnIndex[i] = nb;
   }
}
/*****************************************************************************/
void CGeneticAlgTSP::RunGeneticAlgTSP()
/*****************************************************************************/
//  Execute the traveling salesperson problem using a genetic algorithm
//  according to page 297
{
   int n=0,r1,r2,r3,r4,nIndex;
   double r;
   vector<int> vnAux(m_nPopulations,0),
      vnAuxL(2*m_nPopulations,0);
   vector <double> vdLength(2*m_nPopulations,0.0);
   matrix2d <int> mnIndividuals2(2*m_nPopulations,m_nNumberCities,0);

   for (register int k=1; k<=m_nGenerations; ++k) {

//  mutation: for each individual we introduce a single random local modification

      for (register int i=1; i<m_nPopulations; ++i) {
         r = RandInt();
         if (r < m_dMutp) {
            r1 = static_cast<int>(floor(RandInt()*(m_nNumberCities-1)))+1;
            r2 = static_cast<int>(floor(RandInt()*(r1+1)));   //  lattice position
            m_mnIndividuals[i][r1] = r2;
         }
      }

//  reproduction:  additional individuals are produced

       for (register int i=0; i<m_nPopulations; ++i) {
          memcpy(mnIndividuals2[i],m_mnIndividuals[i],m_nNumberCities*sizeof(int));
          vnAux[i] = i;
       }
       while (vnAux.size() > 0) {   //  point (a)
          nIndex = static_cast<int>(RandInt()*(m_nNumberCities-1));
          r1 = static_cast<int>(RandInt()*vnAux.size());
          r2 = vnAux[r1];
          vnAux.erase(vnAux.begin()+r1);
          r3 = static_cast<int>(RandInt()*vnAux.size());
          r4 = vnAux[r3];
          vnAux.erase(vnAux.begin()+r3);

          for (register int i=0; i<=nIndex; ++i) { //  point (b)
             mnIndividuals2[m_nPopulations+r2][i] = m_mnIndividuals[r2][i];
             mnIndividuals2[m_nPopulations+r4][i] = m_mnIndividuals[r4][i];
          }
          for (register int i=nIndex+1; i<m_nNumberCities; ++i) {
             mnIndividuals2[m_nPopulations+r2][i] = m_mnIndividuals[r4][i];
             mnIndividuals2[m_nPopulations+r4][i] = m_mnIndividuals[r2][i];
          }
       }
       vnAux.resize(m_nPopulations,0);

//  selection: keep the individuals with the highest fitness

       for (register int i=0; i<2*m_nPopulations; ++i) {
          vdLength[i] = GetLength(mnIndividuals2[i]);
          vnAuxL[i] = i;
       }
       SortLength(vdLength,vnAuxL);
      
       for (register int i=0; i<m_nPopulations; ++i) { //  copy to working space
          r1 = vnAuxL[i];
          memcpy(m_mnIndividuals[i],mnIndividuals2[r1],m_nNumberCities*sizeof(int));
       }
       if (!(++n%100))
          cout << "Length = " << fixed << setw(15) << setprecision(5) << vdLength[0]
               << endl;
          //  just so you know that the program is working
       if (fabs(vdLength[0]-m_nNumberCities)<=1.0e-6 ||   //  accuracy achieved?
           fabs(vdLength[0]-(sqrt(2.0)-1.0+m_nNumberCities))<=1.0e-6)
          break;   //  quit
   }   //  end loop over generations
   m_dLength = GetLength(m_mnIndividuals[0]);
   m_nFinGeneration = n;
}
