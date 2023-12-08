//  Realization of the CMolDyn class, Chapter 7

#include "stdafx.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include "MolDyn.h"

using namespace std;

/*****************************************************************************/
CMolDyn::CMolDyn(void)
/*****************************************************************************/
//  empty constructor, please use the property assign to set the
//  problem parameters parameters
{
   m_nMax = m_nNum = 0;
}
/*****************************************************************************/
CMolDyn::CMolDyn(int n_max, int num)
/*****************************************************************************/
//  input parameters:
//  n_max ... # of particles
//  num ..... # of particles per row in a square lattice
{
   m_nMax = n_max;
   m_nNum = num;
}
/*****************************************************************************/
CMolDyn::~CMolDyn(void)
/*****************************************************************************/
{
}
/*****************************************************************************/
void CMolDyn::assign(int n_max, int num)
/*****************************************************************************/
//  use this property to initialze the molecular dynamic parameters if an
//  empty constructor has been used to instantiate the CMolDyn class

//  input parameters:
//  n_max ... # of particles
//  num ..... # of particles per row in a square lattice
{
   m_nMax = n_max;
   m_nNum = num;
}
/*****************************************************************************/
void CMolDyn::SetParameter(double mass, double g, double sigma, double eps,
   double length)
/*****************************************************************************/
//  set specific parameters:
//  mass .... particle mass
//  g ....... acceleration due to gravity
//  sigma ... parameter \sigma of the Lennard-Jones potential, Eq. (7.4)
//  eps ..... parameter \epsilon of the Lennard-Jones potential, Eq. (7.4)
//  length .. box length (x-axis extension)
{
   m_dMass = mass;
   m_dG = g;
   m_dSigma = sigma;
   m_dEps = eps;
   m_dLength = length;
   m_vXi.assign(m_nMax,0.0);   //  assign required space for vectors and arrays
   m_vYi.assign(m_nMax,0.0);
   m_vVxi.assign(m_nMax,0.0);
   m_vVyi.assign(m_nMax,0.0);
   m_vFxi.assign(m_nMax,0.0);
   m_vFyi.assign(m_nMax,0.0);
   m_vV.assign(m_nMax,0.0);
   m_mFxij.assign(m_nMax,m_nMax,0.0);
   m_mFyij.assign(m_nMax,m_nMax,0.0);
}
/*****************************************************************************/
void CMolDyn::InitSqareLattice(double xstart, double ystart, double deltax,
   double deltay)
/*****************************************************************************/
//  place the particles into a square lattice
//  xstart .... x-ccordinate of the bottom left particle
//  ystart .... y-ccordinate of the bottom left particle
//  deltax .... x-axis distance between two particles
//  deltay .... y-axis distance between two particles
{
   m_dXstart = m_vXi[0] = xstart;
   m_dH0 = m_vYi[0] = ystart;
   for (register int i=1; i<m_nMax; ++i) {
      if ((i+1)%m_nNum == 1) {
         m_vXi[i] = xstart;
         m_vYi[i] = m_vYi[i-1]+deltay;
      } else {
         m_vXi[i] = m_vXi[i-1]+deltax;
         m_vYi[i] = m_vYi[i-1];
      }
      m_vV[i] = sqrt(m_vVxi[i]*m_vVxi[i]+m_vVyi[i]*m_vVyi[i]);
   }
}
/*****************************************************************************/
double CMolDyn::LeapFrog(double t)
/*****************************************************************************/
//  execute the leap-frog algorithm for time instance t
{
   double xdiff,xdiff2,ydiff,ydiff2,rij,rij2,fact,fact1;

   for (register int i=0; i<m_nMax; ++i) {
      for (register int j=0; j<m_nMax; ++j) {
         xdiff = m_vXi[i]-m_vXi[j];
         xdiff2 = xdiff*xdiff;
         ydiff = m_vYi[i]-m_vYi[j];
         ydiff2 = ydiff*ydiff;
         rij = sqrt(xdiff2+ydiff2);
         rij2 = rij*rij;
         fact = m_dEps/rij;
         fact1 = 24.0*(m_dSigma/rij2)*(2.0*pow(fact,12.0)-pow(fact,6.0));
                   // Eq. (7.7)
         if (i == j) {   //  calculate f_{ij} in Eq. (7.8)
            m_mFxij[i][j] = 0.0;
            m_mFyij[i][j] = 0.0;
         } else {
            m_mFxij[i][j] = fact1*xdiff;
            m_mFyij[i][j] = fact1*ydiff;
         }
      }
      m_vFxi[i] = m_mFxij.sum(i);  //  Eq. (7.8), x-component
      m_vFyi[i] = m_mFyij.sum(i)-m_dMass*m_dG; //  Eq. (7.8), y-component
   }
   if (t >= m_dTau) {
      for (register int j=0; j<m_nMax; ++j) {   // Eq. (7.20), second line
         m_vVxi[j] += m_vFxi[j]*m_dTau;
         m_vVyi[j] += m_vFyi[j]*m_dTau;
      }
   } else {   //  Eq. (7.20), third line
      for (register int j=0; j<m_nMax; ++j) {
         m_vVxi[j] += 0.5*m_vFxi[j]*m_dTau;
         m_vVyi[j] += 0.5*m_vFyi[j]*m_dTau;
      }
   }
   for (register int j=0; j<m_nMax; ++j) {  //  Eq. (7.20), first line
      m_vXi[j] += m_dTau*m_vVxi[j];
      m_vYi[j] += m_dTau*m_vVyi[j];
   }
   t += m_dTau;  //  set next time instance
   for (register int j=0; j<m_nMax; ++j) {  //  handle boundary conditions
                                            //  Eqs. (7.24) and (7.25)
      if (m_vXi[j] < 0.0) {   //  left-hand wall?
         m_vXi[j] = -m_vXi[j];
         m_vVxi[j] = -m_vVxi[j];
      } else if (m_vXi[j] > m_dLength) {  //  right-hand wall?
         m_vXi[j] = 2.0*m_dLength-m_vXi[j];
         m_vVxi[j] = - m_vVxi[j];
      }
      if (m_vYi[j] < 0.0) {   //  bottom?
         m_vYi[j] = -m_vYi[j];
         m_vVyi[j] = -m_vVyi[j];
      }
      m_vV[j] = sqrt(m_vVxi[j]*m_vVxi[j]+m_vVyi[j]*m_vVyi[j]);
           //  calculate velocity
   }

   return t;
}
