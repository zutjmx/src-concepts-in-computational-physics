//  Realization the CPottsModel class, Chapter 18: Potts Model

#include "stdafx.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include "matrix2d.h"
#include "PottsModel.h"

using namespace std;

/*****************************************************************************/
CPottsModel::CPottsModel(void)
/*****************************************************************************/
//  empty constructor, please use the property assign to define the
//  system's parameters
{
}
/*****************************************************************************/
CPottsModel::CPottsModel(int n, double J, double h, double ts, double te,
                         int nTsteps)
/*****************************************************************************/
//  input parameter:
//  n ........... system size (n*n)
//  J ........... exchange parameter
//  h ........... external field
//  ts .......... starting temperature
//  te .......... final temperature
//  nTsteps ..... number of temperature steps
{
   m_nN = n;
   m_nN2 = m_nN*m_nN;
   m_dKbTstart = ts;
   m_dKbTend = te;
   m_nTsteps = nTsteps;
   m_dKbTstep = (m_dKbTstart-m_dKbTend)/nTsteps;
   m_dJ = J;
   m_dH = h;
   m_mnNeighbors.assign(m_nN2,4,0);
}
/*****************************************************************************/
CPottsModel::~CPottsModel(void)
/*****************************************************************************/
{
}
/*****************************************************************************/
void CPottsModel::assign(int n, double J, double h, double ts, double te,
                         int nTsteps)
/*****************************************************************************/
//  use this property to define the parameters of the Potts model if the empty
//  constructor has been used to instantiate the CPottsModel class.

//  input parameter:
//  n ........... system size (n*n)
//  J ........... exchange parameter
//  h ........... external field
//  ts .......... starting temperature
//  te .......... final temperature
//  nTsteps ..... number of temperature steps
{
   m_nN = n;
   m_nN2 = m_nN*m_nN;
   m_dKbTstart = ts;
   m_dKbTend = te;
   m_nTsteps = nTsteps;
   m_dKbTstep = (m_dKbTstart-m_dKbTend)/nTsteps;
   m_dJ = J;
   m_dH = h;
   m_mnNeighbors.assign(m_nN2,4,0);
}
/*****************************************************************************/
void CPottsModel::GetNeighbor(void)
/*****************************************************************************/
//  define the array neighbors following the recipe on page 210 
{
   for (register int i=1, j=0; i<=m_nN2; ++i, ++j) {
      if (i+m_nN <= m_nN2)            //  up
         m_mnNeighbors[j][0] = i+m_nN;
      else
         m_mnNeighbors[j][0] = i-m_nN*(m_nN-1);
      if (i%m_nN)                  //  right
         m_mnNeighbors[j][1] = i+1;
      else
         m_mnNeighbors[j][1] = i-m_nN+1;
      if (i-m_nN >= 1)                //  down
         m_mnNeighbors[j][2] = i-m_nN;
      else
         m_mnNeighbors[j][2] = i+m_nN*(m_nN-1);
      if ((i-1)%m_nN)              //  left
         m_mnNeighbors[j][3] = i-1;
      else
         m_mnNeighbors[j][3] = i+m_nN-1;
   }
}
/*****************************************************************************/
void CPottsModel::HotStart(void)
/*****************************************************************************/
//  execute hot start
{
   for (register int i=0; i<m_nN2; ++i)
      m_vnSigma[i] = static_cast<int>(m_nQ*RandInt())+1;
}
/*****************************************************************************/
void CPottsModel::SweepR(double dKbT)
/*****************************************************************************/
//  perform a random sweep of the lattice points according to the
//  recipe (3) execution of the code, page 212
//  input parameter:
//  dKbT ... actual temperature
{
   int nIndex,nDelta11,nDelta12,nDelta13,nDelta14,nSigAux,nDelta21,nDelta22,
      nDelta23,nDelta24,nNext[4];
   double dDeltaE,dZ,dPr;

   for (register int iter=0; iter<m_nN2; ++iter) {
      nIndex = static_cast<int>(floor(m_nN2*RandInt()));
      for (register int i=0; i<4; ++i)
         nNext[i] = m_mnNeighbors[nIndex][i];
      nDelta11 = Delta(m_vnSigma[nIndex],m_vnSigma[nNext[0]-1]);
      nDelta12 = Delta(m_vnSigma[nIndex],m_vnSigma[nNext[1]-1]);
      nDelta13 = Delta(m_vnSigma[nIndex],m_vnSigma[nNext[2]-1]);
      nDelta14 = Delta(m_vnSigma[nIndex],m_vnSigma[nNext[3]-1]);
      nSigAux = static_cast<int>((m_nQ)*RandInt())+1;
      if (m_vnSigma[nIndex] <= nSigAux)
         ++nSigAux;
      nDelta21 = Delta(nSigAux,m_vnSigma[nNext[0]-1]);
      nDelta22 = Delta(nSigAux,m_vnSigma[nNext[1]-1]);
      nDelta23 = Delta(nSigAux,m_vnSigma[nNext[2]-1]);
      nDelta24 = Delta(nSigAux,m_vnSigma[nNext[3]-1]);
	   dDeltaE = m_dJ*(nDelta11+nDelta12+nDelta13+nDelta14) //  energy difference
                -m_dJ*(nDelta21+nDelta22+nDelta23+nDelta24);
      dZ = exp(-dDeltaE/dKbT);                //  \exp{-\frac{\Delta e_{ij}}{k_BT}
      dPr = min(dZ,1.0);                      //  Eq. (15.54)
      if (dPr == 1.0 || RandInt() <= dPr)     //  Metropolis
         m_vnSigma[nIndex] = nSigAux;         //  accept new configuration
   }
}
/*****************************************************************************/
void CPottsModel::SweepS(double dKbT)
/*****************************************************************************/
//  perform a sequential sweep of the lattice points according to the
//  recipe (3) execution of the code, page 212
//  input parameter:
//  dKbT ... actual temperature
{
   int nDelta11,nDelta12,nDelta13,nDelta14,nSigAux,nDelta21,nDelta22,
      nDelta23,nDelta24,nNext[4];
   double dDeltaE,dZ,dPr;

   for (register int iter=0; iter<m_nN2; ++iter) {
      for (register int i=0; i<4; ++i)
         nNext[i] = m_mnNeighbors[iter][i];
      nDelta11 = Delta(m_vnSigma[iter],m_vnSigma[nNext[0]-1]);
      nDelta12 = Delta(m_vnSigma[iter],m_vnSigma[nNext[1]-1]);
      nDelta13 = Delta(m_vnSigma[iter],m_vnSigma[nNext[2]-1]);
      nDelta14 = Delta(m_vnSigma[iter],m_vnSigma[nNext[3]-1]);
      nSigAux = static_cast<int>((m_nQ)*RandInt())+1;
      if (m_vnSigma[iter] <= nSigAux)
         ++nSigAux;
      nDelta21 = Delta(nSigAux,m_vnSigma[nNext[0]-1]);
      nDelta22 = Delta(nSigAux,m_vnSigma[nNext[1]-1]);
      nDelta23 = Delta(nSigAux,m_vnSigma[nNext[2]-1]);
      nDelta24 = Delta(nSigAux,m_vnSigma[nNext[3]-1]);
	   dDeltaE = m_dJ*(nDelta11+nDelta12+nDelta13+nDelta14) //  energy difference
                -m_dJ*(nDelta21+nDelta22+nDelta23+nDelta24);
      dZ = exp(-dDeltaE/dKbT);                //  \exp{-\frac{\Delta e_{ij}}{k_BT}
      dPr = min(dZ,1.0);                      //  Eq. (15.54)
      if (dPr == 1.0 || RandInt() <= dPr)     //  Metropolis
         m_vnSigma[iter] = nSigAux;           //  accept new configuration
   }
}
/*****************************************************************************/
void CPottsModel::ErrVar()
/*****************************************************************************/
//  error analysis using the statistical bootstrap algorithm,
//  Section 19.2
{
   int nMeasure=static_cast<int>(m_vdMn.size()), //  # of measurements
      nSample=100,  // number of data points to condider
      nIndex;
   double dEn,
      dMn;
   vector <double> vdMeanEn(nSample,0.0),
      vdMeanMn(nSample,0.0),
      vdMeanEn2(nSample,0.0),
      vdMeanMn2(nSample,0.0),
      vdVarEn(nSample,0.0),
      vdVarMn(nSample,0.0);

   for (register int i=0; i<nSample; ++i) {  //  generate a sample
      for (register int j=0; j<nMeasure; ++j) { //  choose data randomly
         nIndex = static_cast<int>((nMeasure)*RandInt());
         vdMeanEn[i] += (dEn=m_vdEn[nIndex]);
         vdMeanEn2[i] += dEn*dEn;
         vdMeanMn[i] += (dMn=m_vdMn[nIndex]);
         vdMeanMn2[i] += dMn*dMn;
      }
      vdMeanEn[i] /= nMeasure;   //  normalize data
      vdMeanEn2[i] /= nMeasure;
      vdMeanMn[i] /= nMeasure;
      vdMeanMn2[i] /= nMeasure;
      vdVarEn[i] = vdMeanEn2[i]-vdMeanEn[i]*vdMeanEn[i];  //  calculate variances
      vdVarMn[i] = vdMeanMn2[i]-vdMeanMn[i]*vdMeanMn[i];
   }
   m_dVarEnEnd = m_dVarMnEnd = 0.0;
   for (register int i=0; i<nSample; ++i) {
      m_dVarEnEnd += vdVarEn[i];
      m_dVarMnEnd += vdVarMn[i];
   }
   m_dVarEnEnd /= nSample;  //  averages over # of samples
   m_dVarMnEnd /= nSample;
   m_dErrVarEn = m_dErrVarMn = 0.0;
   for (register int i=0; i<nSample; ++i) {
      dEn = vdVarEn[i];
      m_dErrVarEn += dEn*dEn;
      dMn = vdVarMn[i];
      m_dErrVarMn += dMn*dMn;
   }
   m_dErrVarEn -= m_dVarEnEnd*m_dVarEnEnd; //  calculate remaining errors
   m_dErrVarEn = sqrt(m_dErrVarEn)/nSample;
   m_dErrVarMn -= m_dVarMnEnd*m_dVarMnEnd;
   m_dErrVarMn = sqrt(m_dErrVarMn)/nSample;
}
/*****************************************************************************/
void CPottsModel::RunPottsModel(bool rs, int nMeasure, int nTherm, int nStep,
                                int q)
/*****************************************************************************/
//  run the Potts model calculation to generate data for Figs. 18.2, 18.3, 18.4,
//  18.5, and 18.7
//  input parameters:
//  rs ..........  = true random sweep/ = false sequential sweep
//  nMeasure ....  number of measurements
//  nTherm ......  number of thermalization steps
//  nStep .......  number of steps between measurements
//  q ...........  number of q-states
{
   int nNsweep=nTherm+(nMeasure-1)*nStep,k,l=0,nDeltaVec1,nDeltaVec2,
      nDeltaVec3,nDeltaVec4,nMag,nE;
   double dKbT=m_dKbTstart,
      dMag,
      dMeanMag,
      dMeanMag2,
      dEn,
      dMeanEn,
      dMeanEn2;
   vector<int> vnSigma1(m_nN2,0);
   void (CPottsModel::*Sw)(double) = rs ? &CPottsModel::SweepR : &CPottsModel::SweepS;

   GetNeighbor();    //  get neighbor configuration
   m_vnSigma.assign(m_nN2,1);     //  initial configuration (cold start)
   m_nQ = q;
   m_vdTemp.assign(m_nTsteps,0.0);     //  prepare vectors
   m_vdMn.assign(nMeasure,0.0);
   m_vdVarM.assign(m_nTsteps,0.0);
   m_vdVarMn.assign(m_nTsteps,0.0);
   m_vdEn.assign(nMeasure,0.0);
   m_vdVarEn.assign(m_nTsteps,0.0);
   m_vdM.assign(m_nTsteps,0.0);
   m_vdE.assign(m_nTsteps,0.0);
   m_vdVarE.assign(m_nTsteps,0.0);
   m_vdErrVarMn.assign(m_nTsteps,0.0);
   m_vdErrVarEn.assign(m_nTsteps,0.0);

   HotStart();
   while (dKbT >= m_dKbTend) {
      m_vdTemp[l] = dKbT;
      k = 0;
      dMag = dMeanMag = dMeanMag2 = dMeanEn = dMeanEn2 = 0.0;
      for (register int sweep=1; sweep<=nNsweep; ++sweep) {
         (this->*Sw)(dKbT);
         if (sweep>=nTherm && !(sweep%nStep)) {
            nE = nMag = 0;
            for (register int i=0; i<m_nN2; ++i) {
               vnSigma1[i] = m_vnSigma[i];
               if (vnSigma1[i] != 1)  //  calculate for spin configuration Q = 1
                  vnSigma1[i] = 0;    //  see footnote on page 268
               nDeltaVec1 = Delta(m_vnSigma[i],m_vnSigma[m_mnNeighbors[i][0]-1]);
               nDeltaVec2 = Delta(m_vnSigma[i],m_vnSigma[m_mnNeighbors[i][1]-1]);
               nDeltaVec3 = Delta(m_vnSigma[i],m_vnSigma[m_mnNeighbors[i][2]-1]);
               nDeltaVec4 = Delta(m_vnSigma[i],m_vnSigma[m_mnNeighbors[i][3]-1]);
               nMag += vnSigma1[i];            //  magnetization
               nE += nDeltaVec1+nDeltaVec2+nDeltaVec3+nDeltaVec4;
            }
            m_vdMn[k] = dMag = static_cast<double>(nMag)/m_nN2;
                                         //  magnetization per particle
            dMeanMag += dMag;            //  accumulate mean magn. per partcicle
            dMeanMag2 += dMag*dMag;      //  accumulate square of mean magn. per paticle
            dEn = -m_dJ*nE+m_dH*dMag;    //  energy per particle
            m_vdEn[k++] = (dEn /= m_nN2);//  energy per particle, step k
            dMeanEn += dEn;              //  accululate mean energy per particle
            dMeanEn2 += dEn*dEn;         //  accululate square of mean energy per particle
         }
      }
      cout << setprecision(5) << setw(15) << dKbT << setw(10) << k << setw(10) << nMeasure << endl;
      //  just to inform you the program is working
      dMeanMag /= nMeasure;        //  avg magnetization per particel
      dMeanMag2 /= nMeasure;       //  avg square of the magnetization per particle
      dMeanEn /= nMeasure;         //  avg energy per particle
      dMeanEn2 /= nMeasure;        //  avg square energy per particle
      m_vdM[l] = fabs(dMeanMag);   //  store magnetization for temperature dKbT
      m_vdE[l] = dMeanEn;          //  store energy
      m_vdVarM[l] = dMeanMag2-dMeanMag*dMeanMag;  //  variance of <m_1> ~ magnetic susceptibility
      m_vdVarE[l] = dMeanEn2-dMeanEn*dMeanEn;     //  variance of <\epsilon> ~ heat capacity
      ErrVar();                    //  error analysis
      m_vdVarMn[l] = m_dVarMnEnd;  //  save error information for measurements
      m_vdVarEn[l] = m_dVarEnEnd;
      m_vdErrVarMn[l] = m_dErrVarMn;
      m_vdErrVarEn[l++] = m_dErrVarEn;
      dKbT -= m_dKbTstep;          //  reduce temperature for next step
      if (fabs(m_dKbTstep) < 1.0e-5)
         dKbT -= 1.0;              //  reduce temperature for next step
   }
}
