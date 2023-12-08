//  Realization of the CIsing class, Chapter 15: Ising Model

#include "stdafx.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include "matrix2d.h"
#include "Ising.h"

using namespace std;

/*****************************************************************************/
CIsing::CIsing(void)
/*****************************************************************************/
//  empty constructor, please use the property assign to define the
//  system's parameters
{
}
/*****************************************************************************/
CIsing::CIsing(int n, double J, double h, double ts, double te, int nTsteps)
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
   m_dKbTstep = (m_dKbTstart-m_dKbTend)/nTsteps; //  temperature step
   m_dJ = J;
   m_dH = h;
   m_mnNeighbors.assign(m_nN2,4,0);  //  define array for neighbor information
}
/*****************************************************************************/
CIsing::~CIsing(void)
/*****************************************************************************/
{
}
/*****************************************************************************/
void CIsing::assign(int n, double J, double h, double ts, double te, int nTsteps)
/*****************************************************************************/
//  use this property to define the parameters of the Ising model if the empty
//  constructor has been used to instantiate the CIsing class.

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
void CIsing::GetNeighbor(void)
/*****************************************************************************/
//  define the array neighbors following the recipe on page 210 
{
   for (register int i=1, j=0; i<=m_nN2; ++i, ++j) {
      if (i+m_nN <= m_nN2)            //  up
         m_mnNeighbors[j][0] = i+m_nN;
      else
         m_mnNeighbors[j][0] = i-m_nN*(m_nN-1);
      if (i%m_nN)                     //  right
         m_mnNeighbors[j][1] = i+1;
      else
         m_mnNeighbors[j][1] = i-m_nN+1;
      if (i-m_nN >= 1)                //  down
         m_mnNeighbors[j][2] = i-m_nN;
      else
         m_mnNeighbors[j][2] = i+m_nN*(m_nN-1);
      if ((i-1)%m_nN)                 //  left
         m_mnNeighbors[j][3] = i-1;
      else
         m_mnNeighbors[j][3] = i+m_nN-1;
   }
}
/*****************************************************************************/
void CIsing::SweepS(void)
/*****************************************************************************/
//  perform a sequential sweep of the lattice points according to the
//  recipe (3) execution of the code, page 212 required to generate Fig. 15.3
{
   double dCountAcc=0.0,
      dDeltaE,
      dE,
      dMag,
      dPr,
      dZ;
   int nMag,
      nE,
      nS;

   for (register int iter=0; iter<m_nN2; ++iter) {  //  sweep over lattice points
      nE = 0;
      nS = m_vnSigma[iter];   //  actual lattice point
      for (register int i=0; i<4; ++i)  //  sweep over neighbors
         nE += m_vnSigma[m_mnNeighbors[iter][i]-1];
      dDeltaE = 2.0*m_dJ*nS*nE+2.0*m_dH*nS;  //  \Delta E_{ij}, Eq. (15.51)
      dZ = exp(-dDeltaE/m_dKbT);             //  \exp{-\frac{\Delta e_{ij}}{k_BT}
      dPr = min(dZ,1.0);                     //  Eq. (15.54)
      if (dPr == 1.0 || RandInt() <= dPr)    //  Metropolis
         m_vnSigma[iter] = -nS;              //  accept new configuration: flip spin
      nMag = 0;
      dE = 0.0;
      for (register int j=0; j<m_nN2; ++j) {
         nS = m_vnSigma[j];
         nMag += nS;                         //  magnetization
         nE = 0;
         for (register int k=0; k<4; ++k)
            nE += m_vnSigma[m_mnNeighbors[j][k]-1];
         dE += -m_dJ*nS*nE-m_dH*nS;          //  energy
      }
      m_vdM[m_nKindex] = static_cast<double>(nMag)/m_nN2;  // magnetization
      m_vdE[m_nKindex] = dE/m_nN2;           //  energy
      dE = dMag = 0.0;
      for (register int j=0; j<=m_nKindex; ++j) {
         dE += m_vdE[j];      //  accumulate energies
         dMag += m_vdM[j];    //  accumulate magnetization
      }
      m_vdAvgE[m_nKindex] = dE/(m_nKindex+1);  //  avg energy after m_nKindex steps
      m_vdAvgM[m_nKindex] = dMag/(m_nKindex+1);  //  avg magnetization after
                                                 //  m_nKindex steps
      ++m_nKindex;
   }
}
/*****************************************************************************/
void CIsing::SweepR(double dKbT)
/*****************************************************************************/
//  perform a random sweep of the lattice points according to the
//  recipe (3) execution of the code, page 212
//  input parameter:
//  dKbT ... actual temperature
{
   double dDeltaE,
      dPr,
      dZ;
   int nE,
      nS,
      nIndex;

   for (register int iter=0; iter<m_nN2; ++iter) {
      nIndex = static_cast<int>(floor(m_nN2*RandInt())); //  Eq. (15.55)
      nE = 0;
      nS = m_vnSigma[nIndex];
      for (register int i=0; i<4; ++i)
         nE += m_vnSigma[m_mnNeighbors[nIndex][i]-1];
      dDeltaE = 2.0*m_dJ*nS*nE+2.0*m_dH*nS; //  energy change, Eq. (15.51)
      dZ = exp(-dDeltaE/dKbT);              //  \exp{-\frac{\Delta e_{ij}}{k_BT}
      dPr = min(dZ,1.0);                    //  Eq. (15.54)
      if (dPr == 1.0 || RandInt() <= dPr)   //  Metropolis
         m_vnSigma[nIndex] = -nS;           //  accept new configuration,
                                            //  flip spin
   }
}
/*****************************************************************************/
void CIsing::Sweep(double dKbT)
/*****************************************************************************/
//  perform a sequential sweep of the lattice points according to the
//  recipe (3) execution of the code, page 212
//  input parameter:
//  dKbT ... actual temperature
{
   double dDeltaE,
      dPr,
      dZ;
   int nS,
      nE;

   for (register int iter=0; iter<m_nN2; ++iter) {  //  sequential sweep
      nS = m_vnSigma[iter];
      nE = 0;
      for (register int i=0; i<4; ++i)
         nE += m_vnSigma[m_mnNeighbors[iter][i]-1];
      dDeltaE = 2.0*m_dJ*nS*nE+2.0*m_dH*nS; //  energy change, Eq. (15.51)
      dZ = exp(-dDeltaE/dKbT);              //  \exp{-\frac{\Delta e_{ij}}{k_BT}
      dPr = min(dZ,1.0);                    //  Eq. (15.54)
      if (dPr == 1.0 || RandInt() <= dPr)   //  Metropolis
         m_vnSigma[iter] = -nS;             //  accept new configuration,
                                            //  flip spin
   }
}
/*****************************************************************************/
void CIsing::ErrVar()
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
         nIndex = static_cast<int>(floor(nMeasure*RandInt()));
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
void CIsing::RunIsingSweeps(bool cw, int nNsweep)
/*****************************************************************************/
//  run the Ising model calculation to generate data for Fig. 15.3
//  input parameters:
//  cw ..........  = true warm start/ = false cold start
//  nNsweep .....  total number of sweeps
{
   double dMag=0.0,
      dEn=0.0;
   int nE,
      nS;

   GetNeighbor();    //  get neighbor configuration
   m_vnSigma.assign(m_nN2,1);     //  initial configuration (cold start)
   m_vdM.assign(nNsweep*m_nN2+1,0.0);     //  prepare vectors
   m_vdE.assign(nNsweep*m_nN2+1,0.0);
   m_vdAvgM.assign(nNsweep*m_nN2+1,0.0);
   m_vdAvgE.assign(nNsweep*m_nN2+1,0.0);
   if (cw) {                    //  hot start
      for (register int i=0; i<m_nN2; ++i)
         if (RandInt() < 0.5)
            m_vnSigma[i] = -1;
   }
   for (register int i=0; i<m_nN2; ++i) {  // initial configuration
      dMag += m_vnSigma[i];            //  magnetization
      nE = 0;
      nS = m_vnSigma[i];
      for (register int j=0; j<4; ++j)
         nE += m_vnSigma[m_mnNeighbors[i][j]-1];
      dEn += -m_dJ*nS*nE-m_dH*nS;      //  energy, Eq. (15.42)
   }
   m_vdAvgM[0] = m_vdM[0] = dMag/m_nN2;  //  magnetization per particle
   m_vdAvgE[0] = m_vdE[0] = dEn/m_nN2;   //  energy per particle
   m_nKindex = 1;                        //  set index for sweeps
   m_dKbT = m_dKbTstart;                 //  starting temperature
   for (register int i=0, k=0; i<nNsweep; ++i, ++k) { //  perform sweeps
      SweepS();     //  sweep lattice sequentially
//      cout << setw(5) << i << endl;      //  control output to see the
                                         //  program is working
   }
}
/*****************************************************************************/
void CIsing::RunIsingModel(bool cw, bool rs, int nMeasure, int nTherm, int nStep)
/*****************************************************************************/
//  run the Ising model calculation to generate data for all other figures
//  input parameters:
//  cw ..........  = true: warm start/ = false cold start
//  rs ..........  = true: random sweep/ = false: squetial sweep
//  nMeasure ....  number of measurements
//  nTherm ......  number of thermalization steps
//  nStep .......  number of steps between measurements
{
   int nNsweep=nTherm+(nMeasure-1)*nStep;  //  total number of sweeps
   double dKbT=m_dKbTstart,   // set temperature to start temperature
      dMag,
      dMeanMag,
      dMeanMag2,
      dEn,
      dMeanEn,
      dMeanEn2;
   int nMag,
      nE,
      nS;
   register int i=0,
      k=0;
   void (CIsing::*Sw)(double) = rs ? &CIsing::SweepR : &CIsing::Sweep;
        //  decide: random or sequential sweep

   GetNeighbor();    //  get neighbor configuration
   m_vnSigma.assign(m_nN2,1);  //  initial configuration (cold start)
   if (cw) {                   //  hot start?
      for (register int i=0; i<m_nN2; ++i)
         if (RandInt() < 0.5)
            m_vnSigma[i] = -1;
   }
   m_vdTemp.assign(m_nTsteps,0.0);  //  prepare vectors
   m_vdM.assign(m_nTsteps,0.0);     //  magetization per particle <m>
   m_vdE.assign(m_nTsteps,0.0);     //  energy per paricle <\epsilon>
   m_vdVarM.assign(m_nTsteps,0.0);  //  variance of <m>
   m_vdVarMn.assign(m_nTsteps,0.0);
   m_vdErrVarMn.assign(m_nTsteps,0.0);  //  error on the variance <m>
   m_vdVarE.assign(m_nTsteps,0.0);  //  variance of <\epsilon>
   m_vdVarEn.assign(m_nTsteps,0.0);
   m_vdErrVarEn.assign(m_nTsteps,0.0);  //  error on the variance <\epsilon> 
   m_vdMn.assign(nMeasure,0.0);
   m_vdEn.assign(nMeasure,0.0);         //  energy per particle
   while (dKbT > m_dKbTend) {           //  loop over temperatures
      m_vdTemp[i] = dKbT;               //  save actual temperature
      dMag = dMeanMag = dMeanMag2 = dMeanEn = dMeanEn2 = 0.0;
      k = 0;
      for (register int sweep=1; sweep<=nNsweep; ++sweep) {
         (this->*Sw)(dKbT);    //  perform sweep
         if (sweep>=nTherm && !(sweep%nStep)) {  // measurement?
            nMag = 0;                            // yes
            dEn = 0.0;
            for (register int j=0; j<m_nN2; ++j) {
               nS = m_vnSigma[j];
               nMag += nS;                    // magnetization
               nE = 0;
               for (register int k=0; k<4; ++k)
                  nE += m_vnSigma[m_mnNeighbors[j][k]-1];
               dEn += -m_dJ*nS*nE-m_dH*nS;    //  energy, Eq. (15.42)
            }
            m_vdMn[k] = dMag = static_cast<double>(nMag)/m_nN2;  // magnetization per particle,
                                                                 // step k
            dMeanMag += dMag;                 //  accumulate mean magn. per partcicle
            dMeanMag2 += dMag*dMag;           //  accumulate square of mean magn. per paticle
            m_vdEn[k++] = dEn /= m_nN2;       //  energy per particle, step k
            dMeanEn += dEn;                   //  accululate mean energy per particle
            dMeanEn2 += dEn*dEn;              //  accululate square of mean energy per particle
         }
      }
      cout << setprecision(5) << setw(15) << dKbT << setw(10) << k << setw(10) << nMeasure << endl;
      //  just to inform you the program is working
      dMeanMag /= nMeasure;         //  avg magnetization per particel
      dMeanMag2 /= nMeasure;        //  avg square of the magnetization per particle
      dMeanEn /= nMeasure;          //  avg energy per particle
      dMeanEn2 /= nMeasure;         //  avg square energy per particle
      m_vdM[i] = fabs(dMeanMag);    //  store magnetization for temperature dKbT
      m_vdE[i] = dMeanEn;           //  store energy
      m_vdVarM[i] = dMeanMag2-dMeanMag*dMeanMag;  //  variance of <m> ~ magnetic susceptibility
      m_vdVarE[i] = dMeanEn2-dMeanEn*dMeanEn;     //  variance of <\epsilon> ~ heat capacity
      ErrVar();                     //  error analysis
      m_vdVarMn[i] = m_dVarMnEnd;   //  save error information for measurements
      m_vdVarEn[i] = m_dVarEnEnd;
      m_vdErrVarMn[i] = m_dErrVarMn;
      m_vdErrVarEn[i] = m_dErrVarEn;
      dKbT -= m_dKbTstep;           //  reduce temperature for next step
      ++i;
   }
}
