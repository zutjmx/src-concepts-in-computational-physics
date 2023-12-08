//  Realization of the CRandWalk class, Chapter 17, Random Walk

#include "stdafx.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include "RandWalk.h"

using namespace std;

/*****************************************************************************/
CRandWalk::CRandWalk(void)
/*****************************************************************************/
//  empty constructor, please use the property assign to define the
//  system's parameters
{
}
/*****************************************************************************/
CRandWalk::CRandWalk(int n, double p)
/*****************************************************************************/
//  input parameters:
//  n ... number of time steps
//  p ... probability p, step 1, page 245
{
   m_nSteps = n;
   m_dP = p;
   m_dQ = 1.0-m_dP;  //  probability q
   m_vnX.assign(m_nSteps,0);  //  start at x(0) = 0
}
/*****************************************************************************/
CRandWalk::~CRandWalk(void)
/*****************************************************************************/
{
}
/*****************************************************************************/
void CRandWalk::assign(int n, double p)
/*****************************************************************************/
//  use this property to define the parameters of the random walk if the empty
//  constructor has been used to instantiate the CRandWalk class.

//  input parameters:
//  n ... number of time steps
//  p ... probability p, step 1, page 245
{
   m_nSteps = n;
   m_dP = p;
   m_dQ = 1.0-m_dP;  //  probability q
   m_vnX.assign(m_nSteps,0);  //  start at x(0) = 0
}
/*****************************************************************************/
void CRandWalk::RunRandWalk()
/*****************************************************************************/
{
   int r;

   for (register int i=1; i<m_nSteps; ++i) {
      r = RandInt() < m_dP ? +1 : -1;  //  step 2 and 3, page 245
      m_vnX[i] = m_vnX[i-1]+r;         //  step 4, page 245
   }
}