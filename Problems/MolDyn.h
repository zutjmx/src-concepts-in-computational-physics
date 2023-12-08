//  Definition of the CMolDyn class

#pragma once

#include "matrix2d.h"

class CMolDyn
{
public:
   CMolDyn(void);
   CMolDyn(int,int);
   ~CMolDyn(void);
   void assign(int,int);
   void SetParameter(double,double,double,double,double);
   void SetTimeParam(double t0, double tau) {
   //  input parameters:
   //  t0 .... starting time
   //  tau ... time step
      m_dT0 = t0;
      m_dTau = tau;
   }
   void InitSqareLattice(double,double,double,double);
   double LeapFrog(double);

private:
   int m_nMax,    //  # of particles
      m_nNum;     //  # of particles per row in a square lattice
   double m_dMass,//  particle mass
      m_dG,       //  acceleration due to gravity
      m_dSigma,   //  parameter \sigma of the Lennard-Jones potential, Eq. (7.4)
      m_dEps,     //  parameter \epsilon of the Lennard-Jones potential, Eq. (7.4)
      m_dLength,  //  box length (x-axis extension)
      m_dH0,      //  y-ccordinate of the bottom left particle
      m_dXstart,  //  x-ccordinate of the bottom left particle
      m_dT0,      //  starting time
      m_dTau;     //  time step
   matrix2d <double> m_mFxij,  //  two-particle force, x-component
      m_mFyij;                 //  two-particle force, y-component
   std::vector<double> m_vFxi, //  external force, x-component
      m_vFyi;                  //  external force, y-component

public:    //  open acces to results
   std::vector<double> m_vXi,  //  x coordinate of particle i
      m_vYi,                   //  y coordinate of particle i
      m_vVxi,                  //  velocity in x direction of particle i
      m_vVyi,                  //  velocity in y direction of particle i
      m_vV;                    //  velocity of particle i
};

