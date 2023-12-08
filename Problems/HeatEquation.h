//  Definition of the CHeatEquation class, Chapter 9

#pragma once

#include <vector>

class CHeatEquation
{
public:
   CHeatEquation(void);
   CHeatEquation(int,double,double,double,double);
   ~CHeatEquation(void);
   void assign(int,double,double,double,double);
   void SetHeatSource(double,double,double);
   int SolveEquation(void);

private:
   int TriDiag(void);

//  Variables

private:
   int m_nN;               //  # of temperature steps
   double m_dT0,           //  temperature T_0
      m_dTN,               //  temperature T_N
      m_dKappa,            //  kappa, thermal diffusivity
      m_dL,                //  Length
      m_dH,                //  sub-spacing
      m_dH2,               //  h^2
      m_dTheta,            //  Theta, maximum height
      m_dEll;              //  half width
   std::vector<double>m_vDiag,  //  diagonal of the A-Matrix
      m_vSuperd,                //  super-diagonal of the A-Matrix
      m_vSubd,                  //  sub-diagonal of the A-Matrix
      m_vGamma,                 //  Heat Sink/Source
      m_vRhs;                   //  right hand side

public:   //  give access to the results
   std::vector<double> m_vTemp, //  temperatures at the grid-points
      m_vGridPoint;             //  Grid-Points
};

