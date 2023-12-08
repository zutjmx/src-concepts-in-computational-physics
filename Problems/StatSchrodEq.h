//  Definition of the CStatSchrodEq class, Chapter 10

#pragma once

class CStatSchrodEq
{
public:
   CStatSchrodEq(void);
   CStatSchrodEq(int,int);
   ~CStatSchrodEq(void);
   void assign(int,int);
   bool Numerov(double&,double,double);
   double ExpValue(std::vector<double>& v);
   void Potential1(void);
   void Potential2(void);
   void Potential3(void);

private:
   double SolveEquation(double);
   void Normalize(void);

private:
   int m_nMax,              //  number of grid-points N
      m_nMaxIter;           //  max # of iterations
   double m_dH,             //  distance between two grid-points
      m_dH2,                //  m_dH*m_dH
      m_d5ov6,              //  5/6
      m_d1ov6,              //  1/6
      m_dPi;                //  number Pi in machine accuracy

public:  //  open access to results
   std::vector<double> m_vPhil, //  \phi_n(s_l) 
      m_vSl,                    //  grid-points s_l
      m_vSl2,                   //  grid-points s_l, squared
      m_vVl;                    //  potential v
};

