//  Definition of the CFractRandWalk class, Section 17.4, fractal random walk

#pragma once
class CFractRandWalk
{
public:
   CFractRandWalk(void);
   CFractRandWalk(int,double,double);
   ~CFractRandWalk(void);
   void seed(int s) {  //  set seed of the system's random number generator
      srand(s);
   }
   void assign(int,double,double);
   void RunFractRandWalk(void);

private:
   double RandInt(void) {  //  use sytem random number generator to generate
                           //  uniformliy distributed random numbers \in [0,1)
      return static_cast<double>(rand())/(static_cast<double>(RAND_MAX)+1.0);
   }

private:
   int m_nN;        //  number of time steps
   double m_dBeta,  //  parameter \beta, Eq. (17.81)
      m_dInvBeta,   //  1/\beta
      m_dTau;       //  minimal waiting time

public:   //  allow access to results
   std::vector<double> m_vdX,  //  x(t)
      m_vdT;                   //  t
};

