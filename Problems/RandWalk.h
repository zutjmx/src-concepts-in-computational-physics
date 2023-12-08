//  Definition of the CRandWalk class, Chapter 17, Random Walk

#pragma once
class CRandWalk
{
public:
   CRandWalk(void);
   CRandWalk(int,double);
   ~CRandWalk(void);
   void seed(int s) {  //  set seed of the system's random number generator
      srand(s);
   }
   void assign(int,double);
   void RunRandWalk(void);

private:
   double RandInt(void) {  //  use sytem random number generator to generate
                           //  uniformliy distributed random numbers \in [0,1)
      return static_cast<double>(rand())/(static_cast<double>(RAND_MAX)+1.0);
   }

private:
   int m_nSteps;   //  number of time steps
   double m_dP,    //  probability p
      m_dQ;        //  probability q

public:  //  allow access to results
   std::vector<int> m_vnX;   //  x(t)
};

