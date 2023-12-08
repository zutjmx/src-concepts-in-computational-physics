// Definition of the CNormDist class
//
// CNormDist generates randum numbers according to a normal distribution
// using direct sampling described in Chapter 13, Section 13.1 Introduction
#pragma once
class CNormDist
{
public:
   CNormDist(void);
   ~CNormDist(void);
   void seed(int s) {      //  set the seed of the system random number generator
      srand(s);
   }

private:
   double RandInt(void) {  //  generate uniformly distributed random numbers in [0,1)
                           //  using the system random number generator
      return static_cast<double>(rand())/(static_cast<double>(RAND_MAX)+1);
   }

public:
   void sample(double& z1, double& z2);  //  sample the normal distribution
                                         //  according to Eq. (13.8)

private:
   double m_dPi;           //  number Pi in system accuracy
};

