//  Definition of the CLevy class, Section 17.4, Lévy flight

#pragma once
class CLevy
{
public:
   CLevy(void);
   CLevy(int,double,double);
   ~CLevy(void);
   void seed(int s) {  //  set seed of the system's random number generator
      srand(s);
   }
   void assign(int,double,double);
   void RunLevy(void);
   void RunLevy2d(void);
   void RunWiener2d(void);

private:
   double RandInt(void) {  //  use sytem random number generator to generate
                           //  uniformliy distributed random numbers \in [0,1)
      return static_cast<double>(rand())/(static_cast<double>(RAND_MAX)+1.0);
   }

private:
   int m_nN;         //  number of time steps
   double m_dAlpha,  //  Lévy index \alpha
      m_dL,          //  minimal flight length
      m_dInvAlpha,   //  1/\alpha
      m_dPi;         //  number Pi in machine accuracy

public:   //  allow access to results
   std::vector<double> m_vdX,  //  contains x(t)
      m_vdT;                   //  contains t
};
