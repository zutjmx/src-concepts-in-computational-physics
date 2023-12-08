//  Definition of the CFractLevy class, Section 17.4, fractal time Lévy flight

#pragma once
class CFractLevy
{
public:
   CFractLevy(void);
   CFractLevy(int,double,double,double,double);
   ~CFractLevy(void);
   void seed(int s) {  //  set seed of the system's random number generator
      srand(s);
   }
   void RunFractLevy(void);
   void assign(int,double,double,double,double);

private:
   double RandInt(void) {  //  use sytem random number generator to generate
                           //  uniformliy distributed random numbers \in [0,1)
      return static_cast<double>(rand())/(static_cast<double>(RAND_MAX)+1.0);
   }

private:
   int m_nN;         //  number of time steps
   double m_dAlpha,  //  Lévy index \alpha
      m_dInvAlpha,   //  1/\alpha
      m_dL,          //  minimal flight length
      m_dBeta,       //  parameter \beta, Eq. (17.81)
      m_dInvBeta,    //  1/\beta
      m_dTau;        //  minimal waiting time

public:   //  allow access to results
   std::vector<double> m_vdT,  //  t
      m_vdX;                   //  x(t)
};

