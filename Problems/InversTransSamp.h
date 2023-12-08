//  Definition of the CInversTransSamp class, Section 13.2

#pragma once
class CInversTransSamp
{
public:
   CInversTransSamp(void);
   CInversTransSamp(double,double(*)(double));
   ~CInversTransSamp(void);
   void seed(int s) {
      srand(s);
   }
   void Histogram(int,int);
   void assign(double,double(*)(double));

private:
   double RandInt(void) {  //  use sytem random number generator to generate
                           //  uniformliy distributed random numbers \in [0,1)
      return static_cast<double>(rand())/(static_cast<double>(RAND_MAX)+1);
   }

private:
   double (*m_Fx)(double);  //  pointer to inverse function
   double m_dLambda;        //  parameter \lambda, Eq. 13.24
   std::vector<double> m_vdX;  //  generated random numbers

public:  //  allow acces to results
   double m_dHistScale;  //  histogram scale
   std::vector<int> m_vnCount;  //  counts per bin
   std::vector<double> m_vdRange; //  histogram bins
};

