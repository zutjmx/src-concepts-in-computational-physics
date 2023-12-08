//  Definition of the CRejectSamp class, Section 13.3 Rejection metod

#pragma once
class CRejectSamp
{
public:
   CRejectSamp(void);
   CRejectSamp(double,double,double,double(*)(double),double(*)(double),
      double(*)(double));
   ~CRejectSamp(void);
   void assign(double,double,double,double(*)(double),double(*)(double),
      double(*)(double));
   void seed(int s) {
      srand(s);
   }
   void Histogram(int,int);

private:
   double RandInt(void) {  //  use sytem random number generator to generate
                           //  uniformliy distributed random numbers \in [0,1)
      return static_cast<double>(rand())/(static_cast<double>(RAND_MAX)+1);
   }

private:
   double (*m_Qx)(double),  //  pointer to the function which defines the required pdf
      (*m_Hx)(double),      //  pointer to the envelope function
      (*m_InvHx)(double);   //  pointer to the inverse envelope function
   double m_dLambda,        //  envelop function: parameter \lambda, Eq. (13.24)
      m_dSigma,             //  required pdf: parameter \sigma, Eq. (13.33)
      m_dCmin;              //  c_min, Eq. (13.41)
   std::vector<double> m_vdX;  //  collect random numbers

public:   //  enable access to results
   double m_dHistScale;           //  histogram scale
   std::vector<int> m_vnCount;    //  histogram: counts per bin
   std::vector<double> m_vdRange; //  histogram: bin range
};

