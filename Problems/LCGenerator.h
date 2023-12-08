// Declaration of the CLCGenerator class, linear congruential generator
// Section 12.2

#pragma once

class CLCGenerator
{
public:
   CLCGenerator(void);
   CLCGenerator(unsigned long long,unsigned long long,
      unsigned long long,unsigned long long);
   ~CLCGenerator(void);
   void assign(unsigned long long,unsigned long long,
      unsigned long long,unsigned long long);
   double rand(void);
   void seed(unsigned long long s) {
      m_nX0 = s;
   }
   void SpectralTest(int);
   void Histogram(int,int);

private:
   unsigned long long m_nA,  //  Park-Miller parameter a
      m_nC,                  //  Park-Miller parameter c
      m_nM,                  //  Park-Miller parameter m
      m_nX0;                 //  seed
   double m_dM;              //  Park-Miller parameter m

public:    //  allow access to results
   double m_dHistScale;      //  scale of histogram
   std::vector<double> m_vdX1,  // spectral test
      m_vdX2,
      m_vdRange;             //  histrogram range
   std::vector<unsigned int> m_vnCount; // vector with counts
};
