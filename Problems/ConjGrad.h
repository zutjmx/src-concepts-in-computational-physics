//  Definition of the CConjGrad class, Appendix H: Deterministic Optimization

#pragma once
class CConjGrad
{
public:
   CConjGrad(void);
   CConjGrad(double(*)(double,double),double(*)(double,double),
      double(*)(double,double));
   ~CConjGrad(void);
   int Iterate(double x, double y);
   void assign(double(*)(double,double),double(*)(double,double),
      double(*)(double,double));

private:
   double (*m_Fxy)(double,double),  //  pointer to function f(x,y)
      (*m_GradFx)(double,double),   //  pointer to a function which returns
                                    //  the first partial derivative of f(x,y) with
                                    //  respect to x
      (*m_GradFy)(double,double);   //  pointer to a function which returns
                                    //  the first partial derivative of f(x,y) with
                                    //  respect to y

public:  //  allow access to results
   std::vector<double>m_vdXiter, //  x-coordinate of the n-th iteration step
      m_vdYiter;   //  y-coordinate of the n-th iteration step
   int m_nImax;    //  number of iteration steps
};

