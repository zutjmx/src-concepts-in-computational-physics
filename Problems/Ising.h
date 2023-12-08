//  Definition of the CIsing class, Chapter 15: Ising Model

#pragma once
class CIsing
{
public:
   CIsing(void);
   CIsing(int,double,double,double,double=0.0,int=0);
   ~CIsing(void);
   void seed(int s) {  //  set seed of the system's random number generator
      srand(s);
   }
   void assign(int,double,double,double,double=0.0,int=0);
   void RunIsingSweeps(bool,int);
   void RunIsingModel(bool,bool,int,int,int);

private:
   double RandInt(void) {  //  use sytem random number generator to generate
                           //  uniformliy distributed random numbers \in [0,1)
      return static_cast<double>(rand())/(static_cast<double>(RAND_MAX)+1);
   }
   void GetNeighbor(void);
   void Sweep(double);
   void SweepS(void);
   void SweepR(double);
   void ErrVar(void);

public:   //  alloow access to results
   std::vector<int> m_vnSigma;  //  spin configuration
   std::vector<double> m_vdM,   //  magnetization
      m_vdAvgM,     //  magnetization per particle, Fig. 15.3
      m_vdVarM,     //  error on <m>
      m_vdVarMn,    //  variance <m> ~ susceptibility
      m_vdErrVarMn, //  error on variance <m> ~ error on susceptibility
      m_vdE,        //  <\epsilon>
      m_vdAvgE,     //  energy per particle, Fig. 15.3
      m_vdVarE,     //  error on <\epsilon>
      m_vdVarEn,    //  variance <\epsilon> ~ heat capacity
      m_vdErrVarEn, //  error of var <\epsilon> ~ error on heat capacity
      m_vdTemp;     //  temperatures
   int m_nKindex;

private:
   int m_nN,             //  system size
      m_nN2,             //  system size m_nN*m_nN
      m_nTsteps;         //  number  of temperature steps
   double m_dKbTstart,   //  starting temperature
      m_dKbTend,         //  final temperature
      m_dKbTstep,        //  temperature step
      m_dJ,              //  exchange parameter
      m_dH,              //  external field
      m_dKbT,            //  actual temperature
      m_dVarEnEnd,       //  vectors for error analysis
      m_dErrVarEn,
      m_dVarMnEnd,
      m_dErrVarMn;
   std::vector<double> m_vdMn,       //  <m> per step
      m_vdEn;                        //  <\epsilon> per step
   matrix2d<int> m_mnNeighbors;      //  array for neighbor indices
};
