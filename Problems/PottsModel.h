//  Definition the CPottsModel class, Chapter 18: Potts Model

#pragma once
class CPottsModel
{
public:
   CPottsModel(void);
   CPottsModel(int,double,double,double,double=0.0,int=0);
   ~CPottsModel(void);
   void seed(int s) {  //  set seed of the system's random number generator
      srand(s);
   }
   void assign(int,double,double,double,double=0.0,int=0);
   void RunPottsModel(bool,int,int,int,int);

private:
   double RandInt(void) {  //  use sytem random number generator to generate
                           //  uniformliy distributed random numbers \in [0,1)
      return static_cast<double>(rand())/(static_cast<double>(RAND_MAX)+1.0);
   }
   void GetNeighbor(void);
   int Delta(int a, int b) {  //  calculate for spin configuration Q = 1
      return a==b ? 1 : 0;    //  see footnote on page 268
   }
   void HotStart(void);
   void SweepR(double);
   void SweepS(double);
   void ErrVar(void);

public:   //  alloow access to results
   std::vector<int> m_vnSigma;  //  spin configuration
   std::vector<double> m_vdM,   //  magnetization
      m_vdAvgM,     //  magnetization per particle, Fig. 18.2
      m_vdVarM,     //  error on <m>
      m_vdMn,
      m_vdVarMn,    //  variance <m> ~ susceptibility
      m_vdErrVarMn, //  error on variance <m> ~ error on susceptibility
      m_vdE,        //  <\epsilon>
      m_vdEn,
      m_vdAvgE,     //  energy per particle, Fig. 18.3
      m_vdVarE,     //  error on <\epsilon>
      m_vdVarEn,    //  variance <\epsilon> ~ heat capacity
      m_vdErrVarEn, //  error of var <\epsilon> ~ error on heat capacity
      m_vdTemp;     //  temperatures
   int m_nKindex;

private:
   int m_nN,             //  system size
      m_nN2,             //  system size m_nN*m_nN
      m_nTsteps,         //  number  of temperature steps
      m_nQ;              //  number of q-states
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
   matrix2d<int> m_mnNeighbors;      //  array for neighbor indices
};

