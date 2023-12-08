//  Definition of the COrnsteinUhlenbeck class, Section 17.3

#pragma once
class COrnsteinUhlenbeck
{
public:
   COrnsteinUhlenbeck(void);
   ~COrnsteinUhlenbeck(void);
   COrnsteinUhlenbeck(int,double,double,double);
   void assign(int,double,double,double);
   void RunOrnUhl(double,double);

private:
   int m_nSteps;      //  number of time steps
   double m_dBeta,    //  parameter \beta of Eq. (17.57)
      m_dA,           //  parameter A of Eq. (17.57)
      m_dTstep,       //  time step
      m_dVar;

public:  //  alow access to results
   std::vector<double> m_vdV,  //  velocities
      m_vdX,                   //  trajectory
      m_vdTime;                //  time
};

