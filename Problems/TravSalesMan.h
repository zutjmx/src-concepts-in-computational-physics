//  Definition the CTravSalesMan class, Section 20.3: Simulated Annealing

#pragma once
class CTravSalesMan
{
public:
   CTravSalesMan(void);
   CTravSalesMan(int,int);
   ~CTravSalesMan(void);
   void assign(int,int);
   void RunTravSalesMan(int,double,double);
   void RegularLattice(void);
   void ClusterCities(int);
   void Initialize(void);
   void seed(int s) {  //  set seed of the system's random number generator
      srand(s);
   }

private:
   double RandInt(void) {  //  use sytem random number generator to generate
                           //  uniformliy distributed random numbers \in [0,1)
      return static_cast<double>(rand())/(static_cast<unsigned int>(RAND_MAX)+1);
   }
   double LengthT(std::vector<int>&);
   double CalcT0(void);

private:
   int m_nCitiesX,    //  number of lattice points in x-direction
      m_nCitiesY,     //  number of lattice points in y-direction
      m_nNumberCities;//  total # of cities
   std::vector<int> m_vnX1,  //  position of city in x-direction
      m_vnY1,                //  position of city in y-direction
      m_vnInArr;             //  array with city number

public:  //  allow access to results
   double m_dT0;             //  starting temperature
   std::vector<int> m_vnX,   //  x-positions
      m_vnY,                 //  y-positions
      m_vnRoute;             //  route by city numbers
   std::vector<double> m_vdTemp,  //  temperature T for Figs. 20.3 and 20.4
      m_vdMitl,                   //  <H>_T(T)
      m_vdVar;                    //  C_h(T)
};

