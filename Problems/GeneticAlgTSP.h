//  Definition of the CGeneticAlgTSP class, Section 20.4: genetic algorithms
//  traveling salesperson problem

#pragma once
class CGeneticAlgTSP
{
public:
   CGeneticAlgTSP(void);
   CGeneticAlgTSP(int,int,int,int,double);
   ~CGeneticAlgTSP(void);
   void assign(int,int,int,int,double);
   void RegularGrid(void);
   void Initialize(void);
   void RunGeneticAlgTSP();
   double GetLength(int*);
   void seed(int s) {  //  set seed of the system's random number generator
      srand(s);
   }

private:
   double RandInt(void) {  //  use sytem random number generator to generate
                           //  uniformliy distributed random numbers \in [0,1)
      return static_cast<double>(rand())/(static_cast<double>(RAND_MAX)+1.0);
   }
   std::vector<int>& FlipLR(int *,int);
   void GetTour(std::vector<int>&);
   double LengthT(std::vector<int>& o);
   void SortLength(std::vector<double>&,std::vector<int>&);

private:
   int m_nCitiesX,    //  number of lattice points in x-direction
      m_nCitiesY,     //  number of lattice points in y-direction
      m_nPopulations, //  number of populations
      m_nGenerations, //  number of generations
      m_nNumberCities;//  total # of cities
   double m_dMutp;    //  mutation probability
   std::vector<int>  m_vnX1,  //  position of city in x-direction
      m_vnY1;                 //  position of city in y-direction
   matrix2d<int> m_mnIndividuals;  //  array to store the individuals
                                   //  together with their routes

public:  //  allow access to results
   double m_dLength;               //  optimum length
   int m_nFinGeneration;           //  final number of generations
   std::vector<int> m_vnX,         //  final x-positions
      m_vnY,                       //  final y-postions
      m_vnRoute;                   //  final route
};

