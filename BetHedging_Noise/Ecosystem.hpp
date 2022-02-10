//
//  Ecosystem.hpp
//  BetHedging_Noise
//
//  Created by Florian on 01/07/2019.
//  Copyright © 2019 Florian. All rights reserved.
//

#ifndef Ecosystem_hpp
#define Ecosystem_hpp

#include <stdio.h>
#include "Compartment.hpp"
#include "Cell.hpp"
#include <vector>
#include <deque>

class Ecosystem{
    
public:
    
    /* ###Methods### */
    
    
    /* ##Constructors## */
    
    Ecosystem();
    
    // #Total-free constructor# //
    Ecosystem(gsl_rng *random_Simulation, double switchProbability, double delta_t, double delta_l, int timeStepRatio, int lengthPuddle, int widthPuddle, bool adjacencyNeed, double Nut1Input,
              double Nut2Input,Compartment insideCompartment, Compartment outsideCompartment, int initialNumberCell, Cell cell,std::ofstream* _fileoutcompartments,std::ofstream* _fileoutcells);
    
    Ecosystem(const Ecosystem& _ecosystem);
    
    void operator=(const Ecosystem& _ecosystem);
    
    /* ##Destructor## */
    
    ~Ecosystem(); //Deletion of a puddle, useful because two pointers are used here
    
    
    /* ##Other methods## */
    
    /* #Print the state of the puddle# */
    void PrintEcosystem(); //adjacent mesh contains lateral outside compartments influencing the inside of the studied puddle.
    
    /* #Processing in the puddle# */
    
    /* #1.Nutrient diffuson towards cells# */
    void NutrientDiffusion(int cell, double* concNutEnvPre,double freeCompartmentVolume);
    
    
    /* #2.Birth-death definition# */
    //i.Birth of a cell
    void CellReproduction(int cell); //Not yet defined
    //ii.Death of a cell
    bool CellDeath(int cell,double delta_t); //Not yet defined
    //iii.Cell dynamics
    void CellBiomoleculesDynamics(int cell,double delta_t_cellDyn);
    void CellEnergeticDynamics(int cell,double delta_t);
    /* #3.Processing of a time-step in a compartement# */
    void TimeStepInsideEcosystem(int t_current, double delta_t); //Compartment time step proceedings
    void OutputGetting(int simulationNumber,int t_output);
    
    /* ## Final.Processing of a time-step in the puddle */
    //void TimeStepPuddle(); //Puddle time step proceedings need to be defined
    
protected:
    
    /* ##Attributes## */
    
    /* 1.Dimensions of the puddle */
    int  m_lengthInside;
    int  m_widthInside;
    double  m_compartmentVolume;
    // Need for an adjacent mesh
    bool m_adjacencyNeed;
    
    /* 2.Creation of a puddle as a double pointer, each pointer of pointer containing a compartment as if it was an array */
    Compartment **m_puddle;//Pointer on a compartment of the environment//
    
    /* 3.Initialization of the neighbours for each of the compartment in the puddle */
    int ****m_neighbours;
    
    /* 4.Definition of the parameters for diffusion, function of time step, length of a compartment and diffusion coefficients */
    double m_alphaDF_empty[4][2]; //m_alpha[x][0] represents diffusion rate in water when there is no excluded volume
    double m_alphaDF_actual[4][2]; //represents actual diffusion rate accounting for excluded volume
    double m_cellDiffusionRate;
    double m_nutrientSwitchProbability;
    double m_Nutrient1Input;
    double m_Nutrient2Input;
    
    /* 5.Creation of a dynamic vector containing the cells living in the ecosystem at a given timestep */
    std::vector <Cell> m_cell;
    double m_permeabilityCell; //Defined at the level of the puddle where diffusion takes place, and once and for all because this does not evolve in the model
    
    /* 6.Creation of a dynmic vector of dynamic vector containing any of the living cell compartment position */
    std::vector < std::vector < int > > m_positionCell;
    
    /* 7.Creation of an integer storing the ratio between gene products dynamics and other processes*/
    int m_timeStepRatio;//between gene network dynamics and other functions
    double m_dt_NutUptake;
    
    /* 8.Creation of a pointer containing a random number*/
    gsl_rng * _SimuRandom_;
    
    /* 9.Pointer on files */
    std::ofstream* fileoutcompartments;
    std::ofstream* fileoutcells;
};

/* ###Generic functions### */

/* #Definition of the Explicit Finite Differences operator# */
//double OperatorFiniteDifferencesExplicit(double alpha, double c, int nbNeighbours, double n1,double n2,double n3,double n4,double n5,double n6); Why is this here?

#endif /* Ecosystem_hpp */
