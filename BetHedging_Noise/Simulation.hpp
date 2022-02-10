//
//  Simulation.hpp
//  BetHedging_Noise
//
//  Created by Florian on 02/07/2019.
//  Copyright © 2019 Florian. All rights reserved.
//

#ifndef Simulation_hpp
#define Simulation_hpp

#include <stdio.h>
#include <fstream>
#include <sstream>
#include "Puddle.hpp"

class Simulation
{
    
    
public:
    
    Simulation();
    
    //Constructor of Etienne Rajon's bewildering (or wild!) donut: but is it hotter in the outside inside the hole or in the fryer?//
    Simulation(int _simulationNumber,
               int _initialCellAmount,
               int _nutrientTypeAmount,
               int _TFamount,
               int _compartmentAmountLinear,
               double _switchProbability,
               double _nutrientConcentrationRatio,
               double _mutationRate,
               double _volumeMutationRate,
               double _aPropDesc,
               std::ofstream* fileoutcompartments,
               std::ofstream* fileoutcells);
    
    Simulation(const Simulation &_simulation);
    
    void operator=(const Simulation &_simulation);
    
    //Destructor//
    ~Simulation();
    
    //Diffusion between compartments
    void SimulationProcessing(int _simulationDuration);
    
private:
    
    int m_delta_t;
    int m_simulationTimeStep;
    int m_simulationNumber;
    double m_nutrientConcentrationRatio;
    DonutPuddle m_donutPuddle;
    
    gsl_rng *RandomSimulationNumber;
    
    
};

#endif /* Simulation_hpp */


