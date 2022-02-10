//
//  Puddle.hpp
//  BetHedging_Noise
//
//  Created by Florian on 02/07/2019.
//  Copyright © 2019 Florian. All rights reserved.
//

#ifndef Puddle_hpp
#define Puddle_hpp

#include <stdio.h>
#include <cstdio>
#include <string>
#include <math.h>
#include "Ecosystem.hpp"

class DonutPuddle : public Ecosystem
{
public:
    
    DonutPuddle();
    
    //Constructor of Etienne Rajon's bewildering (or wild!) donut: but is it hotter in the outside inside the hole or in the fryer?//
    DonutPuddle(gsl_rng *random_Simulation,
                double switchProbability,
                double Nut1Input,
                double Nut2Input,
                double delta_t,
                double delta_l,
                int timeStepRatio,
                int lengthPuddle,
                int widthPuddle,
                Compartment insideCompartment,
                Compartment outsideCompartment,
                int initialNumberCell,
                Cell cell,
                std::ofstream* fileoutcompartments,
                std::ofstream* fileoutcells);
    
    DonutPuddle(const DonutPuddle &_donutPuddle);
    
    void operator=(const DonutPuddle &_donutPuddle);
    
    //Destructor//
    ~DonutPuddle();
    
    //Diffusion between compartments
    void DiffusionBetweenCompartments();
    
private:
    
};

#endif /* Puddle_hpp */
