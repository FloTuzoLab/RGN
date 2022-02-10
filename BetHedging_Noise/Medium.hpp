//
//  Medium.hpp
//  BetHedging
//
//  Created by florian labourel on 30/06/2019.
//  Copyright © 2019 florian labourel. All rights reserved.
//

#ifndef Medium_hpp
#define Medium_hpp

#include <stdio.h>
#include <cstdio>
#include <iostream>
#include <string>
#include <math.h>
#include <map>
#include <vector>

using namespace std;

class Medium{
    
public:
    
    /* ###Methods### */
    
    
    /* ##Constructors## */
    Medium();
    // #Total constructor# //
    Medium(double diffusionNutrients, double diffusionCells);
    
    // #Constructor using a pre-defined medium# //
    Medium(const Medium& m_medium);
    
    //No need for a destructor: Nevertheless, need to ask!//
    // #Operator = # //
    void operator=(const Medium& _medium);
    
    /* ##Other methods## */
    
    
    // #Setting the settings# //
    void SetDiffusionNutrients(double a);
    void SetDiffusionCommonGoods(double a);
    void SetDiffusionCells(double a);
    
    // #Getting the settings# //
    double GetDiffusionNutrients() const;
    double GetDiffusionCommonGoods() const;
    double GetDiffusionCells() const;
    
    // #Printing of the medium features# //
    void Print() const;
    
    
private:
    
    /* ##Attributes## */
    double m_diffusionNutrients;    //The diffusion coefficient of nutrients between the compartment and its outside environment; can possibly be modified to account for different coefficients regarding the type of the molecules concerned.
    double m_diffusionCommonGoods;  //The diffusion coefficient of Common Goods between the compartment and its outside environment.
    double m_diffusionCells;        //The diffusion oefficient of cells between the compartment and its outside environment. Can possibly be set to 0.
};

#endif /* Medium_hpp */
