//
//  Medium.hpp
//  Cell_Evolution
//
//  Created by Florian on 21/02/2018.
//  Copyright © 2018 Florian. All rights reserved.
//

#ifndef Medium_hpp
#define Medium_hpp

#include <cstdio>
#include <iostream>
#include <string>
#include <math.h>
#include <map>
#include <vector>

class Medium{
    
public:
    
    /* ###Methods### */
    
    
        /* ##Constructors## */
    
    
            // #Default constructor# //
            Medium();
    
            // #Total constructor# //
            Medium(double diffusionNutrients,double diffusionCommonGoods, double diffusionCells);
    
            // #Constructor using a pre-defined medium# //
            Medium(const Medium& m_medium);
    
        //No need for a destructor: Nevertheless, need to ask!//
    
    
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
            void Afficher() const;
    
    
private:
    
    /* ##Attributes## */
        double m_diffusionNutrients;    //The diffusion coefficient of nutrients between the compartment and its outside environment; can possibly be modified to account for different coefficients regarding the type of the molecules concerned.
        double m_diffusionCommonGoods;  //The diffusion coefficient of Common Goods between the compartment and its outside environment.
        double m_diffusionCells;        //The diffusion oefficient of cells between the compartment and its outside environment. Can possibly be set to 0.
};

#endif /* Medium_hpp */
