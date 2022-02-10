//
//  Compartment.hpp
//  BetHedging
//
//  Created by florian labourel on 01/07/2019.
//  Copyright © 2019 florian labourel. All rights reserved.
//

#ifndef Compartment_hpp
#define Compartment_hpp

#include <stdio.h>
#include "Medium.hpp"

class Compartment{
    
public:
    
    /* ###Methods### */

    /* ##Constructors## */
    
    // #Half constructor# //
    //Compartment(double length, double height,double concentrationNutrientStandard,double concentrationNutrientHidden, double concentrationNutrientRevealed, double concentrationCommonGoods, Medium medium,double KcommonGoodsKinetics);
    Compartment();
    
    // #Unlocated compartment# //
    Compartment(double const length, double height, Medium medium, double concentrationNut1, double concentrationNut2);
    // #Total free constructor# //
    Compartment(int reference[2], double const length, double height, Medium medium, double concentrationNut1, double concentrationNut2);
    
    // #Constructor using a pre-defined compartment# //
    Compartment(const Compartment& m_compartment);
    
    //#Operator= # //
    void operator=(const Compartment& _compartment);
    
    /* ##Destructor## */
    
    ~Compartment(); //Deletion of compartment, useful because of the two pointers used
    
    
    /* ##Other methods## */
    
    
    // #Setting the settings# //
    void SetReference(int r1, int r2);
    void SetNeighboursRelationships(int numNeighbour,int r1, int r2);
    void SetMedium(double diffusionNutrients, double diffusionCells);
    void SetConcentrationNut1(double a);
    void SetConcentrationNut2(double a);
    void SetActualVolume(double actualVolume);
    
    // #Getting the settings# //
    int GetNeighbour(int numNeighbour,int numCoordinate) const;
    double GetDiffusionNutrients() const;
    //double GetDiffusionCommonGoods() const;
    double GetDiffusionCells() const;
    double GetConcentrationNut1() const;
    double GetConcentrationNut2() const;
    double GetHeight() const;
    double GetActualVolume() const;
    double* ModifyConcentrationNut1();
    double* ModifyConcentrationNut2();
    
    // #Printing the state of the compartment# //
    void PrintState();//Print the state of the compartment in order to test the different methods.
    //void Outputfilling(); //unnecessary
    
    // #Operations on the compartment# //
    //void DepleteCompartment(std::map <std::string,double> &totalNutrientBalance, double membranePermeability,double totalMembraneSurface, double delta_t); //Depletion of a compartment by a given cell.
    /*void HiddenNutrientBalance(double delta_t); //Calculus of the balance of nutrient following the revelation of hidden nutrients.
     */ //useless for this instance of the model
    double Volume(); //Volume of the compartment.
    double m_ActualVolume;
    
    /*##public variables*/
    
    //ii.Compartment concentrations
    double m_concentrationNut1; //The concentration of the standard nutrient, i.e. the nutrient always available.
    double m_concentrationNut2; //The concentration of the hidden nutrient, available only after its revelation by a "common good" molecule.
    //Removal of interactions within the compartment
    
private:
    
    /* ##Attributes## */
    
    /* 1.References of the compartment and its neighbours */
    int *m_labelCompartment; //Pointer vector containing the compartment label.
    int **m_labelNeighbours; //Pointer vector containing the compartment neighbours labels (as referenced above, i.e. in a vector): 4 neighbours for each compartment, including possible referencing of the outside environment.
    
    /* 2.Compartment shape characterization */
    //The compartment has a cubic shape.
    double m_length; //Compartement length of a compartment
    double m_height;  //Compartment width, i.e. width of the puddle
    
    /* 3.Compartment chemical properties */
    
    //i.Diffusion properties
    Medium m_medium;
    
   
};

/* ###Generic functions### */

/* ## Calculus of the volume of a cube## */
double VolumeCalculus(double length,double height); //Calculus of the volume generated (cubic shape as mentionned above)

#endif /* Compartment_hpp */
