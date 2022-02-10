//
//  Compartiment.hpp
//  Cell_Evolution
//
//  Created by Florian on 16/02/2018.
//  Copyright © 2018 Florian. All rights reserved.
//

//A unit of a 2D space
#ifndef Compartiment_hpp
#define Compartiment_hpp

#include "Medium.hpp"

class Compartiment{
    
public:
    
    /* ###Methods### */
    
    
        /* ##Constructors## */
    
    
            // #Default constructor# //
            Compartiment();
    
            // #Half constructor# //
            Compartiment(double length, double height,double concentrationNutrientStandard,double concentrationNutrientHidden, double concentrationNutrientRevealed, double concentrationCommonGoods, Medium medium,double KcommonGoodsKinetics);
    
            // #Total free constructor# //
            Compartiment(int reference[2], double const length, double height, Medium medium, double KcommonGoodsKinetics,double concentrationNutrientStandard, double concentrationNutrientHidden, double concentrationNutrientRevealed, double concentrationCommonGoods);
    
            // #Constructor using a pre-defined compartment# //
            Compartiment(const Compartiment& m_compartiment);
    
            //#Operator= # //
            void operator=(const Compartiment& _compartment);
    
        /* ##Destructor## */
    
        ~Compartiment(); //Deletion of compartment, useful because of the two pointers used
    
    
        /* ##Other methods## */
    
    
            // #Setting the settings# //
            void SetReference(int r1, int r2);
            void SetNeighboursRelationships(int numNeighbour,int r1, int r2);
            void SetMedium(double diffusionNutrients, double diffusionCommonGoods, double diffusionCells);
            void SetConcentrationNutrientStandard(double a);
            void SetConcentrationNutrientHidden(double a);
            void SetConcentrationNutrientRevealed(double a);
            void SetConcentrationCommonGoods(double a);
    
            // #Getting the settings# //
            int GetNeighbour(int numNeighbour,int numCoordinate) const;
            double GetDiffusionNutrients() const;
            double GetDiffusionCommonGoods() const;
            double GetDiffusionCells() const;
            double GetConcentrationNutrientStandard() const;
            double GetConcentrationNutrientHidden() const;
            double GetConcentrationNutrientRevealed() const;
            double GetConcentrationCommonGoods() const;
            double GetHeight() const;
    
            // #Printing the state of the compartment# //
            void PrintState();//Print the state of the compartment in order to test the different methods.
            void Outputfilling();

            // #Operations on the compartment# //
            void DepleteCompartment(std::map <std::string,double> &totalNutrientBalance, double membranePermeability,double totalMembraneSurface, double delta_t); //Depletion of a compartment by a given cell.
            void EnrichCompartment(double numberPublicGoodsReleased,double numberStandardNutrientsReleased, double numberRevealedNutrientsReleased);
            void HiddenNutrientBalance(double delta_t); //Calculus of the balance of nutrient following the revelation of hidden nutrients.
            double Volume(); //Volume of the compartment.
    
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
        double m_KcommonGoodsKinetics; //Reaction rate between common Goods and hidden nutrients (lacks information on the unit)
    
        //ii.Compartment concentrations
        double m_concentrationNutrientStandard; //The concentration of the standard nutrient, i.e. the nutrient always available.
        double m_concentrationNutrientHidden; //The concentration of the hidden nutrient, available only after its revelation by a "common good" molecule.
        double m_concentrationNutrientRevealed; //The total concentration of nutrients in the compartment including the nutrients revealed by the common goods.
        double m_concentrationCommonGoods; //The concentration in Common Goods produced by cells, which can interact wth hidden nutrients thus revealing them.
    
};

    /* ###Generic functions### */

        /* ## Calculus of the volume of a cube## */
        double VolumeCalculus(double length,double height); //Calculus of the volume generated (cubic shape as mentionned above)

#endif /* Compartiment_hpp */
