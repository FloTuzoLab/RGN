//
//  Cell.hpp
//  BetHedging_Noise
//
//  Created by Florian on 28/05/2019.
//  Copyright © 2019 Florian. All rights reserved.
//


#ifndef Cell_hpp
#define Cell_hpp
#include "Allele.hpp"

#include <stdio.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <vector>
#include <map>

#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <random>
using namespace std;


class Cell
{
    
    
public:
    
    /* ###Methods### */
    
    
    /* ##Constructors## */
    
    Cell();
    
    // #Partially configurable constructor# //
    Cell(int nutrientTypes,
         double Kon,
         double Koff,
         double kcat,
         double KM,
         double permeability,
         double tauTranscriptionOncePolymeraseBound,
         double basalBindingProbabilityPolymerase,
         double enhancedBindingProbabilityPolymerase,
         double tauTraduction,
         double tauDegradationTranscript,
         double tauDegradationProtein,
         double cellVolume,
         double volumeRatioDescendant,
         double volumeCellMitosis,
         double volumeMutationRate,
         std::map <std::string,PhysiologicalAllele> physiologicalAlleles,
         std::vector <RegulatoryAllele> positiveRegulatoryAlleles,
         std::vector <RegulatoryAllele> negativeRegulatoryAlleles,
         gsl_rng *random
         );
    
    // #Constructor using a full pre-defined cell# //
    /*Cell(const Cell &m_cell, bool cellSizerCycle,  gsl_rng *randomSimulation); useless for BH model */
    Cell(const Cell &_cell);
    
    void operator=(const Cell &m_cell);
    
    // #Outside/pseudo constructor using some of the features of a pre-defined cell, i.e. those which are heritable by a daughter, and modifying those of the mother cell which concern products apportioned with respect to the post-mitosis volumes.# //
    void BuildDaughterCell(Cell& _motherCell);
    
    /* ##Destructor## */
    
    ~Cell();// At least needed to delete the pointer of random number.
    
    /* ##Other methods## */
    
    // #Printing the state of the cell# //
    void PrintCell();
    std::string printCellMatrix();
    void Outputfilling(std::ofstream* &_fileoutcell); //need to see what outputs will be needed
    
    // #Operations on the cell# //
    
    /* #0.Reevaluate connectivity vector# */
    //i)Evaluate connectivity matrix
    void EvaluateBSTFConnectivity();
    void ReevaluateBSTFConnectivityAfterFactorSequenceMutation(int transcriptionFactorNumber);
    void ReevaluateBSTFConnectivityAfterRegulatoryFactorSequenceMutation(int transcriptionFactorNumber);
    void ReevaluateBSTFConnectivityAfterRegulatoryPhysiologicalSequenceMutation(std::string function);
    //void ReconstructBSTFConnectivityMatrixAfterEfficiencyPhysiologicalMutation(std::string function);
    
    //ii)Reevaluate connectivity after mutation on a specific gene
    //void EvaluateGeneConnectivity();
    void BasalRateDynamics(int timeStepRatio, double delta_t);
    /* #1.Transcripts production# */
    void TranscriptsDynamics(double delta_t);// In principle, need for a general function with balance of transcript
    /* #2.Proteins production# */
    void ProteinsDynamics(double delta_t);// In principle, need for a general function with balance of proteins
    /* #3.Nutrient uptake# */
    double NutrientUptake(double concNutEnv, double concNutInt, double dt); // change here in that nutrients are counted in energy units.
    void SetNutrientsAfterDiffusion(double V_comp,double* &concNutEnv, double *concNutEnvInit, double dt);
    //bool EnergyBudgetUsefulnessTest(double delta_t);All the energy can be converted to growth
    /* #4.Energy production# */
    void EnergyProduction(double delta_t);
    /* #5.Energy use# */
    void BasalMetabolismCalculus(double delta_t); // contingent to the choice of modelization procedures
    // #i.Withdrawal by the basal metabolism
    void BasalMetabolismWithdrawal(double delta_t); //also decreases efficiency of basal mechanisms
    // #ii.a)Cell growth
    void Growth(); //via phospholipids production
    /*// #ii.b)Cell Nut1 molecule production
    void Nut1EnzymesProduction(int totalNumberProteins,bool energyWholeBudgetUsefulness,double delta_t);
    // #ii.c)Cell Public Goods production
    void Nut2EnzymesProduction(int totalNumberProteins,bool energyWholeBudgetUsefulness,double delta_t);*/ //useless
    // #ii. Whole energy allocation process
    //void EnergyBalance(bool energyWholeBudgetUsefulness,double delta_t); //needed only if model relies on emptying the cell at the end
    
    /* #6.Phenotype effects */
    // #i.Maintenance metabolites effect# //
    //void GarbageNeutralization(double K_garbageNeutralization,double delta_t); //No garbage in bet hedging models
    /*// #iii.Nutrient release if remainder# // //needed only if model relies on emptying the cell at the end
    double StandardNutrientRelease();
    double RevealedNutrientRelease();*/
    /* #7.Mitosis checking */
    // #i.Checking if replication checkpoint reached#
    bool ReplicationCheckpoint(); //should be called mitosis
    /* #8.Death checking */
    bool DeathTest(double delta_t); // need to remove internal death causation
    /* #9.Mutation checking and outcome */
    void MutationProcessing();
    
    // #Getting the settings# //
    
    double GetConcentrationNut1();
    double GetPermeabilityCell();
    double GetAllometryBasalAlpha();
    double GetTauTranscriptionOncePolymeraseBound();
    double GetBasalBindingProbabilityPolymerase();
    double GetEnhancedBindingProbabilityPolymerase();
    double GetTauDegradationTranscript();
    double GetTauDegradationProtein();
    double GetVolumeCell();
    double GetSACell();
    int GetTranscripts(std::string t);
    int GetTFTranscripts(int g);//t-ieth type of transcript
    int GetProteins(std::string t);
    int GetTranscriptionFactors(int g);//p-ieth type of protein
    int GetGenerationNumber();
    double GetMitosisSize();
    
    // #Setting the settings# //
    void MatrixConversion(std::vector <int> &matrixCode);
    //Preparing results //
    
    void GeneTabularConversion(std::string &geneCode);
    
private:
    
    /* ##Attributes## */
    /* 0.Generaion number */
    int m_generationNumber;
    
    /* 1.Chemico-physical properties of the cell */
    //o. Afinities between Transcription factors and Binding-sites
    double m_K__on;                       //binding reaction constant
    double m_K__off;                      //dissociation reaction constant
    //i.Membrane permeability
    double m_permeabilityCell;            //Parameter: as such, it is not subject to evolution
    //ii.Allometric coefficients for the efficiency of production (allometry production) and the greediness of basal needs
    /*double m_allometryProductionAlpha;    //Parameter: as such, it is not subject to evolution*/
    //double m_allometryBasalAlpha;         //Parameter: as such, it is not subject to evolution
    //iii.Conversion rate between Nutrient and Product or By-product
    double m_k__cat;
    double m_K__M;
    /*double m_etaNutrientToGarbage; Garbage needed if intrinsic death rate*/
    //iv.Cost of protein production
    /*std::map <std::string,double> m_metaboliteCost;*/
    //v.Duration of DNA replication
    /*double m_synthesisDuration;*/
    
    /* 2.Transcripts dynamics */
    /*int m_connectivityMatrix[14][10];*/
    int **m_connectivityMatrix;                     //Connectivity between TF and BS
    int *m_totalBindingSites;                       //Number of each type binding-sites in the whole genome
    //i.Transcription rates
    double m_tauTranscriptionOncePolymeraseBound;   //Need to be defined specifically by the intertwining between TF and BS
    double m_basalBindingProbabilityPolymerase;
    double m_enhancedBindingProbabilityPolymerase;  //+ function to pick randomly degradation
    //ii.Degradation rate
    double m_tauDegradationTranscript;
    
    /* 3.Proteins and metabolites dynamics */
    //i.Traduction rate
    double m_tauTraduction;
    /*double m_etaTranslationActivityPostCheckpoint;*/
    //ii.Degradation rate
    double m_tauDegradationProtein;
    //iii.Metabolite and Enzymes
    //double m_enzymeProductionRate;
    /*double m_metaboliteMaintenanceProductionRate;
    double m_metaboliteCommonGoodsProductionRate;
    double m_metaboliteDegradationRate;*/
    
    /* 4.Dynamical physical features */
    //i.Size of the cell
    double m_volumeCell;
    //ii.Current time since DNA synthesis beginning
    /*double m_timeSinceSynthesisBeginning;*/
    
    /* 5.Content of the cell */
    //o. Calibration of basal metabolism constraints
    double m_currentBasalMetabolism;
    double m_synthesisRateProp;
    int m_ribosomesAmount;
    int m_mRNA_ribosomes;
    //i.Energetic content
    int *m_numberNutUnits;
    //int m_basalEnergyUnits;
    double m_numberEnergyUnits;//needed to distinguish converted nutrients from that not available to growth
    //ii.By-product content: garbage
    /*double m_numberGarbageUnits;*/
    //iii.Transcripts content
    //int m_numberGrowthTranscripts;
    /*int m_numberMaintenanceTranscripts;*/
    int m_numberNut1Transcripts;
    int m_numberNut2Transcripts;
    //int m_numberSizerTranscripts;
    int *m_numberTFTranscripts; //Transcripts of transcription factors coding genes
    //iv.Proteins content
    //int m_numberGrowthProteins;
    /*int m_numberMaintenanceProteins;*/
    int m_numberNut1Proteins; //Enzymes in the model
    int m_numberNut2Proteins;
    /*int m_numberSizerProteins; Sizer useful when dealing with cell size evolution that we think is involved in transition towards multicellularity */
    int *m_numberTranscriptionFactors;
    //iv.Maintenance metabolite content
    //double m_numberNut1Enzymes;
    //double m_numberNut2Enzymes;
    int m_numberNutrientTypes;
    
    /* 6.Genetics */
    double m_volumeMutationRate;
    double m_volumeRatioDescendant; //Conserved to be able to influence the ability to differentiate
    /*double m_garbageProportionDescendant;*/
    /*bool m_cellCycleSizer;*/
    double m_volumeCellMitosis; //Provisory: need to be defined by genes
    /*bool m_sizerConcentration;*/
    /*double m_sizerThresholdMitosis;*/ //Featuring the concentration threshold of sizer at which DNA synthesis starts in order to reproduce.
    /*std::string m_cellCycle; //type of cell cycle, including G2 or not, with simplification made: the cell is either in G1 either in S/G2/M*/
    
    /*Genes*/
    
    //Physiological genes//
    std::vector <std::string> m_physiologicalFunctions;
    std::map <std::string,PhysiologicalAllele> m_physiologicalAlleles;
    
    //Regulatory genes//
    std::vector <RegulatoryAllele> m_regulatoryAlleles;
    
    /* 7.Genomics */
    
    /* double m_sizeMinimalGenome; //A cell cannot be smaller than the size of the minimal genome allowing life; coherent with the use of a sizer because if mitosis is determined by a threshold, cells having no sizer will divide continously and reach a size of zero which means nothing at al inside the cell. */
    
    /* 8.Random generator */
    gsl_rng *_Random;//(pseudo) random number used in the processing for the diferent draws.
    
};

/* ###Generic functions### */

/* #Conversion of a volume into a surface area# */
double VtoSAconversion(double V_cell);

/* #Conversion of a surface area into a volume# */
double SAtoVconversion(double SA_cell);

/* #Calculus of a concentration for a molecule in a given cell# */
long double ConcentrationCalculus(double MoleculeNumber, double CellVolume);

double Molecule_micro_toMole_L_Concentration(double MoleculeConc);

/* #Proportion calculus of a given type entity among all the similar entities regardless of their type# */
double Proportion(double sub, double total);

/* Product of elements contained in a vector, raised to the power contained in another vector */
double WeightProduct(const int nb_elements, double* elements,int* repeat_numbers);

void GetKey_valuesMapofPhysiologicalAllele(std::vector<std::string> &keys, std::map<std::string,PhysiologicalAllele> Map);

#endif /* Cell_hpp */
