//
//  Cell.cpp
//  Cell_Evolution
//
//  Created by Florian on 02/03/2018.
//  Copyright © 2018 Florian. All rights reserved.
//

#include "Cell.hpp"
#include <sstream>

using namespace std;

/*  ##Constructors##    */



    /*  #Default Constructor: only some basic physical properties defined#  */

Cell::Cell() /*These two features are not subject to evolution and make no sense when fixed to 0; hence, it is better to define them even in the default construtor.*/ :
/* 0.Generaion number */
m_generationNumber(0),

/*1.Properties of the cell*/

//Fixed:
//i.Chemico-physical properties
m_K__on(pow(10,-1)),                                        //Binding reaction constant
m_K__off(pow(10,-2)),                                       //Dissciation reaction constant
m_permeabilityCell(pow(10,-1)),                             //Permeability of the cell membrane
m_allometryProductionAlpha(1),         //Allometry of energy production: how scale the production of energy per nutrient unit with cell volume
m_allometryBasalAlpha(1),                   //Allometry of basal metabolism: how scale the cost of basal metabolism with cell volume
m_etaNutrientToEnergy(30),                                          //Conversion rate between a nutrient and energy
m_etaNutrientToGarbage(1),                       //Conversion rate between a nutrient and a garbage
m_synthesisDuration(100),

//Variable:
//ii.Transcription rate of genes                                    //Further, these rates will rely on BSTF dynamics
m_tauTranscriptionOncePolymeraseBound(0.1),
m_basalBindingProbabilityPolymerase(0.1),
m_enhancedBindingProbabilityPolymerase(1),
//iii.Traduction rate of transcripts
m_tauTraduction(5),
m_etaTranslationActivityPostCheckpoint(0.5),
//iv.Degradation rate of gene products
m_tauDegradationTranscript(0.01),               //Same degradation rate for each transcript
m_tauDegradationProtein(0.01),                  //Same degradation rate for each protein
//v.Dynamics rate of metabolites aside from chemical interactions of interest
m_metaboliteGrowthProductionRate(10),
m_metaboliteMaintenanceProductionRate(1000),
m_metaboliteCommonGoodsProductionRate(100),
m_metaboliteDegradationRate(0.01),              //Same degradation rate as proteins, without further knowledge

//v.Size of the cell
m_volumeCell(100),
//vi.Current time since DNA synthesis beginning
m_timeSinceSynthesisBeginning(0),

/*2.Content of the cell*/
//o.Calibration of basal metabolism constraints
m_SA_calibration(1),
m_basalMetabolismSAcalibration(pow(10,6)),
m_currentBasalMetabolism(0),
//i.Energetic content
m_numberStandardNutrientUnits(0),                                   //In a first time, nutrient units are immediately consumed
m_numberRevealedNutrientUnits(0),
//m_basalEnergyUnits(0),
m_numberEnergyUnits(0),                                             //In a first time, energy units are immediately consumed
//ii.By-product content: garbage
m_numberGarbageUnits(0),
//iii.Transcripts content
m_numberGrowthTranscripts(50),
m_numberMaintenanceTranscripts(50),
m_numberPublicGoodsTranscripts(0),
m_numberSizerTranscripts(50),
//iv.Proteins content
m_numberGrowthProteins(0),
m_numberMaintenanceProteins(0),
m_numberPublicGoodsProteins(0),
m_numberSizerProteins(0),
//v.Metabolites content
m_numberMaintenanceMetabolites(0),
m_numberPublicGoodsMetabolites(0),

/*3.Genetics*/

//i.Cell volume at mitosis
m_volumeProportionDescendant(0.5),
m_garbageProportionDescendant(0.5),
m_cellCycleSizer(true),
m_volumeCellMitosis(500),
m_volumeMutationRate(500),
m_sizerConcentration(false),
m_sizerThresholdMitosis(1000),
m_sizeMinimalGenome(pow(10,-3)),
//2.Genes definition
m_cellCycle("GSM")

{
    m_metaboliteCost.insert(pair<string,double>("Growth",1));
    m_metaboliteCost.insert(pair<string,double>("Maintenance",1));
    m_metaboliteCost.insert(pair<string,double>("Public Goods",1));              //Cost for :{"Growth","Maintenance","Public goods"}
    //rnd=randomSimulation;
    
    int numberTFgenes=10;
    int numberPhysiologicalgenes=4;
    m_numberTFTranscripts=new int[numberTFgenes];
    m_numberTranscriptionFactors=new int[numberTFgenes];
    m_totalBindingSites=new int[numberTFgenes];
    
    for (int g=0;g<numberTFgenes;g++){
        m_numberTFTranscripts[g]=0;
        m_numberTranscriptionFactors[g]=0;
        m_totalBindingSites[g]=0;
    }
    
    m_connectivityMatrix=new int*[numberTFgenes+numberPhysiologicalgenes];
    for (int Gene=0;Gene<(numberTFgenes+numberPhysiologicalgenes);Gene++){
        m_connectivityMatrix[Gene]=new int[numberTFgenes];
    }
    
    for (int Gene=0;Gene<numberTFgenes+numberPhysiologicalgenes;Gene++){
        for (int typeOfTFNumber=0;typeOfTFNumber<numberTFgenes;typeOfTFNumber++){
            m_connectivityMatrix[Gene][typeOfTFNumber]=0;
        }
    }
    for (int typeOfTFNumber=0;typeOfTFNumber<numberTFgenes;typeOfTFNumber++){
        m_totalBindingSites[typeOfTFNumber]=0;
    }
    
    int lengthRegulatorySequences(1000);
    int lengthBindingSites(8);
    string genericSequence;
    string genericBindingSite;
    int simulationNumber(1);
    gsl_rng *rnd=gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rnd, simulationNumber);
    Random=rnd;
    
    for(int pos=0;pos<lengthRegulatorySequences;pos++){
        int drawResult=gsl_ran_flat(Random, 0, 4);
        if (drawResult<1){
            genericSequence.push_back('A');
        }
        else if (drawResult<2){
            genericSequence.push_back('C');
        }
        else if (drawResult<3){
            genericSequence.push_back('G');
        }
        else{
            genericSequence.push_back('T');
        }
    }
    
    for(int pos=0;pos<lengthBindingSites;pos++){
        int drawResult=gsl_ran_flat(Random, 0, 4);
        if (drawResult<1){
            genericBindingSite.push_back('A');
        }
        else if (drawResult<2){
            genericBindingSite.push_back('C');
        }
        else if (drawResult<3){
            genericBindingSite.push_back('G');
        }
        else{
            genericBindingSite.push_back('T');
        }
    }
    
    PhysiologicalAllele genericPhysiologicalAllele(genericSequence,pow(10,-6),true,true,rnd);
    
    
    m_physiologicalAlleles.insert(pair<string,PhysiologicalAllele>("Growth",genericPhysiologicalAllele));
    m_physiologicalAlleles.insert (pair<string,PhysiologicalAllele>("Maintenance",genericPhysiologicalAllele));
    m_physiologicalAlleles.insert(pair<string,PhysiologicalAllele>("Public Goods",genericPhysiologicalAllele));
    m_physiologicalAlleles.insert(pair<string,PhysiologicalAllele>("Cell sizer",genericPhysiologicalAllele));
    
    GetKey_valuesMapofPhysiologicalAllele(m_physiologicalFunctions, m_physiologicalAlleles);
    
    RegulatoryAllele posreg(genericSequence,pow(10,-6),rnd,genericBindingSite,1);
    RegulatoryAllele negreg(genericSequence,pow(10,-6),rnd,genericBindingSite,-1);
    for (int typeTFnumber=0;typeTFnumber<numberTFgenes/2;typeTFnumber++){
        m_regulatoryAlleles.push_back(posreg);
        m_regulatoryAlleles.push_back(negreg);
    }
    gsl_rng_free(rnd);
}

    /* #Partially configurable constructor# */

Cell::Cell(double K__on,
           double K__off,
           double permeability,
           double allometryProductionCoefficient,
           double allometryBasalCoefficient,
           double etaNutrientToGarbage,
           int synthesisDuration,
           double tauTranscriptionOncePolymeraseBound,
           double basalBindingProbabilityPolymerase,
           double enhancedBindingProbabilityPolymerase,
           double tauTraduction,
           double etaTranslationActivityPostCheckpoint,
           double tauDegradationTranscript,
           double tauDegradationProtein,
           map <string,double> metaboliteCost,
           double SA_calibration,
           double basalMetabolismSAcalibration,
           double cellVolume,
           double volumeProportionDescendant,
           double garbageProportionDescendant,
           double volumeCellMitosis,
           double volumeMutationRate,
           bool sizerConcentration,
           double sizerThresholdMitosis,
           string cellCycle,
           map <string,PhysiologicalAllele> physiologicalAlleles,
           vector <RegulatoryAllele> positiveRegulatoryAlleles,
           vector <RegulatoryAllele> negativeRegulatoryAlleles
           ):
    /* 0.Generaion number */
    m_generationNumber(0),

    /*1.Properties of the cell*/

//Fixed:
    //i.Chemico-physical properties
    m_K__on(K__on),                                                     //Binding reaction constant
    m_K__off(K__off),                                                   //Dissociation reaction constant
    m_permeabilityCell(permeability),                                   //Permeability of the cell membrane
    m_allometryProductionAlpha(allometryProductionCoefficient),         //Allometry of energy production: how scale the production of energy per nutrient unit with cell volume
    m_allometryBasalAlpha(allometryBasalCoefficient),                   //Allometry of basal metabolism: how scale the cost of basal metabolism with cell volume
    m_etaNutrientToEnergy(30),                                          //Conversion rate between a nutrient and energy
    m_etaNutrientToGarbage(etaNutrientToGarbage),                       //Conversion rate between a nutrient and a garbage
    m_synthesisDuration(synthesisDuration),

//Variable:
    //ii.Transcription rate of genes                                    //Further, these rates will rely on BSTF dynamics
    m_tauTranscriptionOncePolymeraseBound(tauTranscriptionOncePolymeraseBound),
    m_basalBindingProbabilityPolymerase(basalBindingProbabilityPolymerase),
    m_enhancedBindingProbabilityPolymerase(enhancedBindingProbabilityPolymerase),
    //iii.Traduction rate of transcripts
    m_tauTraduction(tauTraduction),
    m_etaTranslationActivityPostCheckpoint(etaTranslationActivityPostCheckpoint),
    //iv.Degradation rate of gene products
    m_tauDegradationTranscript(tauDegradationTranscript),               //Same degradation rate for each transcript
    m_tauDegradationProtein(tauDegradationProtein),                     //Same degradation rate for each protein
    //v.Dynamics rate of metabolites aside from chemical interactions of interest
    m_metaboliteGrowthProductionRate(10),
    m_metaboliteMaintenanceProductionRate(10000),
    m_metaboliteCommonGoodsProductionRate(1000),
    m_metaboliteDegradationRate(0.01),              //Same degradation rate as proteins, without further knowledge
    //vi.Size of the cell
    m_volumeCell(cellVolume),
    m_timeSinceSynthesisBeginning(0),

    /*2.Content of the cell*/
    //o.Calibration of basal metabolism constraints
    m_SA_calibration(SA_calibration),
    m_basalMetabolismSAcalibration(basalMetabolismSAcalibration),
    m_currentBasalMetabolism(0),
    //i.Energetic content
    m_numberStandardNutrientUnits(0),                                   //In a first time, nutrient units are immediately consumed
    m_numberRevealedNutrientUnits(0),                                          //In a first time, nutrient units are immediately consumed
    //m_basalEnergyUnits(0),
    m_numberEnergyUnits(0),                                             //In a first time, energy units are immediately consumed
    //ii.By-product content: garbage
    m_numberGarbageUnits(0),//numberGarbageUnits),
    //iii.Transcripts content
    m_numberGrowthTranscripts(5),
    m_numberMaintenanceTranscripts(5),
    m_numberPublicGoodsTranscripts(0),
    m_numberSizerTranscripts(0),
    //iv.Proteins content
    m_numberGrowthProteins(5000),
    m_numberMaintenanceProteins(5000),
    m_numberPublicGoodsProteins(0),
    m_numberSizerProteins(0),
    //v.Maintenance metabolites content
    m_numberMaintenanceMetabolites(0),
    m_numberPublicGoodsMetabolites(0),

    /*3.Genetics*/

    //i.Cell volume at mitosis
    m_volumeProportionDescendant(volumeProportionDescendant),
    m_garbageProportionDescendant(garbageProportionDescendant),
    m_volumeCellMitosis(volumeCellMitosis),
    m_volumeMutationRate(volumeMutationRate),
    m_sizerConcentration(sizerConcentration),
    m_sizerThresholdMitosis(sizerThresholdMitosis),
    //ii.Genes
    m_cellCycle(cellCycle)
{
    m_sizeMinimalGenome=1;
    int simulationNumber(1);
    gsl_rng *rnd=gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rnd, simulationNumber);
    Random=rnd;
    
    m_metaboliteCost=metaboliteCost; //Cost for :{"Growth","Maintenance","Public goods"}
    
    int numberTFgenes= (int) positiveRegulatoryAlleles.size()+ (int) negativeRegulatoryAlleles.size();
    m_numberTFTranscripts=new int[numberTFgenes];
    m_numberTranscriptionFactors=new int[numberTFgenes];
    m_totalBindingSites=new int[numberTFgenes];
    for (int g=0;g<numberTFgenes;g++){
        m_numberTFTranscripts[g]=0;//numberOtherTranscripts[g];
        m_numberTranscriptionFactors[g]=0;//numberOtherProteins[g];
        m_totalBindingSites[g]=0;
    }
    
    m_physiologicalAlleles=physiologicalAlleles;
    GetKey_valuesMapofPhysiologicalAllele(m_physiologicalFunctions, m_physiologicalAlleles);
    
    int totalNumberGenes=numberTFgenes+ (int) physiologicalAlleles.size();
    
    m_connectivityMatrix=new int*[totalNumberGenes];
    for (int Gene=0;Gene<totalNumberGenes;Gene++){
        m_connectivityMatrix[Gene]=new int[numberTFgenes];
    }
    
    for (int Gene=0;Gene<totalNumberGenes;Gene++){
        for (int typeOfTFNumber=0;typeOfTFNumber<numberTFgenes;typeOfTFNumber++){
            m_connectivityMatrix[Gene][typeOfTFNumber]=0;
        }
    }
    
    for (int typeOfTFNumber=0;typeOfTFNumber<numberTFgenes;typeOfTFNumber++){
        m_totalBindingSites[typeOfTFNumber]=0;
    }
    for (int typeTFnumber=0;typeTFnumber<positiveRegulatoryAlleles.size();typeTFnumber++){
        m_regulatoryAlleles.push_back(positiveRegulatoryAlleles[typeTFnumber]);
        
    }
    for (int typeTFnumber=0;typeTFnumber<negativeRegulatoryAlleles.size();typeTFnumber++){
        m_regulatoryAlleles.push_back(negativeRegulatoryAlleles[typeTFnumber]);
    }
    gsl_rng_free(rnd);
}



// #Constructor using a full pre-defined cell# //

Cell::Cell(const Cell &m_cell, bool cellSizerCycle, gsl_rng *randomSimulation):
    /* 0.Generaion number */
    m_generationNumber(m_cell.m_generationNumber),

    /*1.Properties of the cell*/

//Fixed:
    //i.Chemico-physical properties
    m_K__on(m_cell.m_K__on),                                                                       //Binding reaction constant
    m_K__off(m_cell.m_K__off),                                                                     //Dissciation reaction constant
    m_permeabilityCell(m_cell.m_permeabilityCell),                                                 //Permeability of the cell membrane
    m_allometryProductionAlpha(m_cell.m_allometryProductionAlpha),                                 //Allometry of energy production: how scale the production of energy per nutrient unit with cell volume
    m_allometryBasalAlpha(m_cell.m_allometryBasalAlpha),                                           //Allometry of basal metabolism: how scale the cost of basal metabolism with cell volume
    m_etaNutrientToEnergy(m_cell.m_etaNutrientToEnergy),                                           //Conversion rate between a nutrient and energy
    m_etaNutrientToGarbage(m_cell.m_etaNutrientToGarbage),                                         //Conversion rate between a nutrient and a garbage
    m_metaboliteCost(m_cell.m_metaboliteCost),        //Cost for :{"Growth","Maintenance","Public goods"}

    m_synthesisDuration(m_cell.m_synthesisDuration),

//Variable:
    //ii.Transcription rate of genes                                                               //Further, these rates will rely on BSTF dynamics
    m_tauTranscriptionOncePolymeraseBound(m_cell.m_tauTranscriptionOncePolymeraseBound),
    m_basalBindingProbabilityPolymerase(m_cell.m_basalBindingProbabilityPolymerase),
    m_enhancedBindingProbabilityPolymerase(m_cell.m_enhancedBindingProbabilityPolymerase),
    //iii.Traduction rate of transcripts
    m_tauTraduction(m_cell.m_tauTraduction),
    m_etaTranslationActivityPostCheckpoint(m_cell.m_etaTranslationActivityPostCheckpoint),
    //iv.Degradation rate of gene products
    m_tauDegradationTranscript(m_cell.m_tauDegradationTranscript),                                 //Same degradation rate for each transcript
    m_tauDegradationProtein(m_cell.m_tauDegradationProtein),                                       //Same degradation rate for each protein
    //v.Metabolite dynamics rates
    m_metaboliteMaintenanceProductionRate(m_cell.m_metaboliteMaintenanceProductionRate),
    m_metaboliteCommonGoodsProductionRate(m_cell.m_metaboliteCommonGoodsProductionRate),
    m_metaboliteGrowthProductionRate(m_cell.m_metaboliteGrowthProductionRate),
    m_metaboliteDegradationRate(m_cell.m_metaboliteDegradationRate),
    //vi.Size of the cell
    m_volumeCell(m_cell.m_volumeCell),
    m_timeSinceSynthesisBeginning(m_cell.m_timeSinceSynthesisBeginning),

    /*2.Content of the cell*/

    //o.Calibration of basal metabolism constraints
    m_SA_calibration(m_cell.m_SA_calibration),
    m_basalMetabolismSAcalibration(m_cell.m_basalMetabolismSAcalibration),
    m_currentBasalMetabolism(m_cell.m_currentBasalMetabolism),
    //i.Energetic content
    m_numberStandardNutrientUnits(m_cell.m_numberStandardNutrientUnits),                           //In a first time, nutrient units are immediately consumed
    m_numberRevealedNutrientUnits(m_cell.m_numberRevealedNutrientUnits),                           //In a first time, nutrient units are immediately consumed
    //m_basalEnergyUnits(m_cell.m_basalEnergyUnits),
    m_numberEnergyUnits(m_cell.m_numberEnergyUnits),                                               //In a first time, energy units are immediately consumed
    //ii.By-product content: garbage
    m_numberGarbageUnits(m_cell.m_numberGarbageUnits),
    //iii.Transcripts content
    m_numberGrowthTranscripts(m_cell.m_numberGrowthTranscripts),
    m_numberMaintenanceTranscripts(m_cell.m_numberMaintenanceTranscripts),
    m_numberPublicGoodsTranscripts(m_cell.m_numberPublicGoodsTranscripts),
    m_numberSizerTranscripts(m_cell.m_numberSizerTranscripts),
    //iv.Proteins content
    m_numberGrowthProteins(m_cell.m_numberGrowthProteins),
    m_numberMaintenanceProteins(m_cell.m_numberMaintenanceProteins),
    m_numberPublicGoodsProteins(m_cell.m_numberPublicGoodsProteins),
    m_numberSizerProteins(m_cell.m_numberSizerProteins),
    //v.Metabolites content
    m_numberMaintenanceMetabolites(m_cell.m_numberMaintenanceMetabolites),
    m_numberPublicGoodsMetabolites(m_cell.m_numberPublicGoodsMetabolites),

    /*3.Genetics*/

    //i.Cell volume at mitosis
    m_volumeProportionDescendant(m_cell.m_volumeProportionDescendant),
    m_garbageProportionDescendant(m_cell.m_garbageProportionDescendant),
    m_cellCycleSizer(cellSizerCycle),
    m_volumeCellMitosis(m_cell.m_volumeCellMitosis),
    m_volumeMutationRate(m_cell.m_volumeMutationRate),
    m_sizerConcentration(m_cell.m_sizerConcentration),
    m_sizerThresholdMitosis(m_cell.m_sizerThresholdMitosis),
    m_sizeMinimalGenome(m_cell.m_sizeMinimalGenome),

    m_cellCycle(m_cell.m_cellCycle),
    m_physiologicalAlleles(m_cell.m_physiologicalAlleles)

{
    Random=randomSimulation;
    
    GetKey_valuesMapofPhysiologicalAllele(m_physiologicalFunctions, m_physiologicalAlleles);
    
    int numberTFgenes=(int) m_cell.m_regulatoryAlleles.size();
    int totalGeneNumber=numberTFgenes+(int) m_cell.m_physiologicalAlleles.size();
    
    m_numberTFTranscripts=new int[numberTFgenes];
    m_numberTranscriptionFactors=new int[numberTFgenes];
    m_totalBindingSites=new int[numberTFgenes];
    m_connectivityMatrix=new int*[totalGeneNumber];
    for (int Gene=0;Gene<totalGeneNumber;Gene++){
        m_connectivityMatrix[Gene]=new int[numberTFgenes];
    }
    
    for (int g=0;g<numberTFgenes;g++){
        m_numberTFTranscripts[g]=m_cell.m_numberTFTranscripts[g];
        m_numberTranscriptionFactors[g]=m_cell.m_numberTranscriptionFactors[g];
        m_totalBindingSites[g]=m_cell.m_totalBindingSites[g];
    }
    for (int Gene=0;Gene<totalGeneNumber;Gene++){
        for (int typeOfTFNumber=0;typeOfTFNumber<numberTFgenes;typeOfTFNumber++){
            m_connectivityMatrix[Gene][typeOfTFNumber]=m_cell.m_connectivityMatrix[Gene][typeOfTFNumber];
        }
    }
    for (int typeOfTFNumber=0;typeOfTFNumber<numberTFgenes;typeOfTFNumber++){
        m_totalBindingSites[typeOfTFNumber]=m_cell.m_totalBindingSites[typeOfTFNumber];
    }
    for (int typeTFnumber=0;typeTFnumber<numberTFgenes;typeTFnumber++){
        m_regulatoryAlleles.push_back(m_cell.m_regulatoryAlleles[typeTFnumber]);
    }
}

Cell::Cell(const Cell &m_cell):
/* 0.Generaion number */
m_generationNumber(m_cell.m_generationNumber),

/*1.Properties of the cell*/

//Fixed:
//i.Chemico-physical properties
m_K__on(m_cell.m_K__on),                                                                       //Binding reaction constant
m_K__off(m_cell.m_K__off),                                                                     //Dissciation reaction constant
m_permeabilityCell(m_cell.m_permeabilityCell),                                                 //Permeability of the cell membrane
m_allometryProductionAlpha(m_cell.m_allometryProductionAlpha),                                 //Allometry of energy production: how scale the production of energy per nutrient unit with cell volume
m_allometryBasalAlpha(m_cell.m_allometryBasalAlpha),                                           //Allometry of basal metabolism: how scale the cost of basal metabolism with cell volume
m_etaNutrientToEnergy(m_cell.m_etaNutrientToEnergy),                                           //Conversion rate between a nutrient and energy
m_etaNutrientToGarbage(m_cell.m_etaNutrientToGarbage),                                         //Conversion rate between a nutrient and a garbage
m_metaboliteCost(m_cell.m_metaboliteCost),        //Cost for :{"Growth","Maintenance","Public goods"}
m_synthesisDuration(m_cell.m_synthesisDuration),

//Variable:
//ii.Transcription rate of genes                                                               //Further, these rates will rely on BSTF dynamics
m_tauTranscriptionOncePolymeraseBound(m_cell.m_tauTranscriptionOncePolymeraseBound),
m_basalBindingProbabilityPolymerase(m_cell.m_basalBindingProbabilityPolymerase),
m_enhancedBindingProbabilityPolymerase(m_cell.m_enhancedBindingProbabilityPolymerase),
//iii.Traduction rate of transcripts
m_tauTraduction(m_cell.m_tauTraduction),
m_etaTranslationActivityPostCheckpoint(m_cell.m_etaTranslationActivityPostCheckpoint),
//iv.Degradation rate of gene products
m_tauDegradationTranscript(m_cell.m_tauDegradationTranscript),                                 //Same degradation rate for each transcript
m_tauDegradationProtein(m_cell.m_tauDegradationProtein),                                       //Same degradation rate for each protein
//v.Metabolite dynamics rates
m_metaboliteMaintenanceProductionRate(m_cell.m_metaboliteMaintenanceProductionRate),
m_metaboliteCommonGoodsProductionRate(m_cell.m_metaboliteCommonGoodsProductionRate),
m_metaboliteGrowthProductionRate(m_cell.m_metaboliteGrowthProductionRate),
m_metaboliteDegradationRate(m_cell.m_metaboliteDegradationRate),
//vi.Size of the cell
m_volumeCell(m_cell.m_volumeCell),
m_timeSinceSynthesisBeginning(m_cell.m_timeSinceSynthesisBeginning),

/*2.Content of the cell*/

//o.Calibration of basal metabolism constraints
m_SA_calibration(m_cell.m_SA_calibration),
m_basalMetabolismSAcalibration(m_cell.m_basalMetabolismSAcalibration),
m_currentBasalMetabolism(m_cell.m_currentBasalMetabolism),
//i.Energetic content
m_numberStandardNutrientUnits(m_cell.m_numberStandardNutrientUnits),                           //In a first time, nutrient units are immediately consumed
m_numberRevealedNutrientUnits(m_cell.m_numberRevealedNutrientUnits),                           //In a first time, nutrient units are immediately consumed
//m_basalEnergyUnits(m_cell.m_basalEnergyUnits),
m_numberEnergyUnits(m_cell.m_numberEnergyUnits),                                               //In a first time, energy units are immediately consumed
//ii.By-product content: garbage
m_numberGarbageUnits(m_cell.m_numberGarbageUnits),
//iii.Transcripts content
m_numberGrowthTranscripts(m_cell.m_numberGrowthTranscripts),
m_numberMaintenanceTranscripts(m_cell.m_numberMaintenanceTranscripts),
m_numberPublicGoodsTranscripts(m_cell.m_numberPublicGoodsTranscripts),
m_numberSizerTranscripts(m_cell.m_numberSizerTranscripts),
//iv.Proteins content
m_numberGrowthProteins(m_cell.m_numberGrowthProteins),
m_numberMaintenanceProteins(m_cell.m_numberMaintenanceProteins),
m_numberPublicGoodsProteins(m_cell.m_numberPublicGoodsProteins),
m_numberSizerProteins(m_cell.m_numberSizerProteins),
//v.Metabolites content
m_numberMaintenanceMetabolites(m_cell.m_numberMaintenanceMetabolites),
m_numberPublicGoodsMetabolites(m_cell.m_numberPublicGoodsMetabolites),

/*3.Genetics*/

//i.Cell volume at mitosis
m_volumeProportionDescendant(m_cell.m_volumeProportionDescendant),
m_garbageProportionDescendant(m_cell.m_garbageProportionDescendant),
m_cellCycleSizer(m_cell.m_cellCycleSizer),
m_volumeCellMitosis(m_cell.m_volumeCellMitosis),
m_volumeMutationRate(m_cell.m_volumeMutationRate),
m_sizerConcentration(m_cell.m_sizerConcentration),
m_sizerThresholdMitosis(m_cell.m_sizerThresholdMitosis),
m_sizeMinimalGenome(m_cell.m_sizeMinimalGenome),
m_cellCycle(m_cell.m_cellCycle),

m_physiologicalAlleles(m_cell.m_physiologicalAlleles),
m_physiologicalFunctions(m_cell.m_physiologicalFunctions)

{
    //gsl_rng *rnd=gsl_rng_alloc(gsl_rng_default);
    Random=m_cell.Random;
    
    int numberTFgenes= (int) m_cell.m_regulatoryAlleles.size();
    int totalGeneNumber= numberTFgenes+ (int) m_cell.m_physiologicalAlleles.size();
    
    
    m_numberTFTranscripts=new int[numberTFgenes];
    m_numberTranscriptionFactors=new int[numberTFgenes];
    m_totalBindingSites=new int[numberTFgenes];
    m_connectivityMatrix=new int*[totalGeneNumber];
    for (int Gene=0;Gene<totalGeneNumber;Gene++){
        m_connectivityMatrix[Gene]=new int[numberTFgenes];
    }
    
    for (int g=0;g<numberTFgenes;g++){
        m_numberTFTranscripts[g]=m_cell.m_numberTFTranscripts[g];
        m_numberTranscriptionFactors[g]=m_cell.m_numberTranscriptionFactors[g];
        m_totalBindingSites[g]=m_cell.m_totalBindingSites[g];
    }
    for (int Gene=0;Gene<totalGeneNumber;Gene++){
        for (int typeOfTFNumber=0;typeOfTFNumber<numberTFgenes;typeOfTFNumber++){
            m_connectivityMatrix[Gene][typeOfTFNumber]=m_cell.m_connectivityMatrix[Gene][typeOfTFNumber];
        }
    }
    for (int typeOfTFNumber=0;typeOfTFNumber<numberTFgenes;typeOfTFNumber++){
        m_totalBindingSites[typeOfTFNumber]=m_cell.m_totalBindingSites[typeOfTFNumber];
    }
    for (int typeTFnumber=0;typeTFnumber<numberTFgenes;typeTFnumber++){
        m_regulatoryAlleles.push_back(m_cell.m_regulatoryAlleles[typeTFnumber]);
    }
}


void Cell::operator=(const Cell &m_cell)
{
    Random=m_cell.Random;
    
    m_generationNumber=(m_cell.m_generationNumber);
    
    m_K__on=(m_cell.m_K__on);                                                                       //Binding reaction constant
    m_K__off=(m_cell.m_K__off);                                                                     //Dissciation reaction constant
    m_permeabilityCell=(m_cell.m_permeabilityCell);                                                 //Permeability of the cell membrane
    m_allometryProductionAlpha=(m_cell.m_allometryProductionAlpha);                                 //Allometry of energy production: how scale the production of energy per nutrient unit with cell volume
    m_allometryBasalAlpha=(m_cell.m_allometryBasalAlpha);                                           //Allometry of basal metabolism: how scale the cost of basal metabolism with cell volume
    m_etaNutrientToEnergy=(m_cell.m_etaNutrientToEnergy);                                           //Conversion rate between a nutrient and energy
    m_etaNutrientToGarbage=(m_cell.m_etaNutrientToGarbage);                                         //Conversion rate between a nutrient and a garbage
    m_metaboliteCost=(m_cell.m_metaboliteCost);       //Cost for :{"Growth","Maintenance","Public goods"}
    m_synthesisDuration=(m_cell.m_synthesisDuration);
    
    //Variable:
    //ii.Transcription rate of genes                                                               //Further, these rates will rely on BSTF dynamics
    m_tauTranscriptionOncePolymeraseBound=
        (m_cell.m_tauTranscriptionOncePolymeraseBound);
    m_basalBindingProbabilityPolymerase=
        (m_cell.m_basalBindingProbabilityPolymerase);
    m_enhancedBindingProbabilityPolymerase=
        (m_cell.m_enhancedBindingProbabilityPolymerase);
    //iii.Traduction rate of transcripts
    m_tauTraduction=(m_cell.m_tauTraduction);
    m_etaTranslationActivityPostCheckpoint=(m_cell.m_etaTranslationActivityPostCheckpoint);
    //iv.Degradation rate of gene products
    m_tauDegradationTranscript=(m_cell.m_tauDegradationTranscript);                                 //Same degradation rate for each transcript
    m_tauDegradationProtein=(m_cell.m_tauDegradationProtein);                                       //Same degradation rate for each protein
    //v.Metabolite dynamics rates
    m_metaboliteMaintenanceProductionRate=(m_cell.m_metaboliteMaintenanceProductionRate);
    m_metaboliteCommonGoodsProductionRate=(m_cell.m_metaboliteCommonGoodsProductionRate);
    m_metaboliteGrowthProductionRate=(m_cell.m_metaboliteGrowthProductionRate);
    m_metaboliteDegradationRate=(m_cell.m_metaboliteDegradationRate);
    //vi.Size of the cell
    m_volumeCell=(m_cell.m_volumeCell);
    m_timeSinceSynthesisBeginning=(m_cell.m_timeSinceSynthesisBeginning);
    
    /*2.Content of the cell*/
    
    //o.Calibration of basal metabolism constraints
    m_SA_calibration=(m_cell.m_SA_calibration);
    m_basalMetabolismSAcalibration=(m_cell.m_basalMetabolismSAcalibration);
    m_currentBasalMetabolism=(m_cell.m_currentBasalMetabolism);
    //i.Energetic content
    m_numberStandardNutrientUnits=(m_cell.m_numberStandardNutrientUnits);                           //In a first time, nutrient units are immediately consumed
    m_numberRevealedNutrientUnits=(m_cell.m_numberRevealedNutrientUnits);
    //m_basalEnergyUnits=(m_cell.m_basalEnergyUnits);
    m_numberEnergyUnits=(m_cell.m_numberEnergyUnits);                                               //In a first time, energy units are immediately consumed
    //ii.By-product content: garbage
    m_numberGarbageUnits=(m_cell.m_numberGarbageUnits);
    //iii.Transcripts content
    m_numberGrowthTranscripts=(m_cell.m_numberGrowthTranscripts);
    m_numberMaintenanceTranscripts=(m_cell.m_numberMaintenanceTranscripts);
    m_numberPublicGoodsTranscripts=(m_cell.m_numberPublicGoodsTranscripts);
    m_numberSizerTranscripts=(m_cell.m_numberSizerTranscripts);
    //iv.Proteins content
    m_numberGrowthProteins=(m_cell.m_numberGrowthProteins);
    m_numberMaintenanceProteins=(m_cell.m_numberMaintenanceProteins);
    m_numberPublicGoodsProteins=(m_cell.m_numberPublicGoodsProteins);
    m_numberSizerProteins=(m_cell.m_numberSizerProteins);
    //v.Maintenance metabolites content
    m_numberMaintenanceMetabolites=(m_cell.m_numberMaintenanceMetabolites);
    m_numberPublicGoodsMetabolites=(m_cell.m_numberPublicGoodsMetabolites);
    
    /*3.Genetics*/
    
    //i.Cell volume at mitosis
    m_volumeProportionDescendant=(m_cell.m_volumeProportionDescendant);
    m_garbageProportionDescendant=(m_cell.m_garbageProportionDescendant);
    m_cellCycleSizer=(m_cell.m_cellCycleSizer);
    m_volumeCellMitosis=(m_cell.m_volumeCellMitosis);
    m_volumeMutationRate=(m_cell.m_volumeMutationRate);
    m_sizerConcentration=(m_cell.m_sizerConcentration);
    m_sizerThresholdMitosis=(m_cell.m_sizerThresholdMitosis);
    m_sizeMinimalGenome=(m_cell.m_sizeMinimalGenome);
    
    //ii.Genes
    m_cellCycle=(m_cell.m_cellCycle);
    
    //PrintCell();
    
    int numberTFgenes= (int) m_cell.m_regulatoryAlleles.size();
    int totalGeneNumber= numberTFgenes + (int) m_physiologicalAlleles.size();
    
    if (m_cell.m_physiologicalAlleles.size()!=m_physiologicalAlleles.size()){
        for (int Gene=0;Gene<totalGeneNumber;Gene++){
            delete[] m_connectivityMatrix[Gene];
        }
        delete[] m_connectivityMatrix;
        
        totalGeneNumber= numberTFgenes + (int) m_cell.m_physiologicalAlleles.size();
        m_connectivityMatrix=new int*[totalGeneNumber];
        for (int Gene=0;Gene<totalGeneNumber;Gene++){
            m_connectivityMatrix[Gene]=new int[numberTFgenes];
        }
   }
    
    
    m_physiologicalAlleles=(m_cell.m_physiologicalAlleles);
    m_physiologicalFunctions=(m_cell.m_physiologicalFunctions);
    //GetKey_valuesMapofPhysiologicalAllele(m_physiologicalFunctions, m_physiologicalAlleles);
    
    for (int g=0;g<numberTFgenes;g++){
        m_numberTFTranscripts[g]=m_cell.m_numberTFTranscripts[g];
        m_numberTranscriptionFactors[g]=m_cell.m_numberTranscriptionFactors[g];
        m_totalBindingSites[g]=m_cell.m_totalBindingSites[g];
    }
    int compteur=0;
    
    for (int Gene=0;Gene<totalGeneNumber;Gene++){
        for (int typeOfTFNumber=0;typeOfTFNumber<numberTFgenes;typeOfTFNumber++){
            compteur+=1;
            
            m_connectivityMatrix[Gene][typeOfTFNumber]=m_cell.m_connectivityMatrix[Gene][typeOfTFNumber];
        }
    }
    for (int typeOfTFNumber=0;typeOfTFNumber<numberTFgenes;typeOfTFNumber++){
        m_totalBindingSites[typeOfTFNumber]=m_cell.m_totalBindingSites[typeOfTFNumber];
    }
    for (int typeTFnumber=0;typeTFnumber<numberTFgenes;typeTFnumber++){
        m_regulatoryAlleles[typeTFnumber]=m_cell.m_regulatoryAlleles[typeTFnumber];
    }
}

Cell::~Cell()
{
    delete[] m_numberTFTranscripts;
    delete[] m_numberTranscriptionFactors;
    delete[] m_totalBindingSites;
    
    const int totalGeneNumber= (int) m_regulatoryAlleles.size() + (int) m_physiologicalAlleles.size();
    for (int Gene=0;Gene<totalGeneNumber;Gene++){
        delete[] m_connectivityMatrix[Gene];
    }
    delete[] m_connectivityMatrix;
}

void Cell::BuildDaughterCell(Cell& _cell)
{
    m_generationNumber+=1;
    
    int numberTFgenes= (int) _cell.m_regulatoryAlleles.size();
    m_volumeCell=SAtoVconversion(VtoSAconversion(_cell.m_volumeCell)*m_volumeProportionDescendant);
    
    _cell.m_volumeCell=SAtoVconversion(VtoSAconversion(_cell.m_volumeCell)*(1-m_volumeProportionDescendant));
    
    m_timeSinceSynthesisBeginning=0;
    _cell.m_timeSinceSynthesisBeginning=0;
    
    /*2.Content of the cell*/
    
    //i.Energetic content
    m_numberStandardNutrientUnits=0;                           //In a first time, nutrient units are immediately consumed
    _cell.m_numberStandardNutrientUnits=0;
    m_numberRevealedNutrientUnits=0;                           //In a first time, nutrient units are immediately consumed
    _cell.m_numberRevealedNutrientUnits=0;
    m_numberEnergyUnits=0;                             //In a first time, energy units are immediately consumed
    //ii.By-product content: garbage
    m_numberGarbageUnits=gsl_ran_binomial(Random,m_volumeProportionDescendant,(int)_cell.m_numberGarbageUnits);
    _cell.m_numberGarbageUnits-=m_numberGarbageUnits;
    //iii.Transcripts content
    m_numberGrowthTranscripts=gsl_ran_binomial(Random,m_volumeProportionDescendant,_cell.m_numberGrowthTranscripts);
    _cell.m_numberGrowthTranscripts-=m_numberGrowthTranscripts;
    m_numberMaintenanceTranscripts=gsl_ran_binomial(Random,m_volumeProportionDescendant,_cell.m_numberMaintenanceTranscripts);
    _cell.m_numberMaintenanceTranscripts-=m_numberMaintenanceTranscripts;
    m_numberPublicGoodsTranscripts=gsl_ran_binomial(Random,m_volumeProportionDescendant,_cell.m_numberPublicGoodsTranscripts);
    _cell.m_numberPublicGoodsTranscripts-=m_numberPublicGoodsTranscripts;
    m_numberSizerTranscripts=gsl_ran_binomial(Random,m_volumeProportionDescendant,_cell.m_numberSizerTranscripts);
    _cell.m_numberSizerTranscripts-=m_numberSizerTranscripts;
    for (int TFnumber=0;TFnumber<numberTFgenes;TFnumber++){
        m_numberTFTranscripts[TFnumber]=gsl_ran_binomial(Random,m_volumeProportionDescendant,_cell.m_numberTFTranscripts[TFnumber]);
        _cell.m_numberTFTranscripts[TFnumber]-=m_numberTFTranscripts[TFnumber];
    }
    //iv.Proteins content
    m_numberGrowthProteins=gsl_ran_binomial(Random,m_volumeProportionDescendant,_cell.m_numberGrowthProteins);
    _cell.m_numberGrowthProteins-=m_numberGrowthProteins;
    m_numberMaintenanceProteins=gsl_ran_binomial(Random,m_volumeProportionDescendant,_cell.m_numberMaintenanceProteins);
    _cell.m_numberMaintenanceProteins-=m_numberMaintenanceProteins;
    m_numberPublicGoodsProteins=gsl_ran_binomial(Random,m_volumeProportionDescendant,_cell.m_numberPublicGoodsProteins);
    _cell.m_numberPublicGoodsProteins-=m_numberPublicGoodsProteins;
    m_numberSizerProteins=gsl_ran_binomial(Random,m_volumeProportionDescendant,_cell.m_numberSizerProteins);
    _cell.m_numberSizerProteins-=m_numberSizerProteins;
    for (int TFnumber=0;TFnumber<numberTFgenes;TFnumber++){
        m_numberTranscriptionFactors[TFnumber]=gsl_ran_binomial(Random,m_volumeProportionDescendant,_cell.m_numberTranscriptionFactors[TFnumber]);
            _cell.m_numberTranscriptionFactors[TFnumber]-=m_numberTranscriptionFactors[TFnumber];
    }
    //v.Maintenance metabolites content
    m_numberMaintenanceMetabolites=gsl_ran_binomial(Random,m_volumeProportionDescendant,_cell.m_numberMaintenanceMetabolites);
    _cell.m_numberMaintenanceMetabolites-=m_numberMaintenanceMetabolites;
    m_numberPublicGoodsMetabolites=gsl_ran_binomial(Random,m_volumeProportionDescendant,_cell.m_numberPublicGoodsMetabolites);
    _cell.m_numberPublicGoodsMetabolites-=m_numberPublicGoodsMetabolites;
    
    /*3.Genetics*/
    
    //i.Cell volume at mitosis
    m_volumeCellMitosis=(_cell.m_volumeCellMitosis);
    m_sizerThresholdMitosis=(_cell.m_sizerThresholdMitosis);
    
}

/* ##Other methods## */


    // #Printing the state of the cell# //

std::string Cell::printCellMatrix(){
    stringstream tmp("");
    for (int typeOfTFNumber=0;typeOfTFNumber<m_regulatoryAlleles.size();typeOfTFNumber++){
        for (int gene=0;gene<(m_regulatoryAlleles.size()+m_physiologicalAlleles.size());gene++){
            tmp << m_connectivityMatrix[gene][typeOfTFNumber];
            tmp << "_";
        }
    }
    const std::string tmpstr = tmp.str() ;
    const std::size_t sz = tmpstr.size() ;
    return (tmpstr.substr(0,sz-1));
}

void Cell::PrintCell(){
    
    cout << "Generation number of the cell: " << m_generationNumber << endl;
    cout << "BS number replic: " << m_sizerThresholdMitosis << endl;
    //cout << "Time since synthesis began: " << m_timeSinceSynthesisBeginning << endl;
    cout << "Growth proteins: " << m_numberGrowthProteins << endl;
    cout << "Maintenance proteins: " << m_numberMaintenanceProteins << endl;
    cout << "Public Goods proteins: " << m_numberPublicGoodsProteins << endl;
    cout << "Sizer proteins: " << m_numberSizerProteins << endl;
    cout << "Volume at mitosis: " << m_volumeCellMitosis << endl;
    cout << "Volume: " << m_volumeCell << endl;
    cout << "Garbage: " << m_numberGarbageUnits << endl;
    cout <<"Maintenance metabolites: " << m_numberMaintenanceMetabolites << endl;
    cout << "Growth transcripts:" << m_numberGrowthTranscripts << endl;
    cout << "Maintenance transcripts: " << m_numberMaintenanceTranscripts << endl;
    cout << "Public Goods transcripts: " << m_numberPublicGoodsTranscripts << endl;
    cout << "Sizer transcripts: " << m_numberSizerTranscripts << endl;
    cout << "Connectivity matrix:" << endl;
    for (int typeOfTFNumber=0;typeOfTFNumber<m_regulatoryAlleles.size();typeOfTFNumber++){
        for (int gene=0;gene<(m_regulatoryAlleles.size()+m_physiologicalAlleles.size());gene++){
            cout << m_connectivityMatrix[gene][typeOfTFNumber] << " ";
        }
        cout << endl;
    }
}

    /* #Operations on the cell# */

        /* #0.Binding site transcription factor connectivity evaluation# */
/*
void Cell::GeneEvaluateConnectivity(char geneType, int ) //method evaluating the connevtivity between one specific gene regulatory sequence and any transcription factor
*/

   
void Cell::EvaluateBSTFConnectivity(){//As a by-product, the method calculate the total number of binding-sites for each transcription factor, which is stored in the constant vector m_totalBindingSites to avoid too many calculus at each timestep
    int _TFGeneNumber=(int) m_regulatoryAlleles.size();
    int _PhysioGenenumber=(int) m_physiologicalAlleles.size();
    int totalGeneNumber= _TFGeneNumber + _PhysioGenenumber;
    for (int typeOfTFNumber=0;typeOfTFNumber<_TFGeneNumber;typeOfTFNumber++){
        for (int physiologicalGeneNumber=0;physiologicalGeneNumber<m_physiologicalAlleles.size();physiologicalGeneNumber++){
            m_connectivityMatrix[physiologicalGeneNumber][typeOfTFNumber]= m_physiologicalAlleles[m_physiologicalFunctions[physiologicalGeneNumber]].BSTFMatchCount(m_regulatoryAlleles[typeOfTFNumber]);
        }
        for (int typeOfTFNumber_tested=0;typeOfTFNumber_tested<_TFGeneNumber;typeOfTFNumber_tested++){
            m_connectivityMatrix[_PhysioGenenumber+typeOfTFNumber_tested][typeOfTFNumber]=m_regulatoryAlleles[typeOfTFNumber_tested].BSTFMatchCount(m_regulatoryAlleles[typeOfTFNumber]);
        }
    }
    for (int typeofTFNumber=0;typeofTFNumber<_TFGeneNumber;typeofTFNumber++){
        for (int gene_number=0;gene_number<totalGeneNumber;gene_number++){
            m_totalBindingSites[typeofTFNumber]+=m_connectivityMatrix[gene_number][typeofTFNumber];
        }
    }
}

void Cell::ReevaluateBSTFConnectivityAfterFactorSequenceMutation(int transcriptionFactorNumber){
    int _TFGeneNumber=(int) m_regulatoryAlleles.size();
    int _PhysioGenenumber=(int) m_physiologicalAlleles.size();
    int totalGeneNumber= _TFGeneNumber + _PhysioGenenumber;
    
    for (int physiologicalGeneNumber=0;physiologicalGeneNumber<_PhysioGenenumber;physiologicalGeneNumber++){
        m_connectivityMatrix[physiologicalGeneNumber][transcriptionFactorNumber]= m_physiologicalAlleles[m_physiologicalFunctions[physiologicalGeneNumber]].BSTFMatchCount(m_regulatoryAlleles[transcriptionFactorNumber]);
    }
    
    for (int typeOfTFNumber_tested=0;typeOfTFNumber_tested<_TFGeneNumber;typeOfTFNumber_tested++){
        m_connectivityMatrix[_PhysioGenenumber+typeOfTFNumber_tested][transcriptionFactorNumber]=m_regulatoryAlleles[typeOfTFNumber_tested].BSTFMatchCount(m_regulatoryAlleles[transcriptionFactorNumber]);
    }
    
    m_totalBindingSites[transcriptionFactorNumber]=0;
    for (int gene_number=0;gene_number<totalGeneNumber;gene_number++){
        m_totalBindingSites[transcriptionFactorNumber]+=m_connectivityMatrix[gene_number][transcriptionFactorNumber];
    }
}

void Cell::ReevaluateBSTFConnectivityAfterRegulatoryFactorSequenceMutation(int transcriptionFactorNumber){
    int _TFGeneNumber=(int) m_regulatoryAlleles.size();
    int _PhysioGenenumber=(int) m_physiologicalAlleles.size();
    int totalGeneNumber= _TFGeneNumber + _PhysioGenenumber;
    
    for (int typeOfTFNumber=0;typeOfTFNumber<_TFGeneNumber;typeOfTFNumber++){
        m_connectivityMatrix[_PhysioGenenumber+transcriptionFactorNumber][typeOfTFNumber] =m_regulatoryAlleles[transcriptionFactorNumber].BSTFMatchCount(m_regulatoryAlleles[typeOfTFNumber]);
    }
    for (int typeofTFNumber=0;typeofTFNumber<_TFGeneNumber;typeofTFNumber++){
        m_totalBindingSites[typeofTFNumber]=0;
        for (int gene_number=0;gene_number<totalGeneNumber;gene_number++){
            m_totalBindingSites[typeofTFNumber]+=m_connectivityMatrix[gene_number][typeofTFNumber];
        }
    }
}

void Cell::ReevaluateBSTFConnectivityAfterRegulatoryPhysiologicalSequenceMutation(string function){
    int _TFGeneNumber=(int) m_regulatoryAlleles.size();
    int _PhysioGenenumber=(int) m_physiologicalAlleles.size();
    int totalGeneNumber= _TFGeneNumber + _PhysioGenenumber;
    
    auto find1 = find(begin(m_physiologicalFunctions), end(m_physiologicalFunctions), function);
    int functionPosition=(int) distance(begin(m_physiologicalFunctions),find1);
    for (int typeOfTFNumber=0;typeOfTFNumber<_TFGeneNumber;typeOfTFNumber++){
        m_connectivityMatrix[functionPosition][typeOfTFNumber]=m_physiologicalAlleles[function].BSTFMatchCount(m_regulatoryAlleles[typeOfTFNumber]);
    }
    for (int typeofTFNumber=0;typeofTFNumber<_TFGeneNumber;typeofTFNumber++){
        m_totalBindingSites[typeofTFNumber]=0;
        for (int gene_number=0;gene_number<totalGeneNumber;gene_number++){
            m_totalBindingSites[typeofTFNumber]+=m_connectivityMatrix[gene_number][typeofTFNumber];
        }
    }
}

void Cell::ReconstructBSTFConnectivityMatrixAfterEfficiencyPhysiologicalMutation(string function){
    
    //cout << "taill vector fonction2: " << m_physiologicalFunctions.size() << endl;
    //cout << "taill vector alleles2: " << m_physiologicalAlleles.size() << endl;
    int _PhysioGenenumber=(int) m_physiologicalAlleles.size();
    int _TFGeneNumber=(int) m_regulatoryAlleles.size();
    int totalGeneNumber= _TFGeneNumber + _PhysioGenenumber;
    for (int Gene=0;Gene<totalGeneNumber;Gene++){
        delete[] m_connectivityMatrix[Gene];
    }
    delete[] m_connectivityMatrix;
    
    
    auto find1 = find(begin(m_physiologicalFunctions), end(m_physiologicalFunctions), function);
    int functionPosition=(int) distance(begin(m_physiologicalFunctions),find1);
    //cout << function << endl;
    m_physiologicalAlleles.erase(function);
    m_physiologicalFunctions.erase(m_physiologicalFunctions.begin()+functionPosition);
    
    if (m_physiologicalAlleles.size()!=m_physiologicalFunctions.size()){
        cout << "Problem!!!!!" << endl;
    }
    
    _PhysioGenenumber=(int) m_physiologicalAlleles.size();
    _TFGeneNumber=(int) m_regulatoryAlleles.size();
    totalGeneNumber= _TFGeneNumber + _PhysioGenenumber;
    m_connectivityMatrix=new int*[totalGeneNumber];
    for (int Gene=0;Gene<totalGeneNumber;Gene++){
        m_connectivityMatrix[Gene]=new int[_TFGeneNumber];
    }
    
    //cout << "taill vector fonction: " << m_physiologicalFunctions.size() << endl;
    //cout << "taill vector alleles: " << m_physiologicalAlleles.size() << endl;
    for (int typeOfTFNumber=0;typeOfTFNumber<_TFGeneNumber;typeOfTFNumber++){
        for (int physiologicalGeneNumber=0;physiologicalGeneNumber<m_physiologicalAlleles.size();physiologicalGeneNumber++){
            m_connectivityMatrix[physiologicalGeneNumber][typeOfTFNumber]= m_physiologicalAlleles[m_physiologicalFunctions[physiologicalGeneNumber]].BSTFMatchCount(m_regulatoryAlleles[typeOfTFNumber]);
        }
        for (int typeOfTFNumber_tested=0;typeOfTFNumber_tested<_TFGeneNumber;typeOfTFNumber_tested++){
            m_connectivityMatrix[_PhysioGenenumber+typeOfTFNumber_tested][typeOfTFNumber]=m_regulatoryAlleles[typeOfTFNumber_tested].BSTFMatchCount(m_regulatoryAlleles[typeOfTFNumber]);
        }
    }
    for (int typeofTFNumber=0;typeofTFNumber<_TFGeneNumber;typeofTFNumber++){
        for (int gene_number=0;gene_number<totalGeneNumber;gene_number++){
            m_totalBindingSites[typeofTFNumber]+=m_connectivityMatrix[gene_number][typeofTFNumber];
        }
    }
}

        /* #1.Transcripts Production# */

void Cell::TranscriptsDynamics(double delta_t){
    int _TFGeneNumber=(int) m_regulatoryAlleles.size();
    int _PhysioGenenumber=(int) m_physiologicalAlleles.size();
    
    //o.Calculus of the BS and TF total concentrations in the cell, i.e. including those bound and those free
    double cBStot[_TFGeneNumber];
    double cTFtot[_TFGeneNumber];
    for (int TF_number=0;TF_number<_TFGeneNumber;TF_number++){
        if(m_totalBindingSites[TF_number]>0){
            cBStot[TF_number]=ConcentrationCalculus(m_totalBindingSites[TF_number], m_volumeCell);
        }
        else{
            cBStot[TF_number]=0;
        }
        if(m_numberTranscriptionFactors[TF_number]>0){
            cTFtot[TF_number]=ConcentrationCalculus(m_numberTranscriptionFactors[TF_number], m_volumeCell);
        }
        else{
            cTFtot[TF_number]=0;
        }
    }
    //i.Calculus of the equilibrium concentration of each transcription factor - binding site complex, keeping in mind that a transcription factor can only recognise one binding-site: the case in which two transcription factors have become identical is not considered for the moment but has to be in order to test the importance of the phenomenon.
    double cBSTFequilibrium[_TFGeneNumber];
    for (int TF_number=0;TF_number<_TFGeneNumber;TF_number++){
        if (cTFtot[TF_number]>0){
        cBSTFequilibrium[TF_number]=
            (cBStot[TF_number]+cTFtot[TF_number])/2+
        (m_K__off-pow(pow(m_K__off+m_K__on*(cBStot[TF_number]+cTFtot[TF_number]),2)-4*pow(m_K__on,2)*cBStot[TF_number]*cTFtot[TF_number],0.5))/(2*m_K__on);
        }
        else{
            cBSTFequilibrium[TF_number]=0;
        }
    }
    //cout << "Concentration TF1: " <<cBSTFequilibrium[0]<< endl;
    double* probability_nonBSTF=new double[_TFGeneNumber];
    for (int TF_number=0;TF_number<_TFGeneNumber;TF_number++){
        probability_nonBSTF[TF_number]=1-cBSTFequilibrium[TF_number]/cBStot[TF_number];
    }
    //cout << "prob TF1: " <<probability_nonBSTF[0]<< endl;
    //iii.Calculus of the transcription rate of each gene
    //#Initialization#//
    map <string,double> tauTranscriptionPhysiologicalGenes;
    
    for (int geneNumber=0;geneNumber<_PhysioGenenumber;geneNumber++){
        tauTranscriptionPhysiologicalGenes.insert(pair<string,double>(m_physiologicalFunctions[geneNumber],m_tauTranscriptionOncePolymeraseBound));
    }
    double tauTranscriptionFactorsGeneTranscription[_TFGeneNumber];
    //#Calculus#//
    
    double bindingPolymeraseProbabilityWhileNoHelp[_PhysioGenenumber];
    double positiveTranscriptionFactorBindingProbability[_PhysioGenenumber];
    double bindingPolymeraseProbabilityWhileHelpAndNoRepression[_PhysioGenenumber];
    
    double TFbindingPolymeraseProbabilityWhileNoHelp[_TFGeneNumber];
    double TFpositiveTranscriptionFactorBindingProbability[_TFGeneNumber];
    double TFbindingPolymeraseProbabilityWhileHelpAndNoRepression[_TFGeneNumber];
    
        //a.Physiological genes transcription
    for (int geneNumber=0;geneNumber<_PhysioGenenumber;geneNumber++){
    bindingPolymeraseProbabilityWhileNoHelp[geneNumber]=m_basalBindingProbabilityPolymerase*WeightProduct(_TFGeneNumber,probability_nonBSTF,m_connectivityMatrix[geneNumber]);
    
        positiveTranscriptionFactorBindingProbability[geneNumber]=1;
        bindingPolymeraseProbabilityWhileHelpAndNoRepression[geneNumber]=m_enhancedBindingProbabilityPolymerase;
            for (int TF_number=0;TF_number<_TFGeneNumber;TF_number++){
                if(m_regulatoryAlleles[TF_number].NatureOfEffect()=='A'){
                    positiveTranscriptionFactorBindingProbability[geneNumber]*=pow((probability_nonBSTF[TF_number]),m_connectivityMatrix[geneNumber][TF_number]); //Probability of at least one activator bound to one of each BS
                }
                else if(m_regulatoryAlleles[TF_number].NatureOfEffect()=='R'){
                    bindingPolymeraseProbabilityWhileHelpAndNoRepression[geneNumber]*=(pow(probability_nonBSTF[TF_number],m_connectivityMatrix[geneNumber][TF_number]));
                }
            }
        bindingPolymeraseProbabilityWhileHelpAndNoRepression[geneNumber]*=1-positiveTranscriptionFactorBindingProbability[geneNumber];
        
        if (m_physiologicalFunctions[geneNumber]=="Cell sizer" && ReplicationCheckpoint()){
            tauTranscriptionPhysiologicalGenes[m_physiologicalFunctions[geneNumber]]=0;     //Redundant because all transrciption processes should be stopped after synthesis started
        }
        else{
            tauTranscriptionPhysiologicalGenes[m_physiologicalFunctions[geneNumber]]*=
            (bindingPolymeraseProbabilityWhileNoHelp[geneNumber]+bindingPolymeraseProbabilityWhileHelpAndNoRepression[geneNumber]);
        }
    }
    //cout << "bindingprob sans aide: "<<positiveTranscriptionFactorBindingProbability[1] << endl;
    //cout << "prob bind sans aide: " << bindingPolymeraseProbabilityWhileNoHelp[1] << endl;
    //cout << "prob bind avec aide: " << bindingPolymeraseProbabilityWhileHelpAndNoRepression[1] << endl;
    
        //e.Transcription factors transcription
    
    for (int TF_gene_number=0;TF_gene_number<_TFGeneNumber;TF_gene_number++){
        tauTranscriptionFactorsGeneTranscription[TF_gene_number]=m_tauTranscriptionOncePolymeraseBound;
        TFpositiveTranscriptionFactorBindingProbability[TF_gene_number]=1;
        TFbindingPolymeraseProbabilityWhileNoHelp[TF_gene_number]= m_basalBindingProbabilityPolymerase*WeightProduct(_TFGeneNumber,probability_nonBSTF,m_connectivityMatrix[_PhysioGenenumber+TF_gene_number]);
    
        TFbindingPolymeraseProbabilityWhileHelpAndNoRepression[TF_gene_number]=m_enhancedBindingProbabilityPolymerase;
        for (int TF_number=0;TF_number<_TFGeneNumber;TF_number++){
            if(m_regulatoryAlleles[TF_number].NatureOfEffect()=='A'){
                TFpositiveTranscriptionFactorBindingProbability[TF_gene_number]*=
                    pow((probability_nonBSTF[TF_number]),m_connectivityMatrix[_PhysioGenenumber+TF_gene_number][TF_number]);
            }
            else if(m_regulatoryAlleles[TF_number].NatureOfEffect()=='R'){
                TFbindingPolymeraseProbabilityWhileHelpAndNoRepression[TF_gene_number]*=
                    pow(probability_nonBSTF[TF_number],m_connectivityMatrix[_PhysioGenenumber+TF_gene_number][TF_number]);
            }
        }
        TFbindingPolymeraseProbabilityWhileHelpAndNoRepression[TF_gene_number]*=1-TFpositiveTranscriptionFactorBindingProbability[TF_gene_number];
        
        tauTranscriptionFactorsGeneTranscription[TF_gene_number]*=(TFbindingPolymeraseProbabilityWhileNoHelp[TF_gene_number]+TFbindingPolymeraseProbabilityWhileHelpAndNoRepression[TF_gene_number]);
    }
    //cout << tauTranscriptionPhysiologicalGenes["Maintenance"] << endl;
    //iv.Production and degradation of transcripts during a timestep
    
    if (ReplicationCheckpoint()){
        if(m_numberGrowthTranscripts>0){
            m_numberGrowthTranscripts-=gsl_ran_binomial(Random,m_tauDegradationTranscript*delta_t,m_numberGrowthTranscripts);
        }
        
        if(m_numberMaintenanceTranscripts>0){
            m_numberMaintenanceTranscripts-=gsl_ran_binomial(Random,m_tauDegradationTranscript*delta_t,m_numberMaintenanceTranscripts);
        }
        
        if(m_numberPublicGoodsTranscripts>0){
            m_numberPublicGoodsTranscripts-=gsl_ran_binomial(Random,m_tauDegradationTranscript*delta_t,m_numberPublicGoodsTranscripts);
        }
        
        if(m_numberSizerTranscripts>0){
            m_numberSizerTranscripts-=gsl_ran_binomial(Random,m_tauDegradationTranscript*delta_t,m_numberSizerTranscripts);
        }
        for (int TF_number=0;TF_number<_TFGeneNumber;TF_number++){
            if(m_numberTFTranscripts[TF_number]>0){
                m_numberTFTranscripts[TF_number]-=gsl_ran_binomial(Random,m_tauDegradationTranscript*delta_t,m_numberTFTranscripts[TF_number]);
            }
            
        }
    }
    else{
        if(m_numberGrowthTranscripts>0){
            m_numberGrowthTranscripts-=gsl_ran_binomial(Random,m_tauDegradationTranscript*delta_t,m_numberGrowthTranscripts);
        }
        
        m_numberGrowthTranscripts+=gsl_ran_poisson(Random,tauTranscriptionPhysiologicalGenes["Growth"]*delta_t);
    
        if(m_numberMaintenanceTranscripts>0){
            m_numberMaintenanceTranscripts-=gsl_ran_binomial(Random,m_tauDegradationTranscript*delta_t,m_numberMaintenanceTranscripts);
        }
        
        m_numberMaintenanceTranscripts+=gsl_ran_poisson(Random,tauTranscriptionPhysiologicalGenes["Maintenance"]*delta_t);
    
        if(m_numberPublicGoodsTranscripts>0){
            m_numberPublicGoodsTranscripts-=gsl_ran_binomial(Random,m_tauDegradationTranscript*delta_t,m_numberPublicGoodsTranscripts);
        }
        
        m_numberPublicGoodsTranscripts+=gsl_ran_poisson(Random,tauTranscriptionPhysiologicalGenes["Public Goods"]*delta_t);
        
        if (m_cellCycleSizer){
            if(m_numberSizerTranscripts>0){
                m_numberSizerTranscripts-=gsl_ran_binomial(Random,m_tauDegradationTranscript*delta_t,m_numberSizerTranscripts);
            }
            m_numberSizerTranscripts+=gsl_ran_poisson(Random,tauTranscriptionPhysiologicalGenes["Cell sizer"]*delta_t);
            }
    
        for (int TF_number=0;TF_number<_TFGeneNumber;TF_number++){
            if(m_numberTFTranscripts[TF_number]>0){
                m_numberTFTranscripts[TF_number]-=gsl_ran_binomial(Random,m_tauDegradationTranscript*delta_t,m_numberTFTranscripts[TF_number]);
            }
            
            m_numberTFTranscripts[TF_number]+=gsl_ran_poisson(Random,tauTranscriptionFactorsGeneTranscription[TF_number]*delta_t);
        }
    }
    delete[] probability_nonBSTF;
}

        /* #2.Proteins Production# */

void Cell::ProteinsDynamics(double delta_t){
        double concomitantRibosomeTranslating=5;
        map <string, map <string, double > > tauTraduction;
        tauTraduction["Non cycle-dependent Genes"]["G"]=m_tauTraduction;
        tauTraduction["Non cycle-dependent Genes"]["SM"]=m_tauTraduction*m_etaTranslationActivityPostCheckpoint;
        if (m_physiologicalAlleles.count("Cell sizer")>0){
        tauTraduction["Cell Sizer"]["G"]=m_tauTraduction;
        tauTraduction["Cell Sizer"]["SM"]=0;
        }
        string cyclePhase;
        if (ReplicationCheckpoint()){
            cyclePhase="SM";
        }
        else{
            cyclePhase="G";
        }
        if(m_numberGrowthProteins>0){
            m_numberGrowthProteins-=gsl_ran_binomial(Random,m_tauDegradationProtein*delta_t,m_numberGrowthProteins);
        }
    
        if(m_numberGrowthTranscripts>0){
            m_numberGrowthProteins+=gsl_ran_poisson(Random,tauTraduction["Non cycle-dependent Genes"][cyclePhase]*concomitantRibosomeTranslating*delta_t*m_numberGrowthTranscripts);
        }
    
        //ii.Maintenance proteins dynamics
        if(m_numberMaintenanceProteins>0){
            m_numberMaintenanceProteins-=gsl_ran_binomial(Random,m_tauDegradationProtein*delta_t,m_numberMaintenanceProteins);
        }
    
        if(m_numberMaintenanceTranscripts>0){
            m_numberMaintenanceProteins+=gsl_ran_poisson(Random,tauTraduction["Non cycle-dependent Genes"][cyclePhase]*concomitantRibosomeTranslating*delta_t*m_numberMaintenanceTranscripts);
        }
    
        //iii.Public goods proteins dynamics
        if(m_numberPublicGoodsProteins>0){
            m_numberPublicGoodsProteins-=gsl_ran_binomial(Random,m_tauDegradationProtein*delta_t,m_numberPublicGoodsProteins);
        }
    
        if(m_numberPublicGoodsTranscripts>0){
            m_numberPublicGoodsProteins+=gsl_ran_poisson(Random,tauTraduction["Non cycle-dependent Genes"][cyclePhase]*concomitantRibosomeTranslating*delta_t*m_numberPublicGoodsTranscripts);
        }
    
        //iv.Sizer proteins dynamics
        //While synthesis, translation still goes on, but while G2 ansd mitosis, translation decreases a lot and eventually stops: considering that synthesis time encompasses those three steps and that G2/mitosis lasts as long as synthesis (according to litterature, reasonable statement), explaining the ratio of 2. Notice that here, we have assumed that growh only occurs in a G1 manner, assuming the one growth step model does not distort evolution outcomes.
        if (m_cellCycleSizer){
            double tauConsumptionSizer=0.005;
            if(m_numberSizerProteins>0){
                m_numberSizerProteins-=gsl_ran_binomial(Random,(m_tauDegradationProtein+tauConsumptionSizer)*delta_t,m_numberSizerProteins); //Constant consumption such as in the Eukaryotic system cyclin/CDK
            }
    
            if(m_numberSizerTranscripts>0){
                m_numberSizerProteins+=gsl_ran_poisson(Random,tauTraduction["Cell Sizer"][cyclePhase]*concomitantRibosomeTranslating*delta_t*m_numberSizerTranscripts*m_currentBasalMetabolism/(1.5*pow(10,8)));
            }
        }
        //v.Transcription factors dynamics
        int numberTFgenes= (int) m_regulatoryAlleles.size();
        
        for (int TFnumber=0;TFnumber<numberTFgenes;TFnumber++){
            if(m_numberTranscriptionFactors[TFnumber]>0){
                m_numberTranscriptionFactors[TFnumber]-=gsl_ran_binomial(Random,m_tauDegradationProtein*delta_t,m_numberTranscriptionFactors[TFnumber]);
            }
            
            if(m_numberTFTranscripts[TFnumber]>0){
                m_numberTranscriptionFactors[TFnumber]+=gsl_ran_poisson(Random,tauTraduction["Non cycle-dependent Genes"][cyclePhase]*concomitantRibosomeTranslating*delta_t*m_numberTFTranscripts[TFnumber]);
            }
            
        }
    }


        /* #3.Nutrient uptake# */

void Cell::SetNutrientsAfterNutrientDiffusion(double totalMembraneSurface, double totalStandardNutrientUptake,double totalRevealedNutrientUptake){
    m_numberStandardNutrientUnits+= totalStandardNutrientUptake*((VtoSAconversion(m_volumeCell)*pow(10,-10))/totalMembraneSurface);
    m_numberRevealedNutrientUnits+= totalRevealedNutrientUptake*((VtoSAconversion(m_volumeCell)*pow(10,-10))/totalMembraneSurface);
    }

        /* 3bis. Testing possible expenses */
bool Cell::EnergyBudgetUsefulnessTest(double delta_t){
    double totalBudget=pow(m_etaNutrientToEnergy,m_allometryProductionAlpha)*(m_numberStandardNutrientUnits+m_numberRevealedNutrientUnits);
    //cout<< "Energy: " << totalBudget << endl;
    BasalMetabolismCalculus(delta_t);
    //cout <<"basal: " << m_currentBasalMetabolism;
    double maximumConsumption=(m_numberGrowthProteins*m_metaboliteCost["Growth"]+(m_numberMaintenanceProteins*m_metaboliteCost["Maintenance"]+m_numberPublicGoodsProteins*m_metaboliteCost["Public Goods"]))*delta_t+m_currentBasalMetabolism;
    return (maximumConsumption>totalBudget);
}

        /* #4.Energy Production and collateral Garbage production# */
void Cell::EnergyProduction(bool energyWholeBudgetUsefulness,double delta_t){
    if (energyWholeBudgetUsefulness){
        if(pow(m_etaNutrientToEnergy,m_allometryProductionAlpha)*(m_numberStandardNutrientUnits+m_numberRevealedNutrientUnits)>0){
                        m_numberGarbageUnits+= pow(m_etaNutrientToGarbage,m_allometryProductionAlpha)*(m_numberStandardNutrientUnits+m_numberRevealedNutrientUnits);
            //Note: Here, there's an allometric coefficient representing the mitochondrial efficiency in the nutrient treatment which may also be linked to their concentration.
            m_numberEnergyUnits+=pow(m_etaNutrientToEnergy,m_allometryProductionAlpha)*(m_numberStandardNutrientUnits+m_numberRevealedNutrientUnits);   //Same allometric coefficient.
        }
    m_numberStandardNutrientUnits=0;
    m_numberRevealedNutrientUnits=0;
    }
    else{
        double maximumConsumption=(m_numberGrowthProteins*m_metaboliteCost["Growth"]+m_numberMaintenanceProteins*m_metaboliteCost["Maintenance"]+ m_numberPublicGoodsProteins*m_metaboliteCost["Public Goods"])*delta_t+m_currentBasalMetabolism;
        m_numberEnergyUnits+=maximumConsumption;
        m_numberGarbageUnits+=maximumConsumption/m_etaNutrientToEnergy*m_etaNutrientToGarbage;
        double currentDepletedNutrient=maximumConsumption/30;
        double proportionStandardNutrient=m_numberStandardNutrientUnits/(m_numberStandardNutrientUnits+m_numberRevealedNutrientUnits);
        m_numberStandardNutrientUnits-=currentDepletedNutrient*proportionStandardNutrient;
        m_numberRevealedNutrientUnits-=currentDepletedNutrient*(1-proportionStandardNutrient);
    }
}

        /* #5.Energy use# */

            // #i.Withdrawal by the basal metabolism
void Cell::BasalMetabolismCalculus(double delta_t){
    m_currentBasalMetabolism=m_basalMetabolismSAcalibration*delta_t*pow(VtoSAconversion(m_volumeCell)/m_SA_calibration,m_allometryBasalAlpha);
}

void Cell::BasalMetabolismWithdrawal(){
    if (m_numberEnergyUnits>0){
        if (m_currentBasalMetabolism < m_numberEnergyUnits){
            m_numberEnergyUnits-=m_currentBasalMetabolism;
            //m_basalEnergyUnits+=m_basalMetabolismSAcalibration*delta_t*pow(VtoSAconversion(m_volumeCell)/SAcalibration,m_allometryBasalAlpha);
            //cout << "Unités d'énergie:  " <<  m_basalEnergyUnits << endl;
        }
        else{
            m_numberEnergyUnits=0;
        }
    }
}

                // #ii.a)Cell growth
void Cell::Growth(int totalNumberProteins,bool energyWholeBudgetUsefulness,double delta_t){
    if (totalNumberProteins>0){
        double growthFractionProteins=Proportion(m_numberGrowthProteins,totalNumberProteins);
        int numberPhospholipids(0);
        if (energyWholeBudgetUsefulness){
            numberPhospholipids=(int) ((m_numberEnergyUnits*m_metaboliteGrowthProductionRate/m_metaboliteCost["Growth"])*growthFractionProteins);
        }
        else{
            numberPhospholipids=(int) (m_numberGrowthProteins*m_metaboliteGrowthProductionRate*delta_t);
        }
        if (numberPhospholipids>0){
                m_volumeCell=SAtoVconversion(VtoSAconversion(m_volumeCell)+numberPhospholipids*2*pow(10,-7));
        }
    }
}


                // #ii.b)Cell maintenance molecule production
void Cell::MaintenanceMetabolitesProduction(int totalNumberProteins,bool energyWholeBudgetUsefulness,double delta_t){
    m_numberMaintenanceMetabolites-=gsl_ran_binomial(Random,m_tauDegradationProtein*delta_t,m_numberMaintenanceMetabolites);
    if(totalNumberProteins>0){
        double maintenanceFractionProteins=Proportion(m_numberMaintenanceProteins,totalNumberProteins);
        if (energyWholeBudgetUsefulness){
            //cout << "Metabolic Production Rate: " << m_metaboliteGenericProductionRate << endl;
                    m_numberMaintenanceMetabolites+=(m_numberEnergyUnits*m_metaboliteMaintenanceProductionRate/m_metaboliteCost["Maintenance"])*maintenanceFractionProteins;
            }
        else{
            m_numberMaintenanceMetabolites+=(m_numberMaintenanceProteins*m_metaboliteMaintenanceProductionRate);
        }
    }
}

                // #ii.c)Cell Public Goods production
void Cell::PublicGoodsMetabolitesProduction(int totalNumberProteins,bool energyWholeBudgetUsefulness,double delta_t){
    if(totalNumberProteins>0){
        double publicGoodsFractionProteins=Proportion(m_numberPublicGoodsProteins,totalNumberProteins);
        if (energyWholeBudgetUsefulness){
            m_numberPublicGoodsMetabolites+=(m_numberEnergyUnits*m_metaboliteCommonGoodsProductionRate/m_metaboliteCost["Public Goods"])*publicGoodsFractionProteins;
            }
        else{
            m_numberPublicGoodsMetabolites+=(m_numberPublicGoodsProteins*m_metaboliteCommonGoodsProductionRate);
        }
    }
}

            // #ii)Whole energy allocation process summarized
void Cell::EnergyBalance(bool energyWholeBudgetUsefulness,double delta_t){
    //i)Energy production
    EnergyProduction(energyWholeBudgetUsefulness,delta_t);
    //ii)Basal metabolism withdrawal
    BasalMetabolismWithdrawal();
    //iii)Energy allocation
    int totalNumberTraitProteins=m_numberGrowthProteins+m_numberMaintenanceProteins+m_numberPublicGoodsProteins;
    if (m_timeSinceSynthesisBeginning==0 || m_cellCycle=="G1SG2M"){
    Growth(totalNumberTraitProteins,energyWholeBudgetUsefulness,delta_t);
    }
    MaintenanceMetabolitesProduction(totalNumberTraitProteins,energyWholeBudgetUsefulness,delta_t);
    PublicGoodsMetabolitesProduction(totalNumberTraitProteins,energyWholeBudgetUsefulness,delta_t);
    m_numberEnergyUnits=0;
}

        /* #6.Phenotype effects */

            // #i.Maintenance metabolites effect# //
//Analytical expression should by far be preferred in order to minimize modelling time.
void Cell::GarbageNeutralization(double K_garbageNeutralization,double delta_t){
    double _AvogadroConstant=6.02*pow(10.0,23);
    double timeStepNeutralizedGarbage(0);
    double maintenanceMetabolitesConcentration=ConcentrationCalculus(m_numberMaintenanceMetabolites/_AvogadroConstant, m_volumeCell/pow(10,15));
    double garbageUnitsConcentration=ConcentrationCalculus(m_numberGarbageUnits/_AvogadroConstant, m_volumeCell/pow(10,15));
    //cout << "Garbage concentration: " << garbageUnitsConcentration-maintenanceMetabolitesConcentration << endl;
    m_numberMaintenanceMetabolites-=m_numberMaintenanceMetabolites*m_metaboliteDegradationRate;
    if (m_numberMaintenanceMetabolites!=0 && m_numberGarbageUnits!=0){
            timeStepNeutralizedGarbage = maintenanceMetabolitesConcentration*garbageUnitsConcentration*(exp(K_garbageNeutralization*delta_t*(garbageUnitsConcentration-maintenanceMetabolitesConcentration))-1)/(garbageUnitsConcentration*exp(K_garbageNeutralization*delta_t*(garbageUnitsConcentration-maintenanceMetabolitesConcentration))-maintenanceMetabolitesConcentration);
        if (m_numberMaintenanceMetabolites>timeStepNeutralizedGarbage)
            m_numberMaintenanceMetabolites-=timeStepNeutralizedGarbage*m_volumeCell*_AvogadroConstant/pow(10,15);
        else{
            m_numberMaintenanceMetabolites=0;
        }
        if (m_numberMaintenanceMetabolites<0){
            m_numberMaintenanceMetabolites=0;
        }
    
        if (m_numberGarbageUnits>timeStepNeutralizedGarbage){
            m_numberGarbageUnits-=timeStepNeutralizedGarbage*m_volumeCell*_AvogadroConstant/pow(10,15);}
        else{
            m_numberGarbageUnits=0;
        }
        if (m_numberGarbageUnits<0){
            m_numberGarbageUnits=0;
        }
    }
    }

            // #ii.Public goods release# //
double Cell::PublicgoodsRelease(){
    double numberPublicGoodsReleased=m_numberPublicGoodsMetabolites;
    m_numberPublicGoodsMetabolites=0;
    return(numberPublicGoodsReleased);
}

            // #iii.Nutrient release if remainder# //

double Cell::StandardNutrientRelease(){
    double numberExcessStandardNutrient=m_numberStandardNutrientUnits;
    m_numberStandardNutrientUnits=0;
    return(numberExcessStandardNutrient);
}

double Cell::RevealedNutrientRelease(){
    double numberExcessRevealedNutrient=m_numberRevealedNutrientUnits;
    m_numberRevealedNutrientUnits=0;
    return(numberExcessRevealedNutrient);
}

        /* #7.Mitosis checking */

            // #i.Checking if mitosis checkpoint reached#
bool Cell::ReplicationCheckpoint(){
    if (m_cellCycleSizer){
        if(m_sizerConcentration){
            double sizerConcTempo=ConcentrationCalculus(m_numberSizerProteins, m_volumeCell);
            if (m_timeSinceSynthesisBeginning>0){
                if (m_timeSinceSynthesisBeginning==1 && m_cellCycle=="GSM"){
                    m_volumeCellMitosis=m_volumeCell;
                }
                return(true);
            }
            else{
                return(sizerConcTempo>m_sizerThresholdMitosis);
            }
        }
        else{
            if (m_timeSinceSynthesisBeginning>0){
                if (m_timeSinceSynthesisBeginning==1 && m_cellCycle=="GSM"){
                    m_volumeCellMitosis=m_volumeCell;
                }
                return(true);
            }
            else{
                return(m_numberSizerProteins>m_sizerThresholdMitosis);
            }
        }
    }
    else{
        return(m_volumeCell>m_volumeCellMitosis);
    }
}

            // #ii.Processing of DNA synthesis#
void Cell::SynthesisProcessing(double delta_t){
    if (ReplicationCheckpoint()){
        m_timeSinceSynthesisBeginning+=delta_t;
    }
}

            // #iii.Checking if mitosis checkpoint reached#
bool Cell::MitosisCheckpoint(){
    if (m_cellCycleSizer){
        if(ReplicationCheckpoint()){
            double sizerCurrentConcentration=ConcentrationCalculus(m_numberSizerProteins, m_volumeCell);
            return (m_timeSinceSynthesisBeginning>m_synthesisDuration && sizerCurrentConcentration<m_sizerThresholdMitosis*pow(10,-3));
        }
        else{
            return(false);
        }
    }
    else{
        if(ReplicationCheckpoint()){
            return (m_timeSinceSynthesisBeginning>m_synthesisDuration);
        }
        else{
            return(false);
        }
    }
}

        /* #8.Death checking */
bool Cell::DeathTest(double delta_t){
    long double concentrationGarbageUnits=ConcentrationCalculus(m_numberGarbageUnits,m_volumeCell);
    double alpha=pow(10,5);
    double beta=5;
    double basalDeathRate=0.05;
    //cout << 1/(1+pow((concentrationGarbageUnits/alpha),-beta)) << endl;
    //cout << (1.0-(basalDeathRate+(1-basalDeathRate)/(1+pow((concentrationGarbageUnits/alpha),-beta)))) << endl;
    double survivalProbability(1);
    
    if ( (1.0-(basalDeathRate+(1-basalDeathRate)/(1+pow((concentrationGarbageUnits/alpha),-beta)))) < 0 ){
        survivalProbability=0;
    }
    else if ( (1.0-(basalDeathRate+(1-basalDeathRate)/(1+pow((concentrationGarbageUnits/alpha),-beta)))) > 1 ){
        survivalProbability=1;
    }
    else{
        survivalProbability=pow(1.0-(basalDeathRate+(1-basalDeathRate)/(1+pow((concentrationGarbageUnits/alpha),-beta))),delta_t/3600);
        //cout << "Survie: " << survivalProbability<< endl;
    }
    int deathTest=gsl_ran_binomial(Random, survivalProbability, 1);
    return(deathTest==0);
}

        /* #9.Mutation processing */
void Cell::MutationProcessing(){
    int _TFGeneNumber=(int) m_regulatoryAlleles.size();
    int _PhysioGenenumber=(int) m_physiologicalAlleles.size();
    // #Regulatory sequences and binding-site recognition sequences# //
    if (_PhysioGenenumber>0){
        for (int _physiologicalGeneNumber=0;_physiologicalGeneNumber<_PhysioGenenumber;_physiologicalGeneNumber++){
            //cout << m_physiologicalFunctions[_physiologicalGeneNumber] << endl;
            if(m_physiologicalAlleles[m_physiologicalFunctions[_physiologicalGeneNumber]].EfficiencyMutation()){
                //cout << m_physiologicalFunctions[_physiologicalGeneNumber] << endl;
                ReconstructBSTFConnectivityMatrixAfterEfficiencyPhysiologicalMutation(m_physiologicalFunctions[_physiologicalGeneNumber]);
                //PrintCell();
                _physiologicalGeneNumber-=1;
                _PhysioGenenumber-=1;
            }
            else if(m_physiologicalAlleles[m_physiologicalFunctions[_physiologicalGeneNumber]].RegulatorySequenceMutation()){
                ReevaluateBSTFConnectivityAfterRegulatoryPhysiologicalSequenceMutation(m_physiologicalFunctions[_physiologicalGeneNumber]);
            }
        }
    }
    for (int _regulatoryGeneNumber=0;_regulatoryGeneNumber<_TFGeneNumber;_regulatoryGeneNumber++){
        if (m_regulatoryAlleles[_regulatoryGeneNumber].BindingSiteRecognitionMutation()){
            ReevaluateBSTFConnectivityAfterFactorSequenceMutation(_regulatoryGeneNumber);
        }
        if(m_regulatoryAlleles[_regulatoryGeneNumber].RegulatorySequenceMutation()){
            ReevaluateBSTFConnectivityAfterRegulatoryFactorSequenceMutation(_regulatoryGeneNumber);
        }
    }
    
    if (m_cellCycleSizer==false){
        int volumeMitosisMutation;
        volumeMitosisMutation=gsl_ran_binomial(Random,m_volumeMutationRate, 1);
        if (volumeMitosisMutation==1){
            m_volumeCellMitosis+=gsl_ran_gaussian(Random,2);
        }
    }
    else{
        int sizerThresholdMutation;
        sizerThresholdMutation=gsl_ran_binomial(Random, m_volumeMutationRate, 1);
        int positiveMutation;
        positiveMutation=gsl_ran_binomial(Random, 0.5, 1)*sizerThresholdMutation;
        if (positiveMutation==1){
            m_sizerThresholdMitosis+=gsl_ran_poisson(Random, 2);
        }
        else{
            if(sizerThresholdMutation==1){
                m_sizerThresholdMitosis-=gsl_ran_poisson(Random, 2);
            }
        }
    }
}

bool Cell::MitosisFailure(){
    return(m_volumeCell<m_sizeMinimalGenome);
}

    /* #Getting the settings# */

double Cell::GetPermeabilityCell(){
    return(m_permeabilityCell);
}

double Cell::GetAllometryProductionAlpha(){
    return(m_allometryProductionAlpha);
}

double Cell::GetAllometryBasalAlpha(){
    return(m_allometryBasalAlpha);
}

double Cell::GetEtaNutrientToGarbage(){
    return(m_etaNutrientToGarbage);
}

double Cell::GetTauTranscriptionOncePolymeraseBound(){
    return(m_tauTranscriptionOncePolymeraseBound);
}

double Cell::GetBasalBindingProbabilityPolymerase(){
    return(m_basalBindingProbabilityPolymerase);
}

double Cell::GetEnhancedBindingProbabilityPolymerase(){
    return(m_enhancedBindingProbabilityPolymerase);
}

double Cell::GetTauDegradationTranscript(){
    return(m_tauDegradationTranscript);
}

double Cell::GetTauDegradationProtein(){
    return(m_tauDegradationProtein);
}

int Cell::GetMetaboliteCost(string physiologicalFunction){
    return(m_metaboliteCost[physiologicalFunction]);
}

double Cell::GetVolumeCell(){
    return(m_volumeCell);
}

double Cell::GetSACell(){
    return(VtoSAconversion(m_volumeCell));
}

long int Cell::GetNumberGarbageUnits(){
    return(m_numberGarbageUnits);
}

int Cell::GetTranscripts(char t){
    if (t=='g'){
        return(m_numberGrowthTranscripts);
    }
    else if (t=='m'){
        return(m_numberMaintenanceTranscripts);
    }
    else if (t=='p'){
        return(m_numberPublicGoodsTranscripts);
    }
    else if (t=='s'){
        return(m_numberSizerTranscripts);
    }
    else{
        return(0);
    }
}

int Cell::GetTFTranscripts(int g){
        return(m_numberTFTranscripts[g]);
}

int Cell::GetProteins(char t){
    if (t=='g'){
        return(m_numberGrowthProteins);
    }
    else if (t=='m'){
        return(m_numberMaintenanceProteins);
    }
    else if (t=='p'){
        return(m_numberPublicGoodsProteins);
    }
    else if (t=='s'){
        return(m_numberSizerProteins);
    }
    else{
        return(0);
    }
}

int Cell::GetTranscriptionFactors(int g){
        return(m_numberTranscriptionFactors[g]);
}

int Cell::GetGenerationNumber(){
    return(m_generationNumber);
}

//May be provisory getting the size of the cell at mitosis

double Cell::GetMitosisSize(){
    return(m_volumeCellMitosis);
}

double Cell::GetMitosisThreshold(){
    return(m_sizerThresholdMitosis);
}

void Cell::SetGarbageUnits(double numberGarbageunits){
    m_numberGarbageUnits+=numberGarbageunits;
}

//Connectivity matrix conversion
void Cell::MatrixConversion(vector < int > &matrixCode){
    for (int gene=0; gene< m_regulatoryAlleles.size()+m_physiologicalFunctions.size();gene++){
        for (int TF=0; TF<m_regulatoryAlleles.size();TF++){
            matrixCode.push_back(m_connectivityMatrix[gene][TF]);
        }
    }
}

void Cell::GeneTabularConversion(string &geneCode){
    int count=0;
    while (count<m_physiologicalFunctions.size()){
        geneCode+=m_physiologicalFunctions[count];
        count+=1;
    }
}

//General methods out of the class

double VtoSAconversion(double V_cell){
    return pow((V_cell*V_cell*36*M_PI),1.0/3);
}

double SAtoVconversion(double SA_cell){
    return pow(SA_cell,3.0/2)/(6*pow(M_PI,1.0/2));
}
           
long double ConcentrationCalculus(double MoleculeNumber, double CellVolume){
    return(MoleculeNumber/CellVolume);
}
           
double Proportion(double sub, double total){
    return(sub/total);
}

double WeightProduct(const int nb_elements, double* elements,int* repeat_numbers){
    int element_number(0);
    double product(1);
        while (element_number<nb_elements){
            product*=pow(elements[element_number],repeat_numbers[element_number]);
            element_number++;
        }
    return (product);
}

void GetKey_valuesMapofPhysiologicalAllele(vector<string> &keys, map<string,PhysiologicalAllele> Map){
    for (map<string,PhysiologicalAllele>::iterator it=Map.begin();it!=Map.end();it++){
        string _keyTested=it->first;
        int _vectorElementlementNumber(0);
        if (keys.size()==0){
            keys.push_back(_keyTested);
        }
        else{
            while (_vectorElementlementNumber<keys.size()){
                if (keys[_vectorElementlementNumber]==_keyTested){
                    break;
                }
                else if(_vectorElementlementNumber<keys.size()-1){
                    _vectorElementlementNumber+=1;
                }
                else{
                    keys.push_back(_keyTested);
                }
            }
        }
    }
}

