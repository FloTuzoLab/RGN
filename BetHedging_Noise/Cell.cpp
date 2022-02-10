//
//  Cell.cpp
//  BetHedging_Noise
//
//  Created by Florian on 28/05/2019.
//  Copyright © 2019 Florian. All rights reserved.
//

#include "Cell.hpp"

double const AvogadroConstant=6.02*pow(10.0,23); //unit: molecules number
double const UptakeConversion=AvogadroConstant*pow(10,-15);//uptake in mol.L^-1 translated into molecules.micrometer^-3

/*  ##Constructors##    */

Cell::Cell(){
    
}

/* Constructor mainly defined from the outside */

Cell::Cell(int nutrientTypes,
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
           map <string,PhysiologicalAllele> physiologicalAlleles,
           vector <RegulatoryAllele> positiveRegulatoryAlleles,
           vector <RegulatoryAllele> negativeRegulatoryAlleles,
           gsl_rng *random
           ):
/* 0.Generaion number */
m_generationNumber(0),

/*1.Properties of the cell*/

//Fixed:
//i.Chemico-physical properties
m_k__cat(kcat),                                                     //Binding reaction constant
m_K__M(KM),                                                   //Dissociation reaction constant
m_K__on(Kon),
m_K__off(Koff),
m_permeabilityCell(permeability),                                   //Permeability of the cell membrane
//Allometry of basal metabolism: how scale the cost of basal metabolism with cell volume

//Variable:
//ii.Transcription rate of genes                                    //Further, these rates will rely on BSTF dynamics
m_tauTranscriptionOncePolymeraseBound(tauTranscriptionOncePolymeraseBound),
m_basalBindingProbabilityPolymerase(basalBindingProbabilityPolymerase),
m_enhancedBindingProbabilityPolymerase(enhancedBindingProbabilityPolymerase),
//iii.Traduction rate of transcripts
m_tauTraduction(tauTraduction),
//iv.Degradation rate of gene products
m_tauDegradationTranscript(tauDegradationTranscript),               //Same degradation rate for each transcript
m_tauDegradationProtein(tauDegradationProtein),                     //Same degradation rate for each protein
//m_enzymeProductionRate(enzymeProductionRate),
//vi.Size of the cell
m_volumeCell(cellVolume),

/*2.Content of the cell*/
//o.Calibration of basal metabolism constraints
m_currentBasalMetabolism(0),
m_synthesisRateProp(0),
m_ribosomesAmount(2*pow(10,5)),
m_mRNA_ribosomes(pow(10,2)),
//i.Energetic content                                        //In a first time, nutrient units are immediately consumed
m_numberNutrientTypes(nutrientTypes),
//m_basalEnergyUnits(0),
m_numberEnergyUnits(0),                                             //In a first time, energy units are immediately consumed
//ii.By-product content: garbage
//iii.Transcripts content
m_numberNut1Transcripts(50),
m_numberNut2Transcripts(50),
//iv.Proteins content
m_numberNut1Proteins(100000),
m_numberNut2Proteins(100000),
/*//v.Maintenance metabolites content
m_numberNut1Enzymes(0),
m_numberNut2Enzymes(0),*/

/*3.Genetics*/

//i.Cell volume at mitosis
m_volumeRatioDescendant(volumeRatioDescendant),
m_volumeCellMitosis(volumeCellMitosis),
m_volumeMutationRate(volumeMutationRate)
//ii.Genes
{
    //gsl_rng *rnd=gsl_rng_alloc(gsl_rng_default);
    _Random=random;
    
    int numberTFgenes= (int) positiveRegulatoryAlleles.size()+ (int) negativeRegulatoryAlleles.size();
    m_numberTFTranscripts=new int[numberTFgenes];
    m_numberTranscriptionFactors=new int[numberTFgenes];
    m_totalBindingSites=new int[numberTFgenes];
    for (int g=0;g<numberTFgenes;g++){
        m_numberTFTranscripts[g]=0;//numberOtherTranscripts[g];
        m_numberTranscriptionFactors[g]=0;//numberOtherProteins[g];
        m_totalBindingSites[g]=0;
    }
    
    m_numberNutUnits=new int[m_numberNutrientTypes];
    for (int _nutType=0;_nutType<m_numberNutrientTypes;_nutType++){
        m_numberNutUnits[_nutType]=0;
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
}

// #Constructor using a full pre-defined cell# //
/*Cell(const Cell &m_cell, bool cellSizerCycle,  gsl_rng *randomSimulation); useless for BH model */
Cell::Cell(const Cell &_cell):
    /* 0.Generaion number */
    m_generationNumber(_cell.m_generationNumber),
    
    /*1.Properties of the cell*/
    
    //Fixed:
    //i.Chemico-physical properties
    m_k__cat(_cell.m_k__cat),                                                                       //Binding reaction constant
    m_K__M(_cell.m_K__M),
    m_K__on(_cell.m_K__on),
    m_K__off(_cell.m_K__off),

    //Dissciation reaction constant
    m_permeabilityCell(_cell.m_permeabilityCell),                                                 //Permeability of the cell membrane
    
    //Variable:
    //ii.Transcription rate of genes                                                               //Further, these rates will rely on BSTF dynamics
    m_tauTranscriptionOncePolymeraseBound(_cell.m_tauTranscriptionOncePolymeraseBound),
    m_basalBindingProbabilityPolymerase(_cell.m_basalBindingProbabilityPolymerase),
    m_enhancedBindingProbabilityPolymerase(_cell.m_enhancedBindingProbabilityPolymerase),
    //iii.Traduction rate of transcripts
    m_tauTraduction(_cell.m_tauTraduction),
    //iv.Degradation rate of gene products
    m_tauDegradationTranscript(_cell.m_tauDegradationTranscript),                                 //Same degradation rate for each transcript
    m_tauDegradationProtein(_cell.m_tauDegradationProtein),                                       //Same degradation rate for each protein
    //v.Metabolite dynamics rates
    //m_enzymeProductionRate(_cell.m_enzymeProductionRate),
    //vi.Size of the cell
    m_volumeCell(_cell.m_volumeCell),
    
    /*2.Content of the cell*/

    m_currentBasalMetabolism(_cell.m_currentBasalMetabolism),
    m_synthesisRateProp(_cell.m_synthesisRateProp),
    m_ribosomesAmount(_cell.m_ribosomesAmount),
    m_mRNA_ribosomes(_cell.m_mRNA_ribosomes),
    //i.Energetic content //In a first time, nutrient units are immediately consumed
    //m_basalEnergyUnits(m_cell.m_basalEnergyUnits),
    m_numberNutrientTypes(_cell.m_numberNutrientTypes),
    m_numberEnergyUnits(_cell.m_numberEnergyUnits),                                               //In a first time, energy units are immediately consumed
    //iii.Transcripts content
    m_numberNut1Transcripts(_cell.m_numberNut1Transcripts),
    m_numberNut2Transcripts(_cell.m_numberNut2Transcripts),
    //iv.Proteins content
    m_numberNut1Proteins(_cell.m_numberNut1Proteins),
    m_numberNut2Proteins(_cell.m_numberNut2Proteins),
    
    /*3.Genetics*/
    
    //i.Cell volume at mitosis
    m_volumeRatioDescendant(_cell.m_volumeRatioDescendant),
    m_volumeCellMitosis(_cell.m_volumeCellMitosis),
    m_volumeMutationRate(_cell.m_volumeMutationRate),
    
    m_physiologicalAlleles(_cell.m_physiologicalAlleles),
    m_physiologicalFunctions(_cell.m_physiologicalFunctions)
    
    {
        //gsl_rng *rnd=gsl_rng_alloc(gsl_rng_default);
        _Random=_cell._Random;
        
        int numberTFgenes= (int) _cell.m_regulatoryAlleles.size();
        int totalGeneNumber= numberTFgenes+ (int) _cell.m_physiologicalAlleles.size();
        
        m_numberNutUnits=new int[m_numberNutrientTypes];
        for (int _nutType=0;_nutType<m_numberNutrientTypes;_nutType++){
            m_numberNutUnits[_nutType]=_cell.m_numberNutUnits[_nutType];
        }
        
        m_numberTFTranscripts=new int[numberTFgenes];
        m_numberTranscriptionFactors=new int[numberTFgenes];
        m_totalBindingSites=new int[numberTFgenes];
        m_connectivityMatrix=new int*[totalGeneNumber];
        for (int Gene=0;Gene<totalGeneNumber;Gene++){
            m_connectivityMatrix[Gene]=new int[numberTFgenes];
        }
        
        for (int g=0;g<numberTFgenes;g++){
            m_numberTFTranscripts[g]=_cell.m_numberTFTranscripts[g];
            m_numberTranscriptionFactors[g]=_cell.m_numberTranscriptionFactors[g];
            m_totalBindingSites[g]=_cell.m_totalBindingSites[g];
        }
        for (int Gene=0;Gene<totalGeneNumber;Gene++){
            for (int typeOfTFNumber=0;typeOfTFNumber<numberTFgenes;typeOfTFNumber++){
                m_connectivityMatrix[Gene][typeOfTFNumber]=_cell.m_connectivityMatrix[Gene][typeOfTFNumber];
            }
        }
        for (int typeOfTFNumber=0;typeOfTFNumber<numberTFgenes;typeOfTFNumber++){
            m_totalBindingSites[typeOfTFNumber]=_cell.m_totalBindingSites[typeOfTFNumber];
        }
        for (int typeTFnumber=0;typeTFnumber<numberTFgenes;typeTFnumber++){
            m_regulatoryAlleles.push_back(_cell.m_regulatoryAlleles[typeTFnumber]);
        }
}

void Cell::operator=(const Cell &m_cell)
{
    _Random=m_cell._Random;
    
    m_generationNumber=(m_cell.m_generationNumber);
    //Speed reaction constant
    m_k__cat=(m_cell.m_k__cat);
    m_K__M=(m_cell.m_K__M);
    m_K__on=(m_cell.m_K__on);
    m_K__off=(m_cell.m_K__off);
    
    m_permeabilityCell=(m_cell.m_permeabilityCell);                                                 //Permeability of the cell membrane
    
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
    //iv.Degradation rate of gene products
    m_tauDegradationTranscript=(m_cell.m_tauDegradationTranscript);                                 //Same degradation rate for each transcript
    m_tauDegradationProtein=(m_cell.m_tauDegradationProtein);                                       //Same degradation rate for each protein
    //v.Metabolite dynamics rates
    //m_enzymeProductionRate=(m_cell.m_enzymeProductionRate);
    
    //vi.Size of the cell
    m_volumeCell=(m_cell.m_volumeCell);
    
    /*2.Content of the cell*/
    
    //o.Calibration of basal metabolism constraints
    m_currentBasalMetabolism=(m_cell.m_currentBasalMetabolism);
    m_synthesisRateProp=(m_cell.m_synthesisRateProp);
    m_ribosomesAmount=(m_cell.m_ribosomesAmount);
    m_mRNA_ribosomes=(m_cell.m_mRNA_ribosomes);
    //i.Energetic content
    //m_basalEnergyUnits=(m_cell.m_basalEnergyUnits);
    m_numberNutrientTypes=(m_cell.m_numberNutrientTypes);
    m_numberEnergyUnits=(m_cell.m_numberEnergyUnits);                                               //In a first time, energy units are immediately consumed
    //iii.Transcripts content
    m_numberNut1Transcripts=(m_cell.m_numberNut1Transcripts);
    m_numberNut2Transcripts=(m_cell.m_numberNut2Transcripts);
    //iv.Proteins content
    m_numberNut1Proteins=(m_cell.m_numberNut1Proteins);
    m_numberNut2Proteins=(m_cell.m_numberNut2Proteins);
    /*//v.Maintenance enzymes content
    m_numberNut1enzymes=(m_cell.m_numberNut1enzymes);
    m_numberNut2enzymes=(m_cell.m_numberNut2enzymes);*/
    
    /*3.Genetics*/
    
    //i.Cell volume at mitosis
    m_volumeRatioDescendant=(m_cell.m_volumeRatioDescendant);
    m_volumeCellMitosis=(m_cell.m_volumeCellMitosis);
    m_volumeMutationRate=(m_cell.m_volumeMutationRate);
    
    
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
    //m_numberNutUnits=new int[m_numberNutrientTypes];
    for (int _nutType=0;_nutType<m_numberNutrientTypes;_nutType++){
        m_numberNutUnits[_nutType]=0;
    }
    
    m_physiologicalAlleles=(m_cell.m_physiologicalAlleles);
    m_physiologicalFunctions=(m_cell.m_physiologicalFunctions);
    //GetKey_valuesMapofPhysiologicalAllele(m_physiologicalFunctions, m_physiologicalAlleles);
    
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
        m_regulatoryAlleles[typeTFnumber]=m_cell.m_regulatoryAlleles[typeTFnumber];
    }
}

void Cell::BuildDaughterCell(Cell& _motherCell)
{
    m_generationNumber+=1;
    _motherCell.m_generationNumber+=1;
    int numberTFgenes= (int) _motherCell.m_regulatoryAlleles.size();
    m_volumeCell=pow((1/(1+pow(m_volumeRatioDescendant,2.0/3.0))),3.0/2.0)*_motherCell.m_volumeCell;
    
    _motherCell.m_volumeCell=pow((pow(m_volumeRatioDescendant,2.0/3.0)/(1+pow(m_volumeRatioDescendant,2.0/3.0))),3.0/2.0)*_motherCell.m_volumeCell;
    //std::cout <<m_volumeCell << endl;
    //std::cout << _motherCell.m_volumeCell <<endl;
    
    /*2.Content of the cell*/
    double motherCellProp=m_volumeRatioDescendant/(1+m_volumeRatioDescendant);
    double daughterCellProp=1-motherCellProp;
    
    //i.Energetic content
    for (int _nutType=0;_nutType<m_numberNutrientTypes;_nutType++){
        m_numberNutUnits[_nutType]=daughterCellProp*_motherCell.m_numberNutUnits[_nutType];
                        //In a first time, nutrient units are immediately consumed
        _motherCell.m_numberNutUnits[_nutType]=motherCellProp*_motherCell.m_numberNutUnits[_nutType];
    }
    
    m_numberEnergyUnits=daughterCellProp*_motherCell.m_numberEnergyUnits;                             //In a first time, energy units are immediately consumed
    _motherCell.m_numberEnergyUnits=motherCellProp*_motherCell.m_numberEnergyUnits;
    
    //iii.Transcripts content
    m_numberNut1Transcripts=gsl_ran_binomial(_Random,daughterCellProp,_motherCell.m_numberNut1Transcripts);
    _motherCell.m_numberNut1Transcripts-=m_numberNut1Transcripts;
    
    m_numberNut2Transcripts=gsl_ran_binomial(_Random,daughterCellProp,_motherCell.m_numberNut2Transcripts);
    _motherCell.m_numberNut2Transcripts-=m_numberNut2Transcripts;
    
    for (int TFnumber=0;TFnumber<numberTFgenes;TFnumber++){
        m_numberTFTranscripts[TFnumber]=gsl_ran_binomial(_Random,daughterCellProp,_motherCell.m_numberTFTranscripts[TFnumber]);
        _motherCell.m_numberTFTranscripts[TFnumber]-=m_numberTFTranscripts[TFnumber];
    }
    //iv.Proteins content
    m_numberNut1Proteins=gsl_ran_binomial(_Random,daughterCellProp,_motherCell.m_numberNut1Proteins);
    _motherCell.m_numberNut1Proteins-=m_numberNut1Proteins;
    
    m_numberNut2Proteins=gsl_ran_binomial(_Random,daughterCellProp,_motherCell.m_numberNut2Proteins);
    _motherCell.m_numberNut2Proteins-=m_numberNut2Proteins;
    
    for (int TFnumber=0;TFnumber<numberTFgenes;TFnumber++){
        m_numberTranscriptionFactors[TFnumber]=gsl_ran_binomial(_Random,daughterCellProp,_motherCell.m_numberTranscriptionFactors[TFnumber]);
        _motherCell.m_numberTranscriptionFactors[TFnumber]-=m_numberTranscriptionFactors[TFnumber];
    }
    
    /*
    //v.Maintenance metabolites content
    m_numberMaintenanceMetabolites=gsl_ran_binomial(Random,m_volumeProportionDescendant,_cell.m_numberMaintenanceMetabolites);
    _cell.m_numberMaintenanceMetabolites-=m_numberMaintenanceMetabolites;
    m_numberPublicGoodsMetabolites=gsl_ran_binomial(Random,m_volumeProportionDescendant,_cell.m_numberPublicGoodsMetabolites);
    _cell.m_numberPublicGoodsMetabolites-=m_numberPublicGoodsMetabolites;*/
    
    /*3.Genetics*/
    
    //i.Cell volume at mitosis
    m_volumeCellMitosis=(_motherCell.m_volumeCellMitosis); //Questionned: is it necessary?
    
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
    delete[] m_numberNutUnits;
}

/* ##Other methods## */

// #Printing the state of the cell# //

void Cell::PrintCell(){
    
    cout << "Generation number of the cell: " << m_generationNumber << endl;
    cout << "Nutrient 1 proteins: " << m_numberNut1Proteins << endl;
    cout << "Nutrient 2 proteins: " << m_numberNut2Proteins << endl;
    cout << "Volume at mitosis: " << m_volumeCellMitosis << endl;
    cout << "Volume: " << m_volumeCell << endl;
    cout << "Nutrient 1 transcripts:" << m_numberNut1Transcripts << endl;
    cout << "Nutrient 2 transcripts: " << m_numberNut2Transcripts << endl;
    
    cout << "Connectivity matrix:" << endl;
    for (int typeOfTFNumber=0;typeOfTFNumber<m_regulatoryAlleles.size();typeOfTFNumber++){
        for (int gene=0;gene<(m_regulatoryAlleles.size()+m_physiologicalAlleles.size());gene++){
            cout << m_connectivityMatrix[gene][typeOfTFNumber] << " ";
        }
        cout << endl;
    }
}

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

/* #Operations on the cell# */

/* #0.Reevaluate connectivity vector# */

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

void Cell::ReevaluateBSTFConnectivityAfterFactorSequenceMutation(int transcriptionFactorNumber){//Mutation on recognition sequence
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

void Cell::ReevaluateBSTFConnectivityAfterRegulatoryFactorSequenceMutation(int transcriptionFactorNumber){//Mutation on TF regulatory sequence
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

/*
void Cell::ReconstructBSTFConnectivityMatrixAfterEfficiencyPhysiologicalMutation(string function){
    int _PhysioGenenumber=(int) m_physiologicalAlleles.size();
    int _TFGeneNumber=(int) m_regulatoryAlleles.size();
    int totalGeneNumber= _TFGeneNumber + _PhysioGenenumber;
    for (int Gene=0;Gene<totalGeneNumber;Gene++){
        delete[] m_connectivityMatrix[Gene];
    }
    delete[] m_connectivityMatrix;
    
    auto find1 = find(begin(m_physiologicalFunctions), end(m_physiologicalFunctions), function);
    int functionPosition=(int) distance(begin(m_physiologicalFunctions),find1);
    m_physiologicalAlleles.erase(function);
    m_physiologicalFunctions.erase(m_physiologicalFunctions.begin()+functionPosition);
    
    _PhysioGenenumber=(int) m_physiologicalAlleles.size();
    _TFGeneNumber=(int) m_regulatoryAlleles.size();
    totalGeneNumber= _TFGeneNumber + _PhysioGenenumber;
    m_connectivityMatrix=new int*[totalGeneNumber];
    for (int Gene=0;Gene<totalGeneNumber;Gene++){
        m_connectivityMatrix[Gene]=new int[_TFGeneNumber];
    }
    
    for (int typeOfTFNumber=0;typeOfTFNumber<_TFGeneNumber;typeOfTFNumber++){
        for (int physiologicalGeneNumber=0;physiologicalGeneNumber<m_physiologicalAlleles.size();physiologicalGeneNumber++){
            m_connectivityMatrix[physiologicalGeneNumber][typeOfTFNumber]= m_physiologicalAlleles[m_physiologicalFunctions[physiologicalGeneNumber]].BSTFMatchCount(m_regulatoryAlleles[typeOfTFNumber]);
        }
        for (int typeOfTFNumber_tested=0;typeOfTFNumber_tested<_TFGeneNumber;typeOfTFNumber_tested++){
            m_connectivityMatrix[_PhysioGenenumber+typeOfTFNumber_tested][typeOfTFNumber]=m_regulatoryAlleles[typeOfTFNumber_tested].BSTFMatchCount(m_regulatoryAlleles[typeOfTFNumber]);
        }
    }
    for (int typeofTFNumber=0;typeofTFNumber<_TFGeneNumber;typeofTFNumber++){
        m_totalBindingSites[typeofTFNumber]=0; //added from the Master model
        for (int gene_number=0;gene_number<totalGeneNumber;gene_number++){
            m_totalBindingSites[typeofTFNumber]+=m_connectivityMatrix[gene_number][typeofTFNumber];
        }
    }
}*/

/* #1.Transcripts production# */
void Cell::TranscriptsDynamics(double delta_t){
    int _TFGeneNumber=(int) m_regulatoryAlleles.size();
    int _PhysioGenenumber=(int) m_physiologicalAlleles.size();
    
    //o.Calculus of the BS and TF total concentrations in the cell, i.e. including those bound and those free
    double cBStot[_TFGeneNumber];
    double cTFtot[_TFGeneNumber];
    for (int TF_number=0;TF_number<_TFGeneNumber;TF_number++){
        if(m_totalBindingSites[TF_number]>0){
            cBStot[TF_number]=Molecule_micro_toMole_L_Concentration(ConcentrationCalculus(m_totalBindingSites[TF_number], m_volumeCell)) ;
        }
        else{
            cBStot[TF_number]=0;
        }
        if(m_numberTranscriptionFactors[TF_number]>0){
            cTFtot[TF_number]=Molecule_micro_toMole_L_Concentration(ConcentrationCalculus(m_numberTranscriptionFactors[TF_number], m_volumeCell));
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
            (m_K__off-pow(pow(m_K__off+m_K__on*(cBStot[TF_number]+cTFtot[TF_number]),2.)-4*pow(m_K__on,2)*cBStot[TF_number]*cTFtot[TF_number],0.5))/(2*m_K__on);
        }
        else{
            cBSTFequilibrium[TF_number]=0;
        }
    }
    double* probability_nonBSTF=new double[_TFGeneNumber];
    for (int TF_number=0;TF_number<_TFGeneNumber;TF_number++){
        probability_nonBSTF[TF_number]=1-cBSTFequilibrium[TF_number]/cBStot[TF_number];//(1-P_BSTF) in theoretical writings
    }
    
    //iii.Calculus of the transcription rate of each gene
    //#Initialization#//
    map <string,double> tauTranscriptionPhysiologicalGenes;
    
    for (int geneNumber=0;geneNumber<_PhysioGenenumber;geneNumber++){
        tauTranscriptionPhysiologicalGenes.insert(pair<string,double>(m_physiologicalFunctions[geneNumber],m_tauTranscriptionOncePolymeraseBound));
    }
    double tauTranscriptionFactorsGeneTranscription[_TFGeneNumber];
    //#Calculus#//
    
    double bindingPolymeraseProbabilityWhileNoHelp[_PhysioGenenumber];
    double positiveTranscriptionFactorNonBindingProbability[_PhysioGenenumber];
    double bindingPolymeraseProbabilityWhileHelpAndNoRepression[_PhysioGenenumber];
    
    //a.Physiological genes transcription //Need to be completed if genes are to be added for each function should be individualized through its name
    for (int geneNumber=0;geneNumber<_PhysioGenenumber;geneNumber++){
        bindingPolymeraseProbabilityWhileNoHelp[geneNumber]=m_basalBindingProbabilityPolymerase*WeightProduct(_TFGeneNumber,probability_nonBSTF,m_connectivityMatrix[geneNumber]);
        
        positiveTranscriptionFactorNonBindingProbability[geneNumber]=1;
        bindingPolymeraseProbabilityWhileHelpAndNoRepression[geneNumber]=m_enhancedBindingProbabilityPolymerase;
        for (int TF_number=0;TF_number<_TFGeneNumber;TF_number++){
            if(m_regulatoryAlleles[TF_number].NatureOfEffect()=='A'){
                positiveTranscriptionFactorNonBindingProbability[geneNumber]*=pow((probability_nonBSTF[TF_number]),m_connectivityMatrix[geneNumber][TF_number]); //Probability of at least one activator bound to one of each BS
            }
            else if(m_regulatoryAlleles[TF_number].NatureOfEffect()=='R'){
                bindingPolymeraseProbabilityWhileHelpAndNoRepression[geneNumber]*=(pow(probability_nonBSTF[TF_number],m_connectivityMatrix[geneNumber][TF_number]));
            }
        }
        bindingPolymeraseProbabilityWhileHelpAndNoRepression[geneNumber]*=1-positiveTranscriptionFactorNonBindingProbability[geneNumber];
        
        tauTranscriptionPhysiologicalGenes[m_physiologicalFunctions[geneNumber]]*=
            (bindingPolymeraseProbabilityWhileNoHelp[geneNumber]+bindingPolymeraseProbabilityWhileHelpAndNoRepression[geneNumber]);
    }
    
    /*Checking the probabilities of being bound*/
    //cout << "bindingprob sans aide: "<<positiveTranscriptionFactorBindingProbability[1] << endl;
    //cout << "prob bind sans aide: " << bindingPolymeraseProbabilityWhileNoHelp[1] << endl;
    //cout << "prob bind avec aide: " << bindingPolymeraseProbabilityWhileHelpAndNoRepression[1] << endl;
    
    //e.Transcription factors transcription
    
    double TFbindingPolymeraseProbabilityWhileNoHelp[_TFGeneNumber];
    double TFpositiveTranscriptionFactorBindingProbability[_TFGeneNumber];
    double TFbindingPolymeraseProbabilityWhileHelpAndNoRepression[_TFGeneNumber];
    
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
    
    //a.Silencing the switch between replication or no replication
    /*if (ReplicationCheckpoint()){
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
    else{*/
        if(m_numberNut1Transcripts>0){
            m_numberNut1Transcripts-=gsl_ran_binomial(_Random,m_tauDegradationTranscript*delta_t,m_numberNut1Transcripts);
        }
        
        m_numberNut1Transcripts+=gsl_ran_poisson(_Random,tauTranscriptionPhysiologicalGenes["Nut1"]*delta_t);
        
        if(m_numberNut2Transcripts>0){
            m_numberNut2Transcripts-=gsl_ran_binomial(_Random,m_tauDegradationTranscript*delta_t,m_numberNut2Transcripts);
        }
        
        m_numberNut2Transcripts+=gsl_ran_poisson(_Random,tauTranscriptionPhysiologicalGenes["Nut2"]*delta_t);
        
        for (int TF_number=0;TF_number<_TFGeneNumber;TF_number++){
            if(m_numberTFTranscripts[TF_number]>0){
                m_numberTFTranscripts[TF_number]-=gsl_ran_binomial(_Random,m_tauDegradationTranscript*delta_t,m_numberTFTranscripts[TF_number]);
            }
            
            m_numberTFTranscripts[TF_number]+=gsl_ran_poisson(_Random,tauTranscriptionFactorsGeneTranscription[TF_number]*delta_t);
        }
    //}
    delete[] probability_nonBSTF;
}// In principle, need for a general function with balance of transcript


/* #2.Proteins Production# */

void Cell::ProteinsDynamics(double delta_t){
    
    //i.Nutrient proteins dynamics
    //Nutrient 1
    if(m_numberNut1Proteins>0){
        m_numberNut1Proteins-=gsl_ran_binomial(_Random,m_tauDegradationProtein*delta_t,m_numberNut1Proteins);
    }
    
    if(m_numberNut1Transcripts>0){
        m_numberNut1Proteins+=gsl_ran_poisson(_Random,m_tauTraduction*delta_t*m_numberNut1Transcripts);
    }
    //Nutrient 2
    if(m_numberNut2Proteins>0){
        m_numberNut2Proteins-=gsl_ran_binomial(_Random,m_tauDegradationProtein*delta_t,m_numberNut2Proteins);
    }
    
    if(m_numberNut2Transcripts>0){
        m_numberNut2Proteins+=gsl_ran_poisson(_Random,m_tauTraduction*delta_t*m_numberNut2Transcripts);
    }
    
    
    //While synthesis, translation still goes on, but while G2 ansd mitosis, translation decreases a lot and eventually stops: considering that synthesis time encompasses those three steps and that G2/mitosis lasts as long as synthesis (according to litterature, reasonable statement), explaining the ratio of 2. Notice that here, we have assumed that growh only occurs in a G1 manner, assuming the one growth step model does not distort evolution outcomes.
    
    //ii.Transcription factors dynamics
    int amountTFgenes= (int) m_regulatoryAlleles.size();
    
    for (int TFnumber=0;TFnumber<amountTFgenes;TFnumber++){
        if(m_numberTranscriptionFactors[TFnumber]>0){
            m_numberTranscriptionFactors[TFnumber]-=gsl_ran_binomial(_Random,m_tauDegradationProtein*delta_t,m_numberTranscriptionFactors[TFnumber]);
        }
        
        if(m_numberTFTranscripts[TFnumber]>0){
            m_numberTranscriptionFactors[TFnumber]+=gsl_ran_poisson(_Random,m_tauTraduction*delta_t*m_numberTFTranscripts[TFnumber]);
        }
    }
}

/* #3.Nutrient uptake# */
double Cell::NutrientUptake(double concNutEnv, double concNutInt, double dt){
    double dNutInt_dt=(-m_permeabilityCell*VtoSAconversion(m_volumeCell)*(concNutInt-concNutEnv)*dt);
    return(dNutInt_dt);
} // change here in that nutrients are counted in energy units.

void Cell::SetNutrientsAfterDiffusion(double V_comp,double* &concNutEnv,double* concNutEnvInit, double dt){
    double NutUptake[m_numberNutrientTypes];
    //std::cout << "m_NutCellpre: " <<  m_numberNutUnits[0] << std::endl;
    //std::cout << "concEnvInit:" <<concNutEnvInit[0] << std::endl;
    for (int nutTypeAmount=0;nutTypeAmount<m_numberNutrientTypes;nutTypeAmount++){
        NutUptake[nutTypeAmount]=NutrientUptake(concNutEnvInit[nutTypeAmount], Molecule_micro_toMole_L_Concentration(ConcentrationCalculus(m_numberNutUnits[nutTypeAmount], m_volumeCell)) , dt)*UptakeConversion;//Need to initalize a pointer
        if (NutUptake[nutTypeAmount]>0){
            m_numberNutUnits[nutTypeAmount]+=NutUptake[nutTypeAmount];
            //std::cout << "m_NutCell: " <<  m_numberNutUnits[nutTypeAmount] << std::endl;
            concNutEnv[nutTypeAmount]-=Molecule_micro_toMole_L_Concentration(NutUptake[nutTypeAmount]/V_comp);
            //std::cout << "concEnv: " <<  concNutEnv[nutTypeAmount] << std::endl;
            //std::cout << "concNutEnv" << concNutEnv[nutTypeAmount] << std::endl;
        }
        else{}
        
    }
    //std::cout << "Nutrient upatke: " << m_numberNutUnits[0] << std::endl;
}



//bool EnergyBudgetUsefulnessTest(double delta_t); /not needed for the moment
/* #4.Energy production# */
void Cell::EnergyProduction(double delta_t){
    double Copt=5;
    //std::cout << "Nut1prot:" << m_numberNutUnits[0] << std::endl;
    //std::cout << "Vol:" << m_volumeCell << std::endl;
    double Prot=4+(m_numberNut1Proteins+m_numberNut2Proteins)/(m_volumeCell*1000); //Necessary thinking about the relevance to define it through parameters of transcription and translation (same for variation in the rates according to allocation to basal metabolism)
    //std::cout << "Prot" << Prot << std::endl;
    double K__M=m_K__M*pow(10,+abs(Prot-Copt)/Copt);//recalculus of KM
    //std::cout << "KM" << K__M << std::endl;
    //std::cout << "Conc Nut 1:" << Molecule_micro_toMole_L_Concentration(m_numberNutUnits[0]/m_volumeCell) << std::endl;
    double energyAmountProduced[2]={
        m_k__cat*m_numberNut1Proteins*Molecule_micro_toMole_L_Concentration(m_numberNutUnits[0]/m_volumeCell)/(K__M+Molecule_micro_toMole_L_Concentration(m_numberNutUnits[0]/m_volumeCell))
                                                                                ,
        m_k__cat*m_numberNut2Proteins*Molecule_micro_toMole_L_Concentration(m_numberNutUnits[1]/m_volumeCell)/(K__M+Molecule_micro_toMole_L_Concentration(m_numberNutUnits[1]/m_volumeCell))};
    //std::cout << "Energy: " << Molecule_micro_toMole_L_Concentration(m_numberNutUnits[0]/m_volumeCell) << std::endl;
    for (int _nut=0;_nut<2;_nut++){
        if(energyAmountProduced[_nut]>m_numberNutUnits[_nut]){
            m_numberEnergyUnits+=m_numberNutUnits[_nut];
            m_numberNutUnits[_nut]=0;
        }
        else{
            m_numberEnergyUnits+=energyAmountProduced[_nut];
            m_numberNutUnits[_nut]-=energyAmountProduced[_nut];
        }
    }
    //std::cout << "Energy produced: " << m_numberEnergyUnits << endl;
    //std::cout << "Energy: " << m_numberEnergyUnits << std::endl;
}
/* #5.Energy use# */
void Cell::BasalMetabolismCalculus(double delta_t){
    m_currentBasalMetabolism=5*pow(10,3)*pow(m_volumeCell,0.88)*delta_t; //Lynch and Marinov 2015, with ATPtoNut=26
} // contingent to the choice of modelization procedures
// #i.Withdrawal by the basal metabolism
void Cell::BasalMetabolismWithdrawal(double delta_t){
    //std::cout << "Energy pré basal: " << m_numberEnergyUnits << std::endl;
    BasalMetabolismCalculus(delta_t);
    double propTotalProtNut=(m_numberNut1Proteins+m_numberNut2Proteins)/pow(10,6);
    //std::cout << propTotalProtNut << endl;
    //5 denoting the amount of ribosomes concomitantly translating
    if(m_numberEnergyUnits>m_currentBasalMetabolism){
        m_synthesisRateProp+=1;
        m_numberEnergyUnits-=m_currentBasalMetabolism;
        //std::cout << "Energy before: " << m_numberEnergyUnits << endl;
        m_numberEnergyUnits*=(1-propTotalProtNut);
        //std::cout << "Energy after: " << m_numberEnergyUnits << endl;
    }
    else{
        m_synthesisRateProp+=(m_numberEnergyUnits/m_currentBasalMetabolism);
        m_numberEnergyUnits=0;
    }
    //std::cout << "synthèse basale:" << m_synthesisRateProp << endl;
         //std::cout << "Energy post basal: " << m_numberEnergyUnits << std::endl;
}

void Cell::BasalRateDynamics(int timeStepRatio, double delta_t){
    
    m_synthesisRateProp/=timeStepRatio;
    //std::cout << "synthèse basale:" << m_synthesisRateProp << endl;
    //double degPropBasalMet=exp(-m_tauDegradationProtein*delta_t);
    
    //if (m_mRNA_ribosomes>10){
    m_mRNA_ribosomes-=gsl_ran_binomial(_Random, m_tauDegradationTranscript*delta_t,m_mRNA_ribosomes);
    //}
    //if (m_ribosomesAmount>5*pow(10,3)){
    m_ribosomesAmount-=gsl_ran_binomial(_Random, m_tauDegradationProtein*delta_t,m_ribosomesAmount);
    //}
    
    //cout << m_ribosomesAmount << endl;
    
    double tauMaxTranscription=5*pow(10,-2);
    double tauMaxTranslation=5*pow(10,-1);
    
    double tauMaxDegradation_mRNA=5*pow(10,-4);
    double tauMaxDegradation_proteins=2.5*pow(10,-4);
    
    m_mRNA_ribosomes+=gsl_ran_poisson(_Random,m_synthesisRateProp*tauMaxTranscription*delta_t);
    //std::cout << "ribotr:" << m_mRNA_ribosomes << std::endl;
    m_ribosomesAmount+=gsl_ran_poisson(_Random,m_synthesisRateProp*tauMaxTranslation*delta_t*m_mRNA_ribosomes);
    //std::cout << "ribo:" << m_ribosomesAmount << std::endl;
    
    m_tauDegradationTranscript=tauMaxDegradation_mRNA*m_ribosomesAmount/(2*pow(10,5));
    m_tauDegradationProtein=tauMaxDegradation_proteins*m_ribosomesAmount/(2*pow(10,5));
    
    m_tauTraduction=tauMaxTranslation*m_ribosomesAmount/(2*pow(10,5));
    m_tauTranscriptionOncePolymeraseBound=tauMaxTranscription*m_ribosomesAmount/(2*pow(10,5)); //polymeras dynamics considered as identical as ribosomes one
    /*
    std::cout << "protein degradation rate: " << m_tauDegradationProtein << std::endl;
    std::cout << "transcript degradation rate: " << m_tauDegradationTranscript << std::endl;
    std::cout << "translation rate: " << m_tauTraduction << std::endl;
    std::cout << "transcription rate: " << m_tauTranscriptionOncePolymeraseBound << std::endl;
    */
     
        //translation rates rely on the amount of ribosomes (that we may quantify to 10^5, may be linked to cell sizes through variations in degradations rates of mRNAs and proteins)
    //cout << m_tauTraduction << endl;
    
    //cout << m_numberNut1Proteins << endl;
   
    m_synthesisRateProp=0;
}

// #ii.a)Cell growth
void Cell::Growth(){
    if(m_numberEnergyUnits>0){
        double dS_dt=m_numberEnergyUnits/(10*pow(10,8));
        m_volumeCell=VtoSAconversion(SAtoVconversion(m_volumeCell)+dS_dt);
        m_numberEnergyUnits=0;
    }
    //std::cout << "Cell size: " << m_volumeCell << endl;
} //via phospholipids production
// #ii.b)Cell maintenance molecule production

/* #7.Mitosis checking */

// #i.Checking if mitosis checkpoint reached#
bool Cell::ReplicationCheckpoint(){
    return(m_volumeCell>m_volumeCellMitosis);
}

/* #8.Death checking */
bool Cell::DeathTest(double delta_t){
    double basalDeathRate=0.1;
    double survivalProbability;
    survivalProbability=pow((1-basalDeathRate),(delta_t/5000.0));
    int deathTest=gsl_ran_binomial(_Random, survivalProbability, 1);
    return(deathTest==0);
}

/* #9.Mutation processing */
void Cell::MutationProcessing(){
    int _TFGeneNumber=(int) m_regulatoryAlleles.size();
    int _PhysioGenenumber=(int) m_physiologicalAlleles.size();
    // #Regulatory sequences and binding-site recognition sequences# //
    if (_PhysioGenenumber>0){
        for (int _physiologicalGeneNumber=0;_physiologicalGeneNumber<_PhysioGenenumber;_physiologicalGeneNumber++){
        if(m_physiologicalAlleles[m_physiologicalFunctions[_physiologicalGeneNumber]].RegulatorySequenceMutation()){
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
    int volumeMitosisMutation;
    volumeMitosisMutation=gsl_ran_binomial(_Random,m_volumeMutationRate, 1);
    if (volumeMitosisMutation==1){
        m_volumeCellMitosis+=gsl_ran_gaussian(_Random,2);
    }
}

/* #Getting the settings# */

double Cell::GetConcentrationNut1(){
    return(Molecule_micro_toMole_L_Concentration(m_numberNutUnits[0]/m_volumeCell));
}

double Cell::GetPermeabilityCell(){
    return(m_permeabilityCell);
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

double Cell::GetVolumeCell(){
    return(m_volumeCell);
}

double Cell::GetSACell(){
    return(VtoSAconversion(m_volumeCell));
}


int Cell::GetTranscripts(std::string t){
    if (t=="Nut1"){
        return(m_numberNut1Transcripts);
    }
    else if (t=="Nut2"){
        return(m_numberNut2Transcripts);
    }
    else{
        return(0);
    }
}

int Cell::GetTFTranscripts(int g){
    return(m_numberTFTranscripts[g]);
}

int Cell::GetProteins(std::string t){
    if (t=="Nut1"){
        return(m_numberNut1Proteins);
    }
    else if (t=="Nut2"){
        return(m_numberNut2Proteins);
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

//Connectivity matrix conversion
void Cell::MatrixConversion(std::vector <int> &matrixCode){
    for (int gene=0; gene< m_regulatoryAlleles.size()+m_physiologicalFunctions.size();gene++){
        for (int TF=0; TF<m_regulatoryAlleles.size();TF++){
            matrixCode.push_back(m_connectivityMatrix[gene][TF]);
        }
    }
}
/*
void Cell::GeneTabularConversion(string &geneCode){
    int count=0;
    while (count<m_physiologicalFunctions.size()){
        geneCode+=m_physiologicalFunctions[count];
        count+=1;
    }
}*/

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

double Molecule_micro_toMole_L_Concentration(double MoleculeConc){
    return(MoleculeConc/(6.02*pow(10,8)));
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
