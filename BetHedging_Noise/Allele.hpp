//
//  Allele.hpp
//  BetHedging
//
//  Created by florian labourel on 30/06/2019.
//  Copyright © 2019 florian labourel. All rights reserved.
//

#ifndef Allele_hpp
#define Allele_hpp

#include <iostream>
#include <cstdio>
#include <string>

#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <random>

using namespace std;

/*## 1.Generic allele definition ##*/
class Allele
{
    
public:
    /* ###Methods### */
    
    /* ##Constructors## */
    Allele();
    
    // #Total free constructor# //
    Allele(std::string regulatorySequence,double mutationRate,gsl_rng *random_Simulation);
    
    // #Constructor using a full pre-defined cell#
    Allele(const Allele& m_allele);
    
    /* ##Destructor## */
    
    ~Allele();
    
    // #Operator = # //
    void operator=(const Allele& _allele);
    
    /* ##Other methods## */
    
    // #Printing the state of the allele# //
    void PrintAllele();
    bool RegulatorySequenceMutation();
    
protected:
    std::string m_regulatorySequence;
    double m_mutationRate;
    gsl_rng *_random_Simulation;//(pseudo) random number used in the processing for the diferent draws.
};

/*## 2.Pre-initalization of a regulatory allele##*/


class RegulatoryAllele;

/*## 3.Specification of physiological allele definition ##*/

class PhysiologicalAllele : public Allele
{
public:
    PhysiologicalAllele();
    PhysiologicalAllele(std::string regulatorySequence,double mutationRate,gsl_rng *random_Simulation,std::string nameFunction);
    PhysiologicalAllele(const PhysiologicalAllele& m_physiologicalAllele);
    void operator=(const PhysiologicalAllele& m_physiologicalAllele);
    int BSTFMatchCount(const RegulatoryAllele& _regulatoryAllele);
    //bool EfficiencyMutation();
    
private:
    //bool m_efficiencyGene;
    //bool m_efficiencyMutation;
    std::string m_functionGene;
};

/*## 4.Specification of regulatory allele definition ##*/
class RegulatoryAllele : public Allele
{
public:
    RegulatoryAllele();
    RegulatoryAllele(std::string regulatorySequence,double mutationRate,gsl_rng *random_Simulation, std::string bindingSequence,char regulationEffect);
    RegulatoryAllele(const RegulatoryAllele& m_regulatoryallele);
    void operator=(const RegulatoryAllele& m_regulatoryAllele);
    int BSTFMatchCount(const RegulatoryAllele& m_regulatoryAllele);
    friend int PhysiologicalAllele::BSTFMatchCount(const RegulatoryAllele& _regulatoryAllele);
    bool BindingSiteRecognitionMutation();
    char NatureOfEffect();
    void PrintRegulatoryAllele();
    
private:
    std::string m_bindingSequence;
    char m_regulationEffect;//take {'A' for activator; 'R' for repressor}
};




#endif /* Allele_hpp */
