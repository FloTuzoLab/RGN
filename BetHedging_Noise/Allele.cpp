//
//  Allele.cpp
//  BetHedging
//
//  Created by florian labourel on 30/06/2019.
//  Copyright © 2019 florian labourel. All rights reserved.
//

#include "Allele.hpp"
using namespace std;
Allele::Allele():
m_mutationRate(pow(10,-6))
{
    gsl_rng *rnd=gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rnd, 1);
    _random_Simulation=rnd;
    for (int position=0;position<=999;position++){
        m_regulatorySequence.push_back('A');
    }
    gsl_rng_free(rnd);
}

Allele::Allele(std::string regulatorySequence,double mutationRate,gsl_rng *random_Simulation):
m_mutationRate(mutationRate)
{
    _random_Simulation=random_Simulation;
    for (int position=0;position<=regulatorySequence.length();position++){
        m_regulatorySequence.push_back(regulatorySequence[position]);
    }
}
Allele::Allele(const Allele& m_allele):
m_mutationRate(m_allele.m_mutationRate),
_random_Simulation(m_allele._random_Simulation)
{
    for (int position=0;position<=m_allele.m_regulatorySequence.length();position++){
        m_regulatorySequence[position]=m_allele.m_regulatorySequence[position];
    }
}

// #Destructor# //
Allele::~Allele()
{
    //Useless here
}

// Operator = //
void Allele::operator=(const Allele& _allele)
{
    m_mutationRate=_allele.m_mutationRate;
    _random_Simulation=_allele._random_Simulation;
    for (int position=0;position<=_allele.m_regulatorySequence.length();position++){
        m_regulatorySequence[position]=_allele.m_regulatorySequence[position];
    }
}

// #Other methods# //

// Printing function //
void Allele::PrintAllele()
{
    cout << "Mutation rate:" << m_mutationRate << endl;
    cout << "Regulatory sequence:" << endl << m_regulatorySequence << endl;
}

bool Allele::RegulatorySequenceMutation()
{
    vector <int> _mutationPosition;
    //char mutationValue[1];
    int mutation=0;
    for (int sitePosition=0;sitePosition<m_regulatorySequence.length();sitePosition++){
        mutation=gsl_ran_binomial(_random_Simulation, m_mutationRate, 1);
        if (mutation==1){
            _mutationPosition.push_back(sitePosition);
        }
    }
    for (int mutationPosition=0;mutationPosition<_mutationPosition.size();mutationPosition++){
        
        switch(m_regulatorySequence[_mutationPosition[mutationPosition]]){
            case 'A':{
                char nucleotideAlphabet[3]={'T','G','C'}; gsl_ran_choose(_random_Simulation,&m_regulatorySequence[_mutationPosition[mutationPosition]],1,nucleotideAlphabet,3,sizeof(char));
            }
                break;
            case 'C':{
                char nucleotideAlphabet[3]={'T','G','A'}; gsl_ran_choose(_random_Simulation,&m_regulatorySequence[_mutationPosition[mutationPosition]],1,nucleotideAlphabet,3,sizeof(char));
            }
                break;
            case 'G':{
                char nucleotideAlphabet[3]={'T','C','A'}; gsl_ran_choose(_random_Simulation,&m_regulatorySequence[_mutationPosition[mutationPosition]],1,nucleotideAlphabet,3,sizeof(char));
            }
                break;
            case 'T':{
               char nucleotideAlphabet[3]={'G','C','A'}; gsl_ran_choose(_random_Simulation,&m_regulatorySequence[_mutationPosition[mutationPosition]],1,nucleotideAlphabet,3,sizeof(char));
            }
                break;
        }
        
    }
    return (_mutationPosition.size()>0);
}


//Class PhysiologicalAllele//
/* ##Specific features of a physiological allele object## */

/* Constructor */

PhysiologicalAllele::PhysiologicalAllele()
{
}

PhysiologicalAllele::PhysiologicalAllele(std::string regulatorySequence,
                                         double mutationRate,
                                         gsl_rng *random_Simulation,
                                         std::string nameFunction
                                         ):
Allele(regulatorySequence,mutationRate,random_Simulation)
{
    m_functionGene=nameFunction;
}


PhysiologicalAllele::PhysiologicalAllele(const PhysiologicalAllele& m_physiologicalAllele):
Allele(m_physiologicalAllele.m_regulatorySequence,m_physiologicalAllele.m_mutationRate,m_physiologicalAllele._random_Simulation)
{
    m_functionGene=m_physiologicalAllele.m_functionGene;
}

void PhysiologicalAllele::operator=(const PhysiologicalAllele& m_physiologicalAllele)
{
    _random_Simulation=m_physiologicalAllele._random_Simulation;
    m_mutationRate=m_physiologicalAllele.m_mutationRate;
    m_regulatorySequence=m_physiologicalAllele.m_regulatorySequence;
    m_functionGene=m_physiologicalAllele.m_functionGene;
}

//Class RegulatoryAllele//
/* ##Specific features of a regulatory allele object## */

/* Constructor */

//Removal of default constructor (if problem)

RegulatoryAllele::RegulatoryAllele(std::string regulatorySequence,
                                   double mutationRate,
                                   gsl_rng *random_Simulation,
                                   std::string bindingSequence,
                                   char regulationEffect):
Allele(regulatorySequence,mutationRate,random_Simulation),
m_regulationEffect(regulationEffect)
{
    m_bindingSequence=bindingSequence;
}

RegulatoryAllele::RegulatoryAllele(const RegulatoryAllele& m_regulatoryAllele):
Allele(m_regulatoryAllele.m_regulatorySequence,m_regulatoryAllele.m_mutationRate,m_regulatoryAllele._random_Simulation)
{
    m_regulationEffect=m_regulatoryAllele.m_regulationEffect;
    m_bindingSequence=m_regulatoryAllele.m_bindingSequence;
}

void RegulatoryAllele::operator=(const RegulatoryAllele& m_regulatoryAllele)
{
    _random_Simulation=m_regulatoryAllele._random_Simulation;
    m_regulationEffect=m_regulatoryAllele.m_regulationEffect;
    m_mutationRate=m_regulatoryAllele.m_mutationRate;
    m_regulatorySequence=m_regulatoryAllele.m_regulatorySequence;
    for (int position=0;position<8;position++){
        m_bindingSequence[position]=m_regulatoryAllele.m_bindingSequence[position];
    }
}

//BSTFMatchCount needs to be better re-understood

int RegulatoryAllele::BSTFMatchCount(const RegulatoryAllele& m_regulatoryAllele)
{
    int count = 0;
    size_t positionNumber = m_regulatorySequence.find(m_regulatoryAllele.m_bindingSequence, 0); // first occurrence
    while(positionNumber != string::npos)
    {
        count++;
        positionNumber = m_regulatorySequence.find(m_regulatoryAllele.m_bindingSequence, positionNumber + 1);
    }
    return (count);
}

//regulatoryallele passed for reference from outside the object being physiological
int PhysiologicalAllele::BSTFMatchCount(const RegulatoryAllele& _regulatoryAllele){
    int count = 0;
    size_t positionNumber = m_regulatorySequence.find(_regulatoryAllele.m_bindingSequence, 0); // first occurrence
    while(positionNumber != string::npos)
    {
        count++;
        positionNumber = m_regulatorySequence.find(_regulatoryAllele.m_bindingSequence, positionNumber + 1);
    }
    return (count);
}

bool RegulatoryAllele::BindingSiteRecognitionMutation()
{
    vector <int> _mutationPosition;
    for (int sitePosition=0;sitePosition<m_bindingSequence.length();sitePosition++){
        if (gsl_ran_binomial(_random_Simulation,m_mutationRate,1)==1){
            _mutationPosition.push_back(sitePosition);
        }
    }
    for (int mutationPosition=0;mutationPosition<_mutationPosition.size();mutationPosition++){
        
        switch(m_bindingSequence[_mutationPosition[mutationPosition]]){
            case 'A':{
                char nucleotideAlphabet[3]={'T','G','C'}; gsl_ran_choose(_random_Simulation,&m_bindingSequence[_mutationPosition[mutationPosition]],1,nucleotideAlphabet,3,sizeof(char));
            }
                break;
            case 'C':
            {
                char nucleotideAlphabet[3]={'T','G','A'}; gsl_ran_choose(_random_Simulation,&m_bindingSequence[_mutationPosition[mutationPosition]],1,nucleotideAlphabet,3,sizeof(char));
            }
                break;
            case 'G':
            {
                char nucleotideAlphabet[3]={'T','C','A'}; gsl_ran_choose(_random_Simulation,&m_bindingSequence[_mutationPosition[mutationPosition]],1,nucleotideAlphabet,3,sizeof(char));
            }
                break;
            case 'T':
            {
                char nucleotideAlphabet[3]={'G','C','A'}; gsl_ran_choose(_random_Simulation,&m_bindingSequence[_mutationPosition[mutationPosition]],1,nucleotideAlphabet,3,sizeof(char));
                break;
            }
        }
    }
    return (_mutationPosition.size()>0);
    
}

char RegulatoryAllele::NatureOfEffect()
{
    return (m_regulationEffect);
}

void RegulatoryAllele::PrintRegulatoryAllele()
{
    PrintAllele();
    cout << "Binding-sequence recognised:" << endl << m_bindingSequence << endl;
    cout << "Nature de l'effet: " << m_regulationEffect << endl;
}
//Weird: what happens here?
