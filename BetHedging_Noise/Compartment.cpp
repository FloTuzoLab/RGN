//
//  Compartment.cpp
//  BetHedging
//
//  Created by florian labourel on 01/07/2019.
//  Copyright © 2019 florian labourel. All rights reserved.
//

#include "Compartment.hpp"

using namespace std;

/* ###Generic functions### */

double const AvogadroConstant=6.02*pow(10.0,23); //unit: molecules number

/* ## Calculus of the volume of a cube## */

double VolumeCalculus(double length, double height){
    return pow(length,2)*height;
    }

/* ###Methods### */
Compartment::Compartment(){
   
    m_labelCompartment=new int[2];
    m_labelCompartment[0]=0;
    m_labelCompartment[1]=0;
    
    m_labelNeighbours=new int*[4];
    m_labelNeighbours[0]=new int[2];
    m_labelNeighbours[0][0]=0;
    m_labelNeighbours[0][1]=0;
    
    m_labelNeighbours[1]=new int[2];
    m_labelNeighbours[1][0]=0;
    m_labelNeighbours[1][1]=0;
    
    m_labelNeighbours[2]=new int[2];
    m_labelNeighbours[2][0]=0;
    m_labelNeighbours[2][1]=0;
    
    m_labelNeighbours[3]=new int[2];
    m_labelNeighbours[3][0]=0;
    m_labelNeighbours[3][1]=0;
}

Compartment::Compartment(double const length,
                         double height,
                         Medium medium,
                         double concentrationNut1,
                         double concentrationNut2):
/* 2.Compartment shape characterization */
m_length(length), //same as above but defined by the user
m_height(height),
m_ActualVolume(Volume()),
/* 3.Compartment chemical properties */

//i.Diffusion properties
m_medium(medium),

//ii.Compartment concentrations
m_concentrationNut1(concentrationNut1),
m_concentrationNut2(concentrationNut2)
{
    
    /* 1.References of the compartment and its neighbours */
    //i.References of the compartment
    m_labelCompartment=new int[2];
    m_labelCompartment[0]=0;
    m_labelCompartment[1]=0;
    
    m_labelNeighbours=new int*[4];
    m_labelNeighbours[0]=new int[2];
    m_labelNeighbours[0][0]=0;
    m_labelNeighbours[0][1]=0;
    
    m_labelNeighbours[1]=new int[2];
    m_labelNeighbours[1][0]=0;
    m_labelNeighbours[1][1]=0;
    
    m_labelNeighbours[2]=new int[2];
    m_labelNeighbours[2][0]=0;
    m_labelNeighbours[2][1]=0;
    
    m_labelNeighbours[3]=new int[2];
    m_labelNeighbours[3][0]=0;
    m_labelNeighbours[3][1]=0;
}

/* ##Constructors## */
// #III Complete definition constructor: no pre-specification, even for the type of neighbouring relationships# //
Compartment::Compartment(int reference[2],
                           double const length,
                           double height,
                           Medium medium,
                           double concentrationNut1,
                           double concentrationNut2):
/* 2.Compartment shape characterization */
m_length(length), //same as above but defined by the user
m_height(height),
m_ActualVolume(Volume()),

/* 3.Compartment chemical properties */

//i.Diffusion properties
m_medium(medium),

//ii.Compartment concentrations
m_concentrationNut1(concentrationNut1),
m_concentrationNut2(concentrationNut2)
{
    /* 1.References of the compartment and its neighbours */
    //i.References of the compartment
    int r1 = reference[0];
    int r2 = reference[1];
    
    m_labelCompartment=new int[2];
    m_labelCompartment[0]=r1;
    m_labelCompartment[1]=r2;
    
    //ii.References of the neighbours
        m_labelNeighbours=new int*[4];
        m_labelNeighbours[0]=new int[2];
        m_labelNeighbours[0][0]=r1-1;
        m_labelNeighbours[0][1]=r2;
        
        m_labelNeighbours[1]=new int[2];
        m_labelNeighbours[1][0]=r1+1;
        m_labelNeighbours[1][1]=r2;
        
        m_labelNeighbours[2]=new int[2];
        m_labelNeighbours[2][0]=r1;
        m_labelNeighbours[2][1]=r2-1;
        
        m_labelNeighbours[3]=new int[2];
        m_labelNeighbours[3][0]=r1;
        m_labelNeighbours[3][1]=r2+1;
}


// #IV Constructor using a pre-defined compartment# //
Compartment::Compartment(const Compartment& m_compartment):

/* 2.Compartment shape characterization */
m_length(m_compartment.m_length), //same as above but defined by the user
m_height(m_compartment.m_height),
m_ActualVolume(m_compartment.m_ActualVolume),

/* 3.Compartment chemical properties */

//i.Diffusion properties
m_medium(m_compartment.m_medium),

//ii.Compartment concentrations
m_concentrationNut1(m_compartment.m_concentrationNut1),
m_concentrationNut2(m_compartment.m_concentrationNut2)

{
    /* 1.References of the compartment and its neighbours */
    //i.References of the compartment
    
    m_labelCompartment=new int[2];
    m_labelCompartment[0]=m_compartment.m_labelCompartment[0];
    m_labelCompartment[1]=m_compartment.m_labelCompartment[1];
    
   
        m_labelNeighbours=new int*[4];
        m_labelNeighbours[0]=new int[2];
        m_labelNeighbours[1]=new int[2];
        m_labelNeighbours[2]=new int[2];
        m_labelNeighbours[3]=new int[2];
        for (int i=0;i<=3;i++){
            for (int k=0;k<=1;k++){
                m_labelNeighbours[i][k]=m_compartment.m_labelNeighbours[i][k];
            }
        }
}


void Compartment::operator=(const Compartment& _compartment)
{
    if (this != & _compartment){
    delete[] m_labelCompartment;
        for (int i=0;i<=3;i++){
            delete[] m_labelNeighbours[i];
        }
        delete[] m_labelNeighbours;
    }
    m_length=(_compartment.m_length); //same as above but defined by the user
    m_height=(_compartment.m_height);
    m_ActualVolume=(_compartment.m_ActualVolume);
    
    /* 3.Compartment chemical properties */
    
    //i.Diffusion properties
    m_medium=(_compartment.m_medium);
    
    //ii.Compartment concentrations
    m_concentrationNut1=(_compartment.m_concentrationNut1);
    m_concentrationNut2=(_compartment.m_concentrationNut2);
    
    /* 1.References of the compartment and its neighbours */
    //i.References of the compartment
        m_labelCompartment=new int[2];
    m_labelCompartment[0]=_compartment.m_labelCompartment[0];
    m_labelCompartment[1]=_compartment.m_labelCompartment[1];
    
    //ii.References of the neighbours
    
        m_labelNeighbours=new int*[4];
        m_labelNeighbours[0]=new int[2];
        m_labelNeighbours[1]=new int[2];
        m_labelNeighbours[2]=new int[2];
        m_labelNeighbours[3]=new int[2];
        for (int i=0;i<=3;i++){
            for (int k=0;k<=1;k++){
                m_labelNeighbours[i][k]=_compartment.m_labelNeighbours[i][k];
            }
        }
    
}

/* ##Destructor## */
Compartment::~Compartment()
{
    for (int i=0;i<=3;i++){
        delete[] m_labelNeighbours[i];
    }
    delete[] m_labelNeighbours;
    delete[] m_labelCompartment;
}

/* ##Other methods## */

// #Setting the settings# //

void Compartment::SetReference(int r1, int r2){
    m_labelCompartment[0]=r1;
    m_labelCompartment[1]=r2;
}

void Compartment::SetNeighboursRelationships(int numNeighbour,int r1, int r2){
    m_labelNeighbours[numNeighbour][0]=r1;
    m_labelNeighbours[numNeighbour][1]=r2;
}

void Compartment::SetMedium(double diffusionNutrients, double diffusionCells){
    m_medium.SetDiffusionNutrients(diffusionNutrients);
    m_medium.SetDiffusionCells(diffusionCells);
}


void Compartment::SetConcentrationNut1(double a){
    m_concentrationNut1=a;
}

void Compartment::SetConcentrationNut2(double a){
    m_concentrationNut2=a;
}

void Compartment::SetActualVolume(double actualVolume){
    m_ActualVolume=actualVolume;
}


/* #Getting the settings# */
int Compartment::GetNeighbour(int numNeighbour, int numCoordinate) const{
    return (m_labelNeighbours[numNeighbour][numCoordinate]);
}

double Compartment::GetDiffusionNutrients() const{
    return(m_medium.GetDiffusionNutrients());
}

double Compartment::GetDiffusionCells() const{
    return(m_medium.GetDiffusionCells());
}


double Compartment::GetConcentrationNut1() const{
    return (m_concentrationNut1);
}

double Compartment::GetConcentrationNut2() const{
    return(m_concentrationNut2);
}


double Compartment::GetHeight() const{
    return(m_height);
}

double Compartment::GetActualVolume() const{
    return(m_ActualVolume);
}
/*
double* Compartment::ModifyConcentrationNut1() {
    return(&m_concentrationNut1);
}

double* Compartment::ModifyConcentrationNut2() {
    return(&m_concentrationNut2);
}
 */


// #Printing of the compartment state# //
void Compartment::PrintState(){
    cout<<"Reference:(";
    cout<<m_labelCompartment[0]<<","<<m_labelCompartment[1]<<")"<<endl;
    if (m_labelCompartment[0]==0 && m_labelCompartment[1]==0){
        cout<<"Featuring the outside environment, the outside compartment is neighbouring any other compartment. Isn't that obvious, poor fellow?"<<endl;
    }
    
    else if(m_labelCompartment[0]==0 || m_labelCompartment[1]==0){
        for (int l=0; l<3;l++){
            cout<<"Neighbours references:(";
            cout<<m_labelNeighbours[l][0]<<","<<m_labelNeighbours[l][1]<<")"<<endl;
        }
        cout<<"Only given as a piece of information and, for some kind of elegance, I guess. Absolutely useless for the simulations, isn't it?"<<endl;
    }
    
    else{
        for (int l=0; l<4;l++){
            cout<<"Neighbours references:(";
            cout<<m_labelNeighbours[l][0]<<","<<m_labelNeighbours[l][1]<<")"<<endl;
        }
    }
    
    cout<<"Length:"<<m_length<<endl;
    m_medium.Print();
    cout<<"Nutrient 1 concentration:"<<m_concentrationNut1<<endl;
    cout<<"Nutrient 2 concentration:"<<m_concentrationNut2<<endl;
    cout<<"Volume:"<<Volume()<<endl;
}

// #Operations on the compartment# //

//i.Depletion of a compartment by one of the cell located in it
/*
void Compartiment::DepleteCompartment(map <string,double> &totalNutrientBalance, double membranePermeability,double totalMembraneSurface, double delta_t){
    totalNutrientBalance["Standard Nutrient"]=(1-exp(-(membranePermeability*totalMembraneSurface/Volume())*delta_t))*m_concentrationNutrientStandard*Volume()*AvogadroConstant;
    totalNutrientBalance["Revealed Nutrient"]=(1-exp(-(membranePermeability*totalMembraneSurface/Volume())*delta_t))*m_concentrationNutrientRevealed*Volume()*AvogadroConstant;
    m_concentrationNutrientStandard=exp(-(membranePermeability*totalMembraneSurface/Volume())*delta_t)*m_concentrationNutrientStandard;
    m_concentrationNutrientRevealed=exp(-(membranePermeability*totalMembraneSurface/Volume())*delta_t)*m_concentrationNutrientRevealed;
} //une map permet d'associer des éléments n'étant pas de même nature: which unit is used here?
 */ //possibly useless: need to see if really the case.

//Depletion should be done directly via Get and Set functions

/* ###Generic functions### */

/* ## Calculus of the volume of a cube## */
double Compartment::Volume(){
    return VolumeCalculus(m_length,m_height);
}
