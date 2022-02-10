//
//  Medium.cpp
//  Cell_Evolution
//
//  Created by Florian on 21/02/2018.
//  Copyright © 2018 Florian. All rights reserved.
//

#include "Medium.hpp"

using namespace std;

/* ###Methods### */


/* ##Constructors## */


// #Default constructor# //
Medium::Medium():
    m_diffusionNutrients(pow(10,2)),//Need to find the unit, in order to homgenize
    m_diffusionCommonGoods(pow(10,1)),//Same as below
    m_diffusionCells(0) //Number of cells that leave the compartment: need for a conversion towards the same unit as below molecules
{
}


// #Total constructor# //
Medium::Medium(double diffusionNutrients,double diffusionCommonGoods,double diffusionCells):
    m_diffusionNutrients(diffusionNutrients),
    m_diffusionCommonGoods(diffusionCommonGoods),
    m_diffusionCells(diffusionCells)
{
}


// #Constructor using a pre-defined medium#/
Medium::Medium(const Medium& m_medium):
    m_diffusionNutrients(m_medium.m_diffusionNutrients),
    m_diffusionCommonGoods(m_medium.m_diffusionCommonGoods),
    m_diffusionCells(m_medium.m_diffusionCells)
{
}

//Lack an operator = which could allow for simpler parametrization
//Question: do accessors cost time when not used?

/* ##Other methods## */


//  #Setting the settings# //
void Medium::SetDiffusionNutrients(double a){
    m_diffusionNutrients=a;
}

void Medium::SetDiffusionCommonGoods(double a){
    m_diffusionCommonGoods=a;
}

void Medium::SetDiffusionCells(double a){
    m_diffusionCells=a;
}


// #Getting the settings# //
double Medium::GetDiffusionNutrients() const{
    return(m_diffusionNutrients);
}

double Medium::GetDiffusionCommonGoods() const{
    return(m_diffusionCommonGoods);
}

double Medium::GetDiffusionCells() const{
    return(m_diffusionCommonGoods);
}


// #Printing of the medium features# //
void Medium::Afficher() const{
    cout<<"Diffusion coefficient of molecules:"<<m_diffusionNutrients<<endl;
    cout<<"Diffusion coefficient of common goods:"<<m_diffusionCommonGoods<<endl;
    cout<<"Viscosity of the compartment for cell diffusion:"<<m_diffusionCells<<endl;
}
