//
//  Medium.cpp
//  BetHedging
//
//  Created by florian labourel on 30/06/2019.
//  Copyright © 2019 florian labourel. All rights reserved.
//

#include "Medium.hpp"

Medium::Medium(){
    
}

Medium::Medium(double diffusionNutrients,double diffusionCells):
m_diffusionNutrients(diffusionNutrients),
m_diffusionCells(diffusionCells)
{
}


// #Constructor using a pre-defined medium#/
Medium::Medium(const Medium& m_medium):
m_diffusionNutrients(m_medium.m_diffusionNutrients),
m_diffusionCells(m_medium.m_diffusionCells)
{
}

//Lack an operator = which could allow for simpler parametrization
void Medium::operator=(const Medium& _medium){
    m_diffusionNutrients=(_medium.m_diffusionNutrients);
    m_diffusionCells=(_medium.m_diffusionCells);
}
//Question: do accessors cost time when not used?

/* ##Other methods## */


//  #Setting the settings# //
void Medium::SetDiffusionNutrients(double a){
    m_diffusionNutrients=a;
}

void Medium::SetDiffusionCells(double a){
    m_diffusionCells=a;
}


// #Getting the settings# //
double Medium::GetDiffusionNutrients() const{
    return(m_diffusionNutrients);
}

double Medium::GetDiffusionCells() const{
    return(m_diffusionCommonGoods);
}


// #Printing of the medium features# //
void Medium::Print() const{
    cout<<"Diffusion coefficient of molecules:"<<m_diffusionNutrients<<endl;
    cout<<"Viscosity of the compartment for cell diffusion:"<<m_diffusionCells<<endl;
}
