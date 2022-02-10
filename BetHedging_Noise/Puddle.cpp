//
//  Puddle.cpp
//  BetHedging_Noise
//
//  Created by Florian on 02/07/2019.
//  Copyright © 2019 Florian. All rights reserved.
//

#include "Puddle.hpp"

using namespace std;

double OperatorFiniteDifferencesExplicit(double alpha_int,double alpha_out, double c, int nbInsideNeighbours, int nbOutsideNeighbours, double n1,double n2,double n3,double n4,double c_extup,double c_extdown){
    return(alpha_int*(n1+n2+n3+n4)+alpha_out*(c_extup+c_extdown)+(1-nbInsideNeighbours*alpha_int-nbOutsideNeighbours*alpha_out)*c);
}

/* ##Methods for the DonutPuddle class## */

DonutPuddle::DonutPuddle(){
    
    Ecosystem();
        m_neighbours=new int***[m_lengthInside+1];
        m_neighbours[0]=new int**[1];
        m_neighbours[0][0]=new int*[1];
        m_neighbours[0][0][0]=0;
        
        for (int i=1;i<=m_lengthInside;i++){
            m_neighbours[i]=new int**[m_widthInside+1];
            for (int j=1;j<=m_widthInside;j++){
                m_neighbours[i][j]=new int*[4];
                for (int k=0;k<=3;k++){
                    m_neighbours[i][j][k]=new int[2];
                }
            }
        }
        
        for (int i=1;i<=m_lengthInside;i++){
            
            for (int j=1;j<=m_widthInside;j++){
                
                //First neighbour
                if (i==1){
                    m_neighbours[i][j][0][0]=m_lengthInside;
                    m_neighbours[i][j][0][1]=j;
                }
                else{
                    m_neighbours[i][j][0][0]=i-1;
                    m_neighbours[i][j][0][1]=j;
                }
                
                //Second neighbour
                if (i==m_lengthInside){
                    m_neighbours[i][j][1][0]=1;
                    m_neighbours[i][j][1][1]=j;
                }
                else{
                    m_neighbours[i][j][1][0]=i+1;
                    m_neighbours[i][j][1][1]=j;
                }
                
                //Third neighbour
                if (j==1){
                    m_neighbours[i][j][2][0]=i;
                    m_neighbours[i][j][2][1]=m_widthInside;
                }
                else{
                    m_neighbours[i][j][2][0]=i;
                    m_neighbours[i][j][2][1]=j-1;
                }
                
                //Fourth neighbour
                if (j==m_widthInside){
                    m_neighbours[i][j][3][0]=i;
                    m_neighbours[i][j][3][1]=1;
                }
                else{
                    m_neighbours[i][j][3][0]=i;
                    m_neighbours[i][j][3][1]=j+1;
                }
                //is it reallyuseful? How much time is lost here?
                m_puddle[i][j].SetNeighboursRelationships(0, m_neighbours[i][j][0][0], m_neighbours[i][j][0][1]);
                m_puddle[i][j].SetNeighboursRelationships(1, m_neighbours[i][j][1][0], m_neighbours[i][j][1][1]);
                m_puddle[i][j].SetNeighboursRelationships(2, m_neighbours[i][j][2][0], m_neighbours[i][j][2][1]);
                m_puddle[i][j].SetNeighboursRelationships(3, m_neighbours[i][j][3][0], m_neighbours[i][j][3][1]);
                
                //Setting the spatio temporal diffusion coefficient for each of the molecules in each compartment.
                
                }
            }
}

// #Constructor# //

DonutPuddle::DonutPuddle(gsl_rng *random_Simulation,
                         double switchProbability,
                         double Nut1Input,
                         double Nut2Input,
                         double delta_t,
                         double delta_l,
                         int timeStepRatio,                  //Ratio between timestep for proteins and transcripts dynamics (quite slow) and other processes especially diffusion (which is highly time sensitive)
                         int lengthPuddle,                      //Length of the puddle (in compartment units)
                         int widthPuddle,                       //Width of the puddle (in copartment units)
                         Compartment insideCompartment,        //Properties of a compartment in the puddle
                         Compartment outsideCompartment,       //Properties of a compartment
                         int initialNumberCell,
                         Cell cell,
                         ofstream* _fileoutcompartments,
                         ofstream* _fileoutcells):

Ecosystem(random_Simulation,
          switchProbability,
          delta_t,
          delta_l,
          timeStepRatio,
          lengthPuddle,
          widthPuddle,
          false,
          Nut1Input,
          Nut2Input,
          insideCompartment,
          outsideCompartment,
          initialNumberCell,
          cell,
          _fileoutcompartments,
          _fileoutcells)
{
    m_neighbours=new int***[m_lengthInside+1];
    m_neighbours[0]=new int**[1];
    m_neighbours[0][0]=new int*[1];
    m_neighbours[0][0][0]=0;
    
    for (int i=1;i<=m_lengthInside;i++){
        m_neighbours[i]=new int**[m_widthInside+1];
        for (int j=1;j<=m_widthInside;j++){
            m_neighbours[i][j]=new int*[4];
            for (int k=0;k<=3;k++){
                m_neighbours[i][j][k]=new int[2];
            }
        }
    }
    
    for (int i=1;i<=m_lengthInside;i++){
        
        for (int j=1;j<=m_widthInside;j++){
            
            //First neighbour
            if (i==1){
                m_neighbours[i][j][0][0]=m_lengthInside;
                m_neighbours[i][j][0][1]=j;
            }
            else{
                m_neighbours[i][j][0][0]=i-1;
                m_neighbours[i][j][0][1]=j;
            }
            
            //Second neighbour
            if (i==m_lengthInside){
                m_neighbours[i][j][1][0]=1;
                m_neighbours[i][j][1][1]=j;
            }
            else{
                m_neighbours[i][j][1][0]=i+1;
                m_neighbours[i][j][1][1]=j;
            }
            
            //Third neighbour
            if (j==1){
                m_neighbours[i][j][2][0]=i;
                m_neighbours[i][j][2][1]=m_widthInside;
            }
            else{
                m_neighbours[i][j][2][0]=i;
                m_neighbours[i][j][2][1]=j-1;
            }
            
            //Fourth neighbour
            if (j==m_widthInside){
                m_neighbours[i][j][3][0]=i;
                m_neighbours[i][j][3][1]=1;
            }
            else{
                m_neighbours[i][j][3][0]=i;
                m_neighbours[i][j][3][1]=j+1;
            }
            //is it reallyuseful? How much time is lost here?
            m_puddle[i][j].SetNeighboursRelationships(0, m_neighbours[i][j][0][0], m_neighbours[i][j][0][1]);
            m_puddle[i][j].SetNeighboursRelationships(1, m_neighbours[i][j][1][0], m_neighbours[i][j][1][1]);
            m_puddle[i][j].SetNeighboursRelationships(2, m_neighbours[i][j][2][0], m_neighbours[i][j][2][1]);
            m_puddle[i][j].SetNeighboursRelationships(3, m_neighbours[i][j][3][0], m_neighbours[i][j][3][1]);
            
            //Setting the spatio temporal diffusion coefficient for each of the molecules in each compartment.
            
        }
    }
}

DonutPuddle::DonutPuddle(const DonutPuddle &_donutPuddle)

{
    m_lengthInside=(_donutPuddle.m_lengthInside);
    m_widthInside=(_donutPuddle.m_lengthInside);
    m_compartmentVolume=(_donutPuddle.m_compartmentVolume);
    m_adjacencyNeed=(_donutPuddle.m_adjacencyNeed);
    m_cellDiffusionRate=(_donutPuddle.m_cellDiffusionRate);
    m_nutrientSwitchProbability=(_donutPuddle.m_nutrientSwitchProbability);
    m_Nutrient1Input=(_donutPuddle.m_Nutrient1Input);
    m_Nutrient2Input=(_donutPuddle.m_Nutrient2Input);
    m_permeabilityCell=(_donutPuddle.m_permeabilityCell);
    m_timeStepRatio=(_donutPuddle.m_timeStepRatio);
    m_dt_NutUptake=(_donutPuddle.m_dt_NutUptake);
    
    _SimuRandom_=_donutPuddle._SimuRandom_;
    fileoutcompartments=_donutPuddle.fileoutcompartments;
    fileoutcells=_donutPuddle.fileoutcells;
    if (_donutPuddle.m_adjacencyNeed){
        m_puddle=new Compartment*[m_lengthInside+2];
        for (int i=0;i<=m_lengthInside+1;i++){
            m_puddle[i]=new Compartment[m_widthInside+2];
        }
    }
    
    else{
        m_puddle=new Compartment*[m_lengthInside+1];
        m_puddle[0]=new Compartment[1];
        for (int i=1;i<=m_lengthInside;i++){
            m_puddle[i]=new Compartment[m_widthInside+1];
        }
    }
    
    /* 2.Inital definition of compartments */
    if (_donutPuddle.m_adjacencyNeed){// Why is there no outside compartment here?
        for (int i=0;i<=m_lengthInside+1;i++){
            for (int j=0;j<=m_widthInside+1;j++){
                m_puddle[i][j]=_donutPuddle.m_puddle[i][j];
                m_puddle[i][j].SetReference(i,j);
            }
        }
    }
    else{
        m_puddle[0][0]=Compartment(); //is this outside compartment very different than the inside compartment
        for (int i=1;i<=m_lengthInside;i++){
            for (int j=1;j<=m_widthInside;j++){
                m_puddle[i][j]=_donutPuddle.m_puddle[i][j];
                m_puddle[i][j].SetReference(i,j);
            }
        }
    }
    /* 3.Definition of the first cells in the puddle */
    for (int n_c=0;n_c<_donutPuddle.m_cell.size();n_c++){
        m_cell.push_back(_donutPuddle.m_cell[n_c]);
        std::vector <int> newPosition;
        newPosition.push_back(_donutPuddle.m_positionCell[n_c][0]);//
        newPosition.push_back(_donutPuddle.m_positionCell[n_c][1]);
        m_positionCell.push_back(newPosition);
    }
    
    /*for (int n_c=0;n_c<initialNumberCell;n_c++){
     m_cell[n_c].EvaluateBSTFConnectivity();//Need to test this at each beginning of life and nowhere in between
     }*/
    
    /* 4.Definition of the spatio-temporal diffusivity coefficients for one timestep */
    
    //i)Inside diffusivity //same diffusivity for each nutrient
    m_alphaDF_empty[0][0]=_donutPuddle.m_alphaDF_empty[0][0];     //Diffusivity of standard nutrient of any compartment; could have been defined ahead of the loops on the compartments because the compartments have the same physical properties.
    m_alphaDF_empty[1][0]=_donutPuddle.m_alphaDF_empty[1][0] ;    //Diffusivity of hidden nutrient of any compartment
    
    //i)Outside diffusivity // via the surface of the squares separated by their height
    m_alphaDF_empty[0][1]=_donutPuddle.m_alphaDF_empty[0][1] ;    //Diffusivity of standard nutrient of any compartment; could have been defined ahead of the loops on the compartments because the compartments have the same physical properties.
    m_alphaDF_empty[1][1]=_donutPuddle.m_alphaDF_empty[1][1];
    
    
    
    m_neighbours=new int***[_donutPuddle.m_lengthInside+1];
    m_neighbours[0]=new int**[1];
    m_neighbours[0][0]=new int*[1];
    m_neighbours[0][0][0]=0;
    
    for (int i=1;i<=_donutPuddle.m_lengthInside;i++){
        m_neighbours[i]=new int**[_donutPuddle.m_widthInside+1];
        for (int j=1;j<=_donutPuddle.m_widthInside;j++){
            m_neighbours[i][j]=new int*[4];
            for (int k=0;k<=3;k++){
                m_neighbours[i][j][k]=new int[2];
            }
        }
    }
    
    for (int i=1;i<=_donutPuddle.m_lengthInside;i++){
        
        for (int j=1;j<=_donutPuddle.m_widthInside;j++){
            
            //First neighbour
            m_neighbours[i][j][0][0]=_donutPuddle.m_neighbours[i][j][0][0];
            m_neighbours[i][j][0][1]=_donutPuddle.m_neighbours[i][j][0][0];
            
            //Second neighbour
            m_neighbours[i][j][1][0]=_donutPuddle.m_neighbours[i][j][1][0];
            m_neighbours[i][j][1][1]=_donutPuddle.m_neighbours[i][j][1][1];
            
            //Third neighbour
            m_neighbours[i][j][2][0]=_donutPuddle.m_neighbours[i][j][2][0];
            m_neighbours[i][j][2][1]=_donutPuddle.m_neighbours[i][j][2][1];
            
            
            //Fourth neighbour
            m_neighbours[i][j][3][0]=_donutPuddle.m_neighbours[i][j][3][0];
            m_neighbours[i][j][3][1]=_donutPuddle.m_neighbours[i][j][3][1];
            
            //is it reallyuseful? How much time is lost here?
            m_puddle[i][j].SetNeighboursRelationships(0, m_neighbours[i][j][0][0], m_neighbours[i][j][0][1]);
            m_puddle[i][j].SetNeighboursRelationships(1, m_neighbours[i][j][1][0], m_neighbours[i][j][1][1]);
            m_puddle[i][j].SetNeighboursRelationships(2, m_neighbours[i][j][2][0], m_neighbours[i][j][2][1]);
            m_puddle[i][j].SetNeighboursRelationships(3, m_neighbours[i][j][3][0], m_neighbours[i][j][3][1]);
        }
            //Setting the spatio temporal diffusion coefficient for each of the molecules in each compartment.
    }
}

void DonutPuddle::operator=(const DonutPuddle &_donutPuddle)
{
    m_lengthInside=(_donutPuddle.m_lengthInside);
    m_widthInside=(_donutPuddle.m_lengthInside);
    m_compartmentVolume=(_donutPuddle.m_compartmentVolume);
    m_adjacencyNeed=(_donutPuddle.m_adjacencyNeed);
    m_cellDiffusionRate=(_donutPuddle.m_cellDiffusionRate);
    m_nutrientSwitchProbability=(_donutPuddle.m_nutrientSwitchProbability);
    m_Nutrient1Input=(_donutPuddle.m_Nutrient1Input);
    m_Nutrient2Input=(_donutPuddle.m_Nutrient2Input);
    m_permeabilityCell=(_donutPuddle.m_permeabilityCell);
    m_timeStepRatio=(_donutPuddle.m_timeStepRatio);
    m_dt_NutUptake=(_donutPuddle.m_dt_NutUptake);
    
    _SimuRandom_=_donutPuddle._SimuRandom_;
    fileoutcompartments=_donutPuddle.fileoutcompartments;
    fileoutcells=_donutPuddle.fileoutcells;
    if (_donutPuddle.m_adjacencyNeed){
        m_puddle=new Compartment*[m_lengthInside+2];
        for (int i=0;i<=m_lengthInside+1;i++){
            m_puddle[i]=new Compartment[m_widthInside+2];
        }
    }
    
    else{
        m_puddle=new Compartment*[m_lengthInside+1];
        m_puddle[0]=new Compartment[1];
        for (int i=1;i<=m_lengthInside;i++){
            m_puddle[i]=new Compartment[m_widthInside+1];
        }
    }
    
    /* 2.Inital definition of compartments */
    if (_donutPuddle.m_adjacencyNeed){// Why is there no outside compartment here?
        for (int i=0;i<=m_lengthInside+1;i++){
            for (int j=0;j<=m_widthInside+1;j++){
                m_puddle[i][j]=_donutPuddle.m_puddle[i][j];
                m_puddle[i][j].SetReference(i,j);
            }
        }
    }
    else{
        m_puddle[0][0]=_donutPuddle.m_puddle[0][0]; //is this outside compartment very different than the inside compartment
        for (int i=1;i<=m_lengthInside;i++){
            for (int j=1;j<=m_widthInside;j++){
                m_puddle[i][j]=_donutPuddle.m_puddle[i][j];
                m_puddle[i][j].SetReference(i,j);
            }
        }
    }
    /* 3.Definition of the first cells in the puddle */
    for (int n_c=0;n_c<_donutPuddle.m_cell.size();n_c++){
        m_cell.push_back(_donutPuddle.m_cell[n_c]);
        std::vector <int> newPosition;
        newPosition.push_back(_donutPuddle.m_positionCell[n_c][0]);//
        newPosition.push_back(_donutPuddle.m_positionCell[n_c][1]);
        m_positionCell.push_back(newPosition);
    }
    
    /*for (int n_c=0;n_c<initialNumberCell;n_c++){
     m_cell[n_c].EvaluateBSTFConnectivity();//Need to test this at each beginning of life and nowhere in between
     }*/
    
    /* 4.Definition of the spatio-temporal diffusivity coefficients for one timestep */
    
    //i)Inside diffusivity //same diffusivity for each nutrient
    m_alphaDF_empty[0][0]=_donutPuddle.m_alphaDF_empty[0][0];     //Diffusivity of standard nutrient of any compartment; could have been defined ahead of the loops on the compartments because the compartments have the same physical properties.
    m_alphaDF_empty[1][0]=_donutPuddle.m_alphaDF_empty[1][0] ;    //Diffusivity of hidden nutrient of any compartment
    
    //i)Outside diffusivity // via the surface of the squares separated by their height
    m_alphaDF_empty[0][1]=_donutPuddle.m_alphaDF_empty[0][1] ;    //Diffusivity of standard nutrient of any compartment; could have been defined ahead of the loops on the compartments because the compartments have the same physical properties.
    m_alphaDF_empty[1][1]=_donutPuddle.m_alphaDF_empty[1][1];

    m_neighbours=new int***[_donutPuddle.m_lengthInside+1];
    m_neighbours[0]=new int**[1];
    m_neighbours[0][0]=new int*[1];
    m_neighbours[0][0][0]=0;
    
    for (int i=1;i<=_donutPuddle.m_lengthInside;i++){
        m_neighbours[i]=new int**[_donutPuddle.m_widthInside+1];
        for (int j=1;j<=_donutPuddle.m_widthInside;j++){
            m_neighbours[i][j]=new int*[4];
            for (int k=0;k<=3;k++){
                m_neighbours[i][j][k]=new int[2];
            }
        }
    }
    
    for (int i=1;i<=_donutPuddle.m_lengthInside;i++){
        
        for (int j=1;j<=_donutPuddle.m_widthInside;j++){
            
            //First neighbour
            m_neighbours[i][j][0][0]=_donutPuddle.m_neighbours[i][j][0][0];
            m_neighbours[i][j][0][1]=_donutPuddle.m_neighbours[i][j][0][0];
            
            //Second neighbour
            m_neighbours[i][j][1][0]=_donutPuddle.m_neighbours[i][j][1][0];
            m_neighbours[i][j][1][1]=_donutPuddle.m_neighbours[i][j][1][1];
            
            //Third neighbour
            m_neighbours[i][j][2][0]=_donutPuddle.m_neighbours[i][j][2][0];
            m_neighbours[i][j][2][1]=_donutPuddle.m_neighbours[i][j][2][1];
            
            
            //Fourth neighbour
            m_neighbours[i][j][3][0]=_donutPuddle.m_neighbours[i][j][3][0];
            m_neighbours[i][j][3][1]=_donutPuddle.m_neighbours[i][j][3][1];
            
            //is it reallyuseful? How much time is lost here?
            m_puddle[i][j].SetNeighboursRelationships(0, m_neighbours[i][j][0][0], m_neighbours[i][j][0][1]);
            m_puddle[i][j].SetNeighboursRelationships(1, m_neighbours[i][j][1][0], m_neighbours[i][j][1][1]);
            m_puddle[i][j].SetNeighboursRelationships(2, m_neighbours[i][j][2][0], m_neighbours[i][j][2][1]);
            m_puddle[i][j].SetNeighboursRelationships(3, m_neighbours[i][j][3][0], m_neighbours[i][j][3][1]);
        }
        //Setting the spatio temporal diffusion coefficient for each of the molecules in each compartment.
    }
}

DonutPuddle::~DonutPuddle()
{//std::cout << "TET3" << endl;
    
    delete m_neighbours[0][0][0];
    delete m_neighbours[0][0];
    delete m_neighbours[0];
    for (int i=1;i<=m_lengthInside;i++){
        for (int j=1;j<=m_widthInside;j++){
            for (int k=0;k<=3;k++){
                delete[] m_neighbours[i][j][k];
            }
            delete[] m_neighbours[i][j];
        }
        delete[] m_neighbours[i];
    }
    delete[] m_neighbours;
    
}

//RAJON'S DONUT SEEMS TO BE CORRECTLY DEFINED. OH DEAR...//

void DonutPuddle::DiffusionBetweenCompartments()
{
    //#I. Molecules diffusion between compartments# //
    // 1.Initialization //
    double conc[2][m_lengthInside+1][m_widthInside+1];//here, may need a pointer for the dynamic assignation of environment molecules concentrations
    
    //int neigh[4][2];//Only useful in the test version, based on the neighbours recovery//

    
    //First kind of switch in the environment: the sources always bring the same amount of nutrients, albeit different from one another
    int switchTest=gsl_ran_bernoulli(_SimuRandom_, m_nutrientSwitchProbability);
    //std::cout << "nutinput: " << m_Nutrient1Input << endl;
    if (switchTest){
        //std::cout << "YES" << std::endl;
        //std::cout << "Concentration inside puddle test: " <<m_puddle[0][0].m_concentrationNut2 << std::endl;
        /*conc[0][0][0]=m_puddle[0][0].m_concentrationNut2;
        conc[1][0][0]=m_puddle[0][0].m_concentrationNut1;
        m_puddle[0][0].m_concentrationNut1=conc[0][0][0];
        m_puddle[0][0].m_concentrationNut2=conc[1][0][0];*/
        double Nutrient1_Input=m_Nutrient2Input;
        double Nutrient2_Input=m_Nutrient1Input;
        m_Nutrient1Input=Nutrient1_Input;
        m_Nutrient2Input=Nutrient2_Input;
    }
    else{
        //std::cout << "Concentration inside puddle test: " <<m_puddle[0][0].m_concentrationNut2 << std::endl;
        /*conc[0][0][0]=m_puddle[0][0].m_concentrationNut1;
        conc[1][0][0]=m_puddle[0][0].m_concentrationNut2;*/
    }
    
    /*
    for (int i=1;i<=m_lengthInside;i++)
    {
        for (int j=1;j<=m_widthInside;j++)
            
        {
            //Setting the prior concentration for each of the molecule type.
            
            conc[0][i][j]=m_puddle[i][j].m_concentrationNut1;
            
            conc[1][i][j]=m_puddle[i][j].m_concentrationNut2;
        }
    }*/
    
    //std::cout << m_puddle[0][0].m_concentrationNut2 << std::endl;
    //std::cout << conc[1][0][0] << std::endl;
    
    //cout << "conc pré: " << conc [0][1][1] << endl;
    // 2.Calculus of the new concentrations after diffusion inside //
    if(m_lengthInside>1 || m_widthInside>1){
        for (int i=1;i<=m_lengthInside;i++)
        {
            for (int j=1;j<=m_widthInside;j++)
            
            {
                double concNut1(OperatorFiniteDifferencesExplicit(m_alphaDF_actual[0][0],m_alphaDF_actual[0][1],conc[0][i][j],4,1,conc[0][m_neighbours[i][j][0][0]][m_neighbours[i][j][0][1]],conc[0][m_neighbours[i][j][1][0]][m_neighbours[i][j][1][1]],conc[0][m_neighbours[i][j][2][0]][m_neighbours[i][j][2][1]],conc[0][m_neighbours[i][j][3][0]][m_neighbours[i][j][3][1]],conc[0][0][0],0));
            
                double concNut2(OperatorFiniteDifferencesExplicit(m_alphaDF_actual[1][0],m_alphaDF_actual[1][1],conc[1][i][j],4,1,conc[1][m_neighbours[i][j][0][0]][m_neighbours[i][j][0][1]],conc[1][m_neighbours[i][j][1][0]][m_neighbours[i][j][1][1]],conc[1][m_neighbours[i][j][2][0]][m_neighbours[i][j][2][1]],conc[1][m_neighbours[i][j][3][0]][m_neighbours[i][j][3][1]],conc[1][0][0],0));
            
            
                m_puddle[i][j].SetConcentrationNut1(concNut1);
                m_puddle[i][j].SetConcentrationNut2(concNut2);
            }
        }
        //#I. Cells diffusion between compartments# //
        vector <int> cellDiffusionTemporaryArray;
        
        for (int cell=1; cell<=m_cell.size();cell++){
            if (gsl_ran_bernoulli(_SimuRandom_,m_cellDiffusionRate)==1){
                cellDiffusionTemporaryArray.push_back(cell);
            //cout << cellDiffusionTemporaryArray[0];
            }
        }
        for (int c=0; c<cellDiffusionTemporaryArray.size();c++){
            int cellDiffused=cellDiffusionTemporaryArray[c]-1;
            int cellNeighbourChoice;
            cellNeighbourChoice=gsl_ran_flat(_SimuRandom_,0,4);
            for (int coordinate=0;coordinate<=1;coordinate++){
                        m_positionCell[cellDiffused][coordinate]=m_neighbours[m_positionCell[cellDiffused][0]][m_positionCell[cellDiffused][1]][cellNeighbourChoice][coordinate];
            }
        }
    }
    else{
        /*
        //std::cout << "alpha empty: " << m_alphaDF_empty[0][1] << endl;
        //std::cout << "alpha full: " << m_alphaDF_actual[0][1] << endl;
        double concNut1(OperatorFiniteDifferencesExplicit(0,m_alphaDF_actual[0][1],conc[0][1][1],0,1,0,0,0,0,conc[0][0][0],0));
        double concNutTest(OperatorFiniteDifferencesExplicit(0,m_alphaDF_empty[0][1],conc[0][1][1],0,1,0,0,0,0,conc[0][0][0],0));
        //std::cout << "conc full:" << concNut1 << endl;
        //std::cout << "conc test: " << concNutTest << endl;
        
        double concNut2(OperatorFiniteDifferencesExplicit(0,m_alphaDF_empty[0][1],conc[1][1][1],0,1,0,0,0,0,conc[1][0][0],0));
        
        m_puddle[1][1].SetConcentrationNut1(concNut1);
        m_puddle[1][1].SetConcentrationNut2(concNut2);
         */ //In any case, there is a need to introduce a compratment in which the input occurs through net rate, for instance net release of glucose by another organism.

        //std::cout << "volume:" << m_puddle[1][1].GetActualVolume() << endl;
        //std::cout << "alpha full: " << m_alphaDF_actual[1][1] << endl;
        //std::cout << "conc1: " << m_puddle[1][1].m_concentrationNut1 << endl;
        //std::cout << "Concentration Nutrient 1: " << m_puddle[1][1].m_concentrationNut1 << endl;
        //std::cout << "Concentration Output:" << m_alphaDF_empty[1][1]*m_puddle[1][1].m_concentrationNut1 << endl;
        m_puddle[1][1].m_concentrationNut1-=m_alphaDF_empty[1][1]*m_puddle[1][1].m_concentrationNut1;//this rate should be constant because random walks should trigger increase in pressure towards the outside as the concentration increases.
        m_puddle[1][1].m_concentrationNut2-=m_alphaDF_empty[1][1]*m_puddle[1][1].m_concentrationNut2;
        //std::cout <<"alpha_prop:" << (m_alphaDF_actual[1][1]/m_alphaDF_empty[1][1]) << endl;
        m_puddle[1][1].m_concentrationNut1+=(m_alphaDF_actual[1][1]/m_alphaDF_empty[1][1])*m_Nutrient1Input/(m_puddle[1][1].GetActualVolume());
        m_puddle[1][1].m_concentrationNut2+=(m_alphaDF_actual[1][1]/m_alphaDF_empty[1][1])*m_Nutrient2Input/(m_puddle[1][1].GetActualVolume());
        
        //std:: cout << "Actual volume: " << m_puddle[1][1].GetActualVolume() << endl;
        //std::cout << "Concentration Input:" << (m_alphaDF_actual[1][1]/m_alphaDF_empty[1][1])*m_Nutrient1Input/(m_puddle[1][1].GetActualVolume()) << std::endl;
        //slowing down of each of the input and output when overcrowding
    }
    
    /*
     for(double t_nut=0;t_nut<0.999;t_nut+=m_dt_NutUptake){
     concNutEnvPre[0]=m_puddle[compartmentLengthPosition][compartmentWidthPosition].m_concentrationNut1;
     concNutEnvPre[1]=m_puddle[compartmentLengthPosition][compartmentWidthPosition].m_concentrationNut2;
     
     for (int cell=0; cell<numberOfCellsBeforeTimestep;cell++){
     if (m_positionCell[cell][0]==compartmentLengthPosition && m_positionCell[cell][1]==compartmentWidthPosition){
     //2.Nutrient uptake//
     //std::cout << "Cell n°: " << cell << std::endl;
     //std::cout << "t=:" << t_current+t_nut << std::endl;
     if (m_compartmentVolume>0){
     NutrientDiffusion(cell, concNutEnvPre,m_compartmentVolume);
     }
     else{
     //cout << "ERROR!!!!!" << endl;
     
     }
     }
     }
     //cout << "concentration within cell 1 " << m_cell[1].GetConcentrationNut1() << endl;
     //std::cout << "Concentration Nutrient 1: " << m_puddle[1][1].m_concentrationNut1 << endl;
     }*/
}
