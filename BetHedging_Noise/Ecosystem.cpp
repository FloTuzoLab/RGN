//
//  Ecosystem.cpp
//  BetHedging_Noise
//
//  Created by Florian on 01/07/2019.
//  Copyright © 2019 Florian. All rights reserved.
//

#include "Ecosystem.hpp"

using namespace std;

/* ###Methods### */


/* ##Constructors## */

Ecosystem::Ecosystem():m_lengthInside(0),m_widthInside(0),m_compartmentVolume(0){
    
    m_adjacencyNeed=false;
    if (m_adjacencyNeed){
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
    if (m_adjacencyNeed){// Why is there no outside compartment here?
        for (int i=0;i<=m_lengthInside+1;i++){
            for (int j=0;j<=m_widthInside+1;j++){
                m_puddle[i][j]=Compartment(Compartment());
                m_puddle[i][j].SetReference(i,j);
            }
        }
    }
    else{
        m_puddle[0][0]=Compartment(Compartment()); //is this outside compartment very different than the inside compartment
        for (int i=1;i<=m_lengthInside;i++){
            for (int j=1;j<=m_widthInside;j++){
                m_puddle[i][j]=Compartment(Compartment());
                m_puddle[i][j].SetReference(i,j);
            }
        }
    }
}

// #Total-free constructor# //

Ecosystem::Ecosystem(gsl_rng *random_Simulation,
                     double switchProbability,
                     double delta_t,
                     double delta_l,
                     int timeStepRatio,
                     int lengthPuddle,
                     int widthPuddle,
                     bool adjacencyNeed,
                     double Nut1Input,
                     double Nut2Input,
                     Compartment insideCompartment,
                     Compartment outsideCompartment,
                     int initialNumberCell,
                     Cell cell,
                     ofstream* _fileoutcompartments,
                     ofstream* _fileoutcells):
m_nutrientSwitchProbability(switchProbability),m_Nutrient1Input(Nut1Input),m_Nutrient2Input(Nut2Input),m_lengthInside(lengthPuddle),m_widthInside(widthPuddle),m_adjacencyNeed(adjacencyNeed),fileoutcompartments(_fileoutcompartments),fileoutcells(_fileoutcells),m_compartmentVolume(insideCompartment.Volume()),m_dt_NutUptake(0.25)
{
    /* Incipit.Initialization of the random number */
    
    _SimuRandom_=random_Simulation;
    m_timeStepRatio=timeStepRatio;
    double heightPuddle=insideCompartment.GetHeight();//Why would one use an accessor while it can be defined directly in the puddle?
    
    /* 1.Need or need not for a mesh with adjacent squares */
    if (m_adjacencyNeed){
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
    if (m_adjacencyNeed){// Why is there no outside compartment here?
        for (int i=0;i<=m_lengthInside+1;i++){
            for (int j=0;j<=m_widthInside+1;j++){
                m_puddle[i][j]=Compartment(insideCompartment);
                m_puddle[i][j].SetReference(i,j);
            }
        }
    }
    else{
        m_puddle[0][0]=Compartment(outsideCompartment); //is this outside compartment very different than the inside compartment
        for (int i=1;i<=m_lengthInside;i++){
            for (int j=1;j<=m_widthInside;j++){
                m_puddle[i][j]=Compartment(insideCompartment);
                m_puddle[i][j].SetReference(i,j);
            }
        }
    }
    /* 3.Definition of the first cells in the puddle */
    Cell initCell(cell);
    for (int n_c=0;n_c<initialNumberCell;n_c++){
        m_cell.push_back(initCell);
        std::vector <int> newPosition;
        newPosition.push_back(gsl_ran_flat(_SimuRandom_, 1, m_lengthInside+1));// 
        newPosition.push_back(gsl_ran_flat(_SimuRandom_,1,m_widthInside+1));
        m_positionCell.push_back(newPosition);
    }
    m_permeabilityCell=initCell.GetPermeabilityCell();
    
    
    
    for (int n_c=0;n_c<initialNumberCell;n_c++){
        m_cell[n_c].EvaluateBSTFConnectivity();//Need to test this at each beginning of life and nowhere in between
    }
    
    /* 4.Definition of the spatio-temporal diffusivity coefficients for one timestep */
    
    //i)Inside diffusivity //same diffusivity for each nutrient
    m_alphaDF_empty[0][0]=delta_t/pow(delta_l,2)*insideCompartment.GetDiffusionNutrients();      //Diffusivity of standard nutrient of any compartment; could have been defined ahead of the loops on the compartments because the compartments have the same physical properties.
    m_alphaDF_empty[1][0]=delta_t/pow(delta_l,2)*insideCompartment.GetDiffusionNutrients();      //Diffusivity of hidden nutrient of any compartment
    
    //i)Outside diffusivity // via the surface of the squares separated by their height
    m_alphaDF_empty[0][1]=delta_t/pow(heightPuddle,2)*insideCompartment.GetDiffusionNutrients();      //Diffusivity of standard nutrient of any compartment; could have been defined ahead of the loops on the compartments because the compartments have the same physical properties.
    m_alphaDF_empty[1][1]=delta_t/pow(heightPuddle,2)*insideCompartment.GetDiffusionNutrients();      //Diffusivity of hidden nutrient of any compartment
    
    m_cellDiffusionRate=delta_t/pow(delta_l,2)*insideCompartment.GetDiffusionCells();
}

Ecosystem::Ecosystem(const Ecosystem& _ecosystem):
    m_lengthInside(_ecosystem.m_lengthInside),
    m_widthInside(_ecosystem.m_lengthInside),
    m_compartmentVolume(_ecosystem.m_compartmentVolume),
    m_adjacencyNeed(_ecosystem.m_adjacencyNeed),
    m_cellDiffusionRate(_ecosystem.m_cellDiffusionRate),
    m_nutrientSwitchProbability(_ecosystem.m_nutrientSwitchProbability),
    m_Nutrient1Input(_ecosystem.m_Nutrient1Input),
    m_Nutrient2Input(_ecosystem.m_Nutrient2Input),
    m_permeabilityCell(_ecosystem.m_permeabilityCell),
    m_timeStepRatio(_ecosystem.m_timeStepRatio),
    m_dt_NutUptake(_ecosystem.m_dt_NutUptake)

{
    _SimuRandom_=_ecosystem._SimuRandom_;
    fileoutcompartments=_ecosystem.fileoutcompartments;
    fileoutcells=_ecosystem.fileoutcells;
    if (_ecosystem.m_adjacencyNeed){
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
    if (_ecosystem.m_adjacencyNeed){// Why is there no outside compartment here?
        for (int i=0;i<=m_lengthInside+1;i++){
            for (int j=0;j<=m_widthInside+1;j++){
                m_puddle[i][j]=_ecosystem.m_puddle[i][j];
                m_puddle[i][j].SetReference(i,j);
            }
        }
    }
    else{
        m_puddle[0][0]=Compartment(); //is this outside compartment very different than the inside compartment
        for (int i=1;i<=m_lengthInside;i++){
            for (int j=1;j<=m_widthInside;j++){
                m_puddle[i][j]=_ecosystem.m_puddle[i][j];
                m_puddle[i][j].SetReference(i,j);
            }
        }
    }
    /* 3.Definition of the first cells in the puddle */
    for (int n_c=0;n_c<_ecosystem.m_cell.size();n_c++){
        m_cell.push_back(_ecosystem.m_cell[n_c]);
        std::vector <int> newPosition;
        newPosition.push_back(_ecosystem.m_positionCell[n_c][0]);//
        newPosition.push_back(_ecosystem.m_positionCell[n_c][1]);
        m_positionCell.push_back(newPosition);
    }
    
    /*for (int n_c=0;n_c<initialNumberCell;n_c++){
        m_cell[n_c].EvaluateBSTFConnectivity();//Need to test this at each beginning of life and nowhere in between
    }*/
    
    /* 4.Definition of the spatio-temporal diffusivity coefficients for one timestep */
    
    //i)Inside diffusivity //same diffusivity for each nutrient
    m_alphaDF_empty[0][0]=_ecosystem.m_alphaDF_empty[0][0];     //Diffusivity of standard nutrient of any compartment; could have been defined ahead of the loops on the compartments because the compartments have the same physical properties.
    m_alphaDF_empty[1][0]=_ecosystem.m_alphaDF_empty[1][0] ;    //Diffusivity of hidden nutrient of any compartment
    
    //i)Outside diffusivity // via the surface of the squares separated by their height
    m_alphaDF_empty[0][1]=_ecosystem.m_alphaDF_empty[0][1] ;    //Diffusivity of standard nutrient of any compartment; could have been defined ahead of the loops on the compartments because the compartments have the same physical properties.
    m_alphaDF_empty[1][1]=_ecosystem.m_alphaDF_empty[1][1];    //Diffusivity of hidden nutrient of any compartment
    
}

void Ecosystem::operator=(const Ecosystem& _ecosystem){
    
    m_lengthInside=(_ecosystem.m_lengthInside);
    m_widthInside=(_ecosystem.m_lengthInside);
    m_compartmentVolume=(_ecosystem.m_compartmentVolume);
    m_adjacencyNeed=(_ecosystem.m_adjacencyNeed);
    m_cellDiffusionRate=(_ecosystem.m_cellDiffusionRate);
    m_nutrientSwitchProbability=(_ecosystem.m_nutrientSwitchProbability);
    m_Nutrient1Input=(_ecosystem.m_Nutrient1Input);
    m_Nutrient2Input=(_ecosystem.m_Nutrient2Input);
    m_permeabilityCell=(_ecosystem.m_permeabilityCell);
    m_timeStepRatio=(_ecosystem.m_timeStepRatio);
    m_dt_NutUptake=(_ecosystem.m_dt_NutUptake);
    
    _SimuRandom_=_ecosystem._SimuRandom_;
    fileoutcompartments=_ecosystem.fileoutcompartments;
    fileoutcells=_ecosystem.fileoutcells;
    if (_ecosystem.m_adjacencyNeed){
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
    if (_ecosystem.m_adjacencyNeed){// Why is there no outside compartment here?
        for (int i=0;i<=m_lengthInside+1;i++){
            for (int j=0;j<=m_widthInside+1;j++){
                m_puddle[i][j]=_ecosystem.m_puddle[i][j];
                m_puddle[i][j].SetReference(i,j);
            }
        }
    }
    else{
        m_puddle[0][0]=Compartment(); //is this outside compartment very different than the inside compartment
        for (int i=1;i<=m_lengthInside;i++){
            for (int j=1;j<=m_widthInside;j++){
                m_puddle[i][j]=_ecosystem.m_puddle[i][j];
                m_puddle[i][j].SetReference(i,j);
            }
        }
    }
    /* 3.Definition of the first cells in the puddle */
    for (int n_c=0;n_c<_ecosystem.m_cell.size();n_c++){
        m_cell.push_back(_ecosystem.m_cell[n_c]);
        std::vector <int> newPosition;
        newPosition.push_back(_ecosystem.m_positionCell[n_c][0]);//
        newPosition.push_back(_ecosystem.m_positionCell[n_c][1]);
        m_positionCell.push_back(newPosition);
    }
    
    /*for (int n_c=0;n_c<initialNumberCell;n_c++){
     m_cell[n_c].EvaluateBSTFConnectivity();//Need to test this at each beginning of life and nowhere in between
     }*/
    
    /* 4.Definition of the spatio-temporal diffusivity coefficients for one timestep */
    
    //i)Inside diffusivity //same diffusivity for each nutrient
    m_alphaDF_empty[0][0]=_ecosystem.m_alphaDF_empty[0][0];     //Diffusivity of standard nutrient of any compartment; could have been defined ahead of the loops on the compartments because the compartments have the same physical properties.
    m_alphaDF_empty[1][0]=_ecosystem.m_alphaDF_empty[1][0] ;    //Diffusivity of hidden nutrient of any compartment
    
    //i)Outside diffusivity // via the surface of the squares separated by their height
    m_alphaDF_empty[0][1]=_ecosystem.m_alphaDF_empty[0][1] ;    //Diffusivity of standard nutrient of any compartment; could have been defined ahead of the loops on the compartments because the compartments have the same physical properties.
    m_alphaDF_empty[1][1]=_ecosystem.m_alphaDF_empty[1][1];    //Diffusivity of hidden nutrient of any compartment
    
}

Ecosystem::~Ecosystem(){
    
    if (m_adjacencyNeed){
        for (int i=0;i<=m_lengthInside+1;i++){
            delete[] m_puddle[i];
        }
    }
    
    else{
        for (int i=0;i<=m_lengthInside;i++){
            delete[] m_puddle[i];
        }
    }
    delete[] m_puddle;
    //*/
}

//Printing the puddle features//
void Ecosystem::PrintEcosystem()
{
    /*
     cout<<"Outside puddle definition:"<<endl;
     m_puddle[0][0].PrintState();
     if (m_adjacencyNeed)
     {
     cout<<"Lateral puddle definition:"<<endl;
     m_puddle[0][1].PrintState();
     }
     cout<<"Length Puddle:"<<m_lengthInside<<endl;
     cout<<"Width puddle:"<<m_widthInside<<endl;
     
     for (int i=1;i<=m_lengthInside;i++)
     {
     for (int j=1;j<=m_widthInside;j++)
     {
     //cout<<"Type de compartiment du compartiment: ("<<i<<","<<j<<")"<<endl;
     m_puddle[i][j].PrintState();
     }
     }
     */
    cout << "Number of cells: " << m_cell.size() << endl;
    for (int c=1;c<=m_cell.size();c++)
    {
        cout << "Cell n°" << c << endl;
        m_cell[c-1].PrintCell();
        cout << "Position:" << m_positionCell[c-1][0] << "," << m_positionCell[c-1][1] << endl;
        cout << endl;
    }
}

/* #2.Birth-death definition# */

//i.Replication of a cell //is it supobtimal to do te loop in another function?
void Ecosystem::CellReproduction(int cell){
    if (m_cell[cell].ReplicationCheckpoint())
    {
        // #a. Cell reroduction # //
        // # Duplication of the cell
        Cell daughterCell(m_cell[cell]);
        m_cell.push_back(daughterCell);
        m_positionCell.push_back(m_positionCell[cell]); //The new cell takes the same location than that of its mother.
        // # Segregation of gene products # //
        m_cell[m_cell.size()-1].BuildDaughterCell(m_cell[cell]);
        // # Processing of survival to mitosis and checking mutations effect for both cells: half-conservative reproduction# //
        //No failure mitosis as the volume does not evolve
        m_cell[m_cell.size()-1].MutationProcessing();
        m_cell[cell].MutationProcessing();
        
    }
}

//ii.Death of a cell
bool Ecosystem::CellDeath(int cell,double delta_t){
    if (m_cell[cell].DeathTest(delta_t)){
        m_cell.erase(m_cell.begin()+cell); //cell removal
        m_positionCell.erase(m_positionCell.begin()+cell); // cell position removal
        return(true);
    }
    else{
        return(false);
    }
}

void Ecosystem::NutrientDiffusion(int cell, double* concNutEnvPre,double freeCompartmentVolume){
    double* concNutEnv;
    
    concNutEnv=new double[2]{m_puddle[m_positionCell[cell][0]][m_positionCell[cell][1]].m_concentrationNut1,m_puddle[m_positionCell[cell][0]][m_positionCell[cell][1]].m_concentrationNut2};
    //std::cout << "Compartment volume: " <<m_compartmentVolume << std::endl;
    m_cell[cell].SetNutrientsAfterDiffusion(freeCompartmentVolume,concNutEnv,concNutEnvPre, m_dt_NutUptake);
    m_puddle[m_positionCell[cell][0]][m_positionCell[cell][1]].m_concentrationNut1=concNutEnv[0];
    m_puddle[m_positionCell[cell][0]][m_positionCell[cell][1]].m_concentrationNut2=concNutEnv[1];
    delete[] concNutEnv;
}

void Ecosystem::CellBiomoleculesDynamics(int cell,double delta_t_cellDyn){
    m_cell[cell].BasalRateDynamics(m_timeStepRatio,delta_t_cellDyn*m_timeStepRatio);
    m_cell[cell].TranscriptsDynamics(delta_t_cellDyn*m_timeStepRatio);
    m_cell[cell].ProteinsDynamics(delta_t_cellDyn*m_timeStepRatio); //fixed if need to see the effect of a constant content
}

void Ecosystem::CellEnergeticDynamics(int cell, double delta_t){
    m_cell[cell].EnergyProduction(delta_t);
    m_cell[cell].BasalMetabolismWithdrawal( delta_t);
    m_cell[cell].Growth();
}

void Ecosystem::TimeStepInsideEcosystem(int t_current, double delta_t){
    //std::cout << "Nut1Input: " << m_Nutrient1Input << endl;
    //Would be interesting to test if interesting to split processes instead of grouping them.//
    /*O. Whole compartment depletions calculus */ //Note: allows to use an analytical expression preventing numerical instabilities
    int numberOfCellsBeforeTimestep=(int) m_cell.size();
    double* concNutEnvPre;
    concNutEnvPre=new double[2]{0,0};
    //std::cout << "Concentration début: " << m_puddle[1][1].m_concentrationNut1 << endl;
    for (int cell=0; cell<numberOfCellsBeforeTimestep;cell++){
        /*1.Birth-death process*/
        if (CellDeath(cell,delta_t)){
            cell--;
            numberOfCellsBeforeTimestep--;
        }
        /*else*/ if(m_cell[cell].ReplicationCheckpoint()){
            CellReproduction(cell);
            //cout << "tmin repro:" << t_current << endl; to determine fitness differences
        }
    }
    
    
        // #Total cell membrane surface in a compartment# //
    for (int compartmentLengthPosition=1;compartmentLengthPosition<=m_lengthInside;compartmentLengthPosition++){
        for (int compartmentWidthPosition=1;compartmentWidthPosition<=m_widthInside;compartmentWidthPosition++){
            m_compartmentVolume=m_puddle[compartmentLengthPosition][compartmentWidthPosition].Volume();
            for (int cell=0; cell<numberOfCellsBeforeTimestep;cell++){
                if (m_positionCell[cell][0]==compartmentLengthPosition && m_positionCell[cell][1]==compartmentWidthPosition){
                    m_compartmentVolume-=m_cell[cell].GetVolumeCell();
                }
            }
            m_puddle[compartmentLengthPosition][compartmentWidthPosition].SetActualVolume(m_compartmentVolume);
            if (m_compartmentVolume>0){
                m_alphaDF_actual[0][0]=m_alphaDF_empty[0][0]*(m_compartmentVolume/m_puddle[compartmentLengthPosition][compartmentWidthPosition].Volume());      //Diffusivity of standard nutrient of any compartment; could have been defined ahead of the loops on the compartments because the compartments have the same physical properties.
                m_alphaDF_actual[1][0]=m_alphaDF_empty[1][0]*(m_compartmentVolume/m_puddle[compartmentLengthPosition][compartmentWidthPosition].Volume());      //Diffusivity of hidden nutrient of any compartment
            
                //i)Outside diffusivity // via the surface of the squares separated by their height
                m_alphaDF_actual[0][1]=m_alphaDF_empty[0][1]*(m_compartmentVolume/m_puddle[compartmentLengthPosition][compartmentWidthPosition].Volume());      //Diffusivity of standard nutrient of any compartment; could have been defined ahead of the loops on the compartments because the compartments have the same physical properties.
                m_alphaDF_actual[1][1]=m_alphaDF_empty[1][1]*(m_compartmentVolume/m_puddle[compartmentLengthPosition][compartmentWidthPosition].Volume());      //Di
                //std::cout << "Volume compartiment: " << m_compartmentVolume << endl;
                
            }
            
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
            }
            for (int cell=0; cell<numberOfCellsBeforeTimestep;cell++){
                /*3.Energetics*/
                CellEnergeticDynamics(cell,delta_t);
                if (t_current%m_timeStepRatio==0){
                    //m_cell[cell].PrintCell();
                    CellBiomoleculesDynamics(cell,delta_t*m_timeStepRatio);
                }
            }
        }
    }
    
    delete[] concNutEnvPre;
}



/*Generating the file with results*/
void Ecosystem::OutputGetting(int simulationNumber,int t_output)
{
    for (int i=1;i<=m_lengthInside;i++)
    {
        for (int j=1;j<=m_widthInside;j++)
        {
            (*fileoutcompartments)
            <<simulationNumber<<"\t"
            <<t_output<<"\t"
            <<m_compartmentVolume<<"\t"
            <<m_puddle[i][j].m_concentrationNut1<<"\t"
            <<m_puddle[i][j].m_concentrationNut2<<"\t"
            <<m_cell.size()
            <<endl;
        }
    }
    for (int c=0;c<m_cell.size();c++)
    {
        /*
         vector <int> _numericGenomeBarCode;
         m_cell[c-1].MatrixConversion(_numericGenomeBarCode);
         stringstream ss;
         copy( _numericGenomeBarCode.begin(), _numericGenomeBarCode.end(), ostream_iterator<int>(ss, " "));
         string genomeBarCode = ss.str();
         genomeBarCode = genomeBarCode.substr(0, genomeBarCode.length()-1);
         */
        //cout << genomeBarCode << endl;
    (*fileoutcells)
        <<simulationNumber<<"\t"
        <<t_output<<"\t"
        <<m_cell[c].GetGenerationNumber()<<"\t"
        <<m_cell[c].GetVolumeCell()<<"\t"
        <<m_cell[c].GetMitosisSize()<<"\t"
        <<m_cell[c].GetProteins("Nut1")<<"\t"
        <<m_cell[c].GetProteins("Nut2")<<"\t"
        <<m_cell[c].printCellMatrix()
        << endl;
    }
}
