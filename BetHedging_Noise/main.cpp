//
//  main.cpp
//  BetHedging_Noise
//
//  Created by Florian on 28/05/2019.
//  Copyright © 2019 Florian. All rights reserved.
//

#include <iostream>
#include "Simulation.hpp"
#include <math.h>
using namespace std;

int main(int argc, const char * argv[]) {
    
    /* Initialization: Defining the localization of the directory in which parameters are to be found */
    
    //string cdir="/beegfs/home/flabourel/BetHedging/bin/IO_files/";
    string cdir="/Users/florian/Desktop/p_Modelization/Models/Modele_BH/BetHedging_Noise/Results/";
    string filename;
    filename=cdir+"param.txt";
    ifstream filein(filename.c_str());
    
    if (filein){
        cout << "YES!" << endl;
    } //let us know if the file opening has worked//
    
    filein.ignore(10000000, '\n');
    
    int line_read;
    if(argc>1){
        istringstream iss1(argv[1]);
        iss1>>line_read;
    }
    else line_read=10;
    int current_line=1;
     
    while(current_line<line_read){
        filein.ignore(10000000, '\n');
        current_line++;
    }
    

    int line,simulationNumber;//, simulationNumber; //line for the number of the run, simulationNumber for the number of the simulation within each set of parameters
    
    //double cellPermeability;
    //double NutMaxInput;
    double averageTimeSwitch;
    double NutrientRate;
    
    
    //Communicating the information of the .txt to the simulation about to run//
    filein>>line>>simulationNumber>>averageTimeSwitch>>NutrientRate;
    
    //Printing the information about the running simulation//
    std::cout << "line: " << line << endl;
    cout << "n°: " << simulationNumber << endl;
    //cout << "Permeability: " << cellPermeability << endl;
    //cout << "Nutrient Input: " << NutMaxInput << endl;
    cout << "Switching time: " << averageTimeSwitch << endl;
    cout << "NutrientRatio: " <<  NutrientRate << endl;
    
    //Opening the file out//
    stringstream ss1;
    ss1<<line;
    string string1;
    ss1>>string1;
    
    double switchProbability=1/averageTimeSwitch;
    
    //string string1("12");
    //string filename;
    
    ofstream fileoutcompartments, fileoutcells, flieoutsetup;
    
    /*First output file: compartment*/
    filename=cdir+"line"+string1+"_compartments"+".rj";
    fileoutcompartments.open(filename.c_str());
    fileoutcompartments<<"simul\t t \t FreeVolFrac \t [Nutrient1] \t [Nutrient2] \t Popsize"<<endl;
    
    /*Second output file: cell*/
    filename=cdir+"line"+string1+"_cells.rj";
    fileoutcells.open(filename.c_str());
    fileoutcells<<"simul \t t \t Generation number \t CellSize\t MitosisSize \t Nut1Proteins \t Nut2proteins \t GenomeBarCode"<<endl;
    
    Simulation S(simulationNumber,500,2,10,1,switchProbability,NutrientRate,5*pow(10,-4),0,1,&fileoutcompartments,&fileoutcells);
    S.SimulationProcessing(100*pow(10,6));
}
