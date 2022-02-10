//
//  main.cpp
//  Cell_Evolution
//
//  Created by Florian on 16/02/2018.
//  Copyright © 2018 Florian. All rights reserved.
//

#include <fstream>
#include <sstream>
#include "Ecosystem.hpp"

using namespace std;

int main(int argc, char* argv[]) {
    
    /* Initialization: Defining the localization of the directory in which parameters are to be found */
     string cdir="/Users/florian/Desktop/Modelization/Models/Modele_multi_current/Cell_Evolution/";
     //string cdir="/panhome/flabourel/Cell_Evo_last/bin/IO_files/";
    
    string filename;
     filename=cdir+"param.txt";
     ifstream filein(filename.c_str());
    
    if (filein){
        cout << "YES!" << endl;
    } /*let us know if the file opening has worked*/
    
    filein.ignore(10000000, '\n');
    
    int line_read;
    if(argc>1){
        istringstream iss1(argv[1]);
        iss1>>line_read;
    }
    /* TEST fonction istringstream
    else{
        istringstream iss2("2 3");
        iss2>>line_read;
    }
    */
    else line_read=10;
    int current_line=1;
    
    while(current_line<line_read){
        filein.ignore(10000000, '\n');
        current_line++;
    }
    
    /*Incipit:variables possibly taken in outside file*/
    int line, simulationNumber; //line for the number of the run, simulationNumber for the number of the simulation within each set of parameters
    int sizerConcentrationNum;// Use of a concentration of bound sites to reproduce
    int sizerCycleRegulation;// Presence of a
    double eta_trans;// Seems like a maximum rate of transcription
    double cell_diff_ratio;// Diffusion coefficient ratio for cells: why a ratio
    double CG_diff_ratio;// Diffusion coefficient ratio for Common Goods : ratio to the Nurtients diffusion rate
    double desc_prop;// Asymmetry of division: mistake in the coefficient for the moment and may not be the better info to ask to the simulator

   /*Communicating the information of the .txt to the simulation about to run*/ filein>>line>>simulationNumber>>sizerConcentrationNum>>sizerCycleRegulation>>eta_trans>>cell_diff_ratio>>CG_diff_ratio>>desc_prop;
    
    /*Printing the information about the running simulation*/
    cout << "line: " << line << endl;
    cout << "n°: " << simulationNumber << endl;
    cout << "genepresence: " << sizerCycleRegulation << endl;
    cout << "sizerconcpres: " << sizerConcentrationNum << endl;
    cout << "eta_trans: " << eta_trans << endl;
    cout << "cell_diff_ratio: " << cell_diff_ratio << endl;
    cout << "CG_diff_ratio: " << CG_diff_ratio << endl;
    cout << "proportion descendance: " << desc_prop << endl;
    cout << "sizerconcentration: " << sizerConcentrationNum << endl;
    
    /*Definition of the time allocator*/
    time_t start, end;
    time(&start);
    
    /*Opening the file out*/
    stringstream ss1;
    ss1<<line;
    string string1;
    ss1>>string1;
    ofstream fileoutcompartments, fileoutcells, flieoutsetup;
    
    /*First output file: compartment*/
    filename=cdir+"line"+string1+"_compartments"+".rj";
    fileoutcompartments.open(filename.c_str());
    fileoutcompartments<<"simul\tt\t[nutst]\t[nuthid]\t[nutrev]\t[BC]\tpopsize"<<endl;
    
    /*Second output file: cell*/
    filename=cdir+"line"+string1+"_cells.rj";
    fileoutcells.open(filename.c_str());
    fileoutcells<<"simul\tt\tGeneration\tMitosisThreshold\tSize\tMitosisSize\tSizerproteins\tGrowthproteins\tMaintenanceproteins\tPublicGoodsproteins\tGenomeBarCode\tGeneCode"<<endl;
    
    /*Initializing the parameters of a simulation based on the features of the simulaiton defined in .txt*/
    //Remark: this part may be better understood if compartmentalized in another unit than in the main, for example in a simulaiton unit
    bool sizerConcentration;
    if (sizerConcentrationNum==1){
        sizerConcentration=true;
    }
    else{
        sizerConcentration=false;
    }
    
    double sizerThreshold;
    if (sizerConcentration){
        sizerThreshold=10;
    }
    else{
        sizerThreshold=100;
    }
    /*Time features of the simulation*/
    double T_simu(4*pow(10,7));
    double delta_t(0.5);
    int timeStepRatio(10); //Ratio between diffusion dynamics and genomic dynamics, example: 1s. for diffusion processes while 10s. for transcription
    
    
    int numberTF=14;
    double lengthEnvironment(2*pow(10,-2));
    double widthEnvironment(2*pow(10,-2));
    double publicGoodsDiffusivity (pow(10,-9)*CG_diff_ratio);
    double cellDiffusivity (pow(10,-11)*cell_diff_ratio);
    int lengthNumberCompartment(40);
    double commonGoodsConcentrationRatio(1);
    double mutationRate(pow(10,-4)); //mutation rate: per site per generation for regulatory sequences, per generation for efficiency mutations
    
    double ratioBetweenGrowthandCooperationCost(1); //ratio between the energetic cost of the growth cluster and that of Public Goods
    double ratioBetweenK_offAndK_on(500); //ratio between
    double allometryBasalAlpha(1.25); //allometry of basal metabolism
    double minimumSynthesisDuration(1500); //minimum time for ADN replication
    double ATPexpensesCalibration(pow(10,6)); //calibration of basal metabolism on a certain cell: 1micrometer^3?
    double volumeProportionDescendant(desc_prop); //
    double volumeCellMitosis(40); //Volume at mitosis in the case without Sizer protein
    
    double etaTranslationActivityPostCheckpoint(eta_trans);
    double K_commonGoodsKinetics(pow(10,5)); //Kinetics dynamic parameter: unit forgotten
    bool PGefficiencyMutation(false);// No loss of Common Goods cluster of gene possible
    
    // Generation of 20 regulatory sequences corresponding to 16 TF genes and 4 physiological genes//
    
    int lengthRegulatorySequences(1000);
    int lengthBindingSites(8);
    int numberTrait;
    if (sizerCycleRegulation==1){
        numberTrait=4;
    }
    else{
        numberTrait=3;
    }
    const int numberSequence=numberTF+numberTrait;
    string sequence[numberSequence];
    string bindingsite[numberTF];
    gsl_rng *Random=gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(Random, simulationNumber);
    for (int seq=0;seq<numberSequence;seq++){ //ETIENNE: OK; pourrait étre simplifié avec un char[4]={'A', 'T', ...} et gsl_ran_choose; drawResult pourrait être défini en dehors de la boucle
        for(int pos=0;pos<lengthRegulatorySequences;pos++){
            int drawResult=gsl_ran_flat(Random, 0, 4);
            if (drawResult<1){
                sequence[seq].push_back('A');
            }
            else if (drawResult<2){
                sequence[seq].push_back('C');
            }
            else if (drawResult<3){
                sequence[seq].push_back('G');
            }
            else{
                sequence[seq].push_back('T');
            }
        }
    }
    for (int TFnumber=0;TFnumber<numberTF;TFnumber++){
        for(int pos=0;pos<lengthBindingSites;pos++){
            int drawResult=gsl_ran_flat(Random, 0, 4);
            if (drawResult<1){
                bindingsite[TFnumber].push_back('A');
            }
            else if (drawResult<2){
                bindingsite[TFnumber].push_back('C');
            }
            else if (drawResult<3){
                bindingsite[TFnumber].push_back('G');
            }
            else{
                bindingsite[TFnumber].push_back('T');
            }
        }
    }
    
    //Redefinition of the first gene(Growth) in order to fix a starting genotype from pre-simulation
    
    for (int seq=0;seq<numberSequence;seq++){
        if (seq!=2 && seq!=3){
                sequence[seq].replace(0,bindingsite[1].length(),bindingsite[1]); //General activator
        }
    }
    //sequence[0].replace(150,bindingsite[4].length(),bindingsite[4]); //Extra positive regulation of growth
    sequence[1].replace(150,bindingsite[1].length(),bindingsite[5]); //Positive regulation of maintenance
    //sequence[1].replace(250,bindingsite[2].length(),bindingsite[2]);
    /*if (!sizerConcentration){
        sequence[3].replace(450,bindingsite[5].length(),bindingsite[5]);// Test si cela permet de faire survivre les populations avec sizer en concentration
    }*/
    //sequence[5].replace(250,bindingsite[2].length(),bindingsite[2]);
    //sequence[7].replace(150,bindingsite[2].length(),bindingsite[2]);
    //sequence[9].replace(250,bindingsite[2].length(),bindingsite[2]);
    //sequence[10].replace(250,bindingsite[2].length(),bindingsite[2]);
    //sequence[6].replace(250,bindingsite[2].length(),bindingsite[2]);
    
    
    /* Initialization of the objects in the simulations */
    
    //1.Medium definition//
    double nutrientDiffusivity (5*pow(10,-8));
    
    
    Medium medium(nutrientDiffusivity,publicGoodsDiffusivity,cellDiffusivity); //Diffusion coefficients of the medium; order: Nutrients, Public Goods, Cells
    
    //2.Compartments definition//
    
    double heightEnvironmment(5*pow(10,-4));
    
    double lengthCompartment=lengthEnvironment/lengthNumberCompartment;
    int widthNumberCompartment=widthEnvironment/lengthCompartment;
    double standardNutrientConcentration(2*pow(10,-7));
    if(etaTranslationActivityPostCheckpoint==0.2 || volumeProportionDescendant==0.2 || sizerConcentration){
        standardNutrientConcentration*=1;
    }
    double standardNutrientEquilibriumConcentration(standardNutrientConcentration*nutrientDiffusivity/pow(heightEnvironmment,2));
    double hiddenNutrientConcentration(commonGoodsConcentrationRatio*standardNutrientConcentration);
    double revealedNutrientConcentration(0);
    double publicGoodsConcentration(0);
    
    Compartiment InsideCompartment(lengthCompartment,heightEnvironmment,standardNutrientEquilibriumConcentration,hiddenNutrientConcentration,revealedNutrientConcentration,publicGoodsConcentration,medium,K_commonGoodsKinetics);//Concentrations in mol.L^(-1); order: Standard Nutrient, Hidden Nutrient, Revealed Nutrient, Common Goods
    Compartiment UpperOutsideCompartment(lengthCompartment,heightEnvironmment,standardNutrientConcentration,hiddenNutrientConcentration,revealedNutrientConcentration,publicGoodsConcentration,medium,K_commonGoodsKinetics);//Concentrations in mol.L^(-1) ; order: same as above
    Compartiment LateralOutsideCompartment(lengthCompartment,heightEnvironmment,standardNutrientConcentration,hiddenNutrientConcentration,revealedNutrientConcentration,publicGoodsConcentration,medium,K_commonGoodsKinetics);//Concentrations in mol.L^(-1); order: same as above

    //3.Alleles definition//
    double physiologicalMutationRate(mutationRate);
    double regulatoryMutationRate(mutationRate);
    bool geneEfficiency(true);
    
        //i)Physiological alleles//
    PhysiologicalAllele growthallele(sequence[0],physiologicalMutationRate,geneEfficiency,true,Random);          //sequence should result from pre-simulation, mutation rate as a variable parameter.
    PhysiologicalAllele maintenanceallele(sequence[1],physiologicalMutationRate,geneEfficiency,true,Random);     //sequence should result from pre-simulation, mutation rate as a variable parameter.
    PhysiologicalAllele publicgoodsallele(sequence[2],physiologicalMutationRate,geneEfficiency,PGefficiencyMutation,Random);     //sequence should result from pre-simulation, mutation rate as a variable parameter.
    
    std::vector <RegulatoryAllele> positivetranscriptionfactors;               //vector of positive regulatory alleles
    std::vector <RegulatoryAllele> negativetranscriptionfactors;               //vector of negative regulatory alleles
    
        //ii)Regulatory alleles//
    for (int TFnumber=0;TFnumber<numberTF/2;TFnumber++){
        RegulatoryAllele randompositiveTranscriptionFactor(sequence[TFnumber+numberTrait],regulatoryMutationRate,Random,bindingsite[TFnumber],'A');
        RegulatoryAllele randomnegativeTranscriptionFactor(sequence[TFnumber+numberTrait+numberTF/2],regulatoryMutationRate,Random,bindingsite[TFnumber+numberTF/2],'R');
        positivetranscriptionfactors.push_back(randompositiveTranscriptionFactor);
        negativetranscriptionfactors.push_back(randomnegativeTranscriptionFactor);
    }
    
        //iii)Genome assemblage and sizer mechanisms definition (artificial or in the RGN)//
    map <string,PhysiologicalAllele> physiologicalInitialAlleles;
    physiologicalInitialAlleles.insert(pair<string,PhysiologicalAllele>("Growth",growthallele));               //allele encoding growth protein
    physiologicalInitialAlleles.insert(pair<string,PhysiologicalAllele>("Maintenance",maintenanceallele));     //allele encoding maintenance protein
    physiologicalInitialAlleles.insert(pair<string,PhysiologicalAllele>("Public Goods",publicgoodsallele));    //allele encoding public goods protein
    double volumeMutationRate(0);
    if (sizerCycleRegulation==1){
        PhysiologicalAllele cellsizer(sequence[3],physiologicalMutationRate,geneEfficiency,true,Random);             //sequence should result from pre-simulation, mutation rate as a variable parameter.
        physiologicalInitialAlleles.insert(pair<string,PhysiologicalAllele>("Cell sizer",cellsizer));              //allele encoding sizer protein
        volumeMutationRate=(20*mutationRate);
    }
    else{
        volumeMutationRate=(20*mutationRate);
    }
    
    //4.Trait costs//
    map <string,double> metabolitesCost;
    double growthCost(2000);
    double maintenanceCost(200);
    double cooperationCost(growthCost*ratioBetweenGrowthandCooperationCost);
    metabolitesCost.insert(pair<string,double>("Growth",growthCost));
    metabolitesCost.insert(pair<string,double>("Maintenance",maintenanceCost));
    metabolitesCost.insert(pair<string,double>("Public Goods",cooperationCost));
    
    
    
    //5. Cell definition//
    double K_on(pow(10,2));
    double K_off(ratioBetweenK_offAndK_on*K_on);
    double cellPermeability(pow(10,-1));  // Cell permeability to oxygen (dm.s^-1)
    double allometryProductionAlpha(1);
    double stoechiometryNutrienttoGarbage(1);
    double tauTranscriptionOncePolymeraseBound(0.1);
    double basalBindingProbabilityPolymerase(0.01);
    double enhancedBindingProbabilityPolymerase(0.5);
    double tauTraduction(0.1);
    double tauDegradationTranscript(0.001);
    double tauDegradationProtein(0.0005);
    double SAcalibration(1);
    double initialCellVolume(30);
    double garbageProportionDescendant(0.5);
    string cellCycleType("GSM");
    
    Cell cell(K_on,                     //K_on BSTF
               K_off,                     //K_off BSTF
               cellPermeability,                    //permeabilityCell (unit: micrometer.s^-1, corresponds to estimated oxygen permeability)
               allometryProductionAlpha,                              //allometryProductionAlpha (not subject to sensitivity study: means that energy production is strictly proportional to nutrient uptake
               allometryBasalAlpha,                           //allometryBasalAlpha (subject to sensitivity study, from 1 to 1.5)
               stoechiometryNutrienttoGarbage,                              //etaNutrientToGarbage  (number of garbage produced by the conversion of one nutrient)
               minimumSynthesisDuration,                           //synthesisDuration    (duration of DNA synthesis) useless now
               tauTranscriptionOncePolymeraseBound,                            //tauTranscriptionOncePolymeraseBound  (transcription rate once polymerase is bound: acknowledged as a constant, cf. Bio by the number)
               basalBindingProbabilityPolymerase,                           //basalBindingProbabilityPolymerase   (polymerase binding probability without the help of an enhancer: chosen such that the complete range of values from 0 to maximum transcription rate can occur)
               enhancedBindingProbabilityPolymerase,                            //enhancedBindingProbabilityPolymerase (polymerase binding probability wuth the help of an enhancer: chosen such that the maximum transcription rate described in litterature can occur)
               tauTraduction,                            //tauTraduction (translation rate according to scientific litterature, cf. Bio. by the numbers, for one ribosome: generally more than one can be actively translating a single mRNA, function of the energy as a result of the basal metabolism investment (which is only depending on the volume here, not subject to any evolutionary process))
               etaTranslationActivityPostCheckpoint,
               tauDegradationTranscript,                          //tauDegradationTranscript (mRNA degradation rate according to scientific litterature, cf. Bio. by the numbers)
               tauDegradationProtein,                         //tauDegradationProtein  (protein degradation rate according to scientific litterature, cf. Bio. by the numbers)
               metabolitesCost,                //metabolites cost  (cost of producing each metabolite in terms of energy units: subject to wide sensitivity study because few knowledge in the litterature)
               SAcalibration,                              //Surface Area used for basal metabolism calibration
               ATPexpensesCalibration,                      //ATP expenses for the Surface Area used in calibration
               initialCellVolume,                             //cellVolume  (square micrometers)
               volumeProportionDescendant,                            //volumeProportionDescendant (percentage: subject to sensitivity study because we intuit that it could influence differentiation by increasing noise
               garbageProportionDescendant,                            //garbageProportionDescendant (percentage: not used in the first models)
               volumeCellMitosis,                            //volumeCellMitosis  (unused parameter: served at the beginning to parameterize the models)
               volumeMutationRate,
               sizerConcentration,
               sizerThreshold,                              //sizerThreshold   (concentration threshold allowing to progress in the cycle cycle towards DNA replication)
               cellCycleType,                          //type of cell cycle
               physiologicalInitialAlleles,
               positivetranscriptionfactors,
               negativetranscriptionfactors
               );
    
    //
    
    double initialCellNumbers(500);
    
    DonutPuddle donutpuddle(simulationNumber,     //seed initialization
                            sizerCycleRegulation, //set the mechanism of cell division and the resulting process of cell division size mutation
                            delta_t,              //timestep used in the simulation
                            lengthCompartment,            //length of the compartment (unit decimeters, cubic shape)
                            timeStepRatio,
                            lengthNumberCompartment,                   //length of the ecosystem (unit: number of compartment)
                            widthNumberCompartment,                   //width of the ecosystem (unit: number of compartment)
                            InsideCompartment,             //puddle inside compartment
                            UpperOutsideCompartment,             //puddle outside compartment
                            initialCellNumbers,                    //initial number of cells
                            cell,
                            &fileoutcompartments,
                            &fileoutcells
                            );               //initial cells in the simulation
    
    
    
    
    int simu_stepNumber=T_simu/delta_t;
    int t_output=20000;
    donutpuddle.OutputGetting(simulationNumber,1);
    for (int t=1;t<=simu_stepNumber;t++){
        /*In silico reality biological simulation*/
        donutpuddle.TimeStepInsidePuddle(t,delta_t);//Simulation inside each compartment
        donutpuddle.DiffusionBetweenCompartments(); //Diffusion between compartments
        /*Periodic output getting to follow evolutionary trajectories*/
        if (t==t_output){
            donutpuddle.OutputGetting(simulationNumber,t_output);
            t_output+=T_simu/80;
        }
                        /*Test to see if simulations ongoing does simulate the process we want*/
                        /*
                            if (t==t_int){
                         cout << "Pas de temps n°: " << t << endl;
            
                         t_int+=2000;
         
         
                         if (t==t_int_print){
                         donutpuddle.PrintPuddle();
                t_int_print+=2000;
                         }
         
        
                         cout << endl << endl;
                         cout << "Pas de temps n°: " << t << endl;
                         }
                         */
    }
    /*Final output getting: useless if already included in the regular tracking*/
    //donutpuddle.OutputGetting(simulationNumber,T_simu);

    gsl_rng_free(Random);//freeing memory of the random starting point associated with the simulation number
    time(&end);
    filein.close(); //Closing the output files
    
}
//}


    


