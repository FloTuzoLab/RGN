//
//  Simulation.cpp
//  BetHedging_Noise
//
//  Created by Florian on 02/07/2019.
//  Copyright © 2019 Florian. All rights reserved.
//

#include "Simulation.hpp"

Simulation::Simulation(){
    RandomSimulationNumber=gsl_rng_alloc(gsl_rng_default);
    //gsl_rng_set(RandomSimulationNumber,0);
}

Simulation::Simulation(int _simulationNumber,
                       int _initialCellAmount,
                       int _nutrientTypeAmount,
                       int _TFamount,
                       int _compartmentAmountLinear,
                       double _switchProbability,
                       double _nutrientConcentrationRatio, //from the dominant nutrient compared to the total concentrations
                       double _mutationRate,
                       double _volumeMutationRate,
                       double _aPropDesc,
                       std::ofstream* fileoutcompartments,
                       std::ofstream* fileoutcells):
m_nutrientConcentrationRatio(_nutrientConcentrationRatio),m_simulationNumber(_simulationNumber),m_simulationTimeStep(0),m_delta_t(1)//,m_delta_t_diff(0),

{
    RandomSimulationNumber=gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(RandomSimulationNumber, _simulationNumber);
    int initialCellAmount(_initialCellAmount);
    /*Time features of the simulation*/
    
    double delta_l=pow(250,1.0/3.0);
    std::cout << "DELTAL:" << delta_l << endl;
    
    double lengthEnvironment(delta_l*pow(10,1));//in order to set the environment size to 50micro.meter^3
    double widthEnvironment(delta_l*pow(10,1));
    
    //double mutationRate(_mutationRate); //mutation rate: per site per generation for regulatory sequences, per generation for efficiency mutations
    double K_on(pow(10,6));//mol.s^-1 K_on=pow(10,7) to test
    double K_off(pow(10,-2));//sec^-1 >> conversion needed to be used in the cell object K_off=pow(10,-2) seems OK cf. (Donovan et al., 2018)
    //Decreasing K_on could be a realistic way to slow down the switch between cascades
    double kcat(100);//little overestimated (in litterature: around 30-50)
    double KM(pow(10,-3));//very variable among sugars: weak for glucose (10^-5) and high for sucrose (10^-3>>5*10^-3)
    
    //double volumeProportionDescendant(_aPropDesc);
    double volumeCellMitosis(80);
    
    int lengthRegulatorySequences(1000);
    int lengthBindingSites(8);
    
    int traitAmount(2);//Need to be changed if more traits
    int TFamount(_TFamount);
    const int sequenceAmount=traitAmount+TFamount;
    
    string sequence[sequenceAmount];
    string bindingsite[TFamount];
    
    // Generation of sequenceAmount regulatory sequences corresponding to 16 TF genes and 4 physiological genes//
    char nucleotideResult[1];
    char nucleotideAlphabet[4]={'A','T','G','C'};
    
    for (int seq=0;seq<sequenceAmount;seq++){
        
        for(int pos=0;pos<lengthRegulatorySequences;pos++){
            
            gsl_ran_choose(RandomSimulationNumber,nucleotideResult,1,nucleotideAlphabet,4,sizeof(char));
            sequence[seq].push_back(nucleotideResult[0]);
        }
    }
    
    // Generation of TFamount binding sequences
    for (int TFnumber=0;TFnumber<TFamount;TFnumber++){
        for(int pos=0;pos<lengthBindingSites;pos++){
           
            gsl_ran_choose(RandomSimulationNumber,nucleotideResult,1,nucleotideAlphabet,4,sizeof(char));
            bindingsite[TFnumber].push_back(nucleotideResult[0]);
            
        }
    }
    
    //Starting point: positive regulation of nutrient 1 and 2: generalist phenotype//
    /*for (int i=0; i<TFamount;i++){
        sequence[i].replace(200,bindingsite[3].length(),bindingsite[3]);
        
    }*/
    sequence[0].replace(0,bindingsite[0].length(),bindingsite[0]); //S1 enhanced by A1
    sequence[1].replace(0,bindingsite[1].length(),bindingsite[1]); //S2 enhanced by A2
    /*
    sequence[2].replace(0,bindingsite[5].length(),bindingsite[5]);// A1 inhibited by R1
    sequence[3].replace(0,bindingsite[6].length(),bindingsite[6]);//A2 inhibited by R2
    
    sequence[2].replace(100,bindingsite[0].length(),bindingsite[0]);// A1 inhibited by R1
    sequence[3].replace(100,bindingsite[1].length(),bindingsite[1]);//A2 inhibited by R2
     */
    /*
    sequence[0].replace(100,bindingsite[5].length(),bindingsite[5]);//S1 inhibited by R1
    sequence[1].replace(100,bindingsite[6].length(),bindingsite[6]);//S2 inhibite by R2
    sequence[7].replace(0,bindingsite[1].length(),bindingsite[1]);//R1 enhanced by A2
    sequence[8].replace(0,bindingsite[0].length(),bindingsite[0]);//R2 enhanced by A1
    */
    //Need to test the effect: i there a need for an already existing positive feedback on TF{0] defined hereafter);
    /*sequence[2].replace(0,bindingsite[0].length(),bindingsite[0]);*/
    
    //Another strating point: general activator leaving only Nutrient 2 genes with no enhancing activity
    /*for (int seq=0;seq<sequenceAmount;seq++){
        if (seq!=1){
           sequence[seq].replace(0,bindingsite[0].length(),bindingsite[0]); //General activator
        }
    }*/
    
    
    /* Initialization of the objects in the simulations */
    double NutrientInputInit(25);
    double cellPermeability(1*pow(10,-2));  // Cell permeability to glucose (micrometer.s^-1) -- including some of it through passive diffusion. Above [Nut]=10^-2, passive diffusion is enough for the cell to be fed.
    
    //1.Medium definition//
    double nutrientDiffusivity (5*pow(10,2));//unit:micrometer^2/s.
    double Nut1Input(_nutrientConcentrationRatio*(NutrientInputInit*pow(10,1)));//unit:mol.dm^-3
    double Nut2Input((1-_nutrientConcentrationRatio)*(NutrientInputInit*pow(10,1)));//unit:mol.dm^-3//Equivalent to 2.10^-2mol.L^-1 when nocells if right calculus
    double nutrientInitconcentration(pow(10,-2));//same unit
    double cellDiffusivity (0);
    
    
    if (_compartmentAmountLinear>1){
        cellDiffusivity=(pow(10,-10));
        Medium medium(nutrientDiffusivity,cellDiffusivity); //Diffusion coefficients of the medium; order: Nutrients, Public Goods, Cells
    }
    
    Medium medium(nutrientDiffusivity,cellDiffusivity);
    
    //2.Compartments definition//
    
    double heightEnvironmment(delta_l*pow(10,1));
    double lengthCompartment=lengthEnvironment/_compartmentAmountLinear;
    int widthNumberCompartment=widthEnvironment/lengthCompartment;
    int defaultLocation[2]{1,1};
    
    Compartment InsideCompartment(defaultLocation,lengthCompartment,heightEnvironmment,medium,1*nutrientInitconcentration,0*nutrientInitconcentration);//Concentrations in mol.L^(-1); order: Standard Nutrient, Hidden Nutrient, Revealed Nutrient, Common Goods
    InsideCompartment.PrintState();
    Compartment UpperOutsideCompartment(lengthCompartment,heightEnvironmment,medium,0,0);//Concentrations in mol.L^(-1) ; order: same as above
    
    
    //3.Alleles definition//
    double physiologicalMutationRate(_mutationRate);
    double regulatoryMutationRate(_mutationRate);
    
    //i)Physiological alleles//
    PhysiologicalAllele Nut1allele(sequence[0],physiologicalMutationRate,RandomSimulationNumber,"Nut1");          //sequence should result from pre-simulation, mutation rate as a variable parameter.
    PhysiologicalAllele Nut2allele(sequence[1],physiologicalMutationRate,RandomSimulationNumber,"Nut2");     //sequence should result from pre-simulation, mutation rate as a variable parameter.
    
    std::vector <RegulatoryAllele> positivetranscriptionfactors;               //vector of transcription enhancers
    std::vector <RegulatoryAllele> negativetranscriptionfactors;               //vector of transcription inhbitors
    
    //ii)Regulatory alleles//
    for (int TFnumber=0;TFnumber<int(TFamount/2);TFnumber++){
        RegulatoryAllele randompositiveTranscriptionFactor(sequence[TFnumber+traitAmount],regulatoryMutationRate,RandomSimulationNumber,bindingsite[TFnumber],'A');
        RegulatoryAllele randomnegativeTranscriptionFactor(sequence[TFnumber+traitAmount+int(TFamount/2)],regulatoryMutationRate,RandomSimulationNumber,bindingsite[TFnumber+int(TFamount/2)],'R');
        positivetranscriptionfactors.push_back(randompositiveTranscriptionFactor);
        negativetranscriptionfactors.push_back(randomnegativeTranscriptionFactor);
    }
    
    //iii)Genome assemblage and sizer mechanisms definition (artificial or in the RGN)//
    map <string,PhysiologicalAllele> physiologicalInitialAlleles;
    physiologicalInitialAlleles.insert(pair<string,PhysiologicalAllele>("Nut1",Nut1allele));               //allele encoding Nutrient 1 protein
    physiologicalInitialAlleles.insert(pair<string,PhysiologicalAllele>("Nut2",Nut2allele));     //allele encoding Nutrient 2 protein
    
    //4. Cell definition//
    
    double volumeMutationRate(_volumeMutationRate);
    
    double tauTranscriptionOncePolymeraseBound(0.1);
    double basalBindingProbabilityPolymerase(0.01);
    double enhancedBindingProbabilityPolymerase(0.5);
    double tauTraduction(0.1);
    double tauDegradationTranscript(5*pow(10,-4));
    double tauDegradationProtein(2.5*pow(10,-4));
    double initialCellVolume(30);
    
    Cell cell(_nutrientTypeAmount,
              K_on,                     //K_on BSTF
              K_off,                    //K_off BSTF
              kcat,                     //enzyme max activity
              KM,                       //substrate concentration for half max rate
              cellPermeability,         //permeabilityCell (unit: micrometer.s^-1, corresponds to estimated oxygen permeability)
              tauTranscriptionOncePolymeraseBound,                            //tauTranscriptionOncePolymeraseBound  (transcription rate once polymerase is bound: acknowledged as a constant, cf. Bio by the number)
              basalBindingProbabilityPolymerase,                           //basalBindingProbabilityPolymerase   (polymerase binding probability without the help of an enhancer: chosen such that the complete range of values from 0 to maximum transcription rate can occur)
              enhancedBindingProbabilityPolymerase,                            //enhancedBindingProbabilityPolymerase (polymerase binding probability wuth the help of an enhancer: chosen such that the maximum transcription rate described in litterature can occur)
              tauTraduction,                            //tauTraduction (translation rate according to scientific litterature, cf. Bio. by the numbers, for one ribosome: generally more than one can be actively translating a single mRNA, function of the energy as a result of the basal metabolism investment (which is only depending on the volume here, not subject to any evolutionary process))
              tauDegradationTranscript,                          //tauDegradationTranscript (mRNA degradation rate according to scientific litterature, cf. Bio. by the numbers)
              tauDegradationProtein,                         //tauDegradationProtein  (protein degradation rate according to scientific litterature, cf. Bio. by the numbers)
              //ATP expenses for the Surface Area used in calibration
              initialCellVolume,                             //cellVolume  (square micrometers)
              _aPropDesc,                            //volumeProportionDescendant (percentage: subject to sensitivity study because we intuit that it could influence differentiation by increasing noise
              volumeCellMitosis,                            //volumeCellMitosis  (unused parameter: served at the beginning to parameterize the models)
              volumeMutationRate,
              physiologicalInitialAlleles,
              positivetranscriptionfactors,
              negativetranscriptionfactors,
              RandomSimulationNumber
              );
    
    cell.EvaluateBSTFConnectivity();
    //cell.PrintCell(); Pb. of connectivity matrix
    std::cout << initialCellAmount << endl;
    std::cout << cell.printCellMatrix() << endl;
    //5. Puddle definition//
    
    int timeStepRatio(15); //Ratio between diffusion dynamics and genomic dynamics, example: 1s. for diffusion processes while 10s. for transcription
    
    DonutPuddle donutPuddle1(RandomSimulationNumber,       //seed initialization
                            _switchProbability,
                             Nut1Input,
                             Nut2Input,
                            m_delta_t,                      //timestep used in the simulation
                            lengthCompartment,            //length of the compartment (unit decimeters, cubic shape)
                            timeStepRatio,
                            _compartmentAmountLinear,     //length of the ecosystem (unit: number of compartment)
                            widthNumberCompartment,       //width of the ecosystem (unit: number of compartment)
                            InsideCompartment,            //puddle inside compartment
                            UpperOutsideCompartment,      //puddle outside compartment
                            initialCellAmount,            //initial number of cells
                            cell,
                            fileoutcompartments,
                            fileoutcells
                            );               //initial cells in the simulation
    //here: not passed as a reference: does it matter?
    m_donutPuddle=donutPuddle1;
    //m_donutPuddle.PrintEcosystem();
}

Simulation::Simulation(const Simulation &_simulation){
    
    m_delta_t=_simulation.m_delta_t;
    m_simulationTimeStep=_simulation.m_simulationTimeStep;
    m_simulationNumber=_simulation.m_simulationNumber;
    m_nutrientConcentrationRatio=_simulation.m_nutrientConcentrationRatio;
    m_donutPuddle=_simulation.m_donutPuddle;
    RandomSimulationNumber=_simulation.RandomSimulationNumber;
}

void Simulation::operator=(const Simulation &_simulation){
    
    m_delta_t=_simulation.m_delta_t;
    m_simulationTimeStep=_simulation.m_simulationTimeStep;
    m_simulationNumber=_simulation.m_simulationNumber;
    m_nutrientConcentrationRatio=_simulation.m_nutrientConcentrationRatio;
    m_donutPuddle=_simulation.m_donutPuddle;
    RandomSimulationNumber=_simulation.RandomSimulationNumber;
}

//##Destructor##// void as the random number gsl is freed in the simulation module
Simulation::~Simulation(){
    gsl_rng_free(RandomSimulationNumber);
}

void Simulation::SimulationProcessing(int _simulationDuration){
    
    /*Definition of the time allocator*/
    std::cout << _simulationDuration << endl;
    
    int m_simulationTimeStep(_simulationDuration/m_delta_t);
    int t_output=100;
    m_donutPuddle.OutputGetting(m_simulationNumber,0);
    
    for (int t=1;t<=m_simulationTimeStep;t++){
        
        /*In silico reality biological simulation*/
        m_donutPuddle.TimeStepInsideEcosystem(t,m_delta_t);//Simulation inside each compartment
        m_donutPuddle.DiffusionBetweenCompartments(); //Diffusion between compartments
        /*Periodic output getting to follow evolutionary trajectories*/
        if (t==t_output){
            m_donutPuddle.OutputGetting(m_simulationNumber,t_output);
            t_output+=m_simulationTimeStep/200;
        }
        
    }
    /*Final output getting: useless if already included in the regular tracking*/
    m_donutPuddle.OutputGetting(m_simulationNumber,_simulationDuration);
    
    //freeing memory of the random starting point associated with the simulation number
    //time(&end);
    
    
}
