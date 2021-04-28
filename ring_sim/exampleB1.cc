#include "B1RegularDetectorConstruction.hh"
#include "B1ActionInitialization.hh"
#include "params.hh"
#include "globalFunctions.hh"
#include "G4SystemOfUnits.hh"

//TODO: how to control G4MULTITHREADED in a generic way???
#include "G4MTRunManager.hh"
#include "G4RunManager.hh"

#include "G4GenericBiasingPhysics.hh"
#include "G4UImanager.hh"
#include "G4ScoringManager.hh"
#include "QBBC.hh"
#include "FTFP_BERT.hh"
#include "G4PhysListFactory.hh"
#include "G4EmPenelopePhysics.hh"
#include "PhysListEmStandard.hh"
#include "G4HumanPhantomPhysicsList.hh"
#include "B1ExtraPhysics.hh"
#include "B1EnergyDeposit.hh"


#include "G4Run.hh"

#include "G4NistManager.hh"


#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "Randomize.hh"
#include <iostream>
#include <vector>
#include "Shielding.hh"

//calculate MFP
#include "G4EmCalculator.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4ComptonScattering.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4PEEffectFluoModel.hh"
#include "PhysicsList.hh"
#include "G4StepLimiterPhysics.hh"
#include "B1ModularPhysicsList.hh"

#include <ctime>
#include <iostream>
#include <boost/python.hpp>
//#include <boost/python/dict.hpp>
#include <boost/python/numpy.hpp>
#include "params_dict.hh"
#include "ParamsSingleton.hh"

namespace p = boost::python;
namespace np = boost::python::numpy;

char const* greet()
{
  return "hello, world"; 
}

void print_dict(boost::python::dict x){

}

void printXSPerAtom(G4Element* el); //implemented at the end
void printXSPerAtom_25(G4Element* el, G4double* spectrum); //implemented at the end

//global array
//TODO: seems like this is not in use - consider removing!
//B1EnergyDeposit* detectorsArray[NUM_OF_THREADS];


p::tuple run_main(boost::python::dict x)
{

	p::list elements_py = create_elements_list(x);
	p::list A_py = create_A_list(x);
	p::list Z_py = create_Z_list(x);
	np::ndarray material_py = create_material_array(x);
	np::ndarray ids_py = create_ids_array(x);
	np::ndarray dens_py = create_dens_array(x);
	np::ndarray src_py = create_src_array(x);
	np::ndarray src_orient_py = create_src_orient_array(x);
	np::ndarray det_py = create_det_array(x);
	np::ndarray det_orient_py = create_det_orient_array(x);
	np::ndarray spectrum_py = create_spectrum_array(x);
	np::ndarray projection_code_py = create_projection_code_array(x);
	np::ndarray intens_code_py = create_intens_code_array(x);
	np::ndarray phantome_offset_orient_py = create_phantom_loc_orient_array(x);
	np::ndarray error_py = create_error_array(x);

	ParamsSingleton* pp = ParamsSingleton::Instance();
	pp->parse_dict(x);
	pp->extract_elements_list(elements_py);
	pp->extract_A_list(A_py);
	pp->extract_Z_list(Z_py);
	pp->extract_materials_array(material_py);
	pp->extract_ids_array(ids_py);
	pp->extract_dens_array(dens_py);
	pp->extract_src_array(src_py);
	pp->extract_src_orient_array(src_orient_py);
	pp->extract_det_array(det_py);
	pp->extract_det_orient_array(det_orient_py);
	pp->extract_spectrum_array(spectrum_py);
	pp->extract_projection_code_array(projection_code_py);
	pp->extract_intens_code_array(intens_code_py);
	pp->extract_phantom_offset_orientation_array(phantome_offset_orient_py);
	pp->extract_error_array(error_py);
	pp->create_run_images_arr();
	if (pp->int_map["i_CALC_GRADIENT"] == 0){
		pp->create_images_arr();
		
	}else{
		pp->create_P_arr();
		pp->create_grad_arr();
	}
  // Detect interactive mode (if no arguments) and define UI session
	G4UIExecutive* ui = 0;

	// Choose the Random engine
  	CLHEP::HepJamesRandom* eng = new CLHEP::HepJamesRandom;
      std::time_t result = std::time(nullptr);
      eng->setSeed(result);
  	G4Random::setTheEngine(eng);
//  	G4Random::setTheEngine(new CLHEP::RanecuEngine);

//	// Choose the Random engine
//	CLHEP::RanecuEngine* eng = new CLHEP::RanecuEngine;
//	std::time_t result = std::time(nullptr);
//	eng->setSeed(result);
//	//eng->setIndex(8);
//	//eng->setSeed(8);
//	G4Random::setTheEngine(eng);
	// Construct the default run manager
	G4RunManager* runManager = NULL;
	if (pp->int_map["i_MULTI"] == 1)
	{
		  G4MTRunManager* runManagerMT = new G4MTRunManager;
		  runManagerMT->SetNumberOfThreads(pp->int_map["i_NUM_OF_THREADS"]);
		  std::cout << "numOfThreads: " << runManagerMT-> GetNumberOfThreads() << std::endl;
		  std::cout << "multi!!!! " << std::endl;
		  runManager = (G4RunManager*)runManagerMT;

	}  else {
		  runManager = new G4RunManager;
		  std::cout << "Single-Threaded Mode" <<  std::endl;
	}


	// Activate UI-command base scorer
	G4ScoringManager * scManager = G4ScoringManager::GetScoringManager();
	scManager->SetVerboseLevel(pp->int_map["i_VERBOSE_SCORING"]);

	// Set mandatory initialization classes
	// Detector construction
	B1DetectorConstruction* theGeometry = 0;
	theGeometry = new B1RegularDetectorConstruction();
	runManager->SetUserInitialization(theGeometry);

	// Physics list
	/*
	G4VModularPhysicsList* physicsList;

	if (pp->double_map["f_PHANTOM_PRODUCTION_CUTS"] == 1){
		physicsList = new B1ModularPhysicsList("myList");
	}
	else{
		G4PhysListFactory factory;
		if (pp->int_map["i_LIV_MODE"] == 1){
			physicsList = factory.GetReferencePhysList("FTFP_BERT_LIV");
			// G4StringemName = G4String("local");
			// G4VPhysicsConstructor* emPhysicsList = new PhysListEmStandard(emName);
  
		} else{
			physicsList = factory.GetReferencePhysList("FTFP_BERT_EMV");
		}
			//FTFP_BERT_PEN
			//G4VModularPhysicsList* physicsList = new QBBC;
	}
	//physicsList->SetCutsForRegion(10*km,"waterRegion");

	//biasing
	G4GenericBiasingPhysics* biasingPhysics = new G4GenericBiasingPhysics();
	  if ( pp->int_map["i_BIASING"] == 1 )
		{
		  //TODO:fix
		  //biasingPhysics->Bias("gamma");
		  // -- Create list of physics processes to be biased: only brem. in this case:
		  //TODO: is this necessary??
		  std::vector< G4String > processToBias;
		  processToBias.push_back("compt");
		  processToBias.push_back("Rayl");
		  // -- Pass the list to the G4GenericBiasingPhysics, which will wrap the eBrem
		  // -- process of e- and e+ to activate the biasing of it:
		  biasingPhysics->PhysicsBias("gamma", processToBias);

		  physicsList->RegisterPhysics(biasingPhysics);
		  G4cout << "      ********************************************************* " << G4endl;
		  G4cout << "      ********** processes are wrapped for biasing ************ " << G4endl;
		  G4cout << "      ********************************************************* " << G4endl;
		}
	  else
		{
		  G4cout << "      ************************************************* " << G4endl;
		  G4cout << "      ********** processes are not wrapped ************ " << G4endl;
		  G4cout << "      ************************************************* " << G4endl;
		}

	physicsList->SetVerboseLevel(pp->int_map["i_VERBOSE_PHYSICS_LIST"]);

	if (pp->double_map["f_DETECTOR_SPECIAL_CUTS"] == 1) {
		physicsList->RegisterPhysics(new B1ExtraPhysics());
	}
	
	// runManager->SetUserInitialization(physicsList);
	*/
	PhysicsList* phys = new PhysicsList;
	// phys->AddPhysicsList("local");
	 runManager->SetUserInitialization(phys);

	// User action initialization
	runManager->SetUserInitialization(new B1ActionInitialization(theGeometry));


		// Initialize visualization
	G4VisManager* visManager = new G4VisExecutive;
	// G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
	// G4VisManager* visManager = new G4VisExecutive("Quiet");
	visManager->SetVerboseLevel(pp->int_map["i_VERBOSE_VIS"]);
	visManager->Initialize();

	// Get the pointer to the User Interface manager
	G4UImanager* UImanager = G4UImanager::GetUIpointer();
	 
	 // UImanager->ApplyCommand("/gate/physics/processes/PhotoElectric/setModel StandardModel");
	 UImanager->ApplyCommand("/run/verbose " + IntToString(pp->int_map["i_VERBOSE_RUN"]));
	 UImanager->ApplyCommand("/event/verbose " + IntToString(pp->int_map["i_VERBOSE_EVENT"]));
	 UImanager->ApplyCommand("/tracking/verbose " + IntToString(pp->int_map["i_VERBOSE_TRACK"]));


	// Process macro or start UI session
	if ( ! ui ) {
		// batch mode
		G4String command = "/control/execute ";
//		G4String fileName = argv[1];
		G4String fileName = "ring_sim_run_simple.mac";
		UImanager->ApplyCommand(command+fileName);
	}
	else {
		// interactive mode
		UImanager->ApplyCommand("/control/execute init_vis.mac");

		ui->SessionStart();
		delete ui;
	}

	//******************************************************************
	//print elements cross section
	if (pp->int_map["i_PRINT_ELEMENTS_XS"] == 1){
		G4cout << "      ********************************************************* " << G4endl;
		  G4cout << "      **********  printing cross sections ************ " << G4endl;
		  G4cout << "      ********************************************************* " << G4endl;
		printXSPerAtom_25(new G4Element("Hydrogen","H",1.0,1.008 * g/mole), pp->spectrum_arr);
		// printXSPerAtom_25(new G4Element("Helium","He",2.0,4.0026 * g/mole ));
		// printXSPerAtom_25(new G4Element( "Lithium","Li",3.0, 6.941  * g/mole ));
		// printXSPerAtom_25(new G4Element("Beryllium","Be",4.0, 9.012182  * g/mole ));
		// printXSPerAtom_25(new G4Element("Boron","B",5.0, 10.811  * g/mole ));
		printXSPerAtom_25(new G4Element( "Carbon","C",6.0, 12.011 * g/mole ), pp->spectrum_arr);
		printXSPerAtom_25(new G4Element( "Nitrogen","N",7.0, 14.007 * g/mole ), pp->spectrum_arr);
		printXSPerAtom_25(new G4Element( "Oxygen","O",8.0, 16.00  * g/mole ), pp->spectrum_arr);
		// printXSPerAtom_25(new G4Element( "Fluorine","F",9.0, 18.998404  * g/mole ));
		// printXSPerAtom_25(new G4Element(  "Neon","Ne",10.0, 20.1797  * g/mole ));
		// printXSPerAtom_25(new G4Element( "Sodium","Na",11.0, 22.98977 * g/mole ));
		// printXSPerAtom_25(new G4Element( "Magnesium","Mg",12.0, 24.305 * g/mole ));
		// printXSPerAtom_25(new G4Element( "Aluminum","Al",13.0, 26.981539 * g/mole ));
		printXSPerAtom_25(new G4Element( "Phosphorus","P",15.0, 30.97376 * g/mole ), pp->spectrum_arr);
		// printXSPerAtom_25(new G4Element( "Sulfur","S",16.0,32.065* g/mole ));
		// printXSPerAtom_25(new G4Element( "Chlorine","Cl",17.0, 35.453* g/mole ));
		// printXSPerAtom_25(new G4Element( "Argon","Ar",18.0, 39.948 * g/mole ));
		printXSPerAtom_25(new G4Element( "Potassium","K",19.0, 39.0983* g/mole ), pp->spectrum_arr);
		printXSPerAtom_25(new G4Element("Calcium","Ca",20.0, 40.078* g/mole ), pp->spectrum_arr);
		// printXSPerAtom_25(new G4Element( "Scandium","Sc",21.0, 44.95591 * g/mole ));
		// printXSPerAtom_25(new G4Element( "Titanium","Ti",22.0, 47.867 * g/mole ));
		// printXSPerAtom_25(new G4Element( "Vanadium","V",23.0, 50.9415 * g/mole ));
		// printXSPerAtom_25(new G4Element( "Chromium","Cr",24.0, 51.9961 * g/mole ));
		// printXSPerAtom_25(new G4Element( "Manganese","Mn",25.0, 54.93805 * g/mole ));
		// printXSPerAtom_25(new G4Element( "Iron","Fe",26, 55.845* g/mole ));
		// printXSPerAtom_25(new G4Element( "Iodine","I",53, 126.90447 * g/mole ));
		// printXSPerAtom_25(new G4Element( "Lead","Pb",82, 207.2 * g/mole ));
	}

	//**********************************************************************

	if (pp->int_map["i_WRITE_TO_FILES"] == 1){
		pp->free_run_images_mem();
		if (pp->int_map["i_CALC_GRADIENT"] == 0){
			pp->free_images_mem();
			
			
		}else{
			pp->free_grad_mem();
			pp->free_P_mem();
		}
	
		delete visManager;
		delete runManager;

		return p::make_tuple(0);
	} else {
		if (pp->int_map["i_CALC_GRADIENT"] == 1){
			G4int N_shots = pp->int_map["i_NUM_SHOTS"];
			G4int N_voxels = pp->int_map["i_NUM_OF_VOXELS"];
			G4int N_elements = pp->int_map["i_NUM_OF_ELEMENTS"];

			p::tuple shape_grad = p::make_tuple(N_shots, N_voxels, N_elements);
			np::dtype dtype_grad = np::dtype::get_builtin<G4double>();
			np::ndarray grad_result_arr = np::zeros(shape_grad, dtype_grad);

			p::tuple shape_P = p::make_tuple(N_shots, N_voxels);
			np::dtype dtype_P = np::dtype::get_builtin<G4double>();
			np::ndarray P_result_arr = np::zeros(shape_P, dtype_P);

			for (G4int i=0; i<pp->int_map["i_NUM_SHOTS"]; i++){
				for (G4int j=0; j<pp->int_map["i_NUM_OF_VOXELS"]; j++){
					P_result_arr[i][j] = pp->P_arr[i][j];
					for (G4int k=0; k<pp->int_map["i_NUM_OF_ELEMENTS"]; k++){
						grad_result_arr[i][j][k] = pp->grad_arr[i][j][k];
					}
				}
			}



			pp->free_grad_mem();
			pp->free_P_mem();

			delete visManager;
			delete runManager;

			return p::make_tuple(grad_result_arr, P_result_arr);

		} else {
			G4int N_shots = pp->int_map["i_NUM_SHOTS"];
			G4int N_scorers = pp->int_map["i_NUM_OF_SCORERS"];
			G4int N_cols = pp->int_map["i_NUM_OF_DET_COLS"];
			G4int N_rows = pp->int_map["i_NUM_OF_DET_ROWS"];
			p::tuple shape = p::make_tuple(N_shots, N_scorers, N_rows, N_cols);
			np::dtype dtype = np::dtype::get_builtin<double>();
			np::ndarray images_result_arr = np::zeros(shape, dtype);

			for (G4int i=0; i< N_shots; i++){
				for (G4int j=0; j< N_scorers; j++){
					for (G4int k=0; k< N_rows; k++){
						for (G4int l=0; l< N_cols; l++){
							images_result_arr[i][j][k][l] = pp->images_arr[i][j][k][l];
						} //for l
					} // for k
				} // for j
			} // for i
			pp->free_images_mem();



			delete visManager;
			delete runManager;

			return p::make_tuple(images_result_arr);
		} // if pp->int_map["i_CALC_GRADIENT"] == 1
	}  // if pp->int_map["i_WRITE_TO_FILES"] == 1


}


int main(int argc,char** argv)
{
	namespace p = boost::python;
	namespace np = boost::python::numpy;


	Py_Initialize();
	np::initialize();
	std::string extracted_val = "../run_inputs/dict.txt";

	p::dict di = read_params_dict(extracted_val);

	// p::list elements_py = create_elements_list(di);
	// p::list A_py = create_A_list(di);
	// p::list Z_py = create_Z_list(di);
	// np::ndarray material_py = create_material_array(di);
	// np::ndarray id_array = create_ids_array(di);
	// np::ndarray dens_array = create_dens_array(di);
	// np::ndarray srcs_array = create_src_array(di);
	// np::ndarray srcs_orient_array = create_src_orient_array(di);
	// np::ndarray dets_array = create_det_array(di);
	// np::ndarray dets_orient_array = create_det_orient_array(di);
	// np::ndarray spectrum_array = create_spectrum_array(di);
	// np::ndarray projection_code_array = create_projection_code_array(di);
	// np::ndarray intens_code_array = create_intens_code_array(di);
	// np::ndarray phantom_loc_orient_array = create_phantom_loc_orient_array(di);
	// np::ndarray error_array = create_error_array(di);


	//TODO: add additional arguments to run_main
	run_main(di);

	return 1;
}


void printXSPerAtom(G4Element* el){
	ParamsSingleton* pp = ParamsSingleton::Instance();

	G4EmCalculator emCalculator;
	std::string elementName = el->GetName();
	G4int elementZ = el->GetZ();
	std::ofstream outputElementFile;

	std::string fileName =  "../run_outputs/ElementsXS/" + IntToString(elementZ) + ".csv";
	outputElementFile.open(fileName.c_str());

	G4double highE = pp->double_map["f_PARTICLE_ENERGY"]*keV + 3*keV; //3 is buffer
	G4double DeltaEnergy = 0.1*keV; //delta Energy
	G4double Energy = DeltaEnergy;

	G4double comptXS;
	G4double photXS;

	// print phot
	while (Energy<highE) {
		photXS = emCalculator.ComputeCrossSectionPerAtom(Energy,"gamma","phot",el,0);
		outputElementFile << photXS/cm2 << ",";
		Energy = Energy + DeltaEnergy;
	}
	outputElementFile << "\n";
	// print compt
	Energy = DeltaEnergy; //delta Energy
	while (Energy<highE) {
		comptXS = emCalculator.ComputeCrossSectionPerAtom(Energy,"gamma","compt",el,0);
		outputElementFile << comptXS/cm2 << ",";
		Energy = Energy + DeltaEnergy;
	}
	outputElementFile << "\n";
	outputElementFile.close();

}

void printXSPerAtom_25(G4Element* el, G4double* spectrum){
	G4EmCalculator emCalculator;
	std::string elementName = el->GetSymbol();
	G4int elementZ = el->GetZ();
	G4double elementA = el->GetAtomicMassAmu();
	
	std::vector<G4double> cross_sections_absorb;
	std::vector<G4double> cross_sections_compt;
	std::vector<G4double> cross_sections_rayl;
	G4double u = 1.660539e-24;
	
	std::vector<G4double> energies;
	
  
    for (int i = 1; i <= 150; i++)
        energies.push_back(i);
	
	G4double comptXS;
	G4double RaylXS;
	G4double photXS;
	
	for (G4double energy : energies) {
		photXS = emCalculator.ComputeCrossSectionPerAtom(energy*keV, "gamma","phot",el,0);
		comptXS = emCalculator.ComputeCrossSectionPerAtom(energy*keV, "gamma","compt",el,0);
		RaylXS = emCalculator.ComputeCrossSectionPerAtom(energy*keV, "gamma","Rayl",el,0);
		// G4double cross_section = (photXS + comptXS + RaylXS)/(u*elementA);
		// G4cout << "ELEMENT: " << elementName << " ENERGY: " << energy << " ABSORB: " << photXS/cm2 << " COMPTON: " << comptXS/cm2 << " Rayl: " << RaylXS/cm2 << " A: " << elementA << " u: " << u << G4endl;
		cross_sections_absorb.push_back(photXS);
		cross_sections_compt.push_back(comptXS);
		cross_sections_rayl.push_back(RaylXS);
	}
		
	
	
	std::ofstream outputElementFile;
	
	std::string fileName = "/home/yonatangat/Projects/scattering_tomo/elmental_media_info/" + elementName + ".csv";
	outputElementFile.open(fileName.c_str());

	// print phot

	for( G4int i = 0; i<energies.size(); i++)
		outputElementFile << energies[i] << ',' << spectrum[i] << ','  << cross_sections_absorb[i]/cm2 << ',' << cross_sections_compt[i]/cm2 << ',' << cross_sections_rayl[i]/cm2 << "\n";

	// outputElementFile << "\n";
	outputElementFile.close();

}
BOOST_PYTHON_MODULE(exampleB1_lib)   {

  using namespace boost::python;
  Py_Initialize();
  np::initialize();
  def("greet", greet);
  def("run_main", run_main);
  def("print_dict", print_dict);
}
