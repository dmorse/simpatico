/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MdCommandManager.h"                
#include "MdSimulation.h"  
#include "md_potentials.h"  
#include <mcMd/commands/MdCommandFactory.h> 

#include <util/format/Str.h> 
#include <util/format/Dbl.h> 

namespace McMd
{

   using namespace Util;

   // Constructor.
   MdCommandManager::MdCommandManager(MdSimulation& simulation)
    : CommandManager(),
      simulationPtr_(&simulation),
      systemPtr_(&simulation.system())
   {
      // Note: No command setClassName("MdCommandManager")
      // Retains name CommandManager set by base class.
   }

   // Constructor.
   MdCommandManager::MdCommandManager(MdSimulation& simulation, 
		                            MdSystem& system)
    : CommandManager(),
      simulationPtr_(&simulation),
      systemPtr_(&system)
   {}

   // Destructor.
   MdCommandManager::~MdCommandManager()
   {}

   // Return pointer to a new CommandFactory.
   Factory<Command>* MdCommandManager::newDefaultFactory() const
   {  return new MdCommandFactory(simulation(), system()); }

   // Attempt to read one standard (built-in) command.
   bool 
   MdCommandManager::readStandardCommand(std::string name, std::istream& in)
   {
      std::string filename;
      std::ifstream inputFile;
      std::ofstream outputFile;
      bool success = true;

      if (name == "READ_CONFIG") {
         in >> filename;
         Log::file() << Str(filename, 15) << std::endl;
         simulation().fileMaster().openInputFile(filename, inputFile);
         // Log::file() << "Opened config file" << std::endl;
         system().readConfig(inputFile);
         // Log::file() << "Finished reading config file" << std::endl;
         inputFile.close();
      } else
      if (name == "THERMALIZE") {
         double temperature;
         in >> temperature;
         Log::file() << Dbl(temperature, 15, 6) << std::endl;
         system().setBoltzmannVelocities(temperature);
         system().removeDriftVelocity();
      } else
      if (name == "SIMULATE") {
         int endStep;
         in >> endStep;
         Log::file() << Int(endStep, 15) << std::endl;
         bool isContinuation = false;
         simulation().simulate(endStep, isContinuation);
      } else
      if (name == "CONTINUE") {
         if (simulation().iStep() == 0) {
            UTIL_THROW("Attempt to continue simulation when iStep_ == 0");
         }
         int endStep;
         in >> endStep;
         Log::file() << Int(endStep, 15) << std::endl;
         bool isContinuation = true;
         simulation().simulate(endStep, isContinuation);
      } else
      if (name == "ANALYZE_CONFIGS") {
         int min, max;
         in >> min >> max >> filename;
         Log::file() <<  Int(min, 15) << Int(max, 15)
                     <<  Str(filename, 20) << std::endl;
         simulation().analyzeConfigs(min, max, filename);
      } else
      if (name == "WRITE_CONFIG") {
         in >> filename;
         Log::file() << Str(filename, 15) << std::endl;
         simulation().fileMaster().openOutputFile(filename, outputFile);
         system().writeConfig(outputFile);
         outputFile.close();
      } else
      if (name == "WRITE_PARAM") {
         in >> filename;
         Log::file() << Str(filename, 15) << std::endl;
         simulation().fileMaster().openOutputFile(filename, outputFile);
         simulation().writeParam(outputFile);
         outputFile.close();
      } else
      if (name == "SET_CONFIG_IO") {
         std::string classname;
         in >> classname;
         Log::file() << Str(classname, 15) << std::endl;
         system().setConfigIo(classname);
      } else
      if (name == "ANALYZE_TRAJECTORY") {
         std::string classname;
         std::string filename;
         int min, max;
         in >> min >> max >> classname >> filename;
         Log::file() << " " << Str(classname,15) 
                     << " " << Str(filename, 15)
                     << std::endl;
         simulation().analyzeTrajectory(min, max, classname, filename);
      } else 
      if (name == "GENERATE_MOLECULES") {
         DArray<double> diameters;
         DArray<int> capacities;
         int nAtomType = simulation().nAtomType();
         int nSpecies = simulation().nSpecies();
         diameters.allocate(nAtomType);
         capacities.allocate(nSpecies);

         // Parse name
         in >> system().boundary();
         Log::file() << "  " << system().boundary();
         Label capacityLabel("Capacities:");
         in >> capacityLabel;
         for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
            in >> capacities[iSpecies];
            Log::file() << "  " << capacities[iSpecies];
         }
         Label diameterLabel("Diameters:");
         in >> diameterLabel;
         for (int iType=0; iType < nAtomType; iType++) {
            in >> diameters[iType];
            Log::file() << "  " << diameters[iType];
         }
         Log::file() << std::endl;

         system().generateMolecules(capacities, diameters);
      } else
      #ifndef UTIL_MPI
      #ifndef SIMP_NOPAIR
      if (name == "SET_PAIR") {
         std::string paramName;
         int typeId1, typeId2; 
         double value;
         in >> paramName >> typeId1 >> typeId2 >> value;
         Log::file() << "  " <<  paramName 
                     << "  " <<  typeId1 << "  " <<  typeId2
                     << "  " <<  value << std::endl;
         system().pairPotential()
                 .set(paramName, typeId1, typeId2, value);
      } else
      #endif 
      #ifdef SIMP_BOND
      if (name == "SET_BOND") {
         std::string paramName;
         int typeId; 
         double value;
         in >> paramName >> typeId >> value;
         Log::file() << "  " <<  paramName << "  " <<  typeId 
                     << "  " <<  value << std::endl;
         system().bondPotential().set(paramName, typeId, value);
      } else
      #endif
      #ifdef SIMP_ANGLE
      if (name == "SET_ANGLE") {
         std::string paramName;
         int typeId; 
         double value;
         in >> paramName >> typeId >> value;
         Log::file() << "  " <<  paramName << "  " <<  typeId 
                     << "  " <<  value << std::endl;
         system().anglePotential().set(paramName, typeId, value);
      } else
      #endif 
      #ifdef SIMP_DIHEDRAL
      if (name == "SET_DIHEDRAL") {
         std::string paramName;
         int typeId; 
         double value;
         in >> paramName >> typeId >> value;
         Log::file() << "  " <<  paramName << "  " <<  typeId 
                     << "  " <<  value << std::endl;
         system().dihedralPotential().set(paramName, typeId, value);
      } else
      #endif 
      #endif // ifndef UTIL_MPI
      {
         // Failed to match command name
         success = false;
      }
      return success;
   }

}
