/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McCommandManager.h"                
#include "McSimulation.h"  
#include "mc_potentials.h"  
#include <mcMd/commands/McCommandFactory.h> 

#include <util/format/Str.h> 
#include <util/format/Dbl.h> 

namespace McMd
{

   using namespace Util;

   // Constructor.
   McCommandManager::McCommandManager(McSimulation& simulation)
    : CommandManager(),
      simulationPtr_(&simulation),
      systemPtr_(&simulation.system())
   {
      // Note: No setClassname("MdCommandManager");
      // Retain name CommandManager set by base class.
   }  

   // Constructor.
   McCommandManager::McCommandManager(McSimulation& simulation, 
                                      McSystem& system)
    : CommandManager(),
      simulationPtr_(&simulation),
      systemPtr_(&system)
   {}

   // Destructor.
   McCommandManager::~McCommandManager()
   {}

   /// Return pointer to a new CommandFactory.
   Factory<Command>* McCommandManager::newDefaultFactory() const
   {  return new McCommandFactory(*simulationPtr_, *systemPtr_); }

   // Attempt to read one standard (built-in) command.
   bool 
   McCommandManager::readStandardCommand(std::string name, std::istream& in)
   {
      std::string   filename;
      std::ifstream inputFile;
      std::ofstream outputFile;
      bool success = true;

      if (name == "SET_CONFIG_IO") {
         std::string classname;
         in >> classname;
         Log::file() << Str(classname, 15) << std::endl;
         system().setConfigIo(classname);
      } else
      if (name == "READ_CONFIG") {
         in >> filename;
         Log::file() << Str(filename, 15) << std::endl;
         simulation().fileMaster().openInputFile(filename, inputFile);
         system().readConfig(inputFile);
         inputFile.close();
      } else
      if (name == "SIMULATE") {
         int endStep;
         in >> endStep;
         Log::file() << "  " << endStep << std::endl;
         bool isContinuation = false;
         simulation().simulate(endStep, isContinuation);
      } else
      if (name == "CONTINUE") {
         if (simulation().iStep() == 0) {
            UTIL_THROW("Attempt to continue when iStep() == 0");
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
         Log::file() << "  " <<  min << "  " <<  max
                     << "  " <<  filename << std::endl;
         simulation().analyzeConfigs(min, max, filename);
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
      if (name == "WRITE_CONFIG") {
         in >> filename;
         Log::file() << Str(filename, 15) << std::endl;
         simulation().fileMaster().openOutputFile(filename, outputFile);
         system().writeConfig(outputFile);
         outputFile.close();
      } else
      if (name == "WRITE_PARAM") {
         in >> filename;
         Log::file() << "  " << filename << std::endl;
         simulation().fileMaster().openOutputFile(filename, outputFile);
         writeParam(outputFile);
         outputFile.close();
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
         Log::file() << "\n  Boundary:    " << system().boundary();
         Label capacityLabel("Capacities:");
         in >> capacityLabel;
         Log::file() << "\n  Capacities: ";
         for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
            in >> capacities[iSpecies];
            Log::file() << "  " << capacities[iSpecies];
         }
         Label diameterLabel("Diameters:");
         in >> diameterLabel;
         Log::file() << "\n  Diameters: ";
         for (int iType=0; iType < nAtomType; iType++) {
            in >> diameters[iType];
            Log::file() << "  " << diameters[iType];
         }
         Log::file() << std::endl;

         system().generateMolecules(capacities, diameters);

      } else
      if (name == "DEFORM_CELL") {

         // Read in configuration from file
         in >> filename;
         Log::file() << Str(filename, 15) << std::endl;
         simulation().fileMaster().openInputFile(filename, inputFile);
         system().readConfig(inputFile);
         inputFile.close();

         int nSpecies = simulation().nSpecies();
         System::MoleculeIterator molIter;
         Molecule::AtomIterator atomIter;
         for (int iSpec=0; iSpec < nSpecies; ++iSpec) {
            system().begin(iSpec, molIter);
            for ( ; molIter.notEnd(); ++molIter) {
               molIter->begin(atomIter);
               for ( ; atomIter.notEnd(); ++atomIter) {
                  Vector cartPosition, genPosition;
                  cartPosition = atomIter->position();
                  system().boundary().transformCartToGen(cartPosition, genPosition);
                  atomIter->position() = genPosition;
               }
            }
         }

         // Read in new boundary
         in >> system().boundary();
         Log::file() << "  " << system().boundary();
         Log::file() << std::endl;

         for (int iSpec=0; iSpec < nSpecies; ++iSpec) {
            system().begin(iSpec, molIter);
            for ( ; molIter.notEnd(); ++molIter) {
               molIter->begin(atomIter);
               for ( ; atomIter.notEnd(); ++atomIter) {
                  Vector cartPosition, genPosition;
                  genPosition = atomIter->position();
                  system().boundary().transformGenToCart(genPosition, cartPosition);
                  atomIter->position() = cartPosition;
               }
            }
         }

         // Write out configuration to file
         in >> filename;
         Log::file() << Str(filename, 15) << std::endl;
         simulation().fileMaster().openOutputFile(filename, outputFile);
         system().writeConfig(outputFile);
         outputFile.close();

         #ifndef SIMP_NOPAIR 
         // Generate cell list
         system().pairPotential().buildCellList();
         #endif

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
      #endif // ifdef SIMP_DIHEDRAL
      #endif // ifndef UTIL_MPI
      {
         // Command name not recognized
         success = false;
      }
      return success;
   }

}
