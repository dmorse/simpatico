#ifndef MCMD_SYSTEM_TEST_H
#define MCMD_SYSTEM_TEST_H

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

#include <mcMd/simulation/Simulation.h>
#include <mcMd/simulation/System.h>
#include <simp/species/Species.h>
#include <mcMd/chemistry/Molecule.h>

using namespace Util;
using namespace McMd;

class SystemTest : public ParamFileTest
{

public:

   SystemTest()
   {}  

   virtual void setUp();

   void addMolecules();

   void dump();

   void testAddMolecules();

private:

   Simulation simulation_;

   System system_;

   /// Needed to create a stand-alone System object
   int    nAtomType_;

};



void SystemTest::setUp()
{

   #ifdef SIMP_ANGLE
   #ifdef SIMP_DIHEDRAL
   openFile("in/SimulationAngleDihedral"); 
   #else
   openFile("in/SimulationAngle"); 
   #endif
   #else
   #ifdef SIMP_DIHEDRAL
   openFile("in/SimulationDihedral"); 
   #else
   openFile("in/Simulation"); 
   #endif
   #endif

   simulation_.readParameters(file());
   system_.setId(0);
   nAtomType_ = 2;
   system_.setSimulation(simulation_);
   system_.allocateMoleculeSets();
}

void SystemTest::addMolecules()
{
   Species *speciesPtr;
   int      i, j;
   for (i = 0; i < simulation_.nSpecies() ; ++i ) {
      speciesPtr = &simulation_.species(i);
      for (j = 0; j < speciesPtr->capacity(); ++j) {
         system_.addMolecule(simulation_.getMolecule(i));
      }
   }
}

void SystemTest::dump()
{
   const Species             *speciesPtr;
   //const System::MoleculeSet *moleculeSetPtr;
   const Molecule            *moleculePtr;
   const Atom                *atomPtr;
   int   iSpecies, iMolecule, nMolecule, iAtom, nAtom;

   std::cout << std::endl; 
   for (iSpecies=0; iSpecies < simulation_.nSpecies() ; ++iSpecies ) {
      speciesPtr     = &simulation_.species(iSpecies);
      nAtom          = speciesPtr->nAtom(); 
      //moleculeSetPtr = &(system_.moleculeSet(iSpecies));
      //nMolecule      = moleculeSetPtr->size(); 
      nMolecule      = system_.nMolecule(iSpecies);

      std::cout << "Species Id     = " << speciesPtr->id() << std::endl; 
      std::cout << "# in System    = " << nMolecule << std::endl;
      for (iMolecule = 0; iMolecule < nMolecule; ++iMolecule) {
         moleculePtr = &system_.molecule(iSpecies,iMolecule);
         std::cout << "Molecule Id    = " << moleculePtr->id() << std::endl;
         for (iAtom = 0; iAtom < nAtom; ++iAtom) {
            atomPtr = &moleculePtr->atom(iAtom);
            std::cout << "   Atom # " << atomPtr->id()
                      << "   type "   << atomPtr->typeId() << std::endl;
         }
      }
   }
}

void SystemTest::testAddMolecules()
{
   printMethod(TEST_FUNC);
   using std::cout;

   addMolecules();
  
   Species *speciesPtr;
   int      iSpecies;
   for (iSpecies=0; iSpecies < simulation_.nSpecies() ; ++iSpecies ) {
      speciesPtr = &simulation_.species(iSpecies);
      try {
         speciesPtr->isValid();
      } catch (Exception e) {
         e.write(std::cout);
         TEST_ASSERT(0);
      }
      TEST_ASSERT(system_.nMolecule(iSpecies) == speciesPtr->capacity());
   }

   if (verbose() > 0) dump();

}

#if 0
void SystemTest::testReadConfig()
{
   printMethod(TEST_FUNC);
   using std::cout;

   std::ifstream  configFile("in/config");
   system_.readConfig(configFile);

   if (verbose() > 0) {
      std::cout << std::endl; 
      system_.writeConfig(std::cout);
   }

}


void SystemTest::testInitMoleculeIterator()
{
   printMethod(TEST_FUNC);
   using std::cout;

   // Add all molecules to the System
   addMolecules();

   System::MoleculeIterator  molIter;
   Molecule::AtomIterator    atomIter;
   int    iMolecule;

   // Count molecules of each Species
   for (int i=0; i < simulation_.nSpecies() ; ++i ) {
      iMolecule = 0; 
      for (system_.begin(i, molIter); molIter.notEnd(); ++molIter) {
         ++iMolecule;
      }
      TEST_ASSERT(iMolecule == system_.nMolecule(i));
   }

   // Output atomic positions
   if (verbose() > 1) {
      std::cout << std::endl;
      for (int i=0; i < simulation_.nSpecies() ; ++i ) {
         std::cout << "species Id    = " << i << std::endl;
         for (system_.begin(i, molIter); molIter.notEnd(); ++molIter) {
            std::cout << "Molecule Id    = " << molIter->id() << std::endl;
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               std::cout << "   " << atomIter->id();
               std::cout << "   " << atomIter->typeId() << std::endl;
            }
         }
      }
   }
}
#endif

TEST_BEGIN(SystemTest)
TEST_ADD(SystemTest, testAddMolecules)
TEST_END(SystemTest)

#endif
