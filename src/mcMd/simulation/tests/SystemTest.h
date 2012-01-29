#ifndef SYSTEM_TEST_H
#define SYSTEM_TEST_H

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

#include <mcMd/simulation/Simulation.h>
#include <mcMd/simulation/System.h>
#include <mcMd/species/Species.h>
#include <mcMd/chemistry/Molecule.h>

using namespace Util;
using namespace McMd;

class SystemTest : public ParamFileTest<Simulation>
{

public:

   SystemTest()
   {}  

   virtual void setUp();

   void addMolecules();

   void dump();

   void testAddMolecules();

private:

   // Accessor Simulation& object() inherited from ParamFileTest<Simulation>

   System system_;

   /// Needed to create a stand-alone System object
   int    nAtomType_;

};



void SystemTest::setUp()
{

   #ifdef MCMD_ANGLE
   #ifdef MCMD_DIHEDRAL
   openFile("in/SimulationAngleDihedral"); 
   #else
   openFile("in/SimulationAngle"); 
   #endif
   #else
   #ifdef MCMD_DIHEDRAL
   openFile("in/SimulationDihedral"); 
   #else
   openFile("in/Simulation"); 
   #endif
   #endif

   object().readParam(file());
   system_.setId(0);
   nAtomType_ = 2;
   system_.setSimulation(object());
   system_.allocateMoleculeSets();
}

void SystemTest::addMolecules()
{
   Species *speciesPtr;
   int      i, j;
   for (i = 0; i < object().nSpecies() ; ++i ) {
      speciesPtr = &object().species(i);
      for (j = 0; j < speciesPtr->capacity(); ++j) {
         system_.addMolecule(speciesPtr->reservoir().pop());
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
   for (iSpecies=0; iSpecies < object().nSpecies() ; ++iSpecies ) {
      speciesPtr     = &object().species(iSpecies);
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
   for (iSpecies=0; iSpecies < object().nSpecies() ; ++iSpecies ) {
      speciesPtr = &object().species(iSpecies);
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
   for (int i=0; i < object().nSpecies() ; ++i ) {
      iMolecule = 0; 
      for (system_.begin(i, molIter); !molIter.atEnd(); ++molIter) {
         ++iMolecule;
      }
      TEST_ASSERT(iMolecule == system_.nMolecule(i));
   }

   // Output atomic positions
   if (verbose() > 1) {
      std::cout << std::endl;
      for (int i=0; i < object().nSpecies() ; ++i ) {
         std::cout << "species Id    = " << i << std::endl;
         for (system_.begin(i, molIter); !molIter.atEnd(); ++molIter) {
            std::cout << "Molecule Id    = " << molIter->id() << std::endl;
            for (molIter->begin(atomIter); !atomIter.atEnd(); ++atomIter) {
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
