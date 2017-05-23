#ifndef MCMD_MC_SIMULATION_TEST_H
#define MCMD_MC_SIMULATION_TEST_H

#include <mcMd/mcSimulation/McSimulation.h>
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/mdSimulation/MdSystem.h>
#include <simp/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <mcMd/chemistry/Activate.h>
#include <mcMd/potentials/pair/McPairPotential.h>
#include <mcMd/potentials/bond/BondPotential.h>
#ifdef SIMP_ANGLE
#include <mcMd/potentials/angle/AnglePotential.h>
#endif

#include <util/archives/BinaryFileOArchive.h>
#include <util/archives/BinaryFileIArchive.h>
#include <util/archives/Serializable_includes.h>

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

#include <string>
#include <fstream>

using namespace Util;
using namespace McMd;

class McSimulationTest : public ParamFileTest
{

public:

   McSimulationTest();
   virtual void setUp();

   // Test functions
   void testReadParamBond();
   void testReadConfigBond();
   void testPairEnergy();
   void testBondEnergy();
   void testActivate();
   void testMdSystemCopy();
   void testSimulateBond();
   void testWriteRestartBond();
   void testReadRestart();

   #ifdef SIMP_ANGLE
   void testReadParamAngle();
   void testAngleEnergy();
   void testSimulateAngle();
   #endif

private:

   McSimulation simulation_;
   McSystem& system_;

   // Utility functions
   void readParam(const char* filename);
   void readConfig(const char* filename);

};

McSimulationTest::McSimulationTest()
 : ParamFileTest(),
   system_(simulation_.system())
{ 
   //setVerbose(2); 
   //ParamComposite::setEcho(true);
}

void McSimulationTest::setUp()
{
   simulation_.fileMaster().setRootPrefix(filePrefix());
} 

void McSimulationTest::readParam(const char* filename)
{  
   openFile(filename); 
   simulation_.readParam(file());
   file().close();
}

void McSimulationTest::readConfig(const char* filename)
{  
   openFile(filename); 
   system_.readConfig(file());
   file().close();
}

// Test methods

void McSimulationTest::testReadParamBond()
{
   printMethod(TEST_FUNC);
   readParam("in/McSimulation");

   try {
      simulation_.isValid();
   } catch (Exception e) {
      std::cout << e.message();
      TEST_ASSERT(0);
   }

   if (verbose() > 1) {
      std::cout << std::endl;
      simulation_.writeParam(std::cout);
   }
}

void McSimulationTest::testReadConfigBond()
{
   printMethod(TEST_FUNC);

   readParam("in/McSimulation");
   try {
      readConfig("in/config");
      simulation_.isValid();
   } catch (Exception e) {
      std::cout << e.message();
      TEST_ASSERT(0);
   }
}

#ifdef SIMP_ANGLE
void McSimulationTest::testReadParamAngle()
{
   printMethod(TEST_FUNC);

   readParam("in/McSimulationAngle"); 
   if (verbose() > 1) {
      std::cout << std::endl;
      simulation_.writeParam(std::cout);
   }
}
#endif

void McSimulationTest::testPairEnergy()
{ 
   printMethod(TEST_FUNC);
   std::cout << std::endl;

   readParam("in/McSimulation"); 
   readConfig("in/config");

   // simulation_.simulate(10);

   System::MoleculeIterator molIter;
   Molecule::AtomIterator atomIter;
   double total = system_.pairPotential().energy();
   double de;
   double energy = 0.0;
   for (int is=0; is < simulation_.nSpecies(); ++is) {
      for (system_.begin(is, molIter); molIter.notEnd(); ++molIter) {
         for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
            de = system_.pairPotential().atomEnergy(*atomIter);
            //std::cout.width(5);
            //std::cout << atomIter->id() << "     " << de << std::endl;
            energy += de;
         }
      }
   }
   if (verbose() > 1) {
      std::cout << "Total atomPairEnergy = " << 0.5*energy << std::endl;
      std::cout << "Total PairEnergy     = " << total << std::endl;
   }
   TEST_ASSERT(eq(0.5*energy, total));
}

void McSimulationTest::testBondEnergy()
{ 
   printMethod(TEST_FUNC);
   std::cout << std::endl;

   readParam("in/McSimulation"); 
   readConfig("in/config");

   //simulation_.simulate(10);

   double total = system_.bondPotential().energy();
   System::MoleculeIterator molIter;
   Molecule::AtomIterator atomIter;
   double energy, de;
   energy = 0.0;
   for (int is=0; is < simulation_.nSpecies(); ++is) {
      for (system_.begin(is, molIter); molIter.notEnd(); ++molIter) {
         for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
            de = system_.bondPotential().atomEnergy(*atomIter);
            // std::cout.width(5);
            // std::cout << atomIter->id() << "     " << de << std::endl;
            energy += de;
         }
      }
   }
   if (verbose() > 1) {
      std::cout << "Total atomBondEnergy = " << 0.5*energy << std::endl;
      std::cout << "Total bondEnergy     = " << total << std::endl;
   }
   TEST_ASSERT(eq(0.5*energy, total));
}

#ifdef SIMP_ANGLE
void McSimulationTest::testAngleEnergy()
{ 
   printMethod(TEST_FUNC);
   std::cout << std::endl;

   readParam("in/McSimulationAngle"); 
   readConfig("in/config");

   //simulation_.simulate(10);

   double total = system_.anglePotential().energy();

   System::MoleculeIterator molIter;
   Molecule::AtomIterator atomIter;
   double energy, de;
   energy = 0.0;
   for (int is=0; is < simulation_.nSpecies(); ++is) {
      if (simulation_.species(is).nAngle() > 0) {
         for (system_.begin(is, molIter); molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               de = system_.anglePotential().atomEnergy(*atomIter);
               // std::cout.width(5);
               // std::cout << atomIter->id() << "     " << de << std::endl;
               energy += de;
            }
         }
      }
   }
   if (verbose() > 1) {
      std::cout << "Total angleEnergy     = " << total << std::endl;
      std::cout << "Total atomAngleEnergy = " << energy/3.0 << std::endl;
   }
   TEST_ASSERT(eq(energy/3.0, total));
}
#endif

void McSimulationTest::testActivate()
{ 
   printMethod(TEST_FUNC);
   std::cout << std::endl;

   readParam("in/McSimulation"); 
   readConfig("in/config");

   int speciesId = 1;
   int moleculeId = 1;
   int atomId = 2;
   Species& species = simulation_.species(speciesId);
   Molecule& molecule = system_.molecule(speciesId, moleculeId);
   Atom& atom = molecule.atom(atomId);
   Species::AtomBondIdArray bondIds = species.atomBondIds(atomId);
 
   // Test initial state 
   TEST_ASSERT(atom.isActive());
   for (int i = 0; i < bondIds.size(); ++i) {
      TEST_ASSERT(molecule.bond(bondIds[i]).isActive());
   }
   for (int i = 0; i < molecule.nBond(); ++i) {
      TEST_ASSERT(molecule.bond(i).isActive());
      TEST_ASSERT(molecule.bond(i).checkInactive());
   }

   // Deactivate atom atomId 
   Activate::deactivate(molecule.atom(atomId));
   TEST_ASSERT(!atom.isActive());
   for (int i = 0; i < bondIds.size(); ++i) {
      TEST_ASSERT(!molecule.bond(bondIds[i]).isActive());
      TEST_ASSERT(molecule.bond(bondIds[i]).nInActive() == 1);
   }
   for (int i = 0; i < molecule.nBond(); ++i) {
      TEST_ASSERT(molecule.bond(i).checkInactive());
   }

   // Deactivate atom atomId + 1
   Activate::deactivate(molecule.atom(atomId+1));
   TEST_ASSERT(!atom.isActive());
   TEST_ASSERT(!molecule.atom(atomId+1).isActive());
   for (int i = 0; i < bondIds.size(); ++i) {
      TEST_ASSERT(!molecule.bond(bondIds[i]).isActive());
   }
   for (int i = 0; i < molecule.nBond(); ++i) {
      TEST_ASSERT(molecule.bond(i).checkInactive());
   }

   // Reactivate both atoms
   Activate::reactivate(molecule.atom(atomId));
   Activate::reactivate(molecule.atom(atomId+1));
   TEST_ASSERT(atom.isActive());
   TEST_ASSERT(molecule.atom(atomId+1).isActive());
   for (int i = 0; i < bondIds.size(); ++i) {
      TEST_ASSERT(molecule.bond(bondIds[i]).isActive());
      TEST_ASSERT(molecule.bond(bondIds[i]).nInActive() == 0);
   }
   for (int i = 0; i < molecule.nBond(); ++i) {
      TEST_ASSERT(molecule.bond(i).isActive());
      TEST_ASSERT(molecule.bond(i).checkInactive());
   }
}

void McSimulationTest::testMdSystemCopy()
{
   printMethod(TEST_FUNC);
   std::cout << std::endl;

   readParam("in/McSimulation"); 
   readConfig("in/config");

   simulation_.simulate(10);

   MdSystem mdSystem(system_);

   std::ifstream mdSystemFile;
   openInputFile("in/MdSystemCopy", mdSystemFile); 

   mdSystem.readParam(mdSystemFile);

   std::cout << "MC Potential Energy      = " 
             << system_.potentialEnergy() << std::endl;
   std::cout << "MD Potential Energy      = " 
             << mdSystem.potentialEnergy() << std::endl;

   mdSystemFile.close(); 
}

void McSimulationTest::testSimulateBond()
{
   printMethod(TEST_FUNC);

   readParam("in/McSimulation"); 
   readConfig("in/config"); 

   std::string baseFileName("simulate.0");
   simulation_.save(baseFileName);

   simulation_.simulate(20);

   baseFileName = "simulate.20";
   simulation_.save(baseFileName);

}

#ifdef SIMP_ANGLE
void McSimulationTest::testSimulateAngle()
{
   printMethod(TEST_FUNC);

   readParam("in/McSimulationAngle"); 
   readConfig("in/md.config");

   std::string baseFileName("simulateAngle.0");
   simulation_.save(baseFileName);

   simulation_.simulate(40);

   baseFileName = "simulateAngle.20";
   simulation_.save(baseFileName);
}
#endif

void McSimulationTest::testWriteRestartBond()
{
   printMethod(TEST_FUNC);
   std::cout << std::endl;

   readParam("in/McSimulation"); 
   readConfig("in/config"); 

   std::string baseFileName("writeRestart.0");
   simulation_.save(baseFileName);

   simulation_.simulate(10);
   baseFileName = "writeRestart.10";
   simulation_.save(baseFileName);

   bool isContinuation = true;
   simulation_.simulate(20, isContinuation);

   baseFileName = "writeRestart.20";
   simulation_.save(baseFileName);
}

void McSimulationTest::testReadRestart()
{
   printMethod(TEST_FUNC);
   std::cout << std::endl;

   std::string baseFileName("writeRestart.10");
   simulation_.load(baseFileName);

   baseFileName = "readRestart";
   simulation_.save(baseFileName);
}

TEST_BEGIN(McSimulationTest)
TEST_ADD(McSimulationTest, testReadParamBond)
TEST_ADD(McSimulationTest, testReadConfigBond)
TEST_ADD(McSimulationTest, testPairEnergy)
TEST_ADD(McSimulationTest, testBondEnergy)
TEST_ADD(McSimulationTest, testActivate)
TEST_ADD(McSimulationTest, testMdSystemCopy)
TEST_ADD(McSimulationTest, testSimulateBond)
TEST_ADD(McSimulationTest, testWriteRestartBond)
//TEST_ADD(McSimulationTest, testReadRestart)
#ifdef SIMP_ANGLE
TEST_ADD(McSimulationTest, testReadParamAngle)
TEST_ADD(McSimulationTest, testAngleEnergy)
TEST_ADD(McSimulationTest, testSimulateAngle)
#endif
TEST_END(McSimulationTest)

#endif
