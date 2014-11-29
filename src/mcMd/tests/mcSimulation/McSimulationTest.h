#ifndef MCMD_MC_SIMULATION_TEST_H
#define MCMD_MC_SIMULATION_TEST_H

#include <mcMd/mcSimulation/McSimulation.h>
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/mdSimulation/MdSystem.h>
//#include <mcMd/mcSimulation/serialize.h>
#include <mcMd/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <mcMd/potentials/pair/McPairPotential.h>
#include <mcMd/potentials/bond/BondPotential.h>
#ifdef INTER_ANGLE
#include <mcMd/potentials/angle/AnglePotential.h>
#endif

#include <util/archives/MemoryCounter.h>
#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>
#include <util/archives/TextFileOArchive.h>
#include <util/archives/TextFileIArchive.h>
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
   void readParam();
   void readConfig(const char* filename);

   void testReadParamBond();
   void testReadParam();
   void testPairEnergy();
   void testBondEnergy();
   void testMdSystemCopy();
   void testSimulateBond();
   void testWriteRestartBond();
   void testReadRestart();

   #ifdef INTER_ANGLE
   void testAngleEnergy();
   void testSimulateAngle();
   #endif

private:

   McSimulation simulation_;
   McSystem& system_;

};

McSimulationTest::McSimulationTest()
 : ParamFileTest(),
   system_(simulation_.system())
{ 
   setVerbose(2); 
}

void McSimulationTest::setUp()
{
   simulation_.fileMaster().setRootPrefix(filePrefix());
} 

void McSimulationTest::readParam()
{  
   #ifdef INTER_ANGLE
   openFile("in/McSimulationAngle"); 
   #else
   openFile("in/McSimulation"); 
   #endif
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
   openFile("in/McSimulation");
   simulation_.readParam(file());
   file().close();

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

void McSimulationTest::testReadParam()
{
   printMethod(TEST_FUNC);
   //ParamComposite::setEcho(true);

   readParam();
   if (verbose() > 1) {
      std::cout << std::endl;
      simulation_.writeParam(std::cout);
   }

   #if 1
   simulation_.readCommands();
   try {
      simulation_.isValid();
   } catch (Exception e) {
      std::cout << e.message();
      TEST_ASSERT(0);
   }
   #endif
}

void McSimulationTest::testPairEnergy()
{ 
   printMethod(TEST_FUNC);
   std::cout << std::endl;

   //readParam();
   openFile("in/McSimulation"); 
   simulation_.readParam(file());
   simulation_.readCommands();

   simulation_.simulate(10);

   System::MoleculeIterator molIter;
   Molecule::AtomIterator atomIter;
   double total = system_.pairPotential().energy();

   double energy, de;
   energy = 0.0;
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

   openFile("in/McSimulation"); 
   simulation_.readParam(file());
   file().close();
   simulation_.readCommands();

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

#ifdef INTER_ANGLE
void McSimulationTest::testAngleEnergy()
{ 
   printMethod(TEST_FUNC);
   std::cout << std::endl;

   //readParam(simulation_);
   openFile("in/McSimulationAngle"); 
   simulation_.readParam(file());
   file().close();
   simulation_.readCommands();

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

void McSimulationTest::testMdSystemCopy()
{
   printMethod(TEST_FUNC);
   std::cout << std::endl;

   //readParam(simulation_);
   openFile("in/McSimulation"); 
   simulation_.readParam(file());
   file().close();
   simulation_.readCommands();

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

   std::cout << std::endl;
   openFile("in/McSimulation"); 
   simulation_.readParam(file());
   file().close();
   simulation_.readCommands();

   std::string baseFileName("simulate.0");
   simulation_.save(baseFileName);

   simulation_.simulate(20);

   baseFileName = "simulate.20";
   simulation_.save(baseFileName);

}

void McSimulationTest::testWriteRestartBond()
{
   printMethod(TEST_FUNC);

   std::cout << std::endl;
   //readParam(simulation_);
   //simulation_.readCommands();
   openFile("in/McSimulation"); 
   simulation_.readParam(file());
   file().close();
   simulation_.readCommands();

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

#ifdef INTER_ANGLE
void McSimulationTest::testSimulateAngle()
{
   printMethod(TEST_FUNC);
   std::cout << std::endl;
   openFile("in/McSimulationAngle"); 
   simulation_.readParam(file());
   file().close();

   //simulation_.readCommands();
   readConfig("in/md.config");

   std::cout << std::endl;

   std::string baseFileName("simulateAngle.0");
   simulation_.save(baseFileName);

   simulation_.simulate(40);

   baseFileName = "simulateAngle.20";
   simulation_.save(baseFileName);

}
#endif

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
TEST_ADD(McSimulationTest, testReadParam)
TEST_ADD(McSimulationTest, testPairEnergy)
TEST_ADD(McSimulationTest, testBondEnergy)
TEST_ADD(McSimulationTest, testMdSystemCopy)
TEST_ADD(McSimulationTest, testSimulateBond)
TEST_ADD(McSimulationTest, testWriteRestartBond)
//TEST_ADD(McSimulationTest, testReadRestart)
#ifdef INTER_ANGLE
TEST_ADD(McSimulationTest, testAngleEnergy)
TEST_ADD(McSimulationTest, testSimulateAngle)
#endif
TEST_END(McSimulationTest)

#endif
