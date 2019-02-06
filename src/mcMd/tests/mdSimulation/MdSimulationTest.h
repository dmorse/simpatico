#ifndef MCMD_MD_SIMULATION_TEST_H
#define MCMD_MD_SIMULATION_TEST_H

#include <mcMd/mdSimulation/MdSimulation.h>
#include <mcMd/mdSimulation/MdSystem.h>
#include <mcMd/mdIntegrators/MdIntegrator.h>
#include <mcMd/potentials/pair/MdPairPotential.h>
#include <mcMd/potentials/bond/BondPotential.h>
#include <simp/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>

#include <test/UnitTestRunner.h>
#include <test/UnitTest.h>

using namespace Util;
using namespace McMd;

class MdSimulationTest : public UnitTest
{

public:

   MdSimulationTest();
   virtual void setUp();

   void testReadParam();
   void testSetZeroVelocities();
   void testSetBoltzmannVelocities();
   void testBuildPairList();
   void testPairEnergy();
   void testAddPairForces();
   void testBondEnergy();
   void testAddBondForces();
   void testCalculateForces();
   void testStep();
   void testSimulate();
   void testWriteRestart();
   void testReadRestart();

private:

   MdSimulation simulation_;
   MdSystem&    system_;

};


MdSimulationTest::MdSimulationTest()
 : UnitTest(),
   system_(simulation_.system())
{} 

void MdSimulationTest::setUp()
{  
   simulation_.fileMaster().setRootPrefix(filePrefix()); 
   //setVerbose(2);
}

void MdSimulationTest::testReadParam()
{
   if (isIoProcessor()) {
      printMethod(TEST_FUNC);
      std::cout << std::endl;
   
      std::ifstream paramFile;
      openInputFile("in/MdSimulation", paramFile); 
      //ParamComponent::setEcho(true);
      simulation_.readParam(paramFile);
      //ParamComponent::setEcho(false);
      simulation_.readCommands();
      std::cout << std::endl;
   
      try {
         simulation_.isValid();
      } catch (Exception e) {
         TEST_ASSERT(0);
      }
   
      if (verbose() > 1) {
         std::cout << std::endl;
         simulation_.writeParam(std::cout);
         simulation_.system().writeConfig(std::cout);
      }
   }
}

void MdSimulationTest::testSetZeroVelocities()
{ 
   if (isIoProcessor()) {
      printMethod(TEST_FUNC);
      std::cout << std::endl;
   
      // Read the parameter file
      std::ifstream paramFile;
      openInputFile("in/MdSimulation", paramFile); 
      simulation_.readParam(paramFile);
      paramFile.close();
      simulation_.readCommands();
      std::cout << std::endl;
   
      double energy;
      energy = system_.kineticEnergy();
      std::cout << "kinetic energy = " << energy << std::endl;
   
      system_.setZeroVelocities();
      energy = system_.kineticEnergy();
      std::cout << "kinetic energy = " << energy << std::endl;
   }
}

void MdSimulationTest::testSetBoltzmannVelocities()
{ 
   if (isIoProcessor()) {
      printMethod(TEST_FUNC);
      std::cout << std::endl;
   
      std::ifstream paramFile;
      openInputFile("in/MdSimulation", paramFile); 
      simulation_.readParam(paramFile);
      paramFile.close();
      simulation_.readCommands();
      std::cout << std::endl;
   
      double energy;
      energy = system_.kineticEnergy();
      std::cout << "kinetic energy = " << energy << std::endl;
   
      double temperature = 1.0;
      system_.setBoltzmannVelocities(temperature);
   
      energy = system_.kineticEnergy();
      std::cout << "kinetic energy = " << energy << std::endl;
   }
}

void MdSimulationTest::testBuildPairList()
{ 
   if (isIoProcessor()) {
      printMethod(TEST_FUNC);
      std::cout << std::endl;
   
      std::ifstream paramFile;
      openInputFile("in/MdSimulation", paramFile); 
      simulation_.readParam(paramFile);
      paramFile.close();
      simulation_.readCommands();
      if (isIoProcessor()) std::cout << std::endl;
   
      system_.pairPotential().buildPairList();
   
      //bool isContinuation = false;
      //simulation_.simulate(2, isContinuation);
   
      try {
         simulation_.isValid();
      } catch (Exception e) {
         std::cout << e.message();
         TEST_ASSERT(0);
      }
   }
}


void MdSimulationTest::testPairEnergy()
{ 
   if (isIoProcessor()) {
      printMethod(TEST_FUNC);
      std::cout << std::endl;
   
      std::ifstream paramFile;
      openInputFile("in/MdSimulation", paramFile); 
      simulation_.readParam(paramFile);
      paramFile.close();
      simulation_.readCommands();
      std::cout << std::endl;
   
      system_.pairPotential().buildPairList();
      double energy = system_.pairPotential().energy();
      std::cout << "Pair energy: " << energy << std::endl;

      try {
         simulation_.isValid();
      } catch (Exception e) {
         std::cout << e.message();
         TEST_ASSERT(0);
      }
   }
}

void MdSimulationTest::testAddPairForces()
{ 
   if (isIoProcessor()) {
      printMethod(TEST_FUNC);
      std::cout << std::endl;
   
      std::ifstream paramFile;
      openInputFile("in/MdSimulation", paramFile); 
      simulation_.readParam(paramFile);
      paramFile.close();
      simulation_.readCommands();
      if (isIoProcessor()) std::cout << std::endl;
   
      system_.pairPotential().buildPairList();
      system_.setZeroForces();
      system_.pairPotential().addForces();

      try {
         simulation_.isValid();
      } catch (Exception e) {
         std::cout << e.message();
         TEST_ASSERT(0);
      }
   }
}

void MdSimulationTest::testBondEnergy()
{ 
   if (isIoProcessor()) {
      printMethod(TEST_FUNC);
      std::cout << std::endl;
   
      std::ifstream paramFile;
      openInputFile("in/MdSimulation", paramFile); 
      simulation_.readParam(paramFile);
      paramFile.close();
      simulation_.readCommands();
      std::cout << std::endl;
   
      double energy = system_.bondPotential().energy();
      std::cout << "Bond energy = " << energy << std::endl;
   }
}

void MdSimulationTest::testAddBondForces()
{ 
   if (isIoProcessor()) {
      printMethod(TEST_FUNC);
   
      std::ifstream paramFile;
      openInputFile("in/MdSimulation", paramFile); 
      simulation_.readParam(paramFile);
      paramFile.close();
      simulation_.readCommands();
      if (isIoProcessor()) std::cout << std::endl;
   
      //double temperature = 1.0;
      //system_.setBoltzmannVelocities(temperature);
      //simulation_.simulate(1000);
   
      system_.setZeroForces();
      system_.bondPotential().addForces();
   }
}

void MdSimulationTest::testCalculateForces()
{ 
   if (isIoProcessor()) {
      printMethod(TEST_FUNC);
   
      std::ifstream paramFile;
      openInputFile("in/MdSimulation", paramFile); 
      simulation_.readParam(paramFile);
      paramFile.close();
      simulation_.readCommands();
      std::cout << std::endl;
   
      system_.pairPotential().buildPairList();
      system_.calculateForces();

      try {
         simulation_.isValid();
      } catch (Exception e) {
         std::cout << e.message();
         TEST_ASSERT(0);
      }
   }
}

void MdSimulationTest::testStep()
{
   if (isIoProcessor()) {
      printMethod(TEST_FUNC);
      if (isIoProcessor()) std::cout << std::endl;
   
      std::ifstream paramFile;
      openInputFile("in/MdSimulation", paramFile); 
      simulation_.readParam(paramFile);
      paramFile.close();
      simulation_.readCommands();
   
      std::cout << std::endl;
   
      double kinetic, potential;
      system_.pairPotential().buildPairList();
      system_.calculateForces();
      system_.mdIntegrator().setup();
      for (int i=0; i < 10; ++i) {
   
         kinetic   = system_.kineticEnergy(); 
         potential = system_.potentialEnergy(); 
         std::cout << kinetic << "  " << potential 
                   << "  " << kinetic + potential << std::endl;
   
         system_.mdIntegrator().step();
      }
   
      try {
         simulation_.isValid();
      } catch (Exception e) {
         std::cout << e.message();
         TEST_ASSERT(0);
      }

      kinetic   = system_.kineticEnergy(); 
      potential = system_.potentialEnergy(); 
      std::cout << kinetic << "  " << potential << "  " 
                << kinetic + potential << std::endl;
   }
}

void MdSimulationTest::testSimulate()
{
   if (isIoProcessor()) {
      printMethod(TEST_FUNC);
      std::cout << std::endl;
   
      std::ifstream paramFile;
      openInputFile("in/MdSimulation", paramFile); 
      simulation_.readParam(paramFile);
      paramFile.close();
      simulation_.readCommands();
      std::cout << std::endl;
   
      //simulation_.save("simulate.0");
   
      simulation_.simulate(2000);
   
      //simulation_.save("simulate.2000");
      try {
         simulation_.isValid();
      } catch (Exception e) {
         std::cout << e.message();
         TEST_ASSERT(0);
      }
   }
}

void MdSimulationTest::testWriteRestart()
{
   if (isIoProcessor()) {
      printMethod(TEST_FUNC);
      std::cout << std::endl;
   
      std::ifstream paramFile;
      openInputFile("in/MdSimulation", paramFile); 
      simulation_.readParam(paramFile);
      paramFile.close();
      simulation_.readCommands();
      std::cout << std::endl;
   
      // Save initial state to tmp/begin
      simulation_.save("tmp/begin.rst");
   
      // Run simulation of 10000 steps
      simulation_.simulate(10000);
   
      // Save state at iStep = 10000 to file middle.rst
      simulation_.save("tmp/middle.rst");
   
      // Write configuration at iStep = 1000 to middle.cfg
      std::ofstream configFile("tmp/middle.cfg");
      simulation_.system().writeConfig(configFile);
      configFile.close();
   
      bool isContinuation = true;
      simulation_.simulate(10100, isContinuation);
   
      // Write configuration at iStep = 101000 
      configFile.open("tmp/end.cfg");
      simulation_.system().writeConfig(configFile);
      configFile.close();

      try {
         simulation_.isValid();
      } catch (Exception e) {
         std::cout << e.message();
         TEST_ASSERT(0);
      }
   }

}

void MdSimulationTest::testReadRestart()
{
   if (isIoProcessor()) {
      printMethod(TEST_FUNC);
      std::cout << std::endl;
   
      // Restart after iStep = 10000
      simulation_.load("tmp/middle.rst");
   
      // Save new restart file at iStep = 10000
      simulation_.save("tmp/middle2.rst");
   
      // Save config file at iStep = 10000
      std::ofstream configFile("tmp/middle2.cfg");
      simulation_.system().writeConfig(configFile);
      configFile.close();
   
      // Run to iStep = 10100
      bool isContinuation = true;
      simulation_.simulate(10100, isContinuation);
   
      configFile.open("tmp/end2.cfg");
      simulation_.system().writeConfig(configFile);
      configFile.close();
   }
}

TEST_BEGIN(MdSimulationTest)
TEST_ADD(MdSimulationTest, testReadParam)
TEST_ADD(MdSimulationTest, testSetZeroVelocities)
TEST_ADD(MdSimulationTest, testSetBoltzmannVelocities)
TEST_ADD(MdSimulationTest, testBuildPairList)
TEST_ADD(MdSimulationTest, testPairEnergy)
TEST_ADD(MdSimulationTest, testAddPairForces)
TEST_ADD(MdSimulationTest, testBondEnergy)
TEST_ADD(MdSimulationTest, testAddBondForces)
TEST_ADD(MdSimulationTest, testCalculateForces)
TEST_ADD(MdSimulationTest, testStep)
TEST_ADD(MdSimulationTest, testSimulate)
TEST_ADD(MdSimulationTest, testWriteRestart)
TEST_ADD(MdSimulationTest, testReadRestart)
TEST_END(MdSimulationTest)

#endif
