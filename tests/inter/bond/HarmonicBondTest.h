#ifndef HARMONIC_BOND_TEST_H
#define HARMONIC_BOND_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <inter/bond/HarmonicBond.h>
#include <mcMd/chemistry/Atom.h>
#include <util/random/Random.h>
#include <util/containers/RArray.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Inter;

class HarmonicBondTest : public UnitTest 
{

private:

   HarmonicBond bondPotential;

public:

   void setUp()
   {
      // Set Boundary Lengths
      bondPotential.setNBondType(2);

      // Read parameters from file
      std::ifstream in;
      openInputFile("bond/in/HarmonicBond", in);
      bondPotential.readParameters(in);
      in.close();
   }


   void tearDown()
   {}


   void testSetUp() 
   {
      printMethod(TEST_FUNC);
   }


   void testWrite() {
      printMethod(TEST_FUNC);

      // Verbose output
      if (verbose() > 0) { 
         bondPotential.writeParam(std::cout);
      }
   }


   void testEnergy() 
   {
      int    i;
      double rsq;
      double energy;
     
      printMethod(TEST_FUNC);
      
      rsq = 1.00;
      i   = 0;
      energy = bondPotential.energy(rsq, i);
      if (verbose() > 1) {
         std::cout << "bondPotential.energy(1.00, 0) = " << energy << std::endl;
      }

      rsq = 0.81;
      i   = 0;
      energy = bondPotential.energy(rsq, i);
      if (verbose() > 1) {
         std::cout << "bondPotential.energy(0.81, 0) = " << energy << std::endl;
      }

      rsq = 1.21000;
      i   = 0;
      energy = bondPotential.energy(rsq, i);
      if (verbose() > 1) {
         std::cout << "bondPotential.energy(1.21, 0) = " << energy << std::endl;
      }

      rsq = 1.44;
      i   = 0;
      energy = bondPotential.energy(rsq, i);
      if (verbose() > 1) {
         std::cout << "bondPotential.energy(1.44, 0) = " << energy << std::endl;
      }
   }

   void testForceOverR() 
   {
      int    i;
      double rsq;
      double force;
 
      printMethod(TEST_FUNC);
     
      rsq = 1.00;
      i   = 0;
      force = bondPotential.forceOverR(rsq, i);
      if (verbose() > 1) {
         std::cout << "bondPotential.forceOverR(1.00, 0) = " << force << std::endl;
      }

      rsq = 0.81;
      i   = 0;
      force = bondPotential.forceOverR(rsq, i);
      if (verbose() > 1) {
         std::cout << "bondPotential.forceOverR(0.81, 0) = " << force << std::endl;
      }

      rsq = 1.21; 
      i   = 0;
      force = bondPotential.forceOverR(rsq, i);
      if (verbose() > 1) {
         std::cout << "bondPotential.forceOverR(1.21, 0) = " << force << std::endl;
      }

      rsq = 1.44;
      i   = 0;
      force = bondPotential.forceOverR(rsq, i);
      if (verbose() > 1) {
         std::cout << "bondPotential.forceOverR(1.44, 0) = " << force << std::endl;
      }
   }


   void testRandomBondLength() 
   {
      int    type, i;
      double r, beta;
  
      Random *random;
 
      printMethod(TEST_FUNC);
 
      std::ifstream in;
      openInputFile("potentials/bond/in/Random", in);
      random = new Random;
      random->readParam(in);

      beta = 1.0;
      type = 0;
      for (i=0; i < 20; i++) {
         r = bondPotential.randomBondLength(random, beta, type);
         std::cout << "random bond length = " << r << std::endl;
      }
   }
     

};

TEST_BEGIN(HarmonicBondTest)
TEST_ADD(HarmonicBondTest, testSetUp)
TEST_ADD(HarmonicBondTest, testWrite)
TEST_ADD(HarmonicBondTest, testEnergy)
TEST_ADD(HarmonicBondTest, testForceOverR)
   //TEST_ADD(HarmonicBond, testRandomBondLength);
TEST_END(HarmonicBondTest)


#endif
