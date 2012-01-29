#ifndef HARMONIC_L0_BOND_TEST_H
#define HARMONIC_L0_BOND_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <mcMd/potentials/bond/HarmonicL0Bond.h>
#include <mcMd/chemistry/Atom.h>
#include <util/random/Random.h>
#include <util/containers/RArray.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace McMd;

class HarmonicL0BondTest : public UnitTest 
{

private:

   HarmonicL0Bond bondPotential;

public:

   void setUp()
   {
      // Set Boundary Lengths
      bondPotential.setNBondType(2);

      // Read parameters from file
      std::ifstream in;
      openInputFile("bond/in/HarmonicL0Bond", in);
      bondPotential.readParam(in);
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


//   /*
//   void testEnergy2() 
//   {
//      printMethod(TEST_FUNC);
//      Vector r1, r2;
//      double rsq;
//      double energy;
//      int    i;
//      RArray<Atom> atoms;
//      
//      Atom::allocate(2, atoms);
//      Atom& atom1 = atoms[0];
//      Atom& atom2 = atoms[1];
//      bondPotential.writeParam(std::cout);
//      
//      r1 = Vector(1.0, 0.5, 0.5);
//      r2 = Vector(2.0, 0.5, 0.5);
//      atom1.setPosition(r1);
//      atom1.setTypeId(0);
//      atom2.setPosition(r2);
//      atom2.setTypeId(2);
//      rsq = 1.00;
//      i   = 0;
//      energy = bondPotential.energy(rsq, i);
//      std::cout << "bondPotential.energy(1.00, 0) = %lf\n", energy);
//
//      rsq = 0.81;
//      i   = 0;
//      energy = bondPotential.energy(rsq, i);
//      std::cout << "bondPotential.energy(0.81, 0) = %lf\n", energy);
//
//      rsq = 1.21000;
//      i   = 0;
//      energy = bondPotential.energy(rsq, i);
//      std::cout << "bondPotential.energy(1.21, 0) = %lf\n", energy);
//
//      rsq = 1.44;
//      i   = 0;
//      energy = bondPotential.energy(rsq, i);
//      std::cout << "bondPotential.energy(1.44, 0) = %lf\n", energy);
//   }
//   */


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
 
      random = new Random;
      std::ifstream in;
      openInputFile("bond/in/Random", in);
      random->readParam(in);

      beta = 1.0;
      type = 0;
      for (i=0; i < 20; i++) {
         r = bondPotential.randomBondLength(random, 1.0, 0);
         std::cout << "random bond length = " << r << std::endl;
      }
   }
     

};

TEST_BEGIN(HarmonicL0BondTest)
TEST_ADD(HarmonicL0BondTest, testSetUp)
TEST_ADD(HarmonicL0BondTest, testWrite)
TEST_ADD(HarmonicL0BondTest, testEnergy)
TEST_ADD(HarmonicL0BondTest, testForceOverR)
   //TEST_ADD(HarmonicL0BondTest, testRandomBondLength);
TEST_END(HarmonicL0BondTest)

#endif
