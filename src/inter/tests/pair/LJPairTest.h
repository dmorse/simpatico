#ifndef LJ_PAIR_TEST_H
#define LJ_PAIR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <inter/pair/LJPair.h>
#include <util/containers/RArray.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Inter;

class LJPairTest : public UnitTest 
{

private:

   int     nAtomType;
   LJPair  pairPotential;

   //const static double eps;

public:

   void setUp()
   {

      // Set number of Atom types.
      nAtomType = 2;
      pairPotential.setNAtomType(nAtomType);

      // Read parameters from file
      std::ifstream in;
      openInputFile("pair/in/LJPair", in);
      pairPotential.readParam(in);
      in.close();
   }

   void tearDown(){
   }

   void testSetUp() {
      printMethod(TEST_FUNC);
   }


   void testWrite() {

      printMethod(TEST_FUNC);

      // Verbose output
      if (verbose() > 0) {
         std::cout << std::endl; 
         pairPotential.writeParam(std::cout);
      }

   }


   void testEnergy1() {
      int    i, j;
      double rsq;
      double energy;
     
      printMethod(TEST_FUNC);

      rsq = 1.00;
      i   = 0;
      j   = 1;
      energy = pairPotential.energy(rsq, i, j);

      std::cout << std::endl; 
      std::cout << "pairPotential.energy(1.00, 0, 1) = " << energy << std::endl;

      rsq = 0.81;
      i   = 0;
      j   = 1;
      energy = pairPotential.energy(rsq, i, j);
      std::cout << "pairPotential.energy(0.81, 0, 1) = " << energy << std::endl;

      rsq = 1.25992;
      i   = 0;
      j   = 1;
      energy = pairPotential.energy(rsq, i, j);
      std::cout << "pairPotential.energy(1.15992, 0, 1) = " << energy << std::endl;

      rsq = 1.50;
      i   = 0;
      j   = 1;
      energy = pairPotential.energy(rsq, i, j);
      std::cout << "pairPotential.energy(1.50, 0, 1) = " << energy << std::endl;
   }


   void testForceOverR() {
      int    i, j;
      double rsq;
      double force;
     
      printMethod(TEST_FUNC);
     
      rsq = 1.00;
      i   = 0;
      j   = 1;
      force = pairPotential.forceOverR(rsq, i, j);
      std::cout << std::endl;
      std::cout << "pairPotential.forceOverR(1.00, 0, 1) = " << force << std::endl;

      rsq = 0.81;
      i   = 0;
      j   = 1;
      force = pairPotential.forceOverR(rsq, i, j);
      std::cout << "pairPotential.forceOverR(0.81, 0, 1) = " << force << std::endl;

      rsq = 1.25990; // Slightly less than cutoff**2 = 1.259916452
      i   = 0;
      j   = 1;
      force = pairPotential.forceOverR(rsq, i, j);
      std::cout << "pairPotential.forceOverR(1.15990, 0, 1) = " << force << std::endl;

      rsq = 1.50;
      i   = 0;
      j   = 1;
      force = pairPotential.forceOverR(rsq, i, j);
      std::cout << "pairPotential.forceOverR(1.50, 0, 1) = " << force << std::endl;
   }

   void testGetSet() {
      printMethod(TEST_FUNC);

      TEST_ASSERT(eq(pairPotential.epsilon(0, 0), 1.0));
      TEST_ASSERT(eq(pairPotential.epsilon(1, 1), 1.0));
      TEST_ASSERT(eq(pairPotential.epsilon(0, 1), 2.0));
      TEST_ASSERT(eq(pairPotential.epsilon(1, 0), 2.0));
      TEST_ASSERT(eq(pairPotential.sigma(0, 0), 1.0));
      TEST_ASSERT(eq(pairPotential.sigma(1, 1), 1.0));
      TEST_ASSERT(eq(pairPotential.sigma(0, 1), 1.0));
      TEST_ASSERT(eq(pairPotential.sigma(1, 0), 1.0));

      pairPotential.setEpsilon(0, 1, 1.3);
      TEST_ASSERT(eq(pairPotential.epsilon(0, 0), 1.0));
      TEST_ASSERT(eq(pairPotential.epsilon(1, 1), 1.0));
      TEST_ASSERT(eq(pairPotential.epsilon(0, 1), 1.3));
      TEST_ASSERT(eq(pairPotential.epsilon(1, 0), 1.3));
      
      pairPotential.setEpsilon(0, 0, 1.1);
      TEST_ASSERT(eq(pairPotential.epsilon(0, 0), 1.1));
      TEST_ASSERT(eq(pairPotential.epsilon(1, 1), 1.0));
      TEST_ASSERT(eq(pairPotential.epsilon(0, 1), 1.3));
      TEST_ASSERT(eq(pairPotential.epsilon(1, 0), 1.3));
      
      pairPotential.setSigma(0, 1, 1.05);
      TEST_ASSERT(eq(pairPotential.sigma(0, 0), 1.0));
      TEST_ASSERT(eq(pairPotential.sigma(1, 1), 1.0));
      TEST_ASSERT(eq(pairPotential.sigma(0, 1), 1.05));
      TEST_ASSERT(eq(pairPotential.sigma(1, 0), 1.05));
      

   }

};

TEST_BEGIN(LJPairTest)
TEST_ADD(LJPairTest, testSetUp)
TEST_ADD(LJPairTest, testWrite)
TEST_ADD(LJPairTest, testEnergy1)
TEST_ADD(LJPairTest, testForceOverR)
TEST_ADD(LJPairTest, testGetSet)
TEST_END(LJPairTest)


#endif
