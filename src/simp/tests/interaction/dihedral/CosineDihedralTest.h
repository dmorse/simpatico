#ifndef COSINE_DIHEDRAL_TEST_H
#define COSINE_DIHEDRAL_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include "DihedralTestTemplate.h"
#include <simp/interaction/dihedral/CosineDihedral.h>

// #include <util/random/Random.h>
// #include <util/containers/RArray.h>

// #include <iostream>
// #include <fstream>

using namespace Util;
using namespace Simp;

class CosineDihedralTest : public DihedralTestTemplate<CosineDihedral>
{

protected:

   using DihedralTestTemplate<CosineDihedral>::readParamFile;

public:

   void setUp() {
      eps_ = 1.0E-6;
      setNDihedralType(1);
      readParamFile("in/CosineDihedral");
   }

   void testSetUp() 
   {  
      printMethod(TEST_FUNC); 
      if (verbose() > 0) { 
         interaction_.writeParam(std::cout);
      }
   }

   void testEnergy() 
   {
      printMethod(TEST_FUNC); 
      double energy;

      b1_ = Vector(1.0, 0.0,  0.0);
      b2_ = Vector(0.0, 1.0,  0.0);
      b3_ = Vector(0.0, 0.0,  1.0);
      type_ = 0;
      energy = interaction_.energy(b1_, b2_, b3_, type_);
      TEST_ASSERT(eq(energy, 1.0));
      // std::cout << std::endl;
      // std::cout << energy << std::endl;

      b1_ = Vector( 1.1,  0.2, -0.3);
      b2_ = Vector( 0.1,  0.9,  0.2);
      b3_ = Vector(-0.1,  0.2,  1.0);
      type_ = 0;
      energy = interaction_.energy(b1_, b2_, b3_, type_);
      // std::cout << energy << std::endl;

   }

   void testForce() 
   {
      printMethod(TEST_FUNC); 

      b1_ = Vector(1.0, 0.0,  0.0);
      b2_ = Vector(0.0, 1.0,  0.0);
      b3_ = Vector(0.0, 0.0,  1.0);
      type_ = 0;
      forceTest();

      b1_ = Vector( 1.1,  0.2, -0.3);
      b2_ = Vector( 0.1,  0.9,  0.2);
      b3_ = Vector(-0.1,  0.2,  1.0);
      type_ = 0;
      forceTest();
   }



};

TEST_BEGIN(CosineDihedralTest)
TEST_ADD(CosineDihedralTest, testSetUp)
TEST_ADD(CosineDihedralTest, testEnergy)
TEST_ADD(CosineDihedralTest, testForce)
TEST_END(CosineDihedralTest)

#endif
