#ifndef MULTI_HARMONIC_DIHEDRAL_TEST_H
#define MULTI_HARMONIC_DIHEDRAL_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include "DihedralTestTemplate.h"
#include <inter/dihedral/MultiHarmonicDihedral.h>

// #include <util/random/Random.h>
// #include <util/containers/RArray.h>

// #include <iostream>
// #include <fstream>

using namespace Util;
using namespace Inter;

class MultiHarmonicDihedralTest : public DihedralTestTemplate<MultiHarmonicDihedral>
{

protected:

   using DihedralTestTemplate<MultiHarmonicDihedral>::readParamFile;

public:

   void setUp() {
      eps_ = 1.0E-6;
      setNDihedralType(1);
      readParamFile("in/MultiHarmonicDihedral");
   }

   void testSetUp() 
   {  
      printMethod(TEST_FUNC); 

      std::cout << std::endl;
      interaction_.writeParam(std::cout);

      if (verbose() > 0) { 
         interaction_.writeParam(std::cout);
      }
   }

   void testEnergy() 
   {
      printMethod(TEST_FUNC); 
      Torsion torsion;
      double energy, phi;

      b1_ = Vector(1.0, 0.0,  0.0);
      b2_ = Vector(0.0, 1.0,  0.0);
      b3_ = Vector(0.0, 0.0,  1.0);
      type_ = 0;
      energy = interaction_.energy(b1_, b2_, b3_, type_);
      TEST_ASSERT(eq(energy, 1.0));
      torsion.computeAngle(b1_, b2_, b3_);
      phi = torsion.phi();
      TEST_ASSERT(eq(energy, 1.0 + cos(3.0*phi)));

      // std::cout << std::endl;
      // std::cout << energy << std::endl;

      b1_ = Vector( 1.1,  0.2, -0.3);
      b2_ = Vector( 0.1,  0.9,  0.2);
      b3_ = Vector(-0.1,  0.2,  1.0);
      type_ = 0;
      energy = interaction_.energy(b1_, b2_, b3_, type_);
      torsion.computeAngle(b1_, b2_, b3_);
      phi = torsion.phi();
      TEST_ASSERT(eq(energy, 1.0 + cos(3.0*phi)));
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

TEST_BEGIN(MultiHarmonicDihedralTest)
TEST_ADD(MultiHarmonicDihedralTest, testSetUp)
TEST_ADD(MultiHarmonicDihedralTest, testEnergy)
TEST_ADD(MultiHarmonicDihedralTest, testForce)
TEST_END(MultiHarmonicDihedralTest)

#endif
