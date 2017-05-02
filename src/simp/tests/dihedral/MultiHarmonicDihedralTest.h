#ifndef MULTI_HARMONIC_DIHEDRAL_TEST_H
#define MULTI_HARMONIC_DIHEDRAL_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include "DihedralTestTemplate.h"
#include <simp/interaction/dihedral/MultiHarmonicDihedral.h>

using namespace Util;
using namespace Simp;

class MultiHarmonicDihedralTest 
 : public DihedralTestTemplate<MultiHarmonicDihedral>
{

protected:

   using DihedralTestTemplate<MultiHarmonicDihedral>::setNDihedralType;
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

   void testGet() 
   {
      printMethod(TEST_FUNC); 
      double k0, k1, k2, k3, k4;
      k0 = interaction_.get("k0", type_);
      k1 = interaction_.get("k1", type_);
      k2 = interaction_.get("k2", type_);
      k3 = interaction_.get("k3", type_);
      k4 = interaction_.get("k4", type_);
      TEST_ASSERT(eq(k0, 1.0));
      TEST_ASSERT(eq(k1, 0.0));
      TEST_ASSERT(eq(k2, 0.0));
      TEST_ASSERT(eq(k3, 1.0));
      TEST_ASSERT(eq(k4, 0.0));
   }

   void testSetGet() 
   {
      printMethod(TEST_FUNC); 
      double k0, k1, k2, k3, k4;
      k0 = interaction_.get("k0", type_);
      k1 = interaction_.get("k1", type_);
      k2 = interaction_.get("k2", type_);
      k3 = interaction_.get("k3", type_);
      k4 = interaction_.get("k4", type_);
      TEST_ASSERT(eq(k0, 1.0));
      TEST_ASSERT(eq(k1, 0.0));
      TEST_ASSERT(eq(k2, 0.0));
      TEST_ASSERT(eq(k3, 1.0));
      TEST_ASSERT(eq(k4, 0.0));

      interaction_.set("k0", type_, 0.9);
      interaction_.set("k1", type_, 0.5);
      interaction_.set("k2", type_, 0.2);
      interaction_.set("k4", type_, 0.1);
      k0 = interaction_.get("k0", type_);
      k1 = interaction_.get("k1", type_);
      k2 = interaction_.get("k2", type_);
      k3 = interaction_.get("k3", type_);
      k4 = interaction_.get("k4", type_);
      TEST_ASSERT(eq(k0, 0.9));
      TEST_ASSERT(eq(k1, 0.5));
      TEST_ASSERT(eq(k2, 0.2));
      TEST_ASSERT(eq(k3, 1.0));
      TEST_ASSERT(eq(k4, 0.1));
   }

   void testEnergy() 
   {
      printMethod(TEST_FUNC); 

      b1_ = Vector(1.0, 0.0,  0.0);
      b2_ = Vector(0.0, 1.0,  0.0);
      b3_ = Vector(0.0, 0.0,  1.0);
      type_ = 0;
      energyTest();
      double energy = interaction_.energy(b1_, b2_, b3_, type_);
      TEST_ASSERT(eq(energy, 1.0));
      TEST_ASSERT(energyTest());

      b1_ = Vector( 1.1,  0.2, -0.3);
      b2_ = Vector( 0.1,  0.9,  0.2);
      b3_ = Vector(-0.1,  0.2,  1.0);
      type_ = 0;
      TEST_ASSERT(energyTest());

      interaction_.set("k0", type_, 0.9);
      interaction_.set("k1", type_, 0.5);
      interaction_.set("k2", type_, 0.2);
      interaction_.set("k4", type_, 0.1);
      TEST_ASSERT(energyTest());
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

   bool energyTest()
   {
      Torsion torsion;
      double e, fe, phi;
      double k0, k1, k2, k3, k4;

      torsion.computeAngle(b1_, b2_, b3_);
      phi = torsion.phi();

      k0 = interaction_.get("k0", type_);
      k1 = interaction_.get("k1", type_);
      k2 = interaction_.get("k2", type_);
      k3 = interaction_.get("k3", type_);
      k4 = interaction_.get("k4", type_);
      fe =  k0 + k1*cos(phi) + k2*cos(2.0*phi);
      fe += k3*cos(3.0*phi) + k4*cos(4.0*phi);

      e = interaction_.energy(b1_, b2_, b3_, type_);
      return (eq(e, fe));
   }

};

TEST_BEGIN(MultiHarmonicDihedralTest)
TEST_ADD(MultiHarmonicDihedralTest, testSetUp)
TEST_ADD(MultiHarmonicDihedralTest, testGet)
TEST_ADD(MultiHarmonicDihedralTest, testSetGet)
TEST_ADD(MultiHarmonicDihedralTest, testEnergy)
TEST_ADD(MultiHarmonicDihedralTest, testForce)
TEST_END(MultiHarmonicDihedralTest)

#endif
