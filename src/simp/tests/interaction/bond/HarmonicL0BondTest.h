#ifndef HARMONIC_L0_BOND_TEST_H
#define HARMONIC_L0_BOND_TEST_H

#include <simp/interaction/bond/HarmonicL0Bond.h>
#include <simp/tests/interaction/bond/BondTestTemplate.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Simp;

class HarmonicL0BondTest : public BondTestTemplate<HarmonicL0Bond>
{

protected:

   using BondTestTemplate<HarmonicL0Bond>::setNBondType;
   using BondTestTemplate<HarmonicL0Bond>::readParamFile;
   using BondTestTemplate<HarmonicL0Bond>::forceOverR;
   using BondTestTemplate<HarmonicL0Bond>::energy;

public:

   void setUp()
   {
      eps_ = 1.0E-6;
      setNBondType(1);
      readParamFile("in/HarmonicL0Bond");
   }

   void testSetUp() 
   {
      printMethod(TEST_FUNC);
      if (verbose() > 0) {
         std::cout << std::endl; 
         interaction_.writeParam(std::cout);
      }
   }

   void testEnergy() 
   {
      printMethod(TEST_FUNC);
      double e;

      e = energy(0.0, 0);
      if (verbose() > 0) {
         std::cout << std::endl; 
         std::cout << "energy(0.0, 0) = " << e << std::endl;
      }
      TEST_ASSERT(eq(e, 0.0));

      e = energy(0.49, 0);
      if (verbose() > 0) {
         std::cout << "energy(0.49, 0) = " << e << std::endl;
      }
      TEST_ASSERT(e > 0.0);

      e = energy(1.0, 0);
      if (verbose() > 0) {
         std::cout << std::endl; 
         std::cout << "energy(1.00, 0) = " << e << std::endl;
      }

      e = energy(1.21, 0);
      if (verbose() > 0) {
         std::cout << "energy(1.25, 0) = " << e << std::endl;
      }

   }

   void testForceOverR() 
   {
      printMethod(TEST_FUNC);
      double f;
      type_ = 0;
 
      rsq_ = 1.0;
      f = forceOverR();
      TEST_ASSERT(testForce());
      if (verbose() > 0) {
         std::cout << std::endl;
         std::cout << "forceOverR(1.00, 0) = " <<  f << std::endl;
      }
      // TEST_ASSERT(eq(f, 48.0));

      rsq_ = 0.81;
      f = forceOverR();
      TEST_ASSERT(testForce());
      if (verbose() > 0) {
         std::cout << "forceOverR(0.81, 0) = " <<  f << std::endl;
      }
      // TEST_ASSERT(f > 48.0);

      rsq_ = 1.21;
      f = forceOverR();
      TEST_ASSERT(testForce());
      if (verbose() > 0) {
         std::cout << "forceOverR(1.21, 0) = " <<  f << std::endl;
      }
      // TEST_ASSERT(eq(f, 0.0));

   }

   void testGetSet() {
      printMethod(TEST_FUNC);

      TEST_ASSERT(eq(interaction_.kappa(0), 8.0));
      interaction_.set("kappa", 0, 9.0);
      TEST_ASSERT(eq(interaction_.kappa(0), 9.0));
      TEST_ASSERT(eq(interaction_.get("kappa", 0), 9.0));
   }

   void testSaveLoad() {
      printMethod(TEST_FUNC);

      Serializable::OArchive oar;
      openOutputFile("out/HarmonicL0.rst", oar.file());
      interaction_.save(oar);
      oar.file().close();
      
      Serializable::IArchive iar;
      openInputFile("out/HarmonicL0.rst", iar.file());

      HarmonicBond clone;
      clone.setNBondType(1);
      clone.load(iar);
      iar.file().close();

      TEST_ASSERT(eq(interaction_.kappa(0), clone.kappa(0)));

      TEST_ASSERT( eq(interaction_.energy(0.95, 0), 
                      clone.energy(0.95, 0) ));
      TEST_ASSERT( eq(interaction_.forceOverR(0.95, 0), 
                      clone.forceOverR(0.95, 0) ));
   }

   void testRandomBondLength() 
   {
      int    type, i;
      double r, beta;
  
      Random *random;
 
      printMethod(TEST_FUNC);
 
      std::ifstream in;
      openInputFile("in/Random", in);
      random = new Random;
      random->readParam(in);

      beta = 1.0;
      type = 0;
      for (i=0; i < 20; i++) {
         r = interaction_.randomBondLength(random, beta, type);
         std::cout << "random bond length = " << r << std::endl;
      }
   }

};

TEST_BEGIN(HarmonicL0BondTest)
TEST_ADD(HarmonicL0BondTest, testSetUp)
TEST_ADD(HarmonicL0BondTest, testEnergy)
TEST_ADD(HarmonicL0BondTest, testForceOverR)
TEST_ADD(HarmonicL0BondTest, testGetSet)
TEST_ADD(HarmonicL0BondTest, testSaveLoad)
//TEST_ADD(HarmonicL0BondTest, testRandomBondLength)
TEST_END(HarmonicL0BondTest)

#endif
