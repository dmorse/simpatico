#ifndef LJ_PAIR_TEST_H
#define LJ_PAIR_TEST_H

#include <inter/pair/LJPair.h>
#include <inter/pair/PairTest.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Inter;

class LJPairTest : public PairTest<LJPair>
{

public:

   void setUp()
   {
      nAtomType_ = 2;
      interaction_.setNAtomType(nAtomType_);
      readParamFile("in/LJPair");
   }

   void tearDown(){}

   void testSetUp() {
      printMethod(TEST_FUNC);
   }

   void testWrite() {

      printMethod(TEST_FUNC);

      // Verbose output
      if (verbose() > 0) {
         std::cout << std::endl; 
         interaction_.writeParam(std::cout);
      }

   }

   void testEnergy1() {
     
      printMethod(TEST_FUNC);
      std::cout << std::endl; 

      ParameterSet parameter(interaction_);
      parameter.set(1.00, 0, 1);
      std::cout << "interaction_.energy(1.00, 0, 1) = " << parameter.energy() << std::endl;

      parameter.set(0.81, 0, 1);
      std::cout << "interaction_.energy(0.81, 0, 1) = " << parameter.energy() << std::endl;

      parameter.set(1.25992, 0, 1);
      std::cout << "interaction_.energy(1.15992, 0, 1) = " << parameter.energy() << std::endl;

      parameter.set(1.5, 0, 1);
      std::cout << "interaction_.energy(1.50, 0, 1) = " << parameter.energy() << std::endl;
   }


   void testForceOverR() {
      printMethod(TEST_FUNC);
     
      ParameterSet parameter(interaction_);
      parameter.set(1.00, 0, 1);
      TEST_ASSERT(parameter.testForce());
      std::cout << "interaction_.forceOverR(1.00, 0, 1) = " << parameter.forceOverR() << std::endl;

      parameter.set(0.81, 0, 1);
      TEST_ASSERT(parameter.testForce());
      std::cout << "interaction_.forceOverR(0.81, 0, 1) = " << parameter.forceOverR() << std::endl;

      parameter.set(1.25990, 0, 1);
      TEST_ASSERT(parameter.testForce());
      std::cout << "interaction_.forceOverR(1.25990, 0, 1) = " << parameter.forceOverR() << std::endl;

      parameter.set(1.5, 0, 1);
      TEST_ASSERT(parameter.testForce());
      std::cout << "interaction_.forceOverR(1.50, 0, 1) = " << parameter.forceOverR() << std::endl;
   }

   void testGetSet() {
      printMethod(TEST_FUNC);

      TEST_ASSERT(eq(interaction_.epsilon(0, 0), 1.0));
      TEST_ASSERT(eq(interaction_.epsilon(1, 1), 1.0));
      TEST_ASSERT(eq(interaction_.epsilon(0, 1), 2.0));
      TEST_ASSERT(eq(interaction_.epsilon(1, 0), 2.0));
      TEST_ASSERT(eq(interaction_.sigma(0, 0), 1.0));
      TEST_ASSERT(eq(interaction_.sigma(1, 1), 1.0));
      TEST_ASSERT(eq(interaction_.sigma(0, 1), 1.0));
      TEST_ASSERT(eq(interaction_.sigma(1, 0), 1.0));

      interaction_.setEpsilon(0, 1, 1.3);
      TEST_ASSERT(eq(interaction_.epsilon(0, 0), 1.0));
      TEST_ASSERT(eq(interaction_.epsilon(1, 1), 1.0));
      TEST_ASSERT(eq(interaction_.epsilon(0, 1), 1.3));
      TEST_ASSERT(eq(interaction_.epsilon(1, 0), 1.3));
      
      interaction_.setEpsilon(0, 0, 1.1);
      TEST_ASSERT(eq(interaction_.epsilon(0, 0), 1.1));
      TEST_ASSERT(eq(interaction_.epsilon(1, 1), 1.0));
      TEST_ASSERT(eq(interaction_.epsilon(0, 1), 1.3));
      TEST_ASSERT(eq(interaction_.epsilon(1, 0), 1.3));
      
      interaction_.setSigma(0, 1, 1.05);
      TEST_ASSERT(eq(interaction_.sigma(0, 0), 1.0));
      TEST_ASSERT(eq(interaction_.sigma(1, 1), 1.0));
      TEST_ASSERT(eq(interaction_.sigma(0, 1), 1.05));
      TEST_ASSERT(eq(interaction_.sigma(1, 0), 1.05));

   }

   void testModify() {
      printMethod(TEST_FUNC);

      TEST_ASSERT(eq(interaction_.epsilon(0, 0), 1.0));
      TEST_ASSERT(eq(interaction_.epsilon(1, 1), 1.0));
      TEST_ASSERT(eq(interaction_.epsilon(0, 1), 2.0));
      TEST_ASSERT(eq(interaction_.epsilon(1, 0), 2.0));
      TEST_ASSERT(eq(interaction_.sigma(0, 0), 1.0));
      TEST_ASSERT(eq(interaction_.sigma(1, 1), 1.0));
      TEST_ASSERT(eq(interaction_.sigma(0, 1), 1.0));
      TEST_ASSERT(eq(interaction_.sigma(1, 0), 1.0));

      interaction_.set("epsilon", 0, 1, 1.3);
      TEST_ASSERT(eq(interaction_.epsilon(0, 0), 1.0));
      TEST_ASSERT(eq(interaction_.epsilon(1, 1), 1.0));
      TEST_ASSERT(eq(interaction_.epsilon(0, 1), 1.3));
      TEST_ASSERT(eq(interaction_.epsilon(1, 0), 1.3));
      
      interaction_.set("epsilon", 0, 0, 1.1);
      TEST_ASSERT(eq(interaction_.epsilon(0, 0), 1.1));
      TEST_ASSERT(eq(interaction_.epsilon(1, 1), 1.0));
      TEST_ASSERT(eq(interaction_.epsilon(0, 1), 1.3));
      TEST_ASSERT(eq(interaction_.epsilon(1, 0), 1.3));
      TEST_ASSERT(eq(interaction_.sigma(0, 0), 1.0));
      TEST_ASSERT(eq(interaction_.sigma(1, 1), 1.0));
      TEST_ASSERT(eq(interaction_.sigma(0, 1), 1.0));
      TEST_ASSERT(eq(interaction_.sigma(1, 0), 1.0));
      
      interaction_.set("epsilon", 0, 1, 1.05);
      TEST_ASSERT(eq(interaction_.epsilon(0, 0), 1.1));
      TEST_ASSERT(eq(interaction_.epsilon(1, 1), 1.0));
      TEST_ASSERT(eq(interaction_.epsilon(0, 1), 1.05));
      TEST_ASSERT(eq(interaction_.epsilon(1, 0), 1.05));
      TEST_ASSERT(eq(interaction_.sigma(0, 0), 1.0));
      TEST_ASSERT(eq(interaction_.sigma(1, 1), 1.0));
      TEST_ASSERT(eq(interaction_.sigma(0, 1), 1.0));
      TEST_ASSERT(eq(interaction_.sigma(1, 0), 1.0));

      interaction_.set("sigma", 0, 1, 1.05);
      TEST_ASSERT(eq(interaction_.epsilon(0, 0), 1.1));
      TEST_ASSERT(eq(interaction_.epsilon(1, 1), 1.0));
      TEST_ASSERT(eq(interaction_.epsilon(0, 1), 1.05));
      TEST_ASSERT(eq(interaction_.epsilon(1, 0), 1.05));
      TEST_ASSERT(eq(interaction_.sigma(0, 0), 1.0));
      TEST_ASSERT(eq(interaction_.sigma(1, 1), 1.0));
      TEST_ASSERT(eq(interaction_.sigma(0, 1), 1.05));
      TEST_ASSERT(eq(interaction_.sigma(1, 0), 1.05));

      TEST_ASSERT(eq(interaction_.epsilon(0, 0), interaction_.get("epsilon", 0, 0)));
      TEST_ASSERT(eq(interaction_.epsilon(1, 1), interaction_.get("epsilon", 1, 1)));
      TEST_ASSERT(eq(interaction_.epsilon(1, 0), interaction_.get("epsilon", 1, 0)));
      TEST_ASSERT(eq(interaction_.epsilon(0, 1), interaction_.get("epsilon", 0, 1)));
      TEST_ASSERT(eq(interaction_.sigma(0, 0), interaction_.get("sigma", 0, 0)));
      TEST_ASSERT(eq(interaction_.sigma(1, 1), interaction_.get("sigma", 1, 1)));
      TEST_ASSERT(eq(interaction_.sigma(1, 0), interaction_.get("sigma", 1, 0)));
      TEST_ASSERT(eq(interaction_.sigma(0, 1), interaction_.get("sigma", 0, 1)));

   }

   void testSaveLoad() {
      printMethod(TEST_FUNC);

      Serializable::OArchive oar;
      openOutputFile("out/serial", oar.file());
      interaction_.save(oar);
      oar.file().close();
      
      Serializable::IArchive iar;
      openInputFile("out/serial", iar.file());

      LJPair clone;
      clone.load(iar);

      TEST_ASSERT(eq(interaction_.epsilon(0, 0), clone.epsilon(0,0)));
      TEST_ASSERT(eq(interaction_.epsilon(1, 0), clone.epsilon(1,0)));
      TEST_ASSERT(eq(interaction_.epsilon(0, 1), clone.epsilon(0,1)));
      TEST_ASSERT(eq(interaction_.epsilon(1, 1), clone.epsilon(1,1)));
      TEST_ASSERT(eq(interaction_.sigma(0, 0), clone.sigma(0,0)));
      TEST_ASSERT(eq(interaction_.sigma(1, 0), clone.sigma(1,0)));
      TEST_ASSERT(eq(interaction_.sigma(0, 1), clone.sigma(0,1)));
      TEST_ASSERT(eq(interaction_.sigma(1, 1), clone.sigma(1,1)));

      TEST_ASSERT(eq(interaction_.energy(0.95, 0, 1), clone.energy(0.95, 0, 1)));
      TEST_ASSERT(eq(interaction_.forceOverR(0.95, 0, 1), clone.forceOverR(0.95, 0, 1)));
      TEST_ASSERT(eq(interaction_.energy(1.25, 0, 1), clone.energy(1.25, 0, 1)));
      TEST_ASSERT(eq(interaction_.forceOverR(1.25, 0, 1), clone.forceOverR(1.25, 0, 1)));
      TEST_ASSERT(eq(interaction_.energy(1.26, 0, 1), clone.energy(1.26, 0, 1)));
      TEST_ASSERT(eq(interaction_.forceOverR(1.26, 0, 1), clone.forceOverR(1.26, 0, 1)));
   }

};

TEST_BEGIN(LJPairTest)
TEST_ADD(LJPairTest, testSetUp)
TEST_ADD(LJPairTest, testWrite)
TEST_ADD(LJPairTest, testEnergy1)
TEST_ADD(LJPairTest, testForceOverR)
TEST_ADD(LJPairTest, testGetSet)
TEST_ADD(LJPairTest, testModify)
TEST_ADD(LJPairTest, testSaveLoad)
TEST_END(LJPairTest)

#endif
