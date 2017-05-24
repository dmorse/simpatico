#ifndef LJ_PAIR_TEST_H
#define LJ_PAIR_TEST_H

#include <simp/interaction/pair/LJPair.h>
#include <simp/tests/interaction/pair/PairTestTemplate.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Simp;

class LJPairTest : public PairTestTemplate<LJPair>
{

protected:

   using PairTestTemplate<LJPair>::setNAtomType;
   using PairTestTemplate<LJPair>::readParamFile;
   using PairTestTemplate<LJPair>::forceOverR;
   using PairTestTemplate<LJPair>::energy;

public:

   void setUp()
   {
      eps_ = 1.0E-6;
      setNAtomType(2);
      readParamFile("in/LJPair");
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

      e = energy(1.0, 0, 1);
      if (verbose() > 0) {
         std::cout << std::endl; 
         std::cout << "energy(1.00, 0, 1) = " << e << std::endl;
      }
      TEST_ASSERT(eq(e, 2.0));

      e = energy(0.81, 0, 1);
      if (verbose() > 0) {
         std::cout << "energy(0.81, 0, 1) = " << e << std::endl;
      }
      TEST_ASSERT(e > 2.0);

      e = energy(1.25, 0, 1);
      if (verbose() > 0) {
         std::cout << "energy(1.25, 0, 1) = " << e << std::endl;
      }
      TEST_ASSERT(e > 0.0);

      e = energy(1.25992105, 0, 1);
      if (verbose() > 0) {
         std::cout << "energy(1.25992105, 0, 1) = " << e << std::endl;
      }
      TEST_ASSERT(eq(e, 0.0));

      e = energy(1.5, 0, 1);
      if (verbose() > 0) {
         std::cout << "energy(1.5, 0, 1) = " << e << std::endl;
      }
      TEST_ASSERT(eq(e, 0.0));
   }

   void testForceOverR() 
   {
      printMethod(TEST_FUNC);
 
      type1_ = 0;
      type2_ = 1;
      double f;
 
      rsq_ = 1.0;
      f = forceOverR();
      TEST_ASSERT(testForce());
      if (verbose() > 0) {
         std::cout << std::endl;
         std::cout << "forceOverR(1.00, 0, 1) = " <<  f << std::endl;
      }
      TEST_ASSERT(eq(f, 48.0));

      rsq_ = 0.81;
      f = forceOverR();
      TEST_ASSERT(testForce());
      if (verbose() > 0) {
         std::cout << "forceOverR(0.81, 0, 1) = " <<  f << std::endl;
      }
      TEST_ASSERT(f > 48.0);

      rsq_ = 1.25992105;
      f = forceOverR();
      TEST_ASSERT(testForce());
      if (verbose() > 0) {
         std::cout << "forceOverR(1.25992105, 0, 1) = " <<  f << std::endl;
      }
      TEST_ASSERT(std::abs(f) < eps_);

      //Note: Do not test beyond cutoff: result is undefined.
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
      openOutputFile("out/LJPair.rst", oar.file());
      interaction_.save(oar);
      oar.file().close();

      Serializable::IArchive iar;
      openInputFile("out/LJPair.rst", iar.file());

      LJPair clone;
      clone.setNAtomType(2);
      clone.loadParameters(iar);

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
TEST_ADD(LJPairTest, testEnergy)
TEST_ADD(LJPairTest, testForceOverR)
TEST_ADD(LJPairTest, testGetSet)
TEST_ADD(LJPairTest, testModify)
TEST_ADD(LJPairTest, testSaveLoad)
TEST_END(LJPairTest)

#endif
