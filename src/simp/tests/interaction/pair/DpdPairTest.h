#ifndef DPD_PAIR_TEST_H
#define DPD_PAIR_TEST_H

#include <simp/interaction/pair/DpdPair.h>
#include <simp/tests/interaction/pair/PairTestTemplate.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Simp;

class DpdPairTest : public PairTestTemplate<DpdPair>
{

protected:

   using PairTestTemplate<DpdPair>::setNAtomType;
   using PairTestTemplate<DpdPair>::readParamFile;
   using PairTestTemplate<DpdPair>::forceOverR;
   using PairTestTemplate<DpdPair>::energy;

public:

   void setUp()
   {
      eps_ = 1.0E-6;
      setNAtomType(2);
      readParamFile("in/DpdPair");
      // setVerbose(1);
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
      TEST_ASSERT(eq(e, 0.0));

      e = energy(0.36, 0, 1);
      if (verbose() > 0) {
         std::cout << "energy(0.6, 0, 1) = " << e << std::endl;
      }
      TEST_ASSERT(eq(e, 0.16));

      e = energy(1.21, 0, 1);
      if (verbose() > 0) {
         std::cout << "energy(1.1, 0, 1) = " << e << std::endl;
      }
      TEST_ASSERT(eq(e, 0.0));

   }

   void testForceOverR() 
   {
      printMethod(TEST_FUNC);
 
      type1_ = 0;
      type2_ = 1;
      double f;

      double r = 0.5;
      rsq_ = r*r;
      f = forceOverR();
      TEST_ASSERT(testForce());
      if (verbose() > 0) {
         std::cout << std::endl;
         std::cout << "forceOverR(0.25, 0, 1) = " <<  f << std::endl;
      }
      TEST_ASSERT(eq(f*r, 1.0));

      r = 0.8;
      rsq_ = r*r;
      f = forceOverR();
      TEST_ASSERT(testForce());
      if (verbose() > 0) {
         std::cout << "forceOverR(0.64, 0, 1) = " <<  f << std::endl;
      }
      TEST_ASSERT(eq(f*r, 0.4));

      type1_ = 1;
      type2_ = 1;
      f = forceOverR();
      TEST_ASSERT(testForce());
      if (verbose() > 0) {
         std::cout << "forceOverR(0.81, 1, 1) = " <<  f << std::endl;
      }
      TEST_ASSERT(eq(f*r, 0.2));


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
      openOutputFile("out/serial", oar.file());
      interaction_.save(oar);
      oar.file().close();

      Serializable::IArchive iar;
      openInputFile("out/serial", iar.file());

      DpdPair clone;
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

TEST_BEGIN(DpdPairTest)
TEST_ADD(DpdPairTest, testSetUp)
TEST_ADD(DpdPairTest, testEnergy)
TEST_ADD(DpdPairTest, testForceOverR)
TEST_ADD(DpdPairTest, testGetSet)
TEST_ADD(DpdPairTest, testModify)
TEST_ADD(DpdPairTest, testSaveLoad)
TEST_END(DpdPairTest)

#endif
