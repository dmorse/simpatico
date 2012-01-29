#ifndef MD_INTEGRATOR_TEST_H
#define MD_INTEGRATOR_TEST_H

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

#include <util/format/Dbl.h>
#include <mcMd/mdSimulation/MdSimulation.h>
#include <mcMd/mdSimulation/MdSystem.h>
#include <mcMd/mdIntegrators/MdIntegrator.h>
#include <mcMd/mdIntegrators/NVEIntegrator.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace McMd;

class MdIntegratorTest : public ParamFileTest<MdSimulation>
{

public:

   MdIntegratorTest()
    : ParamFileTest<MdSimulation>(),
      systemPtr_(0)
   {}

   virtual void setUp()
   {
      openFile("in/MdIntegrator");
      object().readParam(file());
      systemPtr_ = &object().system();
   }

   void testReadParam();
   void testIntegrator();

private:

   MdSystem*     systemPtr_;

};



void MdIntegratorTest::testReadParam()
{
   printMethod(TEST_FUNC);
   using std::cout;

   try {
      object().isValid();
   } catch (Exception e) {
      std::cout << e.message();
      TEST_ASSERT(0);
   }

   if (verbose() > 1) {
      object().writeParam(std::cout);
   }

}

void MdIntegratorTest::testIntegrator()
{
   printMethod(TEST_FUNC);
   double kinetic, potential;

   /// BUG:: A configurition needs to be read or generated at this point.

   systemPtr_->buildPairList();
   systemPtr_->calculateForces();
   systemPtr_->setBoltzmannVelocities(1.0);

   for (int i=0; i < 10; ++i) {

      kinetic   = systemPtr_->kineticEnergy(); 
      potential = systemPtr_->potentialEnergy(); 
      std::cout << Dbl(kinetic, 20) << Dbl(potential, 20) 
           << Dbl(kinetic + potential, 20) << std::endl;

      systemPtr_->mdIntegrator().step();
   }

   kinetic   = systemPtr_->kineticEnergy(); 
   potential = systemPtr_->potentialEnergy(); 
   std::cout << Dbl(kinetic, 20) <<  Dbl(potential, 20) 
             << Dbl(kinetic + potential, 20) << std::endl;

}

TEST_BEGIN(MdIntegratorTest)
TEST_ADD(MdIntegratorTest, testIntegrator)
TEST_END(MdIntegratorTest)

#endif
