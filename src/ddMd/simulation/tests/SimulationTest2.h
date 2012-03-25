#ifndef DDMD_SYSTEM_TEST_H
#define DDMD_SYSTEM_TEST_H

#include <ddMd/system/System.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/storage/GhostIterator.h>
#include <ddMd/interaction/Interaction.h>
#include <util/random/Random.h>
#include <util/format/Dbl.h>
#include <util/mpi/MpiLogger.h>

#ifdef UTIL_MPI
#define TEST_MPI
#endif

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace DdMd;

class SystemTest : public ParamFileTest<System>
{

public:

   virtual void setUp()
   {}

   void testIntegrate();

};

inline void SystemTest::testIntegrate()
{
   printMethod(TEST_FUNC); 

   openFile("in/param2"); 
   object().readParam(file()); 

   Domain&  domain  = object().domain();
   AtomStorage& storage = object().atomStorage();
   int myRank = domain.gridRank();

   std::string filename("in/config2");
   object().readConfig(filename);
   object().exchanger().exchangeGhosts();

   //object().interaction().setMethodId(1);
   object().interaction().findNeighbors();

   double temperature = 1.0;
   object().setBoltzmannVelocities(temperature);

   // Exchange ghosts among processors.

   // Calculate energies before integration
   double kinetic   = object().kineticEnergy();
   double potential = object().pairPotentialEnergy();
   if (myRank == 0) {
      std::cout << Dbl(kinetic) << Dbl(potential) 
                << Dbl(kinetic + potential) << std::endl;
   }


   for (int i = 0; i < 5; ++i ) {

      object().integrate(10000);
   
      // Calculate energies after integration
      kinetic   = object().kineticEnergy();
      potential = object().pairPotentialEnergy();
      if (myRank == 0) {
         std::cout << Dbl(kinetic) << Dbl(potential) 
                   << Dbl(kinetic + potential) << std::endl;
      }

      int nAtomAll = object().nAtomTotal();
      //if (myRank == 0) {
      //   TEST_ASSERT(nAtomAll == 100);
      //}
      //TEST_ASSERT(object().isValid());

   }

}


TEST_BEGIN(SystemTest)
TEST_ADD(SystemTest, testIntegrate)
TEST_END(SystemTest)

#endif
