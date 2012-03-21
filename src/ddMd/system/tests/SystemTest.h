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
   {  
   }

   void testReadParam();

   void testReadConfig();

   void testExchangeAtoms();

   void testExchange();

   void testUpdate();

   void testCalculateForces();

   void testIntegrate1();

};

inline void SystemTest::testReadParam()
{  
   printMethod(TEST_FUNC); 

   openFile("in/param1"); 
   object().readParam(file()); 
   if (verbose() > 0) {
      object().writeParam(std::cout);
   }
   // TEST_ASSERT(object().buffer().isAllocated())
   TEST_ASSERT(object().domain().isInitialized());
   TEST_ASSERT(object().domain().hasBoundary());
}

inline void SystemTest::testReadConfig()
{
   printMethod(TEST_FUNC); 

   openFile("in/param1"); 
   object().readParam(file()); 

   Domain&  domain  = object().domain();
   AtomStorage& storage = object().atomStorage();

   std::string filename("in/config1");
   object().readConfig(filename);

   int nAtomAll = object().nAtomTotal();
   int myRank = domain.gridRank();
   if (myRank == 0) {
      //std::cout << "Total atom count = " << nAtomAll << std::endl;
      TEST_ASSERT(nAtomAll == 100);
   }
   std::cout.flush();

   #if 0
   MpiLogger logger;
   logger.begin();

   std::cout << std::endl;
   std::cout << "Processor " << myRank << std::endl;
   std::cout << "Position  ";
   for (int i = 0; i < Dimension; ++i) {
       std::cout << domain.gridCoordinate(i) << "  ";
   }
   std::cout << std::endl;

   for (int i = 0; i < Dimension; ++i) {
      std::cout << "bound(" << i << ", 0) = " 
                << Dbl(domain.domainBound(i, 0));
      std::cout << "  bound(" << i << ", 1) = " 
                << Dbl(domain.domainBound(i, 1));
      std::cout << std::endl;
   }
   std::cout << std::endl;
   std::cout << "Lengths " << object().boundary().lengths()
             << std::endl;
   std::cout << "nAtom  " << storage.nAtom() << std::endl;
   std::cout << "nGhost " << storage.nGhost()<< std::endl;
   std::cout << std::endl;
   logger.end();
   #endif

   // Check that all atoms are within the processor domain.
   int j = 0;
   AtomIterator atomIter;
   storage.begin(atomIter);
   for ( ; !atomIter.atEnd(); ++atomIter) {
      j++;
      TEST_ASSERT( domain.isInDomain( atomIter->position() ) );
   }
   TEST_ASSERT(j == storage.nAtom());

}

inline void SystemTest::testExchangeAtoms()
{
   printMethod(TEST_FUNC); 

   openFile("in/param1"); 
   object().readParam(file()); 

   Domain&  domain  = object().domain();
   AtomStorage& storage = object().atomStorage();
   Random&  random = object().random();

   std::string filename("in/config1");
   object().readConfig(filename);

   // Range of random increments for the atom positions.
   double range1 = double(-0.8);
   double range2 = double(0.8);

   // Add a random increment to atom positions
   AtomIterator atomIter;
   for (int i = 0; i < Dimension; ++i) {
     storage.begin(atomIter);
     for ( ; !atomIter.atEnd(); ++atomIter) {
        atomIter->position()[i] += random.uniform(range1, range2);
     }
   }

   // Exchange atoms among processors
   object().exchanger().exchangeAtoms();

   // Check that all atoms are accounted for after exchange.
   int myRank   = domain.gridRank();
   int nAtomAll = object().nAtomTotal();
   if (myRank == 0) {
      // std::cout << "Total atom count = " << nAtomAll << std::endl;
      TEST_ASSERT(nAtomAll == 100);
   }

   // Check that all atoms are within the processor domain.
   int j = 0;
   storage.begin(atomIter);
   for ( ; !atomIter.atEnd(); ++atomIter) {
      j++;
      TEST_ASSERT( domain.isInDomain( atomIter->position() ) );
   }
   TEST_ASSERT(j == storage.nAtom());

}

inline void SystemTest::testExchange()
{
   printMethod(TEST_FUNC); 

   openFile("in/param1"); 
   object().readParam(file()); 

   Domain&  domain  = object().domain();
   AtomStorage& storage = object().atomStorage();
   Random&  random = object().random();

   std::string filename("in/config1");
   object().readConfig(filename);

   // Range of random increments for the atom positions.
   double range1 = double(-0.8);
   double range2 = double(0.8);

   // Add a random increment to atom positions
   AtomIterator atomIter;
   for (int i = 0; i < Dimension; ++i) {
     storage.begin(atomIter);
     for ( ; !atomIter.atEnd(); ++atomIter) {
        atomIter->position()[i] += random.uniform(range1, range2);
     }
   }

   // Exchange atoms among processors
   object().exchanger().exchange();

   // Check that all atoms are accounted for after exchange.
   int myRank   = domain.gridRank();
   int nAtomAll = object().nAtomTotal();
   if (myRank == 0) {
      TEST_ASSERT(nAtomAll == 100);
   }

   // Check that all atoms are within the processor domain.
   int j = 0;
   storage.begin(atomIter);
   for ( ; !atomIter.atEnd(); ++atomIter) {
      j++;
      TEST_ASSERT( domain.isInDomain( atomIter->position() ) );
   }
   TEST_ASSERT(j == storage.nAtom());
   int nAtom = storage.nAtom();

   // Check that all ghosts are outside the processor domain.
   GhostIterator ghostIter;
   storage.begin(ghostIter);
   for ( ; !ghostIter.atEnd(); ++ghostIter) {
      TEST_ASSERT(!domain.isInDomain(ghostIter->position()));
   }

   int nGhostAll = object().nGhostTotal();
   if (myRank == 0) {
      std::cout << "Total ghost count = " << nGhostAll << std::endl;
   }

}

inline void SystemTest::testUpdate()
{
   printMethod(TEST_FUNC); 

   openFile("in/param1"); 
   object().readParam(file()); 

   Domain&  domain  = object().domain();
   AtomStorage& storage = object().atomStorage();
   Random&  random  = object().random();

   std::string filename("in/config1");
   object().readConfig(filename);

   // Add a random increment to atom positions
   double range1 = double(-0.8);
   double range2 = double(0.8);
   AtomIterator atomIter;
   for (int i = 0; i < Dimension; ++i) {
     storage.begin(atomIter);
     for ( ; !atomIter.atEnd(); ++atomIter) {
        atomIter->position()[i] += random.uniform(range1, range2);
     }
   }

   // Exchange atoms among processors
   object().exchanger().exchange();

   // Check that all atoms are accounted for after exchange.
   int myRank   = domain.gridRank();
   int nAtomAll = object().nAtomTotal();
   if (myRank == 0) {
      TEST_ASSERT(nAtomAll == 100);
   }

   // Check that all atoms are within the processor domain.
   int j = 0;
   storage.begin(atomIter);
   for ( ; !atomIter.atEnd(); ++atomIter) {
      j++;
      TEST_ASSERT( domain.isInDomain( atomIter->position() ) );
   }
   TEST_ASSERT(j == storage.nAtom());
   int nAtom = storage.nAtom();

   // Check that all ghosts are outside the processor domain.
   GhostIterator ghostIter;
   storage.begin(ghostIter);
   for ( ; !ghostIter.atEnd(); ++ghostIter) {
      TEST_ASSERT(!domain.isInDomain(ghostIter->position()));
   }

   int nGhostAll = object().nGhostTotal();
   if (myRank == 0) {
      std::cout << "Total ghost count = " << nGhostAll << std::endl;
   }

   #if 0
   //Print number of atoms on each processor after the ghost exchange.
   MpiLogger logger;
   logger.begin();
   std::cout << "Processor " << myRank 
             << " : Post-ghost exchange Atoms count = "
             << storage.nAtom() << std::endl;
   logger.end();

   // Print number of ghosts on each processor after the exchange.
   logger.begin();
   std::cout << "Processor " << myRank 
             << " : Post-ghost exchange Ghost count = "
             << storage.nGhost() << std::endl;
   logger.end();
   #endif

   // Add a random increment to atom positions
   range1 = double(-0.1);
   range2 = double(+0.1);
   storage.begin(atomIter);
   for ( ; !atomIter.atEnd(); ++atomIter) {
      for (int i = 0; i < Dimension; ++i) {
         atomIter->position()[i] += random.uniform(range1, range2);
      }
   }

   object().exchanger().update();

}

inline void SystemTest::testCalculateForces()
{
   printMethod(TEST_FUNC); 

   openFile("in/param1"); 
   object().readParam(file()); 

   Domain&  domain  = object().domain();
   AtomStorage& storage = object().atomStorage();
   int myRank = domain.gridRank();

   std::string filename("in/config1");
   object().readConfig(filename);

   // Exchange ghosts among processsors
   object().exchanger().exchange();
   object().interaction().findNeighbors();

   object().interaction().calculateForces();

   Vector f(0.0); // total force on processor
   Vector t(0.0); // total force on all processors
   AtomIterator iter;
   storage.begin(iter);
   for ( ; !iter.atEnd(); ++iter) {
      f += iter->force();
   }

   #ifdef UTIL_MPI
   domain.communicator().Reduce(&f[0], &t[0], 1, 
                                 MPI::DOUBLE, MPI::SUM, 0);
   domain.communicator().Reduce(&f[1], &t[1], 1, 
                                 MPI::DOUBLE, MPI::SUM, 0);
   domain.communicator().Reduce(&f[2], &t[2], 1, 
                                 MPI::DOUBLE, MPI::SUM, 0);

   if (myRank == 0) {
     std::cout << t << std::endl;
   }
   #endif

}

inline void SystemTest::testIntegrate1()
{
   printMethod(TEST_FUNC); 

   openFile("in/param2"); 
   object().readParam(file()); 

   Domain&  domain  = object().domain();
   AtomStorage& storage = object().atomStorage();
   int myRank = domain.gridRank();

   std::string filename("in/config2");
   object().readConfig(filename);
   object().exchanger().exchange();

   //object().interaction().setMethodId(0);
   object().interaction().findNeighbors();

   double temperature = 1.0;
   object().setBoltzmannVelocities(temperature);

   // Calculate energies before integration
   double kinetic   = object().kineticEnergy();
   double potential = object().pairPotentialEnergy();
   if (myRank == 0) {
      std::cout << Dbl(kinetic) << Dbl(potential) 
                << Dbl(kinetic + potential) << std::endl;
   }

   for (int i = 0; i < 3; ++i ) {

      object().integrate(500);

      // Calculate energies after integration
      kinetic   = object().kineticEnergy();
      potential = object().pairPotentialEnergy();
      if (myRank == 0) {
         std::cout << Dbl(kinetic) << Dbl(potential) 
                   << Dbl(kinetic + potential) << std::endl;
      }
  
      TEST_ASSERT(object().isValid());

   }

}

TEST_BEGIN(SystemTest)
TEST_ADD(SystemTest, testReadParam)
TEST_ADD(SystemTest, testReadConfig)
TEST_ADD(SystemTest, testExchangeAtoms)
TEST_ADD(SystemTest, testExchange)
TEST_ADD(SystemTest, testUpdate)
TEST_ADD(SystemTest, testCalculateForces)
TEST_ADD(SystemTest, testIntegrate1)
TEST_END(SystemTest)

#endif
