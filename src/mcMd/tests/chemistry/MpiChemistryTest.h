#ifdef  UTIL_MPI
#ifndef MCMD_MPI_CHEMISTRY_TEST_H
#define MCMD_MPI_CHEMISTRY_TEST_H

#define TEST_MPI

#include <mpi.h>
#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <mcMd/chemistry/AtomType.h>
#include <mcMd/simulation/McMd_mpi.h>

#include <simp/species/SpeciesGroup.h>

#include <iostream>

using namespace Util;
using namespace Simp;
using namespace McMd;

class MpiChemistryTest : public UnitTest
{

public:

   MpiChemistryTest()
    : UnitTest()
   {  setCommunicator(MPI_COMM_WORLD); }

   void setUp()
   {}

   void tearDown()
   {}

   void testSendRecvAtomType() 
   {
      printMethod(TEST_FUNC);
      AtomType value;
      if (mpiRank() == 1) {
         TEST_ASSERT(communicator() == MPI_COMM_WORLD);
         value.setMass(5.0);
         value.setName("MyType");
         Util::send<AtomType>(communicator(), value, 0, 37);
      } else
      if (mpiRank() == 0) {
         TEST_ASSERT(communicator() == MPI_COMM_WORLD);
         Util::recv<AtomType>(communicator(), value, 1, 37);
         TEST_ASSERT(eq(value.mass(), 5.0));
         std::cout << value << std::endl;
         TEST_ASSERT(!value.name().compare("MyType"));
      }
   }

   void testBcastAtomType() 
   {
      printMethod(TEST_FUNC);
      AtomType value;
      if (mpiRank() == 1) {
         value.setMass(5.0);
         value.setName("MyType");
         Util::bcast<AtomType>(communicator(), value, 1);
      } else
      if (mpiRank() == 0) {
         Util::bcast<AtomType>(communicator(), value, 1);
         TEST_ASSERT(eq(value.mass(), 5.0));
         std::cout << value << std::endl;
         TEST_ASSERT(!value.name().compare("MyType"));
      }
   }

   void testSendRecvSpeciesGroup() 
   {
      printMethod(TEST_FUNC);
      SpeciesGroup<2> value;
      if (mpiRank() == 1) {
         value.setAtomId(0, 5);
         value.setAtomId(1, 3);
         value.setTypeId(8);
         Util::send< SpeciesGroup<2> >(communicator(), value, 0, 37);
      } else
      if (mpiRank() == 0) {
         Util::recv< SpeciesGroup<2> >(communicator(), value, 1, 37);
         TEST_ASSERT(value.atomId(0) == 5);
         TEST_ASSERT(value.atomId(1) == 3);
         TEST_ASSERT(value.typeId() == 8);
      }
   }

   void testBcastSpeciesGroup() 
   {
      printMethod(TEST_FUNC);
      SpeciesGroup<2> value;
      if (mpiRank() == 1) {
         value.setAtomId(0, 5);
         value.setAtomId(1, 3);
         value.setTypeId(8);
         Util::bcast<SpeciesGroup<2> >(communicator(), value, 1);
      } else
      if (mpiRank() == 0) {
         Util::bcast<SpeciesGroup<2> >(communicator(), value, 1);
         TEST_ASSERT(value.atomId(0) == 5);
         TEST_ASSERT(value.atomId(1) == 3);
         TEST_ASSERT(value.typeId() == 8);
      }
   }

};

TEST_BEGIN(MpiChemistryTest)
TEST_ADD(MpiChemistryTest, testSendRecvAtomType)
TEST_ADD(MpiChemistryTest, testBcastAtomType)
TEST_ADD(MpiChemistryTest, testSendRecvSpeciesGroup)
TEST_ADD(MpiChemistryTest, testBcastSpeciesGroup)
TEST_END(MpiChemistryTest)

#endif
#endif // ifdef UTIL_MPI
