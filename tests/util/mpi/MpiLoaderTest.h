#ifndef MPI_LOADER_TEST_H
#define MPI_LOADER_TEST_H

#include <util/mpi/MpiFileIo.h>
#include <util/mpi/MpiLoader.h>
#include <util/archives/Serializable.h>

#ifndef TEST_MPI
#define TEST_MPI
#endif

#include <mpi.h>
#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <iostream>

using namespace Util;

class MpiLoaderTest : public UnitTest
{

public:

   MpiLoaderTest()
    : UnitTest(),
      fileIo_(),
      loader_(fileIo_)
   {}

   void testSetCommunicator() 
   {
      printMethod(TEST_FUNC);
      TEST_ASSERT(!loader().hasIoCommunicator());
      loader().setIoCommunicator(communicator());
      TEST_ASSERT(loader().hasIoCommunicator());
      TEST_ASSERT(&loader().ioCommunicator() == &communicator());
      loader().clearCommunicator();
      TEST_ASSERT(!loader().hasIoCommunicator());
   }

   void testIsIoProcessor1() 
   {
      printMethod(TEST_FUNC);
      if (mpiRank() == 0) {
         TEST_ASSERT(loader().isIoProcessor());
      } else
      if (mpiRank() == 1) {
         TEST_ASSERT(loader().isIoProcessor());
      }
   }

   void testIsIoProcessor2() 
   {
      printMethod(TEST_FUNC);
      loader().setIoCommunicator(communicator());
      if (mpiRank() == 0) {
         TEST_ASSERT(loader().isIoProcessor());
      } else
      if (mpiRank() == 1) {
         TEST_ASSERT(!loader().isIoProcessor());
      }
   }

   MpiFileIo& fileIo()
   { return fileIo_; }

   MpiFileIo& loader()
   { return loader_; }

private:

   MpiFileIo   fileIo_;
   MpiLoader<Serializable::IArchive> loader_;

};

TEST_BEGIN(MpiLoaderTest)
TEST_ADD(MpiLoaderTest, testSetCommunicator)
TEST_ADD(MpiLoaderTest, testIsIoProcessor1)
TEST_ADD(MpiLoaderTest, testIsIoProcessor2)
TEST_END(MpiLoaderTest)

#endif
