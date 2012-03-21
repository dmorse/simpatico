#ifndef DDMD_CONFIG_IO_TEST_H
#define DDMD_CONFIG_IO_TEST_H

#include <ddMd/configIo/ConfigIo.h>
#include <ddMd/communicate/Domain.h>
#include <ddMd/communicate/Buffer.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/storage/BondStorage.h>
#include <util/mpi/MpiLogger.h>

#ifdef UTIL_MPI
#ifndef TEST_MPI
#define TEST_MPI
#endif
#endif

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <test/ParamFileTest.h>

using namespace Util;
using namespace DdMd;

class ConfigIoTest: public ParamFileTest<ConfigIo>
{
public:

   virtual void setUp()
   {}

   void testDistribute()
   {
      printMethod(TEST_FUNC);

      Boundary boundary;
      Domain   domain;
      Buffer   buffer;
      AtomStorage  atomStorage;
      BondStorage  bondStorage;
      std::ifstream configFile;

      // Set connections between objects
      domain.setBoundary(boundary);
      object().associate(domain, boundary, atomStorage, bondStorage, buffer);

      #ifdef UTIL_MPI
      // Set communicators
      domain.setGridCommunicator(communicator());
      domain.setParamCommunicator(communicator());
      atomStorage.setParamCommunicator(communicator());
      bondStorage.setParamCommunicator(communicator());
      buffer.setParamCommunicator(communicator());
      object().setParamCommunicator(communicator());
      #else
      domain.setRank(0);
      #endif

      // Open parameter file
      openFile("in/ConfigIo");

      domain.readParam(file());
      atomStorage.readParam(file());
      bondStorage.readParam(file());
      buffer.readParam(file());
      object().readParam(file());

      // Finish reading parameter file
      closeFile();

      object().readConfig("in/config");

   }
};

TEST_BEGIN(ConfigIoTest)
TEST_ADD(ConfigIoTest, testDistribute)
TEST_END(ConfigIoTest)

#endif /* CONFIG_IO_TEST_H */
