#ifndef DDMD_CONFIG_IO_TEST_H
#define DDMD_CONFIG_IO_TEST_H

#include <ddMd/configIos/DdMdConfigIo.h>
#include <ddMd/communicate/Domain.h>
#include <ddMd/communicate/Buffer.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/storage/BondStorage.h>
#include <ddMd/storage/AngleStorage.h>
#include <ddMd/storage/DihedralStorage.h>
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

class ConfigIoTest: public ParamFileTest<DdMdConfigIo>
{
public:

   virtual void setUp()
   {}

   void testReadConfig()
   {
      printMethod(TEST_FUNC);

      DdMdConfigIo configIo;
      Boundary boundary;
      Domain   domain;
      Buffer   buffer;
      AtomStorage  atomStorage;
      BondStorage  bondStorage;
      #ifdef INTER_ANGLE
      AngleStorage  angleStorage;
      #endif
      #ifdef INTER_DIHEDRAL
      DihedralStorage  dihedralStorage;
      #endif

      // Set connections between objects
      domain.setBoundary(boundary);
      configIo.associate(domain, boundary, atomStorage, bondStorage, 
                         #ifdef INTER_ANGLE
                         angleStorage,
                         #endif
                         #ifdef INTER_DIHEDRAL
                         dihedralStorage,
                         #endif
                         buffer);

      #ifdef UTIL_MPI
      // Set communicators
      domain.setGridCommunicator(communicator());
      domain.setParamCommunicator(communicator());
      atomStorage.setParamCommunicator(communicator());
      bondStorage.setParamCommunicator(communicator());
      #ifdef INTER_ANGLE
      angleStorage.setParamCommunicator(communicator());
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralStorage.setParamCommunicator(communicator());
      #endif
      buffer.setParamCommunicator(communicator());
      configIo.setParamCommunicator(communicator());
      #else
      domain.setRank(0);
      #endif

      // Open parameter file
      std::ifstream file;
      openInputFile("in/ConfigIo", file);

      domain.readParam(file);
      atomStorage.readParam(file);
      bondStorage.readParam(file);
      #ifdef INTER_ANGLE
      angleStorage.readParam(file);
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralStorage.readParam(file);
      #endif
      buffer.readParam(file);
      configIo.readParam(file);
      file.close();

      openInputFile("in/config", file);
      configIo.readConfig(file, MaskBonded);

   }

   void testReadWriteConfig()
   {
      printMethod(TEST_FUNC);

      DdMdConfigIo configIo;
      Boundary boundary;
      Domain   domain;
      Buffer   buffer;
      AtomStorage  atomStorage;
      BondStorage  bondStorage;
      #ifdef INTER_ANGLE
      AngleStorage  angleStorage;
      #endif
      #ifdef INTER_DIHEDRAL
      DihedralStorage  dihedralStorage;
      #endif

      // Set connections between objects
      domain.setBoundary(boundary);
      configIo.associate(domain, boundary, atomStorage, bondStorage, 
                         #ifdef INTER_ANGLE
                         angleStorage,
                         #endif
                         #ifdef INTER_DIHEDRAL
                         dihedralStorage,
                         #endif
                         buffer);

      #ifdef UTIL_MPI
      // Set communicators
      domain.setGridCommunicator(communicator());
      domain.setParamCommunicator(communicator());
      atomStorage.setParamCommunicator(communicator());
      bondStorage.setParamCommunicator(communicator());
      #ifdef INTER_ANGLE
      angleStorage.setParamCommunicator(communicator());
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralStorage.setParamCommunicator(communicator());
      #endif
      buffer.setParamCommunicator(communicator());
      configIo.setParamCommunicator(communicator());
      #else
      domain.setRank(0);
      #endif

      // Open parameter file
      std::ifstream file;
      openInputFile("in/ConfigIo", file);

      domain.readParam(file);
      atomStorage.readParam(file);
      bondStorage.readParam(file);
      #ifdef INTER_ANGLE
      angleStorage.readParam(file);
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralStorage.readParam(file);
      #endif
      buffer.readParam(file);
      configIo.readParam(file);
      file.close();

      openInputFile("in/config", file);
      configIo.readConfig(file, MaskBonded);

      std::ofstream out;
      openOutputFile("out", out);
      configIo.writeConfig(out);

   }

};

TEST_BEGIN(ConfigIoTest)
TEST_ADD(ConfigIoTest, testReadConfig)
TEST_ADD(ConfigIoTest, testReadWriteConfig)
TEST_END(ConfigIoTest)

#endif /* CONFIG_IO_TEST_H */
