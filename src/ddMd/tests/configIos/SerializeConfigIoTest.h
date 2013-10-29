#ifndef DDMD_SERIALIZE_CONFIG_IO_TEST_H
#define DDMD_SERIALIZE_CONFIG_IO_TEST_H

#include <ddMd/configIos/DdMdConfigIo.h>
#include <ddMd/configIos/SerializeConfigIo.h>
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

class SerializeConfigIoTest: public ParamFileTest
{
      SerializeConfigIo configIo;
      DdMdConfigIo      ddMdConfigIo;
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

public:

   virtual void setUp()
   {
      // Set connections between objects
      domain.setBoundary(boundary);
      ddMdConfigIo.associate(domain, boundary, atomStorage, bondStorage, 
                         #ifdef INTER_ANGLE
                         angleStorage,
                         #endif
                         #ifdef INTER_DIHEDRAL
                         dihedralStorage,
                         #endif
                         buffer);
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
      domain.setIoCommunicator(communicator());
      atomStorage.setIoCommunicator(communicator());
      bondStorage.setIoCommunicator(communicator());
      #ifdef INTER_ANGLE
      angleStorage.setIoCommunicator(communicator());
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralStorage.setIoCommunicator(communicator());
      #endif
      buffer.setIoCommunicator(communicator());
      ddMdConfigIo.setIoCommunicator(communicator());
      configIo.setIoCommunicator(communicator());
      #else
      domain.setRank(0);
      #endif

      // Open parameter file
      std::ifstream file;
      #ifdef INTER_ANGLE
         #ifdef INTER_DIHEDRAL 
         openInputFile("in/ConfigIo_a_d", file);
         #else  // ifndef INTER_DIHEDRAL
         openInputFile("in/ConfigIo_a", file);
         #endif // INTER_DIHEDRAL
      #else  // inndef INTER_ANGLE
      openInputFile("in/ConfigIo", file);
      #endif // INTER_ANGLE

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
      ddMdConfigIo.readParam(file);
      file.close();

      openInputFile("in/SerializeConfigIo", file);
      configIo.readParam(file);
      file.close();
   }

   void clearStorage() 
   {
      atomStorage.clearAtoms();
      atomStorage.clearGhosts();
      bondStorage.clearGroups();
      #ifdef INTER_ANGLE
      angleStorage.clearGroups();
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralStorage.clearGroups();
      #endif
   }

   void testReadWriteConfig()
   {
      printMethod(TEST_FUNC);

      std::ifstream inFile;
      openInputFile("in/config", inFile);
      ddMdConfigIo.readConfig(inFile, MaskBonded);
      inFile.close();

      std::ofstream outFile;
      openOutputFile("binary", outFile);
      configIo.writeConfig(outFile);
      outFile.close();

      clearStorage();

      openInputFile("binary", inFile);
      configIo.readConfig(inFile, MaskBonded);
      inFile.close();

      openOutputFile("out", outFile);
      ddMdConfigIo.writeConfig(outFile);
      outFile.close();

   }

};

TEST_BEGIN(SerializeConfigIoTest)
TEST_ADD(SerializeConfigIoTest, testReadWriteConfig)
TEST_END(SerializeConfigIoTest)

#endif /* CONFIG_IO_TEST_H */
