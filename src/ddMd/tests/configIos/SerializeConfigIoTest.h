#ifndef DDMD_SERIALIZE_CONFIG_IO_TEST_H
#define DDMD_SERIALIZE_CONFIG_IO_TEST_H

#include <ddMd/configIos/DdMdConfigIo.h>
#include <ddMd/configIos/SerializeConfigIo.h>
#include <ddMd/communicate/Domain.h>
#include <ddMd/communicate/Buffer.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/storage/GroupStorage.tpp>
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
   bool hasAngle;
   bool hasDihedral;

public:

   virtual void setUp()
   {
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
      domain.setIoCommunicator(communicator());
      atomStorage.setIoCommunicator(communicator());
      bondStorage.setIoCommunicator(communicator());
      #ifdef INTER_ANGLE
      angleStorage.setIoCommunicator(communicator());
      #else
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralStorage.setIoCommunicator(communicator());
      #endif
      buffer.setIoCommunicator(communicator());
      configIo.setIoCommunicator(communicator());
      #else
      domain.setRank(0);
      #endif

      hasAngle = false;
      #ifdef INTER_ANGLE
      hasAngle = true;
      #endif
      hasDihedral = false;
      #ifdef INTER_DIHEDRAL
      hasDihedral = true;
      #endif
   }

   void readParam()
   {
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
      configIo.setIoCommunicator(communicator());
      #else
      domain.setRank(0);
      #endif // ifdef UTIL_MPI

      // Open parameter file
      std::ifstream file;
      if (hasDihedral) {
         openInputFile("in/ConfigIo_a_d", file);
      } else {
         if (hasAngle) {
            openInputFile("in/ConfigIo_a", file);
         } else {
            openInputFile("in/ConfigIo", file);
         }
      }

      domain.readParam(file);
      buffer.readParam(file);

      atomStorage.associate(domain, boundary, buffer);
      atomStorage.readParam(file);
   
      #ifdef INTER_BOND
      bondStorage.associate(domain, atomStorage, buffer);
      bondStorage.readParam(file);
      #endif
   
      #ifdef INTER_ANGLE
      if (hasAngle) {
         angleStorage.associate(domain, atomStorage, buffer);
         angleStorage.readParam(file);
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (hasDihedral) {
         dihedralStorage.associate(domain, atomStorage, buffer);
         dihedralStorage.readParam(file);
      }
      #endif

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
      readParam();

      std::ifstream inFile;
      openInputFile("in/config", inFile);
      configIo.readConfig(inFile, MaskBonded);
      inFile.close();

      std::ofstream outFile;
      openOutputFile("out2", outFile);
      configIo.writeConfig(outFile);
      outFile.close();

      clearStorage();

      openInputFile("out2", inFile);
      configIo.readConfig(inFile, MaskBonded);
      inFile.close();

      #if 0
      openOutputFile("out", outFile);
      configIo.writeConfig(outFile);
      outFile.close();
      #endif

   }

};

TEST_BEGIN(SerializeConfigIoTest)
TEST_ADD(SerializeConfigIoTest, testReadWriteConfig)
TEST_END(SerializeConfigIoTest)

#endif /* CONFIG_IO_TEST_H */
