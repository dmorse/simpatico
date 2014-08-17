#ifndef DDMD_CONFIG_IO_TEST_H
#define DDMD_CONFIG_IO_TEST_H

#include <ddMd/configIos/DdMdConfigIo.h>
#include <ddMd/communicate/Domain.h>
#include <ddMd/communicate/Buffer.h>
#include <ddMd/communicate/GroupDistributor.h>
#include <ddMd/communicate/GroupDistributor.tpp>
#include <ddMd/communicate/GroupCollector.h>
#include <ddMd/communicate/GroupCollector.tpp>
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

class ConfigIoTest: public ParamFileTest
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

      // Domain and buffer must be initialized before the Distributor
      // and Collector objects can be associated and allocated
      atomStorage.distributor().associate(domain, boundary, atomStorage, buffer);
      atomStorage.collector().associate(domain, atomStorage, buffer);
      atomStorage.distributor().allocate(100);
      atomStorage.collector().allocate(100);
      atomStorage.readParam(file);
   
      #ifdef INTER_BOND
      bondStorage.distributor().associate(domain, atomStorage, bondStorage, buffer);
      bondStorage.collector().associate(domain, bondStorage, buffer);
      bondStorage.distributor().allocate(100);
      bondStorage.collector().allocate(100);
      bondStorage.readParam(file);
      #endif
   
      #ifdef INTER_ANGLE
      if (hasAngle) {
         angleStorage.distributor().associate(domain, atomStorage, angleStorage, buffer);
         angleStorage.collector().associate(domain, angleStorage, buffer);
         angleStorage.distributor().allocate(100);
         angleStorage.collector().allocate(100);
         angleStorage.readParam(file);
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (hasDihedral) {
         dihedralStorage.distributor().associate(domain, atomStorage, dihedralStorage, buffer);
         dihedralStorage.collector().associate(domain, dihedralStorage, buffer);
         dihedralStorage.distributor().allocate(100);
         dihedralStorage.collector().allocate(100);
         dihedralStorage.readParam(file);
      }
      #endif

      configIo.readParam(file);
      file.close();
   }

   void testReadConfig()
   {
      printMethod(TEST_FUNC);
      readParam();
      std::ifstream file;
      openInputFile("in/config", file);
      configIo.readConfig(file, MaskBonded);
   }

   void testReadWriteConfig()
   {
      printMethod(TEST_FUNC);

      readParam();
      std::ifstream file;
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
