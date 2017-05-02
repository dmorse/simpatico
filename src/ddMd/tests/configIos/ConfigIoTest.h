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
   #ifdef SIMP_ANGLE
   AngleStorage  angleStorage;
   #endif
   #ifdef SIMP_DIHEDRAL
   DihedralStorage  dihedralStorage;
   #endif
   bool hasAngle;
   bool hasDihedral;

public:

   ConfigIoTest()
    : configIo(false) //hasMolecules = false
   {}

   virtual void setUp()
   {
      // Set connections between objects
      domain.setBoundary(boundary);
      configIo.associate(domain, boundary, atomStorage, bondStorage, 
                         #ifdef SIMP_ANGLE
                         angleStorage,
                         #endif
                         #ifdef SIMP_DIHEDRAL
                         dihedralStorage,
                         #endif
                         buffer);

      #ifdef UTIL_MPI
      // Set communicators
      domain.setGridCommunicator(communicator());
      domain.setIoCommunicator(communicator());
      atomStorage.setIoCommunicator(communicator());
      bondStorage.setIoCommunicator(communicator());
      #ifdef SIMP_ANGLE
      angleStorage.setIoCommunicator(communicator());
      #else
      #endif
      #ifdef SIMP_DIHEDRAL
      dihedralStorage.setIoCommunicator(communicator());
      #endif
      buffer.setIoCommunicator(communicator());
      configIo.setIoCommunicator(communicator());
      #else
      domain.setRank(0);
      #endif

      hasAngle = false;
      #ifdef SIMP_ANGLE
      hasAngle = true;
      #endif
      hasDihedral = false;
      #ifdef SIMP_DIHEDRAL
      hasDihedral = true;
      #endif
   }

   void readParam()
   {
      // Set connections between objects
      domain.setBoundary(boundary);
      configIo.associate(domain, boundary, atomStorage, bondStorage, 
                         #ifdef SIMP_ANGLE
                         angleStorage,
                         #endif
                         #ifdef SIMP_DIHEDRAL
                         dihedralStorage,
                         #endif
                         buffer);

      #ifdef UTIL_MPI
      // Set communicators
      domain.setGridCommunicator(communicator());
      domain.setIoCommunicator(communicator());
      atomStorage.setIoCommunicator(communicator());
      bondStorage.setIoCommunicator(communicator());
      #ifdef SIMP_ANGLE
      angleStorage.setIoCommunicator(communicator());
      #endif
      #ifdef SIMP_DIHEDRAL
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
   
      #ifdef SIMP_BOND
      bondStorage.associate(domain, atomStorage, buffer);
      bondStorage.readParam(file);
      #endif
   
      #ifdef SIMP_ANGLE
      if (hasAngle) {
         angleStorage.associate(domain, atomStorage, buffer);
         angleStorage.readParam(file);
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (hasDihedral) {
         dihedralStorage.associate(domain, atomStorage, buffer);
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
