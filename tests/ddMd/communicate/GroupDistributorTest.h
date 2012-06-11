#ifndef DDMD_GROUP_DISTRIBUTOR_TEST_H
#define DDMD_GROUP_DISTRIBUTOR_TEST_H

#include <ddMd/communicate/AtomDistributor.h>
#include <ddMd/communicate/GroupDistributor.tpp>
#include <ddMd/communicate/Domain.h>
#include <ddMd/communicate/Buffer.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/chemistry/Atom.h>
#include <ddMd/storage/BondStorage.h>
#include <ddMd/storage/AngleStorage.h>
#include <ddMd/storage/DihedralStorage.h>
#include <util/param/Label.h>
#include <util/space/Grid.h>
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

class GroupDistributorTest: public ParamFileTest< GroupDistributor<2> >
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
      AtomDistributor  atomDistributor;
      BondStorage  bondStorage;
      GroupDistributor<2>  bondDistributor;
      AngleStorage angleStorage;
      GroupDistributor<3>  angleDistributor;
      //DihedralStorage  dihedralStorage;
      //GroupDistributor<4>  dihedralDistributor;
      std::ifstream configFile;

      // Set connections between atomDistributors
      domain.setBoundary(boundary);
      atomDistributor.associate(domain, boundary, atomStorage, buffer);
      bondDistributor.associate(domain, atomStorage, bondStorage, buffer);
      angleDistributor.associate(domain, atomStorage, angleStorage, buffer);
      //dihedralDistributor.associate(domain, atomStorage, dihedralStorage, buffer);

      #ifdef UTIL_MPI
      // Set communicators
      domain.setGridCommunicator(communicator());
      domain.setParamCommunicator(communicator());
      buffer.setParamCommunicator(communicator());
      atomStorage.setParamCommunicator(communicator());
      atomDistributor.setParamCommunicator(communicator());
      bondStorage.setParamCommunicator(communicator());
      bondDistributor.setParamCommunicator(communicator());
      angleStorage.setParamCommunicator(communicator());
      angleDistributor.setParamCommunicator(communicator());
      //dihedralStorage.setParamCommunicator(communicator());
      //dihedralDistributor.setParamCommunicator(communicator());
      #else
      domain.setRank(0);
      #endif

      // Open parameter file
      #ifdef UTIL_MPI
      openFile("in/GroupDistributor.213");
      #else
      openFile("in/GroupDistributor.111");
      #endif

      domain.readParam(file());
      atomStorage.readParam(file());
      bondStorage.readParam(file());
      angleStorage.readParam(file());
      //dihedralStorage.readParam(file());
      buffer.readParam(file());
      atomDistributor.readParam(file());
      bondDistributor.readParam(file());
      angleDistributor.readParam(file());
      //dihedralDistributor.readParam(file());

      // Finish reading parameter file
      closeFile();

      int atomCount = 0;  // Number of atoms to be distributed by master
      int bondCount = 0;  // Number to bonds be distributed by master
      int angleCount = 0; // Number to angles be distributed by master
      int myRank    = domain.gridRank();

      // If I am the master processor.
      if (myRank == 0) {
         configFile.open("in/config");
         configFile >> Label("BOUNDARY");
         configFile >> boundary;
      }
      bcast(domain.communicator(), boundary, 0);

      if (myRank == 0) {

         // Read Max number of atoms to be distributed by the master processor
         configFile >> Label("ATOMS");
         configFile >> Label("nAtom") >> atomCount;

         //std::cout << std::endl;
         //std::cout << "Num Atoms to be distributed = " 
         //          << atomCount << std::endl;

         // Initialize the sendbuffer.
         atomDistributor.setup();

         // Fill the atom atomDistributors
         Atom*   ptr;
         int     id, typeId;
         for(int i = 0; i < atomCount; ++i) {
            ptr = atomDistributor.newAtomPtr();

            configFile >> id >> typeId;
            ptr->setId(id);
            ptr->setTypeId(typeId);
            configFile >> ptr->position();
            configFile >> ptr->velocity();
            ptr->velocity() = ptr->position();

            atomDistributor.addAtom();

         }
         file().close();

         // Send any atoms not sent previously.
         atomDistributor.send();

      } else { // If I am not the master processor

         atomDistributor.receive();

      }

      int recvCount = atomStorage.nAtom();
      AtomIterator iter;
      atomStorage.begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         TEST_ASSERT(domain.isInDomain(iter->position()));
      }

      #if 0
      #ifdef UTIL_MPI
      MpiLogger logger;
      logger.begin();
      std::cout << "Processor: " << myRank
                << ", recvCount = " << recvCount << std::endl;
      logger.end();
      #endif 
      #endif 


      // Check that all atoms are accounted for after distribution.
      #ifdef UTIL_MPI
      int nRecvAll;
      communicator().Reduce(&recvCount, &nRecvAll, 1, MPI::INT, MPI::SUM, 0);
      if (myRank == 0) {
         //std::cout << "Total atom count = " << nRecvAll << std::endl;
         TEST_ASSERT(nRecvAll == atomCount);
      }
      #else
      //std::cout << "Total atom count = " << recvCount << std::endl;
      TEST_ASSERT(recvCount == atomCount);
      #endif

      // Read bonds
      if (myRank == 0) {

         // Read Max number of atoms to be distributed by the master processor
         configFile >> Label("BONDS");
         configFile >> Label("nBond") >> bondCount;

         //std::cout << std::endl;
         //std::cout << "Num Bonds to be distributed = " 
         //          << bondCount << std::endl;

         #if UTIL_MPI
         // Initialize the sendbuffer.
         bondDistributor.setup();
         #endif

         // Read bonds and add to bondDistributor
         int i;
         Group<2>* ptr;
         for(i = 0; i < bondCount; ++i) {
            ptr = bondDistributor.newPtr();
            configFile >> *ptr;
            bondDistributor.add();
         }

         // Send any bonds not sent previously.
         bondDistributor.send();

      } else { // If I am not the master processor

         bondDistributor.receive();

      }
      #if 0
      #endif

      // Read angles
      if (myRank == 0) {

         // Read Max number of atoms to be distributed by the master processor
         configFile >> Label("ANGLES");
         configFile >> Label("nAngle") >> angleCount;

         //std::cout << std::endl;
         //std::cout << "Num Bonds to be distributed = " 
         //          << angleCount << std::endl;

         #if UTIL_MPI
         // Initialize the sendbuffer.
         angleDistributor.setup();
         #endif

         // Read angles and add to angleDistributor
         int i;
         Group<3>* ptr;
         for(i = 0; i < angleCount; ++i) {
            ptr = angleDistributor.newPtr();
            configFile >> *ptr;
            angleDistributor.add();
         }

         // Send any angles not sent previously.
         angleDistributor.send();

      } else { // If I am not the master processor

         angleDistributor.receive();

      }

   }

   void testDistribute2()
   {
      printMethod(TEST_FUNC);

      Boundary boundary;
      Domain   domain;
      Buffer   buffer;
      AtomStorage atomStorage;
      BondStorage bondStorage;
      AtomDistributor atomDistributor;
      GroupDistributor<2> bondDistributor;
      std::ifstream configFile;

      // Set connections between atomDistributors
      domain.setBoundary(boundary);
      atomDistributor.associate(domain, boundary, atomStorage, buffer);
      bondDistributor.associate(domain, atomStorage, bondStorage, buffer);

      #ifdef UTIL_MPI
      // Set communicators
      domain.setGridCommunicator(communicator());
      domain.setParamCommunicator(communicator());
      atomStorage.setParamCommunicator(communicator());
      bondStorage.setParamCommunicator(communicator());
      buffer.setParamCommunicator(communicator());
      atomDistributor.setParamCommunicator(communicator());
      bondDistributor.setParamCommunicator(communicator());
      #else
      domain.setRank(0);
      #endif

      // Open parameter file
      #ifdef UTIL_MPI
      openFile("in2/BondDistributor.213");
      #else
      openFile("in2/BondDistributor.111");
      #endif

      domain.readParam(file());
      atomStorage.readParam(file());
      bondStorage.readParam(file());
      buffer.readParam(file());
      atomDistributor.readParam(file());
      bondDistributor.readParam(file());

      // Finish reading parameter file
      closeFile();

      int atomCount = 0; // Number of atoms to be distributed by master
      int bondCount = 0; // Number to bonds be distributed by master
      int myRank    = domain.gridRank();

      // If I am the master processor.
      if (myRank == 0) {
         configFile.open("in2/config");
         configFile >> Label("BOUNDARY");
         configFile >> boundary;
      }
      bcast(domain.communicator(), boundary, 0);

      if (myRank == 0) {

         // Read Max number of atoms to be distributed by the master processor
         configFile >> Label("ATOMS");
         configFile >> Label("nAtom") >> atomCount;

         //std::cout << std::endl;
         //std::cout << "Num Atoms to be distributed = " 
         //          << atomCount << std::endl;

         // Initialize the sendbuffer.
         atomDistributor.setup();

         // Fill the atom atomDistributors
         Atom*   ptr;
         int     id, typeId;
         for(int i = 0; i < atomCount; ++i) {
            ptr = atomDistributor.newAtomPtr();

            configFile >> id >> typeId;
            ptr->setId(id);
            ptr->setTypeId(typeId);
            configFile >> ptr->position();
            configFile >> ptr->velocity();
            ptr->velocity() = ptr->position();

            atomDistributor.addAtom();

         }
         file().close();

         // Send any atoms not sent previously.
         atomDistributor.send();

      } else { // If I am not the master processor

         atomDistributor.receive();

      }

      int recvCount = atomStorage.nAtom();
      AtomIterator iter;
      atomStorage.begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         TEST_ASSERT(domain.isInDomain(iter->position()));
      }

      #if 0
      #ifdef UTIL_MPI
      MpiLogger logger;
      logger.begin();
      std::cout << "Processor: " << myRank
                << ", recvCount = " << recvCount << std::endl;
      logger.end();
      #endif 
      #endif 

      // Check that all atoms are accounted for after distribution.
      #ifdef UTIL_MPI
      int nRecvAll;
      communicator().Reduce(&recvCount, &nRecvAll, 1, MPI::INT, MPI::SUM, 0);
      if (myRank == 0) {
         //std::cout << "Total atom count = " << nRecvAll << std::endl;
         TEST_ASSERT(nRecvAll == atomCount);
      }
      #else
      //std::cout << "Total atom count = " << recvCount << std::endl;
      TEST_ASSERT(recvCount == atomCount);
      #endif

      // Read bonds
      if (myRank == 0) {

         // Read Max number of atoms to be distributed by the master processor
         configFile >> Label("BONDS");
         configFile >> Label("nBond") >> bondCount;

         std::cout << std::endl;
         std::cout << "Num Bonds to be distributed = " 
                   << bondCount << std::endl;

         #if UTIL_MPI
         // Initialize the sendbuffer.
         bondDistributor.setup();
         #endif

         // Read bonds and add to bondDistributor
         int i;
         Group<2>* ptr;
         for(i = 0; i < bondCount; ++i) {
            ptr = bondDistributor.newPtr();
            configFile >> *ptr;
            bondDistributor.add();
         }

         // Send any bonds not sent previously.
         bondDistributor.send();

      } else { // If I am not the master processor

         bondDistributor.receive();

      }

      bondStorage.computeNTotal(domain.communicator());
      if (myRank == 0) {
         int nTotal = bondStorage.nTotal();
         //std::cout << "BondStorage.nTotal() =" 
         //          << nTotal << std::endl;
         TEST_ASSERT(bondCount == bondStorage.nTotal());
      }

   }
};

TEST_BEGIN(GroupDistributorTest)
TEST_ADD(GroupDistributorTest, testDistribute)
//TEST_ADD(GroupDistributorTest, testDistribute2)
TEST_END(GroupDistributorTest)

#endif /* BOND_DISTRIBUTOR_TEST_H */
