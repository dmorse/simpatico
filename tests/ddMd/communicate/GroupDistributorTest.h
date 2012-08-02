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
//#include <ddMd/storage/DihedralStorage.h>
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
      #ifdef INTER_ANGLE
      AngleStorage angleStorage;
      GroupDistributor<3>  angleDistributor;
      #endif
      #ifdef INTER_DIHEDRAL
      //DihedralStorage  dihedralStorage;
      //GroupDistributor<4>  dihedralDistributor;
      #endif
      std::ifstream configFile;

      // Create associations for distributors
      domain.setBoundary(boundary);
      atomDistributor.associate(domain, boundary, atomStorage, buffer);
      bondDistributor.associate(domain, atomStorage, bondStorage, buffer);
      #ifdef INTER_ANGLE
      angleDistributor.associate(domain, atomStorage, angleStorage, buffer);
      #endif
      #ifdef INTER_DIHEDRAL
      //dihedralDistributor.associate(domain, atomStorage, dihedralStorage, 
      //                              buffer);
      #endif

      #ifdef UTIL_MPI
      // Set communicators
      domain.setGridCommunicator(communicator());
      domain.setParamCommunicator(communicator());
      buffer.setParamCommunicator(communicator());
      atomStorage.setParamCommunicator(communicator());
      atomDistributor.setParamCommunicator(communicator());
      bondStorage.setParamCommunicator(communicator());
      bondDistributor.setParamCommunicator(communicator());
      #ifdef INTER_ANGLE
      angleStorage.setParamCommunicator(communicator());
      angleDistributor.setParamCommunicator(communicator());
      #endif
      #ifdef INTER_DIHEDRAL
      //dihedralStorage.setParamCommunicator(communicator());
      //dihedralDistributor.setParamCommunicator(communicator());
      #endif
      #else // ifndef UTIL_MPI
      domain.setRank(0);
      #endif

      // Open and read parameter file
      #ifdef UTIL_MPI
      openFile("in/GroupDistributor.213");
      #else
      openFile("in/GroupDistributor.111");
      #endif

      domain.readParam(file());
      buffer.readParam(file());

      atomStorage.readParam(file());
      atomDistributor.readParam(file());
      int atomCount = 0;  // Number of atoms to be distributed by master

      bondStorage.readParam(file());
      bondDistributor.readParam(file());
      int bondCount = 0;  // Number of bonds be distributed by master

      #ifdef INTER_ANGLE
      angleStorage.readParam(file());
      angleDistributor.readParam(file());
      int angleCount = 0; // Number of angles be distributed by master
      #endif

      #ifdef INTER_DIHEDRAL
      //dihedralStorage.readParam(file());
      //dihedralDistributor.readParam(file());
      //int dihedralCount = 0; // Number of dihedrals be distributed by master
      #endif

      // Finish reading parameter file
      closeFile();

      // If I am the master processor.
      int myRank = domain.gridRank();
      if (myRank == 0) {
         configFile.open("in/config");
         configFile >> Label("BOUNDARY");
         configFile >> boundary;
      }
      bcast(domain.communicator(), boundary, 0);

      if (myRank == 0) {

         // Read number of atoms to be distributed 
         configFile >> Label("ATOMS");
         configFile >> Label("nAtom") >> atomCount;

         //std::cout << std::endl;
         //std::cout << "Num Atoms to be distributed = " 
         //          << atomCount << std::endl;

         // Initialize the sendbuffer.
         atomDistributor.setup();

         // Read atoms
         Atom*   ptr;
         int     id, typeId;
         for(int i = 0; i < atomCount; ++i) {
            ptr = atomDistributor.newAtomPtr();

            configFile >> id >> typeId;
            ptr->setId(id);
            ptr->setTypeId(typeId);

            // Read a position from file.
            if (UTIL_ORTHOGONAL) {
               configFile >> ptr->position();
            } else {
               Vector r;
               configFile >> r;
               boundary.transformCartToGen(r, ptr->position());
            }

            configFile >> ptr->velocity();

            atomDistributor.addAtom();

         }
         file().close();

         // Send any atoms not sent previously.
         atomDistributor.send();

      } else { // If I am not the master processor

         atomDistributor.receive();

      }

      // Check that all atoms are in correct processor domain.
      AtomIterator iter;
      atomStorage.begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         TEST_ASSERT(domain.isInDomain(iter->position()));
      }

      // Check that all atoms are accounted for after distribution.
      #ifdef UTIL_MPI
      atomStorage.computeNAtomTotal(communicator());
      if (myRank == 0) {
         TEST_ASSERT(atomStorage.nAtomTotal() == atomCount);
      }
      #else
      TEST_ASSERT(atomStorage.nAtom() == atomCount);
      #endif

      // Read bonds
      if (myRank == 0) {

         // Read number of bonds to be distributed 
         configFile >> Label("BONDS");
         configFile >> Label("nBond") >> bondCount;

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
      // Note: Validation is done inside send and receive methods

      #ifdef INTER_ANGLE
      // Read angles
      if (myRank == 0) {

         // Read number of angles to be distributed 
         configFile >> Label("ANGLES");
         configFile >> Label("nAngle") >> angleCount;

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
      // Note: Validation is done inside send and receive methods
      #endif

   }

   void testDistribute2()
   {
      printMethod(TEST_FUNC);

      Boundary boundary;
      Domain   domain;
      Buffer   buffer;
      AtomStorage  atomStorage;
      AtomDistributor  atomDistributor;
      BondStorage  bondStorage;
      GroupDistributor<2>  bondDistributor;
      #ifdef INTER_ANGLE
      AngleStorage angleStorage;
      GroupDistributor<3>  angleDistributor;
      #endif
      #ifdef INTER_DIHEDRAL
      //DihedralStorage  dihedralStorage;
      //GroupDistributor<4>  dihedralDistributor;
      #endif
      std::ifstream configFile;

      // Create associations for distributors
      domain.setBoundary(boundary);
      atomDistributor.associate(domain, boundary, atomStorage, buffer);
      bondDistributor.associate(domain, atomStorage, bondStorage, buffer);
      #ifdef INTER_ANGLE
      angleDistributor.associate(domain, atomStorage, angleStorage, buffer);
      #endif
      #ifdef INTER_DIHEDRAL
      //dihedralDistributor.associate(domain, atomStorage, dihedralStorage, 
      //                              buffer);
      #endif

      #ifdef UTIL_MPI
      // Set communicators
      domain.setGridCommunicator(communicator());
      domain.setParamCommunicator(communicator());
      buffer.setParamCommunicator(communicator());
      atomStorage.setParamCommunicator(communicator());
      atomDistributor.setParamCommunicator(communicator());
      bondStorage.setParamCommunicator(communicator());
      bondDistributor.setParamCommunicator(communicator());
      #ifdef INTER_ANGLE
      angleStorage.setParamCommunicator(communicator());
      angleDistributor.setParamCommunicator(communicator());
      #endif
      #ifdef INTER_DIHEDRAL
      //dihedralStorage.setParamCommunicator(communicator());
      //dihedralDistributor.setParamCommunicator(communicator());
      #endif
      #else // ifndef UTIL_MPI
      domain.setRank(0);
      #endif

      // Open and read parameter file
      #ifdef UTIL_MPI
      openFile("in2/GroupDistributor.213");
      #else
      openFile("in2/GroupDistributor.111");
      #endif

      domain.readParam(file());
      buffer.readParam(file());

      atomStorage.readParam(file());
      atomDistributor.readParam(file());
      int atomCount = 0;  // Number of atoms to be distributed by master

      bondStorage.readParam(file());
      bondDistributor.readParam(file());
      int bondCount = 0;  // Number of bonds be distributed by master

      #ifdef INTER_ANGLE
      angleStorage.readParam(file());
      angleDistributor.readParam(file());
      int angleCount = 0; // Number of angles be distributed by master
      #endif

      #ifdef INTER_DIHEDRAL
      //dihedralStorage.readParam(file());
      //dihedralDistributor.readParam(file());
      //int dihedralCount = 0; // Number of dihedrals be distributed by master
      #endif

      // Finish reading parameter file
      closeFile();

      // If I am the master processor.
      int myRank = domain.gridRank();
      if (myRank == 0) {
         configFile.open("in2/config");
         configFile >> Label("BOUNDARY");
         configFile >> boundary;
      }
      bcast(domain.communicator(), boundary, 0);

      if (myRank == 0) {

         // Read number of atoms to be distributed 
         configFile >> Label("ATOMS");
         configFile >> Label("nAtom") >> atomCount;

         //std::cout << std::endl;
         //std::cout << "Num Atoms to be distributed = " 
         //          << atomCount << std::endl;

         // Initialize the sendbuffer.
         atomDistributor.setup();

         // Read atoms
         Atom*   ptr;
         int     id, typeId;
         for(int i = 0; i < atomCount; ++i) {
            ptr = atomDistributor.newAtomPtr();

            configFile >> id >> typeId;
            ptr->setId(id);
            ptr->setTypeId(typeId);

            // Read a position from file.
            if (UTIL_ORTHOGONAL) {
               configFile >> ptr->position();
            } else {
               Vector r;
               configFile >> r;
               boundary.transformCartToGen(r, ptr->position());
            }

            configFile >> ptr->velocity();
            //ptr->velocity() = ptr->position();

            atomDistributor.addAtom();

         }
         file().close();

         // Send any atoms not sent previously.
         atomDistributor.send();

      } else { // If I am not the master processor

         atomDistributor.receive();

      }

      // Check that all atoms are in correct processor domain.
      AtomIterator iter;
      atomStorage.begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         TEST_ASSERT(domain.isInDomain(iter->position()));
      }

      // Check that all atoms are accounted for after distribution.
      #ifdef UTIL_MPI
      atomStorage.computeNAtomTotal(communicator());
      if (myRank == 0) {
         TEST_ASSERT(atomStorage.nAtomTotal() == atomCount);
      }
      #else
      TEST_ASSERT(atomStorage.nAtom() == atomCount);
      #endif

      // Read bonds
      if (myRank == 0) {

         // Read number of bonds to be distributed 
         configFile >> Label("BONDS");
         configFile >> Label("nBond") >> bondCount;

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
      // Note: Validation is done inside send and receive methods

      #ifdef INTER_ANGLE
      // Read angles
      if (myRank == 0) {

         // Read number of angles to be distributed 
         configFile >> Label("ANGLES");
         configFile >> Label("nAngle") >> angleCount;

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
      // Note: Validation is done inside send and receive methods
      #endif

   }
};

TEST_BEGIN(GroupDistributorTest)
TEST_ADD(GroupDistributorTest, testDistribute)
TEST_ADD(GroupDistributorTest, testDistribute2)
TEST_END(GroupDistributorTest)

#endif // DDMD_GROUP_DISTRIBUTOR_TEST_H 
