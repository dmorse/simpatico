#ifndef DDMD_BUFFER_TEST_H
#define DDMD_BUFFER_TEST_H

#include <ddMd/communicate/Buffer.h>
#include <ddMd/chemistry/Atom.h>
#include <ddMd/chemistry/AtomArray.h>
#include <ddMd/chemistry/Bond.h>
#include <util/math/feq.h>
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

class BufferTest: public ParamFileTest
{

   Buffer buffer_;

public:

   virtual void setUp()
   {
      // Allocate a Buffer with space
      // for 4 local atoms or 4 ghosts.
      // buffer_.allocate(4, 4);

      #ifdef UTIL_MPI
      buffer_.setIoCommunicator(communicator());
      #endif

      // Read parameter file
      openFile("in/Buffer");
      buffer_.readParam(file());
      closeFile();

   }

   void testCapacities()
   {
      printMethod(TEST_FUNC);
      TEST_ASSERT(buffer_.atomCapacity() == 4);
      TEST_ASSERT(buffer_.ghostCapacity() > 4);
   }

   #ifdef UTIL_MPI
   void testPackAtom()
   {
      printMethod(TEST_FUNC);

      AtomArray atoms;
      Vector    pos, vel;
      atoms.allocate(3);

      // Atom 0
      atoms[0].setId(0);
      atoms[0].setTypeId(0);
      pos[0] = 1.0;
      pos[1] = 2.0;
      pos[2] = 3.0;
      atoms[0].position() = pos;
      vel[0] = 3.0;
      vel[1] = 5.0;
      vel[2] = 7.0;
      atoms[0].velocity() = vel;
      atoms[0].plan().setFlags(7);

      // Atom 1
      atoms[1].setId(1);
      atoms[1].setTypeId(0);
      pos[0] = 1.3;
      pos[1] = 2.5;
      pos[2] = 3.7;
      atoms[1].position() = pos;
      vel[0] = 1.3;
      vel[1] = 3.5;
      vel[2] = 67.7;
      atoms[1].velocity() = vel;
      atoms[1].plan().setFlags(7);

      TEST_ASSERT(buffer_.isAllocated());
      buffer_.clearSendBuffer();

      buffer_.beginSendBlock(Buffer::ATOM);
      //buffer_.packAtom(atoms[0]);
      //buffer_.packAtom(atoms[1]);
      atoms[0].packAtom(buffer_);
      atoms[1].packAtom(buffer_);
      buffer_.endSendBlock(Buffer::ATOM);
   }

   void testPackGhost()
   {
      printMethod(TEST_FUNC);

      AtomArray atoms;
      Vector pos;
      atoms.allocate(3);

      atoms[0].setId(2);
      atoms[0].setTypeId(1);
      pos[0] = 1.9;
      pos[1] = 4.2;
      pos[2] = 5.3;
      atoms[0].position() = pos;
      atoms[0].plan().setFlags(7);

      atoms[1].setId(3);
      atoms[1].setTypeId(1);
      pos[0] = 1.1;
      pos[1] = 4.4;
      pos[2] = 5.5;
      atoms[1].position() = pos;
      atoms[1].plan().setFlags(7);

      TEST_ASSERT(buffer_.isAllocated());
      buffer_.clearSendBuffer();

      buffer_.beginSendBlock(Buffer::GHOST);
      TEST_ASSERT(buffer_.sendSize() == 0);
      //buffer_.packGhost(atoms[0]);
      //buffer_.packGhost(atoms[1]);
      atoms[0].packGhost(buffer_);
      atoms[1].packGhost(buffer_);
      buffer_.endSendBlock();

   }

   //Test Method for MPI Sendrecv method for local atoms
   void testAtomSendRecv()
   {
      printMethod(TEST_FUNC);
      //std::cout << std::endl;

      AtomArray atoms;
      Vector    pos, vel;
      int       myrank, commsize;

      atoms.allocate(4);
      myrank   = MPI::COMM_WORLD.Get_rank();
      commsize = MPI::COMM_WORLD.Get_size();

      // Fill one local atom object. 
      // Add processor's rank to the position and velocity vectors.
      atoms[0].setId(0);
      atoms[0].setTypeId(0);
      pos[0] = myrank  + 1.0;
      pos[1] = myrank  + 2.0;
      pos[2] = myrank  + 3.0;
      atoms[0].position() = pos;
      vel[0] = myrank  + 4.0;
      vel[1] = myrank  + 5.0;
      vel[2] = myrank  + 6.0;
      atoms[0].velocity() = vel;
      atoms[0].plan().setFlags(7);

      // Fill another local atom object
      atoms[1].setId(1);
      atoms[1].setTypeId(0);
      pos[0] = myrank  + 10.3;
      pos[1] = myrank  + 20.5;
      pos[2] = myrank  + 30.7;
      atoms[1].position() = pos;
      vel[0] = myrank  + 40.3;
      vel[1] = myrank  + 50.3;
      vel[2] = myrank  + 60.3;
      atoms[1].velocity() = vel;
      atoms[1].plan().setFlags(7);

      //Initialize the sendbuffer, set atomtype to ATOM
      buffer_.clearSendBuffer();

      // Pack 2 local atoms into the send buffer
      buffer_.beginSendBlock(Buffer::ATOM);
      //buffer_.packAtom(atoms[0]);
      //buffer_.packAtom(atoms[1]);
      atoms[0].packAtom(buffer_);
      atoms[1].packAtom(buffer_);
      buffer_.endSendBlock();

      // Send the sendbuffer to processor to the right. 
      // Receive from the processor on the left.
      int  source = (myrank + commsize - 1) % commsize;
      int  dest   = (myrank + commsize + 1) % commsize;
      buffer_.sendRecv(MPI::COMM_WORLD, source, dest);

      // Unpack atoms
      buffer_.beginRecvBlock();
      TEST_ASSERT(buffer_.recvSize() == 2);
      atoms[2].unpackAtom(buffer_);
      atoms[3].unpackAtom(buffer_);
      TEST_ASSERT(buffer_.recvSize() == 0);
      buffer_.endRecvBlock();

      TEST_ASSERT(atoms[2].id() == 0);
      TEST_ASSERT(atoms[2].typeId() == 0);
      TEST_ASSERT(feq(atoms[2].position()[0], 1.0 + double(source)));
      TEST_ASSERT(feq(atoms[2].position()[1], 2.0 + double(source)));
      TEST_ASSERT(feq(atoms[2].position()[2], 3.0 + double(source)));
      TEST_ASSERT(feq(atoms[2].velocity()[0], 4.0 + double(source)));
      TEST_ASSERT(feq(atoms[2].velocity()[1], 5.0 + double(source)));
      TEST_ASSERT(feq(atoms[2].velocity()[2], 6.0 + double(source)));
      TEST_ASSERT(atoms[2].plan().flags() == 7);

      TEST_ASSERT(atoms[3].id() == 1);
      TEST_ASSERT(atoms[3].typeId() == 0);
      TEST_ASSERT(feq(atoms[3].position()[0], 10.3 + double(source)));
      TEST_ASSERT(feq(atoms[3].position()[1], 20.5 + double(source)));
      TEST_ASSERT(feq(atoms[3].position()[2], 30.7 + double(source)));
      TEST_ASSERT(feq(atoms[3].velocity()[0], 40.3 + double(source)));
      TEST_ASSERT(feq(atoms[3].velocity()[1], 50.3 + double(source)));
      TEST_ASSERT(feq(atoms[3].velocity()[2], 60.3 + double(source)));
      TEST_ASSERT(atoms[3].plan().flags() == 7);

      #if 0
      MpiLogger logger;
      logger.begin();
      int id, typeId;
      id     = atoms[2].id();
      typeId = atoms[2].typeId();
      pos = atoms[2].position();
      vel = atoms[2].velocity();

      std::cout << std::endl;
      std::cout << "proc: " << myrank 
                << ", source " << source << ", dest   " <<  dest
                << ", comm size " << commsize << std::endl
                << "proc: " << myrank 
                << ", id " << id << ", typeId " << typeId << std::endl
                << "proc: " << myrank << ", position = "
                << pos[0] << " " << pos[1] << " " << pos[2] << std::endl
                << "proc " << myrank 
                << "velocity = "
                <<  vel[0] << ", " << vel[1] << ", " << vel[2] << std::endl;


      id     = atoms[3].id();
      typeId = atoms[3].typeId();
      pos = atoms[3].position();
      vel = atoms[3].velocity();

      std::cout << "proc: " << myrank 
                << ", source " << source << ", dest   " <<  dest
                << ", comm size " << commsize << std::endl
                << "proc: " << myrank 
                << ", id " << id << ", typeId " << typeId << std::endl
                << "proc: " << myrank << ", position = "
                << pos[0] << " " << pos[1] << " " << pos[2] << std::endl
                << "proc " << myrank 
                << "velocity = "
                <<  vel[0] << ", " << vel[1] << ", " << vel[2] << std::endl;


      logger.end();
      #endif

   }

   //Test Method for MPI Sendrecv method for ghost atoms
   void testGhostSendRecv()
   {
      printMethod(TEST_FUNC);
      //std::cout << std::endl;

      AtomArray atoms;
      Vector    pos;
      int myrank, commsize;
      atoms.allocate(4);

      myrank   = MPI::COMM_WORLD.Get_rank();
      commsize = MPI::COMM_WORLD.Get_size();

      // Fill one local atom object. Add processor's rank to the position
      // and velocity vectors.
      atoms[0].setId(0);
      atoms[0].setTypeId(0);
      pos[0] = myrank  + 10.0;
      pos[1] = myrank  + 20.0;
      pos[2] = myrank  + 30.0;
      atoms[0].position() = pos;
      atoms[0].plan().setFlags(3);

      //Fill another local atom object
      atoms[1].setId(1);
      atoms[1].setTypeId(0);
      pos[0] = myrank  + 100.0;
      pos[1] = myrank  + 200.0;
      pos[2] = myrank  + 300.0;
      atoms[1].position() = pos;
      atoms[1].plan().setFlags(2);

      //Initialize the sendbuffer, set atomtype to GHOST
      buffer_.clearSendBuffer();

      //Pack 2 local atoms into the send buffer
      buffer_.beginSendBlock(Buffer::GHOST);
      //buffer_.packGhost(atoms[0]);
      //buffer_.packGhost(atoms[1]);
      atoms[0].packGhost(buffer_);
      atoms[1].packGhost(buffer_);
      buffer_.endSendBlock();

      // Send the sendbuffer to processor to the right. Receive from the
      int  source = (myrank + commsize - 1) % commsize;
      int  dest   = (myrank + commsize + 1) % commsize;
      buffer_.sendRecv(MPI::COMM_WORLD, source, dest);

      // Unpack ghost atoms
      buffer_.beginRecvBlock();
      TEST_ASSERT(buffer_.recvSize() == 2);
      atoms[2].unpackGhost(buffer_);
      TEST_ASSERT(buffer_.recvSize() == 1);
      atoms[3].unpackGhost(buffer_);
      TEST_ASSERT(buffer_.recvSize() == 0);
      buffer_.endRecvBlock();

      TEST_ASSERT(atoms[2].id() == 0);
      TEST_ASSERT(atoms[2].typeId() == 0);
      TEST_ASSERT(feq(atoms[2].position()[0], 10.0 + double(source)));
      TEST_ASSERT(feq(atoms[2].position()[1], 20.0 + double(source)));
      TEST_ASSERT(feq(atoms[2].position()[2], 30.0 + double(source)));
      TEST_ASSERT(atoms[2].plan().flags() == 3);

      TEST_ASSERT(atoms[3].id() == 1);
      TEST_ASSERT(atoms[3].typeId() == 0);
      TEST_ASSERT(feq(atoms[3].position()[0], 100.0 + double(source)));
      TEST_ASSERT(feq(atoms[3].position()[1], 200.0 + double(source)));
      TEST_ASSERT(feq(atoms[3].position()[2], 300.0 + double(source)));
      TEST_ASSERT(atoms[3].plan().flags() == 2);

   }

   //Test Method for MPI Send and MPI Recv methods for local atoms
   void testAtomSend_Recv()
   {
      printMethod(TEST_FUNC);
      //std::cout << std::endl;

      AtomArray atoms;
      Vector    pos, vel;
      int myrank, commsize;

      atoms.allocate(4);
      myrank   = MPI::COMM_WORLD.Get_rank();
      commsize = MPI::COMM_WORLD.Get_size();

      // Fill one local atom object. Add processor's rank to the
      // position and velocity vectors.
      atoms[0].setId(0);
      atoms[0].setTypeId(0);
      pos[0] = myrank  + 1.0;
      pos[1] = myrank  + 2.0;
      pos[2] = myrank  + 3.0;
      atoms[0].position() = pos;
      vel[0] = myrank  + 4.0;
      vel[1] = myrank  + 5.0;
      vel[2] = myrank  + 6.0;
      atoms[0].velocity() = vel;

      //Fill another local atom object
      atoms[1].setId(1);
      atoms[1].setTypeId(0);
      pos[0] = myrank  + 10.3;
      pos[1] = myrank  + 20.3;
      pos[2] = myrank  + 30.3;
      atoms[1].position() = pos;
      vel[0] = myrank  + 40.3;
      vel[1] = myrank  + 50.3;
      vel[2] = myrank  + 60.3;
      atoms[1].velocity() = vel;

      // Initialize the sendbuffer, set atomtype to ATOM
      buffer_.clearSendBuffer();

      // Pack 2 local atoms into the send buffer
      // Set isComplete parameter false.
      buffer_.beginSendBlock(Buffer::ATOM);
      atoms[0].packAtom(buffer_);
      atoms[1].packAtom(buffer_);
      buffer_.endSendBlock(false);

      // Send the sendbuffer to processor to the right.
      int  dest   = (myrank + commsize + 1) % commsize;
      buffer_.send(MPI::COMM_WORLD, dest);

      // Receive from the processor on the left.
      int  source = (myrank + commsize - 1) % commsize;
      buffer_.recv(MPI::COMM_WORLD, source);

      // Unpack 2 atoms from receive buffer
      bool isComplete = buffer_.beginRecvBlock();
      TEST_ASSERT(!isComplete);
      atoms[2].unpackAtom(buffer_);
      atoms[3].unpackAtom(buffer_);
      TEST_ASSERT(buffer_.recvSize() == 0);
      buffer_.endRecvBlock();

      TEST_ASSERT(atoms[2].id() == 0);
      TEST_ASSERT(atoms[2].typeId() == 0);
      TEST_ASSERT(feq(atoms[2].position()[0], 1.0 + double(source)));
      TEST_ASSERT(feq(atoms[2].position()[1], 2.0 + double(source)));
      TEST_ASSERT(feq(atoms[2].position()[2], 3.0 + double(source)));
      TEST_ASSERT(feq(atoms[2].velocity()[0], 4.0 + double(source)));
      TEST_ASSERT(feq(atoms[2].velocity()[1], 5.0 + double(source)));
      TEST_ASSERT(feq(atoms[2].velocity()[2], 6.0 + double(source)));

      TEST_ASSERT(atoms[3].id() == 1);
      TEST_ASSERT(atoms[3].typeId() == 0);
      TEST_ASSERT(feq(atoms[3].position()[0], 10.3 + double(source)));
      TEST_ASSERT(feq(atoms[3].position()[1], 20.3 + double(source)));
      TEST_ASSERT(feq(atoms[3].position()[2], 30.3 + double(source)));
      TEST_ASSERT(feq(atoms[3].velocity()[0], 40.3 + double(source)));
      TEST_ASSERT(feq(atoms[3].velocity()[1], 50.3 + double(source)));
      TEST_ASSERT(feq(atoms[3].velocity()[2], 60.3 + double(source)));

      #if 0
      MpiLogger logger;
      logger.begin();
      int id, typeId;
      id     = atoms[2].id();
      typeId = atoms[2].typeId();
      pos = atoms[2].position();
      vel = atoms[2].velocity();

      std::cout << std::endl;
      std::cout << "proc: " << myrank 
                << ", source " << source << ", dest   " <<  dest
                << ", comm size " << commsize << std::endl
                << "proc: " << myrank << ": " 
                << ", id " << id << ", typeId " << typeId << std::endl
                << "proc: " << myrank << ": " 
                << ", position = "
                << pos[0] << " " << pos[1] << " " << pos[2] << std::endl
                << "proc " << myrank << ": " 
                << "velocity = "
                <<  vel[0] << ", " << vel[1] << ", " << vel[2] << std::endl;

      id     = atoms[3].id();
      typeId = atoms[3].typeId();
      pos = atoms[3].position();
      vel = atoms[3].velocity();

      std::cout << "proc: " << myrank 
                << ", source " << source << ", dest   " <<  dest
                << ", comm size " << commsize << std::endl
                << "proc: " << myrank << ": " 
                << ", id " << id << ", typeId " << typeId << std::endl
                << "proc " << myrank << ": " 
                << ", position = "
                << pos[0] << " " << pos[1] << " " << pos[2] << std::endl
                << "proc " << myrank << ": " 
                << "velocity = "
                <<  vel[0] << ", " << vel[1] << ", " << vel[2] << std::endl;

      logger.end();
      #endif

   }

   //Test Method for MPI Send and MPI Recv methods for ghost atoms
   void testGhostSend_Recv()
   {
      printMethod(TEST_FUNC);
      //std::cout << std::endl;

      AtomArray atoms;
      Vector    pos;
      int myrank, commsize;
      atoms.allocate(4);

      myrank   = MPI::COMM_WORLD.Get_rank();
      commsize = MPI::COMM_WORLD.Get_size();

      // Fill one local atom object. Add processor's rank to the position
      // and velocity vectors.
      atoms[0].setId(0);
      atoms[0].setTypeId(0);
      pos[0] = myrank  + 10.0;
      pos[1] = myrank  + 20.0;
      pos[2] = myrank  + 30.0;
      atoms[0].position() = pos;

      //Fill another local atom object
      atoms[1].setId(1);
      atoms[1].setTypeId(0);
      pos[0] = myrank  + 100.3;
      pos[1] = myrank  + 200.3;
      pos[2] = myrank  + 300.3;
      atoms[1].position() = pos;

      //Initialize the sendbuffer, set atomtype to GHOST
      buffer_.clearSendBuffer();
      buffer_.beginSendBlock(Buffer::GHOST);

      //Pack 2 local atoms into the send buffer
      //buffer_.packGhost(atoms[0]);
      //buffer_.packGhost(atoms[1]);
      atoms[0].packGhost(buffer_);
      atoms[1].packGhost(buffer_);
      buffer_.endSendBlock(true);

      // Send the sendbuffer to processor to the right.
      int  dest   = (myrank + commsize + 1) % commsize;
      buffer_.send(MPI::COMM_WORLD, dest);

      // Receive from the processor on the left.
      int  source = (myrank + commsize - 1) % commsize;
      buffer_.recv(MPI::COMM_WORLD, source);

      // Unpack ghost atoms
      bool isComplete = buffer_.beginRecvBlock();
      TEST_ASSERT(isComplete);
      TEST_ASSERT(buffer_.recvSize() == 2);
      atoms[2].unpackGhost(buffer_);
      TEST_ASSERT(buffer_.recvSize() == 1);
      atoms[3].unpackGhost(buffer_);
      TEST_ASSERT(buffer_.recvSize() == 0);
      buffer_.endRecvBlock();

      TEST_ASSERT(atoms[2].id() == 0);
      TEST_ASSERT(atoms[2].typeId() == 0);
      TEST_ASSERT(feq(atoms[2].position()[0], 10.0 + double(source)));
      TEST_ASSERT(feq(atoms[2].position()[1], 20.0 + double(source)));
      TEST_ASSERT(feq(atoms[2].position()[2], 30.0 + double(source)));

      TEST_ASSERT(atoms[3].id() == 1);
      TEST_ASSERT(atoms[3].typeId() == 0);
      TEST_ASSERT(feq(atoms[3].position()[0], 100.3 + double(source)));
      TEST_ASSERT(feq(atoms[3].position()[1], 200.3 + double(source)));
      TEST_ASSERT(feq(atoms[3].position()[2], 300.3 + double(source)));

   }

   // Test Method for MPI Send and MPI Recv methods for local atoms
   void testAtomGhostSend_Recv()
   {
      printMethod(TEST_FUNC);
      //std::cout << std::endl;

      AtomArray atoms;
      Vector    pos, vel;
      int myrank, commsize;

      atoms.allocate(8);
      myrank   = MPI::COMM_WORLD.Get_rank();
      commsize = MPI::COMM_WORLD.Get_size();

      // Atom[0]. 
      // Add processor's rank to the position and velocity vectors.
      atoms[0].setId(0);
      atoms[0].setTypeId(0);
      pos[0] = myrank  + 1.0;
      pos[1] = myrank  + 2.0;
      pos[2] = myrank  + 3.0;
      atoms[0].position() = pos;
      vel[0] = myrank  + 4.0;
      vel[1] = myrank  + 5.0;
      vel[2] = myrank  + 6.0;
      atoms[0].velocity() = vel;

      // Atom[1]. 
      atoms[1].setId(1);
      atoms[1].setTypeId(0);
      pos[0] = myrank  + 10.3;
      pos[1] = myrank  + 20.3;
      pos[2] = myrank  + 30.3;
      atoms[1].position() = pos;
      vel[0] = myrank  + 40.3;
      vel[1] = myrank  + 50.3;
      vel[2] = myrank  + 60.3;
      atoms[1].velocity() = vel;

      // Atom[2]. 
      // Add processor's rank to the position and velocity vectors.
      atoms[2].setId(0);
      atoms[2].setTypeId(0);
      pos[0] = myrank  + 15.0;
      pos[1] = myrank  + 25.0;
      pos[2] = myrank  + 35.0;
      atoms[2].position() = pos;
      vel[0] = myrank  + 45.0;
      vel[1] = myrank  + 55.0;
      vel[2] = myrank  + 65.0;
      atoms[2].velocity() = vel;

      // Atom[3]. 
      atoms[3].setId(1);
      atoms[3].setTypeId(0);
      pos[0] = myrank  + 12.3;
      pos[1] = myrank  + 22.3;
      pos[2] = myrank  + 32.3;
      atoms[3].position() = pos;
      vel[0] = myrank  + 42.3;
      vel[1] = myrank  + 52.3;
      vel[2] = myrank  + 62.3;

      buffer_.clearSendBuffer();

      // Pack block of 2 local atoms into the send buffer
      // Set isComplete parameter false.
      TEST_ASSERT(buffer_.sendSize() == 0);
      buffer_.beginSendBlock(Buffer::ATOM);
      atoms[0].packAtom(buffer_);
      TEST_ASSERT(buffer_.sendSize() == 1);
      atoms[1].packAtom(buffer_);
      TEST_ASSERT(buffer_.sendSize() == 2);
      buffer_.endSendBlock(false);

      // Pack block 2 ghost atoms into the send buffer
      buffer_.beginSendBlock(Buffer::GHOST);
      atoms[2].packGhost(buffer_);
      atoms[3].packGhost(buffer_);
      buffer_.endSendBlock(true); // set isComplete = true

      // Send sendbuffer to next processor in ring.
      int  dest   = (myrank + commsize + 1) % commsize;
      buffer_.send(MPI::COMM_WORLD, dest);

      // Receive from previous processor in ring.
      int  source = (myrank + commsize - 1) % commsize;
      buffer_.recv(MPI::COMM_WORLD, source);

      // Unpack block of 2 local atoms from receive buffer
      bool isComplete;
      TEST_ASSERT(buffer_.recvSize() == 0);
      isComplete = buffer_.beginRecvBlock();
      TEST_ASSERT(!isComplete);
      TEST_ASSERT(buffer_.recvSize() == 2);
      atoms[4].unpackAtom(buffer_);
      TEST_ASSERT(buffer_.recvSize() == 1);
      atoms[5].unpackAtom(buffer_);
      TEST_ASSERT(buffer_.recvSize() == 0);
      buffer_.endRecvBlock();

      TEST_ASSERT(atoms[4].id() == 0);
      TEST_ASSERT(atoms[4].typeId() == 0);
      TEST_ASSERT(feq(atoms[4].position()[0], 1.0 + double(source)));
      TEST_ASSERT(feq(atoms[4].position()[1], 2.0 + double(source)));
      TEST_ASSERT(feq(atoms[4].position()[2], 3.0 + double(source)));
      TEST_ASSERT(feq(atoms[4].velocity()[0], 4.0 + double(source)));
      TEST_ASSERT(feq(atoms[4].velocity()[1], 5.0 + double(source)));
      TEST_ASSERT(feq(atoms[4].velocity()[2], 6.0 + double(source)));

      TEST_ASSERT(atoms[5].id() == 1);
      TEST_ASSERT(atoms[5].typeId() == 0);
      TEST_ASSERT(feq(atoms[5].position()[0], 10.3 + double(source)));
      TEST_ASSERT(feq(atoms[5].position()[1], 20.3 + double(source)));
      TEST_ASSERT(feq(atoms[5].position()[2], 30.3 + double(source)));
      TEST_ASSERT(feq(atoms[5].velocity()[0], 40.3 + double(source)));
      TEST_ASSERT(feq(atoms[5].velocity()[1], 50.3 + double(source)));
      TEST_ASSERT(feq(atoms[5].velocity()[2], 60.3 + double(source)));

      // Unpack block of 2 ghost atoms from receive buffer
      TEST_ASSERT(buffer_.recvSize() == 0);
      isComplete = buffer_.beginRecvBlock();
      TEST_ASSERT(isComplete);
      atoms[6].unpackGhost(buffer_);
      atoms[7].unpackGhost(buffer_);
      TEST_ASSERT(buffer_.recvSize() == 0);
      buffer_.endRecvBlock();

      TEST_ASSERT(atoms[6].id() == 0);
      TEST_ASSERT(atoms[6].typeId() == 0);
      TEST_ASSERT(feq(atoms[6].position()[0], 15.0 + double(source)));
      TEST_ASSERT(feq(atoms[6].position()[1], 25.0 + double(source)));
      TEST_ASSERT(feq(atoms[6].position()[2], 35.0 + double(source)));

      TEST_ASSERT(atoms[7].id() == 1);
      TEST_ASSERT(atoms[7].typeId() == 0);
      TEST_ASSERT(feq(atoms[7].position()[0], 12.3 + double(source)));
      TEST_ASSERT(feq(atoms[7].position()[1], 22.3 + double(source)));
      TEST_ASSERT(feq(atoms[7].position()[2], 32.3 + double(source)));

   }

   // Test Method for MPI Sendrecv method for bonds
   void testBondSendRecv()
   {
      printMethod(TEST_FUNC);
      //std::cout << std::endl;

      DArray<Bond> bonds;
      int          myrank, commsize;

      bonds.allocate(4);
      myrank   = MPI::COMM_WORLD.Get_rank();
      commsize = MPI::COMM_WORLD.Get_size();

      // Fill one local atom object. 
      bonds[0].setId(0);
      bonds[0].setTypeId(0);
      bonds[0].setAtomId(0, myrank + 34);
      bonds[0].setAtomId(1, myrank + 35);

      // Fill another local atom object
      bonds[1].setId(1);
      bonds[1].setTypeId(0);
      bonds[1].setAtomId(0, myrank + 38);
      bonds[1].setAtomId(1, myrank + 39);

      //Initialize the sendbuffer, set atomtype to ATOM
      buffer_.clearSendBuffer();
      buffer_.beginSendBlock(Buffer::GROUP2);

      // Pack 2 bonds into the send buffer
      //buffer_.packGroup(bonds[0]);
      //buffer_.packGroup(bonds[1]);
      bonds[0].pack(buffer_);
      bonds[1].pack(buffer_);
      buffer_.endSendBlock();

      // Send the sendbuffer to processor to the right. 
      // Receive from the processor on the left.
      int  source = (myrank + commsize - 1) % commsize;
      int  dest   = (myrank + commsize + 1) % commsize;
      buffer_.sendRecv(MPI::COMM_WORLD, source, dest);

      // Unpack bonds
      buffer_.beginRecvBlock();
      TEST_ASSERT(buffer_.recvSize() == 2);
      bonds[2].unpack(buffer_);
      bonds[3].unpack(buffer_);
      TEST_ASSERT(buffer_.recvSize() == 0);
      buffer_.endRecvBlock();

      TEST_ASSERT(bonds[2].id() == 0);
      TEST_ASSERT(bonds[2].typeId() == 0);
      TEST_ASSERT(bonds[2].atomId(0) == 34 + source);
      TEST_ASSERT(bonds[2].atomId(1) == 35 + source);

      TEST_ASSERT(bonds[3].id() == 1);
      TEST_ASSERT(bonds[3].typeId() == 0);
      TEST_ASSERT(bonds[3].atomId(0) == 38 + source);
      TEST_ASSERT(bonds[3].atomId(1) == 39 + source);

      #if 0
      MpiLogger logger;
      logger.begin();
      int id, typeId;
      id     = bonds[2].id();
      typeId = bonds[2].typeId();
      pos = bonds[2].position();
      vel = bonds[2].velocity();

      std::cout << std::endl;
      std::cout << "proc: " << myrank 
                << ", source " << source << ", dest   " <<  dest
                << ", comm size " << commsize << std::endl
                << "proc: " << myrank 
                << ", id " << id << ", typeId " << typeId << std::endl
                << "proc: " << myrank << ", position = "
                << pos[0] << " " << pos[1] << " " << pos[2] << std::endl
                << "proc " << myrank 
                << "velocity = "
                <<  vel[0] << ", " << vel[1] << ", " << vel[2] << std::endl;


      id     = bonds[3].id();
      typeId = bonds[3].typeId();
      pos = bonds[3].position();
      vel = bonds[3].velocity();

      std::cout << "proc: " << myrank 
                << ", source " << source << ", dest   " <<  dest
                << ", comm size " << commsize << std::endl
                << "proc: " << myrank 
                << ", id " << id << ", typeId " << typeId << std::endl
                << "proc: " << myrank << ", position = "
                << pos[0] << " " << pos[1] << " " << pos[2] << std::endl
                << "proc " << myrank 
                << "velocity = "
                <<  vel[0] << ", " << vel[1] << ", " << vel[2] << std::endl;


      logger.end();
      #endif

   }

   // Test Method for MPI Sendrecv method for bonds
   void testBondBcast()
   {
      printMethod(TEST_FUNC);

      DArray<Bond> bonds;
      int          myrank, commsize, source;

      bonds.allocate(4);
      myrank   = MPI::COMM_WORLD.Get_rank();
      commsize = MPI::COMM_WORLD.Get_size();
      source   = 0;

      if (myrank == source) {

         // Fill one bond object. 
         bonds[0].setId(0);
         bonds[0].setTypeId(0);
         bonds[0].setAtomId(0, myrank + 34);
         bonds[0].setAtomId(1, myrank + 35);
   
         // Fill another bondobject
         bonds[1].setId(1);
         bonds[1].setTypeId(0);
         bonds[1].setAtomId(0, myrank + 38);
         bonds[1].setAtomId(1, myrank + 39);
   
         //Initialize the sendbuffer, set atomtype to GROUP
         buffer_.clearSendBuffer();
         buffer_.beginSendBlock(Buffer::GROUP2);
   
         // Pack 2 bonds into the send buffer
         //buffer_.packGroup(bonds[0]);
         //buffer_.packGroup(bonds[1]);
         bonds[0].pack(buffer_);
         bonds[1].pack(buffer_);
         buffer_.endSendBlock();

      }

      buffer_.bcast(MPI::COMM_WORLD, source);

      if (myrank != source) {

         // Unpack bonds
         buffer_.beginRecvBlock();
         TEST_ASSERT(buffer_.recvSize() == 2);
         bonds[2].unpack(buffer_);
         bonds[3].unpack(buffer_);
         TEST_ASSERT(buffer_.recvSize() == 0);
         buffer_.endRecvBlock();
   
         TEST_ASSERT(bonds[2].id() == 0);
         TEST_ASSERT(bonds[2].typeId() == 0);
         TEST_ASSERT(bonds[2].atomId(0) == 34 + source);
         TEST_ASSERT(bonds[2].atomId(1) == 35 + source);
   
         TEST_ASSERT(bonds[3].id() == 1);
         TEST_ASSERT(bonds[3].typeId() == 0);
         TEST_ASSERT(bonds[3].atomId(0) == 38 + source);
         TEST_ASSERT(bonds[3].atomId(1) == 39 + source);

      }

   }
   #endif

};

TEST_BEGIN(BufferTest)
TEST_ADD(BufferTest, testCapacities)
#ifdef UTIL_MPI
TEST_ADD(BufferTest, testPackAtom)
TEST_ADD(BufferTest, testPackGhost)
TEST_ADD(BufferTest, testAtomSendRecv)
TEST_ADD(BufferTest, testGhostSendRecv)
TEST_ADD(BufferTest, testAtomSend_Recv)
TEST_ADD(BufferTest, testGhostSend_Recv)
TEST_ADD(BufferTest, testAtomGhostSend_Recv)
TEST_ADD(BufferTest, testBondSendRecv)
TEST_ADD(BufferTest, testBondBcast)
#endif
TEST_END(BufferTest)

#endif /* BUFFER_TEST_H */
