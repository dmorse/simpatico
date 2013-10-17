#ifndef DDMD_DOMAIN_TEST_H
#define DDMD_DOMAIN_TEST_H

#include <util/global.h>
#include <ddMd/communicate/Domain.h>
#include <util/boundary/Boundary.h>
#include <util/space/Grid.h>
#include <util/format/Int.h>
#include <util/mpi/MpiLogger.h>

#ifdef UTIL_MPI
#ifndef TEST_MPI
#define TEST_MPI
#endif
#endif

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <test/ParamFileTest.h>

#include <iostream>

using namespace Util;
using namespace DdMd;

class DomainTest : public ParamFileTest
{

private:

    Domain domain_;

public:

   virtual void setUp()
   {
      #ifdef UTIL_MPI  
      domain_.setGridCommunicator(communicator()); 
      domain_.setIoCommunicator(communicator()); 
      #endif
   }

   void testReadParam()
   {  
      printMethod(TEST_FUNC); 
   
      Grid grid;
      Boundary  boundary;
      IntVector dimensions;
      IntVector partner;

      #if UTIL_MPI
      openFile("in/Domain"); 
      #else
      openFile("in/Domain.111"); 
      #endif

      #ifdef UTIL_MPI
      #if 0
      MpiLogger logger;
      logger.begin();
      std::cout << std::endl;
      std::cout << "Processor " << mpiRank() 
                << ",  isIoProcessor " << isIoProcessor() 
                << ",  isFileOpen "    << file().is_open()
                << std::endl;
      logger.end();
      #endif
      #else
      int rank = 0;
      domain_.setRank(rank);
      #endif

      domain_.setBoundary(boundary);
      domain_.readParam(file()); 

      #if 0
      #ifdef UTIL_MPI
      logger.begin();
      #endif

      std::cout << "Rank = " << domain_.gridRank() << std::endl;

      int i, j, k, m;
      std::cout << "Dimensions  = ";
      for (int i = 0; i < Dimension; ++i) {
         dimensions[i] = domain_.grid().dimension(i);
         std::cout << dimensions[i] << "   ";
      } 
      std::cout << std::endl;
      grid.setDimensions(dimensions);

      std::cout << "Rank        = " << domain_.gridRank() << std::endl;
      std::cout << "Coordinates = ";
      for (int i = 0; i < Dimension; ++i) {
         std::cout << domain_.gridCoordinate(i) << "   ";
      } 
      std::cout << std::endl;

      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < 2; ++j) {
            std::cout << "   Direction   ";
            std::cout << i << "   " << j << "   ";

            m = domain_.destRank(i, j);
            partner = grid.position(m);
            std::cout << "   dest   ";
            for (k = 0; k < Dimension; ++k) {
               std::cout << partner[k] << "   ";
            }

            m = domain_.sourceRank(i, j);
            partner = grid.position(m);
            std::cout << "   source   ";
            for (k = 0; k < Dimension; ++k) {
               std::cout << partner[k] << "   ";
            }

            std::cout << "shift   =  " << domain_.shift(i, j);
            std::cout << std::endl;
         }
      }
      std::cout << std::endl;

      #ifdef UTIL_MPI
      logger.end();
      #endif

      #endif // if 0

   }

   #if UTIL_MPI
   void testPing()
   {  
      printMethod(TEST_FUNC); 
 
      Boundary boundary; 
      MPI::Request request[2]; 
      MpiLogger logger;

      domain_.setBoundary(boundary);

      openFile("in/Domain"); 
      domain_.readParam(file()); 

      //std::cout << std::endl;
      
      int i, j, d, s, myRank, rr;
      myRank = domain_.gridRank();
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < 2; ++j) {
            s = domain_.sourceRank(i, j);
            d = domain_.destRank(i, j);

            #if 0
            int size = domain_.communicator().Get_size();
            logger.begin();
            std::cout << Int(i) << Int(j);
            std::cout << Int(myRank) << Int(d) << Int(s) << std::endl;
            std::cout << "Size is "<<size<<std::endl;
            logger.end();
            #endif

            if (s != myRank) {
               request[0] = domain_.communicator().Irecv(&rr, 1, MPI::INT, s, 34);
            }
            if (d != myRank) {
               request[1] = domain_.communicator().Isend(&myRank, 1, MPI::INT, d, 34);
            }
            if (s != myRank) {
               request[0].Wait();
               TEST_ASSERT(rr == s);
            }
            if (d != myRank) {
               request[1].Wait();
            }
         }
      }

   }
   #endif

};

TEST_BEGIN(DomainTest)
TEST_ADD(DomainTest, testReadParam)
#ifdef UTIL_MPI
TEST_ADD(DomainTest, testPing)
#endif
TEST_END(DomainTest)

#endif
