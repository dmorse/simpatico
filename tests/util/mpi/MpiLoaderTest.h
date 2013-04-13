#ifndef MPI_LOADER_TEST_H
#define MPI_LOADER_TEST_H

#include <util/mpi/MpiFileIo.h>
#include <util/mpi/MpiLoader.h>
//#include <util/archives/Serializable.h>
#include <util/archives/BinaryFileOArchive.h>
#include <util/archives/BinaryFileIArchive.h>
#include <util/space/Vector.h>
#include <util/containers/DArray.h>
#include <util/containers/FArray.h>

#ifndef TEST_MPI
#define TEST_MPI
#endif

#include <mpi.h>
#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <iostream>
#include <fstream>
#include <complex>

using namespace Util;

class MpiLoaderTest : public UnitTest
{

public:

   MpiLoaderTest();
   // void setUp() {}
   // void tearDown() {}
   void testSetCommunicator();
   void testIsIoProcessor1();
   void testIsIoProcessor2();
   void testOArchiveConstructor1();
   void testOArchiveConstructor2();
   void testPack();

private:

   MpiFileIo  fileIo_;
   MpiLoader<BinaryFileIArchive> loader_;

};

MpiLoaderTest::MpiLoaderTest()
 : UnitTest(),
   fileIo_(),
   loader_(fileIo_)
{}

// void setUp() {}
// void tearDown() {}

void MpiLoaderTest::testSetCommunicator() 
{
   printMethod(TEST_FUNC);
   TEST_ASSERT(!loader_.hasIoCommunicator());
   loader_.setIoCommunicator(communicator());
   TEST_ASSERT(loader_.hasIoCommunicator());
   TEST_ASSERT(&loader_.ioCommunicator() == &communicator());
   loader_.clearCommunicator();
   TEST_ASSERT(!loader_.hasIoCommunicator());
}

void MpiLoaderTest::testIsIoProcessor1() 
{
   printMethod(TEST_FUNC);
   if (mpiRank() == 0) {
      TEST_ASSERT(loader_.isIoProcessor());
   } else
   if (mpiRank() == 1) {
      TEST_ASSERT(loader_.isIoProcessor());
   }
}

void MpiLoaderTest::testIsIoProcessor2() 
{
   printMethod(TEST_FUNC);
   loader_.setIoCommunicator(communicator());
   if (mpiRank() == 0) {
      TEST_ASSERT(loader_.isIoProcessor());
   } else
   if (mpiRank() == 1) {
      TEST_ASSERT(!loader_.isIoProcessor());
   }
}

void MpiLoaderTest::testOArchiveConstructor1()
{
   printMethod(TEST_FUNC);
   BinaryFileOArchive  v;
   if (isIoProcessor()) {
      openOutputFile("binary", v.file());
      v.file().close();
   }
} 

void MpiLoaderTest::testOArchiveConstructor2()
{
   printMethod(TEST_FUNC);
   if (isIoProcessor()) {
     BinaryFileOArchive  v("dummy");
     v.file().close();
   }
} 

void MpiLoaderTest::testPack()
{
   printMethod(TEST_FUNC);
   loader_.setIoCommunicator(communicator());

   // Declare variables
   int i1, i2;
   double d1, d2;
   std::complex<double> c1, c2;
   std::string s1, s2;
   Vector a1, a2;
   double b1[4];
   double b2[4];
   double m1[2][2];
   double m2[2][2];
   DArray<double> e1;
   DArray<double> e2;
   e1.allocate(4);
   e2.allocate(4);
   FArray<double, 4> f1;
   FArray<double, 4> f2;
   DMatrix<double> g1;
   DMatrix<double> g2;
   g1.allocate(2, 2);
   g2.allocate(2, 2);

   // Initialize variables
   i1 = 3;
   d1 = 45.0;
   c1 = std::complex<double>(3.0, 4.0);
   s1 = "My string has spaces";
   a1[0] =  1.3;
   a1[1] = -2.3;
   a1[2] =  3.3;

   b1[0] = 9.0;
   b1[1] = 8.0;
   b1[2] = 7.0;
   b1[3] = 6.0;

   m1[0][0] = 13.0;
   m1[0][1] = 14.0;
   m1[1][0] = 15.0;
   m1[1][1] = 16.0;

   e1[0] = 9.0;
   e1[1] = 3.0;
   e1[2] = 7.0;
   e1[3] = 6.0;

   f1[0] = 9.0;
   f1[1] = 3.0;
   f1[2] = 7.0;
   f1[3] = 6.0;
  
   g1(0, 0) = 12.0;
   g1(0, 1) = 14.0;
   g1(1, 0) = 19.0;
   g1(1, 1) = 16.0;

   if (isIoProcessor()) {
      BinaryFileOArchive  v;
      openOutputFile("binary", v.file());
  
      // Write variables to OArchive v
      v << i1;
      v & d1;
      v << s1;
      v << a1;
      v.pack(b1, 4);
      v.pack(m1[0], 2, 2);
      v << e1;
      v << f1;
      v << g1;
      //v << c1;

      v.file().close();
   }

   // Create IArchive u
   BinaryFileIArchive u;
   if (isIoProcessor()) {
      openInputFile("binary", u.file());
   }

   loader_.load(u, i2);  // int
   TEST_ASSERT(i1 == i2);

   loader_.load(u, d2);  // double
   TEST_ASSERT(d1 == d2);

   loader_.load(u, s2);   // string
   TEST_ASSERT(s1 == s2);

   loader_.load(u, a2);   // Vector
   TEST_ASSERT(a1 == a2);

   loader_.load(u, b2, 4);    // double C array
   for (int j = 0; j < 4; ++j) {
      TEST_ASSERT(b1[j] == b2[j]);
   }

   loader_.load(u, m2[0], 2, 2); // double 2D C array 
   int i, j;
   for (i = 0; i < 2; ++i) {
      for (j = 0; j < 2; ++j) {
         TEST_ASSERT(eq(m1[i][j], m2[i][j]));
      }
   }

   loader_.load(u, e2, 4);  // DArray<double>
   for (int j = 0; j < 4; ++j) {
      TEST_ASSERT(e1[j] == e2[j]);
   }

   loader_.load(u, f2);  // FArray<double>
   for (int j = 0; j < 4; ++j) {
      TEST_ASSERT(f1[j] == f2[j]);
   }

   loader_.load(u, g2, 2, 2); // DMatrix<double>
   for (i = 0; i < 2; ++i) {
      for (j = 0; j < 2; ++j) {
         TEST_ASSERT(eq(g1(i, j), g2(i, j)));
      }
   }

   //loader_.load(u, c2);   // complex
   //TEST_ASSERT(c1 == c2);

}

TEST_BEGIN(MpiLoaderTest)
TEST_ADD(MpiLoaderTest, testSetCommunicator)
TEST_ADD(MpiLoaderTest, testIsIoProcessor1)
TEST_ADD(MpiLoaderTest, testIsIoProcessor2)
TEST_ADD(MpiLoaderTest, testOArchiveConstructor1)
TEST_ADD(MpiLoaderTest, testOArchiveConstructor2)
TEST_ADD(MpiLoaderTest, testPack)
TEST_END(MpiLoaderTest)

#endif
