#ifndef MONOCLINIC_BOUNDARYMI_TEST
#define MONOCLINIC_BOUNDARYMI_TEST

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/boundary/MonoclinicBoundaryMI.h>
#include <util/boundary/serialize.h>
#include <util/space/Vector.h>
#include <util/space/IntVector.h>

#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>

#include <fstream>

using namespace Util;

class MonoclinicBoundaryMITest : public UnitTest 
{

private:

   MonoclinicBoundaryMI boundary;

public:

   void setUp()
   {};

   void tearDown()
   {};

   void testInitialize() 
   {
      printMethod(TEST_FUNC);

      Vector L, Lp;
      double d;

      L[0] = 2.0;
      L[1] = 3.0;
      L[2] = 4.0;
      d = 1.0;

      boundary.setOrthorhombic(L,d);
      Lp = boundary.lengths();


      // Assertions
      TEST_ASSERT(boundary.isValid());

      // Verbose output
/*      if (verbose() > 1) {
         printf("Boundary.minima_: %lf %lf %lf\n", 
                 boundary.minima_[0], boundary.minima_[1], boundary.minima_[2]);
         printf("Boundary.maxima_: %lf %lf %lf\n", 
                 boundary.maxima_[0], boundary.maxima_[1], boundary.maxima_[2]);
         printf("Boundary.L   tilt: %lf %lf %lf %lf\n", 
                 boundary.lengths_[0], boundary.lengths_[1], boundary.lengths_[2], tilt);
      }
*/
      std::cout << std::endl;
      std::cout << "BravaisBasis(1)   " << boundary.bravaisBasisVector(0) << std::endl;
      std::cout << "BravaisBasis(2)   " << boundary.bravaisBasisVector(1) << std::endl;
      std::cout << "BravaisBasis(3)   " << boundary.bravaisBasisVector(2) << std::endl;
      std::cout << "ReciprocalBasis(1)" << boundary.reciprocalBasisVector(0) << std::endl;
      std::cout << "ReciprocalBasis(2)" << boundary.reciprocalBasisVector(1) << std::endl;
      std::cout << "ReciprocalBasis(3)" << boundary.reciprocalBasisVector(2) << std::endl;

   }

   void testStreamIO1() 
   {
      printMethod(TEST_FUNC);

      int i;

      // Read parameters from file
      std::ifstream in;
      openInputFile("in/MonoclinicBoundary", in);
     
      in >> boundary;

      std::cout << std::endl;
      std::cout << boundary << std::endl;

      // Assertions
      TEST_ASSERT(boundary.isValid());
      for (i = 0; i < 3; i++) {
         TEST_ASSERT(eq(boundary.minima_[i], 0.0));
         TEST_ASSERT(eq(boundary.maxima_[i], boundary.lengths_[i]));
         TEST_ASSERT(boundary.lengths_[i] > 1.0E-8);
      }

      // Verbose output
      if (verbose() > 1) {
         std::cout << boundary << std::endl;
      }

   }

   void testStreamIO2() 
   {
      printMethod(TEST_FUNC);

      int i;

      // Read parameters from file
      std::ifstream in;
      openInputFile("in/MonoclinicBoundary", in);
     
      in >> boundary;

      std::cout << std::endl;
      std::cout << boundary << std::endl;

      // Assertions
      TEST_ASSERT(boundary.isValid());
      for (i = 0; i < 3; i++) {
         TEST_ASSERT(eq(boundary.minima_[i], 0.0));
         TEST_ASSERT(eq(boundary.maxima_[i], boundary.lengths_[i]));
         TEST_ASSERT(boundary.lengths_[i] > 1.0E-8);
      }

      // Verbose output
      if (verbose() > 1) {
         std::cout << boundary << std::endl;
      }

   }

   void testStreamIO3() 
   {
      printMethod(TEST_FUNC);

      int i;

      // Read parameters from file
      std::ifstream in;
      openInputFile("in/MonoclinicBoundary", in);
     
      in >> boundary;

      std::cout << std::endl;
      std::cout << boundary << std::endl;

      // Assertions
      for (i = 0; i < 3; i++) {
         TEST_ASSERT(eq(boundary.minima_[i], 0.0));
         TEST_ASSERT(eq(boundary.maxima_[i], boundary.lengths_[i]));
         TEST_ASSERT(boundary.lengths_[i] > 1.0E-8);
      }

      // Verbose output
      if (verbose() > 1) {
         std::cout << boundary << std::endl;
      }
   }
   
   void testSerialize() 
   {
      printMethod(TEST_FUNC);

      MemoryOArchive oar;
      MemoryIArchive iar;

      int i;

      // Read parameters from file
      std::ifstream in;
      openInputFile("in/MonoclinicBoundary", in);
     
      in >> boundary;
      oar.allocate(2000);
      oar << boundary;
      iar = oar;

      MonoclinicBoundaryMI clone;
      iar >> clone;

      std::cout << std::endl;
      std::cout << clone << std::endl;

      // Assertions
      TEST_ASSERT(boundary.isValid());
      for (i = 0; i < Dimension; i++) {
         TEST_ASSERT(eq(clone.minima_[i], 0.0));
         TEST_ASSERT(eq(clone.maxima_[i], clone.lengths_[i]));
         TEST_ASSERT(eq(clone.lengths_[i], boundary.lengths_[i]));
         TEST_ASSERT(eq(clone.volume(), boundary.volume()));
         TEST_ASSERT(clone.lengths_[i] > 1.0E-8);
      }

      // Verbose output
      if (verbose() > 1) {
         std::cout << boundary << std::endl;
      }

   }

   void testShift()
   {
      printMethod(TEST_FUNC);

      Vector R, L;
      double d;

      // Setup Boundary
      L[0] = 2.0;
      L[1] = 3.0;
      L[2] = 4.0;
      d = 1.0;
      boundary.setOrthorhombic(L,d);

      R[0] = 2.6;
      R[1] = -0.4;
      R[2] = 2.1;
      boundary.shift(R);
      TEST_ASSERT(eq(R[0], 0.6));
      TEST_ASSERT(eq(R[1], 2.6));
      TEST_ASSERT(eq(R[2], 3.1));

      R[0] = 1.6;
      R[1] = 2.4;
      R[2] = 4.1;
      boundary.shift(R);
      TEST_ASSERT(eq(R[0], 1.6));
      TEST_ASSERT(eq(R[1], 2.4));
      TEST_ASSERT(eq(R[2], 0.1));

      R[0] =  0.6;
      R[1] = -0.4;
      R[2] = -2.1;
      boundary.shift(R);
      TEST_ASSERT(eq(R[0], 0.6));
      TEST_ASSERT(eq(R[1], 2.6));
      TEST_ASSERT(eq(R[2], 2.9));

   };

   void testDistanceSq1()
   {
      printMethod(TEST_FUNC);

      Vector L;
      double d;

      Vector R1;
      Vector R2;
      double dRSq;

      // Setup Boundary
      L[0] = 2.0;
      L[1] = 3.0;
      L[2] = 4.0;
      d = 1.0;

      boundary.setOrthorhombic(L,d);

      R1[0] =  0.0;
      R1[1] =  0.0;
      R1[2] =  0.0;
      R2[0] =  0.0;
      R2[1] =  0.0;
      R2[2] =  2.1;
      dRSq = boundary.distanceSq(R1, R2);
      TEST_ASSERT(eq(dRSq, 3.61));

      R1[0] =  2.0;
      R1[1] =  0.0;
      R1[2] =  0.0;
      R2[0] =  1.1;
      R2[1] =  3.0;
      R2[2] =  1.0;
      dRSq = boundary.distanceSq(R1, R2);
      TEST_ASSERT(eq(dRSq, 0.81));

      R1[0] =  2.0;
      R1[1] =  3.0;
      R1[2] =  5.0;
      R2[0] =  0.0;
      R2[1] =  0.0;
      R2[2] =  0.1;
      dRSq = boundary.distanceSq(R1, R2);
      TEST_ASSERT(eq(dRSq, 0.01));

   };


   void testDistanceSq2()
   {
      Vector    L;
      double d;

      Vector    R1;
      Vector    R2;
      Vector    dR;
      IntVector shift; 
      double    dRSq1, dRSq2, dRSq3;

      printMethod(TEST_FUNC);

      // Setup Boundary
      L[0] = 2.0;
      L[1] = 3.0;
      L[2] = 4.0;
      d = 1.0;

      boundary.setOrthorhombic(L,d);

      R1[0] =  0.0;
      R1[1] =  0.0;
      R1[2] =  0.0;
      R2[0] =  0.0;
      R2[1] =  0.0;
      R2[2] =  2.1;
      dRSq1 = boundary.distanceSq(R1, R2);
      dRSq2 = boundary.distanceSq(R1, R2, dR);
      dRSq3 = boundary.distanceSq(R1, R2, shift);
 

      TEST_ASSERT(eq(dR[0], 0.0));
      TEST_ASSERT(eq(dR[1], 0.0));
      TEST_ASSERT(eq(dR[2], 1.9));
      TEST_ASSERT(eq(dRSq1, 3.61));
      TEST_ASSERT(eq(dRSq2, 3.61));
      TEST_ASSERT(eq(dRSq3, 3.61));

      R1[0] =  2.0;
      R1[1] =  0.0;
      R1[2] =  0.0;
      R2[0] =  1.1;
      R2[1] =  3.0;
      R2[2] =  1.0;
      dRSq1 = boundary.distanceSq(R1, R2);
      dRSq2 = boundary.distanceSq(R1, R2, dR);
      dRSq3 = boundary.distanceSq(R1, R2, shift);
      TEST_ASSERT(eq(dR[0], 0.9));
      TEST_ASSERT(eq(dR[1], 0.0));
      TEST_ASSERT(eq(dR[2], 0.0));
      TEST_ASSERT(eq(dRSq1, 0.81));
      TEST_ASSERT(eq(dRSq2, 0.81));
      TEST_ASSERT(eq(dRSq3, 0.81));

      R1[0] =  2.2;
      R1[1] =  3.0;
      R1[2] =  5.1;
      R2[0] =  0.2;
      R2[1] =  0.0;
      R2[2] =  0.1;
      dRSq1 = boundary.distanceSq(R1, R2);
      dRSq2 = boundary.distanceSq(R1, R2, dR);
      dRSq3 = boundary.distanceSq(R1, R2, shift);
      TEST_ASSERT(eq(dR[0], 0.0));
      TEST_ASSERT(eq(dR[1], 0.0));
      TEST_ASSERT(eq(dR[2], 0.0));
      TEST_ASSERT(eq(dRSq1, 0.0));
      TEST_ASSERT(eq(dRSq2, 0.0));
      TEST_ASSERT(eq(dRSq3, 0.0));

      R1[0] =  2.0;
      R1[1] = -0.1;
      R1[2] =  2.0;
      R2[0] =  0.0;
      R2[1] =  0.0;
      R2[2] =  0.0;
      dRSq1 = boundary.distanceSq(R1, R2);
      dRSq2 = boundary.distanceSq(R1, R2, dR);
      dRSq3 = boundary.distanceSq(R1, R2, shift);
      TEST_ASSERT(eq(dR[0], 0.0));
      TEST_ASSERT(eq(dR[1],-0.1));
      TEST_ASSERT(eq(dR[2], 2.0));

      R1[0] =  0.0;
      R1[1] =  1.5;
      R1[2] =  0.5;
      R2[0] =  0.0;
      R2[1] =  0.0;
      R2[2] =  0.0;
      dRSq1 = boundary.distanceSq(R1, R2);
      dRSq2 = boundary.distanceSq(R1, R2, dR);
      dRSq3 = boundary.distanceSq(R1, R2, shift);
      TEST_ASSERT(eq(dR[0], 0.0));
      TEST_ASSERT(eq(dR[1], 1.5));
      TEST_ASSERT(eq(dR[2], 0.5));

      R1[0] =  0.0;
      R1[1] =  0.0;
      R1[2] =  2.0;
      R2[0] =  0.0;
      R2[1] =  0.0;
      R2[2] =  0.0;
      dRSq1 = boundary.distanceSq(R1, R2);
      dRSq2 = boundary.distanceSq(R1, R2, dR);
      dRSq3 = boundary.distanceSq(R1, R2, shift);
      TEST_ASSERT(eq(dR[0], 0.0));
      TEST_ASSERT(eq(dR[1], 0.0));
      TEST_ASSERT(eq(dR[2], 2.0));

   };

}; 

TEST_BEGIN(MonoclinicBoundaryMITest)
TEST_ADD(MonoclinicBoundaryMITest, testInitialize)
TEST_ADD(MonoclinicBoundaryMITest, testStreamIO1)
TEST_ADD(MonoclinicBoundaryMITest, testStreamIO2)
TEST_ADD(MonoclinicBoundaryMITest, testStreamIO3)
TEST_ADD(MonoclinicBoundaryMITest, testSerialize)
TEST_ADD(MonoclinicBoundaryMITest, testShift)
TEST_ADD(MonoclinicBoundaryMITest, testDistanceSq1)
TEST_ADD(MonoclinicBoundaryMITest, testDistanceSq2)

TEST_END(MonoclinicBoundaryMITest)

#endif
