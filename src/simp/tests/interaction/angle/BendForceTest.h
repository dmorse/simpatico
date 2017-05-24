#ifndef BEND_FORCE_TEST_H
#define BEND_FORCE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <simp/interaction/angle/BendForce.h>

using namespace Util;
using namespace Simp;

class BendForceTest : public UnitTest 
{

private:

   BendForce  angle_;
   Vector  b1_, b2_;
   double  eps_;

public:

   void setUp()
   { eps_ = 1.0E-6; }

   void testComputeAngle()
   {
      printMethod(TEST_FUNC);
      b1_ = Vector( 0.2, 0.90,  0.3);
      b2_ = Vector(-0.1, 0.85, -0.4);
      angle_.computeAngle(b1_, b2_);
      angleTest();

      double product = -0.2*0.1 + 0.90*0.85 - 0.3*0.4;
      double b1Abs = sqrt(0.2*0.2 + 0.90*0.90 + 0.3*0.3);
      double b2Abs = sqrt(0.1*0.1 + 0.85*0.85 + 0.4*0.4);
      double cos   = product/(b1Abs*b2Abs);
      // std::cout << std::endl;
      // std::cout << angle_.cosTheta << std::endl;
      // std::cout << cos << std::endl;
      TEST_ASSERT(eq(angle_.cosTheta, cos));

   } 


   void testComputeDerivatives()
   {
      printMethod(TEST_FUNC);
      b1_ = Vector( 0.2, 0.90,  0.3);
      b2_ = Vector(-0.1, 0.85, -0.4);
      angle_.computeDerivatives(b1_, b2_);
      angleTest();
      derivativeTest();
   } 

   void angleTest()
   {
      double product, b1Abs, b2Abs, c;
      product = b1_[0]*b2_[0] + b1_[1]*b2_[1] + b1_[2]*b2_[2];
      b1Abs = sqrt(b1_[0]*b1_[0] + b1_[1]*b1_[1] + b1_[2]*b1_[2]);
      b2Abs = sqrt(b2_[0]*b2_[0] + b2_[1]*b2_[1] + b2_[2]*b2_[2]);
      c = product/(b1Abs*b2Abs);

      TEST_ASSERT(eq(angle_.cosTheta, c));

      c = angle_.cosTheta;
      double s = angle_.sinTheta();
      double theta = angle_.theta();
      TEST_ASSERT( eq(c, std::cos(theta)) );
      TEST_ASSERT( eq(s, std::sin(theta)) );
   }

   void derivativeTest()
   {
      double cos0 = angle_.cosTheta;
      Vector b10 = b1_;
      Vector b20 = b2_;
      Vector d10 = angle_.d1;
      Vector d20 = angle_.d2;
      double d, cos1, cos2;

      // Derivative with respect to b1
      // std::cout << std::endl;
      for (int i = 0; i < Dimension; ++i) {
         b1_ = b10;
         b2_ = b20;
         b1_[i] += eps_;
         angle_.computeAngle(b1_, b2_);
         cos1 = angle_.cosTheta;
         b1_[i] += eps_;
         angle_.computeAngle(b1_, b2_);
         cos2 = angle_.cosTheta;
         d = (4.0*cos1 - cos2 - 3.0*cos0)/(2.0*eps_);
         // std::cout << d10[i] << "   " << d << "   "
         //          << d - d10[i] << std::endl;
         TEST_ASSERT(fabs(d10[i] - d) < 1.0E-8);
      }

      // Derivative with respect to b2
      for (int i = 0; i < Dimension; ++i) {
         b1_ = b10;
         b2_ = b20;
         b2_[i] += eps_;
         angle_.computeAngle(b1_, b2_);
         cos1 = angle_.cosTheta;
         b2_[i] += eps_;
         angle_.computeAngle(b1_, b2_);
         cos2 = angle_.cosTheta;
         d = (4.0*cos1 - cos2 - 3.0*cos0)/(2.0*eps_);
         // std::cout << d20[i] << "   " << d << "   "
         //           << d - d20[i] << std::endl;
         TEST_ASSERT(fabs(d20[i] - d) < 1.0E-8);
      }

   }

};

TEST_BEGIN(BendForceTest)
TEST_ADD(BendForceTest, testComputeAngle)
TEST_ADD(BendForceTest, testComputeDerivatives)
TEST_END(BendForceTest)

#endif
