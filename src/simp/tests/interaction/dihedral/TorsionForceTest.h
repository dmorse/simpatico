#ifndef TORSION_FORCE_TEST_H
#define TORSION_FORCE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <simp/interaction/dihedral/TorsionForce.h>
#include <cmath>

using namespace Util;
using namespace Simp;

class TorsionForceTest : public UnitTest 
{

private:

   TorsionForce  torsion_;
   Vector  b1_, b2_, b3_;
   double  eps_;

public:

   void setUp()
   { eps_ = 1.0E-6; }

   void testComputeAngle()
   {
      printMethod(TEST_FUNC);

      // Perpendicular crank (cosPhi = 0)
      b1_ = Vector(1.0, 0.0,  0.0);
      b2_ = Vector(0.0, 1.0,  0.0);
      b3_ = Vector(0.0, 0.0,  1.0);
      torsion_.computeAngle(b1_, b2_, b3_);
      angleTest();
      //std::cout << std::endl;
      //std::cout << torsion_.cosPhi << std::endl;
      TEST_ASSERT(eq(torsion_.cosPhi, 0.0));

      // Syn eclipsed (arc) configuration (cosPhi = 1)
      b1_ = Vector(1.0, 0.0,  0.0);
      b2_ = Vector(0.0, 1.0,  0.0);
      b3_ = Vector(-1.0, 0.0, 0.0);
      torsion_.computeAngle(b1_, b2_, b3_);
      angleTest();
      // std::cout << std::endl;
      // std::cout << torsion_.cosPhi << std::endl;
      TEST_ASSERT(eq(torsion_.cosPhi, 1.0));

      // Syn 45 deg twist configuration (cosPhi = 1/sqrt(2) )
      b1_ = Vector(1.0, 0.0,  0.0);
      b2_ = Vector(0.0, 1.0,  0.0);
      b3_ = Vector(-1.0, 0.0, 1.0);
      torsion_.computeAngle(b1_, b2_, b3_);
      angleTest();
      // std::cout << std::endl;
      // std::cout << torsion_.cosPhi << std::endl;
      TEST_ASSERT(eq(torsion_.cosPhi, 1.0/sqrt(2.0)));

      // Anti/trans zig-zag configuration (cosPhi = -1)
      b1_ = Vector(1.0, 0.0,  0.0);
      b2_ = Vector(0.0, 1.0,  0.0);
      b3_ = Vector(1.0, 0.0, 0.0);
      torsion_.computeAngle(b1_, b2_, b3_);
      angleTest();
      //std::cout << std::endl;
      //std::cout << torsion_.cosPhi << std::endl;
      TEST_ASSERT(eq(torsion_.cosPhi, -1.0));

      // Anti/trans 45 deg twist configuration (cosPhi = 1/sqrt(2) )
      b1_ = Vector(1.0, 0.0, 0.0);
      b2_ = Vector(0.0, 1.0, 0.0);
      b3_ = Vector(1.0, 0.0, 1.0);
      torsion_.computeAngle(b1_, b2_, b3_);
      angleTest();
      // std::cout << std::endl;
      // std::cout << torsion_.cosPhi << std::endl;
      TEST_ASSERT(eq(torsion_.cosPhi, -1.0/sqrt(2.0)));

      // Rescaled anti 45 deg twist configuration (cosPhi = 1/sqrt(2) )
      // Testing that angles don't change when vectors are scaled
      b1_ = Vector(2.0, 0.0, 0.0);
      b2_ = Vector(0.0, 1.5, 0.0);
      b3_ = Vector(1.3, 0.0, -1.3);
      torsion_.computeAngle(b1_, b2_, b3_);
      angleTest();
      // std::cout << std::endl;
      // std::cout << torsion_.cosPhi << std::endl;
      TEST_ASSERT(eq(torsion_.cosPhi, -1.0/sqrt(2.0)));

      // Irregular angle
      b1_ = Vector( 1.1,  0.4, -0.3);
      b2_ = Vector( 0.1,  0.9,  0.2);
      b3_ = Vector(-0.2,  0.4,  0.8);
      torsion_.computeAngle(b1_, b2_, b3_);
      angleTest();
      //std::cout << torsion_.cosPhi << std::endl;
   } 

   void testComputeDerivatives()
   {
      printMethod(TEST_FUNC);

      b1_ = Vector(1.0, 0.0,  0.0);
      b2_ = Vector(0.0, 1.0,  0.0);
      b3_ = Vector(0.0, 0.0,  1.0);
      torsion_.computeDerivatives(b1_, b2_, b3_);
      angleTest();
      derivativeTest();

      b1_ = Vector( 1.1,  0.4, -0.3);
      b2_ = Vector( 0.2,  0.9,  0.3);
      b3_ = Vector(-0.3,  0.2,  1.0);
      torsion_.computeDerivatives(b1_, b2_, b3_);
      angleTest();
      derivativeTest();
   } 

   void angleTest()
   {
      // Precondition: CosPhi was computed by calling either
      // torsion_.computeAngle or torsion_computeDerivatives 

      double c = torsion_.cosPhi;
      double s = torsion_.sinPhi();
      double phi = torsion_.phi();
      TEST_ASSERT( c <= 1.00000000001 && c >= -1.000000000001 );
      TEST_ASSERT( s >= 0.0 && s <= 1.0 );
      TEST_ASSERT( eq(s, std::sin(phi)) );
      TEST_ASSERT( eq(c, std::cos(phi)) );
      TEST_ASSERT( eq(s, std::sin(phi)) );
   }

   void derivativeTest()
   {
      // Preconditions:
      // Vectors b1_, b2_, b3_ were set 
      // torsion_.computeDerivatives(b1_, b2_, b3_) was called

      Vector t;
      Vector d10 = torsion_.d1;
      Vector d20 = torsion_.d2;
      Vector d30 = torsion_.d3;
      double cos0 = torsion_.cosPhi;
      double d, cos1, cos2;

      // Check that derivatives are perpendicular to force vectors
      double dot1 = d10.dot(b1_);
      double dot2 = d20.dot(b2_);
      double dot3 = d30.dot(b3_);
      TEST_ASSERT(eq(dot1, 0.0));
      TEST_ASSERT(eq(dot2, 0.0));
      TEST_ASSERT(eq(dot3, 0.0));

      // Derivative with respect to b1
      //std::cout << std::endl;
      for (int i = 0; i < Dimension; ++i) {
         t = b1_;
         t[i] += eps_;
         torsion_.computeAngle(t, b2_, b3_);
         cos1 = torsion_.cosPhi;
         t[i] += eps_;
         torsion_.computeAngle(t, b2_, b3_);
         cos2 = torsion_.cosPhi;
         d = (4.0*cos1 - cos2 - 3.0*cos0)/(2.0*eps_);
         //std::cout << d10[i] << "   " << d << "   "
         //         << d - d10[i] << std::endl;
         TEST_ASSERT(fabs(d10[i] - d) < 1.0E-8);
      }

      // Derivative with respect to b2
      for (int i = 0; i < Dimension; ++i) {
         t = b2_;
         t[i] += eps_;
         torsion_.computeAngle(b1_, t, b3_);
         cos1 = torsion_.cosPhi;
         t[i] += eps_;
         torsion_.computeAngle(b1_, t, b3_);
         cos2 = torsion_.cosPhi;
         d = (4.0*cos1 - cos2 - 3.0*cos0)/(2.0*eps_);
         // std::cout << d20[i] << "   " << d << "   "
         //           << d - d20[i] << std::endl;
         TEST_ASSERT(fabs(d20[i] - d) < 1.0E-8);
      }

      // Derivative with respect to b3
      for (int i = 0; i < Dimension; ++i) {
         t = b3_;
         t[i] += eps_;
         torsion_.computeAngle(b1_, b2_, t);
         cos1 = torsion_.cosPhi;
         t[i] += eps_;
         torsion_.computeAngle(b1_, b2_, t);
         cos2 = torsion_.cosPhi;
         d = (4.0*cos1 - cos2 - 3.0*cos0)/(2.0*eps_);
         //std::cout << d30[i] << "   " << d << "   "
         //          << d - d30[i] << std::endl;
         TEST_ASSERT(fabs(d30[i] - d) < 1.0E-8);
      }

   }

};

TEST_BEGIN(TorsionForceTest)
TEST_ADD(TorsionForceTest, testComputeAngle)
TEST_ADD(TorsionForceTest, testComputeDerivatives)
TEST_END(TorsionForceTest)

#endif
