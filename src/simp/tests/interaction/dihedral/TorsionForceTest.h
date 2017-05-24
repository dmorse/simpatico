#ifndef TORSION_FORCE_TEST_H
#define TORSION_FORCE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <simp/interaction/dihedral/TorsionForce.h>

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
      b1_ = Vector(1.0, 0.0,  0.0);
      b2_ = Vector(0.0, 1.0,  0.0);
      b3_ = Vector(0.0, 0.0,  1.0);
      torsion_.computeAngle(b1_, b2_, b3_);
      angleTest();
      //std::cout << std::endl;
      //std::cout << torsion_.cosPhi << std::endl;
      TEST_ASSERT(eq(torsion_.cosPhi, 0.0));

      b1_ = Vector( 1.1,  0.2, -0.3);
      b2_ = Vector( 0.1,  0.9,  0.2);
      b3_ = Vector(-0.1,  0.4,  1.0);
      torsion_.computeDerivatives(b1_, b2_, b3_);
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

      b1_ = Vector( 1.1,  0.2, -0.3);
      b2_ = Vector( 0.1,  0.9,  0.2);
      b3_ = Vector(-0.1,  0.2,  1.0);
      torsion_.computeDerivatives(b1_, b2_, b3_);
      angleTest();
      derivativeTest();
   } 

   void angleTest()
   {
      double c = torsion_.cosPhi;
      double s = torsion_.sinPhi();
      double phi = torsion_.phi();
      TEST_ASSERT( eq(c, std::cos(phi)) );
      TEST_ASSERT( eq(s, std::sin(phi)) );
   }

   void derivativeTest()
   {
      Vector t;
      Vector d10 = torsion_.d1;
      Vector d20 = torsion_.d2;
      Vector d30 = torsion_.d3;
      double cos0 = torsion_.cosPhi;
      double d, cos1, cos2;

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
