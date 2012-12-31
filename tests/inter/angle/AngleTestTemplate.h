#ifndef ANGLE_TEST_TEMPLATE_H
#define ANGLE_TEST_TEMPLATE_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <inter/angle/BendForce.h>

#include <fstream>
#include <string>

using namespace Util;
using namespace Inter;

template <class Interaction>
class AngleTestTemplate : public UnitTest 
{

public:

   AngleTestTemplate(const char* paramFileName)
    : paramFileName_(paramFileName)
   {}

   void setUp()
   {
      eps_ = 1.0E-5;
      nAngleType_ = 1;
      interaction_.setNAngleType(nAngleType_);
      std::ifstream paramFile;
      openInputFile(paramFileName_, paramFile);
      interaction_.readParam(paramFile);
      paramFile.close();
   }

   void tearDown()
   {}

protected:

   BendForce    bend_;
   Interaction  interaction_;
   Vector       b1_, b2_;
   double       eps_;
   int          nAngleType_;
   std::string  paramFileName_;

   void forceTest()
   {
      bend_.computeDerivatives(b1_, b2_);
      Vector b10 = b1_;
      Vector b20 = b2_;

      Vector f1, f2;
      int type = 0;
      double e0 = interaction_.energy(bend_.cosTheta, type);
      interaction_.force(b1_, b2_, f1, f2, type);

      // Derivative with respect to b1
      double d, e1, e2;
      //std::cout << std::endl;
      for (int i = 0; i < Dimension; ++i) {
         b1_ = b10;
         b2_ = b20;

         b1_[i] += eps_;
         bend_.computeAngle(b1_, b2_);
         e1 = interaction_.energy(bend_.cosTheta, type);

         b1_[i] += eps_;
         bend_.computeAngle(b1_, b2_);
         e2 = interaction_.energy(bend_.cosTheta, type);

         d = (4.0*e1 - e2 - 3.0*e0)/(2.0*eps_);
         //std::cout << f1[i] << "   " << d << "   "
         //          << d - f1[i] << std::endl;
         TEST_ASSERT(fabs(d - f1[i]) < 1.0E-8);
      }

      // Derivative with respect to b2
      for (int i = 0; i < Dimension; ++i) {
         b1_ = b10;
         b2_ = b20;

         b2_[i] += eps_;
         bend_.computeAngle(b1_, b2_);
         e1 = interaction_.energy(bend_.cosTheta, type);

         b2_[i] += eps_;
         bend_.computeAngle(b1_, b2_);
         e2 = interaction_.energy(bend_.cosTheta, type);

         d = (4.0*e1 - e2 - 3.0*e0)/(2.0*eps_);
         //std::cout << f2[i] << "   " << d << "   "
         //          << d - f2[i] << std::endl;
         TEST_ASSERT(fabs(d - f2[i]) < 1.0E-8);
      }
   }

};

#endif
