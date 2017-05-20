#ifndef ANGLE_TEST_TEMPLATE_H
#define ANGLE_TEST_TEMPLATE_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <simp/interaction/angle/BendForce.h>

#include <fstream>
#include <string>

using namespace Util;
using namespace Simp;

template <class Interaction>
class AngleTestTemplate : public UnitTest 
{

protected:

   BendForce    bend_;
   Interaction  interaction_;
   Vector       b1_, b2_;
   double       eps_;
   int          nAngleType_;
   int          type_;

   AngleTestTemplate()
    : eps_(1.0E-5),
      nAngleType_(0),
      type_(0)
   {}

   void setNAngleType(int nAngleType)
   {
      nAngleType_ = nAngleType;
      interaction_.setNAngleType(nAngleType);
   }

   void readParamFile(std::string fileName)
   {
      std::ifstream paramFile;
      openInputFile(fileName, paramFile);
      interaction_.readParam(paramFile);
      paramFile.close();
   }

   void forceTest()
   {
      bend_.computeDerivatives(b1_, b2_);
      Vector b10 = b1_;
      Vector b20 = b2_;

      Vector f1, f2;
      double e0 = interaction_.energy(bend_.cosTheta, type_);
      interaction_.force(b1_, b2_, f1, f2, type_);

      // Derivative with respect to b1
      double d, e1, e2;
      //std::cout << std::endl;
      for (int i = 0; i < Dimension; ++i) {
         b1_ = b10;
         b2_ = b20;

         b1_[i] += eps_;
         bend_.computeAngle(b1_, b2_);
         e1 = interaction_.energy(bend_.cosTheta, type_);

         b1_[i] += eps_;
         bend_.computeAngle(b1_, b2_);
         e2 = interaction_.energy(bend_.cosTheta, type_);

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
         e1 = interaction_.energy(bend_.cosTheta, type_);

         b2_[i] += eps_;
         bend_.computeAngle(b1_, b2_);
         e2 = interaction_.energy(bend_.cosTheta, type_);

         d = (4.0*e1 - e2 - 3.0*e0)/(2.0*eps_);
         //std::cout << f2[i] << "   " << d << "   "
         //          << d - f2[i] << std::endl;
         TEST_ASSERT(fabs(d - f2[i]) < 1.0E-8);
      }
   }

};

#endif
