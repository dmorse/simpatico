#ifndef DIHEDRAL_TEST_TEMPLATE_H
#define DIHEDRAL_TEST_TEMPLATE_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <simp/interaction/dihedral/TorsionForce.h>

#include <fstream>
#include <string>

using namespace Util;
using namespace Simp;

template <class Interaction>
class DihedralTestTemplate : public UnitTest 
{

protected:

   Interaction  interaction_;
   Vector       b1_, b2_, b3_;
   double       eps_;
   int          nDihedralType_;
   int          type_;

   DihedralTestTemplate()
    : eps_(1.0E-5),
      nDihedralType_(0),
      type_(0)
   {}

   void setNDihedralType(int nDihedralType)
   {
      nDihedralType_ = nDihedralType;
      interaction_.setNDihedralType(nDihedralType_);
   }

   void readParamFile(std::string paramFileName)
   {
      std::ifstream paramFile;
      openInputFile(paramFileName, paramFile);
      interaction_.readParam(paramFile);
      paramFile.close();
   }

   void forceTest()
   {
      Vector f1, f2, f3, t;
      double d, e0, e1, e2;

      e0 = interaction_.energy(b1_, b2_, b3_, type_);
      interaction_.force(b1_, b2_, b3_, f1, f2, f3, type_);

      // Derivative with respect to b1
      // std::cout << std::endl;
      for (int i = 0; i < Dimension; ++i) {
         t = b1_;
         t[i] += eps_;
         e1 = interaction_.energy(t, b2_, b3_, type_);
         t[i] += eps_;
         e2 = interaction_.energy(t, b2_, b3_, type_);
         d = (4.0*e1 - e2 - 3.0*e0)/(2.0*eps_);
         // std::cout << f1[i] << "   " << d << "   "
         //           << d - f1[i] << std::endl;
         TEST_ASSERT(fabs(d - f1[i]) < 1.0E-8);
      }

      // Derivative with respect to b2
      for (int i = 0; i < Dimension; ++i) {
         t = b2_;
         t[i] += eps_;
         e1 = interaction_.energy(b1_, t, b3_, type_);
         t[i] += eps_;
         e2 = interaction_.energy(b1_, t, b3_, type_);
         d = (4.0*e1 - e2 - 3.0*e0)/(2.0*eps_);
         // std::cout << f2[i] << "   " << d << "   "
         //           << d - f2[i] << std::endl;
         TEST_ASSERT(fabs(d - f2[i]) < 1.0E-8);
      }

      // Derivative with respect to b3
      for (int i = 0; i < Dimension; ++i) {
         t = b3_;
         t[i] += eps_;
         e1 = interaction_.energy(b1_, b2_, t, type_);
         t[i] += eps_;
         e2 = interaction_.energy(b1_, b2_, t, type_);
         d = (4.0*e1 - e2 - 3.0*e0)/(2.0*eps_);
         // std::cout << f3[i] << "   " << d << "   "
         //           << d - f3[i] << std::endl;
         TEST_ASSERT(fabs(d - f3[i]) < 1.0E-8);
      }

   }

};

#endif
