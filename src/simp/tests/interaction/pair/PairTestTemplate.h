#ifndef PAIR_TEST_TEMPLATE_H
#define PAIR_TEST_TEMPLATE_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/containers/RArray.h>
#include <util/archives/Serializable.h>
#include <util/archives/Serializable_includes.h>

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace Util;
using namespace Simp;

template <class Interaction>
class PairTestTemplate : public UnitTest
{

protected:

   Interaction  interaction_;
   double  rsq_;
   double  eps_;
   int  nAtomType_;
   int  type1_;
   int  type2_;

   PairTestTemplate()
    : interaction_(),
      rsq_(1.0),
      eps_(1.0E-5),
      nAtomType_(0),
      type1_(0),
      type2_(0)
   {}

   void setNAtomType(int nAtomType)
   {
      nAtomType_ = nAtomType;
      interaction_.setNAtomType(nAtomType);
   }

   void readParamFile(std::string filename)
   {
      std::ifstream in;
      openInputFile(filename, in);
      interaction_.readParameters(in);
      in.close();
   }

   double energy(double rsq, int type1, int type2)
   {  return interaction_.energy(rsq, type1, type2); }

   double energy()
   {  return interaction_.energy(rsq_, type1_, type2_); }

   double forceOverR(double rsq, int type1, int type2)
   {  return interaction_.forceOverR(rsq, type1, type2); }

   double forceOverR()
   {  return interaction_.forceOverR(rsq_, type1_, type2_); }

   bool testForce()
   {
      double fOverR  = interaction_.forceOverR(rsq_, type1_, type2_);
      double energy0 = interaction_.energy(rsq_, type1_, type2_);

      double dRsq = eps_*rsq_;
      double rsq1 = rsq_ + dRsq;
      double energy1 = interaction_.energy(rsq1, type1_, type2_);

      double rsq2 = rsq1 + dRsq;
      double energy2 = interaction_.energy(rsq2, type1_, type2_);

      double dEdRsq = (4.0*energy1 - energy2 - 3.0*energy0)/(2.0*dRsq);
      double diff = fOverR + 2.0*dEdRsq;
      diff = diff/(fOverR + 1.0);
      //std::cout << std::setprecision(12)
      //          << fOverR      << "   "
      //          << -2.0*dEdRsq << "   "
      //          << diff        << "   " << std::endl;
      return (fabs(diff) < 1.0E-6);
   }

};

#endif
