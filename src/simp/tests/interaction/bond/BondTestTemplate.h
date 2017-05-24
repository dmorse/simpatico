#ifndef BOND_TEST_TEMPLATE_H
#define BOND_TEST_TEMPLATE_H

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
class BondTestTemplate : public UnitTest
{

protected:

   Interaction  interaction_;
   double  rsq_;
   double  eps_;
   int  nBondType_;
   int  type_;

   BondTestTemplate()
    : interaction_(),
      rsq_(1.0),
      eps_(1.0E-5),
      nBondType_(0),
      type_(0)
   {}

   void setNBondType(int nBondType)
   {
      nBondType_ = nBondType;
      interaction_.setNBondType(nBondType);
   }

   void readParamFile(std::string filename)
   {
      std::ifstream in;
      openInputFile(filename, in);
      interaction_.readParameters(in);
      in.close();
   }

   double energy(double rsq, int type)
   {  return interaction_.energy(rsq, type); }

   double energy()
   {  return interaction_.energy(rsq_, type_); }

   double forceOverR(double rsq, int type)
   {  return interaction_.forceOverR(rsq, type); }

   double forceOverR()
   {  return interaction_.forceOverR(rsq_, type_); }

   bool testForce()
   {
      double fOverR  = interaction_.forceOverR(rsq_, type_);
      double energy0 = interaction_.energy(rsq_, type_);

      double dRsq = eps_*rsq_;
      double rsq1 = rsq_ + dRsq;
      double energy1 = interaction_.energy(rsq1, type_);

      double rsq2 = rsq1 + dRsq;
      double energy2 = interaction_.energy(rsq2, type_);

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
