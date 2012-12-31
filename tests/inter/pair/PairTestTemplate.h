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
using namespace Inter;

template <class Interaction>
class PairTestTemplate : public UnitTest 
{

protected:

   class Pair;

   int  nAtomType_;
   Interaction  interaction_;
   Pair  pair_;

   PairTestTemplate()
    : nAtomType_(0),
      interaction_(),
      pair_(interaction_)
   {}

   void readParamFile(const char* filename)
   {
      std::ifstream in;
      openInputFile(filename, in);
      interaction_.readParameters(in);
      in.close();
   }

   class Pair 
   {

   private:

      Interaction* interactionPtr_;
      double rsq_;
      int    type1_;
      int    type2_;

   public:

      Pair(Interaction& interaction)
       : interactionPtr_(&interaction),
         rsq_(1.0),
         type1_(0),
         type2_(0)
      {}

      Pair& set(double rsq, int type1, int type2)
      {
         rsq_   = rsq;
         type1_ = type1; 
         type2_ = type2; 
         return *this;
      }

      double energy()
      { return interactionPtr_->energy(rsq_, type1_, type2_); }

      double forceOverR()
      {  return interactionPtr_->forceOverR(rsq_, type1_, type2_); }

      bool testForce()
      {
         double fOverR  = forceOverR();
         double energy0 = energy();

         double eps  = 1.0E-8;
         double dRsq = eps*rsq_;

         double rsq1 = rsq_ + dRsq;
         double energy1 = interactionPtr_->energy(rsq1, type1_, type2_);

         double rsq2 = rsq1 + dRsq;
         double energy2 = interactionPtr_->energy(rsq2, type1_, type2_);

         double dEdRsq = (4.0*energy1 - energy2 - 3.0*energy0)/(2.0*dRsq);
         double diff = fOverR + 2.0*dEdRsq;
         diff = diff/(fOverR + 1.0);
         //std::cout << std::setprecision(12) 
         //          << fOverR      << "   " 
         //          << -2.0*dEdRsq << "   "
         //          << diff        << "   " << std::endl;
         return (fabs(diff) < 1.0E-5);
      }

   }; // class Pair

};

#endif
