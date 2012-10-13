#ifndef PAIR_TEST_H
#define PAIR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/containers/RArray.h>
#include <util/archives/Serializable.h>
#include <util/archives/Serializable_includes.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Inter;

template <class Interaction>
class PairTest : public UnitTest 
{

protected:

   int          nAtomType_;
   Interaction  interaction_;

   class ParameterSet 
   {

   private:

      Interaction* interactionPtr_;
      double rsq_;
      int    type1_;
      int    type2_;

   public:

      ParameterSet()
       : interactionPtr_(0),
         rsq_(1.0),
         type1_(0),
         type2_(0)
      {}

      ParameterSet(Interaction& interaction)
       : interactionPtr_(&interaction),
         rsq_(1.0),
         type1_(0),
         type2_(0)
      {}

      ParameterSet& set(Interaction& interaction, double rsq, int type1, int type2)
      {
         interactionPtr_ = &interaction;
         rsq_   = rsq;
         type1_ = type1; 
         type2_ = type2; 
         return *this;
      }

      ParameterSet& set(double rsq, int type1, int type2)
      {
         rsq_   = rsq;
         type1_ = type1; 
         type2_ = type2; 
         return *this;
      }

      double energy()
      { 
         return interactionPtr_->energy(rsq_, type1_, type2_); 
      }

      double forceOverR()
      {  return interactionPtr_->forceOverR(rsq_, type1_, type2_); }

      bool testForce()
      {
         double eps  = 1.0E-7;
         double dRsq = eps*rsq_;
         double energy0 = interactionPtr_->energy(rsq_, type1_, type2_);
         double energy1 = interactionPtr_->energy(rsq_ + dRsq, type1_, type2_);
         double fOverR  = forceOverR();
         //std::cout << precision(12) << fOverR << " " << -2.0*(energy1 - energy0)/dRsq << std::endl;
         double diff    = fOverR + 2.0*(energy1 - energy0)/dRsq;
         diff    = diff/(fOverR + 1.0);
         return (fabs(diff) < 1.0E-5);
      }

   };

public:

   void readParamFile(const char* filename)
   {
      std::ifstream in;
      openInputFile(filename, in);
      interaction_.readParameters(in);
      in.close();
   }

   void testSetUp() {
      printMethod(TEST_FUNC);
   }

   void testWrite() {

      printMethod(TEST_FUNC);

      // Verbose output
      if (verbose() > 0) {
         std::cout << std::endl; 
         interaction_.writeParam(std::cout);
      }

   }

};

#endif
