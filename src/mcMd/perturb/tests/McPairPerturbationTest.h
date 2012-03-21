#ifndef MCMD_MC_PAIR_PERTURBATION_TEST_H
#define MCMD_MC_PAIR_PERTURBATION_TEST_H

#include <mcMd/perturb/mcSystem/McPairPerturbation.h>
#include <mcMd/mcSimulation/McSimulation.h>
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/mdSimulation/MdSystem.h>
#include <mcMd/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/format/Dbl.h>

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace McMd;

class McPairPerturbationTest : public ParamFileTest<McSimulation>
{

public:

   McPairPerturbationTest()
    : ParamFileTest<McSimulation>(),
      system_(object().system()),
      perturbation_(object().system())
   {}

   virtual void setUp()
   {  
      openFile("in/McSimulation"); 
      object().readParam(file());
   }

   void testReadParam();
   void testPairEnergy();
   void testSetParameter();

private:

   McSystem&          system_;
   McPairPerturbation perturbation_;

};


void McPairPerturbationTest::testReadParam()
{
   printMethod(TEST_FUNC);

   try {
      object().isValid();
   } catch (Exception e) {
      std::cout << e.message();
      TEST_ASSERT(0);
   }

   if (verbose() > 1) {
      std::cout << std::endl;
      object().writeParam(std::cout);
   }

}


void McPairPerturbationTest::testPairEnergy()
{ 
   printMethod(TEST_FUNC);
   std::cout << std::endl;

   object().run();

   System::MoleculeIterator molIter;
   Molecule::AtomIterator   atomIter;
   double energy, de;

   energy = 0.0;
   for (int is=0; is < object().nSpecies(); ++is) {
      for (system_.begin(is, molIter); molIter.notEnd(); ++molIter) {
         for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
            de = system_.atomPairEnergy(*atomIter);
            //std::cout.width(5);
            //std::cout << atomIter->id() << "     " << de << std::endl;
            energy += de;
         }
      }
   }
   std::cout << "Total atomPairEnergy = " << 0.5*energy << std::endl;
   std::cout << "TotalPairEnergy      = " << system_.pairEnergy() << std::endl;

}

void McPairPerturbationTest::testSetParameter()
{ 
   printMethod(TEST_FUNC);
   std::cout << std::endl;

   object().run();

   object().system().perturbation().setParameter(2.0);
   std::cout << Dbl(object().system().perturbation().parameter()) << std::endl;
   TEST_ASSERT( eq(object().system().perturbation().parameter(), 2.0) );
   std::cout << Dbl(object().system().perturbation().derivative()) << std::endl;

   //object().system().perturbation().setParameter(2.0);
   //std::cout << Dbl(object().system().perturbation().parameter()) << std::endl;
   //TEST_ASSERT( eq(object().system().perturbation().parameter(), 2.0) );
   //std::cout << Dbl(object().system().perturbation().derivative()) << std::endl;

   System::MoleculeIterator molIter;
   Molecule::AtomIterator   atomIter;
   double energy, de;

   energy = 0.0;
   for (int is=0; is < object().nSpecies(); ++is) {
      for (system_.begin(is, molIter); molIter.notEnd(); ++molIter) {
         for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
            de = system_.atomPairEnergy(*atomIter);
            //std::cout.width(5);
            //std::cout << atomIter->id() << "     " << de << std::endl;
            energy += de;
         }
      }
   }
   std::cout << "Total atomPairEnergy = " << 0.5*energy << std::endl;
   std::cout << "TotalPairEnergy      = " << system_.pairEnergy() << std::endl;

}



TEST_BEGIN(McPairPerturbationTest)
TEST_ADD(McPairPerturbationTest, testReadParam)
TEST_ADD(McPairPerturbationTest, testPairEnergy)
TEST_ADD(McPairPerturbationTest, testSetParameter)
TEST_END(McPairPerturbationTest)

#endif
