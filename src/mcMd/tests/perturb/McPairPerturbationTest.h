#ifndef MCMD_MC_PAIR_PERTURBATION_TEST_H
#define MCMD_MC_PAIR_PERTURBATION_TEST_H

#include <mcMd/perturb/mcSystem/McPairPerturbation.h>
#include <mcMd/mcSimulation/McSimulation.h>
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/mdSimulation/MdSystem.h>
#include <simp/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/format/Dbl.h>

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace McMd;

class McPairPerturbationTest : public ParamFileTest
{

public:

   McPairPerturbationTest()
    : ParamFileTest(),
      system_(simulation_.system()),
      perturbation_(simulation_.system())
   {}

   virtual void setUp()
   {  
      openFile("in/McSimulation"); 
      simulation_.readParam(file());
   }

   void testReadParam();
   void testPairEnergy();
   void testSetParameter();

private:

   McSimulation       simulation_;
   McSystem&          system_;
   McPairPerturbation perturbation_;

};


void McPairPerturbationTest::testReadParam()
{
   printMethod(TEST_FUNC);

   try {
      simulation_.isValid();
   } catch (Exception e) {
      std::cout << e.message();
      TEST_ASSERT(0);
   }

   if (verbose() > 1) {
      std::cout << std::endl;
      simulation_.writeParam(std::cout);
   }

}


void McPairPerturbationTest::testPairEnergy()
{ 
   printMethod(TEST_FUNC);
   std::cout << std::endl;

   simulation_.run();

   System::MoleculeIterator molIter;
   Molecule::AtomIterator   atomIter;
   double energy, de;

   energy = 0.0;
   for (int is=0; is < simulation_.nSpecies(); ++is) {
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

   simulation_.run();

   simulation_.system().perturbation().setParameter(2.0);
   std::cout << Dbl(simulation_.system().perturbation().parameter()) << std::endl;
   TEST_ASSERT( eq(simulation_.system().perturbation().parameter(), 2.0) );
   std::cout << Dbl(simulation_.system().perturbation().derivative()) << std::endl;

   //simulation_.system().perturbation().setParameter(2.0);
   //std::cout << Dbl(simulation_.system().perturbation().parameter()) << std::endl;
   //TEST_ASSERT( eq(simulation_.system().perturbation().parameter(), 2.0) );
   //std::cout << Dbl(simulation_.system().perturbation().derivative()) << std::endl;

   System::MoleculeIterator molIter;
   Molecule::AtomIterator   atomIter;
   double energy, de;

   energy = 0.0;
   for (int is=0; is < simulation_.nSpecies(); ++is) {
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
