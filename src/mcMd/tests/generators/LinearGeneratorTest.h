#ifndef MCMD_LINEAR_GENERATOR_TEST_H
#define MCMD_LINEAR_GENERATOR_TEST_H

#include "GeneratorTest.h"
#include <mcMd/generators/LinearGenerator.h>

using namespace Util;
using namespace McMd;

class LinearGeneratorTest : public GeneratorTest
{

public:

   // Test functions
   void testReadParamBond();

};

// Test methods

void LinearGeneratorTest::testReadParamBond()
{
   printMethod(TEST_FUNC);
   readParam("in/McSimulation");

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

   Species& species = simulation_.species(1);
   LinearGenerator generator(species, system_);
   generator.setBondPotential(system_.bondPotential());

   system_.boundary().setCubic(10.0);
   DArray<double> diameters;

   diameters.allocate(2);
   diameters[0] = 1.0;
   diameters[1] = 1.0;

   CellList cellList;
   Generator::setupCellList(simulation_.atomCapacity(),
                            system_.boundary(),
                            diameters, cellList);
   
   
   bool success = generator.generate(100, diameters, cellList);
   TEST_ASSERT(success);

   std::ofstream outFile("out/linear.cfg");
   system_.writeConfig(outFile);
   outFile.close();
}


TEST_BEGIN(LinearGeneratorTest)
TEST_ADD(LinearGeneratorTest, testReadParamBond)
TEST_END(LinearGeneratorTest)

#endif
