#ifndef MCMD_POINT_GENERATOR_TEST_H
#define MCMD_POINT_GENERATOR_TEST_H

#include "GeneratorTest.h"
#include <mcMd/generators/PointGenerator.h>

using namespace Util;
using namespace McMd;

class PointGeneratorTest : public GeneratorTest
{

public:

   // Test functions
   void testReadParamBond();

};

// Test methods

void PointGeneratorTest::testReadParamBond()
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

   Species& species = simulation_.species(0);
   PointGenerator generator(species, system_);

   system_.boundary().setCubic(10.0);
   DArray<double> diameters;

   diameters.allocate(2);
   diameters[0] = 1.0;
   diameters[1] = 1.0;

   CellList cellList;
   Generator::setupCellList(simulation_.atomCapacity(),
                            system_.boundary(),
                            diameters, cellList);
   
   TEST_ASSERT(generator.generate(500, diameters, cellList));

   std::ofstream outFile("out/point.cfg");
   system_.writeConfig(outFile);
   outFile.close();
}


TEST_BEGIN(PointGeneratorTest)
TEST_ADD(PointGeneratorTest, testReadParamBond)
TEST_END(PointGeneratorTest)

#endif
