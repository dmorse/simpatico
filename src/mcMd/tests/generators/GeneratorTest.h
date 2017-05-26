#ifndef MCMD_GENERATOR_TEST_H
#define MCMD_GENERATOR_TEST_H

#include <mcMd/mcSimulation/McSimulation.h>
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/mdSimulation/MdSystem.h>
#include <simp/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <mcMd/potentials/bond/BondPotential.h>
#ifdef SIMP_ANGLE
#include <mcMd/potentials/angle/AnglePotential.h>
#endif
#include <util/containers/DArray.h>

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

#include <string>
#include <fstream>
#include <iostream>

using namespace Util;
using namespace McMd;

class GeneratorTest : public ParamFileTest
{

public:

   GeneratorTest();
   virtual void setUp();

   // Test functions
   void testReadParamBond();

protected:

   McSimulation simulation_;
   McSystem& system_;

   // Utility functions
   void readParam(const char* filename);
   void readConfig(const char* filename);

};

GeneratorTest::GeneratorTest()
 : ParamFileTest(),
   system_(simulation_.system())
{ 
   //setVerbose(2); 
   //ParamComposite::setEcho(true);
}

void GeneratorTest::setUp()
{
   simulation_.fileMaster().setRootPrefix(filePrefix());
} 

void GeneratorTest::readParam(const char* filename)
{  
   openFile(filename); 
   simulation_.readParam(file());
   file().close();
}

void GeneratorTest::readConfig(const char* filename)
{  
   openFile(filename); 
   system_.readConfig(file());
   file().close();
}


#endif
