#ifndef SIMP_SPECIES_FINDER_TEST_H
#define SIMP_SPECIES_FINDER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <simp/species/SpeciesFinder.h>
#include <util/containers/DArray.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace McMd;

class SpeciesFinderTest : public UnitTest 
{

private:

   SpeciesFinder finder;
   int           nSpecies;
   DArray<int>   nMolecules;
   DArray<int>   nParts;

public:

   void setUp() 
   { 
      setVerbose(1);
      nSpecies = 5;
      nMolecules.allocate(nSpecies);
      nParts.allocate(nSpecies);
      nMolecules[0] = 12;
      nParts[0]     = 8;
      nMolecules[1] = 0;
      nParts[1]     = 7;
      nMolecules[2] = 9;
      nParts[2]     = 3;
      nMolecules[3] = 6;
      nParts[3]     = 0;
      nMolecules[4] = 5;
      nParts[4]     = 8;

      finder.allocate(nSpecies);
      for (int i = 0; i < nSpecies; ++i) {
         finder.setSpecies(i, nMolecules[i], nParts[i]);
      }
      finder.initialize();
   } 

   void tearDown() 
   {}
  
   void testInitialize()
   {
      printMethod(TEST_FUNC);
   }

   void testFindMolecule()
   {
      printMethod(TEST_FUNC);

      SpeciesFinder::Molecule context;
      int is, im, j;
      j = 0;
      for (is = 0; is < nSpecies; ++is) {
         for (im = 0; im < nMolecules[is]; ++im) {
            finder.findMolecule(j, context);
            TEST_ASSERT(context.speciesId == is);
            TEST_ASSERT(context.moleculeId == im);
            ++j;
         }
      }
   }

   void testFindPart()
   {
      printMethod(TEST_FUNC);

      SpeciesFinder::Part context;
      int is, im, ip, j;
      j = 0;
      for (is = 0; is < nSpecies; ++is) {
         // std::cout << "Species " << is 
         //          << " nMolecule " << nMolecules[is]
         //          << " nPart " << nParts[is] << "\n";
         for (im = 0; im < nMolecules[is]; ++im) {
            for (ip = 0; ip < nParts[is]; ++ip) {
               finder.findPart(j, context);
               TEST_ASSERT(context.speciesId == is);
               TEST_ASSERT(context.moleculeId == im);
               TEST_ASSERT(context.partId == ip);
               // std::cout << j << "  "
               //          << context.speciesId << "  "
               //          << context.moleculeId << "  "
               //          << context.partId << "\n";
               ++j;
            }
         }
      }
   }

};

TEST_BEGIN(SpeciesFinderTest)
TEST_ADD(SpeciesFinderTest, testInitialize)
TEST_ADD(SpeciesFinderTest, testFindMolecule)
TEST_ADD(SpeciesFinderTest, testFindPart)
TEST_END(SpeciesFinderTest)

#endif
