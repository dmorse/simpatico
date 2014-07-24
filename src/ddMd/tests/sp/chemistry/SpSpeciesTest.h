#ifndef DDMD_SP_SPECIES_TEST_H
#define DDMD_SP_SPECIES_TEST_H

#include <ddMd/sp/chemistry/SpSpecies.h>
#include <ddMd/sp/chemistry/SpAtom.h>
#include <ddMd/sp/chemistry/SpMolecule.h>

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace DdMd;

class SpSpeciesTest : public ParamFileTest
{

private:

   SpSpecies species_;

public:

   SpSpeciesTest() 
    : species_()
   {}

   //virtual void setUp(){}

   void testAddAtoms();
   void testRead();

};

inline void SpSpeciesTest::testAddAtoms()
{
   printMethod(TEST_FUNC);

   int nAtom = 2;
   int capacity = 3;
   int speciesId = 4;

   species_.initialize(nAtom, capacity); 
   species_.setId(speciesId);
   TEST_ASSERT(species_.nAtom() == nAtom);
   TEST_ASSERT(species_.capacity() == capacity);

   DArray<SpAtom> atoms;
   atoms.allocate(nAtom*capacity);
   int i, j, k;
   k = 0;
   for (i=0; i < species_.capacity(); ++i) {
      for (j=0; j < species_.nAtom(); ++j) {
         atoms[k].speciesId  = speciesId;
         atoms[k].moleculeId = i;
         atoms[k].atomId = j;
         ++k;
      }
   }
   species_.addAtom(atoms[0]);
   species_.addAtom(atoms[3]);
   species_.addAtom(atoms[2]);
   species_.addAtom(atoms[1]);
   TEST_ASSERT(species_.size() == 2);
   TEST_ASSERT(species_.isValid());

   SpSpecies::SpMoleculeIterator iter;
   i = 0;
   for (species_.begin(iter) ; iter.notEnd(); ++iter) {
      TEST_ASSERT(iter->id() == i);
      for (j = 0; j < species_.nAtom(); ++j) {
         TEST_ASSERT(iter->atom(j).atomId == j);
         TEST_ASSERT(iter->atom(j).moleculeId == i);
         TEST_ASSERT(iter->atom(j).speciesId == speciesId);
      }
      ++i;
   }
   TEST_ASSERT(i == species_.size());

   species_.clear();
   TEST_ASSERT(species_.size() == 0);
}

inline void SpSpeciesTest::testRead()
{
   printMethod(TEST_FUNC);

   int nAtom = 2;
   int capacity = 3;
   int speciesId = 4;

   std::ifstream in;
   openInputFile("in/Species", in);
   in >> species_;
   species_.setId(speciesId);
   TEST_ASSERT(species_.nAtom() == nAtom);
   TEST_ASSERT(species_.capacity() == capacity);
   std::cout << species_;

   DArray<SpAtom> atoms;
   atoms.allocate(nAtom*capacity);
   int i, j, k;
   k = 0;
   for (i=0; i < species_.capacity(); ++i) {
      for (j=0; j < species_.nAtom(); ++j) {
         atoms[k].speciesId  = speciesId;
         atoms[k].moleculeId = i;
         atoms[k].atomId = j;
         ++k;
      }
   }

   species_.addAtom(atoms[0]);
   species_.addAtom(atoms[3]);
   species_.addAtom(atoms[2]);
   species_.addAtom(atoms[1]);
   TEST_ASSERT(species_.size() == 2);
   TEST_ASSERT(species_.id() == speciesId);
   TEST_ASSERT(species_.isValid());

   SpSpecies::SpMoleculeIterator iter;
   i = 0;
   for (species_.begin(iter); iter.notEnd(); ++iter) {
      TEST_ASSERT(iter->id() == i);
      for (j = 0; j < species_.nAtom(); ++j) {
         TEST_ASSERT(iter->atom(j).atomId == j);
         TEST_ASSERT(iter->atom(j).moleculeId == i);
         TEST_ASSERT(iter->atom(j).speciesId == speciesId);
      }
      ++i;
   }
   TEST_ASSERT(i == species_.size());

   species_.clear();
   TEST_ASSERT(species_.size() == 0);
   TEST_ASSERT(species_.id() == speciesId);
}


TEST_BEGIN(SpSpeciesTest)
TEST_ADD(SpSpeciesTest, testAddAtoms)
TEST_ADD(SpSpeciesTest, testRead)
TEST_END(SpSpeciesTest)

#endif
