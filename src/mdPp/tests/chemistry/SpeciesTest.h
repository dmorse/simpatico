#ifndef MDPP_SPECIES_TEST_H
#define MDPP_SPECIES_TEST_H

#include <mdPp/chemistry/Species.h>
#include <mdPp/chemistry/Atom.h>
#include <mdPp/chemistry/Molecule.h>

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace MdPp;

class SpeciesTest : public ParamFileTest
{

private:

   Species species_;

public:

   SpeciesTest() 
    : species_()
   {}

   //virtual void setUp(){}

   void testAddAtoms();
   void testRead();

};

inline void SpeciesTest::testAddAtoms()
{
   printMethod(TEST_FUNC);

   species_.initialize(2, 3); 
   species_.setId(3);
   TEST_ASSERT(species_.nAtom() == 2);
   TEST_ASSERT(species_.capacity() == 3);

   DArray<Atom> atoms;
   atoms.allocate(6);
   int i, j, k;
   k = 0;
   for (i=0; i < species_.capacity(); ++i) {
      for (j=0; j < species_.nAtom(); ++j) {
         atoms[k].speciesId  = 3;
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

   Species::MoleculeIterator iter;
   species_.begin(iter);
   i = 0;
   for ( ; iter.notEnd(); ++iter) {
      TEST_ASSERT(iter->id() == i);
      for (j = 0; j < species_.nAtom(); ++j) {
         TEST_ASSERT(iter->atom(j).atomId == j);
         TEST_ASSERT(iter->atom(j).moleculeId == i);
         TEST_ASSERT(iter->atom(j).speciesId == 3);
      }
      ++i;
   }

   species_.clear();
   TEST_ASSERT(species_.size() == 0);
}

inline void SpeciesTest::testRead()
{
   printMethod(TEST_FUNC);

   std::ifstream in;
   openInputFile("in/Species", in);
   in >> species_;
   species_.setId(3);
   TEST_ASSERT(species_.nAtom() == 2);
   TEST_ASSERT(species_.capacity() == 3);
   std::cout << species_;

   DArray<Atom> atoms;
   atoms.allocate(6);
   int i, j, k;
   k = 0;
   for (i=0; i < species_.capacity(); ++i) {
      for (j=0; j < species_.nAtom(); ++j) {
         atoms[k].speciesId  = 3;
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

   Species::MoleculeIterator iter;
   species_.begin(iter);
   i = 0;
   for ( ; iter.notEnd(); ++iter) {
      TEST_ASSERT(iter->id() == i);
      for (j = 0; j < species_.nAtom(); ++j) {
         TEST_ASSERT(iter->atom(j).atomId == j);
         TEST_ASSERT(iter->atom(j).moleculeId == i);
         TEST_ASSERT(iter->atom(j).speciesId == 3);
      }
      ++i;
   }

   species_.clear();
   TEST_ASSERT(species_.size() == 0);
}


TEST_BEGIN(SpeciesTest)
TEST_ADD(SpeciesTest, testAddAtoms)
TEST_ADD(SpeciesTest, testRead)
TEST_END(SpeciesTest)

#endif
