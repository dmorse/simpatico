#ifndef MDPP_SPECIES_STORAGE_TEST_H
#define MDPP_SPECIES_STORAGE_TEST_H

#include <mdPp/storage/SpeciesStorage.h>
#include <mdPp/chemistry/Atom.h>
#include <mdPp/chemistry/Molecule.h>

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace MdPp;

class SpeciesStorageTest : public ParamFileTest
{

private:

   Species species_;

public:

   SpeciesStorageTest() 
    : species_()
   {}

   //virtual void setUp(){}

   void testAddAtoms();
   void testRead();

};

inline void SpeciesStorageTest::testAddAtoms()
{
   printMethod(TEST_FUNC);

   int nAtom = 2;
   int capacity = 3;
   int speciesId = 4;

   species_.initialize(nAtom, capacity); 
   species_.setId(speciesId);
   TEST_ASSERT(species_.nAtom() == nAtom);
   TEST_ASSERT(species_.capacity() == capacity);

   DArray<Atom> atoms;
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

   Species::Iterator iter;
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

inline void SpeciesStorageTest::testRead()
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

   DArray<Atom> atoms;
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

   Species::Iterator iter;
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


TEST_BEGIN(SpeciesStorageTest)
TEST_ADD(SpeciesStorageTest, testAddAtoms)
TEST_ADD(SpeciesStorageTest, testRead)
TEST_END(SpeciesStorageTest)

#endif
