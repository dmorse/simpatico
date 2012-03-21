#ifndef MCMD_TETHER_MASTER_TEST_H
#define MCMD_TETHER_MASTER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <mcMd/tethers/TetherMaster.h>
#include <mcMd/chemistry/Atom.h>
#include <util/containers/RArray.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace McMd;

class TetherMasterTest : public UnitTest 
{

public:

   void setUp()
   {
      atomCapacity_ = 100;
      Atom::allocate(atomCapacity_, atoms_);

      std::ifstream file("in/TetherMaster");
      tetherMaster_.readParam(file);
      file.close();

   }

   void tearDown()
   {
      Atom::deallocate();
   }

   void testReadParam() 
   {
      printMethod(TEST_FUNC);

      std::cout << std::endl;
      tetherMaster_.writeParam(std::cout);
      tetherMaster_.isValid();
   }

   void testAdd() 
   {
      printMethod(TEST_FUNC);
      Vector        anchor;
      const Tether* tetherPtr;

      anchor[0] = 0.5;
      anchor[1] = 1.0;
      anchor[2] = 1.5;

      tetherMaster_.addTether(atoms_[5], anchor);
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[5]));
      tetherPtr = &tetherMaster_.tether(atoms_[5]);
      TEST_ASSERT(tetherPtr->hasAtom());
      TEST_ASSERT(&tetherPtr->atom() == &atoms_[5]);
      TEST_ASSERT(eq(tetherPtr->anchor()[0], anchor[0]));
      TEST_ASSERT(eq(tetherPtr->anchor()[1], anchor[1]));
      TEST_ASSERT(eq(tetherPtr->anchor()[2], anchor[2]));
      TEST_ASSERT(tetherMaster_.nTether() == 1);
      tetherMaster_.isValid();

      tetherMaster_.addTether(atoms_[9], anchor);
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[9]));
      tetherPtr = &tetherMaster_.tether(atoms_[9]);
      TEST_ASSERT(tetherPtr->hasAtom());
      TEST_ASSERT(&tetherPtr->atom() == &atoms_[9]);
      TEST_ASSERT(tetherMaster_.nTether() == 2);
      tetherMaster_.isValid();

   }

   void testAddRemove() 
   {
      printMethod(TEST_FUNC);
      Vector          anchor;
      //const Tether* tetherPtr;

      anchor[0] = 0.5;
      anchor[1] = 1.0;
      anchor[2] = 1.5;

      tetherMaster_.addTether(atoms_[2], anchor);
      tetherMaster_.addTether(atoms_[5], anchor);
      tetherMaster_.addTether(atoms_[9], anchor);
      tetherMaster_.addTether(atoms_[13], anchor);
      tetherMaster_.addTether(atoms_[39], anchor);
      tetherMaster_.addTether(atoms_[68], anchor);
      tetherMaster_.isValid();

      TEST_ASSERT(tetherMaster_.isTethered(atoms_[2]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[5]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[9]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[13]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[39]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[68]));
      TEST_ASSERT(tetherMaster_.nTether() == 6);

      tetherMaster_.removeTether(atoms_[13]);
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[2]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[5]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[9]));
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[13]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[39]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[68]));
      TEST_ASSERT(tetherMaster_.nTether() == 5);
      tetherMaster_.isValid();

      tetherMaster_.removeTether(atoms_[5]);
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[2]));
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[5]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[9]));
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[13]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[39]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[68]));
      TEST_ASSERT(tetherMaster_.nTether() == 4);
      tetherMaster_.isValid();

      tetherMaster_.addTether(atoms_[84], anchor);
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[2]));
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[5]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[9]));
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[13]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[39]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[68]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[84]));
      TEST_ASSERT(tetherMaster_.nTether() == 5);
      tetherMaster_.isValid();

   }

   void testAddRemovePair() 
   {
      printMethod(TEST_FUNC);
      Vector          anchor;

      anchor[0] = 0.5;
      anchor[1] = 1.0;
      anchor[2] = 1.5;

      tetherMaster_.addTether(atoms_[2], anchor);
      tetherMaster_.addTether(atoms_[5], anchor);
      tetherMaster_.addTether(atoms_[9], anchor);
      tetherMaster_.addTether(atoms_[13], anchor);
      tetherMaster_.addTether(atoms_[39], anchor);
      tetherMaster_.addTether(atoms_[68], anchor);
      TEST_ASSERT(tetherMaster_.nTether() == 6);
      tetherMaster_.isValid();

      tetherMaster_.removeTether(atoms_[13]);
      tetherMaster_.removeTether(atoms_[5]);
      tetherMaster_.addTether(atoms_[84], anchor);
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[2]));
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[5]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[9]));
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[13]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[39]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[68]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[84]));
      TEST_ASSERT(tetherMaster_.nTether() == 5);
      tetherMaster_.isValid();

      tetherMaster_.pairTethers(atoms_[9], atoms_[39]);
      tetherMaster_.isValid();

      tetherMaster_.removeTether(atoms_[9]);
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[2]));
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[5]));
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[9]));
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[13]));
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[39]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[68]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[84]));
      TEST_ASSERT(tetherMaster_.nTether() == 3);
      tetherMaster_.isValid();

      tetherMaster_.addTether(atoms_[13], anchor);
      Tether* tetherPtr1 = &tetherMaster_.tether(atoms_[2]);
      Tether* tetherPtr2 = &tetherMaster_.tether(atoms_[68]);
      tetherMaster_.pairTethers(*tetherPtr1, *tetherPtr2);
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[2]));
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[5]));
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[9]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[13]));
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[39]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[68]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[84]));
      TEST_ASSERT(tetherMaster_.nTether() == 4);
      tetherMaster_.isValid();

      tetherMaster_.removeTether(atoms_[2]);
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[2]));
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[5]));
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[9]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[13]));
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[39]));
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[68]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[84]));
      TEST_ASSERT(tetherMaster_.nTether() == 2);

   }

   void testAddRemoveSwitch1() 
   {
      printMethod(TEST_FUNC);
      Vector anchor;

      anchor[0] = 0.5;
      anchor[1] = 1.0;
      anchor[2] = 1.5;

      tetherMaster_.addTether(atoms_[2], anchor);
      tetherMaster_.addTether(atoms_[5], anchor);
      tetherMaster_.addTether(atoms_[9], anchor);
      tetherMaster_.addTether(atoms_[13], anchor);
      tetherMaster_.addTether(atoms_[39], anchor);
      tetherMaster_.addTether(atoms_[68], anchor);
      TEST_ASSERT(tetherMaster_.nTether() == 6);
      tetherMaster_.isValid();

      tetherMaster_.removeTether(atoms_[13]);
      tetherMaster_.removeTether(atoms_[5]);
      tetherMaster_.addTether(atoms_[84], anchor);
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[2]));
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[5]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[9]));
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[13]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[39]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[68]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[84]));
      TEST_ASSERT(tetherMaster_.nTether() == 5);
      tetherMaster_.isValid();

      tetherMaster_.transferTether(atoms_[39], atoms_[43]);
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[2]));
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[5]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[9]));
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[13]));
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[39]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[43]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[68]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[84]));
      TEST_ASSERT(tetherMaster_.nTether() == 5);
      tetherMaster_.isValid();

   }

   void testAddRemoveSwitch2() 
   {
      printMethod(TEST_FUNC);
      Vector          anchor;

      anchor[0] = 0.5;
      anchor[1] = 1.0;
      anchor[2] = 1.5;

      tetherMaster_.addTether(atoms_[2], anchor);
      tetherMaster_.addTether(atoms_[5], anchor);
      tetherMaster_.addTether(atoms_[9], anchor);
      tetherMaster_.addTether(atoms_[13], anchor);
      tetherMaster_.addTether(atoms_[39], anchor);
      tetherMaster_.addTether(atoms_[68], anchor);
      TEST_ASSERT(tetherMaster_.nTether() == 6);
      tetherMaster_.isValid();

      tetherMaster_.removeTether(atoms_[13]);
      tetherMaster_.removeTether(atoms_[5]);
      tetherMaster_.addTether(atoms_[84], anchor);
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[2]));
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[5]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[9]));
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[13]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[39]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[68]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[84]));
      TEST_ASSERT(tetherMaster_.nTether() == 5);
      tetherMaster_.isValid();

      Tether& tether = tetherMaster_.tether(atoms_[39]);
      tetherMaster_.transferTether(tether, atoms_[43]);
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[2]));
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[5]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[9]));
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[13]));
      TEST_ASSERT(!tetherMaster_.isTethered(atoms_[39]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[43]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[68]));
      TEST_ASSERT(tetherMaster_.isTethered(atoms_[84]));
      TEST_ASSERT(tetherMaster_.nTether() == 5);
      tetherMaster_.isValid();

   }

private:

   TetherMaster tetherMaster_;

   RArray<Atom> atoms_;

   int atomCapacity_;

};

TEST_BEGIN(TetherMasterTest)
TEST_ADD(TetherMasterTest, testReadParam)
TEST_ADD(TetherMasterTest, testAdd)
TEST_ADD(TetherMasterTest, testAddRemove)
TEST_ADD(TetherMasterTest, testAddRemovePair)
TEST_ADD(TetherMasterTest, testAddRemoveSwitch1)
TEST_ADD(TetherMasterTest, testAddRemoveSwitch2)
TEST_END(TetherMasterTest)

#endif
