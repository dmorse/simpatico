#ifndef DDMD_GROUP_TEST_H
#define DDMD_GROUP_TEST_H

#include <ddMd/chemistry/Group.h>
#include <ddMd/chemistry/AtomArray.h>
#include <ddMd/chemistry/Atom.h>
#include <util/space/Vector.h>

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace DdMd;

class GroupTest : public UnitTest 
{

public:

   void setUp() {}

   void tearDown() {}
 
   void testConstructor();

   void testSetGet();

   void testFileIo();

};


void GroupTest::testConstructor()
{
   printMethod(TEST_FUNC);
   Group<2> group;
} 

void GroupTest::testSetGet()
{
   printMethod(TEST_FUNC);
   Group<2> group;
   AtomArray atoms;
   atoms.allocate(5);

   TEST_ASSERT(group.nPtr() == 0);

   group.setId(3);
   group.setTypeId(4);
   group.setAtomId(0, 34);
   group.setAtomId(1, 35);
   group.setAtomPtr(1, &atoms[3]);

   TEST_ASSERT(group.id() == 3);
   TEST_ASSERT(group.typeId() == 4);
   TEST_ASSERT(group.atomId(0) == 34);
   TEST_ASSERT(group.atomId(1) == 35);
   TEST_ASSERT(group.atomPtr(0) == 0);
   TEST_ASSERT(group.atomPtr(1) == &atoms[3]);
   TEST_ASSERT(group.nPtr() == 1);

   group.clearAtomPtr(1);
   TEST_ASSERT(group.nPtr() == 0);
   TEST_ASSERT(group.atomPtr(0) == 0);
   TEST_ASSERT(group.atomPtr(1) == 0);

   group.setAtomPtr(0, &atoms[3]);
   group.setAtomPtr(1, &atoms[4]);
   TEST_ASSERT(group.nPtr() == 2);

   group.clearAtomPtr(0);
   TEST_ASSERT(group.atomPtr(0) == 0);
   TEST_ASSERT(group.atomPtr(1) == &atoms[4]);
   TEST_ASSERT(group.nPtr() == 1);
} 

void GroupTest::testFileIo()
{
   printMethod(TEST_FUNC);
   Group<2> group;

   std::ifstream file;
   file.open("in/Group");

   file >> group;

   std::cout << std::endl;
   std::cout << group;
   std::cout << std::endl;

}

TEST_BEGIN(GroupTest)
TEST_ADD(GroupTest, testConstructor)
TEST_ADD(GroupTest, testSetGet)
TEST_ADD(GroupTest, testFileIo)
TEST_END(GroupTest)

#endif
