#ifndef DDMD_GROUP_TEST_H
#define DDMD_GROUP_TEST_H

#include <ddMd/chemistry/Group.h>
#include <ddMd/chemistry/AtomArray.h>
#include <ddMd/chemistry/Atom.h>
#include <util/space/Vector.h>
#include <util/archives/BinaryFileOArchive.h>
#include <util/archives/BinaryFileIArchive.h>

#ifdef UTIL_MPI
#ifndef TEST_MPI
#define TEST_MPI
#endif
#endif

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
   void testSerialize();
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

   if (isIoProcessor()) {
      std::cout << std::endl;
      std::cout << group;
      std::cout << std::endl;
   }

}

void GroupTest::testSerialize() 
{
   printMethod(TEST_FUNC);
   Group<2> v;
   int i1 = 35;
   int i2 = 43;

   // Read from input file
   std::ifstream in;
   openInputFile("in/Group", in);
   in >> v;

   // Write to binary file archive
   BinaryFileOArchive oa;
   openOutputFile("binary", oa.file());
   oa << i1;
   oa << v;
   oa << i2;
   oa.file().close();

   // Write to binary file archive
   Group<2> u;
   int j1, j2;
   BinaryFileIArchive ia;
   openInputFile("binary", ia.file());
   ia >> j1;
   ia >> u;
   ia >> j2;
   
   TEST_ASSERT(j1 == i1);
   TEST_ASSERT(u.typeId() == v.typeId());
   for (int i = 0; i < 2; ++i) {
      TEST_ASSERT(u.atomId(i) == v.atomId(i));
   }
   TEST_ASSERT(j2 == i2);

   if (isIoProcessor()) {
      std::cout << std::endl;
      std::cout << u;
      std::cout << std::endl;
   }

}

TEST_BEGIN(GroupTest)
TEST_ADD(GroupTest, testConstructor)
TEST_ADD(GroupTest, testSetGet)
TEST_ADD(GroupTest, testFileIo)
TEST_ADD(GroupTest, testSerialize)
TEST_END(GroupTest)

#endif
