#ifndef DDMD_BOND_STORAGE_TEST_H
#define DDMD_BOND_STORAGE_TEST_H

#include <ddMd/storage/BondStorage.h>
#include <ddMd/storage/GroupIterator.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/GroupStorage.tpp>
#include <ddMd/chemistry/Bond.h>
#include <util/containers/DPArray.h>

#ifdef UTIL_MPI
#ifndef TEST_MPI
#define TEST_MPI
#endif
#endif

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace DdMd;

class BondStorageTest : public ParamFileTest
{
private:

     BondStorage bondStorage_;

public:

   virtual void setUp()
   { 
      #ifdef UTIL_MPI 
      bondStorage_.setIoCommunicator(communicator());
      #endif

      openFile("in/BondStorage"); 
      bondStorage_.readParam(file()); 
   }

   void testReadParam();
   void testAdd();
   void testAddRemove();
   void testIterator();
   void testFindBonds();
   void testClear();

};

inline void BondStorageTest::testReadParam()
{  
   printMethod(TEST_FUNC); 
   if (verbose() > 0) {
      bondStorage_.writeParam(std::cout);
   }
}

inline void BondStorageTest::testAdd()
{
   printMethod(TEST_FUNC);

   Bond* ptr53 = bondStorage_.add(53);
   TEST_ASSERT(bondStorage_.find(53) == ptr53);
   TEST_ASSERT(ptr53->id() == 53);
   TEST_ASSERT(bondStorage_.size() == 1);
   TEST_ASSERT(bondStorage_.isValid());

   Bond* ptr35 = bondStorage_.add(35);
   TEST_ASSERT(bondStorage_.size() == 2);
   TEST_ASSERT(bondStorage_.find(53) == ptr53);
   TEST_ASSERT(ptr53->id() == 53);
   TEST_ASSERT(bondStorage_.find(35) == ptr35);
   TEST_ASSERT(ptr35->id() == 35);
   TEST_ASSERT(bondStorage_.isValid());

   Bond* ptr18 = bondStorage_.add(18);
   TEST_ASSERT(bondStorage_.find(18) == ptr18);
   TEST_ASSERT(ptr18->id() == 18);
   TEST_ASSERT(bondStorage_.size() == 3);
   TEST_ASSERT(bondStorage_.isValid());
}

void BondStorageTest::testAddRemove()
{
   printMethod(TEST_FUNC);

   // Add three bonds
   Bond* ptr53 = bondStorage_.add(53);
   Bond* ptr35 = bondStorage_.add(35);
   Bond* ptr18 = bondStorage_.add(18);

   TEST_ASSERT(bondStorage_.find(53) == ptr53);
   TEST_ASSERT(ptr53->id() == 53);
   TEST_ASSERT(bondStorage_.find(35) == ptr35);
   TEST_ASSERT(ptr35->id() == 35);
   TEST_ASSERT(bondStorage_.find(18) == ptr18);
   TEST_ASSERT(ptr18->id() == 18);
   TEST_ASSERT(bondStorage_.size() == 3);
   TEST_ASSERT(bondStorage_.isValid());

   Bond* newPtr = bondStorage_.newPtr();
   bondStorage_.returnPtr();

   TEST_ASSERT(bondStorage_.find(53) == ptr53);
   TEST_ASSERT(ptr53->id() == 53);
   TEST_ASSERT(bondStorage_.find(35) == ptr35);
   TEST_ASSERT(ptr35->id() == 35);
   TEST_ASSERT(bondStorage_.find(18) == ptr18);
   TEST_ASSERT(ptr18->id() == 18);
   TEST_ASSERT(bondStorage_.size() == 3);
   TEST_ASSERT(bondStorage_.isValid());

   bondStorage_.remove(ptr53);
   TEST_ASSERT(bondStorage_.find(53) == 0);
   TEST_ASSERT(ptr53->id() < 0);
   TEST_ASSERT(bondStorage_.find(35) == ptr35);
   TEST_ASSERT(ptr35->id() == 35);
   TEST_ASSERT(bondStorage_.find(18) == ptr18);
   TEST_ASSERT(ptr18->id() == 18);
   TEST_ASSERT(bondStorage_.size() == 2);
   TEST_ASSERT(bondStorage_.isValid());

   Bond* ptr67 = bondStorage_.add(67);
   Bond* ptr82 = bondStorage_.add(82);
   Bond* ptr44 = bondStorage_.add(44);
   TEST_ASSERT(bondStorage_.find(53) == 0);
   TEST_ASSERT(bondStorage_.find(35) == ptr35);
   TEST_ASSERT(ptr35->id() == 35);
   TEST_ASSERT(bondStorage_.find(18) == ptr18);
   TEST_ASSERT(ptr18->id() == 18);
   TEST_ASSERT(bondStorage_.find(67) == ptr67);
   TEST_ASSERT(ptr67->id() == 67);
   TEST_ASSERT(bondStorage_.find(82) == ptr82);
   TEST_ASSERT(ptr82->id() == 82);
   TEST_ASSERT(bondStorage_.find(44) == ptr44);
   TEST_ASSERT(ptr44->id() == 44);
   TEST_ASSERT(bondStorage_.size() == 5);
   TEST_ASSERT(bondStorage_.isValid());

   bondStorage_.remove(ptr35);
   TEST_ASSERT(bondStorage_.find(53) == 0);
   TEST_ASSERT(bondStorage_.find(35) == 0);
   TEST_ASSERT(ptr35->id() < 0);
   TEST_ASSERT(bondStorage_.find(18) == ptr18);
   TEST_ASSERT(ptr18->id() == 18);
   TEST_ASSERT(bondStorage_.find(67) == ptr67);
   TEST_ASSERT(ptr67->id() == 67);
   TEST_ASSERT(bondStorage_.find(82) == ptr82);
   TEST_ASSERT(ptr82->id() == 82);
   TEST_ASSERT(bondStorage_.find(44) == ptr44);
   TEST_ASSERT(ptr44->id() == 44);
   TEST_ASSERT(bondStorage_.size() == 4);
   TEST_ASSERT(bondStorage_.isValid());
   #if 0
   #endif

}
 
void BondStorageTest::testIterator()
{
   printMethod(TEST_FUNC);

   DPArray<Bond> bonds;

   bonds.allocate(bondStorage_.capacity());

   // Add bonds
   bonds.append(*bondStorage_.add(53));
   bonds.append(*bondStorage_.add(35));
   bonds.append(*bondStorage_.add(18));
   bonds.append(*bondStorage_.add(44));
   bonds.append(*bondStorage_.add(17));
   bonds.append(*bondStorage_.add(82));
   bonds.append(*bondStorage_.add(39));
   TEST_ASSERT(bondStorage_.size() == 7);
   TEST_ASSERT(bondStorage_.isValid());
 
   GroupIterator<2> localIter;
   int n = 0; 
   for (bondStorage_.begin(localIter); localIter.notEnd(); ++localIter) {
      ++n;
      //std::cout << localIter->id() << std::endl;
   }
   TEST_ASSERT(n == bondStorage_.size());
   TEST_ASSERT(n == 7);

   bondStorage_.remove(&bonds[1]);
   --n;
   TEST_ASSERT(bondStorage_.isValid());
   TEST_ASSERT(n == bondStorage_.size());
   TEST_ASSERT(n == 6);

}

inline void BondStorageTest::testFindBonds()
{
   printMethod(TEST_FUNC);

   DPArray<Bond> bonds;

   bonds.allocate(bondStorage_.capacity());

   // Add bonds
   Bond* ptr53 = bondStorage_.add(53); // 0
   Bond* ptr35 = bondStorage_.add(35); // 1
   Bond* ptr18 = bondStorage_.add(18); // 2
   Bond* ptr44 = bondStorage_.add(44); // 3
   Bond* ptr17 = bondStorage_.add(17); // 4
   Bond* ptr82 = bondStorage_.add(82); // 5
   Bond* ptr39 = bondStorage_.add(39); // 6

   bonds.append(*ptr53); // 0
   bonds.append(*ptr35); // 1
   bonds.append(*ptr18); // 2
   bonds.append(*ptr44); // 3
   bonds.append(*ptr17); // 4
   bonds.append(*ptr82); // 5
   bonds.append(*ptr39); // 6
  
   TEST_ASSERT(bondStorage_.find(53) == ptr53);
   TEST_ASSERT(bondStorage_.find(82) == ptr82);
   TEST_ASSERT(ptr53->id() == 53);
   TEST_ASSERT(bondStorage_.find(54) == 0);
   bondStorage_.remove(ptr82);
   TEST_ASSERT(bondStorage_.find(53) == ptr53);
   TEST_ASSERT(bondStorage_.find(18) == ptr18);
   TEST_ASSERT(bondStorage_.find(44) == ptr44);
   TEST_ASSERT(bondStorage_.find(17) == ptr17);
   TEST_ASSERT(bondStorage_.find(82) == 0);
   TEST_ASSERT(bondStorage_.find(83) == 0);

}

inline void BondStorageTest::testClear()
{
   printMethod(TEST_FUNC);

   // Add bonds
   Bond* ptr53 = bondStorage_.add(53); // 0
   Bond* ptr35 = bondStorage_.add(35); // 1
   Bond* ptr18 = bondStorage_.add(18); // 2
   Bond* ptr44 = bondStorage_.add(44); // 3
   Bond* ptr17 = bondStorage_.add(17); // 4
   Bond* ptr82 = bondStorage_.add(82); // 5
   Bond* ptr39 = bondStorage_.add(39); // 6
  
   TEST_ASSERT(bondStorage_.size() == 7);
   bondStorage_.remove(ptr82);
   TEST_ASSERT(bondStorage_.size() == 6);
   bondStorage_.clearGroups();
   TEST_ASSERT(bondStorage_.isValid());
   TEST_ASSERT(bondStorage_.size() == 0);

}

TEST_BEGIN(BondStorageTest)
TEST_ADD(BondStorageTest, testReadParam)
TEST_ADD(BondStorageTest, testAdd)
TEST_ADD(BondStorageTest, testAddRemove)
TEST_ADD(BondStorageTest, testIterator)
TEST_ADD(BondStorageTest, testFindBonds)
TEST_ADD(BondStorageTest, testClear)
TEST_END(BondStorageTest)

#endif
