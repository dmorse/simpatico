#ifndef SPAN_TYPE_MAP_TEST_H
#define SPAN_TYPE_MAP_TEST_H

#include <spAn/config/TypeMap.h>

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace SpAn;

class TypeMapTest : public ParamFileTest
{

private:

   TypeMap map_;

public:

   TypeMapTest() 
    : map_()
   {}

   virtual void setUp()
   { 
      setVerbose(2);
   }

   void testRead();

};

inline void TypeMapTest::testRead()
{  
   printMethod(TEST_FUNC); 
   std::ifstream file;
   openInputFile("in/TypeMap", file);
   map_.read(file); 
   file.close(); 

   TEST_ASSERT(map_.size() == 3);
   TEST_ASSERT(2 == map_.id("A"));
   TEST_ASSERT(1 == map_.id("B"));
   TEST_ASSERT(0 == map_.id("C"));
   TEST_ASSERT("C" == map_.name(0));
   TEST_ASSERT("B" == map_.name(1));
   TEST_ASSERT("A" == map_.name(2));
}

TEST_BEGIN(TypeMapTest)
TEST_ADD(TypeMapTest, testRead)
TEST_END(TypeMapTest)

#endif
