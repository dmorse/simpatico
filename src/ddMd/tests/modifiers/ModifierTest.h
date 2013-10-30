#ifndef DDMD_MODIFIER_TEST_H
#define DDMD_MODIFIER_TEST_H

#include <ddMd/modifiers/Modifier.h>
#include "ModifierClasses.h"

#ifdef UTIL_MPI
#ifndef TEST_MPI
#define TEST_MPI
#endif
#endif

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <iostream>

using namespace Util;
using namespace DdMd;

class ModifierTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   void testConstructor()
   {
      printMethod(TEST_FUNC);
      ModifierA modifier;
      TEST_ASSERT(modifier.isSet(Modifier::Flags::PostIntegrate1));
      TEST_ASSERT(!modifier.isSet(Modifier::Flags::PreIntegrate1));
      cout << modifier.mask() << std::endl;
   }

};

TEST_BEGIN(ModifierTest)
TEST_ADD(ModifierTest, testConstructor)
TEST_END(ModifierTest)

#endif //ifndef DDMD_MODIFIER_TEST_H
