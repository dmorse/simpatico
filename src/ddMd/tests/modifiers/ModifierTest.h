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
   {
      Label::clear();
      ParamComponent::setEcho(false);
   }

   void tearDown()
   {}

   void testConstructor()
   {
      printMethod(TEST_FUNC);
      DdMd::ModifierA modifier;
      TEST_ASSERT(modifier.isSet(Modifier::Flags::PostIntegrate1));
      TEST_ASSERT(!modifier.isSet(Modifier::Flags::PreIntegrate1));

      //std::cout << std::endl;
      //std::cout << modifier.flags() << std::endl;
   }

   void testReadParam()
   {
      printMethod(TEST_FUNC);
      ModifierA modifier;
      std::ifstream in;
      openInputFile("in/ModifierA", in);
      modifier.readParam(in);

      TEST_ASSERT(modifier.interval() == 10);
      TEST_ASSERT(modifier.isAtInterval(10));
      TEST_ASSERT(modifier.isAtInterval(20));
      TEST_ASSERT(!modifier.isAtInterval(5));

      std::cout << std::endl;
      modifier.writeParam(std::cout);
   }

};

TEST_BEGIN(ModifierTest)
TEST_ADD(ModifierTest, testConstructor)
TEST_ADD(ModifierTest, testReadParam)
TEST_END(ModifierTest)

#endif //ifndef DDMD_MODIFIER_TEST_H
