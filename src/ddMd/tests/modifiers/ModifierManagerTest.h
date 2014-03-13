#ifndef DDMD_MODIFIER_MANAGER_TEST_H
#define DDMD_MODIFIER_MANAGER_TEST_H

#include <ddMd/modifiers/ModifierManager.h>
#include <ddMd/modifiers/Modifier.h>
#include "ModifierClasses.h"
#include "ModifierSubFactory.h"

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

class ModifierManagerTest : public UnitTest 
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
      ModifierManager manager;
      ModifierSubFactory subfactory;
      manager.addSubfactory(subfactory);
   }

   void testReadParam()
   {
      printMethod(TEST_FUNC);
      ModifierManager manager;
      ModifierSubFactory subfactory;
      manager.addSubfactory(subfactory);

      std::ifstream in;
      openInputFile("in/ModifierManager", in);
      //ParamComponent::setEcho(true);
      manager.readParam(in);
      in.close();

      std::cout << std::endl;
      manager.writeParam(std::cout);

      long iStep = 10;
      manager.postIntegrate1(iStep);
   }

};

TEST_BEGIN(ModifierManagerTest)
TEST_ADD(ModifierManagerTest, testConstructor)
TEST_ADD(ModifierManagerTest, testReadParam)
TEST_END(ModifierManagerTest)

#endif 
