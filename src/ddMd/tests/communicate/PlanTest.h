#ifndef DDMD_PLAN_TEST_H
#define DDMD_PLAN_TEST_H

#include <ddMd/communicate/Plan.h>

#ifdef UTIL_MPI
#ifndef TEST_MPI
#define TEST_MPI
#endif
#endif

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <iostream>

using namespace DdMd;

class PlanTest : public UnitTest
{

public:

   virtual void setUp()
   {}

   void testGhostPlan()
   {
      printMethod(TEST_FUNC);

      Plan plan;
      plan.setGhost(0, 0);
      plan.setGhost(1, 0);
      plan.setGhost(2, 1);
      TEST_ASSERT(plan.ghost(0, 0));
      TEST_ASSERT(!plan.ghost(0, 1));
      TEST_ASSERT(plan.ghost(1, 0));
      TEST_ASSERT(!plan.ghost(1, 1));
      TEST_ASSERT(!plan.ghost(2, 0));
      TEST_ASSERT(plan.ghost(2, 1));

      // Clear a field that is not set.
      plan.clearGhost(1, 1);
      TEST_ASSERT(plan.ghost(0, 0));
      TEST_ASSERT(!plan.ghost(0, 1));
      TEST_ASSERT(plan.ghost(1, 0));
      TEST_ASSERT(!plan.ghost(1, 1));
      TEST_ASSERT(!plan.ghost(2, 0));
      TEST_ASSERT(plan.ghost(2, 1));

      // Clear a field that was set.
      plan.clearGhost(2, 1);
      plan.setGhost(2, 0);
      TEST_ASSERT(plan.ghost(0, 0));
      TEST_ASSERT(!plan.ghost(0, 1));
      TEST_ASSERT(plan.ghost(1, 0));
      TEST_ASSERT(!plan.ghost(1, 1));
      TEST_ASSERT(plan.ghost(2, 0));
      TEST_ASSERT(!plan.ghost(2, 1));

   }

   void testExchangePlan()
   {
      printMethod(TEST_FUNC);

      Plan plan;
      plan.setExchange(0, 0);
      plan.setExchange(1, 0);
      plan.setExchange(2, 1);
      TEST_ASSERT(plan.exchange(0, 0));
      TEST_ASSERT(!plan.exchange(0, 1));
      TEST_ASSERT(plan.exchange(1, 0));
      TEST_ASSERT(!plan.exchange(1, 1));
      TEST_ASSERT(!plan.exchange(2, 0));
      TEST_ASSERT(plan.exchange(2, 1));

      // Clear a field that is not set.
      plan.clearExchange(1, 1);
      TEST_ASSERT(plan.exchange(0, 0));
      TEST_ASSERT(!plan.exchange(0, 1));
      TEST_ASSERT(plan.exchange(1, 0));
      TEST_ASSERT(!plan.exchange(1, 1));
      TEST_ASSERT(!plan.exchange(2, 0));
      TEST_ASSERT(plan.exchange(2, 1));

      // Clear a field that was set.
      plan.clearExchange(2, 1);
      plan.setExchange(2, 0);
      TEST_ASSERT(plan.exchange(0, 0));
      TEST_ASSERT(!plan.exchange(0, 1));
      TEST_ASSERT(plan.exchange(1, 0));
      TEST_ASSERT(!plan.exchange(1, 1));
      TEST_ASSERT(plan.exchange(2, 0));
      TEST_ASSERT(!plan.exchange(2, 1));

   }

   void testClear()
   {
      printMethod(TEST_FUNC);

      Plan plan;
      plan.setExchange(0, 0);
      plan.setExchange(1, 0);
      plan.setExchange(2, 1);
      TEST_ASSERT(plan.exchange(0, 0));
      TEST_ASSERT(!plan.exchange(0, 1));
      TEST_ASSERT(plan.exchange(1, 0));
      TEST_ASSERT(!plan.exchange(1, 1));
      TEST_ASSERT(!plan.exchange(2, 0));
      TEST_ASSERT(plan.exchange(2, 1));

      plan.setGhost(0, 0);
      plan.setGhost(1, 0);
      plan.setGhost(2, 1);
      TEST_ASSERT(plan.ghost(0, 0));
      TEST_ASSERT(!plan.ghost(0, 1));
      TEST_ASSERT(plan.ghost(1, 0));
      TEST_ASSERT(!plan.ghost(1, 1));
      TEST_ASSERT(!plan.ghost(2, 0));
      TEST_ASSERT(plan.ghost(2, 1));

      plan.clearFlags();
      TEST_ASSERT(!plan.exchange(0, 0));
      TEST_ASSERT(!plan.exchange(0, 1));
      TEST_ASSERT(!plan.exchange(1, 0));
      TEST_ASSERT(!plan.exchange(1, 1));
      TEST_ASSERT(!plan.exchange(2, 0));
      TEST_ASSERT(!plan.exchange(2, 1));
      TEST_ASSERT(!plan.ghost(0, 0));
      TEST_ASSERT(!plan.ghost(0, 1));
      TEST_ASSERT(!plan.ghost(1, 0));
      TEST_ASSERT(!plan.ghost(1, 1));
      TEST_ASSERT(!plan.ghost(2, 0));
      TEST_ASSERT(!plan.ghost(2, 1));

   }

   void testInserter()
   {
      printMethod(TEST_FUNC);
      Plan plan;
      plan.clearFlags();

      plan.setExchange(0, 0);
      plan.setExchange(1, 0);
      plan.setExchange(2, 1);

      plan.setGhost(0, 0);
      plan.setGhost(1, 1);
      plan.setGhost(2, 0);
      plan.setGhost(2, 1);

      //std::cout << plan;
   }
};

TEST_BEGIN(PlanTest)
TEST_ADD(PlanTest, testGhostPlan)
TEST_ADD(PlanTest, testExchangePlan)
TEST_ADD(PlanTest, testClear)
TEST_ADD(PlanTest, testInserter)
TEST_END(PlanTest)

#endif
