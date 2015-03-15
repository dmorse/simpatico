#ifndef COMMUNICATE_TEST_COMPOSITE_H
#define COMMUNICATE_TEST_COMPOSITE_H

#include "DomainTest.h"
#include "BufferTest.h"
#include "AtomDistributorTest.h"
#include "GroupDistributorTest.h"
#include "PlanTest.h"
#include "ExchangerTest.h"
#include "ExchangerForceTest.h"
#include "AtomCollectorTest.h"
#include "BondCollectorTest.h"

#include <test/CompositeTestRunner.h>

TEST_COMPOSITE_BEGIN(CommunicateTestComposite)

TEST_COMPOSITE_ADD_UNIT(DomainTest)
TEST_COMPOSITE_ADD_UNIT(BufferTest)
TEST_COMPOSITE_ADD_UNIT(AtomDistributorTest)
TEST_COMPOSITE_ADD_UNIT(GroupDistributorTest)
TEST_COMPOSITE_ADD_UNIT(PlanTest)
TEST_COMPOSITE_ADD_UNIT(ExchangerTest)
TEST_COMPOSITE_ADD_UNIT(ExchangerForceTest)
TEST_COMPOSITE_ADD_UNIT(AtomCollectorTest)
TEST_COMPOSITE_ADD_UNIT(BondCollectorTest)
TEST_COMPOSITE_END

#endif
