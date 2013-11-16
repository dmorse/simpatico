#ifndef DDMD_STORAGE_TEST_COMPOSITE_H
#define DDMD_STORAGE_TEST_COMPOSITE_H

#include "AtomMapTest.h"
#include "AtomStorageTest.h"
#include "BondStorageTest.h"

#include <test/CompositeTestRunner.h>

TEST_COMPOSITE_BEGIN(StorageTestComposite)
TEST_COMPOSITE_ADD_UNIT(AtomMapTest);
TEST_COMPOSITE_ADD_UNIT(AtomStorageTest);
TEST_COMPOSITE_ADD_UNIT(BondStorageTest);
TEST_COMPOSITE_END

#endif
