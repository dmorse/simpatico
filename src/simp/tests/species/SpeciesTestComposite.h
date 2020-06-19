#ifndef SIMP_SPECIES_TEST_COMPOSITE_H
#define SIMP_SPECIES_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "SpeciesGroupTest.h"
#include "SpeciesTest.h"
#include "SpeciesFinderTest.h"
#include "HomopolymerTest.h"
#include "MultiblockTest.h"
#include "HomoRingTest.h"

TEST_COMPOSITE_BEGIN(SpeciesTestComposite)
TEST_COMPOSITE_ADD_UNIT(SpeciesGroupTest);
TEST_COMPOSITE_ADD_UNIT(SpeciesTest);
TEST_COMPOSITE_ADD_UNIT(SpeciesFinderTest);
TEST_COMPOSITE_ADD_UNIT(HomopolymerTest);
TEST_COMPOSITE_ADD_UNIT(MultiblockTest);
TEST_COMPOSITE_ADD_UNIT(HomoRingTest);
TEST_COMPOSITE_END

#endif
