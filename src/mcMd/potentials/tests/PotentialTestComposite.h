#ifndef MCMD_POTENTIAL_TEST_COMPOSITE_H
#define MCMD_POTENTIAL_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "pair/LJPairTest.h"
#include "bond/HarmonicBondTest.h"
#include "bond/HarmonicL0BondTest.h"
#ifdef MCMC_DIHEDRAL
#include "dihedral/CosineDihedralTest.h"
#endif

TEST_COMPOSITE_BEGIN(PotentialTestComposite)
TEST_COMPOSITE_ADD_UNIT(LJPairTest);
TEST_COMPOSITE_ADD_UNIT(HarmonicBondTest);
TEST_COMPOSITE_ADD_UNIT(HarmonicL0BondTest);
#ifdef MCMC_DIHEDRAL
TEST_COMPOSITE_ADD_UNIT(CosineDihedralTest);
#endif
TEST_COMPOSITE_END

#endif
