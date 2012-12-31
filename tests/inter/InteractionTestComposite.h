#ifndef INTERACTION_TEST_COMPOSITE_H
#define INTERACTION_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "pair/PairTestComposite.h"
#include "bond/HarmonicBondTest.h"
#include "bond/HarmonicL0BondTest.h"
#include "angle/AngleTestComposite.h"
#ifdef MCMC_DIHEDRAL
#include "dihedral/CosineDihedralTest.h"
#endif

TEST_COMPOSITE_BEGIN(InteractionTestComposite)
addChild(new PairTestComposite, "pair/");

TEST_COMPOSITE_ADD_UNIT(HarmonicBondTest);
TEST_COMPOSITE_ADD_UNIT(HarmonicL0BondTest);

//TEST_COMPOSITE_ADD_UNIT(AngleTest);
addChild(new AngleTestComposite, "angle/");

#ifdef INTER_DIHEDRAL
TEST_COMPOSITE_ADD_UNIT(CosineDihedralTest);
#endif

TEST_COMPOSITE_END

#endif
