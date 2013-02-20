#ifndef INTERACTION_TEST_COMPOSITE_H
#define INTERACTION_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "pair/PairTestComposite.h"
#include "bond/BondTestComposite.h"
#ifdef INTER_ANGLE
#include "angle/AngleTestComposite.h"
#endif
#ifdef INTER_DIHEDRAL
#include "dihedral/DihedralTestComposite.h"
#endif

TEST_COMPOSITE_BEGIN(InteractionTestComposite)
addChild(new PairTestComposite, "pair/");
addChild(new BondTestComposite, "bond/");

#ifdef INTER_ANGLE
addChild(new AngleTestComposite, "angle/");
#endif

#ifdef INTER_DIHEDRAL
addChild(new DihedralTestComposite, "dihedral/");
#endif

TEST_COMPOSITE_END

#endif
