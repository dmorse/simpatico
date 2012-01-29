#ifndef CHAIN_MAKER_H
#define CHAIN_MAKER_H

#include   <util/param/ParamComposite.h>
#include   <mcMd/boundary/Boundary.h>
#include   <mcMd/potentials/BondPotential.h>
#include   <util/random/Random.h>

#include <iostream>

using namespace McMd;
using namespace Util;

class ChainMaker : public ParamComposite
{

public:

   void readParam(std::istream& in);

   void writeChains(std::ostream& out);

private:

   Boundary       boundary_;
   Random         random_;
   BondPotential  bondPotential_;
   int            nAtomPerMolecule_;
   int            nMolecule_;

};
#endif
