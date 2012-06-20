#ifndef ATOMIC_MAKER_H
#define ATOMIC_MAKER_H

#include   <util/param/ParamComposite.h>
#include   <util/boundary/Boundary.h>
#include   <util/random/Random.h>

#include <iostream>

using namespace Util;

class AtomicMaker : public ParamComposite
{

public:

   void readParam(std::istream& in);

   void writeConfig(std::ostream& out);

private:

   Boundary       boundary_;
   Random         random_;
   int            nMolecule_;

};
#endif
