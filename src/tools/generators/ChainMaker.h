#ifndef TOOLS_CHAIN_MAKER_H
#define TOOLS_CHAIN_MAKER_H

#include   <inter/bond/HarmonicBond.h>
#include   <util/param/ParamComposite.h>
#include   <util/boundary/Boundary.h>
#include   <util/random/Random.h>
#include   <util/containers/DArray.h>

#include <iostream>

namespace Tools 
{

   using namespace Inter;
   using namespace Util;
   
   class ChainMaker : public ParamComposite
   {
   
   public:
  
      ChainMaker();
      ~ChainMaker();

      void readParam(std::istream& in);
   
      void writeChains(std::ostream& out);
   
      void writeChainsMcMd(std::ostream& out);
   
      void writeChainsDdMd(std::ostream& out);
   
      void writeChainsDdMdMole(std::ostream& out);
   
   private:
   
      Boundary       boundary_;
      Random         random_;
      HarmonicBond   bondPotential_;
      int            nAtomPerMolecule_;
      int            nMolecule_;
      std::string    outputStyle_;
   
   };

}
#endif
