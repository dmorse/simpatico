#ifndef TOOLS_CHAIN_MAKER_H
#define TOOLS_CHAIN_MAKER_H

#include   <simp/interaction/bond/HarmonicBond.h>
#include   <util/param/ParamComposite.h>
#include   <util/boundary/Boundary.h>
#include   <util/random/Random.h>
#include   <util/containers/DArray.h>

#include <iostream>

namespace Tools 
{

   using namespace Simp;
   using namespace Util;

   /**
   * Generates a simulation box full of randomly generated linear polymers.
   */   
   class ChainMaker : public ParamComposite
   {
   
   public:
 
      /*
      * Constructor
      */ 
      ChainMaker();

      /*
      * Destructor.
      */ 
      ~ChainMaker();

      /*
      * Read parameters from input file.
      */ 
      void readParam(std::istream& in);
   
      /*
      * Write configuration.
      */ 
      void writeChains(std::ostream& out);
   
      /*
      * Write configuration in default McMd format.
      */ 
      void writeChainsMcMd(std::ostream& out);
   
      /*
      * Write configuration in default DdMd format.
      */ 
      void writeChainsDdMd(std::ostream& out);
   
      /*
      * Write configuration in default DdMd molecular format.
      */ 
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
