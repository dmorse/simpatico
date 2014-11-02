#ifndef UTIL_AUTOCORRELATION_TPP
#define UTIL_AUTOCORRELATION_TPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AutoCorrelation.h"  

namespace Util
{

   /*
   * Constructor
   *
   * \param blockFactor ratio of block sizes of subsequent stages
   */
   AutoCorrelation(int blockFactor = 4)
    : AutoCorrStage(blockFactor)
   {}

   /*
   * Set the base output file name.
   */
   void setOutputFileName(std::string outputFileName)
   {  outputFileName_ = outputFileName; }

   /*
   * Register the creation of a descendant stage.
   */
   virtual void registerDescendant(AutoCorrelation* ptr)
   {}

   /*
   * Load state from an archive.
   */
   virtual void loadParameters(Serializable::IArchive& ar)
   {
      AutoCorrStage::load(ar);
   }

   /*
   * Save state to an archive.
   */
   virtual void save(Serializable::OArchive& ar)
   {
      AutoCorrStage::save(ar);
   }

   /*
   * Serialize to/from an archive.
   */
   template <class Archive>
   void serialize(Archive& ar, const unsigned int version)
   {
      AutoCorrStage::serialize(ar, version);
   }

   /*
   * Output the autocorrelation function
   */
   virtual void output()
   {}

}
#endif
