#ifndef UTIL_AUTOCORRELATION_TPP
#define UTIL_AUTOCORRELATION_TPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AutoCorrelation.h"  
#include "AutoCorrStage.tpp"  

#include <string>

namespace Util
{

   /*
   * Constructor
   *
   * \param blockFactor ratio of block sizes of subsequent stages
   */
   template <typename Data, typename Product>
   AutoCorrelation<Data, Product>::AutoCorrelation(int blockFactor = 4)
    : AutoCorrStage(blockFactor)
   {}

   /*
   * Set the base output file name.
   */
   template <typename Data, typename Product>
   void 
   AutoCorrelation<Data, Product>::setOutputFileName(std::string outputFileName)
   {  outputFileName_ = outputFileName; }

   /*
   * Register the creation of a descendant stage.
   */
   template <typename Data, typename Product>
   virtual void 
   AutoCorrelation<Data, Product>::registerDescendant(AutoCorrelation* ptr)
   {}

   /*
   * Load state from an archive.
   */
   template <typename Data, typename Product>
   virtual void 
   AutoCorrelation<Data, Product>::loadParameters(Serializable::IArchive& ar)
   {
      AutoCorrStage::load(ar);
   }

   /*
   * Save state to an archive.
   */
   template <typename Data, typename Product>
   virtual void 
   AutoCorrelation<Data, Product>::save(Serializable::OArchive& ar)
   {
      AutoCorrStage::save(ar);
   }

   /*
   * Serialize to/from an archive.
   */
   template <class Archive>
   template <typename Data, typename Product>
   void 
   AutoCorrelation<Data, Product>::serialize(Archive& ar, 
                                             const unsigned int version)
   {
      AutoCorrStage::serialize(ar, version);
   }

   /*
   * Output the autocorrelation function
   */
   template <typename Data, typename Product>
   virtual void AutoCorrelation<Data, Product>::output()
   {}

}
#endif
