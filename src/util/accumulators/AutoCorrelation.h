#ifndef UTIL_AUTOCORRELATION_H
#define UTIL_AUTOCORRELATION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AutoCorrStage.h"            // base class
#include <util/containers/GArray.h>   // member
#include <util/global.h>

namespace Util
{

   /**
   * Auto-correlation function, using hierarchical algorithm.
   *
   * This class represents the primary stage of a linked list of
   * AutoCorrelation objects.
   *
   * \ingroup Accumulators_Module
   */
   template <typename Data, typename Product>
   class AutoCorrelation : public AutoCorrStage<Data, Product>
   {

   public:

      /**
      * Constructor
      */
      AutoCorrelation();

      /**
      * Output the autocorrelation function
      *
      * \param out output stream.
      */
      virtual void output(std::ostream& out);

      /**
      * Register a descendant stage.
      *
      * This should be called only by a root stage.
      *
      * \param ptr pointer to a descendant AutoCorrelation.
      */
      virtual void registerDescendant(AutoCorrStage<Data, Product>* ptr);

   private:

      // Pointers to descendant AutoCorrStage objects
      GArray< AutoCorrStage<Data, Product>* > descendants_;

   };

}
#endif
