#ifndef MCMD_COLVAR_H
#define MCMD_COLVAR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>
#include <util/misc/Setable.h>

namespace McMd
{

   using namespace Util;

   /**
   * Base class for collective variables.
   *
   * Usage: Use the compute() function to compute and store a value, the
   * the unset() function to unset the stored value, and the value function
   * to retrieve the current value.
   *
   * \ingroup McMd_ColVar_Module
   */
   class ColVar : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      ColVar();

      /**
      * Destructor.
      *
      * Empty default implementation.
      */
      virtual ~ColVar();

      /**
      * Compute and store the value of the collective variable.
      */
      virtual void compute() = 0;

      /**
      * Unset the stored value (mark as unknown or obsolete).
      */
      void unset();

      /**
      * Return current value, compute if not already set.
      */
      virtual double value();

   protected:

      Setable<double> value_;

   };

}
#endif
