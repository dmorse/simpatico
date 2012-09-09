#ifdef INTER_EXTERNAL
#ifndef HOOMD_LAMELLAR_EXTERNAL_H
#define HOOMD_LAMELLAR_EXTERNAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdExternal.h"

#include <hoomd/PotentialExternal.h>
#include <hoomd/EvaluatorExternalLamellar.h>
#include <hoomd/AllDriverPotentialExternalGPU.cuh>

namespace McMd
{
   extern char classNameHoomdLamellar[];

   /**
   * A potential encapsulating the HOOMD Lamellar evaluator
   *
   * \ingroup External_Module
   */
   class HoomdLamellarExternal : public HoomdExternal< EvaluatorExternalLamellar,
      gpu_compute_lamellar_forces, classNameHoomdLamellar >
   {
   
   public:

      /**
      * Default constructor.
      */
      HoomdLamellarExternal();

      /**
      * Copy constructor
      */
      HoomdLamellarExternal(const HoomdLamellarExternal& other);

      /**
      * read parameters from file
      *
      * \param in input stream
      */
      void readParam(std::istream &in);

      /**
      * Set external potential parameter
      *
      */
      void setExternalParameter(double externalParameter);

      /**
      * returns external potential parameter
      */
      double externalParameter() const;

      /**
      * returns the class name
      */
      std::string className() const;

   private:

      /// Index representing the direction perpendicular to the lamellae.
      int perpDirection_;

      /// Interfacial width in lamellar phase.
      double width_;

      /// per-type prefactor of potential
      DArray<double> prefactor_;

      /// external parameter
      double externalParameter_;

      /// Number of periods in a cell
      int periodicity_;

   };
  
}

#endif
#endif
