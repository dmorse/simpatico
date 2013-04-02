#ifdef INTER_EXTERNAL
#ifndef HOOMD_PERIODIC_EXTERNAL_H
#define HOOMD_PERIODIC_EXTERNAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdExternal.h"

#include <hoomd/PotentialExternal.h>
#include <hoomd/EvaluatorExternalPeriodic.h>
#include <hoomd/AllDriverPotentialExternalGPU.cuh>

namespace McMd
{
   extern char classNameHoomdLamellarOrdering[];

   /**
   * A potential encapsulating the HOOMD Periodic evaluator
   *
   * \ingroup External_Module
   */
   class HoomdLamellarOrderingExternal : public HoomdExternal< EvaluatorExternalPeriodic,
      gpu_compute_periodic_forces, classNameHoomdLamellarOrdering >
   {
   
   public:

      /**
      * Default constructor.
      */
      HoomdLamellarOrderingExternal();

      /**
      * Copy constructor
      */
      HoomdLamellarOrderingExternal(const HoomdLamellarOrderingExternal& other);

      /**
      * read parameters from file
      *
      * \param in input stream
      */
      void readParameters(std::istream &in);

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

      /// Interfacial width in periodic phase.
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
