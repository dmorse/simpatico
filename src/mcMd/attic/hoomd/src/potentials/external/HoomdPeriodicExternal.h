#ifdef SIMP_EXTERNAL
#ifndef HOOMD_PERIODIC_EXTERNAL_H
#define HOOMD_PERIODIC_EXTERNAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdExternal.h"

#include <hoomd/PotentialExternal.h>
#include <hoomd/EvaluatorPeriodicExternal.h>
#include <hoomd/AllDriverPotentialExternalGPU.cuh>

namespace McMd
{
   extern char classNameHoomdPeriodic[];

   /**
   * A potential encapsulating the HOOMD Periodic evaluator
   *
   * \ingroup External_Module
   */
   class HoomdPeriodicExternal : public HoomdExternal< EvaluatorPeriodicExternal,
      gpu_compute_periodic_forces, classNameHoomdPeriodic >
   {
   
   public:

      /**
      * Default constructor.
      */
      HoomdPeriodicExternal();

      /**
      * Copy constructor
      */
      HoomdPeriodicExternal(const HoomdPeriodicExternal& other);

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

      /// per-type prefactor of potential
      DArray<double> prefactor_;

      /// external parameter
      double externalParameter_;

      /// Array of Miller index IntVectors for the reciprocal lattice vectors.
      DArray<IntVector>  waveIntVectors_;

      double interfaceWidth_;

      int periodicity_;
   };
  
}

#endif
#endif
