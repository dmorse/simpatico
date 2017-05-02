#ifdef SIMP_EXTERNAL
#ifndef HOOMD_LOCAL_LAMELLAR_ORDERING_EXTERNAL_H
#define HOOMD_LOCAL_LAMELLAR_ORDERING_EXTERNAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdExternal.h"

#include <hoomd/PotentialExternal.h>
#include <hoomd/LocalExternalParams.h>
#include <hoomd/EvaluatorLocalExternal.h>
#include <hoomd/AllDriverPotentialExternalGPU.cuh>

namespace McMd
{
   extern char classNameHoomdLocalLamellarOrdering[];

   /**
   * A potential encapsulating the HOOMD Periodic evaluator
   *
   * \ingroup External_Module
   */
   class HoomdLocalLamellarOrderingExternal : public HoomdExternal< EvaluatorLocalExternal,
      gpu_compute_local_forces, classNameHoomdLocalLamellarOrdering >
   {
   
   public:

      /**
      * Default constructor.
      */
      HoomdLocalLamellarOrderingExternal();

      /**
      * Copy constructor
      */
      HoomdLocalLamellarOrderingExternal(const HoomdLocalLamellarOrderingExternal& other);

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

      /// Index representing the direction parallel to the lamellae.
      int parallelDirection_;

      /// Fraction of length in parallel direction for which lamellar ordering is imposed.
      double fraction_;

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
