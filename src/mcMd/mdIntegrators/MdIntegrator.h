#ifndef MCMD_MD_INTEGRATOR_H
#define MCMD_MD_INTEGRATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <util/boundary/Boundary.h>       // typedef

#include <iostream>

namespace McMd
{

   using namespace Util;

   class Simulation;
   class MdSystem;
   
   /**
   * Abstract base for molecular dynamics integrators.
   *
   * \ingroup McMd_MdIntegrator_Module
   */
   class MdIntegrator : public ParamComposite
   {
   
   public:

      /**
      * Constructor. 
      *
      * \param system parent MdSystem object
      */
      MdIntegrator(MdSystem& system);

      /** 
      * Destructor.
      */
      virtual ~MdIntegrator();

      /**
      * Initialize internal state, if any.
      *
      * This method should be called just before entering
      * the main MD loop. Empty default implementation.
      */
      virtual void setup()
      {}

      /**
      * Take a complete MD integration step.
      */
      virtual void step() = 0;

      /**
      * Get Boundary of parent System by reference.
      */
      Boundary& boundary();

      /**
      * Get parent MdSystem by reference.
      */
      MdSystem& system();

      /**
      * Get parent Simulation by reference.
      */
      Simulation& simulation();

   protected:

      /// Integrator time step 
      double dt_;

   private:

      /// Pointer to a Boundary.
      Boundary* boundaryPtr_;
   
      /// Pointer to parent MdSystem
      MdSystem* systemPtr_;

      /// Pointer to parent Simulation.
      Simulation* simulationPtr_;

   }; 

   // Inline methods

   /*
   * Get Boundary of parent System.
   */
   inline Boundary&   MdIntegrator::boundary()
   {  return *boundaryPtr_; }

   /*
   * Get parent MdSystem.
   */
   inline MdSystem&   MdIntegrator::system()
   {  return *systemPtr_; }

   /*
   * Get parent Simulation.
   */
   inline Simulation& MdIntegrator::simulation()
   {  return *simulationPtr_; }

} 
#endif
