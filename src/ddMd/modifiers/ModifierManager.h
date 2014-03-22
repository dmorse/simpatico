#ifndef DDMD_MODIFIER_MANAGER_H
#define DDMD_MODIFIER_MANAGER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Manager.h>        // base class template
#include <util/containers/GPArray.h>   // member template
#include "Modifier.h"                  // template parameter

namespace DdMd
{

   using namespace Util;

   class Simulation;

   /**
   * Manager for a set of Modifier objects.
   *
   * An ModifierManager maintains a list of Modifier objects and 
   * provides methods to execute specified actions at various 
   * points before and during the main integration loop. 
   *
   * \ingroup DdMd_Manager_Module
   * \ingroup DdMd_Modifier_Module
   */
   class ModifierManager : public Manager<Modifier>
   {

   public:

      /**
      * Default constructor (for unit testing).
      */
      ModifierManager();

      /**
      * Constructor (for use in Simulation).
      */
      ModifierManager(Simulation& simulation);

      /**
      * Destructor.
      */
      ~ModifierManager();

      /**
      * Read main block of parameter file. 
      *
      * \param in input parameter file stream.
      */
      void readParameters(std::istream &in);

      /**
      * Return pointer to a new default factory.
      *
      * Virtual, inherited from Manager<Modifier>.
      */
      Factory<Modifier>* newDefaultFactory() const;

      /// \name Integrator actions 
      //@{ 
   
      /**
      * Setup before entering the main loop.
      */ 
      void setup();

      /** 
      * Call just before the first step of velocity-Verlet algorithm. 
      */
      void preIntegrate1(long iStep);

      /** 
      * Call just after the first step of velocity-Verlet algorithm. 
      */
      void postIntegrate1(long iStep);

      /** 
      * Call on exchange steps before transforming to scaled coordinates.
      */
      void preTransform(long iStep);

      /** 
      * Call on exchange steps after transforming, before exchanging.
      */
      void preExchange(long iStep);

      /** 
      * Call on exchange steps after atom exchange, before reneighboring
      */
      void postExchange(long iStep);

      /** 
      * Call on exchange steps after re-building the neighbor list.
      */
      void postNeighbor(long iStep);

      /** 
      * Call on update steps before updating ghost positions.
      */
      void preUpdate(long iStep);

      /** 
      * Call on update steps after updating ghost positions.
      */
      void postUpdate(long iStep);

      /** 
      * Call after updating but before calculating forces.
      */
      void preForce(long iStep);

      /** 
      * Call after calculating forces
      */
      void postForce(long iStep);

      /** 
      * Call after 2nd integration step, at end of the time step.
      */
      void endOfStep(long iStep);

      //@} 
      /// \name Communication functions
      //@{ 

      /**
      * Pack data into buffer used to exchange atoms.
      */
      void packExchange(long iStep);

      /**
      * Unpack data from buffer used to exchange atoms.
      */
      void unpackExchange(long iStep);

      /**
      * Pack data into buffer used to update ghost positions.
      */
      void packUpdate(long iStep);

      /**
      * Unpack data from buffer used to update ghost positions.
      */
      void unpackUpdate(long iStep);

      /**
      * Pack data into buffer used to reverse update forces.
      *
      * Used only if reverse communication is enabled.
      */
      void packReverseUpdate(long iStep);

      /**
      * Unpack data from the buffer used to reverse update forces.
      *
      * Used only if reverse communication is enabled.
      */
      void unpackReverseUpdate(long iStep);

      //@} 
      
   private:

      /// Pointer to parent Simulation.
      Simulation* simulationPtr_;

      /// Arrays of pointers to modifiers for specific actions.
      GPArray<Modifier> setupModifiers_;
      GPArray<Modifier> preIntegrate1Modifiers_;
      GPArray<Modifier> postIntegrate1Modifiers_;
      GPArray<Modifier> preTransformModifiers_;
      GPArray<Modifier> preExchangeModifiers_;
      GPArray<Modifier> postExchangeModifiers_;
      GPArray<Modifier> postNeighborModifiers_;
      GPArray<Modifier> preUpdateModifiers_;
      GPArray<Modifier> postUpdateModifiers_;
      GPArray<Modifier> preForceModifiers_;
      GPArray<Modifier> postForceModifiers_;
      GPArray<Modifier> endOfStepModifiers_;
      GPArray<Modifier> exchangeModifiers_;
      GPArray<Modifier> updateModifiers_;
      GPArray<Modifier> reverseUpdateModifiers_;

   };

}
#endif
