#ifndef DDMD_ACTOR_MANAGER_H
#define DDMD_ACTOR_MANAGER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Manager.h>        // base class template
#include <util/containers/GPArray.h>   // member template
#include "Actor.h"                     // template parameter

namespace DdMd
{

   using namespace Util;

   class Simulation;

   /**
   * Manager for a set of Actor objects.
   *
   * An ActorManager maintains a list of Actor objects and 
   * provides methods to execute specified actions at various 
   * points before and during the main integration loop. 
   *
   * \ingroup DdMd_Manager_Module
   */
   class ActorManager : public Manager<Actor>
   {

   public:

      /**
      * Default constructor.
      */
      ActorManager();

      /**
      * Constructor.
      */
      ActorManager(Simulation& simulation);

      /**
      * Destructor.
      */
      ~ActorManager();

      /**
      * Read parameter file. 
      *
      * \param in input parameter file stream.
      */
      void readParam(std::istream &in);

      /**
      * Return pointer to a new default factory.
      *
      * Virtual, inherited from Manager<Actor>.
      */
      Factory<Actor>* newDefaultFactory() const;

      /// \name Setup actions 
      //@{ 
   
      void setupPostExchange();
      void setupPostNeighbor();
      void setupPostForce();
  
      //@} 
      /// \name Integration actions (within main loop)
      //@{ 

      void preIntegrate1(long iStep);
      void postIntegrate1(long iStep);
      void preTransform(long iStep);
      void preExchange(long iStep);
      void postExchange(long iStep);
      void postNeighbor(long iStep);
      void preUpdate(long iStep);
      void postUpdate(long iStep);
      void preForce(long iStep);
      void postForce(long iStep);
      void endOfStep(long iStep);

      //@} 
      /// \name Interprocessor communication actions 
      //@{ 

      void packExchange(long iStep);
      void unpackExchange(long iStep);
      void packUpdate(long iStep);
      void unpackUpdate(long iStep);
      void packReverseUpdate(long iStep);
      void unpackReverseUpdate(long iStep);

      //@} 
      
   private:

      /// Pointer to parent Simulation.
      Simulation* simulationPtr_;

      // Arrays of actors for specific actions.
      GPArray<Actor> setupPostExchangeActors_;
      GPArray<Actor> setupPostNeighborActors_;
      GPArray<Actor> setupPostForceActors_;
      GPArray<Actor> preIntegrate1Actors_;
      GPArray<Actor> postIntegrate1Actors_;
      GPArray<Actor> preTransformActors_;
      GPArray<Actor> preExchangeActors_;
      GPArray<Actor> postExchangeActors_;
      GPArray<Actor> postNeighborActors_;
      GPArray<Actor> preUpdateActors_;
      GPArray<Actor> postUpdateActors_;
      GPArray<Actor> preForceActors_;
      GPArray<Actor> postForceActors_;
      GPArray<Actor> endOfStepActors_;
      GPArray<Actor> exchangeActors_;
      GPArray<Actor> updateActors_;
      GPArray<Actor> reverseUpdateActors_;

   };

}
#endif
