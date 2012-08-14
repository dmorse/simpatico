#ifndef DDMD_ACTOR_MANAGER_H
#define DDMD_ACTOR_MANAGER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Actor.h"               // template parameter
#include <util/param/Manager.h>   // base class template

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

      // Setup
 
      /// Execute all re-implemented setupPostExchange() methods.
      void setupPostExchange();

      /// Execute all re-implemented setupPostNeighbor() methods.
      void setupPostNeighbor();

      /// Execute all re-implemented setupPostForce() methods.
      void setupPostForce();
   
      // Integration
      void preIntegrate(long iStep);
      void postIntegrate(long iStep);
   
      void preTransform(long iStep);
      void preExchange(long iStep);
      void postExchange(long iStep);
      void postNeighbor(long iStep);
   
      void preUpdate(long iStep);
      void postUpdate(long iStep);
   
      void preForce(long iStep);
      void postForce(long iStep);
      void endOfStep(long iStep);
   
      // Communication
      void packExchange(long iStep);
      void unpackExchange(long iStep);

      void packUpdate(long iStep);
      void unpackUpdate(long iStep);

      void packReverseUpdate(long iStep);
      void unpackReverseUpdate(long iStep);

      /**
      * Return pointer to a new default factory.
      *
      * Virtual, inherited from Manager<Actor>.
      */
      Factory<Actor>* newDefaultFactory() const;

   private:

      /// Pointer to parent Simulation.
      Simulation* simulationPtr_;

      // Arrays of actions with setup methods 
      std::vector<Actor*> actionsSetupPostExchange_;
      std::vector<Actor*> actionsSetupPostNeighbor_;
      std::vector<Actor*> actionsSetupPostForce_;

      // Arrays of actions with integrate methods 
      std::vector<Actor*> actionsPreIntegrate_;
      std::vector<Actor*> actionsPostIntegrate_;
      std::vector<Actor*> actionsPreTransform_;
      std::vector<Actor*> actionsPreExchange_;
      std::vector<Actor*> actionsPostExchange_;
      std::vector<Actor*> actionsPostNeighbor_;
      std::vector<Actor*> actionsPreUpdate_;
      std::vector<Actor*> actionsPostUpdate_;
      std::vector<Actor*> actionsPreForce_;
      std::vector<Actor*> actionsPostForce_;
      std::vector<Actor*> actionsEndOfStep_;

      // Arrays of actions with communication methods 
      std::vector<Actor*> actionsPackExchange_;
      std::vector<Actor*> actionsUnpackExchange_;
      std::vector<Actor*> actionsPackUpdate_;
      std::vector<Actor*> actionsUnpackUpdate_;
      std::vector<Actor*> actionsPackReverseUpdate_;
      std::vector<Actor*> actionsUnpackReverseUpdate_;

   };

}
#endif
