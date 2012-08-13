#ifndef DDMD_ACTION_MANAGER_H
#define DDMD_ACTION_MANAGER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Action.h"               // template parameter
#include <util/param/Manager.h>   // base class template

namespace DdMd
{

   using namespace Util;

   class Simulation;

   /**
   * Manager for a list of Action objects.
   *
   * \ingroup DdMd_Manager_Module
   */
   class ActionManager : public Manager<Action>
   {

   public:

      /**
      * Constructor.
      */
      ActionManager(Simulation& simulation);

      /**
      * Destructor.
      */
      ~ActionManager();

      /**
      * Read parameter file. 
      *
      * \param in input parameter file stream.
      */
      void readParam(std::istream &in);

      // Setup
   
      void setupPostExchange();
      void setupPostNeighbor();
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
      void packExchange(Buffer& buffer);
      void unpackExchange(Buffer& buffer);
      void packUpdate(Buffer& buffer);
      void unpackUpdate(Buffer& buffer);
      void packReverseUpdate(Buffer& buffer);
      void unpackReverseUpdate(Buffer& buffer);

      /**
      * Return pointer to a new default factory.
      */
      Factory<Action>* newDefaultFactory() const;

   private:

      /// Pointer to parent Simulation.
      Simulation* simulationPtr_;

      // Arrays of actions with setup methods 
      std::vector<Action*> actionsSetupPostExchange_;
      std::vector<Action*> actionsSetupPostNeighbor_;
      std::vector<Action*> actionsSetupPostForce_;

      // Arrays of actions with integrate methods 
      std::vector<Action*> actionsPreIntegrate_;
      std::vector<Action*> actionsPostIntegrate_;
      std::vector<Action*> actionsPreTransform_;
      std::vector<Action*> actionsPreExchange_;
      std::vector<Action*> actionsPostExchange_;
      std::vector<Action*> actionsPostNeighbor_;
      std::vector<Action*> actionsPreUpdate_;
      std::vector<Action*> actionsPostUpdate_;
      std::vector<Action*> actionsPreForce_;
      std::vector<Action*> actionsPostForce_;
      std::vector<Action*> actionsEndOfStep_;

      // Arrays of actions with communication methods 
      std::vector<Action*> actionsPackExchange_;
      std::vector<Action*> actionsUnpackExchange_;
      std::vector<Action*> actionsPackUpdate_;
      std::vector<Action*> actionsUnpackUpdate_;
      std::vector<Action*> actionsPackreverseUpdate_;
      std::vector<Action*> actionsUnpackreverseUpdate_;

   };

}
#endif
