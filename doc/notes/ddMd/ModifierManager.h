#ifndef DDMD_MODIFIER_MANAGER_H
#define DDMD_MODIFIER_MANAGER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Modifier.h"               // template parameter
#include <util/param/Manager.h>     // base class template

namespace DdMd
{

   using namespace Util;

   class Simulation;

   /**
   * Manager for a list of Modifier objects.
   *
   * \ingroup DdMd_Manager_Module
   */
   class ModifierManager : public Manager<Modifier>
   {

   public:

      /**
      * Constructor.
      */
      ModifierManager(Simulation& simulation);

      /**
      * Destructor.
      */
      ~ModifierManager();

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
      void setupEnd();
   
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
      void postFinalIntegrate(long iStep);
   
      // Communication
      void pack_exchange(Buffer& buffer);
      void unpack_exchange(Buffer& buffer);
      void pack_update(Buffer& buffer);
      void unpack_update(Buffer& buffer);
      void pack_reverseUpdate(Buffer& buffer);
      void unpack_reverseUpdate(Buffer& buffer);

      /**
      * Return pointer to a new default factory.
      */
      Factory<Modifier>* newDefaultFactory() const;

   private:

      /// Pointer to parent Simulation.
      Simulation* simulationPtr_;
 
   };

}
#endif
