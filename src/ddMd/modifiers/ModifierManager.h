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
#include "Modifier.h"                     // template parameter

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
      * Default constructor.
      */
      ModifierManager();

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

      /**
      * Return pointer to a new default factory.
      *
      * Virtual, inherited from Manager<Modifier>.
      */
      Factory<Modifier>* newDefaultFactory() const;

      /// \name Integrator actions 
      //@{ 
   
      void setup();
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

      // Arrays of modifiers for specific actions.
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
