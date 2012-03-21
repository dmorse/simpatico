#ifndef MCMD_MC_MOVE_FACTORY_H
#define MCMD_MC_MOVE_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>
#include <mcMd/mcMoves/McMove.h>

namespace McMd
{

   using namespace Util;

   class McSimulation;
   class McSystem;

   /**
   * McMoveFactory for an McSimulation
   *
   * \ingroup Factory_Module
   * \ingroup McMove_Module
   */
   class McMoveFactory : public Factory<McMove>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent simulation
      * \param system     parent system
      */
      McMoveFactory(McSimulation& simulation, McSystem& system)
       : simulationPtr_(&simulation),
         systemPtr_(&system)
      {}

      /** 
      * Return pointer to a new McMove object.
      *
      * \param  className name of a subclass of McMove.
      * \return base class pointer to a new instance of className.
      */
      virtual McMove* factory(const std::string& className) const;

   protected:

      /**
      * Return reference to parent Simulation.
      */
      McSimulation& simulation() const
      { return *simulationPtr_; }

      /**
      * Return reference to parent McSystem.
      */
      McSystem& system() const
      { return *systemPtr_; }

   private:

      McSimulation *simulationPtr_;
      McSystem     *systemPtr_;

   };

}
#endif
