#ifndef SLIPLINK_MC_MOVE_FACTORY_H
#define SLIPLINK_MC_MOVE_FACTORY_H

#include <mcMd/mcMoves/McMoveFactory.h>

namespace McMd
{

   using namespace Util;

   class McSimulation;
   class McSystem;

   /**
   * Custom McMoveFactory.
   */
   class SliplinkMcMoveFactory : public McMoveFactory
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent simulation
      * \param system     parent system
      */
      SliplinkMcMoveFactory(McSimulation& simulation, McSystem& system);

      /** 
      * Return pointer to a new McMove object.
      *
      * \param  className name of a subclass of McMove.
      * \return base class pointer to a new instance of className.
      */
      virtual McMove* factory(const std::string& className) const;

   };

}
#endif
