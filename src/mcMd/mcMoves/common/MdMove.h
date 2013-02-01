#ifndef MCMD_MD_MOVE_H
#define MCMD_MD_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcMoves/SystemMove.h>  // base class
#include <util/containers/DArray.h>   // member template
#include <util/space/Vector.h>         // member template parameter

namespace McMd
{

   using namespace Util;

   class McSystem;
   class MdSystem;
   

   /**
   * MdMove is a simple NVE molecular Dynamics MC move.
   *
   * An MdMove simply runs a (typically) short NVE simulation with initial
   * velocities chosen from a Boltzmann distribution. All attempted moves
   * are accepted.
   *
   * MdMove rigorously satisfies detailed balance only in the limit of
   * vanishing step size dt. The closely related HybridMdMove does satisfy
   * detailed balance, for all dt, but generally accepts a fraction of 
   * attempted moves. 
   *
   * \ingroup McMd_McMove_Module MD_Module
   */
   class MdMove : public SystemMove 
   {
   
   public:
   
      /**
      * Constructor. 
      *
      * Constructs a component MdSystem object.
      */
      MdMove(McSystem& system);
   
      /**
      * Destructor.
      */
      ~MdMove();
   
      /**
      * Read nStep, dt, skin, maxNPair from file.
      *
      * \param in input parameter stream.
      */
      virtual void readParameters(std::istream& in);
   
      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

      /**
      * Generate, attempt and accept or reject a move.
      */
      bool move();
   
   private:
  
      /// MdSystem object used for MD integration
      MdSystem  *mdSystemPtr_;  
   
      /// Number of Md steps per Hybrid MD move
      int  nStep_;

   };

}      
#endif
