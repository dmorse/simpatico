#ifndef MCMD_SYSTEM_MOVE_H
#define MCMD_SYSTEM_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcMoves/McMove.h>      // base class
#include <util/ensembles/EnergyEnsemble.h> // inline function
#include <util/boundary/Boundary.h>        // typedef

namespace McMd
{

   using namespace Util;

   class McSystem;

   /**
   * An McMove that acts on one McSystem.
   *
   * \ingroup McMd_McMove_Module
   */
   class SystemMove : public McMove
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param system parent McSystem
      */
      SystemMove(McSystem& system);
   
      /**
      * Destructor.
      */
      virtual ~SystemMove();
   
   protected:

      /// Get parent McSystem.
      McSystem& system();

      /// Get Boundary object of parent McSystem.
      Boundary& boundary();

      /// Get EnergyEnsemble object of parent McSystem.
      EnergyEnsemble& energyEnsemble();

      /**
      * Boltzmann weight associated with an energy difference.
      *
      * \param energy energy or energy difference.
      * \return Boltzmann weight exp(-beta*energy)
      */
      double boltzmann(double energy);
  
   private:
 
      /// Pointer to parent McSystem object.
      McSystem  *systemPtr_;

      /// Pointer to Boundary of parent System.
      Boundary  *boundaryPtr_;

      /// Pointer to EnergyEnsemble of parent System.
      EnergyEnsemble*  isothermalPtr_;
 
   };

   // Inline methods

   /*
   * Get parent McSystem.
   */
   inline McSystem& SystemMove::system()
   {  return *systemPtr_; }

   /*
   * Get Boundary object of parent McSystem.
   */
   inline Boundary& SystemMove::boundary()
   {  return *boundaryPtr_; }

   /*
   * Get EnergyEnsemble object of parent McSystem.
   */
   inline EnergyEnsemble& SystemMove::energyEnsemble()
   {  return *isothermalPtr_; }

   /*
   * Boltzmann weight associated with an energy difference.
   */
   inline double SystemMove::boltzmann(double energy)
   {  return exp(-isothermalPtr_->beta()*energy); }
   
}
#endif
