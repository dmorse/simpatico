#ifndef MCMD_HYBRID_NPH_MD_MOVE_H
#define MCMD_HYBRID_NPH_MD_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcMoves/SystemMove.h>  // base class
#include <util/containers/DArray.h>   // member template
#include <util/space/Vector.h>         // member template parameter
#include <mcMd/mdIntegrators/NphIntegrator.h>         

namespace McMd
{

   using namespace Util;

   class McSystem;
   class MdSystem;
   

   /**
   * HybridMdMove is a hybrid Molecular Dynamics MC move.
   *
   * \ingroup McMd_McMove_Module MD_Module
   */
   class HybridNphMdMove : public SystemMove 
   {
   
   public:
   
      /**
      * Constructor. 
      *
      * \param system parent McSystem.
      */
      HybridNphMdMove(McSystem& system);
   
      /**
      * Destructor.
      */
      ~HybridNphMdMove();
   
      /**
      * Read nStep and MdSystem parameters from file.
      *
      * \param in input parameter stream
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
      MdSystem      *mdSystemPtr_;  
   
      /// Array to store old atomic positions.
      DArray<Vector> oldPositions_; 

      /// Number of Md steps per Hybrid MD move
      int            nStep_;
      
      /// Pointer to NphIntegrator
      NphIntegrator* nphIntegratorPtr_;

      /// Mass of integrator barostat 
      double barostatMass_;

      /// Integration mode
      LatticeSystem mode_;

   };

}      
#endif
