#ifndef MCMD_CFB_REPTATION_MOVE_H
#define MCMD_CFB_REPTATION_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcMoves/base/CfbEndBase.h>   // base class
#include <mcMd/chemistry/MaskPolicy.h>      // member
#include <util/containers/DArray.h>         // members
#include <util/accumulators/AutoCorr.h>

namespace McMd
{

   using namespace Util;

   class McSystem;
   class Molecule;

   /**
   * Configuration bias reptation move for a Linear species.
   *
   * A reptation move uses a configuration bias algorithm to add a segment 
   * to one end of a chain (the head) and delete one from the other (the 
   * tail).  The direction of the move (i.e., which end is the head) is 
   * chosen at random. The only parameter is nTrial, the number of trial 
   * positions for the new position of the head monomer.
   *
   * \ingroup McMd_McMove_Module
   */
   class CfbReptationMove : public CfbEndBase
   {
   
   public:
   
      /**
      * Constructor. 
      */
      CfbReptationMove(McSystem& system);
   
      /**
      * Read species to which displacement is applied.
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
      * Generate and accept or reject configuration bias move
      */
      virtual bool move();

      /**
      * Output statistics about accepted reptation steps
      */
      virtual void output();
   
   private:
   
      /// Integer index for molecular species.
      int speciesId_;
  
      /// Type id for bonds (must be the same for all bonds).
      int bondTypeId_;
  
      /// Number of junctions.
      int nJunction_;

      /// Array of ids for lower junction atoms.
      DArray<int>  junctions_;

      /// Array of types for lower junction atoms.
      DArray<int>  lTypes_;

      /// Array of types for upper junction atoms.
      DArray<int>  uTypes_;

      /// Policy for masking of nonbonded interactions between bonded atoms.
      MaskPolicy maskPolicy_;

      /// Do we calculate the autocorrelation of accepted MC steps?
      int hasAutoCorr_;

      /// Capacity of accepted step autocorrelation
      int autoCorrCapacity_;

      /// Filename for autocorrelation output
      std::string outputFileName_;

      /// Array of autocorrelators for accepted steps
      DArray<AutoCorr<double, double> > acceptedStepsAccumulators_;

      double junctionFactor(Molecule* molPtr, int sign);

   };

}      
#endif
