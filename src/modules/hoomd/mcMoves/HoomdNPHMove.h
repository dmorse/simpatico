#ifndef HOOMD_NPH_MOVE_H
#define HOOMD_NPH_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdMove.h"

#include <hoomd/TwoStepNPHGPU.h>

namespace McMd
{

   using namespace Util;

   /**
   * HoomdNPHMove is a hybrid Molecular Dynamics MC move using Hoomd-blue
   * 
   * It samples the NPT (isobaric-isothermal) ensemble by integrating the NPH (isoenthalpic-isobaric,
   * Andersen barostat) equations of motion. The box dimensions are updated after an accepted HoomdNPHMove.
   *
   * \ingroup McMove_Module
   */
   class HoomdNPHMove : public HoomdMove
   {

   public:

      /**
      * Constructor.
      *
      * Constructs a component MdSystem object.
      */
      HoomdNPHMove(McSystem& system);

      /**
      * Destructor.
      */
      ~HoomdNPHMove();

      /**
      * Read nStep, dt, skin, maxNPair from file.
      */
      virtual void readParameters(std::istream& in);

      /**
      * Generate, attempt and accept or reject a move.
      */
      bool move();

   protected:
      /// Create NPH integrator
      void createIntegrator();

      /// Andersen barostat mass
      double W_;
     
      /// anisotropic integration mode (can be cubic, orthorhombic or tetragonal)
      TwoStepNPHGPU::integrationMode integrationMode_; 

      /// integration method
      boost::shared_ptr<TwoStepNPHGPU> twoStepNPHGPUSPtr_;

   private:
      
      /// Geometry of simulation cell
      std::string modeIn_;

      bool toImposeConstrain_;
      bool setConstrain_;
      Vector constrainLengths_;
   };

}
#endif
