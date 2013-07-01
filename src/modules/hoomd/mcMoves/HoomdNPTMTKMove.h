#ifndef HOOMD_NPT_MTK_MOVE_H
#define HOOMD_NPT_MTK_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdMove.h"

#include <hoomd/TwoStepNPTMTKGPU.h>

namespace McMd
{

   using namespace Util;

   /**
   * HoomdNPTMTKMove is a hybrid Molecular Dynamics MC move using Hoomd-blue
   * 
   * It samples the NPT (isobaric-isothermal) ensemble by integrating the NPH (isoenthalpic-isobaric,
   * Andersen barostat) equations of motion. The box dimensions are updated after an accepted HoomdNPTMTKMove.
   *
   * \ingroup McMove_Module
   */
   class HoomdNPTMTKMove : public HoomdMove
   {

   public:

      /**
      * Constructor.
      *
      * Constructs a component MdSystem object.
      */
      HoomdNPTMTKMove(McSystem& system);

      /**
      * Destructor.
      */
      ~HoomdNPTMTKMove();

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
      double tauP_;
     
      /// anisotropic integration mode (can be cubic, orthorhombic or tetragonal)
      TwoStepNPTMTKGPU::couplingMode couplingMode_; 

      /// integration method
      boost::shared_ptr<TwoStepNPTMTKGPU> twoStepNPTMTKGPUSPtr_;

   private:
      
      /// Geometry of simulation cell
      std::string couplingModeIn_;

      /// Geometry of simulation cell
      bool imposeConstrain_;

      /// Geometry of simulation cell
      double minAspectRatio_;
     
      /// Geometry of simulation cell
      double maxAspectRatio_;
   };

}
#endif
