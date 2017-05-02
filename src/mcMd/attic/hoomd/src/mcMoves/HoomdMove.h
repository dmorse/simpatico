#ifndef HOOMD_MOVE_H
#define HOOMD_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/simulation/System.h>   // defines MoleculeSetObserver
#include <mcMd/mcMoves/SystemMove.h>  // base class
#include <util/containers/DArray.h>   // member template
#include <util/space/Vector.h>        // member template parameter

#if defined(HOOMD_DEVEL) && defined(MCMD_LINK)
#include <mcMd/links/LinkMaster.h>
#include <mcMd/links/LinkEvents.h>
#include <util/misc/Observer.h>
#include <boost/unordered_map.hpp>
#endif

#include <boost/shared_ptr.hpp>

#include <boost/hoomd_config.h>       // HOOMD
#include <hoomd/SystemDefinition.h>
#include <hoomd/System.h>
#include <hoomd/ExecutionConfiguration.h>
#include <hoomd/ParticleData.h>
#include <hoomd/IntegratorTwoStep.h>
#include <hoomd/ParticleGroup.h>
#include <hoomd/ComputeThermoGPU.h>
#include <hoomd/SFCPackUpdater.h>

namespace McMd
{

   using namespace Util;

   class McSystem;
   class MdSystem;
   class HoomdPairPotential;
   #ifdef HOOMD_DEVEL
   class HoomdBondPotential;
   #endif
   #ifdef SIMP_EXTERNAL
   class HoomdExternalPotential;
   #endif

   /**
   * HoomdMove is a hybrid Molecular Dynamics MC move that uses HOOMD
   *
   * \ingroup McMove_Module
   */
   class HoomdMove : public SystemMove,
   #if defined(MCMD_LINK) && defined(HOOMD_DEVEL)
                     public Observer<LinkAddEvent>,
                     public Observer<LinkResetEvent>,
                     public Observer<LinkRemoveEvent>,
   #endif
                     public MoleculeSetObserver
   {

   public:

      /**
      * Constructor.
      *
      * Constructs a component MdSystem object.
      */
      HoomdMove(McSystem& system);

      /**
      * Destructor.
      */
      ~HoomdMove();

      /**
      * Read nStep, dt, skin, maxNPair from file.
      */
      virtual void readParameters(std::istream& in);

      /**
      * Generate, attempt and accept or reject a move.
      */
      bool move();

      /**
      * connection to moleculeSet change signal
      */
      virtual void notifyMoleculeSetChanged();

      #if defined(HOOMD_DEVEL) && defined(MCMD_LINK)
      /**
      * Connection to LinkAddEvent
      *
      * \param event description of event
      */
      virtual void update(const LinkAddEvent& event);

      /**
      * Connection to LinkResetEvent
      *
      * \param event description of event
      */
      virtual void update(const LinkResetEvent& event);

      /**
      * Connection to LinkRemoveEvent
      *
      * \param event description of event
      */
      virtual void update(const LinkRemoveEvent& event);
      #endif

   protected:

      /// HOOMD SystemDefinition object used for MD integration
      boost::shared_ptr<SystemDefinition>  systemDefinitionSPtr_;

      /// HOOMD SystemDefinition object used for MD integration
      boost::shared_ptr<ExecutionConfiguration>  executionConfigurationSPtr_;

      /// HOOMD particle data
      boost::shared_ptr< ParticleData > particleDataSPtr_;

      /// HOOMD bond data
      boost::shared_ptr< BondData > bondDataSPtr_;

      /// HOOMD integrator
      boost::shared_ptr< IntegratorTwoStep > integratorSPtr_;

      /// HOOMD ComputeThermo
      boost::shared_ptr< ComputeThermoGPU > thermoSPtr_;

      /// pair potential
      HoomdPairPotential *hoomdPairPtr_ ;

      #ifdef HOOMD_DEVEL
      /// bond potential
      HoomdBondPotential *hoomdBondPtr_;
      #endif

      /// HOOMD pair potential
      boost::shared_ptr< ForceCompute > pairForceSPtr_;

      #ifdef SIMP_EXTERNAL
      bool implementExternalPotential_;
      /// HOOMD external potential
      boost::shared_ptr< ForceCompute > externalForceSPtr_;
      #endif 

      /// HOOMD neighbor list
      boost::shared_ptr< NeighborList > nListSPtr_;

      /// HOOMD bond potential
      boost::shared_ptr< ForceCompute > bondForceSPtr_;

      #ifndef HOOMD_DEVEL
      /// internal name of HOOMD bond potential
      std::string bondName_;
      #endif

      /// HOOMD integration method
      boost::shared_ptr< ParticleGroup > groupAllSPtr_;

      /// HOOMD sorter
      boost::shared_ptr< SFCPackUpdater > sorterSPtr_;

      /// GPU to execute on
      int GPUId_ ;

      /// Number of Md steps per Hybrid MD move
      int nStep_;

      /// Time step for HOOMD integrator
      double dt_;

      /// Skin length
      double skin_;

      /// has the MoleculeSet changed since last move()?
      bool moleculeSetHasChanged_;

      /// is HOOMD initialized?
      bool HoomdIsInitialized_;

      /// the period with which the sorter is called
      unsigned int sorterPeriod_;

      /// Boundary lengths
      Vector lengths_;

      /// map BondPotential onto ForceCompute
      boost::shared_ptr<ForceCompute> mapBondPotential();

      /// register bonds
      void addBonds();

      #if defined(HOOMD_DEVEL) && defined(MCMD_LINK)

      /**
      * Register links as bonds.
      *
      * Mutable links are added as regular HOOMD bonds.
      *
      * \note
      * The HOOMD bond table is kept in sync with the links registerd in the
      * LinkMaster via Link{Add,Reset,Remove}Events.
      *
      * - Links are assigned HOOMD bond type ids starting with the first free
      *   bond type id (i.e. nBondType)
      *
      * - logic in \sa HoomdLinkBondPotential::forceCompute() assures that
      *   the HOOMD bond potential is initialized with the parameters of the
      *   link potential for nBondType <= HOOMD bond id <  nBondType+nLinkType
      *
      * This method adds all defined links as bonds at once.
      */
      void addLinks();

      /// Map to translate between link::tag() and HOOMD bond tags
      boost::unordered_map<int, unsigned int> linkBondMap_;
      #endif

      /// Initialize HOOMD simulation
      void initSimulation();

      /// Create integrator
      void createIntegrator();

      /// Generate random velocities
      void generateRandomVelocities(ArrayHandle<Scalar4> h_vel);

   };

}
#endif
