#ifndef HOOMD_MOVE_CPP
#define HOOMD_MOVE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdMove.h"
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/simulation/Simulation.h>
#include <util/ensembles/BoundaryEnsemble.h>
#include <mcMd/chemistry/Atom.h>
#include <simp/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/potentials/pair/McPairPotential.h>
#include <mcMd/potentials/bond/BondFactory.h>
#include <mcMd/potentials/bond/BondPotential.h>
#include <mcMd/potentials/bond/BondPotentialImpl.h>
#include <simp/interaction/bond/HarmonicBond.h>
#include <simp/interaction/bond/HarmonicL0Bond.h>
#include <simp/interaction/bond/FeneBond.h>

#include <modules/hoomd/potentials/pair/HoomdPairFactory.h>
#include <modules/hoomd/potentials/pair/HoomdPairPotential.h>
#include <modules/hoomd/potentials/pair/HoomdPair.h>
#include <modules/hoomd/potentials/bond/HoomdBond.h>
#ifdef SIMP_EXTERNAL
#include <mcMd/potentials/external/ExternalPotential.h>
#include <modules/hoomd/potentials/external/HoomdExternalFactory.h>
#include <modules/hoomd/potentials/external/HoomdExternal.h>
#endif
#if defined(HOOMD_DEVEL) && defined(MCMD_LINK)
#include <mcMd/potentials/link/HoomdLinkBondPotential.h>
#include <mcMd/potentials/link/LinkFactory.cpp>
#endif

#include <hoomd/HOOMDMath.h>
#include <hoomd/ParticleData.h>
#include <hoomd/ForceCompute.h>
#include <hoomd/AllBondPotentials.h>
#include <hoomd/ParticleGroup.h>
#include <hoomd/IntegratorTwoStep.h>
#include <hoomd/TwoStepNVEGPU.h>
#include <hoomd/TwoStepNPTMTKGPU.h>
#include <hoomd/SFCPackUpdater.h>
#include <hoomd/NeighborList.h>


namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /*
   * Constructor
   */
   HoomdMove::HoomdMove(McSystem& system) :
      SystemMove(system),
      GPUId_(0),
      nStep_(0),
      dt_(0),
      skin_(0),
      moleculeSetHasChanged_(false),
      HoomdIsInitialized_(false),
      sorterPeriod_(300)
   {
      #if defined(MCMD_LINK) && !defined(HOOMD_DEVEL)
      UTIL_THROW("Links not supported by HoomdMove (only with HOOMD_DEVEL)");
      #endif
   }

   /*
   * Destructor.
   */
   HoomdMove::~HoomdMove()
   {}

   /*
   * Read parameters
   */
   void HoomdMove::readParameters(std::istream& in)
   {
      if ((! energyEnsemble().isIsothermal()) || (!system().boundaryEnsemble().isRigid()))
         UTIL_THROW("Must be in isothermal-rigid ensemble.");

      readProbability(in);
      read<int>(in, "nStep", nStep_);
      read<double>(in, "dt", dt_);
      read<double>(in, "skin", skin_);
      char* env;
      if ((env = getenv("OMPI_COMM_WORLD_LOCAL_RANK")) != NULL) {
         GPUId_ = atoi(env);
      } else {
         GPUId_ = -1;
      }

      // read<int>(in, "GPUId", GPUId_);
      // create HOOMD execution configuration
      executionConfigurationSPtr_ = 
        boost::shared_ptr<ExecutionConfiguration>(new ExecutionConfiguration(ExecutionConfiguration::GPU,GPUId_));

      // set CUDA error checking
      #ifdef UTIL_DEBUG
      executionConfigurationSPtr_->setCUDAErrorChecking(true);
      #endif

      // subscribe to molecule set notifications
      system().subscribeMoleculeSetChange(*this);

      #if defined(HOOMD_DEVEL) && defined(MCMD_LINK)
      // subscribe to link notifications
      system().linkMaster().Notifier<LinkAddEvent>::registerObserver(*this);
      system().linkMaster().Notifier<LinkResetEvent>::registerObserver(*this);
      system().linkMaster().Notifier<LinkRemoveEvent>::registerObserver(*this);
      #endif 

   }

   #if defined(MCMD_LINK) && defined(HOOMD_DEVEL)
   /**
   * add all active links as bonds
   *
   * \note: this function must only called once, i.e. after addBonds(),
   *        to ensure that the same link is not added multiple times
   *
   * after initialization, insertion/modification/removal of links
   * is taken care of by Link{Add,Reset,Remove}Events
   */
   void HoomdMove::addLinks()
   {
      linkBondMap_.clear();

      // Loop over bonds and initialize them
      int nBondType = simulation().nBondType();
     
      LinkMaster::ConstLinkIterator it;
      system().linkMaster().begin(it);
      for (; it.notEnd(); ++it) {
         // link type ids are mapped according to
         // HOOMD bond id =  nBondType + link type id
         ::Bond bond(nBondType+it->typeId(),
                     it->atom0().id(),
                     it->atom1().id());
         bondDataSPtr_->addBond(bond);
         // keep track of bond
         unsigned int bondTag = bondDataSPtr_->getLastTag();
         int linkTag = it->tag();
         #ifdef UTIL_DEBUG
         std::pair<boost::unordered_map<int, unsigned int>::iterator,bool> ret;
         ret=
         #endif
         linkBondMap_.insert(std::pair<int,unsigned int>(linkTag,bondTag));
         assert(ret.second);
      }

   }
   #endif

   #ifndef HOOMD_DEVEL
   /**
   * return a ForceCompute corresponding to the BondPotential
   */
   // FIXME: uses default GPU block_sizes (not optimized values)
   boost::shared_ptr<ForceCompute> HoomdMove::mapBondPotential()
   {
      if (system().bondPotential().interactionClassName() == "HarmonicBond")
      {
         BondPotentialImpl<HarmonicBond> *harmonicBondPtr =
            dynamic_cast< BondPotentialImpl<HarmonicBond> *>
            (&system().bondPotential());
         boost::shared_ptr<PotentialBondHarmonic>
            force(new PotentialBondHarmonicGPU(systemDefinitionSPtr_,""));
         for (int iType=0; iType<simulation().nBondType(); iType++)
         {
            Scalar2 param = make_scalar2(
                             harmonicBondPtr->interaction().stiffness(iType),
                             harmonicBondPtr->interaction().length(iType));
            force->setParams(iType, param);
         }
         bondName_ = std::string("harmonic");
         return force;
      } else
      if (system().bondPotential().interactionClassName() == "HarmonicL0Bond")
      {
          BondPotentialImpl<HarmonicL0Bond> *harmonicL0BondPtr =
             dynamic_cast< BondPotentialImpl<HarmonicL0Bond> *>
             (&system().bondPotential());
          boost::shared_ptr<PotentialBondHarmonic>
             force(new PotentialBondHarmonicGPU(systemDefinitionSPtr_,""));
          for (int iType=0; iType<simulation().nBondType(); iType++)
          {
             Scalar2 param = make_scalar2(
                             harmonicL0BondPtr->interaction().stiffness(iType),
                             0);
             force->setParams(iType, param);
          }
          bondName_ = std::string("harmonic");
          return force;
       } else
        UTIL_THROW("Bond potential not (yet) supported by HoomdMove.");
   }
   #endif

   #if defined(HOOMD_DEVEL) && defined(MCMD_LINK)
   /**
   * connection to LinkAddEvent
   */ 
   void HoomdMove::update(const LinkAddEvent& event)
   {
      if (!HoomdIsInitialized_) return;

      // Add a link as a HOOMD bond
      Link *linkPtr = event.get();

      ::Bond bond(simulation().nBondType()+linkPtr->typeId(),
                  linkPtr->atom0().id(),
                  linkPtr->atom1().id());
      bondDataSPtr_->addBond(bond);
      // keep track of bond
      unsigned int bondTag = bondDataSPtr_->getLastTag();
      int linkTag = linkPtr->tag();
      #ifdef UTIL_DEBUG
      std::pair<boost::unordered_map<int, unsigned int>::iterator,bool> ret;
      ret=
      #endif
      linkBondMap_.insert(std::pair<int,unsigned int>(linkTag,bondTag));
      assert(ret.second);
   }

   /**
   * connection to LinkResetEvent
   */ 
   void HoomdMove::update(const LinkResetEvent& event)
   {
      if (!HoomdIsInitialized_) return;

      Link *linkPtr = event.get();

      // First remove the link from the HOOMD bond table
      int linkTag = linkPtr->tag();
      boost::unordered_map<int, unsigned>::iterator it;
      it = linkBondMap_.find(linkTag);
      assert(it != linkBondMap_.end());
      unsigned int oldBondTag = it ->second;
      bondDataSPtr_->removeBond(oldBondTag);

      // Add the modified link 
      // link type ids are mapped according to
      // HOOMD bond id =  nBondType + link type id
      ::Bond bond(simulation().nBondType()+linkPtr->typeId(),
         linkPtr->atom0().id(),
         linkPtr->atom1().id());
      bondDataSPtr_->addBond(bond);

      // Update link-bond mapping
      unsigned int newBondTag = bondDataSPtr_->getLastTag();
      it->second = newBondTag;
   }

   /**
   * connection to LinkRemoveEvent
   */ 
   void HoomdMove::update(const LinkRemoveEvent& event)
   {
      if (!HoomdIsInitialized_) return;

      Link *linkPtr = event.get();

      // Remove bond from HOOMD
      int linkTag = linkPtr->tag();
      boost::unordered_map<int, unsigned>::iterator it;
      it = linkBondMap_.find(linkTag);
      assert(it != linkBondMap_.end());
      unsigned int bondTag = it ->second;
      bondDataSPtr_->removeBond(bondTag);

      // Remove link-bond mapping
      linkBondMap_.erase(it); 
   }
   #endif  // HOOMD_DEVEL && MCMD_LINK

   void HoomdMove::initSimulation()
   {
      // set up simulation box
      BoxDim box;
      Boundary boundary = system().boundary();
      lengths_ = boundary.lengths();
      const Scalar3 boundaryLengths = make_scalar3(lengths_[0], lengths_[1], lengths_[2]);
      box.setL(boundaryLengths);
      //box.xlo = -lengths_[0]/2.; box.xhi = lengths_[0]/2.;
      //box.ylo = -lengths_[1]/2.; box.yhi = lengths_[1]/2.;
      //box.zlo = -lengths_[2]/2.; box.zhi = lengths_[2]/2.;

      // create HOOMD system
      systemDefinitionSPtr_ =
         boost::shared_ptr<SystemDefinition>(new SystemDefinition(
         system().nAtom(),             // # atoms
         box,                        // box dimensions
         simulation().nAtomType(),   // # atom types
         #if defined(HOOMD_DEVEL) && defined(MCMD_LINK)
         simulation().nLinkType()+
         #endif
         simulation().nBondType(),   // # bond types
         0,                          // # angle types
         0,                          // # dihedral types
         0,                          // # improper types
         executionConfigurationSPtr_ // execution configuration
         ));

      // access Bond data
      bondDataSPtr_ = systemDefinitionSPtr_->getBondData();
      // add bonds
      addBonds();

      #if defined(HOOMD_DEVEL) && defined(MCMD_LINK)
      // add link
      addLinks();
      #endif

      // Try to obtain ForceCompute from pair potential
      hoomdPairPtr_ = HoomdPairFactory::hoomdPairPotentialPtr(
         system().pairPotential());
      if (!hoomdPairPtr_)
         UTIL_THROW("No Hoomd pair potential defined.");

      pairForceSPtr_ = HoomdPairFactory::hoomdFactory(system().pairPotential(),
         system(), systemDefinitionSPtr_,nListSPtr_,skin_);

      #ifdef SIMP_EXTERNAL
      externalForceSPtr_ = HoomdExternalFactory::hoomdFactory(system().externalPotential(),system(),
         systemDefinitionSPtr_);
      implementExternalPotential_ = system().implementExternalPotential();
      #endif 

      #ifdef HOOMD_DEVEL
      #ifdef MCMD_LINK
      // Create link potential
      hoomdBondPtr_ =
         system().linkFactory().linkBondFactory(
         system().bondPotential().interactionClassName(),
         system().linkPotential().interactionClassName());
      if (!hoomdBondPtr_)
         UTIL_THROW("No (known) combination of HOOMD link and bond potentials defined.");
      bondForceSPtr_ = hoomdBondPtr_->forceCompute(systemDefinitionSPtr_,
         system());
      #else
      // Try to obtain ForceCompute from bond potential
      hoomdBondPtr_ = system().bondFactory().convertPotential(
         system().bondPotential());; 
      if (!hoomdBondPtr_)
         UTIL_THROW("No (known) HOOMD bond potential defined.");

      bondForceSPtr_ = hoomdBondPtr_->forceCompute(
         systemDefinitionSPtr_,system()); 
      #endif // MCMD_LINK
      #else
      bondForceSPtr_ = mapBondPotential();
      #endif // HOOMD_DEVEL

      // group of all particles
      boost::shared_ptr<ParticleSelectorTag> sel(
         new ParticleSelectorTag(systemDefinitionSPtr_, 0, system().nAtom()-1));
      groupAllSPtr_ = boost::shared_ptr<ParticleGroup>(
         new ParticleGroup(systemDefinitionSPtr_, sel));

      // create ComputeThermo
      thermoSPtr_ = boost::shared_ptr<ComputeThermoGPU>(new
         ComputeThermoGPU(systemDefinitionSPtr_, groupAllSPtr_, ""));

      // get ParticleData
      particleDataSPtr_ = systemDefinitionSPtr_->getParticleData();

      // register sorter
      sorterSPtr_ = boost::shared_ptr<SFCPackUpdater>(new SFCPackUpdater(
         systemDefinitionSPtr_));

      HoomdIsInitialized_ = true;
   }
 
   void HoomdMove::createIntegrator()     
   {

      // create NVE Integrator 
      integratorSPtr_ = boost::shared_ptr<IntegratorTwoStep>(new IntegratorTwoStep(systemDefinitionSPtr_,dt_));

      // create IntegrationMethod
      boost::shared_ptr<TwoStepNVEGPU> methodSPtr(
         new TwoStepNVEGPU(systemDefinitionSPtr_,groupAllSPtr_));

      integratorSPtr_->addIntegrationMethod(methodSPtr);

      // register pair and bond forces in Integrator
      integratorSPtr_->addForceCompute(pairForceSPtr_);
      integratorSPtr_->addForceCompute(bondForceSPtr_);
      #ifdef SIMP_EXTERNAL
      if (implementExternalPotential_)
         integratorSPtr_->addForceCompute(externalForceSPtr_);
      #endif


      // set flags for thermodynamic quantities requested by integrator
      // additionally, we require computation of potential energy
      PDataFlags peFlag;
      peFlag[pdata_flag::potential_energy] = 1;
      particleDataSPtr_->setFlags(integratorSPtr_->getRequestedPDataFlags() | peFlag);
   }

   void HoomdMove::addBonds()
   {
      // Loop over bonds and initialize them
      int iSpec;
      System::MoleculeIterator molIter;
      Molecule::BondIterator bondIter;
      for (iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
         for (system().begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
            for (molIter->begin(bondIter); bondIter.notEnd(); ++bondIter) {
               ::Bond bond(bondIter->typeId(),
                           bondIter->atom(0).id(),
                           bondIter->atom(1).id());
               bondDataSPtr_->addBond(bond);
             }
         }
      } 
   }

   /*
    * Generate random velocities
    */
   void HoomdMove::generateRandomVelocities(ArrayHandle<Scalar4> h_vel)
   {
      double temp = energyEnsemble().temperature();
      unsigned int nAtom = (unsigned int) system().nAtom();
      Random& random = simulation().random();

      Scalar3 totalMomentum = make_scalar3(0.0,0.0,0.0);

      for (unsigned int idx=0; idx<nAtom; idx++) {
         Scalar mass = h_vel.data[idx].w;
         Scalar sigma = sqrt(temp / mass);
         Scalar vx = sigma*random.gaussian();
         Scalar vy = sigma*random.gaussian();
         Scalar vz = sigma*random.gaussian();
        
         totalMomentum.x += vx * mass;
         totalMomentum.y += vy * mass;
         totalMomentum.z += vz * mass;

         h_vel.data[idx].x = vx;
         h_vel.data[idx].y = vy;
         h_vel.data[idx].z = vz; 
      }

      // loop through the particles again and remove the system momentum
      totalMomentum.x /= nAtom;
      totalMomentum.y /= nAtom;
      totalMomentum.z /= nAtom;
      for (unsigned int idx = 0; idx < nAtom; idx++)
      {
         Scalar mass = h_vel.data[idx].w;
         h_vel.data[idx].x -= totalMomentum.x / mass;
         h_vel.data[idx].y -= totalMomentum.y / mass;
         h_vel.data[idx].z -= totalMomentum.z / mass;
      }
   }

   /*
    * The molecule set change signal is used to reinitialize the MD simulation
    */
   void HoomdMove::notifyMoleculeSetChanged()
   {
      moleculeSetHasChanged_ = true;
   }

   /*
   * Generate, attempt and accept or reject a Hybrid MD/MC move.
   */
   bool HoomdMove::move()
   {
      if ((!HoomdIsInitialized_) || moleculeSetHasChanged_) {
         initSimulation();
         moleculeSetHasChanged_ = false;
      }

      // We need to create the Integrator every time since we are starting
      // with new coordinates, velocities etc.
      // this does not seem to incur a significant performance decrease
      createIntegrator();

      System::MoleculeIterator molIter;
      Molecule::AtomIterator   atomIter;
      int nSpec = simulation().nSpecies();

      // Increment counter for attempted moves
      incrementNAttempt();

      double oldEnergy, newEnergy;

      {
      // Copy atom coordinates into HOOMD
      ArrayHandle<Scalar4> h_pos(particleDataSPtr_->getPositions(), access_location::host, access_mode::readwrite);
      ArrayHandle<Scalar4> h_vel(particleDataSPtr_->getVelocities(), access_location::host, access_mode::readwrite);
      ArrayHandle<unsigned int> h_tag(particleDataSPtr_->getTags(), access_location::host, access_mode::readwrite);
      ArrayHandle<unsigned int> h_rtag(particleDataSPtr_->getRTags(), access_location::host, access_mode::readwrite);

      for (int iSpec =0; iSpec < nSpec; ++iSpec) {
         system().begin(iSpec, molIter);
         for ( ; molIter.notEnd(); ++ molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               unsigned int idx = (unsigned int) atomIter->id();
               Vector& pos = atomIter->position();
               h_pos.data[idx].x = pos[0] - lengths_[0]/2.;
               h_pos.data[idx].y = pos[1] - lengths_[1]/2.;
               h_pos.data[idx].z = pos[2] - lengths_[2]/2.;
 
               int type = atomIter->typeId();
               h_vel.data[idx].w = simulation().atomType(type).mass();
               h_pos.data[idx].w = __int_as_scalar(type);
               h_tag.data[idx] = idx;
               h_rtag.data[idx] = idx; 
            }
         }
      }
 
      // Generate random velocities
      generateRandomVelocities(h_vel);
      }


      // Notify that the particle order has changed
      particleDataSPtr_->notifyParticleSort();

      // Initialize integrator (calculate forces and potential energy for step 0)
      integratorSPtr_->prepRun(0);

      // Calculate oldEnergy
      thermoSPtr_->compute(0);
      oldEnergy = thermoSPtr_->getLogValue("kinetic_energy",0);
      oldEnergy += thermoSPtr_->getLogValue("potential_energy",0);

      // Integrate nStep_ steps forward
      for (int iStep = 0; iStep < nStep_; ++iStep) {
         integratorSPtr_->update(iStep);

         // do we need to sort the particles?
         // do not sort at time step 0 to speed up short runs
         if (! (iStep % sorterPeriod_) && iStep)
           sorterSPtr_->update(iStep);
      }


      // Calculate new energy
      thermoSPtr_->compute(nStep_);
      newEnergy = thermoSPtr_->getLogValue("kinetic_energy",nStep_);
      newEnergy += thermoSPtr_->getLogValue("potential_energy",nStep_);

      // Decide whether to accept or reject
      bool accept = random().metropolis( boltzmann(newEnergy-oldEnergy) );

      if (accept) {
         // read back integrated positions 
         ArrayHandle<Scalar4> h_pos(particleDataSPtr_->getPositions(), access_location::host, access_mode::read);
         ArrayHandle<Scalar4> h_vel(particleDataSPtr_->getVelocities(), access_location::host, access_mode::read);
         ArrayHandle<unsigned int> h_tag(particleDataSPtr_->getTags(), access_location::host, access_mode::read);
         ArrayHandle<unsigned int> h_rtag(particleDataSPtr_->getRTags(), access_location::host, access_mode::read);
         for (int iSpec = 0; iSpec < nSpec; ++iSpec) {
            system().begin(iSpec, molIter);
            for ( ; molIter.notEnd(); ++molIter) {
               for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
                  unsigned int idx = h_rtag.data[atomIter->id()]; 
                  atomIter->position() = Vector(h_pos.data[idx].x+lengths_[0]/2.,
                                                h_pos.data[idx].y+lengths_[1]/2.,
                                                h_pos.data[idx].z+lengths_[2]/2.);
               }
            }
         }

         system().pairPotential().buildCellList();
         incrementNAccept();
      } else {
         // not accepted, do nothing
      }

      return accept;
   }
 
}
#endif
