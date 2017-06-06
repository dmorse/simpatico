#ifndef HOOMD_NPT_MTK_MOVE_CPP
#define HOOMD_NPT_MTK_MOVE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdNPTMTKMove.h"

#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/simulation/Simulation.h>
#include <simp/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <mcMd/potentials/pair/McPairPotential.h>
#include <util/ensembles/BoundaryEnsemble.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor
   */
   HoomdNPTMTKMove::HoomdNPTMTKMove(McSystem& system) :
      HoomdMove(system), tauP_(0.0)
   {
   }

   /*
   * Destructor.
   */
   HoomdNPTMTKMove::~HoomdNPTMTKMove()
   {}

   /*
   * Read parameters
   */
   void HoomdNPTMTKMove::readParameters(std::istream& in)
   {
      if ((! system().boundaryEnsemble().isIsobaric()) || (!energyEnsemble().isIsothermal()))
         UTIL_THROW("Must be in isothermal-isobaric ensemble.");

      readProbability(in);
      read<int>(in, "nStep", nStep_);
      read<double>(in, "dt", dt_);

      read<std::string>(in, "couplingMode", couplingModeIn_);
      if (couplingModeIn_ == "couple_none")
         couplingMode_ = TwoStepNPTMTKGPU::couple_none;
      else if (couplingModeIn_ == "couple_xy") 
         couplingMode_ = TwoStepNPTMTKGPU::couple_xy;
      else if (couplingModeIn_ == "couple_xz")
         couplingMode_ = TwoStepNPTMTKGPU::couple_xz;
      else if (couplingModeIn_ == "couple_yz")
         couplingMode_ = TwoStepNPTMTKGPU::couple_yz;
      else if (couplingModeIn_ == "couple_xyz")
         couplingMode_ = TwoStepNPTMTKGPU::couple_xyz;
      else
         UTIL_THROW("Unsupported coupling mode.");

      // Should a constraint be applied on aspect ratio in the case of a tetragonal cell
      read<bool>(in, "imposeConstrain", imposeConstrain_);
      if (imposeConstrain_) {
         // Aspect ratio defined as length in unique direction over the other length
         read<double>(in, "minAspectRatio", minAspectRatio_);
         read<double>(in, "maxAspectRatio", maxAspectRatio_);
      }

      // Barostat mass W ~ 3 N T tauP^2
      read<double>(in,"tauP", tauP_);
      char* env;
      if ((env = getenv("OMPI_COMM_WORLD_LOCAL_RANK")) != NULL) {
         GPUId_ = atoi(env);
      } else {
         GPUId_ = -1;
      }

      //toImposeConstrain_ = false;
      //setConstrain_ = false;

      //read<int>(in, "GPUId", GPUId_);
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
 
   void HoomdNPTMTKMove::createIntegrator()     
   {

      // create NVE Integrator 
      integratorSPtr_ = boost::shared_ptr<IntegratorTwoStep>(new IntegratorTwoStep(systemDefinitionSPtr_,dt_));

      boost::shared_ptr<Variant> variantTSPtr(new VariantConst(1.0));
      boost::shared_ptr<Variant> variantPSPtr(new VariantConst(system().boundaryEnsemble().pressure()));

      Species  *speciesPtr;
      int       iSpec, nMolecule;
      int       nAtom = 0;
      for (int iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
         speciesPtr = &simulation().species(iSpec);
         nMolecule  = system().nMolecule(iSpec);
         nAtom += nMolecule*(speciesPtr->nAtom());
      }
      thermoSPtr_->setNDOF(3.0*nAtom);

      // create IntegrationMethod
      twoStepNPTMTKGPUSPtr_ = boost::shared_ptr<TwoStepNPTMTKGPU>(new TwoStepNPTMTKGPU(systemDefinitionSPtr_, groupAllSPtr_,                                                                  thermoSPtr_, 1.0, tauP_, variantTSPtr, variantPSPtr, 
                                                                  couplingMode_, TwoStepNPTMTKGPU::baro_x | TwoStepNPTMTKGPU                                                                  ::baro_y | TwoStepNPTMTKGPU::baro_z, true));

      integratorSPtr_->addIntegrationMethod(twoStepNPTMTKGPUSPtr_);

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

   /*
   * Generate, attempt and accept or reject a Hybrid MD/MC move.
   */
   bool HoomdNPTMTKMove::move()
   {
      if ((!HoomdIsInitialized_) || moleculeSetHasChanged_) {
         initSimulation();
         moleculeSetHasChanged_ = false;
      }

      /*
      if ( !setConstrain_ && toImposeConstrain_) {
         Boundary initBoundary = system().boundary();
         constrainLengths_ = initBoundary.lengths();
         setConstrain_ = true;
      }
      */

      // We need to create the Integrator every time since we are starting
      // with new coordinates, velocities etc.
      // this does not seem to incur a significant performance decrease
      createIntegrator();

      System::MoleculeIterator molIter;
      Molecule::AtomIterator   atomIter;
      int nSpec = simulation().nSpecies();

      // Increment counter for attempted moves
      incrementNAttempt();

      // set hoomd simulation box
      BoxDim box;
      Boundary boundary = system().boundary();
      lengths_ = boundary.lengths();
      const Scalar3 boundaryLengths = make_scalar3(lengths_[0], lengths_[1], lengths_[2]);
      box.setL(boundaryLengths);

      particleDataSPtr_->setGlobalBoxL(boundaryLengths);
      {
      // copy atom coordinates into hoomd
      ArrayHandle<Scalar4> h_pos(particleDataSPtr_->getPositions(), access_location::host, access_mode::readwrite);
      ArrayHandle<Scalar4> h_vel(particleDataSPtr_->getVelocities(), access_location::host, access_mode::readwrite);
      ArrayHandle<unsigned int> h_tag(particleDataSPtr_->getTags(), access_location::host, access_mode::readwrite);
      ArrayHandle<unsigned int> h_rtag(particleDataSPtr_->getRTags(), access_location::host, access_mode::readwrite);
      int nind = 0;
      for (int iSpec =0; iSpec < nSpec; ++iSpec) {
         system().begin(iSpec, molIter);
         for ( ; molIter.notEnd(); ++ molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               nind = nind + 1;
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

      // Done with the data
      // generate integrator variables from a Gaussian distribution
      double etax = 0.0;
      double etay = 0.0;
      double etaz = 0.0;
      Random& random = simulation().random();
      double temp = energyEnsemble().temperature();
      Species  *speciesPtr;
      int       nMolecule;
      int       nAtom = 0;
      for (int iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
         speciesPtr = &simulation().species(iSpec);
         nMolecule  = system().nMolecule(iSpec);
         nAtom += nMolecule*(speciesPtr->nAtom());
      }
      if (couplingMode_ == TwoStepNPTMTKGPU::couple_xyz) {
         // one degree of freedom
         // barostat_energy = 1/2 (1/W) eta_x^2
         double sigma = sqrt( temp/(tauP_*tauP_*temp*3.0*nAtom) );
         etax = sigma*random.gaussian();
         etay = etax;
         etaz = etax;
      } else if (couplingMode_ == TwoStepNPTMTKGPU::couple_none) {
         // three degrees of freedom 
         // barostat_energy = 1/2 (1/W) (eta_x^2 + eta_y^2 + eta_z^2)
         double sigma = sqrt( temp/(tauP_*tauP_*temp*3.0*nAtom) );
         etax = sigma*random.gaussian();
         etay = sigma*random.gaussian();
         etaz = sigma*random.gaussian();
      } else {
         // two degrees of freedom
         // barostat_energy = 1/2 (1/W) eta_x^2 + 1/2 (1/(2W)) eta_y^2
         double sigma1 = sqrt( temp/(tauP_*tauP_*temp*3.0*nAtom) );
         double eta1 = sigma1*random.gaussian();
         double sigma2 = sqrt( temp/(tauP_*tauP_*temp*3.0*nAtom)/2.0 );
         double eta2 = sigma2*random.gaussian();
         if (couplingMode_ == TwoStepNPTMTKGPU::couple_xy) {
            etax = eta2;
            etay = eta2;
            etaz = eta1;
         } else
         if (couplingMode_ == TwoStepNPTMTKGPU::couple_yz) {
            etax = eta1;
            etay = eta2;
            etaz = eta2;
         } else 
         if (couplingMode_ == TwoStepNPTMTKGPU::couple_xz) {
            etax = eta2;
            etay = eta1;
            etaz = eta2;
         } else { 
           UTIL_THROW("Unsupported coupling mode.");
         }
      }
      twoStepNPTMTKGPUSPtr_->setEta(etax,etay,etaz);

      // Notify that the particle order might have changed
      particleDataSPtr_->notifyParticleSort();

      // Calculate oldH (the conserved quantity of the Andersen barostat)
      double volume = lengths_[0]*lengths_[1]*lengths_[2];

      // Initialize integrator (calculate forces and potential energy for step 0)
      integratorSPtr_->prepRun(0);

      // H = U + PV + barostat_energy
      thermoSPtr_->compute(0);
      double oldH = thermoSPtr_->getLogValue("kinetic_energy",0);
      oldH += thermoSPtr_->getLogValue("potential_energy",0);
      oldH += system().boundaryEnsemble().pressure() * volume;
      oldH += integratorSPtr_->getLogValue("npt_barostat_energy",0);
      //std::cout << "oldH kinetic is " << thermoSPtr_->getLogValue("kinetic_energy",0) << std::endl;
      //std::cout << "oldH potential is " << thermoSPtr_->getLogValue("potential_energy",0) << std::endl;
      //std::cout << "oldH pressure is " << system().boundaryEnsemble().pressure() * volume << std::endl;
      //std::cout << "oldH barostat is " << integratorSPtr_->getLogValue("npt_barostat_energy",0) << std::endl;

      // Integrate nStep_ steps forward
      for (int iStep = 0; iStep < nStep_; ++iStep) {
         integratorSPtr_->update(iStep);
         //std::cout << "energy at iStep is " << externalForceSPtr_->getLogValue("external_periodic_energy", iStep) << std::endl; 

         // do we need to sort the particles?
         // do not sort at time step 0 to speed up short runs
         if (! (iStep % sorterPeriod_) && iStep)
           sorterSPtr_->update(iStep);
      }

      Vector newLengths;
      box = particleDataSPtr_->getBox();
      newLengths[0] = (box.getL()).x;
      newLengths[1] = (box.getL()).y;
      newLengths[2] = (box.getL()).z;
      volume = newLengths[0]*newLengths[1]*newLengths[2];

      // Calculate new value of the conserved quantity
      thermoSPtr_->compute(nStep_);
      double newH = thermoSPtr_->getLogValue("kinetic_energy",nStep_);
      newH += thermoSPtr_->getLogValue("potential_energy",nStep_);
      newH += system().boundaryEnsemble().pressure() * volume;
      newH += integratorSPtr_->getLogValue("npt_barostat_energy",nStep_);
      //std::cout << "newH kinetic is " << thermoSPtr_->getLogValue("kinetic_energy",nStep_) << std::endl;
      //std::cout << "newH potential is " << thermoSPtr_->getLogValue("potential_energy",nStep_) << std::endl;
      //std::cout << "newH pressure is " << system().boundaryEnsemble().pressure() * volume << std::endl;
      //std::cout << "newH barostat is " << integratorSPtr_->getLogValue("npt_barostat_energy",nStep_) << std::endl;

      bool accept;

      // Decide whether to accept or reject
      if (imposeConstrain_) {
         double  newAspectRatio;
         if (couplingMode_ == TwoStepNPTMTKGPU::couple_xy) {
            // Aspect ratio defined as length in unique direction over the other length
            newAspectRatio = double(newLengths[2]/newLengths[0]);
            if ( newAspectRatio < minAspectRatio_ || newAspectRatio > maxAspectRatio_) {
               accept = false;
            } else {
               accept = random.metropolis( boltzmann(newH-oldH) );
            }
         } 
         if (couplingMode_ == TwoStepNPTMTKGPU::couple_xz) {
            // Aspect ratio defined as length in unique direction over the other length
            newAspectRatio = double(newLengths[1]/newLengths[0]);
            if ( newAspectRatio < minAspectRatio_ || newAspectRatio > maxAspectRatio_) {
               accept = false;
            } else {
               accept = random.metropolis( boltzmann(newH-oldH) );
            }
         } 
         if (couplingMode_ == TwoStepNPTMTKGPU::couple_yz) {
            // Aspect ratio defined as length in unique direction over the other length
            newAspectRatio = double(newLengths[0]/newLengths[1]);
            if ( newAspectRatio < minAspectRatio_ || newAspectRatio > maxAspectRatio_) {
               accept = false;
            } else {
               accept = random.metropolis( boltzmann(newH-oldH) );
            }
         }
      } else { 
         accept = random.metropolis( boltzmann(newH-oldH) );
      }

      if (accept) {
         // read back new boundary
         box = particleDataSPtr_->getBox();
         lengths_ = newLengths;
         system().boundary().setOrthorhombic(lengths_);

         // read back integrated positions
         ArrayHandle<Scalar4> h_pos(particleDataSPtr_->getPositions(), access_location::host, access_mode::read);
         ArrayHandle<Scalar4> h_vel(particleDataSPtr_->getVelocities(), access_location::host, access_mode::read);
         ArrayHandle<unsigned int> h_tag(particleDataSPtr_->getTags(), access_location::host, access_mode::read);
         ArrayHandle<unsigned int> h_rtag(particleDataSPtr_->getRTags(), access_location::host, access_mode::read);
         //int nind = 0; 
         for (int iSpec = 0; iSpec < nSpec; ++iSpec) {
            system().begin(iSpec, molIter);
            for ( ; molIter.notEnd(); ++molIter) {
               for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
                  //nind = nind + 1;
                  unsigned int idx = h_rtag.data[atomIter->id()]; 
                  atomIter->position() = Vector(h_pos.data[idx].x+lengths_[0]/2.,
                                                h_pos.data[idx].y+lengths_[1]/2.,
                                                h_pos.data[idx].z+lengths_[2]/2.);
                  //if(nind == 85408) {
                    // std::cout << "accept hpos is " << h_pos.data[idx].x << std::endl;
                    // std::cout << "accept position is " << atomIter->position() << std::endl;
                  //}
               }
            }
         }
         // Done with data
         system().pairPotential().buildCellList();
         incrementNAccept();
      } else {
          // not accepted, do nothing
         // read back integrated positions
         int nind = 0;
         for (int iSpec = 0; iSpec < nSpec; ++iSpec) {
            system().begin(iSpec, molIter);
            for ( ; molIter.notEnd(); ++molIter) {
               for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
                  nind = nind + 1;
                  /*
                  if(nind == 1) {
                     std::cout << "reject index is " << atomIter->id() << std::endl;
                     std::cout << "reject position is " << atomIter->position() << std::endl;
                     std::cout << "lengths is " << lengths_ << std::endl;
                  }*/
               }
            }
         }

      }

//      std::cout << "Lx = " << lengths_[0] << " Ly = " << lengths_[1] << " Lz = " << lengths_[2] << std::endl;
      return accept;
   }
 
}
#endif
