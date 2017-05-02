/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "NvtDpdVvIntegrator.h"
#include <mcMd/mdSimulation/MdSystem.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/potentials/pair/MdPairPotential.h>
#include <util/ensembles/EnergyEnsemble.h>
#include <mcMd/neighbor/PairList.h>
#include <mcMd/neighbor/PairIterator.h>
#include <util/boundary/Boundary.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/archives/Serializable_includes.h>
#include <util/space/Vector.h>

#define MCMD_DPD_TYPE 0

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor.   
   */
   NvtDpdVvIntegrator::NvtDpdVvIntegrator(MdSystem& system)
   : MdIntegrator(system),
     dtMinvFactors_(),
     cutoff_(-1.0),
     gamma_(-1.0),
     sigma_(-1.0),
     temperature_(1.0),
     cutoffSq_(-1.0),
     #ifndef SIMP_NOPAIR
     pairListPtr_(&system.pairPotential().pairList()),
     #endif
     boundaryPtr_(&system.boundary()),
     randomPtr_(&system.simulation().random()),
     energyEnsemblePtr_(&system.energyEnsemble()),
     atomCapacity_(system.simulation().atomCapacity()),
     isInitialized_(false)
   {
      // Note: Within the constructor, the method parameter "system" 
      // hides the MdIntegrator::system() method name.

      // Precondition
      if (!system.energyEnsemble().isIsothermal() ) {
         UTIL_THROW("System energy ensemble is not isothermal");
      }
      temperature_  = energyEnsemblePtr_->temperature();

      setClassName("NvtDpdVvIntegrator"); 
   }

   /* 
   * Destructor.   
   */
   NvtDpdVvIntegrator::~NvtDpdVvIntegrator() 
   {}

   /* 
   * Read parameter and configuration files, initialize system.
   */
   void NvtDpdVvIntegrator::readParameters(std::istream &in) 
   {
      // Precondition
      if (isInitialized_) {
         UTIL_THROW("NvtDpdVvIntegrator already initialized");
      }

      read<double>(in, "dt", dt_);
      read<double>(in, "cutoff", cutoff_);
      #ifndef SIMP_NOPAIR
      if (cutoff_ > system().pairPotential().maxPairCutoff()) {
         UTIL_THROW("Error: Dpd pair cutoff > maxCutoff of pair potential");
      }
      #endif
      read<double>(in, "gamma", gamma_);

      // Allocate arrays for internal use
      dissipativeForces_.allocate(atomCapacity_);
      randomForces_.allocate(atomCapacity_);
      dtMinvFactors_.allocate(simulation().nAtomType());

      isInitialized_ = true;
   }

   /*
   * Load the internal state to an archive.
   */
   void NvtDpdVvIntegrator::loadParameters(Serializable::IArchive& ar)
   {  
      // Precondition
      if (isInitialized_) {
         UTIL_THROW("NvtDpdVvIntegrator already initialized");
      }

      loadParameter<double>(ar, "dt", dt_);
      loadParameter<double>(ar, "cutoff", cutoff_);
      #ifndef SIMP_NOPAIR
      if (cutoff_ > system().pairPotential().maxPairCutoff()) {
         UTIL_THROW("Error: Dpd pair cutoff > maxCutoff of pair potential");
      }
      #endif
      loadParameter<double>(ar, "gamma", gamma_);
      ar & temperature_;
      ar & sigma_;
      ar & cutoffSq_;
      ar & atomCapacity_;
      ar & dtMinvFactors_;
      ar & dissipativeForces_;
      ar & randomForces_;

      isInitialized_ = true;
   }

   /*
   * Save the internal state to an archive.
   */
   void NvtDpdVvIntegrator::save(Serializable::OArchive& ar)
   {
      ar & dt_;  
      ar & cutoff_;  
      ar & gamma_;  
      ar & temperature_;
      ar & sigma_;
      ar & cutoffSq_;
      ar & atomCapacity_;
      ar & dtMinvFactors_;
      ar & dissipativeForces_;
      ar & randomForces_;
   }

   /* 
   * Initialize constants and forces.
   */
   void NvtDpdVvIntegrator::setup() 
   {
      if (!isInitialized_) {
         UTIL_THROW("Object must be initialized before use");
      }

      // Calculate dtMinvFactors_ array.
      double mass;
      int nAtomType = simulation().nAtomType();
      for (int i = 0; i < nAtomType; ++i) {
         mass = simulation().atomType(i).mass();
         dtMinvFactors_[i] = 0.5*dt_/mass;
      }

      temperature_ = system().energyEnsemble().temperature();
      sigma_ = sqrt(2.0*gamma_*temperature_/dt_);
      cutoffSq_ = cutoff_*cutoff_;

      #ifndef SIMP_NOPAIR
      system().pairPotential().clearPairListStatistics();
      system().pairPotential().buildPairList();
      #endif

      system().calculateForces();
      computeDpdForces(true);
   }

   /*
   * Compute new dissipative and (optionally) random forces.
   *
   * Dissipative and random forces are stored in separate arrays. 
   * If computeRandom == true, recompute random as well as dissipative 
   * forces.
   */
   void NvtDpdVvIntegrator::computeDpdForces(bool computeRandom)
   {
      Vector e;   // unit vector (r0 - r1)/|r0 - r1|;
      Vector f;   // force vector 
      Vector dv;  // difference in velocities v0 - v1;
      double rsq; // square of distance between atoms.
      double r;   // distance between atoms.
      double fr;  // magnitude of random pair force.
      double fd;  // magnitude of dissipative pair force.
      double wr;  // weighting function for random forces.
      Atom* atom0Ptr;
      Atom* atom1Ptr;
      PairIterator iter;
      int i;

      // Set all dissipative forces to zero.
      for (i=0; i < atomCapacity_; ++i) {
         dissipativeForces_[i].zero();
      }
      if (computeRandom) {
         for (i=0; i < atomCapacity_; ++i) {
            randomForces_[i].zero();
         }
      }

      // Iterator over atom pairs
      for (pairListPtr_->begin(iter); iter.notEnd(); ++iter) {
         iter.getPair(atom0Ptr, atom1Ptr);
         rsq = boundaryPtr_->distanceSq(atom0Ptr->position(), 
                                        atom1Ptr->position(), e);
         if (rsq < cutoffSq_) {
            r = sqrt(rsq);
            e /= r;

            #if MCMD_DPD_TYPE == 0
            wr  = 1.0;
            #endif
            #if MCMD_DPD_TYPE == 1
            wr  = 1.0 - (r/cutoff_);
            #endif
            #if MCMD_DPD_TYPE == 2
            wr  = 1.0 - rsq/cutoffSq_;
            #endif
           
            // Add random forces to atomic vectors.
            if (computeRandom) {
               fr  = sigma_*wr*randomPtr_->gaussian();
               f.multiply(e, fr);
               randomForces_[atom0Ptr->id()] += f;
               randomForces_[atom1Ptr->id()] -= f;
            }

            // Add dissipative forces to array
            dv.subtract(atom0Ptr->velocity(), atom1Ptr->velocity());
            fd = -gamma_*wr*wr*dv.dot(e);
            f.multiply(e, fd);
            dissipativeForces_[atom0Ptr->id()] += f;
            dissipativeForces_[atom1Ptr->id()] -= f;

         }

      }
      system().positionSignal().notify();
      system().velocitySignal().notify();

   }

   /*
   * DPD integrator time step.
   */
   void NvtDpdVvIntegrator::step() 
   {
      Vector dv;
      Vector dr;
      Vector f;
      System::MoleculeIterator molIter;
      int    iSpecies, nSpecies, atomId;

      nSpecies = simulation().nSpecies();

      // 1st half velocity Verlet, loop over atoms 
      Molecule::AtomIterator atomIter;
      for (iSpecies=0; iSpecies < nSpecies; ++iSpecies) {
         for (system().begin(iSpecies, molIter); molIter.notEnd(); ++molIter)
         {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               atomId = atomIter->id();

               // f = (conservative atomic) + random + dissipative force 
               f.add(atomIter->force(), dissipativeForces_[atomId]);
               f += randomForces_[atomId];

               // Velocity increment dv = force*0.5*dt/mass
               dv.multiply(f, dtMinvFactors_[atomIter->typeId()]);
               atomIter->velocity() += dv;

               dr.multiply(atomIter->velocity(), dt_);
               atomIter->position() += dr;

            }
         }
      }
      system().positionSignal().notify();
      system().velocitySignal().notify();

      #if 0
      #ifndef SIMP_NOPAIR
      if (!system().pairPotential().isPairListCurrent()) {
         system().pairPotential().buildPairList();
      }
      #endif
      #endif

      // Note: MdSystem().calculateForces() rebuilds the pair list as needed.

      // Calculate all new forces
      system().calculateForces();
      computeDpdForces(true);

      // 2nd half velocity Verlet, loop over atoms
      for (iSpecies=0; iSpecies < nSpecies; ++iSpecies) {
         system().begin(iSpecies, molIter); 
         for ( ; molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               atomId = atomIter->id();

               // f = conservative + random + dissipative force 
               f.add(atomIter->force(), dissipativeForces_[atomId]);
               f += randomForces_[atomId];

               // Increment velocity by dv = f*0.5*dt/mass
               dv.multiply(f, dtMinvFactors_[atomIter->typeId()]);
               atomIter->velocity() += dv;

            }
         }
      }
      system().velocitySignal().notify();

      // Recompute disipative force (but not random forces).
      computeDpdForces(false);

   }

}
#undef MCMD_DPD_TYPE
