#ifndef MCMD_MD_SYSTEM_CPP
#define MCMD_MD_SYSTEM_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include "MdSystem.h"
#include "MdSimulation.h"
#include <mcMd/configIos/MdConfigIo.h>
#include <mcMd/mcSimulation/McSystem.h>
#include <util/ensembles/BoundaryEnsemble.h>
#include <util/ensembles/EnergyEnsemble.h>
#include <mcMd/neighbor/PairIterator.h>
#include <mcMd/mdIntegrators/MdIntegratorFactory.h>

#include <mcMd/potentials/pair/MdPairPotential.h>
#include <mcMd/potentials/pair/PairFactory.h>

#include <mcMd/potentials/bond/BondPotential.h>
#include <mcMd/potentials/bond/BondFactory.h>

#ifdef INTER_ANGLE
#include <mcMd/potentials/angle/AnglePotential.h>
#include <mcMd/potentials/angle/AngleFactory.h>
#endif
#ifdef INTER_DIHEDRAL
#include <mcMd/potentials/dihedral/DihedralPotential.h>
#endif
#ifdef MCMD_LINK
#include <mcMd/links/LinkMaster.h>
#endif
#ifdef INTER_EXTERNAL
#include <mcMd/potentials/external/ExternalPotential.h>
#endif
#ifdef INTER_TETHER
#include <mcMd/potentials/tether/TetherPotential.h>
#include <mcMd/tethers/TetherMaster.h>
#endif

#include <util/param/Factory.h>
#include <util/space/Vector.h>
#include <util/space/Dimension.h>

#include <mcMd/simulation/stress.h>
#include <util/space/Tensor.h>
#include <util/accumulators/setToZero.h>

#include <fstream>

namespace McMd
{

   using namespace Util;

   /*
   * Default constructor.
   */
   MdSystem::MdSystem()
    : System(),
      #ifndef INTER_NOPAIR
      pairPotentialPtr_(0),
      #endif
      bondPotentialPtr_(0),
      #ifdef INTER_ANGLE
      anglePotentialPtr_(0),
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralPotentialPtr_(0),
      #endif
      #ifdef MCMD_LINK
      linkPotentialPtr_(0),
      #endif
      #ifdef INTER_EXTERNAL
      externalPotentialPtr_(0),
      #endif
      #ifdef INTER_TETHER
      tetherPotentialPtr_(0),
      #endif
      mdIntegratorPtr_(0),
      mdIntegratorFactoryPtr_(0),
      createdMdIntegratorFactory_(false)
   {  setClassName("MdSystem"); }

   /*
   * Constructor, copy of a System.
   */
   MdSystem::MdSystem(McSystem& system)
    : System(system),
      #ifndef INTER_NOPAIR
      pairPotentialPtr_(0),
      #endif
      bondPotentialPtr_(&system.bondPotential()),
      #ifdef INTER_ANGLE
      anglePotentialPtr_(0),
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralPotentialPtr_(0),
      #endif
      #ifdef MCMD_LINK
      linkPotentialPtr_(0),
      #endif
      #ifdef INTER_EXTERNAL
      externalPotentialPtr_(0),
      #endif
      #ifdef INTER_TETHER
      tetherPotentialPtr_(0),
      #endif
      mdIntegratorPtr_(0),
      mdIntegratorFactoryPtr_(0),
      createdMdIntegratorFactory_(false)
   {
      setClassName("MdSystem");

      #ifndef INTER_NOPAIR
      pairPotentialPtr_ = pairFactory().mdFactory(system.pairPotential());
      if (pairPotentialPtr_ == 0) {
         UTIL_THROW("Failed attempt to clone McPairPotential");
      }
      #endif
      #ifdef INTER_ANGLE
      if (system.hasAnglePotential()) {
         anglePotentialPtr_ = &system.anglePotential();
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (system.hasDihedralPotential()) {
         dihedralPotentialPtr_ = &system.dihedralPotential();
      }
      #endif
      #ifdef MCMD_LINK
      if (system.hasLinkPotential()) {
         linkPotentialPtr_ = &system.linkPotential();
      }
      #endif
      #ifdef INTER_EXTERNAL
      if (system.hasExternalPotential()) {
         externalPotentialPtr_ = &system.externalPotential();
      }
      #endif
      #ifdef INTER_TETHER
      if (system.hasTetherPotential()) {
         tetherPotentialPtr_ = &system.tetherPotential();
      }
      #endif
   }

   /*
   * Destructor.
   */
   MdSystem::~MdSystem()
   {
      #ifndef INTER_NOPAIR
      if (pairPotentialPtr_) delete pairPotentialPtr_;
      #endif
      if (!isCopy() && bondPotentialPtr_) delete bondPotentialPtr_;
      #ifdef INTER_ANGLE
      if (!isCopy() && anglePotentialPtr_) delete anglePotentialPtr_;
      #endif
      #ifdef INTER_DIHEDRAL
      if (!isCopy() && dihedralPotentialPtr_) delete dihedralPotentialPtr_;
      #endif
      #ifdef MCMD_LINK
      if (!isCopy() && linkPotentialPtr_) delete linkPotentialPtr_;
      #endif
      #ifdef INTER_EXTERNAL
      if (!isCopy() && externalPotentialPtr_) delete externalPotentialPtr_;
      #endif
      #ifdef INTER_TETHER
      if (!isCopy() && tetherPotentialPtr_) delete tetherPotentialPtr_;
      #endif

      if (mdIntegratorFactoryPtr_ && createdMdIntegratorFactory_) {
         delete mdIntegratorFactoryPtr_;
      }
      if (mdIntegratorPtr_) {
         delete mdIntegratorPtr_;
      }

   }

   /*
   * Get the MdIntegrator factory
   */
   Factory<MdIntegrator>& MdSystem::mdIntegratorFactory()
   {
      if (mdIntegratorFactoryPtr_ == 0) {
         mdIntegratorFactoryPtr_ = new MdIntegratorFactory(*this);
         createdMdIntegratorFactory_ = true;
      }
      return *mdIntegratorFactoryPtr_;
   }




   /*
   * Read parameter and configuration files, initialize system.
   */
   void MdSystem::readParameters(std::istream &in)
   {

      if (!isCopy()) {
         allocateMoleculeSets();
         readFileMaster(in);
         readPotentialStyles(in);
      }

      #ifndef INTER_NOPAIR
      if (!isCopy()) {
         assert(pairPotentialPtr_ == 0);
         pairPotentialPtr_ = pairFactory().mdFactory(pairStyle(), *this);
         if (pairPotentialPtr_ == 0) {
            UTIL_THROW("Failed attempt to create PairPotential");
         }
      }
      readParamComposite(in, *pairPotentialPtr_);
      #endif

      if (!isCopy()) {

         assert(bondPotentialPtr_ == 0);
         if (simulation().nBondType() > 0) {
            bondPotentialPtr_ = bondFactory().factory(bondStyle());
            if (bondPotentialPtr_ == 0) {
               UTIL_THROW("Failed attempt to create bondPotential");
            }
            readParamComposite(in, *bondPotentialPtr_);
         }

         #ifdef INTER_ANGLE
         assert(anglePotentialPtr_ == 0);
         if (simulation().nAngleType() > 0) {
            anglePotentialPtr_ =
                       angleFactory().factory(angleStyle());
            if (anglePotentialPtr_ == 0) {
               UTIL_THROW("Failed attempt to create anglePotential");
            }
            readParamComposite(in, *anglePotentialPtr_);
         }
         #endif

         #ifdef INTER_DIHEDRAL
         assert(dihedralPotentialPtr_ == 0);
         if (simulation().nDihedralType() > 0) {
            dihedralPotentialPtr_ =
                       dihedralFactory().factory(dihedralStyle());
            if (dihedralPotentialPtr_ == 0) {
               UTIL_THROW("Failed attempt to create dihedralPotential");
            }
            readParamComposite(in, *dihedralPotentialPtr_);
         }
         #endif

         #ifdef MCMD_LINK
         assert(linkPotentialPtr_ == 0);
         if (simulation().nLinkType() > 0) {
            readLinkMaster(in);
            linkPotentialPtr_ = linkFactory().factory(linkStyle());
            if (linkPotentialPtr_ == 0) {
               UTIL_THROW("Failed attempt to create linkPotential");
            }
            readParamComposite(in, *linkPotentialPtr_);
         }
         #endif

         #ifdef INTER_EXTERNAL
         assert(externalPotentialPtr_ == 0);
         if (simulation().hasExternal() > 0) {
            externalPotentialPtr_ =
                       externalFactory().factory(externalStyle());
            if (externalPotentialPtr_ == 0) {
               UTIL_THROW("Failed attempt to create externalPotential");
            }
            readParamComposite(in, *externalPotentialPtr_);
         }
         #endif

         #ifdef INTER_TETHER
         if (simulation().hasExternal() > 0) {
            readTetherMaster(in);
            tetherPotentialPtr_ =
                       tetherFactory().factory(tetherStyle(), *this);
            if (tetherPotentialPtr_ == 0) {
               UTIL_THROW("Failed attempt to create tetherPotential");
            }
            readParamComposite(in, *tetherPotentialPtr_);
         }
         #endif

         // Read EnergyEnsemble and BoundaryEnsemble
         readEnsembles(in);
      }

      // Check for MdIntegratorFactory, create default if necessary.
      if (mdIntegratorFactoryPtr_ == 0) {
         mdIntegratorFactoryPtr_ = new MdIntegratorFactory(*this);
         createdMdIntegratorFactory_ = true;
      }
      assert(mdIntegratorFactoryPtr_);
      assert(mdIntegratorPtr_ == 0);

      // Read polymorphic MdIntegrator
      std::string className;
      bool isEnd;
      mdIntegratorPtr_ =
              mdIntegratorFactoryPtr_->readObject(in, *this, className, isEnd);
      if (!mdIntegratorPtr_) {
         std::string msg("Unknown MdIntegrator subclass name: ");
         msg += className;
         UTIL_THROW(msg.c_str());
      }

      #ifdef MCMD_PERTURB
      // Read Perturbation object for free energy perturbation.
      readPerturbation(in);
      #endif
   }

   /*
   * Load parameter and configuration files, initialize system.
   */
   void MdSystem::loadParameters(Serializable::IArchive& ar)
   {
      if (!isCopy()) {
         allocateMoleculeSets();
         loadFileMaster(ar);
         loadPotentialStyles(ar);
      }

      #ifndef INTER_NOPAIR
      if (!isCopy()) {
         assert(pairPotentialPtr_ == 0);
         pairPotentialPtr_ = pairFactory().mdFactory(pairStyle(), *this);
         if (pairPotentialPtr_ == 0) {
            UTIL_THROW("Failed attempt to create PairPotential");
         }
      }
      loadParamComposite(ar, *pairPotentialPtr_);
      #endif

      if (!isCopy()) {

         assert(bondPotentialPtr_ == 0);
         if (simulation().nBondType() > 0) {
            bondPotentialPtr_ = bondFactory().factory(bondStyle());
            if (bondPotentialPtr_ == 0) {
               UTIL_THROW("Failed attempt to create bondPotential");
            }
            loadParamComposite(ar, *bondPotentialPtr_);
         }

         #ifdef INTER_ANGLE
         assert(anglePotentialPtr_ == 0);
         if (simulation().nAngleType() > 0) {
            anglePotentialPtr_ =
                       angleFactory().factory(angleStyle());
            if (anglePotentialPtr_ == 0) {
               UTIL_THROW("Failed attempt to create anglePotential");
            }
            loadParamComposite(ar, *anglePotentialPtr_);
         }
         #endif

         #ifdef INTER_DIHEDRAL
         assert(dihedralPotentialPtr_ == 0);
         if (simulation().nDihedralType() > 0) {
            dihedralPotentialPtr_ =
                       dihedralFactory().factory(dihedralStyle());
            if (dihedralPotentialPtr_ == 0) {
               UTIL_THROW("Failed attempt to create dihedralPotential");
            }
            loadParamComposite(ar, *dihedralPotentialPtr_);
         }
         #endif

         #ifdef MCMD_LINK
         assert(linkPotentialPtr_ == 0);
         if (simulation().nLinkType() > 0) {
            loadLinkMaster(ar);
            linkPotentialPtr_ = linkFactory().factory(linkStyle());
            if (linkPotentialPtr_ == 0) {
               UTIL_THROW("Failed attempt to create linkPotential");
            }
            loadParamComposite(ar, *linkPotentialPtr_);
         }
         #endif

         #ifdef INTER_EXTERNAL
         assert(externalPotentialPtr_ == 0);
         if (simulation().hasExternal() > 0) {
            externalPotentialPtr_ =
                       externalFactory().factory(externalStyle());
            if (externalPotentialPtr_ == 0) {
               UTIL_THROW("Failed attempt to create externalPotential");
            }
            loadParamComposite(ar, *externalPotentialPtr_);
         }
         #endif

         #ifdef INTER_TETHER
         assert(tetherPotentialPtr_ == 0);
         if (simulation().hasTether() > 0) {
            loadTetherMaster(ar);
            tetherPotentialPtr_ =
                       tetherFactory().factory(tetherStyle(), *this);
            if (tetherPotentialPtr_ == 0) {
               UTIL_THROW("Failed attempt to create tetherPotential");
            }
            loadParamComposite(ar, *tetherPotentialPtr_);
         }
         #endif

         // Read EnergyEnsemble and BoundaryEnsemble
         loadEnsembles(ar);
      }

      // Check for MdIntegratorFactory, create default if necessary.
      if (mdIntegratorFactoryPtr_ == 0) {
         mdIntegratorFactoryPtr_ = new MdIntegratorFactory(*this);
         createdMdIntegratorFactory_ = true;
      }
      assert(mdIntegratorFactoryPtr_);
      assert(mdIntegratorPtr_ == 0);

      // Read polymorphic MdIntegrator
      std::string className;
      mdIntegratorPtr_ =
              mdIntegratorFactoryPtr_->loadObject(ar, *this, className);
      if (!mdIntegratorPtr_) {
         std::string msg("Unknown MdIntegrator subclass name: ");
         msg += className;
         UTIL_THROW(msg.c_str());
      }

      #ifdef MCMD_PERTURB
      // Read Perturbation object for free energy perturbation.
      loadPerturbation(ar);
      #endif
   }

   /*
   * Save state, excluding configuration.
   */
   void MdSystem::saveParameters(Serializable::OArchive& ar)
   {
      if (!isCopy()) {
         saveFileMaster(ar);
         savePotentialStyles(ar);
      }
      #ifndef INTER_NOPAIR
      pairPotentialPtr_->save(ar);
      #endif
      if (!isCopy()) {
         if (simulation().nBondType() > 0) {
            bondPotentialPtr_->save(ar);
         }
         #ifdef INTER_ANGLE
         if (simulation().nAngleType() > 0) {
            anglePotentialPtr_->save(ar);
         }
         #endif
         #ifdef INTER_DIHEDRAL
         if (simulation().nDihedralType() > 0) {
            dihedralPotentialPtr_->save(ar);
         }
         #endif
         #ifdef MCMD_LINK
         if (simulation().nLinkType() > 0) {
            saveLinkMaster(ar);
            linkPotentialPtr_->save(ar);
         }
         #endif
         #ifdef INTER_EXTERNAL
         assert(externalPotentialPtr_ == 0);
         if (simulation().hasExternal() > 0) {
            externalPotentialPtr_->save(ar);
         }
         #endif
         #ifdef INTER_TETHER
         if (simulation().hasTether() > 0) {
            saveTetherMaster(ar);
            tetherPotentialPtr_->save(ar);
         }
         #endif
         saveEnsembles(ar);
      }
      std::string className = mdIntegratorPtr_->className();
      ar & className;
      mdIntegratorPtr_->save(ar);
      #ifdef MCMD_PERTURB
      savePerturbation(ar);
      #endif
   }

   /*
   * Read configuration from a specific input stream.
   */
   void MdSystem::readConfig(std::istream &in)
   {
      System::readConfig(in);

      #ifndef INTER_NOPAIR
      pairPotential().clearPairListStatistics();
      pairPotential().buildPairList();
      #endif
      calculateForces();
   }

   /*
   * Load a System configuration from an archive.
   */
   void MdSystem::loadConfig(Serializable::IArchive& ar)
   {
      System::loadConfig(ar);
      #ifndef INTER_NOPAIR
      pairPotential().clearPairListStatistics();
      pairPotential().buildPairList();
      #endif
      calculateForces();
   }

   /*
   * Shift all atoms into primary unit cell.
   */
   void MdSystem::shiftAtoms()
   {
      MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      for (int iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
         for (begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
            molIter->begin(atomIter);
            for ( ; atomIter.notEnd(); ++atomIter) {
               #ifdef MCMD_SHIFT
               boundary().shift(atomIter->position(), atomIter->shift());
               #else
               boundary().shift(atomIter->position());
               #endif
            }
         }
      }
   }

   /*
   * Set force vectors to zero for all atoms in this System.
   */
   void MdSystem::setZeroForces()
   {
      MoleculeIterator       molIter;
      Molecule::AtomIterator atomIter;
      for (int iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
         for (begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               atomIter->force().zero();
            }
         }
      }
   }

   /*
   * Set velocity vectors to zero for all atoms in this System.
   */
   void MdSystem::setZeroVelocities()
   {
      MoleculeIterator       molIter;
      Molecule::AtomIterator atomIter;
      for (int iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
         for (begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               atomIter->velocity().zero();
            }
         }
      }

   }

   /*
   * Set velocities to Boltzmannian random values, for all atoms in System.
   */
   void MdSystem::setBoltzmannVelocities(double temperature)
   {
      Simulation& sim = simulation();
      Random &random = sim.random();
      double scale;
      double mass;
      MoleculeIterator molIter;
      int nSpec = sim.nSpecies();
      int iSpec, j;

      Molecule::AtomIterator atomIter;
      for (iSpec = 0; iSpec < nSpec; ++iSpec) {
         for (begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               mass   = sim.atomType(atomIter->typeId()).mass();
               scale  = sqrt(temperature/mass);
               for (j = 0; j < Dimension; ++j) {
                  atomIter->velocity()[j] = scale*random.gaussian();
               }
            }
         }
      }

   }

   /*
   * Subtract average velocity from all atomic velocities.
   */
   Vector MdSystem::removeDriftVelocity()
   {
      Vector average(0.0);
      MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      int    nSpec = simulation().nSpecies();
      int    nAtom = 0;
      int    iSpec;

      for (iSpec = 0; iSpec < nSpec; ++iSpec) {
         for (begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               average += atomIter->velocity();
               ++nAtom;
            }
         }
      }
      average /= double(nAtom);
      for (iSpec = 0; iSpec < nSpec; ++iSpec) {
         for (begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               atomIter->velocity() -= average;
            }
         }
      }

      return average;
   }

   /*
   * Calculate all forces.
   */
   void MdSystem::calculateForces()
   {
      setZeroForces();
      #ifndef INTER_NOPAIR
      pairPotential().addForces();
      #endif
      if (hasBondPotential()) {
         bondPotential().addForces();
      }
      #ifdef INTER_ANGLE
      if (hasAnglePotential()) {
         anglePotential().addForces();
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (hasDihedralPotential()) {
         dihedralPotential().addForces();
      }
      #endif
      #ifdef MCMD_LINK
      if (hasLinkPotential()) {
         linkPotential().addForces();
      }
      #endif
      #ifdef INTER_EXTERNAL
      if (hasExternalPotential()) {
         externalPotential().addForces();
      }
      #endif
      #ifdef INTER_TETHER
      if (tetherPotentialPtr_) {
         tetherPotential().addForces();
      }
      #endif
   }

   /*
   * Calculate and return total potential energy.
   */
   double MdSystem::potentialEnergy()
   {
      double energy = 0.0;
      #ifndef INTER_NOPAIR
      energy += pairPotential().energy();
      #endif
      if (hasBondPotential()) {
         energy += bondPotential().energy();
      }
      #ifdef INTER_ANGLE
      if (hasAnglePotential()) {
         energy += anglePotential().energy();
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (hasDihedralPotential()) {
         energy += dihedralPotential().energy();
      }
      #endif
      #ifdef MCMD_LINK
      if (hasLinkPotential()) {
         energy += linkPotential().energy();
      }
      #endif
      #ifdef INTER_EXTERNAL
      if (hasExternalPotential()) {
         energy += externalPotential().energy();
      }
      #endif
      #ifdef INTER_TETHER
      if (tetherPotentialPtr_) {
         energy += tetherPotential().energy();
      }
      #endif
      return energy;
   }

   /*
   * Return total kinetic energy.
   */
   double MdSystem::kineticEnergy() const
   {
      const Simulation& sim = simulation();
      double KE, mass;
      ConstMoleculeIterator molIter;
      Molecule::ConstAtomIterator atomIter;
      int iSpec, typeId;

      KE = 0.0;
      for (iSpec=0; iSpec < sim.nSpecies(); ++iSpec) {
         for (begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               typeId = atomIter->typeId();
               mass = sim.atomType(typeId).mass();
               KE += atomIter->velocity().square()*mass;
            }
         }
      }

      return 0.5*KE;
   }

   // Kinetic Stress

   /*
   * Compute the kinetic (velocity) stress.
   */
   template <typename T>
   void MdSystem::computeKineticStressImpl(T& stress) const
   {
      Vector p;
      double mass;
      ConstMoleculeIterator molIter;
      Molecule::ConstAtomIterator atomIter;
      const Vector* vPtr;
      const Simulation& sim = simulation();

      setToZero(stress);
      for (int iSpec=0; iSpec < sim.nSpecies(); ++iSpec) {
         for (begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               vPtr   = &atomIter->velocity();
               mass   = sim.atomType(atomIter->typeId()).mass();
               p.multiply(*vPtr, mass);
               incrementPairStress(*vPtr, p, stress);
            }
         }
      }
      stress /= boundary().volume();
      normalizeStress(stress);
   }

   template <>
   void MdSystem::computeKineticStress<double>(double& stress) const
   {  computeKineticStressImpl(stress); }

   template <>
   void MdSystem::computeKineticStress<Util::Vector>(Util::Vector& stress) const
   {  computeKineticStressImpl(stress); }

   template <>
   void MdSystem::computeKineticStress<Util::Tensor>(Util::Tensor& stress) const
   {  computeKineticStressImpl(stress); }

   // Virial Stress

   template <typename T>
   void MdSystem::computeVirialStressImpl(T& stress) const
   {
      setToZero(stress);

      #ifndef INTER_NOPAIR
      T pairStress;
      pairPotential().computeStress(pairStress);
      stress += pairStress;
      #endif

      if (hasBondPotential()) {
         T bondStress;
         bondPotential().computeStress(bondStress);
         stress += bondStress;
      }

      #ifdef INTER_ANGLE
      if (hasAnglePotential()) {
         T angleStress;
         anglePotential().computeStress(angleStress);
         stress += angleStress;
      }
      #endif

      #ifdef INTER_DIHEDRAL
      if (hasDihedralPotential()) {
         T dihedralStress;
         dihedralPotential().computeStress(dihedralStress);
         stress += dihedralStress;
      }
      #endif

      #ifdef MCMD_LINK
      if (hasLinkPotential()) {
         T linkStress;
         linkPotential().computeStress(linkStress);
         stress += linkStress;
      }
      #endif

   }

   template <>
   void MdSystem::computeVirialStress<double>(double& stress) const
   {  computeVirialStressImpl(stress); }

   template <>
   void MdSystem::computeVirialStress<Util::Vector>(Util::Vector& stress)
      const
   {  computeVirialStressImpl(stress); }

   template <>
   void MdSystem::computeVirialStress<Util::Tensor>(Util::Tensor& stress)
      const
   {  computeVirialStressImpl(stress); }

   // Total Stress

   template <>
   void MdSystem::computeStress<double>(double& stress) const
   {
      computeVirialStress(stress);

      double kineticStress;
      computeKineticStress(kineticStress);
      stress += kineticStress;
   }

   template <>
   void MdSystem::computeStress<Util::Vector>(Util::Vector& stress) const
   {
      Util::Vector kineticStress;
      computeKineticStress(kineticStress);

      computeVirialStress(stress);
      stress += kineticStress;

   }

   template <>
   void MdSystem::computeStress<Util::Tensor>(Util::Tensor& stress) const
   {
      computeVirialStress(stress);

      Util::Tensor kineticStress;
      computeKineticStress(kineticStress);
      stress += kineticStress;
   }

   /*
   * Return true if this MdSystem is valid, or an throw Exception.
   */
   bool MdSystem::isValid() const
   {
      System::isValid();
      #ifndef INTER_NOPAIR
      pairPotentialPtr_->pairList().isValid();
      #endif
      return true;
   }

   /*
   * Return a pointer to a new default ConfigIo.
   */
   ConfigIo* MdSystem::newDefaultConfigIo()
   {  return new MdConfigIo(*this); }

}
#endif
