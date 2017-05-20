/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
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
#include <mcMd/neighbor/CellList.h>
#include <mcMd/mdIntegrators/MdIntegratorFactory.h>
#include <mcMd/generators/Generator.h>
#include <mcMd/generators/generatorFactory.h>

#include <mcMd/potentials/pair/MdPairPotential.h>
#include <mcMd/potentials/pair/PairFactory.h>

#ifdef SIMP_BOND
#include <mcMd/potentials/bond/BondPotential.h>
#include <mcMd/potentials/bond/BondFactory.h>
#endif
#ifdef SIMP_ANGLE
#include <mcMd/potentials/angle/AnglePotential.h>
#include <mcMd/potentials/angle/AngleFactory.h>
#endif
#ifdef SIMP_DIHEDRAL
#include <mcMd/potentials/dihedral/DihedralPotential.h>
#endif
#ifdef SIMP_COULOMB
#include <mcMd/potentials/coulomb/MdCoulombPotential.h>
#include <mcMd/potentials/coulomb/CoulombFactory.h>
#endif
#ifdef SIMP_EXTERNAL
#include <mcMd/potentials/external/ExternalPotential.h>
#endif
#ifdef MCMD_LINK
#include <mcMd/links/LinkMaster.h>
#endif
#ifdef SIMP_TETHER
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
      #ifndef SIMP_NOPAIR
      pairPotentialPtr_(0),
      #endif
      #ifdef SIMP_BOND
      bondPotentialPtr_(0),
      #endif
      #ifdef SIMP_ANGLE
      anglePotentialPtr_(0),
      #endif
      #ifdef SIMP_DIHEDRAL
      dihedralPotentialPtr_(0),
      #endif
      #ifdef SIMP_COULOMB
      coulombPotentialPtr_(0),
      #endif
      #ifdef SIMP_EXTERNAL
      externalPotentialPtr_(0),
      #endif
      #ifdef MCMD_LINK
      linkPotentialPtr_(0),
      #endif
      #ifdef SIMP_TETHER
      tetherPotentialPtr_(0),
      #endif
      mdIntegratorPtr_(0),
      mdIntegratorFactoryPtr_(0),
      createdMdIntegratorFactory_(false)
   {  
      setClassName("MdSystem"); 

      // Set actions taken when particles are moved
      positionSignal().addObserver(*this, &MdSystem::unsetPotentialEnergy);
      positionSignal().addObserver(*this, &MdSystem::unsetVirialStress);
   }

   /*
   * Constructor, copy of a System.
   */
   MdSystem::MdSystem(McSystem& system)
    : System(system),
      #ifndef SIMP_NOPAIR
      pairPotentialPtr_(0),
      #endif
      #ifdef SIMP_BOND
      bondPotentialPtr_(0),
      #endif
      #ifdef SIMP_ANGLE
      anglePotentialPtr_(0),
      #endif
      #ifdef SIMP_DIHEDRAL
      dihedralPotentialPtr_(0),
      #endif
      #ifdef SIMP_COULOMB
      coulombPotentialPtr_(0),
      #endif
      #ifdef SIMP_EXTERNAL
      externalPotentialPtr_(0),
      #endif
      #ifdef MCMD_LINK
      linkPotentialPtr_(0),
      #endif
      #ifdef SIMP_TETHER
      tetherPotentialPtr_(0),
      #endif
      mdIntegratorPtr_(0),
      mdIntegratorFactoryPtr_(0),
      createdMdIntegratorFactory_(false)
   {
      setClassName("MdSystem");

      #ifndef SIMP_NOPAIR
      assert(pairPotentialPtr_ == 0);
      pairPotentialPtr_ = pairFactory().mdFactory(system.pairPotential());
      if (pairPotentialPtr_ == 0) {
         UTIL_THROW("Failed attempt to clone McPairPotential");
      }
      #endif
      #ifdef SIMP_BOND
      assert(bondPotentialPtr_ == 0);
      if (system.hasBondPotential()) {
         bondPotentialPtr_ = &system.bondPotential();
      }
      #endif
      #ifdef SIMP_ANGLE
      assert(anglePotentialPtr_ == 0);
      if (system.hasAnglePotential()) {
         anglePotentialPtr_ = &system.anglePotential();
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      assert(dihedralPotentialPtr_ == 0);
      if (system.hasDihedralPotential()) {
         dihedralPotentialPtr_ = &system.dihedralPotential();
      }
      #endif
      #ifdef SIMP_COULOMB
      #if 0
      assert(coulombPotentialPtr_ == 0);
      if (system.hasCoulombPotential()) {
         coulombPotentialPtr_ = &system.coulombPotential();
      }
      #endif
      #endif
      #ifdef SIMP_EXTERNAL
      if (system.hasExternalPotential()) {
         externalPotentialPtr_ = &system.externalPotential();
      }
      #endif
      #ifdef MCMD_LINK
      if (system.hasLinkPotential()) {
         linkPotentialPtr_ = &system.linkPotential();
      }
      #endif
      #ifdef SIMP_TETHER
      if (system.hasTetherPotential()) {
         tetherPotentialPtr_ = &system.tetherPotential();
      }
      #endif

      // Set actions taken when particles are moved
      positionSignal().addObserver(*this, &MdSystem::unsetPotentialEnergy);
      positionSignal().addObserver(*this, &MdSystem::unsetVirialStress);
   }

   /*
   * Destructor.
   */
   MdSystem::~MdSystem()
   {
      #ifndef SIMP_NOPAIR
      if (pairPotentialPtr_) delete pairPotentialPtr_;
      #endif
      #ifdef SIMP_BOND
      if (!isCopy() && bondPotentialPtr_) delete bondPotentialPtr_;
      #endif
      #ifdef SIMP_ANGLE
      if (!isCopy() && anglePotentialPtr_) delete anglePotentialPtr_;
      #endif
      #ifdef SIMP_DIHEDRAL
      if (!isCopy() && dihedralPotentialPtr_) delete dihedralPotentialPtr_;
      #endif
      #ifdef SIMP_COULOMB
      if (!isCopy() && coulombPotentialPtr_) delete coulombPotentialPtr_;
      #endif
      #ifdef MCMD_LINK
      if (!isCopy() && linkPotentialPtr_) delete linkPotentialPtr_;
      #endif
      #ifdef SIMP_EXTERNAL
      if (!isCopy() && externalPotentialPtr_) delete externalPotentialPtr_;
      #endif
      #ifdef SIMP_TETHER
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

      #ifdef SIMP_COULOMB
      if (!isCopy()) {
         assert(coulombPotentialPtr_ == 0);
         if (simulation().hasCoulomb()) {
            coulombPotentialPtr_ =
                       coulombFactory().factory(coulombStyle());
            if (coulombPotentialPtr_ == 0) {
               UTIL_THROW("Failed attempt to create CoulombPotential");
            }
            readParamComposite(in, *coulombPotentialPtr_);
         }
      }
      #endif

      #ifndef SIMP_NOPAIR
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

         #ifdef SIMP_BOND
         assert(bondPotentialPtr_ == 0);
         if (simulation().nBondType() > 0) {
            bondPotentialPtr_ = bondFactory().factory(bondStyle());
            if (bondPotentialPtr_ == 0) {
               UTIL_THROW("Failed attempt to create bondPotential");
            }
            readParamComposite(in, *bondPotentialPtr_);
         }
         #endif

         #ifdef SIMP_ANGLE
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

         #ifdef SIMP_DIHEDRAL
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

         #ifdef SIMP_EXTERNAL
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

         #ifdef SIMP_TETHER
         if (simulation().hasTether() > 0) {
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

      #ifdef SIMP_COULOMB
      if (!isCopy()) {
         assert(coulombPotentialPtr_ == 0);
         if (simulation().hasCoulomb() > 0) {
            coulombPotentialPtr_ =
                       coulombFactory().factory(coulombStyle());
            if (coulombPotentialPtr_ == 0) {
               UTIL_THROW("Failed attempt to create CoulombPotential");
            }
            loadParamComposite(ar, *coulombPotentialPtr_);
         }
      }
      #endif

      #ifndef SIMP_NOPAIR
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

         #ifdef SIMP_BOND
         assert(bondPotentialPtr_ == 0);
         if (simulation().nBondType() > 0) {
            bondPotentialPtr_ = bondFactory().factory(bondStyle());
            if (bondPotentialPtr_ == 0) {
               UTIL_THROW("Failed attempt to create bondPotential");
            }
            loadParamComposite(ar, *bondPotentialPtr_);
         }
         #endif

         #ifdef SIMP_ANGLE
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

         #ifdef SIMP_DIHEDRAL
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

         #ifdef SIMP_EXTERNAL
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

         #ifdef SIMP_TETHER
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

      #ifdef SIMP_COULOMB
      if (!isCopy()) {
         if (simulation().hasCoulomb()) {
            coulombPotentialPtr_->save(ar);
         }
      }
      #endif

      #ifndef SIMP_NOPAIR
      pairPotential().save(ar);
      #endif

      if (!isCopy()) {
         #ifdef SIMP_BOND
         if (simulation().nBondType() > 0) {
            bondPotential().save(ar);
         }
         #endif
         #ifdef SIMP_ANGLE
         if (simulation().nAngleType() > 0) {
            anglePotential().save(ar);
         }
         #endif
         #ifdef SIMP_DIHEDRAL
         if (simulation().nDihedralType() > 0) {
            dihedralPotential().save(ar);
         }
         #endif
         #ifdef MCMD_LINK
         if (simulation().nLinkType() > 0) {
            saveLinkMaster(ar);
            linkPotential().save(ar);
         }
         #endif
         #ifdef SIMP_EXTERNAL
         if (simulation().hasExternal()) {
            externalPotential().save(ar);
         }
         #endif
         #ifdef SIMP_TETHER
         if (simulation().hasTether()) {
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

      #ifndef SIMP_NOPAIR
      pairPotential().clearPairListStatistics();
      pairPotential().buildPairList();
      #endif
      #ifdef SIMP_COULOMB
      if (hasCoulombPotential()) {
         coulombPotential().makeWaves();
         Log::file() << "Initial coulombPotential nWave = " 
                     << coulombPotential().nWave() << std::endl;
      }
      #endif
      calculateForces();
   }

   /*
   * Load a System configuration from an archive.
   */
   void MdSystem::loadConfig(Serializable::IArchive& ar)
   {
      System::loadConfig(ar);
      #ifndef SIMP_NOPAIR
      pairPotential().clearPairListStatistics();
      pairPotential().buildPairList();
      #endif
      #ifdef SIMP_COULOMB
      if (hasCoulombPotential()) {
         coulombPotential().makeWaves();
      }
      #endif
      calculateForces();
   }


   /*
   * Generate molecules for all species.
   */
   void MdSystem::generateMolecules(
                       Array<int> const & capacities,
                       Array<double> const & diameters)
   {

      // Setup a local cell list
      CellList cellList;
      Generator::setupCellList(simulation().atomCapacity(), 
                               boundary(), diameters, cellList);

      // Generate molecules for each species
      Generator* ptr = 0;
      int nSpecies = simulation().nSpecies();
      bool success = false;
      for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         if (capacities[iSpecies] > 0) {
            ptr = generatorFactory(simulation().species(iSpecies), *this);
            UTIL_CHECK(ptr);
            success = ptr->generate(capacities[iSpecies], 
                                    diameters, cellList);
            delete ptr;
            if (!success) {
               Log::file() << "Failed to complete species " 
                           << iSpecies << "\n";
            }
            UTIL_CHECK(success);
         }
      }

      #ifndef SIMP_NOPAIR
      pairPotential().clearPairListStatistics();
      pairPotential().buildPairList();
      #endif
      #ifdef UTIL_DEBUG
      isValid();
      #endif

      #ifdef SIMP_COULOMB
      if (hasCoulombPotential()) {
         coulombPotential().makeWaves();
      }
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
      #ifndef SIMP_NOPAIR
      // This method builds pair list if needed, and shifts atoms if
      // it builds the pair list.
      pairPotential().addForces();
      #endif
      #ifdef SIMP_BOND
      if (hasBondPotential()) {
         bondPotential().addForces();
      }
      #endif
      #ifdef SIMP_ANGLE
      if (hasAnglePotential()) {
         anglePotential().addForces();
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (hasDihedralPotential()) {
         dihedralPotential().addForces();
      }
      #endif
      #ifdef SIMP_COULOMB
      if (hasCoulombPotential()) {
         coulombPotential().addForces();
      }
      #endif
      #ifdef SIMP_EXTERNAL
      if (hasExternalPotential()) {
         externalPotential().addForces();
      }
      #endif
      #ifdef MCMD_LINK
      if (hasLinkPotential()) {
         linkPotential().addForces();
      }
      #endif
      #ifdef SIMP_TETHER
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
      #ifndef SIMP_NOPAIR
      energy += pairPotential().energy();
      #endif
      #ifdef SIMP_BOND
      if (hasBondPotential()) {
         energy += bondPotential().energy();
      }
      #endif
      #ifdef SIMP_ANGLE
      if (hasAnglePotential()) {
         energy += anglePotential().energy();
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (hasDihedralPotential()) {
         energy += dihedralPotential().energy();
      }
      #endif
      #ifdef SIMP_COULOMB
      if (hasCoulombPotential()) {
         energy += coulombPotential().kSpaceEnergy();
      }
      #endif
      #ifdef SIMP_EXTERNAL
      if (hasExternalPotential()) {
         energy += externalPotential().energy();
      }
      #endif
      #ifdef MCMD_LINK
      if (hasLinkPotential()) {
         energy += linkPotential().energy();
      }
      #endif
      #ifdef SIMP_TETHER
      if (tetherPotentialPtr_) {
         energy += tetherPotential().energy();
      }
      #endif
      return energy;
   }

   /*
   * Unset precomputed potential energy components.
   */
   void MdSystem::unsetPotentialEnergy()
   {
      #ifndef SIMP_NOPAIR
      pairPotential().unsetEnergy();
      #endif
      #ifdef SIMP_BOND
      if (hasBondPotential()) {
         bondPotential().unsetEnergy();
      }
      #endif
      #ifdef SIMP_ANGLE
      if (hasAnglePotential()) {
         anglePotential().unsetEnergy();
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (hasDihedralPotential()) {
         dihedralPotential().unsetEnergy();
      }
      #endif
      #ifdef SIMP_COULOMB
      if (hasCoulombPotential()) {
         coulombPotential().unsetEnergy();
      }
      #endif
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
      T dStress;

      #ifndef SIMP_NOPAIR
      pairPotential().computeStress(dStress);
      stress += dStress;
      #endif
      #ifdef SIMP_BOND
      if (hasBondPotential()) {
         bondPotential().computeStress(dStress);
         stress += dStress;
      }
      #endif
      #ifdef SIMP_ANGLE
      if (hasAnglePotential()) {
         anglePotential().computeStress(dStress);
         stress += dStress;
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (hasDihedralPotential()) {
         dihedralPotential().computeStress(dStress);
         stress += dStress;
      }
      #endif
      #ifdef SIMP_COULOMB
      if (hasCoulombPotential()) {
         coulombPotential().computeStress(dStress);
         stress += dStress;
      }
      #endif
      #ifdef MCMD_LINK
      if (hasLinkPotential()) {
         linkPotential().computeStress(dStress);
         stress += dStress;
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
   * Unset precomputed virial stress components.
   */
   void MdSystem::unsetVirialStress()
   {
      #ifndef SIMP_NOPAIR
      pairPotential().unsetStress();
      #endif
      #ifdef SIMP_BOND
      if (hasBondPotential()) {
         bondPotential().unsetStress();
      }
      #endif
      #ifdef SIMP_ANGLE
      if (hasAnglePotential()) {
         anglePotential().unsetStress();
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (hasDihedralPotential()) {
         dihedralPotential().unsetStress();
      }
      #endif
      #ifdef SIMP_COULOMB
      if (hasCoulombPotential()) {
         coulombPotential().unsetStress();
      }
      #endif
   }

   // Miscellaneous member functions

   /*
   * Return true if this MdSystem is valid, or an throw Exception.
   */
   bool MdSystem::isValid() const
   {
      System::isValid();
      #ifndef SIMP_NOPAIR
      pairPotential().pairList().isValid();
      #endif
      return true;
   }

   /*
   * Return a pointer to a new default ConfigIo.
   */
   ConfigIo* MdSystem::newDefaultConfigIo()
   {  return new MdConfigIo(*this); }

}
