/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McSystem.h"
#include "McSimulation.h"
#include <mcMd/simulation/stress.h>
#include <util/ensembles/BoundaryEnsemble.h>
#include <util/ensembles/EnergyEnsemble.h>
#include <mcMd/generators/Generator.h>
#include <mcMd/generators/generatorFactory.h>

#ifndef SIMP_NOPAIR
#include <mcMd/potentials/pair/McPairPotential.h>
#include <mcMd/potentials/pair/PairFactory.h>
#endif
#ifdef SIMP_BOND
#include <mcMd/potentials/bond/BondPotential.h>
#endif
#ifdef SIMP_ANGLE
#include <mcMd/potentials/angle/AnglePotential.h>
#endif
#ifdef SIMP_DIHEDRAL
#include <mcMd/potentials/dihedral/DihedralPotential.h>
#endif
#ifdef MCMD_LINK
#include <mcMd/links/LinkMaster.h>
#endif
#ifdef SIMP_EXTERNAL
#include <mcMd/potentials/external/ExternalPotential.h>
#include <mcMd/potentials/external/ExternalPotential.h>
#endif
#ifdef SIMP_TETHER
#include <mcMd/tethers/TetherMaster.h>
#include <mcMd/potentials/tether/TetherPotential.h>
#endif

#ifdef MCMD_PERTURB
#include <mcMd/perturb/mcSystem/McPerturbationFactory.h>
#endif

#include <simp/species/Species.h>

#include <util/param/Factory.h>
#include <util/space/Vector.h>
#include <util/space/Tensor.h>
#include <util/space/Dimension.h>
#include <util/archives/Serializable_includes.h>
#include <util/accumulators/setToZero.h>

#include <util/global.h>

#include <fstream>

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /*
   * Constructor.
   */
   McSystem::McSystem()
    :
      System()
      #ifndef SIMP_NOPAIR
      , pairPotentialPtr_(0)
      #endif
      #ifdef SIMP_BOND
      , bondPotentialPtr_(0)
      #endif
      #ifdef SIMP_ANGLE
      , anglePotentialPtr_(0)
      #endif
      #ifdef SIMP_DIHEDRAL
      , dihedralPotentialPtr_(0)
      #endif
      #ifdef MCMD_LINK
      , linkPotentialPtr_(0)
      #endif
      #ifdef SIMP_EXTERNAL
      , externalPotentialPtr_(0)
      #endif
      #ifdef SIMP_TETHER
      , tetherPotentialPtr_(0)
      #endif
   { 
      setClassName("McSystem"); 

      // Actions taken when particles are moved
      positionSignal().addObserver(*this, &McSystem::unsetPotentialEnergies);
      positionSignal().addObserver(*this, &McSystem::unsetVirialStress);
   }


   /*
   * Destructor.
   */
   McSystem::~McSystem()
   {
      #ifndef SIMP_NOPAIR
      if (pairPotentialPtr_) delete pairPotentialPtr_;
      #endif
      #ifdef SIMP_BOND
      if (bondPotentialPtr_) delete bondPotentialPtr_;
      #endif
      #ifdef SIMP_ANGLE
      if (anglePotentialPtr_) delete anglePotentialPtr_;
      #endif
      #ifdef SIMP_DIHEDRAL
      if (dihedralPotentialPtr_) delete dihedralPotentialPtr_;
      #endif
      #ifdef MCMD_LINK
      if (linkPotentialPtr_) delete linkPotentialPtr_;
      #endif
      #ifdef SIMP_EXTERNAL
      if (externalPotentialPtr_) delete externalPotentialPtr_;
      #endif
      #ifdef SIMP_TETHER
      if (tetherPotentialPtr_) delete tetherPotentialPtr_;
      #endif
   }

   // -------------------------------------------------------------
   // Initialization and Parameter I/O 
   
   /*
   * Read parameters from file.
   */
   void McSystem::readParameters(std::istream &in)
   {
      allocateMoleculeSets();
      readFileMaster(in);
      readPotentialStyles(in);

      #ifndef SIMP_NOPAIR
      assert(pairPotentialPtr_ == 0);
      pairPotentialPtr_ = pairFactory().mcFactory(pairStyle(), *this);
      if (pairPotentialPtr_ == 0) {
         UTIL_THROW("Failed attempt to create McPairPotential");
      }
      readParamComposite(in, *pairPotentialPtr_);
      #endif

      #ifdef SIMP_BOND
      assert(bondPotentialPtr_ == 0);
      if (simulation().nBondType() > 0) {
         bondPotentialPtr_ = bondFactory().factory(bondStyle());
         if (bondPotentialPtr_ == 0) {
            UTIL_THROW("Failed attempt to create BondPotential");
         }
         readParamComposite(in, *bondPotentialPtr_);
      }
      #endif

      #ifdef SIMP_ANGLE
      assert(anglePotentialPtr_ == 0);
      if (simulation().nAngleType() > 0) {
         anglePotentialPtr_ = angleFactory().factory(angleStyle());
         if (anglePotentialPtr_ == 0) {
            UTIL_THROW("Failed attempt to create AnglePotential");
         }
         readParamComposite(in, *anglePotentialPtr_);
      }
      #endif

      #ifdef SIMP_DIHEDRAL
      assert(dihedralPotentialPtr_ == 0);
      if (simulation().nDihedralType() > 0) {
         dihedralPotentialPtr_ = dihedralFactory().factory(dihedralStyle());
         if (dihedralPotentialPtr_ == 0) {
            UTIL_THROW("Failed attempt to create DihedralPotential");
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
            UTIL_THROW("Failed attempt to create BondPotential for links");
         }
         readParamComposite(in, *linkPotentialPtr_);
      }
      #endif

      #ifdef SIMP_EXTERNAL
      assert(externalPotentialPtr_ == 0);
      if (simulation().hasExternal()) {
         externalPotentialPtr_ =
            externalFactory().factory(externalStyle());
         if (externalPotentialPtr_ == 0) {
            UTIL_THROW("Failed attempt to create ExternalPotential");
         }
         readParamComposite(in, *externalPotentialPtr_);
      }
      #endif

      #ifdef SIMP_TETHER
      if (simulation().hasTether()) {
         readTetherMaster(in);
         tetherPotentialPtr_ = tetherFactory().factory(tetherStyle(), *this);
         if (tetherPotentialPtr_ == 0) {
            UTIL_THROW("Failed attempt to create TetherPotential");
         }
         readParamComposite(in, *tetherPotentialPtr_);
      }
      #endif

      // Read EnergyEnsemble and BoundaryEnsemble
      readEnsembles(in);

      #ifdef MCMD_PERTURB
      readPerturbation(in);
      #ifdef UTIL_MPI
      readReplicaMove(in);
      #endif
      #endif
   }

   /*
   * Load parameters from an archive.
   */
   void McSystem::loadParameters(Serializable::IArchive& ar)
   {
      allocateMoleculeSets();
      loadFileMaster(ar);
      loadPotentialStyles(ar);

      #ifndef SIMP_NOPAIR
      pairPotentialPtr_ = pairFactory().mcFactory(pairStyle(), *this);
      if (pairPotentialPtr_ == 0) {
         UTIL_THROW("Failed attempt to create McPairPotential");
      }
      loadParamComposite(ar, *pairPotentialPtr_);
      #endif

      #ifdef SIMP_BOND
      assert(bondPotentialPtr_ == 0);
      if (simulation().nBondType() > 0) {
         bondPotentialPtr_ = bondFactory().factory(bondStyle());
         if (bondPotentialPtr_ == 0) {
            UTIL_THROW("Failed attempt to create BondPotential");
         }
         loadParamComposite(ar, *bondPotentialPtr_);
      }
      #endif

      #ifdef SIMP_ANGLE
      assert(anglePotentialPtr_ == 0);
      if (simulation().nAngleType() > 0) {
         anglePotentialPtr_ = angleFactory().factory(angleStyle());
         if (anglePotentialPtr_ == 0) {
            UTIL_THROW("Failed attempt to create AnglePotential");
         }
         loadParamComposite(ar, *anglePotentialPtr_);
      }
      #endif

      #ifdef SIMP_DIHEDRAL
      assert(dihedralPotentialPtr_ == 0);
      if (simulation().nDihedralType() > 0) {
         dihedralPotentialPtr_ = dihedralFactory().factory(dihedralStyle());
         if (dihedralPotentialPtr_ == 0) {
            UTIL_THROW("Failed attempt to create DihedralPotential");
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
            UTIL_THROW("Failed attempt to create BondPotential for links");
         }
         loadParamComposite(ar, *linkPotentialPtr_);
      }
      #endif

      #ifdef SIMP_EXTERNAL
      assert(externalPotentialPtr_ == 0);
      if (simulation().hasExternal()) {
         externalPotentialPtr_ =
            externalFactory().factory(externalStyle());
         if (externalPotentialPtr_ == 0) {
            UTIL_THROW("Failed attempt to create ExternalPotential");
         }
         loadParamComposite(ar, *externalPotentialPtr_);
      }
      #endif

      #ifdef SIMP_TETHER
      if (simulation().hasTether()) {
         loadTetherMaster(ar);
         tetherPotentialPtr_ = tetherFactory().factory(tetherStyle(), *this);
         if (tetherPotentialPtr_ == 0) {
            UTIL_THROW("Failed attempt to create TetherPotential");
         }
         loadParamComposite(ar, *tetherPotentialPtr_);
      }
      #endif

      loadEnsembles(ar);

      #ifdef MCMD_PERTURB
      loadPerturbation(ar);
      #ifdef UTIL_MPI
      loadReplicaMove(ar);
      #endif
      #endif

   }

   /* 
   * Save parameters.
   */
   void McSystem::saveParameters(Serializable::OArchive& ar) 
   {
      saveFileMaster(ar);
      savePotentialStyles(ar);
      #ifndef SIMP_NOPAIR 
      pairPotential().save(ar); 
      #endif
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
         assert(linkPotentialPtr_);
         linkPotentialPtr_->save(ar); 
      }
      #endif
      #ifdef SIMP_EXTERNAL
      assert(externalPotentialPtr_ == 0);
      if (simulation().hasExternal() > 0) {
         assert(externalPotentialPtr_);
         externalPotentialPtr_->save(ar); 
      }
      #endif
      #ifdef SIMP_TETHER
      if (simulation().hasExternal() > 0) {
         saveTetherMaster(ar);
         assert(tetherPotentialPtr_);
         tetherPotentialPtr_->save(ar); 
      }
      #endif

      saveEnsembles(ar);

      #ifdef MCMD_PERTURB
      savePerturbation(ar);
      #ifdef UTIL_MPI
      saveReplicaMove(ar);
      #endif
      #endif
   }
  
   // -------------------------------------------------------------
   // Configuration I/o
   
   /*
   * Read configuration from a specific input stream.
   */
   void McSystem::readConfig(std::istream &in)
   {
      System::readConfig(in);
      #ifndef SIMP_NOPAIR
      pairPotential().buildCellList();
      #endif
   }

   /* 
   * Load a System configuration from an archive.
   */
   void McSystem::loadConfig(Serializable::IArchive& ar)
   {  
      System::loadConfig(ar); 
      #ifndef SIMP_NOPAIR
      pairPotential().buildCellList();
      #endif
   }

   /*
   * Generate molecules for all species.
   */
   void McSystem::generateMolecules(
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
      pairPotential().buildCellList();
      #endif

      #ifdef UTIL_DEBUG
      isValid();
      #endif
   }

   // -------------------------------------------------------------
   // Energy Evaluators (including all components)

   /*
   * Return total potential energy for one Atom.
   */
   double McSystem::atomPotentialEnergy(const Atom &atom) const
   {
      double energy = 0.0;
      #ifndef SIMP_NOPAIR
      energy += pairPotential().atomEnergy(atom);
      #endif
      #ifdef SIMP_BOND
      if (hasBondPotential()) {
         energy += bondPotential().atomEnergy(atom);
      }
      #endif
      #ifdef SIMP_ANGLE
      if (hasAnglePotential()) {
         energy += anglePotential().atomEnergy(atom);
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (hasDihedralPotential()) {
         energy += dihedralPotential().atomEnergy(atom);
      }
      #endif
      #ifdef MCMD_LINK
      if (hasLinkPotential()) {
         energy += linkPotential().atomEnergy(atom);
      }
      #endif
      #ifdef SIMP_EXTERNAL
      if (hasExternalPotential()) {
         energy += externalPotential().atomEnergy(atom);
      }
      #endif
      #ifdef SIMP_TETHER
      if (tetherPotentialPtr_) {
         if (tetherMaster().isTethered(atom)) {
            energy += tetherPotential().atomEnergy(atom);
         }
      }
      #endif
      return energy;
   }

   /*
   * Return total potential energy for this McSystem.
   */
   double McSystem::potentialEnergy() const
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
      #ifdef MCMD_LINK
      if (hasLinkPotential()) {
         energy += linkPotential().energy();
      }
      #endif
      #ifdef SIMP_EXTERNAL
      if (hasExternalPotential()) {
         energy += externalPotential().energy();
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
   * Unset all precomputed potential energy components.
   */
   void McSystem::unsetPotentialEnergies()
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
   }

   // -------------------------------------------------------------
   // Pressure/Stress Evaluators (including all components)

   template <typename T>
   void McSystem::computeVirialStressImpl(T& stress) const
   {
      setToZero(stress);

      #ifndef SIMP_NOPAIR
      T pairStress;
      pairPotential().computeStress(pairStress);
      stress += pairStress;
      #endif
      #ifdef SIMP_BOND
      if (hasBondPotential()) {
         T bondStress;
         bondPotential().computeStress(bondStress);
         stress += bondStress;
      }
      #endif
      #ifdef SIMP_ANGLE
      if (hasAnglePotential()) {
         T angleStress;
         anglePotential().computeStress(angleStress);
         stress += angleStress;
      }
      #endif
      #ifdef SIMP_DIHEDRAL
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

   /*
   * Compute virial pressure (isotropic)
   */
   template <>
   void McSystem::computeVirialStress<double>(double& stress) const
   {  computeVirialStressImpl(stress); }

   /*
   * Compute virial stress (diagonal components)
   */
   template <>
   void McSystem::computeVirialStress<Util::Vector>(Util::Vector& stress) const
   {  computeVirialStressImpl(stress); }

   /*
   * Compute virial stress (full tensor).
   */
   template <>
   void McSystem::computeVirialStress<Util::Tensor>(Util::Tensor& stress) const
   {  computeVirialStressImpl(stress); }

   /*
   * Compute total pressure (isotropic)
   */
   template <>
   void McSystem::computeStress<double>(double& stress) const
   {
      computeVirialStress(stress);
      if (energyEnsemble().isIsothermal()) {
         double rho = nAtom()/boundary().volume();
         double temperature = energyEnsemble().temperature();
         stress += rho*temperature;
      }
   }

   /*
   * Compute total stress (diagonal components)
   */
   template <>
   void McSystem::computeStress<Util::Vector>(Util::Vector& stress) const
   {
      computeVirialStress(stress);
      if (energyEnsemble().isIsothermal()) {
         double rho = nAtom()/boundary().volume();
         double temperature = energyEnsemble().temperature();
         for (int i = 0; i < Dimension; ++i) {
            stress[i] += rho*temperature;
         }
      }
   }

   /*
   * Compute total stress (full tensor)
   */
   template <>
   void McSystem::computeStress<Util::Tensor>(Util::Tensor& stress) const
   {
      computeVirialStress(stress);
      if (energyEnsemble().isIsothermal()) {
         double rho = nAtom()/boundary().volume();
         double temperature = energyEnsemble().temperature();
         for (int i = 0; i < Dimension; ++i) {
            stress(i, i) += rho*temperature;
         }
      }
   }

   /*
   * Unset all precomputed virial stress components.
   */
   void McSystem::unsetVirialStress()
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
   }

   // -------------------------------------------------------------
   // Miscellaneous
   
   #ifdef MCMD_PERTURB
   /*
   * Return a pointer to a new default McPerturbationFactory.
   */
   Factory<Perturbation>* McSystem::newDefaultPerturbationFactory()
   {  return new McPerturbationFactory(*this); }
   #endif

   /*
   * Return true if this McSystem is valid, or an throw Exception.
   */
   bool McSystem::isValid() const
   {
      System::isValid();
      #ifndef SIMP_NOPAIR
      UTIL_CHECK(pairPotentialPtr_);
      pairPotential().cellList().isValid();
      #endif
      return true;
   }

}
