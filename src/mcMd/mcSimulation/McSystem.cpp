#ifndef MCMD_MC_SYSTEM_CPP
#define MCMD_MC_SYSTEM_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "McSystem.h"
#include "McSimulation.h"
#include <mcMd/simulation/stress.h>
#include <mcMd/species/Species.h>
#include <util/ensembles/BoundaryEnsemble.h>
#include <util/ensembles/EnergyEnsemble.h>

#ifndef INTER_NOPAIR
#include <mcMd/potentials/pair/McPairPotential.h>
#include <mcMd/potentials/pair/PairFactory.h>
#endif
#include <mcMd/potentials/bond/BondPotential.h>
#ifdef INTER_ANGLE
#include <mcMd/potentials/angle/AnglePotential.h>
#endif
#ifdef INTER_DIHEDRAL
#include <mcMd/potentials/dihedral/DihedralPotential.h>
#endif
#ifdef MCMD_LINK
#include <mcMd/links/LinkMaster.h>
#endif
#ifdef INTER_EXTERNAL
#include <mcMd/potentials/external/ExternalPotential.h>
#include <mcMd/potentials/external/ExternalPotential.h>
#endif
#ifdef INTER_TETHER
#include <mcMd/tethers/TetherMaster.h>
#include <mcMd/potentials/tether/TetherPotential.h>
#endif

#ifdef MCMD_PERTURB
#include <mcMd/perturb/mcSystem/McPerturbationFactory.h>
#endif

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

   /*
   * Constructor.
   */
   McSystem::McSystem()
    :
      System(),
      #ifndef INTER_NOPAIR
      pairPotentialPtr_(0),
      #endif
      bondPotentialPtr_(0)
      #ifdef INTER_ANGLE
      , anglePotentialPtr_(0)
      #endif
      #ifdef INTER_DIHEDRAL
      , dihedralPotentialPtr_(0)
      #endif
      #ifdef MCMD_LINK
      , linkPotentialPtr_(0)
      #endif
      #ifdef INTER_EXTERNAL
      , externalPotentialPtr_(0)
      #endif
      #ifdef INTER_TETHER
      , tetherPotentialPtr_(0)
      #endif
   { setClassName("McSystem"); }

   /*
   * Destructor.
   */
   McSystem::~McSystem()
   {
      #ifndef INTER_NOPAIR
      if (pairPotentialPtr_) delete pairPotentialPtr_;
      #endif
      if (bondPotentialPtr_) delete bondPotentialPtr_;
      #ifdef INTER_ANGLE
      if (anglePotentialPtr_) delete anglePotentialPtr_;
      #endif
      #ifdef INTER_DIHEDRAL
      if (dihedralPotentialPtr_) delete dihedralPotentialPtr_;
      #endif
      #ifdef MCMD_LINK
      if (linkPotentialPtr_) delete linkPotentialPtr_;
      #endif
      #ifdef INTER_EXTERNAL
      if (externalPotentialPtr_) delete externalPotentialPtr_;
      #endif
      #ifdef INTER_TETHER
      if (tetherPotentialPtr_) delete tetherPotentialPtr_;
      #endif
   }

   /*
   * Read parameters from file.
   */
   void McSystem::readParameters(std::istream &in)
   {
      allocateMoleculeSets();
      readFileMaster(in);
      readPotentialStyles(in);

      #ifndef INTER_NOPAIR
      pairPotentialPtr_ = pairFactory().mcFactory(pairStyle(), *this);
      if (pairPotentialPtr_ == 0) {
         UTIL_THROW("Failed attempt to create McPairPotential");
      }
      readParamComposite(in, *pairPotentialPtr_);
      #endif

      assert(bondPotentialPtr_ == 0);
      if (simulation().nBondType() > 0) {
         bondPotentialPtr_ = bondFactory().factory(bondStyle());
         if (bondPotentialPtr_ == 0) {
            UTIL_THROW("Failed attempt to create BondPotential");
         }
         readParamComposite(in, *bondPotentialPtr_);
      }

      #ifdef INTER_ANGLE
      assert(anglePotentialPtr_ == 0);
      if (simulation().nAngleType() > 0) {
         anglePotentialPtr_ = angleFactory().factory(angleStyle());
         if (anglePotentialPtr_ == 0) {
            UTIL_THROW("Failed attempt to create AnglePotential");
         }
         readParamComposite(in, *anglePotentialPtr_);
      }
      #endif

      #ifdef INTER_DIHEDRAL
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

      #ifdef INTER_EXTERNAL
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

      #ifdef INTER_TETHER
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

      #ifndef INTER_NOPAIR
      pairPotentialPtr_ = pairFactory().mcFactory(pairStyle(), *this);
      if (pairPotentialPtr_ == 0) {
         UTIL_THROW("Failed attempt to create McPairPotential");
      }
      loadParamComposite(ar, *pairPotentialPtr_);
      #endif

      assert(bondPotentialPtr_ == 0);
      if (simulation().nBondType() > 0) {
         bondPotentialPtr_ = bondFactory().factory(bondStyle());
         if (bondPotentialPtr_ == 0) {
            UTIL_THROW("Failed attempt to create BondPotential");
         }
         loadParamComposite(ar, *bondPotentialPtr_);
      }

      #ifdef INTER_ANGLE
      assert(anglePotentialPtr_ == 0);
      if (simulation().nAngleType() > 0) {
         anglePotentialPtr_ = angleFactory().factory(angleStyle());
         if (anglePotentialPtr_ == 0) {
            UTIL_THROW("Failed attempt to create AnglePotential");
         }
         loadParamComposite(ar, *anglePotentialPtr_);
      }
      #endif

      #ifdef INTER_DIHEDRAL
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

      #ifdef INTER_EXTERNAL
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

      #ifdef INTER_TETHER
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
      #ifndef INTER_NOPAIR 
      assert(pairPotentialPtr_);
      pairPotentialPtr_->save(ar); 
      #endif
      if (simulation().nBondType() > 0) {
         assert(bondPotentialPtr_);
         bondPotentialPtr_->save(ar); 
      }
      #ifdef INTER_ANGLE
      if (simulation().nAngleType() > 0) {
         assert(anglePotentialPtr_);
         anglePotentialPtr_->save(ar); 
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (simulation().nDihedralType() > 0) {
         assert(dihedralPotentialPtr_);
         dihedralPotentialPtr_->save(ar); 
      }
      #endif
      #ifdef MCMD_LINK
      if (simulation().nLinkType() > 0) {
         saveLinkMaster(ar);
         assert(linkPotentialPtr_);
         linkPotentialPtr_->save(ar); 
      }
      #endif
      #ifdef INTER_EXTERNAL
      assert(externalPotentialPtr_ == 0);
      if (simulation().hasExternal() > 0) {
         assert(externalPotentialPtr_);
         externalPotentialPtr_->save(ar); 
      }
      #endif
      #ifdef INTER_TETHER
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
  
   /*
   * Read configuration from a specific input stream.
   */
   void McSystem::readConfig(std::istream &in)
   {
      System::readConfig(in);
      #ifndef INTER_NOPAIR
      pairPotential().buildCellList();
      #endif
   }

   /* 
   * Load a System configuration from an archive.
   */
   void McSystem::loadConfig(Serializable::IArchive& ar)
   {  
      System::loadConfig(ar); 
      #ifndef INTER_NOPAIR
      pairPotential().buildCellList();
      #endif
   }

   // Energy Evaluators (including all components)

   /*
   * Return total potential energy for one Atom.
   */
   double McSystem::atomPotentialEnergy(const Atom &atom) const
   {
      double energy = 0;
      #ifndef INTER_NOPAIR
      energy += pairPotential().atomEnergy(atom);
      #endif
      if (hasBondPotential()) {
         energy += bondPotential().atomEnergy(atom);
      }
      #ifdef INTER_ANGLE
      if (hasAnglePotential()) {
         energy += anglePotential().atomEnergy(atom);
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (hasDihedralPotential()) {
         energy += dihedralPotential().atomEnergy(atom);
      }
      #endif
      #ifdef MCMD_LINK
      if (hasLinkPotential()) {
         energy += linkPotential().atomEnergy(atom);
      }
      #endif
      #ifdef INTER_EXTERNAL
      if (hasExternalPotential()) {
         energy += externalPotential().atomEnergy(atom);
      }
      #endif
      #ifdef INTER_TETHER
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

   // Pressure/Stress Evaluators (including all components)

   template <typename T>
   void McSystem::computeVirialStressImpl(T& stress) const
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
   * Return true if this McSystem is valid, or an throw Exception.
   */
   bool McSystem::isValid() const
   {
      System::isValid();
      #ifndef INTER_NOPAIR
      pairPotentialPtr_->cellList().isValid();
      #endif
      return true;
   }

   // -------------------------------------------------------------
   #ifdef MCMD_PERTURB

   /*
   * Return a pointer to a new default McPerturbationFactory.
   */
   Factory<Perturbation>* McSystem::newDefaultPerturbationFactory()
   {  return new McPerturbationFactory(*this); }

   #endif
   // ------------------------------------------------------------------

}
#endif
