#ifndef MCMD_MC_SYSTEM_H
#define MCMD_MC_SYSTEM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

#include <mcMd/simulation/System.h>     // base class
#include <mcMd/neighbor/CellList.h>     // member

namespace McMd
{

   using namespace Util;

   class Atom;
   #ifndef INTER_NOPAIR
   class McPairPotential;
   #endif
   class BondPotential;
   #ifdef INTER_ANGLE
   class AnglePotential;
   #endif
   #ifdef INTER_DIHEDRAL
   class DihedralPotential;
   #endif
   #ifdef INTER_EXTERNAL
   class ExternalPotential;
   #endif
   #ifdef INTER_TETHER
   class TetherPotential;
   #endif

   /**
   * A System for use in a Markov chain Monte Carlo simulation.
   *
   * McSystem provides methods to evaluate energies of individual
   * atoms, and other functions needed in MC simulations.
   *
   * \ingroup McMd_System_Module
   */
   class McSystem : public System
   {

   public:

      /// Constructor.
      McSystem();

      /// Destructor.
      virtual ~McSystem();

      /**
      * Read parameters from input file.
      *
      * This method first calls System::readParameters(in).
      *
      * Finally, it allocates a CellList, and builds the new CellList if an
      * initial configuration has been read previously by System::readParam.
      * The number of cells allocated for the CellList is determined by the
      * dimensions of maxBoundary and the maximum nonbonded cutoff distance.
      *
      * \param in file stream.
      */
      virtual void readParameters(std::istream &in);

      /**
      * Load parameters from an archive, without configuration.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save parameters to an archive, without configuration.
      *
      * \param ar output/saving archive
      */
      virtual void saveParameters(Serializable::OArchive &ar);

      /**
      * Read system configuration from file.
      *
      * Calls System::readConfig() and then builds CellList.
      *
      * \param in configuration file input stream
      */
      virtual void readConfig(std::istream& in);

      /**
      * Load the system configuration from an archive.
      *
      * \param ar input (loading) archive object.
      */
      virtual void loadConfig(Serializable::IArchive& ar);

      /// \name Energy and Stress calculators
      //@{

      /**
      * Calculate the total potential energy for one Atom.
      *
      * \param  atom Atom object of interest
      * \return potential energy of atom
      */
      double atomPotentialEnergy(const Atom& atom) const;

      /**
      * Return total potential energy of this System.
      */
      double potentialEnergy() const;

      /**
      * Compute total virial stress (excludes kinetic contribution).
      *
      * The virial stress is the sum of all pair and bonded stress components.
      * It excludes the kinetic stress.
      *
      * \sa McSystem::computeStress(T&) regarding allowed types for T.
      *
      * \param stress (output) pressure, pressures or Tensor.
      */
      template <typename T>
      void computeVirialStress(T& stress) const;

      /**
      * Compute total pressure (T=double), xyz pressures (T=Vector) or stress
      * (T=Tensor).
      *
      * In this and all other stress / pressure calculators, typename T must be
      * double to obtain total pressure, Util::Vector to obtain the diagonal
      * pressure components, or Util::Tensor to obtain the full stress tensor.
      * These templates work only for these three types. The parameter stress
      * contains the desired value when the function returns.
      *
      * The total stress is the sum of the virial stress arising from
      * interatomic forces and an ideal gas kinetic stress.
      *
      * \param stress (output) pressure, xyz pressures or stress Tensor.
      */
      template <typename T>
      void computeStress(T& stress) const;

      // ----------------------------------------------------------------------

      #ifndef INTER_NOPAIR
      /**
      * Return the McPairPotential by reference.
      */
      McPairPotential& pairPotential() const;
      #endif

      /**
      * Does a bond potential exist?.
      */
      bool hasBondPotential() const;

      /**
      * Return the BondPotential by reference.
      */
      BondPotential& bondPotential() const;

      #ifdef INTER_ANGLE
      /**
      * Does angle potential exist?.
      */
      bool hasAnglePotential() const;

      /**
      * Return AnglePotential by reference.
      */
      AnglePotential& anglePotential() const;
      #endif


      #ifdef INTER_DIHEDRAL
      /**
      * Does a dihedral potential exist?.
      */
      bool hasDihedralPotential() const;

      /**
      * Return the DihedralPotential by reference.
      */
      DihedralPotential& dihedralPotential() const;
      #endif

      #ifdef MCMD_LINK
      /**
      * Does a link potential exist?.
      */
      bool hasLinkPotential() const;

      /**
      * Return the McLinkPotential by reference.
      */
      BondPotential& linkPotential() const;
      #endif

      #ifdef INTER_EXTERNAL
      /**
      * Does an external potential exist?.
      */
      bool hasExternalPotential() const;

      /**
      * Return ExternalPotential by reference.
      */
      ExternalPotential& externalPotential() const;
      #endif

      #ifdef INTER_TETHER
      /**
      * Return the TetherPotential by reference.
      */
      TetherPotential& tetherPotential() const;
      #endif

      //@}

      /**
      * Return true if McSystem is valid, or throw Exception.
      */
      virtual bool isValid() const;

   protected:

      #ifdef MCMD_PERTURB
      /**
      * Return a pointer to a new McPerturbationFactory.
      */
      virtual Factory<Perturbation>* newDefaultPerturbationFactory();
      #endif

      #ifndef INTER_NOPAIR
      /**
      * Set the PairPotential
      *
      * \param pairPotentialPtr pointer to pair potential
      */
      void setPairPotential(McPairPotential* pairPotentialPtr);
      #endif

   private:

      #ifndef INTER_NOPAIR
      /// Array to hold neighbors returned by a CellList.
      mutable CellList::NeighborArray neighbors_;

      McPairPotential* pairPotentialPtr_;
      #endif

      BondPotential* bondPotentialPtr_;

      #ifdef INTER_ANGLE
      /// Pointer to an AnglePotential.
      AnglePotential* anglePotentialPtr_;
      #endif

      #ifdef INTER_DIHEDRAL
      /// Pointer to an DihedralPotential.
      DihedralPotential* dihedralPotentialPtr_;
      #endif

      #ifdef MCMD_LINK
      /// Pointer to the McLinkPotential.
      BondPotential* linkPotentialPtr_;
      #endif

      #ifdef INTER_EXTERNAL
      /// Pointer to an ExternalPotential.
      ExternalPotential* externalPotentialPtr_;
      #endif

      #ifdef INTER_TETHER
      /// Pointer to an TetherPotential.
      TetherPotential* tetherPotentialPtr_;
      #endif

      /*
      * Implementations of the explicit specializations of the public
      * stress calculators computeStress(T& ) etc. for T = double,
      * Util::Vector and Util::Tensor simply call these private method
      * templates. This allows a single template method to be used to
      * compile the 3 relevant explicit specializations of each stress
      * calculator, while keeping the template implementations out of
      * this header file.
      */

      template <typename T>
      void computeVirialStressImpl(T& stress) const;

      #ifdef MCMD_LINK
      template <typename T>
      void computeLinkStressImpl(T& stress) const;
      #endif

   };

   // Inline methods

   #ifndef INTER_NOPAIR
   /*
   * Return the McPairPotential by reference.
   */
   inline McPairPotential& McSystem::pairPotential() const
   {  return *pairPotentialPtr_; }
   #endif

   /*
   * Does a bond potential exist?
   */
   inline bool McSystem::hasBondPotential() const
   {  return bool(bondPotentialPtr_); }

   /*
   * Return the BondPotential by const reference.
   */
   inline BondPotential& McSystem::bondPotential() const
   {  return *bondPotentialPtr_; }

   #ifdef INTER_ANGLE
   /*
   * Does an angle potential exist?
   */
   inline bool McSystem::hasAnglePotential() const
   {  return bool(anglePotentialPtr_); }

   /*
   * Return angle potential by reference.
   */
   inline AnglePotential& McSystem::anglePotential() const
   {  return *anglePotentialPtr_; }
   #endif

   #ifdef INTER_DIHEDRAL
   /*
   * Does a dihedral potential exist?
   */
   inline bool McSystem::hasDihedralPotential() const
   {  return bool(dihedralPotentialPtr_); }

   /*
   * Return dihedral potential by reference.
   */
   inline DihedralPotential& McSystem::dihedralPotential() const
   {  return *dihedralPotentialPtr_; }
   #endif

   #ifdef MCMD_LINK
   /*
   * Does a link potential exist?
   */
   inline bool McSystem::hasLinkPotential() const
   {  return bool(linkPotentialPtr_); }

   /*
   * Return link potential by reference.
   */
   inline BondPotential& McSystem::linkPotential() const
   {  return *linkPotentialPtr_; }
   #endif

   #ifdef INTER_EXTERNAL
   /*
   * Does an external potential exist?
   */
   inline bool McSystem::hasExternalPotential() const
   {  return bool(externalPotentialPtr_); }

   /*
   * Return external potential by reference.
   */
   inline ExternalPotential& McSystem::externalPotential() const
   {  return *externalPotentialPtr_; }
   #endif

   #ifdef INTER_TETHER
   /*
   * Return tether potential by reference.
   */
   inline TetherPotential& McSystem::tetherPotential() const
   {  return *tetherPotentialPtr_; }
   #endif

   #ifndef INTER_NOPAIR
   /*
   * Set the pair potential
   */
   inline void McSystem::setPairPotential(McPairPotential* pairPotentialPtr)
   {  pairPotentialPtr_ = pairPotentialPtr; }
   #endif

}
#endif
