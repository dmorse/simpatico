#ifndef MCMD_MD_SYSTEM_H
#define MCMD_MD_SYSTEM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/simulation/System.h>     // base class

#include <util/global.h>

namespace Util { 
   template <typename T> class Factory; 
}

namespace McMd
{

   using namespace Util;

   // Forward declarations

   #ifndef INTER_NOPAIR
   class MdPairPotential;
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

   class MdIntegrator;
   class McSystem;

   /**
   * A System for Molecular Dynamics simulation.
   *
   * An MdSystem is a System that also has:
   *
   * - Associated Potential objects (MdPairPotential, BondPotential, etc.)
   *
   * - Associated MdIntegrator and MdIntegratorFactory objects.
   *
   * It provides methods to compute forces, total energy, and stress.
   *
   * \ingroup McMd_System_Module
   */
   class MdSystem : public System
   {
   
   public:

      /**   
      * Constructor.
      */
      MdSystem();

      /** 
      * Constructor, copy from McSystem.
      *
      * Used to create child MdSystem for Hybrid MD/MC move.
      *
      * \param system System object to be cloned.
      */
      MdSystem(McSystem& system);

      /** 
      * Destructor.
      */
      virtual ~MdSystem();

      /**
      * Read parameters from input file.
      *
      * This method calls System::readParameters(in), reads potentials, and 
      * allocates the pairList.
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
      * This calls System::readConfig(), followed by
      * pairPotential().buildPairList() and calculateForces().
      * 
      * \param in configuration file input stream
      */
      virtual void readConfig(std::istream& in);

      /**
      * Load the MdSystem configuration from an archive.
      *
      * This calls System::loadConfig(ar), followed by
      * pairPotential().buildPairList() and calculateForces().
      *
      * \param ar input (loading) archive object.
      */
      virtual void loadConfig(Serializable::IArchive& ar);

      /// \name Force, Energy and Stress calculators
      //@{

      /**
      * Compute all forces in this System.
      *
      * This method zero all forces and then calls the addForces()
      * method of each of the potential objects (pair, bond, etc.). 
      * The McPairPotential::addForces() method updates the pair 
      * list if necessary before calculating pair forces. On exit,
      * all atomic forces are updated to values corresponding to
      * current positions.
      */
      void calculateForces();

      /**
      * Compute and return total kinetic energy.
      */
      double kineticEnergy() const;

      /**
      * Compute and return total potential energy.
      */
      double potentialEnergy();

      /**
      * Compute total pressure (T=double), xyz pressures (T=Vector) or stress 
      * (T=Tensor).
      *
      * In this and all other stress / pressure calculator method templates, 
      * typename T must be double to obtain total pressure, Util::Vector to 
      * obtain the x, y, and z diagonal pressure components, or Util::Tensor 
      * to obtain the full stress tensor.  These templates work only for 
      * these three types. The function parameter "stress" is set to the
      * desired value.
      *
      * The total stress is the sum of the virial stress, from interatomic
      * forces, and a kinetic stress arising from the atomic velocities.
      *
      * \param stress (output) pressure, xyz pressures or stress Tensor.
      */
      template <typename T>
      void computeStress(T& stress) const;

      /**
      * Compute total virial stress, from all forces.
      *
      * The virial stress is the sum of all pair and bonded stress components.
      * It excludes the kinetic stress.
      *
      * \sa MdSystem::computeStress(T&) regarding allowed types for T.
      *
      * \param stress (output) pressure, pressures or Tensor.
      */
      template <typename T>
      void computeVirialStress(T& stress) const;

      /**
      * Compute kinetic stress mvv, arising from velocities.
      *
      * \sa MdSystem::computeStress(T&) regarding allowed types for T.
      *
      * \param stress (output) pressure, pressures or Tensor.
      */
      template <typename T>
      void computeKineticStress(T& stress) const;

      //@}
      /// \name Potential energy accessors
      //@{

      #ifndef INTER_NOPAIR
      /**
      * Return MdPairPotential by reference.
      */
      MdPairPotential& pairPotential() const;
      #endif
   
      /**
      * Does a bond potential exist?.
      */
      bool hasBondPotential() const;

      /**
      * Return BondPotential by reference.
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
      * Return DihedralPotential by reference.
      */
      DihedralPotential& dihedralPotential() const;
      #endif

      #ifdef MCMD_LINK
      /**
      * Does a link potential exist?.
      */
      bool hasLinkPotential() const;

      /**
      * Return link potential by reference.
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
      * Return TetherPotential by reference.
      */
      TetherPotential& tetherPotential() const;
      #endif

      //@}
      /// \name Mutators
      //@{

      /** 
      * Set all atomic forces to zero.
      */
      void setZeroForces();

      /**
      * Set all atomic velocities to zero.
      */
      void setZeroVelocities();

      /**
      * Set all velocities to Boltzmann distributed random values.
      *
      * \param temperature temperature used in Boltzmann distribution
      */
      void setBoltzmannVelocities(double temperature);

      /**
      * Subtract average velocity from all atomic velocities.
      *
      * Upon return, the system has a vanishing center-of-mass momentum.
      *
      * \return value of drift velocity before subtraction
      */
      Vector removeDriftVelocity();

      /**
      * Shift all atoms into primary cell.
      *
      * Shifts all atom positions into the primary image of the periodic 
      * boundary. 
      */
      void shiftAtoms();

      //@}
      /// \name Accessors
      //@{

      /**
      * Return the MdIntegrator by reference.
      */
      MdIntegrator& mdIntegrator();

      /**
      * Return the MdIntegrator Factory by reference.
      */
      Factory<MdIntegrator>& mdIntegratorFactory();

      /**
      * Return true if valid, or throw Exception.
      */
      virtual bool isValid() const;

      //@}

   protected:

      /**
      * Return a pointer to a new MdConfigIo object.
      */
      virtual ConfigIo* newDefaultConfigIo();

   private:
      
      #ifndef INTER_NOPAIR 
      /// Pointer to an MdPairPotential. 
      MdPairPotential* pairPotentialPtr_;
      #endif
  
      /// Pointer to an BondPotential. 
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
      /// Pointer to the MdLinkPotential. 
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

      /// Pointer to an MdIntegrator.
      MdIntegrator *mdIntegratorPtr_;

      /// Pointer to an MdIntegratorFactory.
      Factory<MdIntegrator> *mdIntegratorFactoryPtr_;

      /// Did this class create the MdIntegratorFactory?
      bool createdMdIntegratorFactory_;

      /*
      * Implementations of the explicit specializations of the public 
      * stress calculators computeVirialStress(T& ) etc. for T = double,
      * Util::Vector and Util::Tensor simply call these private method 
      * templates. This allows a single private template method to be 
      * used to compile the 3 relevant explicit specializations of each 
      * stress calculator, while allowing the template implementations 
      * to be defined in MdSystem.cpp, rather than in this header file. 
      */

      template <typename T>
      void computeKineticStressImpl(T& stress) const;

      template <typename T>
      void computeVirialStressImpl(T& stress) const;

   };

   // Inline functions
  
   #ifndef INTER_NOPAIR 
   /*
   * Return PairPotential by reference.
   */
   inline MdPairPotential& MdSystem::pairPotential() const
   {
      assert(pairPotentialPtr_);
      return *pairPotentialPtr_; 
   }
   #endif

   /*
   * Does a bond potential exist?
   */
   inline bool MdSystem::hasBondPotential() const
   {
      return bool(bondPotentialPtr_); 
   }

   /*
   * Return bond potential by reference.
   */
   inline BondPotential& MdSystem::bondPotential() const
   {  
      assert(bondPotentialPtr_);  
      return *bondPotentialPtr_; 
   }

   #ifdef INTER_ANGLE
   /*
   * Does an angle potential exist?
   */
   inline bool MdSystem::hasAnglePotential() const
   {
      return bool(anglePotentialPtr_); 
   }

   /*
   * Return angle potential by reference.
   */
   inline AnglePotential& MdSystem::anglePotential() const
   {
      assert(anglePotentialPtr_);  
      return *anglePotentialPtr_; 
   }
   #endif

   #ifdef INTER_DIHEDRAL
   /*
   * Does a dihedral potential exist?
   */
   inline bool MdSystem::hasDihedralPotential() const
   {
      return bool(dihedralPotentialPtr_); 
   }

   /*
   * Return dihedral potential by reference.
   */
   inline DihedralPotential& MdSystem::dihedralPotential() const
   {  
      assert(dihedralPotentialPtr_);  
      return *dihedralPotentialPtr_; 
   }
   #endif

   #ifdef MCMD_LINK
   /*
   * Does a link potential exist?
   */
   inline bool MdSystem::hasLinkPotential() const
   {
      return bool(linkPotentialPtr_); 
   }

   /*
   * Return link potential by reference.
   */
   inline BondPotential& MdSystem::linkPotential() const
   {  
      assert(linkPotentialPtr_);  
      return *linkPotentialPtr_; 
   }
   #endif

   #ifdef INTER_EXTERNAL
   /*
   * Does an external potential exist?
   */
   inline bool MdSystem::hasExternalPotential() const
   {
      return bool(externalPotentialPtr_); 
   }

   /*
   * Return external potential by reference.
   */
   inline ExternalPotential& MdSystem::externalPotential() const
   {  
      assert(externalPotentialPtr_);  
      return *externalPotentialPtr_; 
   }
   #endif

   #ifdef INTER_TETHER
   /*
   * Return tether potential by reference.
   */
   inline TetherPotential& MdSystem::tetherPotential() const
   {  return *tetherPotentialPtr_; }
   #endif

   /*
   * Return the MdIntegrator by reference.
   */
   inline MdIntegrator& MdSystem::mdIntegrator()
   {
      assert(mdIntegratorPtr_);  
      return *mdIntegratorPtr_; 
   }

} 
#endif
