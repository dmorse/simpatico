#ifndef MCMD_MD_SYSTEM_INTERFACE_H
#define MCMD_MD_SYSTEM_INTERFACE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/simulation/SystemInterface.h>  // base class

namespace McMd
{

   using namespace Util;

   class MdSystem;
   class MdPairPotential;
   class BondPotential;
   class AnglePotential;
   class DihedralPotential;
   class ExternalPotential;

   /**
   * An interface to an MdSystem.
   *
   * \ingroup McMd_System_Module
   */
   class MdSystemInterface : public SystemInterface
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param mcSystem parent MdSystem
      */
      MdSystemInterface(MdSystem& mcSystem);
   
      /**
      * Destructor.
      */
      virtual ~MdSystemInterface();
   
   protected:

      #ifndef SIMP_NOPAIR
      /**
      * Get the MdPairPotential.
      */
      MdPairPotential& pairPotential() const;
      #endif
   
      #ifdef SIMP_BOND
      /**
      * Get the BondPotential.
      */
      BondPotential& bondPotential() const;
      #endif
   
      #ifdef SIMP_ANGLE
      /**
      * Get the AnglePotential. 
      */
      AnglePotential& anglePotential() const;
      #endif
   
      #ifdef SIMP_DIHEDRAL
      /**
      * Get the DihedralPotential.
      */
      DihedralPotential& dihedralPotential() const;
      #endif

      #ifdef SIMP_EXTERNAL
      /**
      * Get the ExternalPotential.
      */
      ExternalPotential& externalPotential() const;
      #endif

   private:
 
      #ifndef SIMP_NOPAIR
      MdPairPotential* pairPotentialPtr_;
      #endif

      #ifdef SIMP_BOND
      BondPotential* bondPotentialPtr_;
      #endif
 
      #ifdef SIMP_ANGLE
      AnglePotential* anglePotentialPtr_;
      #endif

      #ifdef SIMP_DIHEDRAL
      DihedralPotential* dihedralPotentialPtr_;
      #endif

      #ifdef SIMP_EXTERNAL
      ExternalPotential* externalPotentialPtr_;
      #endif

   };

   // Inline methods

   #ifndef SIMP_NOPAIR
   /*
   * Get MdPairPotential of MdSystem.
   */
   inline MdPairPotential& MdSystemInterface::pairPotential() const
   {  return *pairPotentialPtr_; }
   #endif

   #ifdef SIMP_BOND
   /*
   * Get the BondPotential.
   */
   inline BondPotential& MdSystemInterface::bondPotential() const
   {  return *bondPotentialPtr_; }
   #endif

   #ifdef SIMP_ANGLE
   /*
   * Get AnglePotential. 
   */
   inline AnglePotential& MdSystemInterface::anglePotential() const
   {  return *anglePotentialPtr_; }
   #endif

   #ifdef SIMP_DIHEDRAL
   /*
   * Get DihedralPotential.
   */
   inline DihedralPotential& MdSystemInterface::dihedralPotential() const
   {  return *dihedralPotentialPtr_; }
   #endif

   #ifdef SIMP_EXTERNAL
   /*
   * Get ExternalPotential.
   */
   inline ExternalPotential& MdSystemInterface::externalPotential() const
   {  return *externalPotentialPtr_; }
   #endif

}
#endif
