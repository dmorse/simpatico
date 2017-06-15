#ifndef MCMD_MC_SYSTEM_INTERFACE_H
#define MCMD_MC_SYSTEM_INTERFACE_H

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

   class McSystem;
   class McPairPotential;
   class BondPotential;
   class AnglePotential;
   class DihedralPotential;
   class ExternalPotential;

   /**
   * An interface to an McSystem.
   *
   * \ingroup McMd_System_Module
   */
   class McSystemInterface : public SystemInterface
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param mcSystem parent McSystem
      */
      McSystemInterface(McSystem& mcSystem);
   
      /**
      * Destructor.
      */
      virtual ~McSystemInterface();
   
   protected:

      #ifndef SIMP_NOPAIR
      /**
      * Get the McPairPotential.
      */
      McPairPotential& pairPotential() const;
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
      McPairPotential* pairPotentialPtr_;
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
   * Get McPairPotential of McSystem.
   */
   inline McPairPotential& McSystemInterface::pairPotential() const
   {  return *pairPotentialPtr_; }
   #endif

   #ifdef SIMP_BOND
   /*
   * Get the BondPotential.
   */
   inline BondPotential& McSystemInterface::bondPotential() const
   {  return *bondPotentialPtr_; }
   #endif

   #ifdef SIMP_ANGLE
   /*
   * Get AnglePotential. 
   */
   inline AnglePotential& McSystemInterface::anglePotential() const
   {  return *anglePotentialPtr_; }
   #endif

   #ifdef SIMP_DIHEDRAL
   /*
   * Get DihedralPotential.
   */
   inline DihedralPotential& McSystemInterface::dihedralPotential() const
   {  return *dihedralPotentialPtr_; }
   #endif

   #ifdef SIMP_EXTERNAL
   /*
   * Get ExternalPotential.
   */
   inline ExternalPotential& McSystemInterface::externalPotential() const
   {  return *externalPotentialPtr_; }
   #endif

}
#endif
