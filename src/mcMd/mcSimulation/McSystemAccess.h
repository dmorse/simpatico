#ifndef MCMD_MC_SYSTEM_ACCESS_H
#define MCMD_MC_SYSTEM_ACCESS_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/simulation/SubSystem.h>  // base class

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
   * Provides faster access to an McSystem.
   *
   * \ingroup McMd_McMove_Module
   */
   class McSystemAccess : public SubSystem
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param mcSystem parent McSystem
      */
      McSystemAccess(McSystem& mcSystem);
   
      /**
      * Destructor.
      */
      virtual ~McSystemAccess();
   
   protected:

      #ifndef INTER_NOPAIR
      /**
      * Get the McPairPotential.
      */
      McPairPotential& pairPotential() const;
      #endif
   
      #ifdef INTER_BOND
      /**
      * Get the BondPotential.
      */
      BondPotential& bondPotential() const;
      #endif
   
      #ifdef INTER_ANGLE
      /**
      * Get the AnglePotential. 
      */
      AnglePotential& anglePotential() const;
      #endif
   
      #ifdef INTER_DIHEDRAL
      /**
      * Get the DihedralPotential.
      */
      DihedralPotential& dihedralPotential() const;
      #endif

      #ifdef INTER_EXTERNAL
      /**
      * Get the ExternalPotential.
      */
      ExternalPotential& externalPotential() const;
      #endif

   private:
 
      #ifndef INTER_NOPAIR
      McPairPotential* pairPotentialPtr_;
      #endif

      #ifdef INTER_BOND
      BondPotential* bondPotentialPtr_;
      #endif
 
      #ifdef INTER_ANGLE
      AnglePotential* anglePotentialPtr_;
      #endif

      #ifdef INTER_DIHEDRAL
      DihedralPotential* dihedralPotentialPtr_;
      #endif

      #ifdef INTER_EXTERNAL
      ExternalPotential* externalPotentialPtr_;
      #endif

   };

   // Inline methods

   #ifndef INTER_NOPAIR
   /*
   * Get McPairPotential of McSystem.
   */
   inline McPairPotential& McSystemAccess::pairPotential() const
   {  return *pairPotentialPtr_; }
   #endif

   #ifdef INTER_BOND
   /*
   * Get the BondPotential.
   */
   inline BondPotential& McSystemAccess::bondPotential() const
   {  return *bondPotentialPtr_; }
   #endif

   #ifdef INTER_ANGLE
   /*
   * Get AnglePotential. 
   */
   inline AnglePotential& McSystemAccess::anglePotential() const
   {  return *anglePotentialPtr_; }
   #endif

   #ifdef INTER_DIHEDRAL
   /*
   * Get DihedralPotential.
   */
   inline DihedralPotential& McSystemAccess::dihedralPotential() const
   {  return *dihedralPotentialPtr_; }
   #endif

   #ifdef INTER_EXTERNAL
   /*
   * Get ExternalPotential.
   */
   inline ExternalPotential& McSystemAccess::externalPotential() const
   {  return *externalPotentialPtr_; }
   #endif

}
#endif
