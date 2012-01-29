#ifndef SYSTEM_DIAGNOSTIC_H
#define SYSTEM_DIAGNOSTIC_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/Diagnostic.h>
#include <mcMd/simulation/Simulation.h>
#include <util/global.h>                   // assert in inline function

namespace McMd
{

   using namespace Util;

   /**
   * Template for Diagnostic associated with one System.
   *
   * The FileMaster associated with a SystemDiagnostic is the one
   * used by the parent System.
   *
   * \ingroup Diagnostic_Module
   */
   template <class SystemType>
   class SystemDiagnostic : public Diagnostic 
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param system parent System
      */
      SystemDiagnostic(SystemType& system);
  
      /**
      * Destructor.
      */
      virtual ~SystemDiagnostic();
  
   protected:
  
      /** 
      * Return reference to parent system.
      */
      SystemType& system();

      #if 0
      bool hasBondPotential() const;

      #ifdef MCMD_ANGLE
      // Does an angle potential exist?
      bool hasAnglePotential() const;
      #endif

      #ifdef MCMD_DIHEDRAL
      // Does a dihedral potential exist?
      bool hasDihedralPotential() const;
      #endif

      #ifdef MCMD_LINK
      // Does a link potential exist?
      bool hasLinkPotential() const;
      #endif

      #ifdef MCMD_EXTERNAL
      // Does an external potential exist?
      bool hasExternalPotential() const;
      #endif

      #ifdef MCMD_TETHER
      // Does a tether potential exist?
      bool hasTetherPotential() const;
      #endif
      #endif

   private:

      /// Pointer to parent system.
      SystemType* systemPtr_;

      #if 0
      // Does a bond potential exist?
      bool hasBondPotential_;

      #ifdef MCMD_ANGLE
      // Does an angle potential exist?
      bool hasAnglePotential_;
      #endif

      #ifdef MCMD_DIHEDRAL
      // Does a dihedral potential exist?
      bool hasDihedralPotential_;
      #endif

      #ifdef MCMD_LINK
      // Does a link potential exist?
      bool hasLinkPotential_;
      #endif

      #ifdef MCMD_EXTERNAL
      // Does an external potential exist?
      bool hasExternalPotential_;
      #endif

      #ifdef MCMD_TETHER
      // Does a tether potential exist?
      bool hasTetherPotential_;
      #endif
      #endif

   };

   /*
   * Constructor.
   */
   template <class SystemType>
   SystemDiagnostic<SystemType>::SystemDiagnostic(SystemType& system) 
    : Diagnostic(),
      systemPtr_(&system)
      #if 0
      , hasBondPotential_(system.simulation().nBondType() > 0)
      #ifdef MCMD_ANGLE
      , hasAnglePotential_(system.simulation().nAngleType() > 0)
      #endif
      #ifdef MCMD_DIHEDRAL
      , hasDihedralPotential_(system.simulation().nDihedralType() > 0)
      #endif
      #ifdef MCMD_LINK
      , hasLinkPotential_(system.simulation().nLinkType() > 0)
      #endif
      #ifdef MCMD_EXTERNAL
      , hasExternalPotential_(system.simulation().hasExternal())
      #endif
      #ifdef MCMD_TETHER
      , hasTetherPotential_(system.simulation().hasTether())
      #endif
      #endif
   {  setFileMaster(system.fileMaster()); }

   /*
   * Destructor.
   */
   template <class SystemType>
   SystemDiagnostic<SystemType>::~SystemDiagnostic() 
   {}

   /*
   * Return reference to parent system.
   */
   template <class SystemType>
   inline SystemType& SystemDiagnostic<SystemType>::system()
   {
      assert(systemPtr_);
      return *systemPtr_; 
   }

   #if 0

   /*
   * Does a bond potential exist?
   */
   template <class SystemType>
   inline bool SystemDiagnostic<SystemType>::hasBondPotential() const
   {  return hasBondPotential_; }

   #ifdef MCMD_ANGLE
   /*
   * Does an angle potential exist?
   */
   template <class SystemType>
   inline bool SystemDiagnostic<SystemType>::hasAnglePotential() const
   {  return hasAnglePotential_; }
   #endif

   #ifdef MCMD_DIHEDRAL
   /// Does a dihedral potential exist?
   template <class SystemType>
   inline bool SystemDiagnostic<SystemType>::hasDihedralPotential() const
   {  return hasDihedralPotential_; }
   #endif

   #ifdef MCMD_LINK
   /// Does a link potential exist?
   template <class SystemType>
   inline bool SystemDiagnostic<SystemType>::hasLinkPotential() const
   { return hasLinkPotential_; }
   #endif

   #ifdef MCMD_EXTERNAL
   /// Does an external potential exist?
   template <class SystemType>
   inline bool SystemDiagnostic<SystemType>::hasExternalPotential() const
   { return hasExternalPotential_; }
   #endif

   #ifdef MCMD_TETHER
   /// Does a tether potential exist?
   template <class SystemType>
   inline bool SystemDiagnostic<SystemType>::hasTetherPotential() const
   { return hasTetherPotential_; }
   #endif

   #endif // if 0

}
#endif
