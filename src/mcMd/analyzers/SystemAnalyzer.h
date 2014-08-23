#ifndef MCMD_SYSTEM_ANALYZER_H
#define MCMD_SYSTEM_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/Analyzer.h>
#include <mcMd/simulation/Simulation.h>
#include <util/global.h>                   // assert in inline function

namespace McMd
{

   using namespace Util;

   /**
   * Template for Analyzer associated with one System.
   *
   * The FileMaster associated with a SystemAnalyzer is the one
   * used by the parent System.
   *
   * \ingroup McMd_Analyzer_Module
   */
   template <class SystemType>
   class SystemAnalyzer : public Analyzer 
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param system parent System
      */
      SystemAnalyzer(SystemType& system);
  
      /**
      * Destructor.
      */
      virtual ~SystemAnalyzer();
  
   protected:
  
      /** 
      * Return reference to parent system.
      */
      SystemType& system();

      #if 0
      bool hasBondPotential() const;

      #ifdef INTER_ANGLE
      // Does an angle potential exist?
      bool hasAnglePotential() const;
      #endif

      #ifdef INTER_DIHEDRAL
      // Does a dihedral potential exist?
      bool hasDihedralPotential() const;
      #endif

      #ifdef MCMD_LINK
      // Does a link potential exist?
      bool hasLinkPotential() const;
      #endif

      #ifdef INTER_EXTERNAL
      // Does an external potential exist?
      bool hasExternalPotential() const;
      #endif

      #ifdef INTER_TETHER
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

      #ifdef INTER_ANGLE
      // Does an angle potential exist?
      bool hasAnglePotential_;
      #endif

      #ifdef INTER_DIHEDRAL
      // Does a dihedral potential exist?
      bool hasDihedralPotential_;
      #endif

      #ifdef MCMD_LINK
      // Does a link potential exist?
      bool hasLinkPotential_;
      #endif

      #ifdef INTER_EXTERNAL
      // Does an external potential exist?
      bool hasExternalPotential_;
      #endif

      #ifdef INTER_TETHER
      // Does a tether potential exist?
      bool hasTetherPotential_;
      #endif
      #endif

   };

   /*
   * Constructor.
   */
   template <class SystemType>
   SystemAnalyzer<SystemType>::SystemAnalyzer(SystemType& system) 
    : Analyzer(),
      systemPtr_(&system)
      #if 0
      , hasBondPotential_(system.simulation().nBondType() > 0)
      #ifdef INTER_ANGLE
      , hasAnglePotential_(system.simulation().nAngleType() > 0)
      #endif
      #ifdef INTER_DIHEDRAL
      , hasDihedralPotential_(system.simulation().nDihedralType() > 0)
      #endif
      #ifdef MCMD_LINK
      , hasLinkPotential_(system.simulation().nLinkType() > 0)
      #endif
      #ifdef INTER_EXTERNAL
      , hasExternalPotential_(system.simulation().hasExternal())
      #endif
      #ifdef INTER_TETHER
      , hasTetherPotential_(system.simulation().hasTether())
      #endif
      #endif
   {  setFileMaster(system.fileMaster()); }

   /*
   * Destructor.
   */
   template <class SystemType>
   SystemAnalyzer<SystemType>::~SystemAnalyzer() 
   {}

   /*
   * Return reference to parent system.
   */
   template <class SystemType>
   inline SystemType& SystemAnalyzer<SystemType>::system()
   {
      assert(systemPtr_);
      return *systemPtr_; 
   }

   #if 0

   /*
   * Does a bond potential exist?
   */
   template <class SystemType>
   inline bool SystemAnalyzer<SystemType>::hasBondPotential() const
   {  return hasBondPotential_; }

   #ifdef INTER_ANGLE
   /*
   * Does an angle potential exist?
   */
   template <class SystemType>
   inline bool SystemAnalyzer<SystemType>::hasAnglePotential() const
   {  return hasAnglePotential_; }
   #endif

   #ifdef INTER_DIHEDRAL
   /// Does a dihedral potential exist?
   template <class SystemType>
   inline bool SystemAnalyzer<SystemType>::hasDihedralPotential() const
   {  return hasDihedralPotential_; }
   #endif

   #ifdef MCMD_LINK
   /// Does a link potential exist?
   template <class SystemType>
   inline bool SystemAnalyzer<SystemType>::hasLinkPotential() const
   { return hasLinkPotential_; }
   #endif

   #ifdef INTER_EXTERNAL
   /// Does an external potential exist?
   template <class SystemType>
   inline bool SystemAnalyzer<SystemType>::hasExternalPotential() const
   { return hasExternalPotential_; }
   #endif

   #ifdef INTER_TETHER
   /// Does a tether potential exist?
   template <class SystemType>
   inline bool SystemAnalyzer<SystemType>::hasTetherPotential() const
   { return hasTetherPotential_; }
   #endif

   #endif // if 0

}
#endif
