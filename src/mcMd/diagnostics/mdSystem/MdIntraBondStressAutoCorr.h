#ifndef MCMD_MD_INTRA_BOND_STRESS_AUTO_CORR_H
#define MCMD_MD_INTRA_BOND_STRESS_AUTO_CORR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/system/IntraBondStressAutoCorr.h> // base template
#include <mcMd/mdSimulation/MdSystem.h>                      // parameter

namespace McMd
{

   using namespace Util;


   /**
   * Autocorrelation for bond stress of a molecule.
   *
   * \ingroup McMd_Diagnostic_Module
   */
   class MdIntraBondStressAutoCorr : public IntraBondStressAutoCorr<MdSystem>
   {
   
   public:
  
      /**
      * Constructor.
      *
      * \param system reference to parent System
      */
      MdIntraBondStressAutoCorr(MdSystem &system);

      /**
      * Destructor.
      */
      virtual ~MdIntraBondStressAutoCorr();
  
      using IntraBondStressAutoCorr<MdSystem>::readParam;
      using IntraBondStressAutoCorr<MdSystem>::setup;
      using IntraBondStressAutoCorr<MdSystem>::sample;
      using IntraBondStressAutoCorr<MdSystem>::output;

      /**
      * Serialize to/from an archive. 
      *
      * \param ar      saving or loading archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void 
   MdIntraBondStressAutoCorr::serialize(Archive& ar, const unsigned int version)
   {  IntraBondStressAutoCorr<MdSystem>::serialize(ar, version); }

}
#endif
