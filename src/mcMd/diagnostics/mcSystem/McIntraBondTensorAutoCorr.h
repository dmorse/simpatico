#ifndef MCMD_MC_INTRA_BOND_TENSOR_AUTO_CORR_H
#define MCMD_MC_INTRA_BOND_TENSOR_AUTO_CORR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/system/IntraBondTensorAutoCorr.h> // base template
#include <mcMd/mcSimulation/McSystem.h>                      // parameter

namespace McMd
{

   using namespace Util;


   /**
   * Autocorrelation for bond stress of a molecule.
   *
   * \ingroup McMd_Diagnostic_Module
   */
   class McIntraBondTensorAutoCorr : public IntraBondTensorAutoCorr<McSystem>
   {
   
   public:
  
      /**
      * Constructor.
      *
      * \param system reference to parent System
      */
      McIntraBondTensorAutoCorr(McSystem &system);

      /**
      * Destructor.
      */
      virtual ~McIntraBondTensorAutoCorr();
  
      using IntraBondTensorAutoCorr<McSystem>::readParameters;
      using IntraBondTensorAutoCorr<McSystem>::loadParameters;
      using IntraBondTensorAutoCorr<McSystem>::save;
      using IntraBondTensorAutoCorr<McSystem>::setup;
      using IntraBondTensorAutoCorr<McSystem>::sample;
      using IntraBondTensorAutoCorr<McSystem>::output;

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
   McIntraBondTensorAutoCorr::serialize(Archive& ar, const unsigned int version)
   {  IntraBondTensorAutoCorr<McSystem>::serialize(ar, version); }

}
#endif
