#ifndef MCMD_MD_INTRA_BOND_TENSOR_AUTO_CORR_H
#define MCMD_MD_INTRA_BOND_TENSOR_AUTO_CORR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/system/IntraBondTensorAutoCorr.h> // base template
#include <mcMd/mdSimulation/MdSystem.h>                      // parameter

namespace McMd
{

   using namespace Util;


   /**
   * Autocorrelation for bond orientation of of a molecule. 
   *
   * See class template IntraBondTensorAutoCorr for further information.
   *
   * \sa McMd::IntraBondTensorAutoCorr
   *
   * \ingroup McMd_Analyzer_Module
   */
   class MdIntraBondTensorAutoCorr : public IntraBondTensorAutoCorr<MdSystem>
   {
   
   public:
  
      /**
      * Constructor.
      *
      * \param system reference to parent System
      */
      MdIntraBondTensorAutoCorr(MdSystem &system);

      /**
      * Destructor.
      */
      virtual ~MdIntraBondTensorAutoCorr();
  
      using IntraBondTensorAutoCorr<MdSystem>::readParam;
      using IntraBondTensorAutoCorr<MdSystem>::setup;
      using IntraBondTensorAutoCorr<MdSystem>::sample;
      using IntraBondTensorAutoCorr<MdSystem>::output;
      using IntraBondTensorAutoCorr<MdSystem>::loadParameters;
      using IntraBondTensorAutoCorr<MdSystem>::save;

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
   MdIntraBondTensorAutoCorr::serialize(Archive& ar, const unsigned int version)
   {  IntraBondTensorAutoCorr<MdSystem>::serialize(ar, version); }

}
#endif
