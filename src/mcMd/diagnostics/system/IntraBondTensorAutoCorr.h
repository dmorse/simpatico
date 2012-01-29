#ifndef INTRA_BOND_TENSOR_AUTO_CORR_H
#define INTRA_BOND_TENSOR_AUTO_CORR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>   // base class template
#include <util/accumulators/AutoCorrArray.h>     // member template 
#include <util/space/Tensor.h>                   // member template parameter
#include <util/containers/DArray.h>              // member template

#include <util/archives/serialize.h>             // used in method template
#include <util/global.h>                         // used in method template


namespace McMd
{

   using namespace Util;

   class Species;


   /**
   * Autocorrelation for bond stress of a molecule.
   *
   * \ingroup Diagnostic_Module
   */
   template <class SystemType>
   class IntraBondTensorAutoCorr : public SystemDiagnostic<SystemType>
   {
   
   public:
  
      /**
      * Constructor.
      *
      * \param system reference to parent System
      */
      IntraBondTensorAutoCorr(SystemType &system);

      /**
      * Destructor.
      */
      virtual ~IntraBondTensorAutoCorr();
  
      /** 
      * Read parameters from file.
      *
      * \param in input parameter stream
      */
      virtual void readParam(std::istream& in);
  
      /** 
      * Set number of molecules and clear accumulator.
      */
      virtual void initialize();
  
      /** 
      * Evaluate end-to-end vectors of all chains, add to ensemble.
      *
      * \param iStep counter for number of steps
      */
      virtual void sample(long iStep);
   
      /**
      * Output results after simulation is completed.
      */
      virtual void output();

      /**
      * Serialize to/from an archive. 
      *
      * \param ar      saving or loading archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);


   protected:

      using SystemDiagnostic<SystemType>::read;
      using SystemDiagnostic<SystemType>::readOutputFileName;
      using SystemDiagnostic<SystemType>::readInterval;
      using SystemDiagnostic<SystemType>::isAtInterval;
      using SystemDiagnostic<SystemType>::writeParam;
      using SystemDiagnostic<SystemType>::system;
      using SystemDiagnostic<SystemType>::outputFileName;
      using SystemDiagnostic<SystemType>::fileMaster;

   private:

      /// Output file stream.
      std::ofstream outputFile_;

      /// Statistical accumulator.
      AutoCorrArray<Tensor, double>  accumulator_;
 
      /// Array of stress values for different molecules. 
      DArray<Tensor> data_;
 
      /// Pointer to relevant Species.
      Species* speciesPtr_;
   
      /// Index of relevant Species.
      int speciesId_;
   
      /// Number of molecules in the species (must not change).
      int nMolecule_;

      /// Number of bonds / molecule in the species (must not change).
      int nBond_;
 
      /// Maximum length of each sequence in AutoCorrArray.
      int capacity_;
   
      /// Has readParam been called?
      bool isInitialized_;

   };

   /*
   * Serialize to/from an archive. 
   */
   template <class SystemType>
   template <class Archive>
   void IntraBondTensorAutoCorr<SystemType>::serialize(Archive& ar, const unsigned int version)
   {
      if (!isInitialized_) {
         UTIL_THROW("Error: Object not initialized.");
      }

      // Data not set by readParam
      ar & accumulator_;
      ar & nMolecule_;

      // Data set by readParam (check consistency).
      serializeCheck(ar, speciesId_, "speciesId");
      serializeCheck(ar, nBond_, "nBond");
      serializeCheck(ar, capacity_, "capacity");
   }

}
#endif
