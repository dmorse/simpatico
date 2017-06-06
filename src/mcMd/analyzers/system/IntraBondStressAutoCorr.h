#ifndef MCMD_INTRA_BOND_STRESS_AUTO_CORR_H
#define MCMD_INTRA_BOND_STRESS_AUTO_CORR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>   // base class template
#include <util/accumulators/AutoCorrArray.h>     // member template 
#include <util/space/Tensor.h>                   // member template parameter
#include <util/containers/DArray.h>              // member template

#include <util/archives/serialize.h>             // used in method template
#include <util/global.h>                         // used in method template

namespace Simp {
   class Species;
}

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /**
   * Autocorrelation for bond stress of a molecule.
   *
   * \ingroup McMd_Analyzer_McMd_Module
   */
   template <class SystemType>
   class IntraBondStressAutoCorr : public SystemAnalyzer<SystemType>
   {
   
   public:
  
      /**
      * Constructor.
      *
      * \param system reference to parent System
      */
      IntraBondStressAutoCorr(SystemType &system);

      /**
      * Destructor.
      */
      virtual ~IntraBondStressAutoCorr();
  
      /** 
      * Read parameters from file.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);
  
      /**
      * Load internal state from archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);
   
      /**
      * Save internal state to archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);
   
      /**
      * Serialize to/from an archive. 
      *
      * \param ar      saving or loading archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /** 
      * Set number of molecules and clear accumulator.
      */
      virtual void setup();
  
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

   protected:

      using SystemAnalyzer<SystemType>::read;
      using SystemAnalyzer<SystemType>::readInterval;
      using SystemAnalyzer<SystemType>::readOutputFileName;
      using SystemAnalyzer<SystemType>::loadParameter;
      using SystemAnalyzer<SystemType>::isAtInterval;
      using SystemAnalyzer<SystemType>::writeParam;
      using SystemAnalyzer<SystemType>::system;
      using SystemAnalyzer<SystemType>::outputFileName;
      using SystemAnalyzer<SystemType>::fileMaster;

   private:

      /// Output file stream.
      std::ofstream outputFile_;

      /// Statistical accumulator.
      AutoCorrArray<Tensor, double>  accumulator_;
 
      /// Array of stress values for different molecules (temporary)
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
   void IntraBondStressAutoCorr<SystemType>::serialize(Archive& ar, const unsigned int version)
   {
      Analyzer::serialize(ar, version);
      ar & speciesId_;
      ar & capacity_;

      ar & nBond_;
      ar & nMolecule_;
      ar & accumulator_;
   }

}
#endif
