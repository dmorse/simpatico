#ifndef MCMD_LINEAR_ROUSE_AUTO_CORR_H
#define MCMD_LINEAR_ROUSE_AUTO_CORR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, Jian Qin and David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>   // base class template
#include <mcMd/simulation/System.h>              // base class template parameter
#include <util/accumulators/AutoCorrArray.h>     // member template 
#include <util/space/Vector.h>                   // member template parameter
#include <util/containers/DArray.h>              // member template

#include <util/archives/serialize.h>             // used in method template
#include <util/global.h>                         // used in method template


namespace McMd
{

   using namespace Util;

   class Species;


   /**
   * Autocorrelation for Rouse mode coefficients of a linear molecule.
   *
   * \ingroup Diagnostic_Module
   */
   class LinearRouseAutoCorr : public SystemDiagnostic<System>
   {
   
   public:
  
      /**
      * Constructor.
      *
      * \param system reference to parent System
      */
      LinearRouseAutoCorr(System &system);

      /**
      * Destructor.
      */
      virtual ~LinearRouseAutoCorr();
  
      /** 
      * Read parameters from file.
      *
      * \param in input parameter stream
      */
      virtual void readParam(std::istream& in);
  
      /** 
      * Set number of molecules, initialize eigenvector, and clear accumulator.
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

   private:

      /// Output file stream.
      std::ofstream outputFile_;

      /// Statistical accumulator.
      AutoCorrArray<Vector, double>  accumulator_;
   
      /// Array of Rouse mode coefficients, one per molecule of species.
      DArray<Vector> data_;
   
      /// Array of coefficients for projecting to the Rouse mode.
      DArray<double> projector_;
   
      /// Pointer to relevant Species.
      Species* speciesPtr_;
   
      /// Index of relevant Species.
      int      speciesId_;
   
      /// Number of molecules in the species (must not change).
      int      nMolecule_;

      /// Number of atoms in the species (must not change).
      int      nAtom_;
 
      /// Index to Rouse mode.
      int      p_;

      /// Maximum length of each sequence in AutoCorrArray.
      int      capacity_;
   
      /// Has readParam been called?
      bool    isInitialized_;

   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void LinearRouseAutoCorr::serialize(Archive& ar, const unsigned int version)
   {
      if (!isInitialized_) {
         UTIL_THROW("Error: Object not initialized.");
      }

      // Data not set by readParam
      ar & accumulator_;
      ar & nMolecule_;

      // Data set by readParam (check consistency).
      serializeCheck(ar, speciesId_, "speciesId");
      serializeCheck(ar, nAtom_, "nAtom");
      serializeCheck(ar, capacity_, "capacity");
   }

}
#endif
